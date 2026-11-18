"""
DE441 access and apparent geocentric coordinates of the Sun.

Ephemeris
---------
The reduced DE441 SPK kernel holds three Chebyshev (type 2) chains,
each split into two segments at JD 2440432.5 (1969-07-30):

    0 -> 3    Earth-Moon barycentre (EMB) w.r.t. Solar System barycentre
    3 -> 399  Earth w.r.t. EMB
    0 -> 10   Sun w.r.t. Solar System barycentre

The segment covering the requested TDB is selected explicitly (jplephem's
``kernel[center, target]`` returns only one segment).  Positions are in
km in the ICRF axes, which DE441 realises; velocities are the analytic
derivatives of the Chebyshev polynomials (km/day).

Apparent place of the Sun (geocentric)
--------------------------------------
1. Earth(t) = EMB(t) + (Earth - EMB)(t), barycentric, at TDB t.
2. Light time: Sun(t - τ) with τ = |Sun(t - τ) - Earth(t)| / c, iterated.
   This yields the astrometric direction in the BCRS/ICRS axes.
3. Gravitational deflection: none (the deflecting body is the source).
4. Annual aberration: ERFA ``ab`` (relativistic, with Earth's
   barycentric velocity), giving the proper direction in the GCRS.
5. Frame bias + precession: ERFA ``ltpb`` — Vondrák, Capitaine &
   Wallace (2011, 2012) long-term model, valid for ±200 000 years and
   agreeing with IAU 2006 near J2000.  IAU 2006 (P03) precession is
   *not* used: its polynomials are only reliable within about ±1000 yr.
6. Nutation: IAU 2000A with IAU 2006 adjustments (``nut06a``), applied
   with the mean obliquity of the same Vondrák model.  See limitations.
7. Apparent ecliptic longitude: the vector is rotated from the true
   equator & equinox of date to the ecliptic by the *true* obliquity.
   Equivalently: longitude on the Vondrák **mean ecliptic of date**,
   measured from the **true equinox of date** (mean longitude + Δψ),
   which is the classical "apparent longitude" of the almanacs.

The diurnal parallax and diurnal aberration are not included (geocentric
place).  At apparent noon the parallax displaces the Sun in declination
only, and diurnal aberration changes the noon instant by ~0.02 s.

Limitations of the reduction (see README for magnitudes)
-------------------------------------------------------
* IAU 2000A nutation is a fit to modern data; its fundamental arguments
  are polynomials whose errors grow away from J2000.  Nutation itself
  is bounded (|Δψ| < 20"), so at ±4000 yr it is still good to a few
  arcseconds (1" of solar longitude ≈ 24 s of time).
* The CIO locator ``s06`` used for the hour angle is an IAU 2006 series;
  far from J2000 its error reaches some tens of arcseconds (a few
  seconds of time in the noon instant).
"""

import bisect
import math

import erfa
import numpy as np
from jplephem.spk import SPK

from .config import get_kernel_path
from .time_scales import tdb_to_tt, tt_to_tdb

AU_KM = 149597870.7
C_KM_PER_S = 299792.458
C_KM_PER_DAY = C_KM_PER_S * 86400.0

_MJD0 = 2400000.5
_REQUIRED_CHAINS = ((0, 3), (3, 399), (0, 10))


# ---------------------------------------------------------------------
# DE441 kernel
# ---------------------------------------------------------------------

class Ephemeris:
    """Explicit-segment access to the reduced DE441 kernel."""

    def __init__(self, path):
        self.path = path
        self._spk = SPK.open(path)
        self._chains = {}

        for segment in self._spk.segments:
            key = (segment.center, segment.target)
            self._chains.setdefault(key, []).append(segment)

        for key in _REQUIRED_CHAINS:
            if key not in self._chains:
                self._spk.close()
                raise ValueError(
                    f"Kernel {path} has no segment {key[0]} -> {key[1]}."
                )
            self._chains[key].sort(key=lambda s: s.start_jd)

        self._starts = {
            key: [s.start_jd for s in segs]
            for key, segs in self._chains.items()
        }

    @property
    def coverage(self):
        """(first, last) TDB Julian Date covered by every required chain."""

        start = max(self._chains[k][0].start_jd for k in _REQUIRED_CHAINS)
        end = min(self._chains[k][-1].end_jd for k in _REQUIRED_CHAINS)
        return start, end

    def segment(self, center, target, jd_tdb):
        """The segment of chain center -> target that covers ``jd_tdb``."""

        key = (center, target)
        segments = self._chains[key]
        i = bisect.bisect_right(self._starts[key], jd_tdb) - 1

        if i >= 0 and segments[i].start_jd <= jd_tdb <= segments[i].end_jd:
            return segments[i]

        raise ValueError(
            f"No DE441 segment {center} -> {target} covers "
            f"JD(TDB) {jd_tdb:.6f}."
        )

    def position(self, center, target, jd_tdb):
        """Position (km) of target w.r.t. center."""

        return self.segment(center, target, jd_tdb).compute(jd_tdb)

    def position_velocity(self, center, target, jd_tdb):
        """Position (km) and velocity (km/day) of target w.r.t. center."""

        seg = self.segment(center, target, jd_tdb)
        return seg.compute_and_differentiate(jd_tdb)

    def earth_barycentric(self, jd_tdb):
        """Barycentric position (km) and velocity (km/day) of the Earth."""

        emb_p, emb_v = self.position_velocity(0, 3, jd_tdb)
        off_p, off_v = self.position_velocity(3, 399, jd_tdb)
        return emb_p + off_p, emb_v + off_v

    def sun_barycentric(self, jd_tdb):
        """Barycentric position (km) of the Sun."""

        return self.position(0, 10, jd_tdb)

    def close(self):
        self._spk.close()


_ephemeris = None


def get_ephemeris():
    """The shared ``Ephemeris`` instance (opened lazily)."""

    global _ephemeris
    if _ephemeris is None:
        _ephemeris = Ephemeris(get_kernel_path())
    return _ephemeris


def close():
    """Close the DE441 kernel (it is reopened on the next use)."""

    global _ephemeris
    if _ephemeris is not None:
        _ephemeris.close()
        _ephemeris = None


# ---------------------------------------------------------------------
# Geocentric Sun
# ---------------------------------------------------------------------

def sun_geocentric_astrometric(jd_tdb):
    """
    Light-time-corrected geocentric vector of the Sun (km, ICRF axes).

    Returns (vector_km, earth_velocity_km_per_day, light_time_days).
    """

    eph = get_ephemeris()
    earth_p, earth_v = eph.earth_barycentric(jd_tdb)

    tau = 0.0
    for _ in range(6):
        vector = eph.sun_barycentric(jd_tdb - tau) - earth_p
        tau_new = math.sqrt(vector @ vector) / C_KM_PER_DAY
        if abs(tau_new - tau) < 1e-12:  # ~0.1 µs
            break
        tau = tau_new

    return vector, earth_v, tau


def sun_apparent_gcrs(jd_tdb):
    """
    Apparent (proper) geocentric direction of the Sun in the GCRS.

    Light time and annual aberration included.  Returns
    (unit_vector, distance_au, light_time_days).
    """

    vector, earth_v, tau = sun_geocentric_astrometric(jd_tdb)
    distance_km = math.sqrt(vector @ vector)
    natural = vector / distance_km

    beta = earth_v / C_KM_PER_DAY
    bm1 = math.sqrt(1.0 - beta @ beta)
    proper = erfa.ab(natural, beta, distance_km / AU_KM, bm1)

    return proper, distance_km / AU_KM, tau


# ---------------------------------------------------------------------
# Precession-nutation (long-term)
# ---------------------------------------------------------------------

def precession_nutation(jd_tt):
    """
    GCRS -> true equator & equinox of date.

    Vondrák long-term frame bias + precession, IAU 2000A/2006 nutation.
    Returns (rnpb, true_obliquity, dpsi, deps), angles in radians.
    """

    epj = erfa.epj(_MJD0, jd_tt - _MJD0)

    rbp = erfa.ltpb(epj)
    eps_mean = math.acos(
        float(np.clip(erfa.ltpequ(epj) @ erfa.ltpecl(epj), -1.0, 1.0))
    )
    dpsi, deps = erfa.nut06a(_MJD0, jd_tt - _MJD0)

    rnpb = erfa.numat(eps_mean, dpsi, deps) @ rbp
    return rnpb, eps_mean + deps, dpsi, deps


# ---------------------------------------------------------------------
# Apparent solar coordinates
# ---------------------------------------------------------------------

def apparent_solar_longitude(jd_tdb):
    """
    Apparent geocentric ecliptic coordinates of the Sun at TDB ``jd_tdb``.

    Longitude is on the mean ecliptic of date from the true equinox of
    date.  Returns dict with longitude (rad, [0, 2π)), latitude (rad),
    distance_au.
    """

    proper, distance_au, _ = sun_apparent_gcrs(jd_tdb)
    rnpb, eps_true, _, _ = precession_nutation(tdb_to_tt(jd_tdb))

    ecl = erfa.rx(eps_true, rnpb) @ proper
    lon, lat = erfa.c2s(ecl)

    return {
        "longitude": float(erfa.anp(lon)),
        "latitude": float(lat),
        "distance_au": distance_au,
    }


def apparent_solar_hour_angle(jd_ut1, jd_tt, longitude_rad):
    """
    Local apparent hour angle of the Sun (rad, in (-π, π]).

    ``jd_ut1`` drives the Earth rotation angle; ``jd_tt`` (the same
    instant in TT, i.e. UT1 + ΔT) drives the Sun's place and
    precession-nutation.  Uses the CIO-based relation

        H = ERA(UT1) + λ_east - α_CIO

    with α_CIO the apparent right ascension measured from the CIO on
    the CIP equator.  Polar motion is neglected (< 1", unknown in the
    past); the place is geocentric.
    """

    proper, _, _ = sun_apparent_gcrs(tt_to_tdb(jd_tt))
    rnpb, _, _, _ = precession_nutation(jd_tt)

    x, y = rnpb[2, 0], rnpb[2, 1]
    s = erfa.s06(_MJD0, jd_tt - _MJD0, x, y)
    cirs = erfa.c2ixys(x, y, s) @ proper
    ra_cio = math.atan2(cirs[1], cirs[0])

    era = erfa.era00(_MJD0, jd_ut1 - _MJD0)
    return float(erfa.anpm(era + longitude_rad - ra_cio))
