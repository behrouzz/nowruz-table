"""
March equinox from DE441.

Definition
----------
The March equinox is the instant at which the Sun's **apparent
geocentric ecliptic longitude** — on the mean ecliptic of date, measured
from the true equinox of date (see ``ephemeris.apparent_solar_longitude``)
— passes through 0° while increasing.

Method
------
1. Approximate instant from the mean tropical year (only a starting
   point; it can be off by a day or more far from J2000).
2. Newton-like refinement using the DE441 longitude and the mean solar
   motion (~0.9856°/day).
3. Final root by Brent's method (``scipy.optimize.brentq``) on a
   bracket of ±0.01 day around the refined value, xtol 1e-9 day (~0.1 ms).

The function evaluated is the longitude wrapped into (-180°, +180°], so
the 360° -> 0° discontinuity becomes a smooth zero crossing.

The independent variable is TDB (the DE441 argument); the result is then
converted TDB -> TT -> UT1 (``time_scales``).  The TT/TDB instants depend
only on the ephemeris and the reduction models; the UT1 instant also
carries the ΔT uncertainty.
"""

import math
from dataclasses import dataclass
from functools import lru_cache

from scipy.optimize import brentq

from .config import (
    MAX_GREGORIAN_YEAR,
    MIN_GREGORIAN_YEAR,
    REFERENCE_MEAN_TIME_OFFSET_DAYS,
)
from .date_conversion import format_jd, jd_to_gregorian_datetime
from .delta_t import sigma_at_jd
from .ephemeris import apparent_solar_longitude
from .time_scales import tdb_to_tt, tt_to_ut1

# March equinox 2000 (TT) and mean tropical year: starting estimate only.
_JD_EQUINOX_2000 = 2451623.8159
_TROPICAL_YEAR = 365.24219
_MEAN_MOTION = 2.0 * math.pi / _TROPICAL_YEAR  # rad / day


def _wrapped_longitude(jd_tdb):
    """Apparent solar longitude in (-π, π] radians."""

    lon = apparent_solar_longitude(jd_tdb)["longitude"]
    return lon - 2.0 * math.pi if lon > math.pi else lon


def check_year(year):
    if not MIN_GREGORIAN_YEAR <= year <= MAX_GREGORIAN_YEAR:
        raise ValueError(
            f"Year {year} is outside the supported range "
            f"{MIN_GREGORIAN_YEAR} .. {MAX_GREGORIAN_YEAR} (astronomical)."
        )


@lru_cache(maxsize=None)
def march_equinox_tdb(year):
    """March equinox of astronomical Gregorian year ``year``, JD(TDB)."""

    check_year(year)

    t = _JD_EQUINOX_2000 + (year - 2000) * _TROPICAL_YEAR
    for _ in range(10):
        step = _wrapped_longitude(t) / _MEAN_MOTION
        t -= step
        if abs(step) < 1e-4:
            break

    half = 0.01
    while True:
        a, b = t - half, t + half
        fa, fb = _wrapped_longitude(a), _wrapped_longitude(b)
        if fa < 0.0 < fb:
            break
        half *= 2.0
        if half > 5.0:
            raise RuntimeError(f"Cannot bracket the equinox of {year}.")

    return brentq(_wrapped_longitude, a, b, xtol=1e-9, rtol=1e-15)


@dataclass(frozen=True)
class Equinox:
    """March equinox instant in the three time scales."""

    year: int              # astronomical Gregorian year
    jd_tdb: float
    jd_tt: float
    jd_ut1: float          # with the model ΔT (+ any offset requested)
    delta_t: float         # ΔT used, seconds
    delta_t_sigma: float   # estimated 1σ of the ΔT model, seconds

    @property
    def local_mean_time(self):
        """(GregorianDate, h, m, s) in local mean time at 52.5° E."""

        return jd_to_gregorian_datetime(
            self.jd_ut1, REFERENCE_MEAN_TIME_OFFSET_DAYS
        )

    def describe(self):
        return (
            f"March equinox {self.year}: "
            f"{format_jd(self.jd_tt, scale='TT')} | "
            f"{format_jd(self.jd_ut1, scale='UT1')} "
            f"(ΔT = {self.delta_t:.1f} s ± {self.delta_t_sigma:.0f} s) | "
            f"{format_jd(self.jd_ut1, REFERENCE_MEAN_TIME_OFFSET_DAYS, 'LMT 52.5E')}"
        )

    def __str__(self):
        return self.describe()


def march_equinox(year, delta_t_offset=0.0):
    """
    March equinox of astronomical Gregorian ``year`` as an ``Equinox``.

    ``delta_t_offset`` (s) is added to the model ΔT for the UT1 value.
    """

    jd_tdb = march_equinox_tdb(year)
    jd_tt = tdb_to_tt(jd_tdb)
    jd_ut1, dt = tt_to_ut1(jd_tt, delta_t_offset)

    return Equinox(
        year=year,
        jd_tdb=jd_tdb,
        jd_tt=jd_tt,
        jd_ut1=jd_ut1,
        delta_t=dt,
        delta_t_sigma=sigma_at_jd(jd_tt),
    )
