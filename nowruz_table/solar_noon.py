"""
Apparent solar noon at the reference meridian (52.5° E).

Apparent noon is the instant when the **local apparent hour angle** of
the Sun's geocentric apparent place is zero:

    H(t) = ERA(UT1) + λ - α_CIO(TT) = 0

* ERA — Earth Rotation Angle, a linear function of UT1 (IAU 2000).
* α_CIO — apparent right ascension of the Sun (light time + aberration,
  Vondrák precession, IAU 2000A nutation) measured from the CIO.
* TT = UT1 + ΔT for the Sun's place; the result is in UT1.

This includes the equation of time implicitly (it is *not* the mean-time
expression UT1 + λ/360).  The root is found with Brent's method in a
±1.8 h window around local mean noon (|equation of time| < 20 min over
the whole range); xtol 1e-9 day.

The civil date is the date in local mean time of the meridian, i.e.
UT1 + λ/360°; for 52.5° E that is UT1 + 3 h 30 m, which coincides with
Iran Standard Time (UTC + 3:30) to within |UT1 - UTC| < 0.9 s.
"""

import math
from dataclasses import dataclass
from functools import lru_cache

from scipy.optimize import brentq

from .config import REFERENCE_LONGITUDE_DEG
from .date_conversion import (
    GregorianDate,
    format_jd,
    gregorian_to_jdn,
    jdn_to_gregorian,
)
from .ephemeris import apparent_solar_hour_angle
from .time_scales import ut1_to_tt


@dataclass(frozen=True)
class SolarNoon:
    date: GregorianDate      # civil date (local mean time at the meridian)
    longitude_deg: float
    jd_ut1: float
    jd_tt: float
    delta_t: float           # seconds

    @property
    def equation_of_time_seconds(self):
        """Apparent minus mean solar time (s); positive = Sun early."""

        mean_noon = self.date.jdn - self.longitude_deg / 360.0
        return (mean_noon - self.jd_ut1) * 86400.0

    def describe(self):
        offset = self.longitude_deg / 360.0
        return (
            f"Apparent noon {self.date} at {self.longitude_deg}° E: "
            f"{format_jd(self.jd_ut1, scale='UT1')} = "
            f"{format_jd(self.jd_ut1, offset, 'LMT')} "
            f"(equation of time {self.equation_of_time_seconds / 60:+.2f} min)"
        )

    def __str__(self):
        return self.describe()


@lru_cache(maxsize=4096)
def _noon_ut1(jdn, longitude_deg, delta_t_offset):
    lon = math.radians(longitude_deg)

    def hour_angle(jd_ut1):
        jd_tt, _ = ut1_to_tt(jd_ut1, delta_t_offset)
        return apparent_solar_hour_angle(jd_ut1, jd_tt, lon)

    mean_noon = jdn - longitude_deg / 360.0
    return brentq(
        hour_angle, mean_noon - 0.075, mean_noon + 0.075,
        xtol=1e-9, rtol=1e-15,
    )


def apparent_solar_noon_jdn(jdn, longitude_deg=REFERENCE_LONGITUDE_DEG,
                            delta_t_offset=0.0):
    """Apparent noon of the civil day ``jdn`` (local mean time date)."""

    jd_ut1 = _noon_ut1(int(jdn), float(longitude_deg), float(delta_t_offset))
    jd_tt, dt = ut1_to_tt(jd_ut1, delta_t_offset)
    return SolarNoon(
        date=jdn_to_gregorian(jdn),
        longitude_deg=longitude_deg,
        jd_ut1=jd_ut1,
        jd_tt=jd_tt,
        delta_t=dt,
    )


def apparent_solar_noon(year, month, day,
                        longitude_deg=REFERENCE_LONGITUDE_DEG,
                        delta_t_offset=0.0):
    """
    Apparent solar noon on proleptic Gregorian date (astronomical year).
    """

    return apparent_solar_noon_jdn(
        gregorian_to_jdn(year, month, day), longitude_deg, delta_t_offset
    )
