"""
ΔT = TT - UT1: Espenak & Meeus (NASA) polynomial model, with uncertainty.

Model
-----
The polynomial expressions are those published by F. Espenak and J. Meeus
for the *Five Millennium Canon of Solar Eclipses* (NASA/TP-2006-214141),
https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html .  They are
valid for astronomical years -1999 ... +3000.  Before -500 and after
+2150 they reduce to the long-term parabola -20 + 32 u^2 (u in centuries
from 1820), which describes the mean tidal deceleration of Earth's
rotation but none of its decade-to-millennium fluctuations.

The argument is a decimal year computed continuously from the instant:
``y = 2000 + (JD - JD(2000-01-01)) / 365.2425``.  (NASA uses
``year + (month - 0.5)/12``; the difference is < 1 day and immaterial.)
The model's piecewise branches are not exactly continuous at their joins
(e.g. ~14 s at -500); this is a property of the published model.

The lunar secular-acceleration correction ``-0.000012932 (y-1955)^2`` s,
which adapts the Canon's ΔT to a lunar ephemeris with a different
tidal term, is **not** applied: this project only computes the Sun, and
that correction is not a universal ΔT correction for DE441.  Its size
(~200 s at -2000) is well inside the ΔT uncertainty there.

Uncertainty (1σ, seconds) — the model documented by Espenak & Meeus
---------------------------------------------------------------------
(https://eclipsewise.com/help/uncertainty.html):

* y < -999 (before 1000 BC): Huber (2000), calibration year -500;
* -999 <= y <= 1200: 0.8 t^2 with t = (y - 1820)/100 (Morrison &
  Stephenson 2004): 431 s at -500, 265 s at 0, 54 s at 1000, 31 s at 1200;
* 1200 < y <= 2005: linear interpolation between 31 s at 1200 and the
  telescopic-era standard errors 5 s at 1700, 1 s at 1800, 0.1 s at 1900
  (Stephenson & Houlden 1986), held at 0.1 s to 2005.  The documentation
  gives no values between 1200 and 1700; the interpolation there is ours;
* y > 2005 (prediction): Huber (2000), calibration year 2005, with a
  floor of 0.1 s.

Huber's estimate of the total standard error (fluctuations plus drift)
N years from the calibration year is

    σ = 365.25 N sqrt[(N Q / 3)(1 + N / M)] / 1000   seconds,

with Q = 0.058 ms²/yr and M = 2500 yr.  It gives 3732 s (62 min) at
-2000, 1900 s at -1500, 612 s at +2500 and 1885 s (31 min) at +3000.  The
model has small jumps at its joins (e.g. ~15 s at -1000).

The table stores the critical ΔT offset of every Nowruz, so the
calendar package can apply any multiple of σ (``sigma_k``,
``sigma_scale``) afterwards.
"""

from .config import MAX_GREGORIAN_YEAR, MIN_GREGORIAN_YEAR

# JD (UT) of 2000-01-01 00:00 and mean Gregorian year.
_JD_2000 = 2451544.5
_GREGORIAN_YEAR_DAYS = 365.2425

# Hard limits of the model (one year of slack around the supported
# range; both outer branches are the same parabola).
_Y_MIN = MIN_GREGORIAN_YEAR - 1.0
_Y_MAX = MAX_GREGORIAN_YEAR + 1.0


def decimal_year(jd):
    """Decimal (astronomical) year of a Julian Date."""

    return 2000.0 + (jd - _JD_2000) / _GREGORIAN_YEAR_DAYS


def _check_range(y):
    if not _Y_MIN <= y <= _Y_MAX:
        raise ValueError(
            f"Year {y:.2f} is outside the ΔT model range "
            f"({MIN_GREGORIAN_YEAR} .. {MAX_GREGORIAN_YEAR})."
        )


def delta_t_seconds(y):
    """
    Espenak–Meeus ΔT = TT - UT1 in seconds for decimal year ``y``.
    """

    _check_range(y)

    if y < -500:
        u = (y - 1820) / 100
        return -20 + 32 * u**2

    if y < 500:
        u = y / 100
        return (
            10583.6
            - 1014.41 * u
            + 33.78311 * u**2
            - 5.952053 * u**3
            - 0.1798452 * u**4
            + 0.022174192 * u**5
            + 0.0090316521 * u**6
        )

    if y < 1600:
        u = (y - 1000) / 100
        return (
            1574.2
            - 556.01 * u
            + 71.23472 * u**2
            + 0.319781 * u**3
            - 0.8503463 * u**4
            - 0.005050998 * u**5
            + 0.0083572073 * u**6
        )

    if y < 1700:
        t = y - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129

    if y < 1800:
        t = y - 1700
        return (
            8.83
            + 0.1603 * t
            - 0.0059285 * t**2
            + 0.00013336 * t**3
            - t**4 / 1174000
        )

    if y < 1860:
        t = y - 1800
        return (
            13.72
            - 0.332447 * t
            + 0.0068612 * t**2
            + 0.0041116 * t**3
            - 0.00037436 * t**4
            + 0.0000121272 * t**5
            - 0.0000001699 * t**6
            + 0.000000000875 * t**7
        )

    if y < 1900:
        t = y - 1860
        return (
            7.62
            + 0.5737 * t
            - 0.251754 * t**2
            + 0.01680668 * t**3
            - 0.0004473624 * t**4
            + t**5 / 233174
        )

    if y < 1920:
        t = y - 1900
        return (
            -2.79
            + 1.494119 * t
            - 0.0598939 * t**2
            + 0.0061966 * t**3
            - 0.000197 * t**4
        )

    if y < 1941:
        t = y - 1920
        return 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3

    if y < 1961:
        t = y - 1950
        return 29.07 + 0.407 * t - t**2 / 233 + t**3 / 2547

    if y < 1986:
        t = y - 1975
        return 45.45 + 1.067 * t - t**2 / 260 - t**3 / 718

    if y < 2005:
        t = y - 2000
        return (
            63.86
            + 0.3345 * t
            - 0.060374 * t**2
            + 0.0017275 * t**3
            + 0.000651814 * t**4
            + 0.00002373599 * t**5
        )

    if y < 2050:
        t = y - 2000
        return 62.92 + 0.32217 * t + 0.005589 * t**2

    if y < 2150:
        u = (y - 1820) / 100
        return -20 + 32 * u**2 - 0.5628 * (2150 - y)

    u = (y - 1820) / 100
    return -20 + 32 * u**2


# (year, 1σ seconds) between 1200 and 2005.
_SIGMA_TABLE = (
    (1200.0, 0.8 * ((1200 - 1820) / 100) ** 2),   # Morrison & Stephenson
    (1700.0, 5.0),
    (1800.0, 1.0),
    (1900.0, 0.1),
    (2005.0, 0.1),
)

# Huber (2000) parameters as used by Espenak & Meeus.
_HUBER_Q = 0.058     # ms^2 / yr
_HUBER_M = 2500.0    # yr, span of observed ΔT (-500 .. +2000)


def huber_sigma_seconds(n_years):
    """Huber (2000) total standard error n_years from the calibration year."""

    n = abs(n_years)
    return 365.25 * n * (n * _HUBER_Q / 3 * (1 + n / _HUBER_M)) ** 0.5 / 1000


def delta_t_sigma_seconds(y):
    """Estimated 1σ uncertainty of ΔT at decimal year ``y`` (seconds)."""

    _check_range(y)

    if y < -999:
        return huber_sigma_seconds(-500 - y)

    if y <= 1200:
        t = (y - 1820) / 100
        return 0.8 * t**2

    if y <= 2005:
        for (y0, s0), (y1, s1) in zip(_SIGMA_TABLE, _SIGMA_TABLE[1:]):
            if y0 <= y <= y1:
                return s0 + (s1 - s0) * (y - y0) / (y1 - y0)

    return max(0.1, huber_sigma_seconds(y - 2005))


def delta_t_at_jd(jd):
    """ΔT (s) at an instant given as a JD (UT1 or TT; the difference is
    irrelevant for this slowly varying model)."""

    return delta_t_seconds(decimal_year(jd))


def sigma_at_jd(jd):
    """1σ ΔT uncertainty (s) at an instant given as a JD."""

    return delta_t_sigma_seconds(decimal_year(jd))
