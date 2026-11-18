"""
Time-scale conversions: TDB <-> TT <-> UT1.

* TDB — argument of the DE441 ephemeris.
* TT  — argument of precession, nutation and the ΔT model.
* UT1 — Earth-rotation angle, apparent solar noon, civil date.

TDB - TT is the periodic relativistic term (|TDB - TT| < 1.7 ms) from
ERFA ``dtdb`` (Fairhead & Bretagnon 1990 series) for a geocentric
observer.  The series' accuracy degrades far from the present, but its
amplitude is bounded, so the error remains at the millisecond level —
eight orders of magnitude below the ΔT uncertainty that dominates this
project.  The conversion is done properly anyway, so that ΔT is never
subtracted from TDB directly.

TT - UT1 = ΔT comes from ``delta_t`` (Espenak–Meeus model).  UTC is not
used anywhere: before 1960 it does not exist, and since 1972
|UT1 - UTC| < 0.9 s.
"""

import erfa

from .delta_t import delta_t_at_jd

_MJD0 = 2400000.5


def tdb_minus_tt_seconds(jd):
    """TDB - TT in seconds at a geocentric observer."""

    return erfa.dtdb(_MJD0, jd - _MJD0, 0.0, 0.0, 0.0, 0.0)


def tt_to_tdb(jd_tt):
    return jd_tt + tdb_minus_tt_seconds(jd_tt) / 86400.0


def tdb_to_tt(jd_tdb):
    # dtdb varies by < 1e-11 s over the 1.7 ms difference of argument.
    return jd_tdb - tdb_minus_tt_seconds(jd_tdb) / 86400.0


def tt_to_ut1(jd_tt, delta_t_offset=0.0):
    """
    TT -> UT1.  ``delta_t_offset`` (s) is added to the model ΔT, for
    sensitivity tests.  Returns (jd_ut1, delta_t_used_seconds).
    """

    # ΔT is evaluated at the UT1 instant (one fixed-point step is ample:
    # the model changes by < 1 s per year).
    dt = delta_t_at_jd(jd_tt) + delta_t_offset
    dt = delta_t_at_jd(jd_tt - dt / 86400.0) + delta_t_offset
    return jd_tt - dt / 86400.0, dt


def ut1_to_tt(jd_ut1, delta_t_offset=0.0):
    """UT1 -> TT.  Returns (jd_tt, delta_t_used_seconds)."""

    dt = delta_t_at_jd(jd_ut1) + delta_t_offset
    return jd_ut1 + dt / 86400.0, dt
