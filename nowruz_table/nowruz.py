"""
Nowruz determination for one Gregorian year.

Rule: let E be the March equinox and D the civil date containing E at
52.5° E (local mean time, UT1 + 3h30m).  Nowruz is D if E is before the
apparent solar noon of D at 52.5° E, otherwise D + 1.

ΔT sensitivity
--------------
In addition to the decision with the model ΔT, the generator stores the
**critical ΔT offset** δc: the change of ΔT (seconds, added to the model)
at which the equinox would coincide with the apparent noon of D, i.e.
at which the decision flips.  It is solved exactly (Newton iteration;
the margin noon - equinox grows with ΔT at a rate of ~0.997, because the
equinox moves 1 s earlier in UT1 per second of ΔT while the noon moves
by only ~0.3 % of that).

* δc > 0: a larger ΔT would move Nowruz one day **earlier**.
* δc < 0: a smaller ΔT would move Nowruz one day **later**.

With σ the 1σ ΔT uncertainty, a Nowruz is robust at k·σ when
|δc| > k·σ.  The calendar package uses this to evaluate robustness for
any k without the ephemeris.  (A change of the civil date D itself never
alters the result: an equinox near midnight is ~12 h from both noons.)
"""

from dataclasses import dataclass

from .config import REFERENCE_LONGITUDE_DEG, REFERENCE_MEAN_TIME_OFFSET_DAYS
from .date_conversion import jd_to_jdn, jdn_to_gregorian
from .equinox import Equinox, march_equinox
from .solar_noon import SolarNoon, apparent_solar_noon_jdn

YEAR_OFFSET = 621  # Gregorian year of Nowruz = Persian year + 621


@dataclass(frozen=True)
class NowruzRecord:
    gregorian_year: int
    nowruz_jdn: int
    equinox: Equinox
    noon: SolarNoon               # apparent noon of the equinox's civil date
    margin_seconds: float         # noon - equinox (UT1); > 0 -> same day
    critical_delta_t_offset: float  # seconds; see module docstring

    @property
    def persian_year(self):
        return self.gregorian_year - YEAR_OFFSET

    @property
    def date(self):
        return jdn_to_gregorian(self.nowruz_jdn)

    def robust(self, sigma_k=2.0):
        return abs(self.critical_delta_t_offset) > (
            sigma_k * self.equinox.delta_t_sigma
        )

    def describe(self, sigma_k=2.0):
        side = "before" if self.margin_seconds > 0 else "at/after"
        status = "robust" if self.robust(sigma_k) else "ΔT-SENSITIVE"
        return "\n".join([
            f"Nowruz of Persian year {self.persian_year}: {self.date}",
            f"  {self.equinox}",
            f"  {self.noon}",
            f"  Equinox {side} noon by {abs(self.margin_seconds):.1f} s; "
            f"decision flips at ΔT offset "
            f"{self.critical_delta_t_offset:+.1f} s -> {status} "
            f"at ±{sigma_k:g}σ",
        ])

    def __str__(self):
        return self.describe()


def _margin(gregorian_year, civil_jdn, delta_t_offset):
    """noon(D) - equinox in UT1 seconds, for model ΔT + offset."""

    eq = march_equinox(gregorian_year, delta_t_offset)
    noon = apparent_solar_noon_jdn(
        civil_jdn, REFERENCE_LONGITUDE_DEG, delta_t_offset
    )
    return (noon.jd_ut1 - eq.jd_ut1) * 86400.0


def critical_delta_t_offset(gregorian_year, civil_jdn, margin_seconds):
    """ΔT offset (s) making the equinox coincide with noon of D."""

    offset = -margin_seconds
    for _ in range(8):
        m = _margin(gregorian_year, civil_jdn, offset)
        offset -= m  # d(margin)/d(offset) ≈ 0.997
        if abs(m) < 1e-4:
            break
    return offset


def compute_nowruz(gregorian_year):
    """Apply the rule with the model ΔT and compute the sensitivity."""

    eq = march_equinox(gregorian_year)
    civil_jdn = jd_to_jdn(eq.jd_ut1, REFERENCE_MEAN_TIME_OFFSET_DAYS)
    noon = apparent_solar_noon_jdn(civil_jdn, REFERENCE_LONGITUDE_DEG)

    margin = (noon.jd_ut1 - eq.jd_ut1) * 86400.0
    jdn = civil_jdn if margin > 0.0 else civil_jdn + 1

    return NowruzRecord(
        gregorian_year=gregorian_year,
        nowruz_jdn=jdn,
        equinox=eq,
        noon=noon,
        margin_seconds=margin,
        critical_delta_t_offset=critical_delta_t_offset(
            gregorian_year, civil_jdn, margin
        ),
    )
