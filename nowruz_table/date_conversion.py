"""
Julian Date and proleptic Gregorian calendar conversion.

Conventions
-----------
* Years are **astronomical** integers: 1 = 1 AD, 0 = 1 BC, -1 = 2 BC,
  -1999 = 2000 BC.  Use ``historical_year`` / ``parse_year`` to convert
  to and from BC/AD notation.
* The calendar is the **proleptic Gregorian** calendar for every date,
  including those before the 1582 reform (no Julian-calendar switch).
* A Julian Day Number (JDN) is the integer JD of the day's noon; the
  civil day with JDN n runs from JD n - 0.5 to JD n + 0.5 (in whatever
  time scale the JD refers to).

Everything uses integer floor division, so the algorithms are valid for
negative years and negative JDNs alike.  Python's ``datetime`` is not
used because it cannot represent years < 1.
"""

import math
import re
from typing import NamedTuple

# JDN of 1970-01-01 (proleptic Gregorian).
_JDN_UNIX_EPOCH = 2440588

GREGORIAN_MONTH_NAMES = (
    "January", "February", "March", "April", "May", "June", "July",
    "August", "September", "October", "November", "December",
)

WEEKDAY_NAMES = (
    "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday",
    "Sunday",
)


def is_gregorian_leap_year(year):
    """Proleptic Gregorian leap-year rule (astronomical year)."""

    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def gregorian_month_length(year, month):
    if month == 2:
        return 29 if is_gregorian_leap_year(year) else 28
    return (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month - 1]


def historical_year(year):
    """Astronomical year -> 'NNNN BC' / 'NNNN AD' string."""

    if year <= 0:
        return f"{1 - year} BC"
    return f"{year} AD"


def parse_year(text):
    """
    Parse a Gregorian year as astronomical number.

    Accepts '2026', '-1999', '2000BC', '2000 BC', '500 AD', '500CE',
    '44 BCE'.  A bare integer is taken as an astronomical year.
    """

    s = str(text).strip().upper().replace(" ", "")
    m = re.fullmatch(r"([+-]?\d+)(BCE|BC|AD|CE)?", s)
    if not m:
        raise ValueError(f"Cannot parse year: {text!r}")

    value = int(m.group(1))
    era = m.group(2)

    if era in ("BC", "BCE"):
        if value < 1:
            raise ValueError("BC years start at 1 BC.")
        return 1 - value
    if era in ("AD", "CE") and value < 1:
        raise ValueError("AD years start at 1 AD.")
    return value


class GregorianDate(NamedTuple):
    """Proleptic Gregorian date with astronomical year numbering."""

    year: int
    month: int
    day: int

    def iso(self):
        """ISO-8601-like string with astronomical year (e.g. -1999-03-21)."""

        sign = "-" if self.year < 0 else ""
        return f"{sign}{abs(self.year):04d}-{self.month:02d}-{self.day:02d}"

    def historical(self):
        """e.g. '21 March 2000 BC'."""

        return (
            f"{self.day} {GREGORIAN_MONTH_NAMES[self.month - 1]} "
            f"{historical_year(self.year)}"
        )

    @property
    def jdn(self):
        return gregorian_to_jdn(self.year, self.month, self.day)

    @property
    def weekday(self):
        return WEEKDAY_NAMES[jdn_weekday(self.jdn)]

    def __str__(self):
        if self.year <= 0:
            return f"{self.iso()} ({self.historical()})"
        return self.iso()


def validate_gregorian(year, month, day):
    if not 1 <= month <= 12:
        raise ValueError(f"Invalid Gregorian month: {month}")
    if not 1 <= day <= gregorian_month_length(year, month):
        raise ValueError(f"Invalid Gregorian day: {year}-{month}-{day}")


def gregorian_to_jdn(year, month, day):
    """Proleptic Gregorian date -> Julian Day Number (integer)."""

    validate_gregorian(year, month, day)

    # Days-from-civil (H. Hinnant), year starting on 1 March.
    y = year - (1 if month <= 2 else 0)
    era = y // 400
    yoe = y - era * 400
    mp = (month + 9) % 12
    doy = (153 * mp + 2) // 5 + day - 1
    doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
    return era * 146097 + doe - 719468 + _JDN_UNIX_EPOCH


def jdn_to_gregorian(jdn):
    """Julian Day Number -> proleptic ``GregorianDate``."""

    z = int(jdn) - _JDN_UNIX_EPOCH + 719468
    era = z // 146097
    doe = z - era * 146097
    yoe = (doe - doe // 1460 + doe // 36524 - doe // 146096) // 365
    doy = doe - (365 * yoe + yoe // 4 - yoe // 100)
    mp = (5 * doy + 2) // 153
    day = doy - (153 * mp + 2) // 5 + 1
    month = mp + 3 if mp < 10 else mp - 9
    year = yoe + era * 400 + (1 if month <= 2 else 0)
    return GregorianDate(year, month, day)


def jdn_weekday(jdn):
    """0 = Monday ... 6 = Sunday (JDN 0 was a Monday)."""

    return int(jdn) % 7


def jd_to_jdn(jd, offset_days=0.0):
    """
    Civil day (JDN) containing instant ``jd``.

    ``offset_days`` shifts from the JD's time scale to local time, e.g.
    52.5/360 for local mean time at 52.5° E when ``jd`` is UT1.
    """

    return math.floor(jd + 0.5 + offset_days)


def jd_to_gregorian_datetime(jd, offset_days=0.0, round_seconds=None):
    """
    Instant -> (GregorianDate, hour, minute, second as float).

    The time scale is that of ``jd`` shifted by ``offset_days``.
    ``round_seconds`` (number of decimals) rounds before splitting, so
    that the carry into minutes/hours/days is handled correctly.
    """

    x = jd + 0.5 + offset_days
    jdn = math.floor(x)
    seconds = (x - jdn) * 86400.0
    if round_seconds is not None:
        seconds = round(seconds, round_seconds)

    if seconds >= 86400.0:
        jdn += 1
        seconds -= 86400.0

    hour = int(seconds // 3600)
    minute = int((seconds - hour * 3600) // 60)
    second = seconds - hour * 3600 - minute * 60
    return jdn_to_gregorian(jdn), hour, minute, second


def format_jd(jd, offset_days=0.0, scale=""):
    """Format an instant as 'YYYY-MM-DD hh:mm:ss.s SCALE'."""

    date, h, m, s = jd_to_gregorian_datetime(jd, offset_days, round_seconds=1)
    return f"{date.iso()} {h:02d}:{m:02d}:{s:04.1f} {scale}".rstrip()
