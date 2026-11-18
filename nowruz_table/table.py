"""
Writing the Nowruz table (CSV).

File layout
-----------
* A header of lines starting with ``#``, grouped into titled sections
  (``# [Provenance]``, ``# [Coverage]``, ...).  Metadata entries have the
  form ``# key: value`` (keys are lower-case identifiers); section titles,
  rules and the ``[Columns]`` guide are plain comments.
* Then a standard CSV header and one row per Gregorian year.

Columns (one row per Gregorian astronomical year G):

=========================  ===============================================
gregorian_year             G (astronomical: 0 = 1 BC, -1999 = 2000 BC)
persian_year               G - 621 (astronomical: 0 = 1 BH)
nowruz_jdn                 Julian Day Number of Farvardin 1
nowruz_date                same, proleptic Gregorian ISO (astronomical year)
equinox_jd_tdb             March equinox, JD in TDB
equinox_jd_tt              March equinox, JD in TT
equinox_jd_ut1             March equinox, JD in UT1 (model ΔT)
equinox_ut1                same, human readable
delta_t_s                  ΔT = TT - UT1 used (Espenak–Meeus), seconds
delta_t_sigma_s            estimated 1σ uncertainty of ΔT, seconds
noon_jd_ut1                apparent noon at 52.5° E on the equinox's civil
                           date, JD in UT1
equation_of_time_s         apparent - mean solar time at that noon, s
margin_s                   noon - equinox (UT1), s; > 0 -> Nowruz = that date
critical_delta_t_offset_s  ΔT change (s) that would flip the decision;
                           > 0 -> Nowruz one day earlier, < 0 -> later
robust                     |critical offset| > k·σ with k = robust_sigma_k
=========================  ===============================================

JDs are written with 9 decimals (~0.1 ms).
"""

import csv
import datetime
import os

from . import __version__
from .config import (
    DEFAULT_SIGMA_K,
    MAX_GREGORIAN_YEAR,
    MIN_GREGORIAN_YEAR,
    AUTHOR,
    GENERATOR_URL,
    REFERENCE_LONGITUDE_DEG,
    get_kernel_path,
)
from .date_conversion import format_jd, jdn_to_gregorian
from .ephemeris import get_ephemeris
from .nowruz import compute_nowruz

FORMAT_VERSION = 1

COLUMNS = (
    "gregorian_year", "persian_year", "nowruz_jdn", "nowruz_date",
    "equinox_jd_tdb", "equinox_jd_tt", "equinox_jd_ut1", "equinox_ut1",
    "delta_t_s", "delta_t_sigma_s", "noon_jd_ut1", "equation_of_time_s",
    "margin_s", "critical_delta_t_offset_s", "robust",
)


def _row(rec, sigma_k):
    eq = rec.equinox
    return {
        "gregorian_year": rec.gregorian_year,
        "persian_year": rec.persian_year,
        "nowruz_jdn": rec.nowruz_jdn,
        "nowruz_date": jdn_to_gregorian(rec.nowruz_jdn).iso(),
        "equinox_jd_tdb": f"{eq.jd_tdb:.9f}",
        "equinox_jd_tt": f"{eq.jd_tt:.9f}",
        "equinox_jd_ut1": f"{eq.jd_ut1:.9f}",
        "equinox_ut1": format_jd(eq.jd_ut1),
        "delta_t_s": f"{eq.delta_t:.3f}",
        "delta_t_sigma_s": f"{eq.delta_t_sigma:.3f}",
        "noon_jd_ut1": f"{rec.noon.jd_ut1:.9f}",
        "equation_of_time_s": f"{rec.noon.equation_of_time_seconds:.3f}",
        "margin_s": f"{rec.margin_seconds:.4f}",
        "critical_delta_t_offset_s": f"{rec.critical_delta_t_offset:.4f}",
        "robust": "1" if rec.robust(sigma_k) else "0",
    }


_RULE_WIDTH = 78

# Short description of every column, written into the header.
_COLUMN_GUIDE = (
    ("gregorian_year", "Gregorian year (astronomical numbering)"),
    ("persian_year", "Persian year = gregorian_year - 621 (0 = 1 BH)"),
    ("nowruz_jdn", "Julian Day Number of Nowruz (1 Farvardin)"),
    ("nowruz_date", "Nowruz as a proleptic Gregorian ISO date"),
    ("equinox_jd_tdb", "March equinox, Julian Date, TDB"),
    ("equinox_jd_tt", "March equinox, Julian Date, TT"),
    ("equinox_jd_ut1", "March equinox, Julian Date, UT1 (model Delta T)"),
    ("equinox_ut1", "March equinox, UT1, readable"),
    ("delta_t_s", "Delta T = TT - UT1 used, seconds"),
    ("delta_t_sigma_s", "estimated 1-sigma uncertainty of Delta T, seconds"),
    ("noon_jd_ut1", "apparent noon at 52.5 E on the equinox's civil date, UT1"),
    ("equation_of_time_s", "apparent minus mean solar time at that noon, s"),
    ("margin_s", "noon - equinox in seconds (> 0: Nowruz is that date)"),
    ("critical_delta_t_offset_s",
     "Delta T change that flips the date (> 0: earlier, < 0: later)"),
    ("robust", "1 if |critical offset| > robust_sigma_k * sigma, else 0"),
)


def _metadata_sections(first, last, sigma_k):
    """Header content as (section title, [(key, value), ...])."""

    kernel = get_kernel_path()
    start, end = get_ephemeris().coverage
    now = datetime.datetime.now(datetime.timezone.utc)
    return [
        ("Provenance", [
            ("created_by", f"{AUTHOR} ({GENERATOR_URL})"),
            ("generator", f"nowruz_table {__version__}"),
            ("generated_utc", now.strftime("%Y-%m-%dT%H:%M:%SZ")),
            ("format_version", FORMAT_VERSION),
        ]),
        ("Coverage", [
            ("first_gregorian_year", first),
            ("last_gregorian_year", last),
            ("year_numbering", "astronomical (0 = 1 BC, -1999 = 2000 BC)"),
            ("calendar", "proleptic Gregorian for all dates"),
        ]),
        ("Calendar rule", [
            ("reference_longitude_deg", REFERENCE_LONGITUDE_DEG),
            ("civil_date", "local mean time at 52.5 E = UT1 + 3h30m"),
            ("rule", "equinox before apparent noon -> Nowruz on its civil "
                     "date, else the next day"),
        ]),
        ("Ephemeris", [
            ("ephemeris", "JPL DE441"),
            ("kernel", os.path.basename(kernel)),
            ("kernel_bytes", os.path.getsize(kernel)),
            ("kernel_coverage_jd_tdb", f"{start} .. {end}"),
        ]),
        ("Astronomical models", [
            ("equinox", "apparent solar longitude = 0 (mean ecliptic, "
                        "true equinox of date)"),
            ("apparent_place", "light time + annual aberration (ERFA ab)"),
            ("precession", "Vondrak et al. 2011 long-term model (ERFA ltpb)"),
            ("nutation", "IAU 2000A with IAU 2006 adjustments (ERFA nut06a)"),
            ("earth_rotation", "Earth Rotation Angle, CIO-based hour angle "
                               "(ERFA era00, s06)"),
            ("tdb_minus_tt", "ERFA dtdb at the geocentre"),
        ]),
        ("Delta T", [
            ("delta_t_model", "Espenak & Meeus 2006 polynomials"),
            ("delta_t_sigma_model", "Huber 2000 before 1000 BC and after "
                                    "2005; 0.8 t^2 to 1200; table to 2005"),
            ("robust_sigma_k", sigma_k),
        ]),
    ]


def _header_lines(first, last, sigma_k):
    """The '#' comment header of the table file."""

    rule = "# " + "=" * (_RULE_WIDTH - 2)
    sections = _metadata_sections(first, last, sigma_k)
    width = max(len(k) for _, items in sections for k, _ in items) + 2

    lines = [
        rule,
        "# Astronomical Nowruz table - Persian (Solar Hijri) calendar",
        f"# Gregorian years {first} .. {last} (astronomical numbering)",
        rule,
    ]
    for title, items in sections:
        lines += ["#", f"# [{title}]"]
        lines += [f"# {key + ':':<{width}}{value}" for key, value in items]

    col_width = max(len(c) for c, _ in _COLUMN_GUIDE) + 2
    lines += ["#", "# [Columns]"]
    lines += [f"#   {c:<{col_width}}{d}" for c, d in _COLUMN_GUIDE]
    lines.append(rule)
    return lines


def write_table(path, first=MIN_GREGORIAN_YEAR, last=MAX_GREGORIAN_YEAR,
                sigma_k=DEFAULT_SIGMA_K, progress=None):
    """
    Compute Nowruz for Gregorian years ``first`` .. ``last`` and write
    the CSV table to ``path``.  ``progress(done, total)`` is called
    periodically if given.  Returns the number of rows written.
    """

    if first > last:
        raise ValueError("Empty year range.")

    header = _header_lines(first, last, sigma_k)
    total = last - first + 1
    rows = []
    previous_jdn = None

    for i, year in enumerate(range(first, last + 1)):
        rec = compute_nowruz(year)
        if previous_jdn is not None and (
            rec.nowruz_jdn - previous_jdn not in (365, 366)
        ):
            raise RuntimeError(
                f"Year length {rec.nowruz_jdn - previous_jdn} before {year}."
            )
        previous_jdn = rec.nowruz_jdn
        rows.append(_row(rec, sigma_k))
        if progress and (i + 1) % 100 == 0:
            progress(i + 1, total)

    tmp = f"{path}.tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        for line in header:
            fh.write(line + "\n")
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)

    if progress:
        progress(total, total)
    return len(rows)


def print_progress(done, total):
    """A ``progress`` callback for ``write_table`` that prints to stdout."""

    end = "\n" if done == total else ""
    print(f"\r  {done}/{total} years", end=end, flush=True)
