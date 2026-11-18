"""
Essential checks of the generator.  Needs the DE441 kernel (see README).
Run from the project folder:  python -m unittest discover -s tests -v
"""

import csv
import math
import os
import re
import sys
import tempfile
import unittest

sys.path.insert(
    0,
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
)

from nowruz_table import (  # noqa: E402
    apparent_solar_longitude,
    apparent_solar_noon,
    compute_nowruz,
    march_equinox,
    write_table,
)
from nowruz_table.date_conversion import (  # noqa: E402
    GregorianDate,
    format_jd,
    gregorian_to_jdn,
    jdn_to_gregorian,
    parse_year,
)
from nowruz_table.delta_t import delta_t_seconds, delta_t_sigma_seconds  # noqa: E402
from nowruz_table.ephemeris import get_ephemeris  # noqa: E402
from nowruz_table.nowruz import _margin  # noqa: E402
from nowruz_table.config import REFERENCE_MEAN_TIME_OFFSET_DAYS  # noqa: E402
from nowruz_table.date_conversion import jd_to_jdn  # noqa: E402

class DateConversion(unittest.TestCase):
    def test_known_jdn(self):
        self.assertEqual(gregorian_to_jdn(2000, 1, 1), 2451545)
        self.assertEqual(gregorian_to_jdn(-4713, 11, 24), 0)
        self.assertEqual(gregorian_to_jdn(1582, 10, 15), 2299161)
        self.assertEqual(jdn_to_gregorian(2451545), GregorianDate(2000, 1, 1))

    def test_round_trip_full_range(self):
        start = gregorian_to_jdn(-2001, 1, 1)
        end = gregorian_to_jdn(3002, 12, 31)
        for jdn in range(start, end + 1, 97):
            self.assertEqual(jdn_to_gregorian(jdn).jdn, jdn)

    def test_weekday_and_years(self):
        self.assertEqual(GregorianDate(2000, 1, 1).weekday, "Saturday")
        self.assertEqual(parse_year("2000 BC"), -1999)
        self.assertEqual(parse_year("1BC"), 0)
        self.assertEqual(parse_year("-1999"), -1999)
        with self.assertRaises(ValueError):
            gregorian_to_jdn(1900, 2, 29)   # proleptic Gregorian rule
        gregorian_to_jdn(-4, 2, 29)         # 5 BC is a leap year


class DeltaT(unittest.TestCase):
    def test_values(self):
        self.assertAlmostEqual(delta_t_seconds(2000.0), 63.86, places=2)
        self.assertAlmostEqual(delta_t_seconds(1900.0), -2.79, places=2)
        self.assertAlmostEqual(delta_t_seconds(0.0), 10583.6, places=1)
        self.assertAlmostEqual(delta_t_sigma_seconds(-500.0), 430.6, delta=1)

    def test_sigma_matches_espenak_tables(self):
        # https://eclipsewise.com/help/uncertainty.html, Tables 1-3
        documented = {-2000: 3732, -1500: 1900, -1000: 622, 0: 265,
                      1000: 54, 1700: 5, 1800: 1, 2500: 612, 3000: 1885}
        for year, sigma in documented.items():
            self.assertAlmostEqual(delta_t_sigma_seconds(year), sigma,
                                   delta=max(1.0, 0.002 * sigma), msg=year)
        with self.assertRaises(ValueError):
            delta_t_seconds(-2500.0)


class Ephemeris(unittest.TestCase):
    def test_segment_selection(self):
        eph = get_ephemeris()
        boundary = 2440432.5
        early = eph.segment(0, 3, boundary - 1)
        late = eph.segment(0, 3, boundary + 1)
        self.assertNotEqual(early.start_jd, late.start_jd)
        # Positions are continuous across the boundary.
        a = early.compute(boundary)
        b = late.compute(boundary)
        self.assertLess(float(sum((a - b) ** 2)) ** 0.5, 1e-3)  # km

    def test_earth_sun_distance(self):
        # Near perihelion (early January) the distance is ~0.9833 au.
        d = apparent_solar_longitude(2460680.5)["distance_au"]
        self.assertAlmostEqual(d, 0.9833, delta=0.0005)


class Equinox(unittest.TestCase):
    # Published March equinoxes (UT, minute precision).
    PUBLISHED_UT = {
        2000: "2000-03-20 07:35",
        2024: "2024-03-20 03:06",
        2025: "2025-03-20 09:01",
        2026: "2026-03-20 14:46",
    }

    def test_published_instants(self):
        for year, text in self.PUBLISHED_UT.items():
            eq = march_equinox(year)
            # Model ΔT is up to ~6 s larger than observed in the 2020s,
            # so compare to within a minute.
            got = format_jd(eq.jd_ut1)[:16]
            minute_after = format_jd(eq.jd_ut1 + 30 / 86400)[:16]
            minute_before = format_jd(eq.jd_ut1 - 30 / 86400)[:16]
            self.assertIn(text, (got, minute_before, minute_after), year)

    def test_longitude_zero_at_root(self):
        for year in (-1999, 1, 2026, 3001):
            eq = march_equinox(year)
            lon = apparent_solar_longitude(eq.jd_tdb)["longitude"]
            lon = lon - 2 * math.pi if lon > math.pi else lon
            self.assertLess(abs(lon), 1e-10)

    def test_time_scale_chain(self):
        eq = march_equinox(2026)
        self.assertLess(abs(eq.jd_tdb - eq.jd_tt) * 86400, 0.002)
        self.assertAlmostEqual(
            (eq.jd_tt - eq.jd_ut1) * 86400, eq.delta_t, places=3
        )


class SolarNoon(unittest.TestCase):
    def test_equation_of_time(self):
        # Equation of time extremes: ~-14.2 min (11 Feb), ~+16.4 min (3 Nov).
        self.assertAlmostEqual(
            apparent_solar_noon(2025, 2, 11).equation_of_time_seconds / 60,
            -14.2, delta=0.1)
        self.assertAlmostEqual(
            apparent_solar_noon(2025, 11, 3).equation_of_time_seconds / 60,
            16.4, delta=0.1)



class Nowruz(unittest.TestCase):
    OFFICIAL = {2016: 20, 2020: 20, 2024: 20, 2025: 21, 2026: 21, 2029: 20}

    def test_official_dates(self):
        for g, day in self.OFFICIAL.items():
            rec = compute_nowruz(g)
            self.assertEqual(jdn_to_gregorian(rec.nowruz_jdn),
                             GregorianDate(g, 3, day))

    def test_critical_offset_flips_decision(self):
        for g in (2223, 2025, -1999):
            rec = compute_nowruz(g)
            civil = jd_to_jdn(rec.equinox.jd_ut1, REFERENCE_MEAN_TIME_OFFSET_DAYS)
            c = rec.critical_delta_t_offset
            self.assertLess(abs(_margin(g, civil, c)), 1e-3)
            before = _margin(g, civil, c - 1.0)
            after = _margin(g, civil, c + 1.0)
            self.assertLess(before, 0.0)
            self.assertGreater(after, 0.0)

    def test_describe(self):
        text = compute_nowruz(2223).describe()
        self.assertIn("2223-03-21", text)
        self.assertIn("ΔT-SENSITIVE", text)

    def test_write_table(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.csv")
            self.assertEqual(write_table(path, 2024, 2026), 3)
            with open(path, encoding="utf-8") as fh:
                lines = fh.read().splitlines()
            meta = dict(
                re.match(r"#\s*([a-z_]+):\s*(.*)$", l).groups()
                for l in lines if re.match(r"#\s*[a-z_]+:", l)
            )
            self.assertEqual(meta["format_version"], "1")
            self.assertIn("Behrouz Safari", meta["created_by"])
            rows = list(csv.DictReader(l for l in lines if not l.startswith("#")))
            self.assertEqual([r["nowruz_date"] for r in rows],
                             ["2024-03-20", "2025-03-21", "2026-03-21"])


if __name__ == "__main__":
    unittest.main()
