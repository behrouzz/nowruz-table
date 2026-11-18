"""
Examples of the astronomical computations behind the Nowruz table.

Run from the project folder with its venv active:

    python example.py
"""

from nowruz_table import (
    apparent_solar_noon,
    close,
    compute_nowruz,
    march_equinox,
    write_table,
)

# March equinox of a Gregorian (astronomical) year: TT, UT1 and local time.
print(march_equinox(2026))

# Apparent solar noon at 52.5° E on a date (includes the equation of time).
print(apparent_solar_noon(2026, 3, 20))

# Full Nowruz decision for a year, with its ΔT sensitivity.
print(compute_nowruz(2026))
print(compute_nowruz(2223))     # equinox only 5.7 s before noon
print(compute_nowruz(-1999))    # 2000 BC

# Individual values.
rec = compute_nowruz(2025)
print(rec.date, rec.margin_seconds, rec.critical_delta_t_offset, rec.robust())

# A small table for a few years (the full one: see generate_table.py).
write_table("nowruz_2020_2030.csv", first=2020, last=2030)

close()   # release the DE441 kernel
