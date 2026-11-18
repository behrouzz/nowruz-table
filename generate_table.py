"""
Generate the full Nowruz table (Gregorian years -1999 .. 3001).
Takes about 1.5 minutes.

Run from the project folder with its venv active:

    python generate_table.py

The table is written to ``nowruz_table.csv`` in this folder.  Copy it to
``persian_calendar/data/nowruz_table.csv`` in the persian-calendar project
(https://github.com/behrouzz/persian-calendar) to use it there.
"""

import os

from nowruz_table import close, print_progress, write_table

OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nowruz_table.csv")

try:
    rows = write_table(OUTPUT, progress=print_progress)
finally:
    close()

print(f"Wrote {rows} rows to {OUTPUT}")
