"""
nowruz_table — one-off generator of the astronomical Nowruz table.

Computes, from JPL DE441, the March equinox, the apparent solar noon at
52.5° E and the resulting Nowruz date for every Gregorian year
-1999 .. 3001, with ΔT and its sensitivity, and writes them to a CSV
file used by the ``persian_calendar`` package
(https://github.com/behrouzz/persian-calendar).

    from nowruz_table import write_table, print_progress
    write_table("nowruz_table.csv", progress=print_progress)
"""

__version__ = "1.0.0"

from .config import get_kernel_path, set_kernel_path  # noqa: E402
from .ephemeris import apparent_solar_longitude, close  # noqa: E402
from .equinox import Equinox, march_equinox  # noqa: E402
from .nowruz import NowruzRecord, compute_nowruz  # noqa: E402
from .solar_noon import SolarNoon, apparent_solar_noon  # noqa: E402
from .table import (  # noqa: E402
    COLUMNS,
    FORMAT_VERSION,
    print_progress,
    write_table,
)
