"""
Configuration of the table generator.

Constants of the calendar definition and the location of the DE441 kernel.

Kernel path resolution order:

1. ``set_kernel_path(path)``;
2. the environment variable ``NOWRUZ_TABLE_KERNEL``;
3. ``kernel/de441_sun_earth_2000BC_3000AD.bsp`` inside the project folder
   (the folder that contains ``nowruz_table/``; git-ignored).
"""

import math
import os

KERNEL_FILENAME = "de441_sun_earth_2000BC_3000AD.bsp"
# This file is <project>/nowruz_table/config.py
PROJECT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
)
DEFAULT_KERNEL_PATH = os.path.join(PROJECT_DIR, "kernel", KERNEL_FILENAME)

AUTHOR = "Behrouz Safari"
GENERATOR_URL = "https://github.com/behrouzz/nowruz-table"
CALENDAR_URL = "https://github.com/behrouzz/persian-calendar"
KERNEL_ENV_VAR = "NOWRUZ_TABLE_KERNEL"

# Reference meridian of the calendar: 52°30'00" E (Iran Standard Time).
REFERENCE_LONGITUDE_DEG = 52.5
REFERENCE_LONGITUDE_RAD = math.radians(REFERENCE_LONGITUDE_DEG)

# Offset of local mean time at the reference meridian from UT1, in days
# (52.5° / 360° = 3 h 30 min).  Defines the civil date boundary.
REFERENCE_MEAN_TIME_OFFSET_DAYS = REFERENCE_LONGITUDE_DEG / 360.0

# Supported Gregorian (astronomical) years for March equinox / Nowruz.
# 3001 is included because the length of the year starting in 3000
# needs the Nowruz of 3001; the kernel covers it.
MIN_GREGORIAN_YEAR = -1999
MAX_GREGORIAN_YEAR = 3001

# Number of standard deviations of ΔT for the robustness flag in the table.
DEFAULT_SIGMA_K = 2.0

_kernel_path = None


def set_kernel_path(path):
    """Set the DE441 kernel path (closes any kernel already opened)."""

    global _kernel_path
    from . import ephemeris

    ephemeris.close()
    _kernel_path = os.path.abspath(path)


def get_kernel_path():
    """Return the configured DE441 kernel path."""

    if _kernel_path is not None:
        return _kernel_path

    return os.path.abspath(
        os.environ.get(KERNEL_ENV_VAR, DEFAULT_KERNEL_PATH)
    )
