# Preparing the DE441 kernel

`nowruz-table` reads a **reduced** copy of the JPL DE441 planetary ephemeris,
`kernel/de441_sun_earth_2000BC_3000AD.bsp` (219 MB). It is not stored in the
repository. This document explains how it was made, so that you can rebuild
it.

## Why a reduced kernel?

The full DE441 is distributed as two SPK files totalling about 3 GB, and it
covers 13 200 BC to 17 191 AD for all planets, the Moon and Pluto. This project
needs only:

| NAIF ID | Body | Why it is needed |
|---:|---|---|
| 10 | Sun | position of the Sun relative to the Solar System barycentre (0) |
| 3 | Earth-Moon barycentre (EMB) | DE441 gives the Earth only relative to the EMB, so the EMB is the link to the barycentre |
| 399 | Earth | position of the Earth relative to the EMB (3) |

The Earth's barycentric position is then *Earth = EMB (0→3) + Earth relative
to EMB (3→399)*.

The time range is cut to 1 January 2001 BC – 1 January 3002 AD. The
calendar needs the March equinoxes of 2000 BC (astronomical year −1999)
through 3001 AD, so this leaves about a year of margin at each end.

NASA's NAIF utility **SPKMERGE** extracts these bodies and this time range
into a new SPK file. The extraction copies the original DE441 Chebyshev
coefficients unchanged, so the reduced kernel has exactly the precision of
DE441.

## Step 1: Download the NAIF utilities and the leap-seconds kernel

Windows 64-bit:

| File | URL |
|---|---|
| `spkmerge.exe` | https://naif.jpl.nasa.gov/pub/naif/utilities/PC_Windows_64bit/spkmerge.exe |
| `brief.exe` | https://naif.jpl.nasa.gov/pub/naif/utilities/PC_Windows_64bit/brief.exe |
| `latest_leapseconds.tls` | https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls |

For Linux or macOS, take `spkmerge` and `brief` from the folder for your
platform under https://naif.jpl.nasa.gov/pub/naif/utilities/. On those
systems, make them executable with `chmod +x spkmerge brief`.

SPKMERGE always needs a leap-seconds kernel to convert time strings, even
though the times below are given in TDB.

## Step 2: Download the two DE441 files

| File | URL |
|---|---|
| `de441_part-1.bsp` | https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp |
| `de441_part-2.bsp` | https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp |

Each file is roughly 1.6 GB. DE441 is split in time: part 1 ends and part 2
begins on **1969-07-30** (JD 2440432.5).

## Step 3: Create a working folder

The examples use `C:\Astronomy\DE441`; any folder works if you adapt the
paths.

```
C:\Astronomy\DE441\
    tools\
        spkmerge.exe
        brief.exe
    kernels\
        latest_leapseconds.tls
    original\
        de441_part-1.bsp
        de441_part-2.bsp
    output\
    subset_de441.txt          (created in step 4)
```

## Step 4: Write the SPKMERGE command file

Create `C:\Astronomy\DE441\subset_de441.txt` with this content:

```
LEAPSECONDS_KERNEL = C:\Astronomy\DE441\kernels\latest_leapseconds.tls

SPK_KERNEL         = C:\Astronomy\DE441\output\de441_sun_earth_2000BC_3000AD.bsp
LOG_FILE           = C:\Astronomy\DE441\output\spkmerge.log

BODIES             = 3, 10, 399
BEGIN_TIME         = 1 JAN 2001 B.C. 00:00:00 TDB
END_TIME           = 1 JAN 3002 A.D. 00:00:00 TDB

SOURCE_SPK_KERNEL  = C:\Astronomy\DE441\original\de441_part-1.bsp
SOURCE_SPK_KERNEL  = C:\Astronomy\DE441\original\de441_part-2.bsp
```

What each line does:

| Keyword | Meaning |
|---|---|
| `LEAPSECONDS_KERNEL` | the leap-seconds file from step 1 |
| `SPK_KERNEL` | the output file to create (it must not already exist) |
| `LOG_FILE` | a log of what was done; SPKMERGE also stores it in the comment area of the new kernel |
| `BODIES` | the NAIF IDs to extract: EMB (3), Sun (10), Earth (399) |
| `BEGIN_TIME`, `END_TIME` | the time range to keep, in TDB, the time scale of DE441 |
| `SOURCE_SPK_KERNEL` | the input files, both DE441 parts |

`BODIES`, `BEGIN_TIME` and `END_TIME` come before the `SOURCE_SPK_KERNEL`
lines, so they apply to both source files.

## Step 5: Run SPKMERGE

In a Windows command prompt (cmd):

```
cd /d C:\Astronomy\DE441\tools
spkmerge.exe C:\Astronomy\DE441\subset_de441.txt
```

This creates `C:\Astronomy\DE441\output\de441_sun_earth_2000BC_3000AD.bsp`
(about 219 MB) and `spkmerge.log`.

The log, which is also embedded in the kernel, shows the commands as
SPKMERGE applied them: split per source file, with the times converted to
UTC. For the kernel used by this project it reads:

```
SOURCE_SPK_KERNEL = ...\de441_part-1.bsp
    INCLUDE_COMMENTS = NO
    BODIES           = 3, 10, 399
    BEGIN_TIME       = 2002 B.C. DEC 31 23:59:18.814
    END_TIME         = 1969 JUL 29 23:59:18.817
SOURCE_SPK_KERNEL = ...\de441_part-2.bsp
    INCLUDE_COMMENTS = NO
    BODIES           = 3, 10, 399
    BEGIN_TIME       = 1969 JUL 29 23:59:18.817
    END_TIME         = 3001 DEC 31 23:58:50.817
```

The "23:59:18" style times are UTC equivalents of the TDB limits (TDB is
ahead of UTC by about 41 s at the start and 69 s at the end). They are not a different
range.

## Step 6: Verify the result

### With BRIEF

```
brief.exe -c C:\Astronomy\DE441\output\de441_sun_earth_2000BC_3000AD.bsp
```

BRIEF lists the bodies in the file, with `-c` showing the centre each one is
given relative to. You should see exactly the Sun (10) and the EMB (3)
relative to the Solar System barycentre (0), and the Earth (399) relative to
the EMB (3), all covering the requested range.

### With Python (this project's venv)

Copy the kernel to `kernel/` in this project (see the next step), activate
the venv, and run from the project folder:

```python
from jplephem.spk import SPK

kernel = SPK.open("kernel/de441_sun_earth_2000BC_3000AD.bsp")
print(kernel)
for s in kernel.segments:
    print(s.center, s.target, s.start_jd, s.end_jd)
```

Expected output:

```
File type DAF/SPK and format LTL-IEEE with 6 segments:
-2000-01-01..1969-07-30  Type 2  Earth Barycenter (3) -> Earth (399)
-2000-01-01..1969-07-30  Type 2  Solar System Barycenter (0) -> Sun (10)
-2000-01-01..1969-07-30  Type 2  Solar System Barycenter (0) -> Earth Barycenter (3)
1969-07-30..3002-01-01  Type 2  Earth Barycenter (3) -> Earth (399)
1969-07-30..3002-01-01  Type 2  Solar System Barycenter (0) -> Sun (10)
1969-07-30..3002-01-01  Type 2  Solar System Barycenter (0) -> Earth Barycenter (3)
3 399 990574.5 2440432.5
0 10 990574.5 2440432.5
0 3 990574.5 2440432.5
3 399 2440432.5 2817517.5
0 10 2440432.5 2817517.5
0 3 2440432.5 2817517.5
```

Each body has **two segments**, one from each DE441 part, meeting at
JD 2440432.5 (1969-07-30). This is why `nowruz_table.ephemeris` selects the
segment covering each requested date explicitly instead of using jplephem's
`kernel[center, target]`, which returns only one segment.

The kernel covers JD(TDB) 990574.5 – 2817517.5. jplephem prints these
limits as dates in the proleptic Gregorian calendar, with astronomical
years: −2000 = 2001 BC.

## Step 7: Install it in this project

Copy the file into `kernel/` inside the `nowruz-table` project folder:

```
nowruz-table/
    kernel/
        de441_sun_earth_2000BC_3000AD.bsp
```

`kernel/` is git-ignored. To keep the kernel somewhere else, set the
environment variable `NOWRUZ_TABLE_KERNEL` to its path, or call
`nowruz_table.set_kernel_path(path)` (see the [README](README.md)).

The original `de441_part-1.bsp` and `de441_part-2.bsp` are not needed by this
project once the reduced kernel exists.

## Reference copy

The kernel used to generate the bundled Nowruz table has:

| Property | Value |
|---|---|
| Size | 219 239 424 bytes |
| SHA-256 | `dfa8fc6f6569e67dc66f8f76a2fb72dbca37f90d240fb091a713f1c4df9a1dbd` |
| Created by SPKMERGE | 2026-11-17 16:32:35 |

A rebuilt kernel contains the same ephemeris data. Its checksum will
nevertheless differ, because the comment area embeds the creation time and
the file paths. Compare the segment listing above, not the checksum.
