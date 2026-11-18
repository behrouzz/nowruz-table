# nowruz-table — Nowruz table generator (DE441)

A one-off generator. It computes the astronomical Nowruz (Farvardin 1) for
every Gregorian year from **2000 BC to 3001 AD** using JPL **DE441**, and
writes the results to a CSV table. The companion package
[persian-calendar](https://github.com/behrouzz/persian-calendar) ships that
table and needs neither the 219 MB kernel nor any scientific library.

Once the table exists, you need this project again only to regenerate it
(for example, after changing the ΔT model).

## Setup

Python ≥ 3.10. The code runs from the project folder inside a virtual
environment that holds the dependencies listed in `requirements.txt`
(`numpy`, `scipy`, `jplephem`, `pyerfa`).

```bash
git clone https://github.com/behrouzz/nowruz-table.git
cd nowruz-table
python -m venv venv
venv\Scripts\python -m pip install -r requirements.txt        # Windows
# venv/bin/python -m pip install -r requirements.txt          # Linux / macOS
```

Activate the venv from inside the project folder:

```bash
venv\Scripts\activate           # Windows, cmd
.\venv\Scripts\Activate.ps1     # Windows, PowerShell
source venv/bin/activate        # Linux / macOS
```

Use `deactivate` to leave it.

### The DE441 kernel

The reduced kernel `de441_sun_earth_2000BC_3000AD.bsp` (219 MB) is not part
of the repository; `kernel/` is git-ignored. How it was built from the
official DE441 files with NAIF's SPKMERGE, and how to verify it, is described
in [kernel_preparation.md](kernel_preparation.md).

Put it in `kernel/` inside the project folder, or point to it elsewhere. It
is located as follows:

1. `nowruz_table.set_kernel_path(PATH)`, if called;
2. the environment variable `NOWRUZ_TABLE_KERNEL`;
3. by default, `kernel/de441_sun_earth_2000BC_3000AD.bsp` in the project folder.

## Usage

The package is used from Python. Run your scripts (or an interactive
`python` session) from inside the project folder with the venv active, so
that `import nowruz_table` finds the package folder.

### Generating the table

`generate_table.py` computes all 5001 years (about 1.5 min) and writes
`nowruz_table.csv` in the project folder:

```bash
python generate_table.py
```

It does nothing more than this:

```python
from nowruz_table import close, print_progress, write_table

write_table("nowruz_table.csv", progress=print_progress)
close()   # release the DE441 kernel
```

`write_table(path, first=-1999, last=3001, sigma_k=2.0, progress=None)` can
also write a partial table; `sigma_k` sets the k of the informational
`robust` column.

To use a new table in the calendar, copy it to
`persian_calendar/data/nowruz_table.csv` in the
[persian-calendar](https://github.com/behrouzz/persian-calendar) project, or
load it from anywhere with `persian_calendar.load_table(path)`.

### Inspecting the computations

```python
from nowruz_table import march_equinox, apparent_solar_noon, compute_nowruz

print(march_equinox(2026))          # equinox in TT, UT1 and local mean time
print(apparent_solar_noon(2026, 3, 20))   # apparent noon at 52.5° E, equation of time
print(compute_nowruz(2026))         # full Nowruz decision with ΔT sensitivity
print(compute_nowruz(-1999))        # 2000 BC (years are astronomical: -1999 = 2000 BC)

rec = compute_nowruz(2223)          # equinox only 5.7 s before noon
print(rec.date)                     # 2223-03-21
print(rec.margin_seconds)           # noon - equinox, s
print(rec.critical_delta_t_offset)  # ΔT change that flips the date, s
print(rec.robust(sigma_k=2))        # False

eq = march_equinox(2026)
print(eq.jd_tdb, eq.jd_tt, eq.jd_ut1, eq.delta_t, eq.delta_t_sigma)
```

`example.py` contains these examples; run it with `python example.py`.

Tests (these need the kernel): `python -m unittest discover -s tests -v`

## The table

The file starts with a header of `#` comment lines, grouped into sections:

* **Provenance:** author, generator version, generation time, format version.
* **Coverage:** year range and numbering conventions.
* **Calendar rule:** reference meridian, civil date and the Nowruz rule.
* **Ephemeris:** kernel name, size and coverage.
* **Astronomical models:** equinox definition, apparent place, precession,
  nutation, Earth rotation, TDB − TT.
* **Delta T:** ΔT and uncertainty models, and the k of the `robust` column.
* **Columns:** a short description of every column.

Machine-readable entries have the form `# key: value`; section titles and the
column guide are plain comments. After the header comes a normal CSV header
line, then one row per Gregorian year (astronomical numbering: 0 = 1 BC,
−1999 = 2000 BC).

| Column | Meaning |
|---|---|
| `gregorian_year` | G (astronomical) |
| `persian_year` | G − 621 (astronomical: 0 = 1 BH) |
| `nowruz_jdn` | Julian Day Number of Farvardin 1 |
| `nowruz_date` | same, proleptic Gregorian ISO date |
| `equinox_jd_tdb`, `equinox_jd_tt`, `equinox_jd_ut1` | March equinox in TDB, TT, UT1 (9 decimals ≈ 0.1 ms) |
| `equinox_ut1` | human-readable UT1 |
| `delta_t_s`, `delta_t_sigma_s` | ΔT = TT − UT1 used, and its estimated 1σ |
| `noon_jd_ut1` | apparent noon at 52.5° E on the civil date containing the equinox |
| `equation_of_time_s` | apparent − mean solar time at that noon |
| `margin_s` | noon − equinox (UT1). If > 0, Nowruz is that date; otherwise it is the next date |
| `critical_delta_t_offset_s` | the change of ΔT that would flip the decision. If > 0, Nowruz would move one day earlier; if < 0, one day later |
| `robust` | 1 if \|critical offset\| > k·σ (k = `robust_sigma_k`) |

**Why store the critical offset?** In UT1, the equinox moves 1 s earlier for
every second added to ΔT, while the noon moves by only about 0.3 % of that.
The offset at which equinox and noon coincide is solved exactly (Newton
iteration on the full calculation). With that offset stored, the calendar
package can test robustness for *any* k without the ephemeris: a Nowruz is
robust at kσ when |δc| > kσ. A change of the civil date itself, when the
equinox is near midnight, never changes the result, because such an equinox
is about 12 h from both noons.

## Calendar rule

For Gregorian year G:

1. E is the March equinox.
2. D is the civil date containing E at **52.5° E**, in local mean time
   (UT1 + 3 h 30 m, i.e. Iran Standard Time).
3. Nowruz of Persian year G − 621 is D if E is *before* apparent solar noon
   of D at 52.5° E, and D + 1 otherwise.

No arithmetic leap cycle is used.

## Astronomical method

| Step | Model | Frame / time scale | Validity |
|---|---|---|---|
| Barycentric Earth & Sun | DE441 Chebyshev, explicit segment choice (split at JD 2440432.5); Earth = EMB + (Earth−EMB) | ICRF, TDB | −2000 … 3002 |
| Light time | iterated Sun(t−τ) − Earth(t) | BCRS, TDB | — |
| Annual aberration | ERFA `ab` (relativistic); Earth velocity from the Chebyshev derivative | GCRS | — |
| Bias + precession | **Vondrák et al. 2011** (`ltpb`, `ltpecl`, `ltpequ`) | GCRS → mean equator/equinox of date | ±200 000 yr |
| Nutation | IAU 2000A + IAU 2006 adjustments (`nut06a`) | → true equator/equinox | modern fit (see limitations) |
| Apparent longitude | rotation by the true obliquity | **mean ecliptic of date, true equinox of date** | — |
| Equinox | Newton start + Brent root (xtol 1e-9 d) on longitude wrapped to (−180°, 180°] | TDB → TT (`dtdb`) → UT1 (ΔT) | — |
| ΔT | Espenak–Meeus 2006 polynomials | TT − UT1 | −1999 … 3000 |
| Apparent noon | root of H = ERA(UT1) + λ − α_CIO (`era00`, `s06`, `c2ixys`) | UT1 for rotation; TT/TDB for the Sun | — |

**Equinox definition.** The Sun's apparent geocentric longitude, measured from
the *true equinox of date* along the Vondrák *mean ecliptic of date*, passes
through 0°. If the IAU 2006 precession (valid for about ±1000 yr) were used
instead of Vondrák, the equinox would shift by 0 s at 2000, by −8 s at year 0
and by −101 s at 2000 BC.

**Apparent noon** is the instant when the local apparent hour angle of the
Sun's apparent place is zero, so the equation of time is included. It is not
the mean-time expression UT1 + λ/360.

## Validation

* Equinoxes match published instants to the minute (2000-03-20 07:35 UT,
  2024 03:06, 2025 09:01, 2026 14:46). With the observed ΔT, 2025 becomes
  09:01:29 UT1, which matches the official 12:31:30 IRST.
* Equation of time: −14.2 min (11 Feb) and +16.4 min (3 Nov).
* Official Iranian Nowruz dates and leap years (1391, 1395, 1399, 1403,
  1408) match.
* The Borkowski (1996) / jalaali-js algorithm gives identical dates for
  1200–1499 AP. Across 1–2379 AP the two differ in 8 years, all with the
  equinox within 6 s – 17 min of noon.
* Meeus mean-equinox formula: within 20 min over −1000 … 3000.
* The full table has 1211 leap years in 5000 (rate 0.2422). 150 Nowruz
  dates are ΔT-sensitive at ±2σ (70 at ±1σ, 235 at ±3σ).

## Limitations and uncertainties

The root finders resolve instants to about 0.1 ms. That is numerical
precision, not accuracy:

1. **ΔT dominates outside roughly 1700–2025.**
   * The 1σ uncertainty follows the model documented by Espenak & Meeus
     ([eclipsewise.com/help/uncertainty.html](https://eclipsewise.com/help/uncertainty.html)):
     - **1000 BC – 1200 AD:** 0.8 t² s (t in centuries from 1820; Morrison &
       Stephenson 2004). That gives 431 s at 500 BC, 265 s at year 0 and
       54 s at 1000.
     - **1200 – 2005:** linear interpolation from 31 s at 1200 to the
       telescopic-era errors of 5 s at 1700, 1 s at 1800 and 0.1 s at 1900.
       The documentation gives no values between 1200 and 1700, so the
       interpolation there is ours.
     - **Before 1000 BC and after 2005:** Huber (2000),
       σ = 365.25·N·√[(N·Q/3)(1 + N/M)]/1000 s, with N years from the
       calibration year (−500 or 2005), Q = 0.058 ms²/yr and M = 2500 yr.
       That gives 62 min at 2000 BC, 10 min in 2500 and 31 min in 3000.
   * The Espenak–Meeus prediction is already about 5 s too large in 2025,
     consistent with σ ≈ 4.6 s there.
   * The ΔT polynomial branches have small jumps at their joins (e.g. ~14 s
     at 500 BC), and so does the σ model (~15 s at 1000 BC).
   * The lunar secular-acceleration correction is deliberately **not**
     applied. It belongs to lunar ephemeris comparisons, not to the Sun from
     DE441.
2. **Nutation** (IAU 2000A) is a fit to modern data. Because nutation is
   bounded, its error at ±4000 yr is estimated at a few arcseconds, i.e.
   about a minute in the equinox or less.
3. **Hour angle**: the error of `s06` far from J2000 amounts to at most a few
   seconds in the noon instant. Polar motion, diurnal aberration (~0.02 s)
   and parallax (no effect on the noon instant) are neglected.
4. **TDB − TT** (`dtdb`) is bounded below 1.7 ms.
5. The civil date uses local mean time at 52.5° E. For the past this is a
   convention, and all dates are proleptic.

## Project layout

```
nowruz_table/
  config.py           constants, kernel path
  date_conversion.py  JDN <-> proleptic Gregorian, BC/AD, formatting
  delta_t.py          Espenak–Meeus ΔT + uncertainty model
  time_scales.py      TDB <-> TT <-> UT1
  ephemeris.py        DE441 segments, apparent Sun, long-term precession-nutation
  equinox.py          March equinox root finding
  solar_noon.py       apparent solar noon
  nowruz.py           Nowruz decision + critical ΔT offset
  table.py            CSV writer
generate_table.py     writes the full table to nowruz_table.csv
example.py            usage examples
tests/test_generator.py
requirements.txt      numpy, scipy, jplephem, pyerfa
kernel_preparation.md how the reduced DE441 kernel is built
kernel/               the DE441 kernel (git-ignored, not in the repository)
```
