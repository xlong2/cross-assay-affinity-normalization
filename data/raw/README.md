# Raw inputs

Neither file here is produced by any script in this repository. They are
supplied inputs, which is why they live in `raw/`.

## `skempi_v2.csv`

The upstream [SKEMPI v2](https://life.bsc.es/pid/skempi2) release, unmodified.
No code in this repository reads it; it is kept as the provenance anchor for
the file below.

## `skempi_v2_temperature_parsed.csv`

The table the pipeline actually reads (`src/calibrate.py`). It is
`skempi_v2.csv` plus four derived columns:

| Column | Meaning |
| --- | --- |
| `Temp parsed ` | temperature parsed out of SKEMPI's free-text `Temperature` field |
| `Delta_G_wt` | wild-type binding free energy |
| `Delta_G_mut ` | mutant binding free energy |
| `Delta_Delta_G` | `Delta_G_mut ` − `Delta_G_wt` |

**Provenance gap:** the script that produced these four columns is not part of
this repository and does not appear anywhere in its git history. The derivation
was done outside the project, so this file cannot be regenerated from
`skempi_v2.csv` with the code provided here. Treat it as a primary input.

Note the trailing spaces in `Delta_G_mut ` and `Temp parsed ` — they are part of
the real column names and `src/calibrate.py` depends on them.
