# User manual

What each script does, what it reads and writes, and where its settings live.
For installation and the one-command run, see [README.md](README.md).

## Running a single script

`make` is the normal entry point. To run one script by hand, stay in the repo root (paths are relative to it):

```bash
export PYTHONPATH=$PWD
python scripts/fig1.py
```

Figure scripts read generated data, so run `make` once first. `make fig1` (or `fig2`, `fig3`, `fig4`, `si`)
builds one figure plus whatever it needs.

## How the pipeline fits together

```
data/raw/*.rdata --> convert_data.py --> data/raw/*.csv --+--> run_pipeline.py ---> fig1, si_filtered, fig2cd
                                                          +--> compute_pairs.py --> figures_pairs
                                                          +--> soil.py, supp/*
neutral_model/v2d --> data/senm/ (500 x 500, 100 runs) -----> run_pipeline.py, supp/*
                  --> neutral_fr.py (1500 x 1500, 3 runs) --> figures_pairs
```

`fig2ab.py` stands alone: it simulates its own small lattice.

## Data preparation

### `scripts/convert_data.py`

Converts the Dryad census files into the tables everything else reads.

- **Reads:** `data/raw/bci.tree<c>.rdata`
- **Writes:** `data/raw/bci.tree<c>.csv`, columns `x, y, name, dbh`
- **Does:** keeps alive stems (`status == 'A'`) and drops stems without coordinates.

This is the slowest small step (about 4 minutes on a laptop).

### `neutral_model/` (`v2d`)

C simulator of the voter model with speciation, sampled at stationarity through its coalescent dual, so it needs
no burn-in. `build.sh` compiles it with `-DDUMP`, which makes it write configurations.

```bash
./v2d <NX> <NY> <NU> <NUMREL> <kernel>
```

- `kernel` is the half-width of a square neighbourhood, in sites.
- It writes `AConf_*` (lines `kx ky species`) and `Abund_*` into the **current directory**.
- It exits with status 1 even on success, so check for the output files instead.
- `V2D_SEED` fixes the seed (the Makefile sets 2026). Without it the seed comes from the clock.

The Makefile runs it inside `data/senm/` as `v2d 500 500 2e-5 100 14`, with the values taken from `config.yaml`.

### `scripts/run_pipeline.py`

Computes everything the spectral figures need, for each scale ℓ = 5–49 m.

- **Reads:** census CSVs, `data/senm/AConf_*`
- **Writes:** in `data/processed/scale_<ℓ>/`:
  - `forest_spectra.npy`: shape (8, 100), one spectrum per census
  - `senm_spectra.npy`: shape (100, 100), one spectrum per realization
  - `forest_correlation_mean.npy`, `forest_correlation_std.npy`
- **Does:** bins each species on a grid of side ℓ, z-scores each row, and takes the Pearson matrix and its eigenvalues.
  BCI species are the top 100 of `reference_census`, tracked across censuses. SENM species are re-selected in
  each realization, because they have no identity across runs.

Spectra are stored in descending order. Rows of the matrices are in **alphabetical** species order, not abundance order.

### `scripts/compute_pairs.py`

Conspecific pair probability F(r) for BCI.

```bash
python scripts/compute_pairs.py {all|species|size} [--raw data/raw] [--out data/processed/pairs]
```

| Mode | Writes | Content |
|---|---|---|
| `all` | `fofr_bci.pkl` | F(r) for the eight censuses |
| `species` | `perspecies_c4.pkl`, `perspecies_c8.pkl` | pair counts of each species with ≥ 200 stems, censuses 4 and 8 |
| `size` | `fofr_size.pkl` | F(r) by dbh class (10–49, 50–99, ≥ 100 mm), censuses 4 and 8 |

Only stems with dbh > 0 inside the plot are used. The method is described under `src/pairs.py`.

### `scripts/neutral_fr.py`

F(r) of the neutral model, for Fig. 3b.

```bash
python scripts/neutral_fr.py --side 1500 --spacing 1 --K 14 --nu 2e-5 [--nrel 3] [--v2d PATH] [--out DIR]
```

- `--side`: lattice side in sites. `--spacing`: metres per site. Keep the spacing that defines the metre scale of
  the spectra (1 m per site), otherwise K means a different distance in Fig. 2 and Fig. 3.
- **Does:** runs `v2d` in a temporary folder, cuts the central 1000 × 500 m window, and computes F(r) per realization.
  It stops with an error if the window doesn't fit in the lattice or the `AConf_*` files are missing.
- **Writes:** `data/processed/pairs/cal_K<K>_nu<nu>.pkl`, with keys `rc`, `F`, `S` (species in the window), `K`, `nu`, `spacing_m`.
- **Prints:** the species count, the large-scale slope b and F(25 m).

## Main-text figures

Outputs go to `figures/`.

### `scripts/fig1.py`: Fig. 1

- **Panels:** (a) map of three example species on the 100 m grid, (b) their occupation vectors on that grid,
  (c) the correlation matrix at ℓ = 20 m, ordered by community.
- **Reads:** `scale_20` correlation matrix (through `src/figdata.py`), census CSV for the map.
- **Writes:** `fig1.{pdf,svg,png}`
- **Settings at the top of the script:** `SP` (species and colours), `MAP_C` (census shown, 8 = 2015), `L_MAP`, `L_C`.

### `scripts/si_filtered.py`: SI filtered matrix

- **Panels:** (a) raw C, (b) MP-filtered C, in the same order as Fig. 1c, (c) spectrum with the MP edge.
- **Writes:** `si_filtered.{pdf,png}`
- **Prints:** the caption numbers: modes kept, variance fraction, raw-vs-filtered correlation, and mean ρ within and between communities.

### `scripts/fig2ab.py`: Fig. 2a–b cartoon

- **Panels:** (a) the update rule on a neighbourhood with a square kernel, (b) one patch coarse-grained at five block sizes.
- **Does:** simulates its own 120 × 60 lattice with a fixed seed. These numbers are for illustration and are separate from
  the calibrated SENM. Settings are in the `CONFIG` block at the top.
- **Writes:** `fig2_ab.{pdf,svg,png}`

### `scripts/fig2cd.py`: Fig. 2c–d

- **Panels:** (c) spectra at ℓ = 5 and 28 m, each divided by the MP edge of its own domain,
  (d) modes above noise against ℓ (5–45 m), BCI vs the neutral model.
- **Reads:** saved spectra only. Counts are recomputed from the spectra and the MP edge.
- **Writes:** `fig2_cd.{pdf,png}` and `fig2d_numbers.csv` (modes and variance above noise per scale, mean and SD).
- **Prints:** ℓ_N, the smallest ℓ from which the BCI mean stays within 2 SD of the neutral ensemble.

Check by hand when the data change:

- the shaded crossover band in panel d is fixed at 22–42 m and doesn't follow ℓ_N;
- the "below/above the crossover" labels, axis limits and annotation positions are fixed too.

### `scripts/fig2.tex`: Fig. 2 assembly

Overlays `fig2_ab.pdf` on `fig2_cd.pdf`. Compile it from `figures/` (the Makefile does) to get `figures/fig2.pdf`.
`fig2cd.py` saves an exact 180 × 112 mm canvas so the overlay lines up.

### `scripts/figures_pairs.py`: Fig. 3 and SI pairs robustness

```bash
python scripts/figures_pairs.py [--K 14] [--nu 2e-5] [--cross 20 35] [--inp DIR] [--out DIR]
```

- **Fig. 3 (`panel_pairs`):** (a) F(r) of BCI with the log law a − b ln r fitted for r ≥ 25 m,
  (b) excess over that law with the k most abundant species removed, plus the neutral model.
- **SI (`pairs_robustness`):** (a) excess by size class, (b) the distance below which each species leaves its own log law
  (a 3-SD departure below 25 m, fit over 25–250 m, census 4).
- **Writes:** both figures as pdf, svg and png.
- **Prints:** how many species depart, the median and the IQR.

`--K` and `--nu` pick the neutral-model file, so they must match what `neutral_fr.py` produced. `--cross` sets the shaded
crossover band, which is a manual choice and not computed from ℓ_N.

### `scripts/soil.py`: Fig. 4

Tests how much of the correlation structure soil explains.

- **Reads:** census CSVs, `data/soil/barro_soil_data.xls`, `data/soil/level.png` (map background)
- **Does:**
  1. Bins 13 soil variables on the grid, fills empty cells with the mean, and z-scores them.
  2. Regresses each species' mean abundance on every variable and its square.
  3. Compares the MP-filtered correlations before and after the regression.
- **Writes:** in `figures/panel_2/`: `scatter.svg` and `panel_{b,c,d}_species.svg` (panels only, not assembled)

Unlike the other figures, this script recomputes its matrices from the CSVs. Soil data are on a 20 m grid, so keep
`scale: 20` in `config.yaml`. The example pairs in `pairs = [...]` are row indices into the alphabetical matrix:
check them if `n_species` or the reference census changes.

## Supplementary figures

Outputs go to `figures_supp/`. These scripts recompute from raw data, and their settings are in `main()`, not `config.yaml`.

| Script | Writes | Shows |
|---|---|---|
| `supp/si_figures.py` | `100_species_rank.pdf`, `mp_Q_validity.pdf` | rank of the top 100 species across censuses; Q = T/N against ℓ |
| `supp/distributions.py` | `fig_eig_distributions_{forest,senm}.svg` | eigenvalue distributions at ℓ = 7 and 45 m with the MP edge |
| `supp/empirical_null.py` | `spin_null_vs_mp_{forest,senm}.svg` | spin nulls against the MP density at ℓ = 7 and 42 m |
| `supp/market_mode_check.py` | `market_mode_check.pdf` | whether the top eigenvector at ℓ = 40 m has uniform sign |

- **`empirical_null.py`** runs 1000 spins (seed 0): "rows" shifts each species circularly, keeping its autocorrelation and
  destroying alignment; "both" also shifts columns. It uses census 4 and SENM realization 4.
- **`market_mode_check.py`:** the rank axis of its second panel is **wrong for now**. It assumes rows follow abundance, but
  they are alphabetical (open issue).

## Library (`src/`)

| Module | Provides |
|---|---|
| `config.py` | `load_config()` reads `config.yaml` |
| `utils.py` | `raw2df` (load census or SENM run), `top_species`, `get_dimensions`, `load_spectra`, `load_C`, `ell_N` |
| `compute.py` | `compute_spatial_series`, `compute_forest_C`, `compute_forest_spectrum`, `compute_senm_spectrum`, `marchenko_pastur_bounds`, `marchenko_pastur_filter` |
| `communities.py` | `L_detect_communities`: Laplacian diffusion distances, Ward clustering |
| `figdata.py` | `forest_matrix(scale)`: C, filtered C, MP edge, community order, for Fig. 1c and the SI figure |
| `pairs.py` | `Domain` (F(r) computation), `logfit`, the distance bins `EDGES` and centres `RC` |
| `plotstyle.py` | paper style (`apply`), colours, `diverging()` colormap, `mm_canvas` |

**`pairs.py` method.** Stems go on a raster, and pair counts at distance r come from FFT autocorrelations with zero padding,
so the window isn't wrapped. F(r) is conspecific pairs divided by all pairs at the same r, which corrects for edges.
Bins are 1 m wide up to 10 m and logarithmic up to 500 m.

**MP edge.** `marchenko_pastur_bounds` returns (λ₊, λ₋) with λ± = (1 ± √(N/T))², using continuous T = (W/ℓ)(H/ℓ).
The data use floored bins `int(W/ℓ)·int(H/ℓ)`, so the two differ slightly for ℓ that don't divide the plot (open issue).

## Changing settings

| To change | Edit |
|---|---|
| Species count, censuses, reference census, number of communities | `config.yaml` |
| SENM and neutral F(r) parameters, number of realizations | `senm` block of `config.yaml`: the `Makefile` reads it, and any edit to `config.yaml` reruns the simulation |
| Scales | `range(5, 50)` in `run_pipeline.py`; `fig2cd.py` plots 5–45 m |
| Soil scale | `scale` in `config.yaml` (only `soil.py` uses it) |
| Seed of the simulator | `V2D_SEED` in the `Makefile` |

## Rerunning

`make` tracks data steps through their outputs and the stamp files `data/raw/.csv`, `data/senm/.done`,
`data/processed/.spectra` and `data/processed/pairs/.bci`. Delete one to force that step and everything after it.
The census CSVs aren't linked to the `.rdata` files, so delete `data/raw/.csv` if those change.

Figure targets always rerun, so `make` repeats the slower `fig4` and `si` steps too. Use `make fig1` and the like to avoid that.

## Troubleshooting

| Symptom | Cause |
|---|---|
| `ModuleNotFoundError: src` | `PYTHONPATH` isn't set, or you aren't in the repo root |
| `FileNotFoundError` in `data/processed/` | the data step hasn't run: use `make <figure>` |
| `Warning: ... file not found ... Returning empty DataFrame`, then a crash | `raw2df` found no input: check the CSVs or `data/senm/`, and the `senm` block of `config.yaml` |
| `findfont: Font family ['TeX Gyre Heros'] not found` | harmless: install the font to match the paper |
| `pdflatex: command not found` at `fig2` | install TeX with TikZ; `fig2_ab` and `fig2_cd` are still written |
