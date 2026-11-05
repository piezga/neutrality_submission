# Neutrality across scales in the BCI forest

Code for the paper's figures: species correlation spectra in the Barro Colorado Island 50-ha plot
compared with a spatially explicit neutral model (SENM).

## Quick start

```bash
git clone <repo-url> neutrality && cd neutrality
conda env create -f environment.yml
conda activate neutrality
make -j4
```

Figures land in `figures/` (main text) and `figures_supp/` (SI).

## On a workstation (no terminal attached)

`make` needs no interaction: plots use the `Agg` backend, progress bars are off and output is unbuffered.
Point `PY` at the environment's interpreter instead of activating it:

```bash
nohup make -j4 PY=$(conda run -n neutrality which python) > make.log 2>&1 &
```

A failed or interrupted step leaves no stamp, so rerunning `make` resumes from that step.

## Requirements

- conda (or mamba): Python packages are in `environment.yml`
- a C compiler (`cc`): builds the neutral-model simulator
- `pdflatex` with TikZ: assembles Fig. 2
- optional: the TeX Gyre Heros font (otherwise DejaVu Sans is used)

## Data

| Path | Content | Source |
|---|---|---|
| `data/raw/bci.tree<c>.rdata` | BCI censuses 1–8 (1982–2015) | Condit et al., Dryad, https://doi.org/10.15146/5xcp-0d46 |
| `data/soil/barro_soil_data.xls` | soil nutrients on a 20 m grid | https://forestgeo.si.edu/sites/neotropics/barro-colorado-island/bci-maps/bci-soil-maps |

If the `.rdata` files are missing, download them from source into `data/raw/` and data/soil/.
Everything else (`data/raw/*.csv`, `data/senm/`, `data/processed/`) is generated.

## What `make` does

| Step | Command | Output |
|---|---|---|
| Census tables, alive stems only | `scripts/convert_data.py` | `data/raw/bci.tree<c>.csv` |
| Build simulator | `neutral_model/build.sh` | `neutral_model/v2d` |
| SENM, 100 realizations (500², K = 14, ν = 2e-5) | `v2d` | `data/senm/` |
| Spectra and BCI correlation matrix, ℓ = 5–49 m | `scripts/run_pipeline.py` | `data/processed/scale_<ℓ>/` |
| Pair probability F(r), BCI | `scripts/compute_pairs.py` | `data/processed/pairs/` |
| F(r), neutral model (1500², 1 m/site) | `scripts/neutral_fr.py` | `data/processed/pairs/cal_K14_nu2e-05.pkl` |
| Fig. 1, SI filtered matrix | `scripts/fig1.py`, `si_filtered.py` | `figures/fig1.*`, `si_filtered.*` |
| Fig. 2 | `scripts/fig2ab.py`, `fig2cd.py`, `fig2.tex` | `figures/fig2.pdf`, `fig2d_numbers.csv` |
| Fig. 3, SI pairs robustness | `scripts/figures_pairs.py` | `figures/panel_pairs.*`, `pairs_robustness.*` |
| Fig. 4 (soil) | `scripts/soil.py` | `figures/panel_2/` |
| SI figures | `supp/*.py` | `figures_supp/` |

Data steps rerun only when their inputs change; figure targets always rerun. See [documentation.md](documentation.md) for what each script does.

## Useful targets

```bash
make fig2              # one figure (fig1, fig2, fig3, fig4, si)
make PY=/path/python   # use a specific interpreter instead of the active one
```

To recompute a step, delete its output or its stamp file (`data/raw/.csv`, `data/senm/.done`,
`data/processed/.spectra`, `data/processed/pairs/.bci`).

## Configuration

Analysis parameters (species count, censuses, SENM size and parameters) are in `config.yaml`.
The simulator seed is set by `V2D_SEED` in the `Makefile` (2026).

## Layout

```
src/            analysis library (spectra, MP bounds, communities, F(r), plot style)
scripts/        pipeline steps and main-text figures
supp/           SI figures
neutral_model/  C code of the voter model with speciation (v2d)
config.yaml     parameters
Makefile        the pipeline
```

## Use of Large Language Models
P.Z. made use of Claude Sonnet 5.5 to troubleshoot and optimize individual scripts, and of Claude Opus 5.5 to write the documentation, improve code readability and assemble the pipeline into a single Makefile.
