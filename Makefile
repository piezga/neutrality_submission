# Every figure of the paper:  make        (make -j4 runs independent steps in parallel)
PY   ?= python
PROC := data/processed
export PYTHONPATH := $(CURDIR)
export V2D_SEED   := 2026
export MPLBACKEND := Agg
export PYTHONUNBUFFERED := 1
export TQDM_DISABLE := 1

# SENM parameters live in config.yaml only
cfg  = $(shell sed -n 's/^ *$(1): *"\?\([^"]*\)"\?$$/\1/p' config.yaml)
K    := $(call cfg,kernel)
NU   := $(call cfg,nu)
NREL := $(call cfg,num_realizations)
FR   := $(PROC)/pairs/cal_K$(K)_nu$(shell printf %g $(NU)).pkl

all: fig1 fig2 fig3 fig4 si
.PHONY: all fig1 fig2 fig3 fig4 si

# ---- data
data/raw/.csv: scripts/convert_data.py
	$(PY) $< && touch $@

neutral_model/v2d: neutral_model/*.c neutral_model/*.h
	bash neutral_model/build.sh

# SENM realizations for the spectra (config.yaml: senm). v2d exits 1 on success, so check its output instead.
data/senm/.done: neutral_model/v2d config.yaml
	rm -rf data/senm && mkdir -p data/senm && cd data/senm && (../../neutral_model/v2d 500 500 $(NU) $(NREL) $(K) 2>/dev/null || true)
	test $$(ls data/senm/AConf_* | wc -l) -eq $(NREL) && touch $@

$(PROC)/.spectra: data/raw/.csv data/senm/.done scripts/run_pipeline.py src/compute.py src/utils.py config.yaml
	$(PY) scripts/run_pipeline.py && touch $@

$(PROC)/pairs/.bci: data/raw/.csv scripts/compute_pairs.py src/pairs.py
	for w in all species size; do $(PY) scripts/compute_pairs.py $$w || exit 1; done && touch $@

$(FR): neutral_model/v2d scripts/neutral_fr.py src/pairs.py
	$(PY) scripts/neutral_fr.py --side 1500 --spacing 1 --K $(K) --nu $(NU)

# ---- figures (figures/ and figures_supp/)
fig1: $(PROC)/.spectra
	mkdir -p figures && $(PY) scripts/fig1.py && $(PY) scripts/si_filtered.py

fig2: $(PROC)/.spectra
	mkdir -p figures && $(PY) scripts/fig2ab.py && $(PY) scripts/fig2cd.py
	cd figures && pdflatex -interaction=batchmode ../scripts/fig2.tex >/dev/null && rm -f fig2.aux fig2.log

fig3: $(PROC)/pairs/.bci $(FR)
	$(PY) scripts/figures_pairs.py --K $(K) --nu $(NU)

fig4: data/raw/.csv
	$(PY) scripts/soil.py

si: data/raw/.csv data/senm/.done
	mkdir -p figures_supp && for s in si_figures distributions empirical_null market_mode_check; do $(PY) supp/$$s.py || exit 1; done
