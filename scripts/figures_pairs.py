"""Fig. 3 (F(r) in BCI and the neutral model) and the SI pairs-robustness figure.

    python scripts/figures_pairs.py [--K 14] [--nu 2e-5] [--cross 20 35]

Inputs in data/processed/pairs: fofr_bci.pkl, perspecies_c{4,8}.pkl, fofr_size.pkl (compute_pairs.py)
and cal_K<K>_nu<nu>.pkl (neutral_fr.py). Outputs in figures/: panel_pairs.*, pairs_robustness.*
"""
import os, sys, pickle, argparse
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src import plotstyle as ps
from src.pairs import logfit

ap = argparse.ArgumentParser()
ap.add_argument('--K', type=int, default=14); ap.add_argument('--nu', type=float, default=2e-5)
ap.add_argument('--cross', type=float, nargs=2, default=(20, 35), help='crossover band shaded in r (m)')
ap.add_argument('--inp', default='data/processed/pairs'); ap.add_argument('--out', default='figures')
args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
ps.apply(); INK, MUT, BAND, BCI, NEU = ps.INK, ps.MUT, ps.BAND, ps.BCI, ps.NEU
CROSS = tuple(args.cross); XT = [1, 2, 5, 10, 20, 50, 100, 200, 500]
load = lambda f: pickle.load(open(os.path.join(args.inp, f), 'rb'))


def xaxis(x):
    x.set_xscale('log'); x.set_xlim(0.8, 560); x.set_xticks(XT); x.set_xticklabels([str(t) for t in XT])
    x.set_xlabel('Distance between two stems, $r$ (m)')


def save(fig, name):
    for e in ('pdf', 'svg'): fig.savefig(f'{args.out}/{name}.{e}')
    fig.savefig(f'{args.out}/{name}.png', dpi=400)


def fig3():
    fig, ax = plt.subplots(1, 2, figsize=(180 / 25.4, 62 / 25.4), gridspec_kw=dict(wspace=0.32))
    o = load('fofr_bci.pkl'); rc = o['rc']; F = np.array([o[c] for c in range(1, 9)]); m = F.mean(0)
    x = ax[0]; x.axvspan(*CROSS, color=BAND, lw=0, zorder=0)
    x.fill_between(rc, F.min(0), F.max(0), color=BCI, alpha=0.22, lw=0)
    x.plot(rc, m, color=BCI, lw=1.2, marker='o', ms=2.2, mew=0, label='BCI, eight censuses (mean and range)')
    fit = rc >= 25; b, a = np.polyfit(np.log(rc[fit]), m[fit], 1); rr = np.geomspace(0.8, 560, 100)
    x.plot(rr, a + b * np.log(rr), color=INK, lw=0.8, ls=(0, (3, 2)), label=r'$a-b\,\ln r$, fitted for $r\geq 25$ m')
    xaxis(x); x.set_ylim(0.04, 0.1); x.set_ylabel('Probability that two stems are conspecific, $F(r)$')
    x.legend(loc='upper right', fontsize=6.5, handlelength=2.2, borderaxespad=0.2)

    x = ax[1]; x.axvspan(*CROSS, color=BAND, lw=0, zorder=0); x.axhline(0, color=MUT, lw=0.5, zorder=1)
    P = [load(f'perspecies_c{c}.pkl') for c in (4, 8)]
    ks = [0, 2, 5, 10, 20]; cmap = LinearSegmentedColormap.from_list('g', ['#00452f', '#009E73', '#9fdcc4'])
    for k, col in zip(ks, [cmap(t) for t in np.linspace(0, 0.72, len(ks))]):
        R = []
        for p in P:
            sps = sorted(p['n'], key=lambda s: -p['n'][s]); Fk = sum(p['pairs'][s] for s in sps[k:]) / p['all']
            R.append(100 * (Fk / logfit(p['rc'], Fk) - 1))
        R = np.array(R); mR = R.mean(0)
        x.fill_between(rc, R.min(0), R.max(0), color=col, alpha=0.25, lw=0, zorder=2)
        x.plot(rc, mR, color=col, lw=1.1, zorder=3); x.text(0.62, mR[0], f'{k}', ha='left', va='center', fontsize=6.5, color=col)
    x.text(0.62, 186, '$k$', ha='left', va='center', fontsize=6.5, color=INK)
    n = load(f'cal_K{args.K}_nu{args.nu:g}.pkl')
    R = np.array([100 * (Fn / logfit(n['rc'], Fn) - 1) for Fn in n['F']])
    x.fill_between(n['rc'], np.nanmin(R, 0), np.nanmax(R, 0), color=NEU, alpha=0.22, lw=0, zorder=2)
    x.plot(n['rc'], np.nanmean(R, 0), color=NEU, lw=1.1, ls=(0, (4, 2)), zorder=3)
    x.text(4, -32, 'neutral model', color=NEU, fontsize=6.5, ha='left', va='center')
    x.text(55, 120, 'BCI without its $k$\nmost abundant species', color='#00704f', fontsize=6.5, ha='left', va='center')
    x.text(np.sqrt(CROSS[0] * CROSS[1]), 172, 'crossover', ha='center', va='center', fontsize=6.5, color=MUT)
    xaxis(x); x.set_xlim(0.5, 560); x.set_ylim(-40, 190); x.set_ylabel('Excess over the logarithmic law (%)')
    for i, xx in enumerate(ax): xx.text(-0.19, 1.04, 'ab'[i], transform=xx.transAxes, fontsize=8, fontweight='bold', va='bottom')
    save(fig, 'panel_pairs')


def si_pairs():
    fig, ax = plt.subplots(1, 2, figsize=(180 / 25.4, 62 / 25.4), gridspec_kw=dict(wspace=0.32))
    s = load('fofr_size.pkl'); rc = s['rc']
    x = ax[0]; x.axvspan(*CROSS, color=BAND, lw=0, zorder=0); x.axhline(0, color=MUT, lw=0.5)
    for lab, col, name, ty in [('10-49mm', BCI, 'saplings, 1–5 cm dbh', 40), ('>=100mm', '#0072B2', 'trees, ≥10 cm dbh', 27)]:
        R = np.array([100 * (s[(c, lab)] / logfit(rc, s[(c, lab)]) - 1) for c in (4, 8)])
        x.fill_between(rc, R.min(0), R.max(0), color=col, alpha=0.22, lw=0); x.plot(rc, R.mean(0), color=col, lw=1.1)
        n = int(np.mean([s[(c, lab, 'n')] for c in (4, 8)])); x.text(45, ty, f'{name}\n({n / 1000:.0f},000 stems)', color=col, fontsize=6.5, va='center')
    xaxis(x); x.set_ylabel('Excess over the logarithmic law (%)'); x.set_ylim(-15, 48)

    o = load('perspecies_c4.pkl'); rc = o['rc']; al = o['all']; fit = (rc >= 25) & (rc <= 250); deps = []
    for sp in o['n']:
        F = o['pairs'][sp] / al; b, a = np.polyfit(np.log(rc[fit]), F[fit], 1); res = F - (a + b * np.log(rc)); sd = res[fit].std()
        ex = np.where((res > 3 * sd) & (rc < 25))[0]; deps.append(rc[ex].max() if len(ex) else np.nan)
    deps = np.array(deps); x = ax[1]
    x.hist(deps[np.isfinite(deps)], bins=np.arange(0, 27.5, 2.5), color=BCI, alpha=0.85, edgecolor='white', lw=0.5)
    md = np.nanmedian(deps); x.axvline(md, color=INK, lw=0.8, ls=(0, (3, 2)))
    x.text(md + 0.6, x.get_ylim()[1] * 0.97, f'median {md:.0f} m', ha='left', va='top', fontsize=6.5)
    x.set_xlabel('Distance below which a species departs\nfrom its own logarithmic law (m)'); x.set_ylabel('Number of species'); x.set_xlim(0, 26)
    for i, xx in enumerate(ax): xx.text(-0.19, 1.04, 'ab'[i], transform=xx.transAxes, fontsize=8, fontweight='bold', va='bottom')
    save(fig, 'pairs_robustness')
    print('species with departure', np.isfinite(deps).sum(), 'of', len(deps), 'median', md, 'IQR', np.nanpercentile(deps, [25, 75]))


fig3(); si_pairs()
