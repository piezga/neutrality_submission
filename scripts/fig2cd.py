"""Fig. 2c-d: BCI vs neutral-model spectra at two scales (c) and modes above the noise edge vs scale (d).
Replaces panel_1. Everything is computed from the spectra saved by the pipeline, with the MP edge from
src.compute.marchenko_pastur_bounds, so counts and edge always refer to the same threshold.
Also writes fig2d_numbers.csv (modes above noise, and the variance they carry, vs scale) and prints ell_N.
Panel positions leave the top-left corner empty for panels a-b.
"""
import os, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src import plotstyle as ps
from src.config import load_config
from src.compute import marchenko_pastur_bounds
from src.utils import load_spectra, ell_N

ps.apply(); matplotlib.rcParams['savefig.bbox'] = 'standard'    # exact 180 x 112 mm canvas: fig2.tex overlays a-b
INK, MUT, NOISE, BCI, NEU = ps.INK, ps.MUT, ps.NOISE, ps.BCI, ps.NEU
BCI_L, NEU_L, NEU_LL = ps.BCI_L, ps.NEU_L, ps.NEU_LL
NS = load_config()['n_species']; OUT = 'figures'; os.makedirs(OUT, exist_ok=True)
edge = lambda kind, l: marchenko_pastur_bounds(kind, NS, l)[0]            # lambda_+
have = lambda l: all(os.path.exists(f'data/processed/scale_{l}/{k}_spectra.npy') for k in ('forest', 'senm'))
Ls = np.array([l for l in range(5, 46) if have(l)])

cnt_b, cnt_n, var_b, var_n = [], [], [], []
for l in Ls:
    F, N = load_spectra('forest', l), load_spectra('senm', l)            # (8, NS), (n_real, NS)
    lb, ln = edge('forest', l), edge('senm', l)
    cnt_b.append((F > lb).sum(1)); cnt_n.append((N > ln).sum(1))
    var_b.append(100 * np.where(F > lb, F, 0).sum(1) / NS); var_n.append(100 * np.where(N > ln, N, 0).sum(1) / NS)
cnt_b, cnt_n, var_b, var_n = map(np.array, (cnt_b, cnt_n, var_b, var_n))
lN = ell_N(Ls, cnt_b, cnt_n); print('ell_N from the count of non-noise modes (2 SD):', lN)
np.savetxt(f'{OUT}/fig2d_numbers.csv', np.c_[Ls, cnt_b.mean(1), cnt_b.std(1, ddof=1), cnt_n.mean(1), cnt_n.std(1, ddof=1),
           var_b.mean(1), var_b.std(1, ddof=1), var_n.mean(1), var_n.std(1, ddof=1)], delimiter=',', fmt='%.4f',
           header='scale_m,bci_modes_mean,bci_modes_sd,neutral_modes_mean,neutral_modes_sd,'
                  'bci_var_pct_mean,bci_var_pct_sd,neutral_var_pct_mean,neutral_var_pct_sd')

fig, ax_mm, plabel = ps.mm_canvas(180, 112); W, H = 180, 112
# ---- c: spectra, each divided by the MP edge of its own domain (noise region is lambda/lambda_+ < 1)
plabel(84, H - 1.5, 'c')
c1 = ax_mm(95.5, 63, 39, 45.5); c2 = ax_mm(140, 63, 39, 45.5, sharey=c1)
R = 30; x = np.arange(1, R + 1)
for a, l, tag in ((c1, 5, 'below the crossover'), (c2, 28, 'above the crossover')):
    F = load_spectra('forest', l)[:, :R] / edge('forest', l)               # spectra are saved in descending order
    N = load_spectra('senm', l)[:, :R] / edge('senm', l)
    a.axhspan(1e-3, 1, color=NOISE, lw=0, zorder=0); a.axhline(1, color=MUT, lw=0.6, zorder=1)
    nm, ns = N.mean(0), N.std(0, ddof=1)
    a.fill_between(x, np.clip(nm - 2 * ns, 1e-3, None), nm + 2 * ns, color=NEU_L, lw=0, alpha=0.85, zorder=2)
    a.plot(x, nm, color=NEU, lw=0.9, zorder=3)
    a.fill_between(x, F.min(0), F.max(0), color=BCI_L, lw=0, alpha=0.9, zorder=3)
    fm = F.mean(0); a.plot(x, fm, color=BCI, lw=0.9, zorder=4)
    up = nm > 1; a.plot(x[up], nm[up], 's', ms=2.4, mfc=NEU, mec='white', mew=0.3, zorder=5)
    up = fm > 1; a.plot(x[up], fm[up], 'o', ms=2.6, mfc=BCI, mec='white', mew=0.3, zorder=6)
    j = int(np.where(Ls == l)[0][0]); kb, kn = cnt_b[j].mean(), cnt_n[j].mean(); close = abs(kb - kn) < 5
    for k, col, ha, dx in ((kb, BCI, 'right' if close else 'center', -0.4 if close else 0),
                           (kn, NEU, 'left' if close else 'center', 0.4 if close else 0)):
        a.plot([k + 0.5, k + 0.5], [0.62, 1.0], color=col, lw=0.9, zorder=6)
        a.text(k + 0.5 + dx, 0.58, f'{k:.1f}', color=col, fontsize=6.5, fontweight='bold', ha=ha, va='top')
    a.set_yscale('log'); a.set_xlim(0.3, R + 0.5); a.set_ylim(0.42, 7); a.set_xticks([1, 10, 20, 30]); a.set_xlabel('Eigenvalue rank')
    a.text(0.97, 0.985, f'$\\ell$ = {l} m', transform=a.transAxes, ha='right', va='top', fontsize=7, fontweight='bold')
    a.text(0.97, 0.905, tag, transform=a.transAxes, ha='right', va='top', fontsize=6, color=MUT)
c1.set_ylabel('Eigenvalue / noise edge, $\\lambda/\\lambda_+$'); c1.set_yticks([0.5, 1, 2, 5])
c1.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}')); c1.yaxis.set_minor_formatter(FuncFormatter(lambda v, p: ''))
plt.setp(c2.get_yticklabels(), visible=False); c1.text(R, 0.45, 'noise', fontsize=6, color=MUT, ha='right', va='bottom')

# ---- d: modes above noise vs scale; the crossover is a fading band starting at X0
plabel(36, 52, 'd'); d = ax_mm(48, 9, 96, 40); X0, X1, LMAX = 22, 42, Ls.max()
xs = np.linspace(X0, LMAX + 0.6, 400); al = np.clip((xs - X0) / (X1 - X0), 0, 1) ** 1.3
grad = np.ones((1, len(xs), 4)); grad[..., :3] = matplotlib.colors.to_rgb('#e2e2e2'); grad[..., 3] = al
d.imshow(grad, extent=[X0, LMAX + 0.6, -1e3, 1e3], aspect='auto', zorder=0, interpolation='bilinear')
nm, ns = cnt_n.mean(1), cnt_n.std(1, ddof=1)
d.fill_between(Ls, nm - 2 * ns, nm + 2 * ns, color=NEU_LL, lw=0, zorder=1); d.fill_between(Ls, nm - ns, nm + ns, color=NEU_L, lw=0, zorder=1)
d.plot(Ls, nm, color=NEU, lw=1.0, zorder=3); rng = np.random.default_rng(3)
for k in range(cnt_b.shape[1]):
    d.plot(Ls + rng.uniform(-0.2, 0.2, len(Ls)), cnt_b[:, k], 'o', ms=1.4, mfc=BCI, mec='none', alpha=0.35, zorder=2)
d.plot(Ls, cnt_b.mean(1), '-', color=BCI, lw=1.0, zorder=4)
d.set_xlim(4.2, LMAX + 0.6); d.set_xticks([5, 15, 25, 35, 45]); d.set_xticks(range(5, 46, 5), minor=True)
d.set_xlabel('Scale, $\\ell$ (m)'); d.set_ylabel('Modes above noise'); d.set_ylim(4.5, 24); d.set_yticks([5, 10, 15, 20])
d.text(7.5, 21.8, 'neutral model', color=NEU, fontsize=6.5, ha='left', va='center'); d.text(6, 7.6, 'BCI', color=BCI, fontsize=6.5, ha='left', va='center')
d.annotate('', xy=(44.5, 21.0), xytext=(26, 21.0), arrowprops=dict(arrowstyle='-|>', color=MUT, lw=0.7, mutation_scale=7))
d.text(35.2, 22.3, 'crossover to a neutral description', fontsize=6, ha='center', va='bottom', color=MUT)
fig.savefig(f'{OUT}/fig2_cd.pdf'); fig.savefig(f'{OUT}/fig2_cd.png', dpi=300)
