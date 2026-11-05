"""SI figure: raw vs Marchenko-Pastur-filtered correlation matrix at l = 20 m (same order as Fig. 1c), and the
spectrum of C with the MP edge. Prints the numbers quoted in the caption."""
import os, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src import plotstyle as ps
from src.figdata import forest_matrix

ps.apply(); INK, MUT, NOISE, BCI, DIV = ps.INK, ps.MUT, ps.NOISE, ps.BCI, ps.diverging()
L_C, OUT = 20, 'figures/si_filtered'

M = forest_matrix(L_C)
C, CS, lp, w, keep, K, eig, T, o, CM, bounds = (M[k] for k in
    ('C', 'CMP', 'lam_plus', 'w', 'keep', 'K', 'eig', 'T', 'order', 'labels', 'bounds'))
NSP = len(M['names'])

# numbers for the caption / text
iu = np.triu_indices(NSP, 1)
r_raw_filt = np.corrcoef(C[iu], CS[iu])[0, 1]
var_frac = w[keep].sum() / NSP
Cr_b = C[np.ix_(o, o)]; same = CM[o][:, None] == CM[o][None, :]
off = ~np.eye(NSP, dtype=bool)
print(f'K = {K} modes above lambda_+ = {lp:.3f}; fraction of total variance {var_frac:.3f}')
print(f'Pearson r between raw and filtered off-diagonal entries: {r_raw_filt:.3f}')
print(f'mean raw rho within clusters {Cr_b[same & off].mean():.3f}, between {Cr_b[~same].mean():.3f}')
CSo = CS[np.ix_(o, o)]
print(f'mean filtered rho within {CSo[same & off].mean():.3f}, between {CSo[~same].mean():.3f}')
print(f'raw off-diag sd {C[iu].std():.3f}, filtered off-diag sd {CS[iu].std():.3f}')

FW, FH = 180, 64
fig, ax_mm, plabel = ps.mm_canvas(FW, FH)

S_, Y0, vmax = 46, 12, 0.3
for j, (M, ttl) in enumerate(((C, 'Raw, $C$'), (CS, f'Filtered, $C^{{(s)}}$ ({K} modes)'))):
    x0 = 6 + j * (S_ + 9)
    plabel(x0 - 5.5, FH - 0.8, 'ab'[j])
    a = ax_mm(x0, Y0, S_, S_)
    Mo = M[np.ix_(o, o)].copy(); np.fill_diagonal(Mo, np.nan)
    im = a.imshow(Mo, cmap=DIV, vmin=-vmax, vmax=vmax, interpolation='nearest', rasterized=True)
    for b_ in bounds: a.axhline(b_, color='white', lw=0.45); a.axvline(b_, color='white', lw=0.45)
    for s_ in a.spines.values(): s_.set_visible(False)
    a.set_xticks([]); a.set_yticks([])
    a.text(0.0, 1.025, ttl, transform=a.transAxes, ha='left', va='bottom', fontsize=7)
    a.text(1.0, 1.025, f'$\\ell$ = {L_C} m', transform=a.transAxes, ha='right', va='bottom', fontsize=6.5, color=MUT)
    a.text(-0.025, 0.5, 'Species $i$', transform=a.transAxes, rotation=90, ha='right', va='center', fontsize=6.5)
    a.text(0.5, -0.025, 'Species $j$', transform=a.transAxes, ha='center', va='top', fontsize=6.5)
cb = ax_mm(6 + S_ + 4.5 - S_ * 0.35, 3.0, S_ * 0.7, 1.5)
bar = fig.colorbar(im, cax=cb, orientation='horizontal', ticks=[-0.3, 0, 0.3])
bar.outline.set_visible(False); cb.tick_params(labelsize=6, length=1.8, pad=1, width=0.4)
cb.set_xticklabels(['−0.3', '0', '0.3'])
cb.text(-0.04, 0.5, '$\\rho_{ij}$', transform=cb.transAxes, ha='right', va='center', fontsize=7)

# c: spectrum and MP edge
plabel(121, FH - 0.8, 'c')
c = ax_mm(132, Y0 + 2, 45, S_ - 2)
rk = np.arange(1, NSP + 1)
c.axhspan(1e-3, lp, color=NOISE, lw=0, zorder=0)
c.axhline(lp, color=MUT, lw=0.6, zorder=1)
c.plot(rk, eig, color=BCI, lw=0.9, zorder=2)
c.plot(rk[:K], eig[:K], 'o', ms=2.6, mfc=BCI, mec='white', mew=0.3, zorder=3)
c.plot(rk[K:], eig[K:], 'o', ms=1.8, mfc='white', mec=BCI, mew=0.4, zorder=3)
c.set_xscale('log'); c.set_yscale('log'); c.set_xlim(0.8, 110)
c.set_ylim(eig.min() * 0.8, eig.max() * 1.5)
c.xaxis.set_major_locator(FixedLocator([1, 10, 100])); c.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
c.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5])); c.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
c.xaxis.set_minor_formatter(NullFormatter()); c.yaxis.set_minor_formatter(NullFormatter())
c.set_xlabel('Eigenvalue rank'); c.set_ylabel('Eigenvalue, $\\lambda$', labelpad=1)
c.text(1.05, eig.min() * 0.92, f'noise: $\\lambda<\\lambda_+$ = {lp:.2f}', ha='left', va='bottom', fontsize=6.2, color=MUT)
c.text(K * 1.25, eig[K - 1] * 1.12, f'{K} modes kept', ha='left', va='bottom', fontsize=6.3, color=BCI)

fig.savefig(f'{OUT}.pdf', dpi=600); fig.savefig(f'{OUT}.png', dpi=450)
print('saved', OUT)
