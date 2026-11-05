"""Fig. 1: a) map of three example species with the l = 100 m grid, b) their occupation vectors on that grid,
c) correlation matrix of the 100 most abundant species at l = 20 m, ordered by the communities of the MP-filtered matrix."""
import os, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Rectangle, FancyArrowPatch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src import plotstyle as ps
from src.figdata import forest_matrix
from src.utils import raw2df, get_dimensions

ps.apply(); INK, MUT, PLOT_BG, YEAR, DIV = ps.INK, ps.MUT, '#f3f3f0', ps.YEAR, ps.diverging()
SP = [('tet2pa', 'Tetragastris panamensis', '#3B6FB6'),      # example species: code, label, colour
      ('cou2cu', 'Coussarea curvigemmia', '#B8336A'),
      ('taline', 'Talisia nervosa', '#D9A21B')]
MAP_C, L_MAP, L_C = 8, 100, 20                                 # census shown in a-b, grid of a-b, scale of c
(W, H), OUT = get_dimensions('forest'), 'figures/fig1'

M = forest_matrix(L_C); NSP = len(M['names'])
Cr = M['C'][np.ix_(M['order'], M['order'])]; np.fill_diagonal(Cr, np.nan); bounds = M['bounds']
df = raw2df('forest', MAP_C, [c for c, _, _ in SP])


def series(df, l):
    nx_, ny_ = W // l, H // l
    xb = (df.x // l).astype(int).clip(0, nx_ - 1); yb = (df.y // l).astype(int).clip(0, ny_ - 1)
    return xb * ny_ + yb, nx_ * ny_


FW, FH = 180, 60
fig, ax_mm, plabel = ps.mm_canvas(FW, FH)
def flow(x0, x1, y, txt):
    fig.patches.append(FancyArrowPatch((x0 / FW, y / FH), (x1 / FW, y / FH), transform=fig.transFigure,
                       arrowstyle='-|>,head_length=2.2,head_width=1.2', lw=0.6, color=MUT, figure=fig))
    fig.text((x0 + x1) / 2 / FW, (y + 1.4) / FH, txt, ha='center', va='bottom', fontsize=5.6, color=MUT,
             linespacing=1.1)

# a: map
plabel(0.5, FH - 0.8, 'a')
AX, AY, AW = 8, 9.5, 74
a = ax_mm(AX, AY, AW, AW / 2)
a.add_patch(Rectangle((0, 0), W, H, color=PLOT_BG, lw=0, zorder=0))
for code, lab, col in SP:
    p = df[df.name == code]
    a.scatter(p.x, p.y, s=1.0, c=col, lw=0, alpha=0.9, zorder=2, rasterized=True)
for x in range(0, W + 1, L_MAP): a.axvline(x, color='white', lw=0.9, zorder=3)
for y in range(0, H + 1, L_MAP): a.axhline(y, color='white', lw=0.9, zorder=3)
cx, cy = W - L_MAP, H - L_MAP          # top-right cell, outlined
a.add_patch(Rectangle((cx, cy), L_MAP, L_MAP, fill=False, ec='white', lw=2.4, zorder=4))
a.add_patch(Rectangle((cx, cy), L_MAP, L_MAP, fill=False, ec=INK, lw=1.1, zorder=5))
dim = dict(arrowstyle='<->', lw=0.5, color=INK, shrinkA=0, shrinkB=0)
a.annotate('', xy=(cx, H + 22), xytext=(W, H + 22), arrowprops=dim, annotation_clip=False)
a.text(cx + L_MAP / 2, H + 30, r'$\ell$', ha='center', va='bottom', fontsize=7, clip_on=False)
a.annotate('', xy=(W + 22, cy), xytext=(W + 22, H), arrowprops=dim, annotation_clip=False)
a.text(W + 32, cy + L_MAP / 2, r'$\ell$', ha='left', va='center', fontsize=7, clip_on=False)
a.set_xlim(0, W); a.set_ylim(0, H); a.set_aspect('equal')
a.set_xticks(range(0, 1001, 200)); a.set_yticks(range(0, 501, 100))
a.set_xlabel('x (m)', labelpad=1.5); a.set_ylabel('y (m)', labelpad=1.5)
for sp in ('left', 'bottom'): a.spines[sp].set_visible(False)
a.tick_params(length=0, pad=2)
lg = ax_mm(AX, AY + AW / 2 + 3.0, AW, 3.5); lg.axis('off'); lg.set_xlim(0, AW); lg.set_ylim(0, 1)
lx = 0.8
for code, lab, col in SP:
    lg.scatter([lx], [0.5], s=9, c=col, lw=0)
    g, sp_ = lab.split()
    lg.text(lx + 1.6, 0.5, f'$\\it{{{g[0]}.\\ {sp_}}}$', fontsize=6.5, va='center', color=INK)
    lx += 22.5
lg.text(0.0, 1.75, f'{YEAR[MAP_C]} census', ha='left', va='center', fontsize=6.5, color=MUT)

flow(83.5, 90, 28, 'count\nper cell')

# b: occupation vectors on the l = 100 m grid
plabel(90.5, FH - 0.8, 'b')
k4, T4 = series(df, L_MAP)
bx, bw, bh, gap = 100, 28.5, 10.6, 3.2
for j, (code, lab, col) in enumerate(SP):
    y0 = AY + (2 - j) * (bh + gap)
    b = ax_mm(bx, y0, bw, bh)
    n = np.bincount(k4[df.name == code], minlength=T4)
    kc = (cx // L_MAP) * (H // L_MAP) + cy // L_MAP + 1   # cell outlined in a
    b.axvspan(kc - 0.6, kc + 0.6, color='#e4e4e4', lw=0, zorder=0)
    b.bar(np.arange(1, T4 + 1), n, width=0.78, color=col, lw=0, zorder=2)
    b.set_xlim(0.2, T4 + 0.8); top = int(n.max()); b.set_ylim(0, top * 1.05)
    b.set_yticks([0, top]); b.tick_params(axis='y', labelsize=5.8, pad=1, length=1.8)
    b.set_xticks([1, 10, 20, 30, 40, 50]); b.tick_params(axis='x', length=1.8, pad=1.5)
    if j < 2: b.set_xticklabels([])
    b.text(0.03, 1.0, f'$n^{{{j + 1}}}_k$', transform=b.transAxes, fontsize=6.5, color=col, va='top', ha='left')
    if j == 0: b.text(kc + 0.6, top * 1.05, 'cell outlined in a', ha='right', va='bottom', fontsize=5.4, color=MUT)
    if j == 2: b.set_xlabel('Cell, $k$', labelpad=1)
fig.text(92.6 / FW, (AY + (3 * bh + 2 * gap) / 2) / FH, 'Stems per cell', rotation=90,
         ha='center', va='center', fontsize=7)

flow(130, 136.5, 28, 'Pearson\n$\\rho_{ij}$')

# c: correlation matrix
plabel(137, FH - 0.8, 'c')
S_ = 38.5
cm = ax_mm(141.3, 13.5, S_, S_)
vmax = 0.3
M_ = Cr.copy()
M_[np.triu_indices(NSP, 1)] = np.nan
im = cm.imshow(M_, cmap=DIV, vmin=-vmax, vmax=vmax, interpolation='nearest', rasterized=True)
for s_ in cm.spines.values(): s_.set_visible(False)
cm.set_xticks([]); cm.set_yticks([])
for b_ in bounds:      # community blocks: thin white separators inside the lower triangle
    cm.plot([-0.5, b_], [b_, b_], color='white', lw=0.45); cm.plot([b_, b_], [b_, NSP - 0.5], color='white', lw=0.45)
cm.set_xlim(-0.5, NSP - 0.5); cm.set_ylim(NSP - 0.5, -0.5)
cm.set_position([141.3 / FW, AY / FH, S_ / FW, S_ / FH])
cm.text(0.98, 0.985, f'{NSP} species', transform=cm.transAxes, ha='right', va='top', fontsize=6.5)
cm.text(0.98, 0.91, f'$\\ell$ = {L_C} m', transform=cm.transAxes, ha='right', va='top', fontsize=6.5, color=MUT)
cm.text(-0.025, 0.5, 'Species $i$', transform=cm.transAxes, rotation=90, ha='right', va='center', fontsize=6.5)
cm.text(0.5, -0.025, 'Species $j$', transform=cm.transAxes, ha='center', va='top', fontsize=6.5)
cb = cm.inset_axes([0.50, 0.72, 0.46, 0.035])
bar = fig.colorbar(im, cax=cb, orientation='horizontal', ticks=[-0.3, 0, 0.3])
bar.outline.set_visible(False); cb.tick_params(labelsize=6, length=1.8, pad=1, width=0.4)
cb.set_xticklabels(['\u22120.3', '0', '0.3'])
cb.text(0.5, 1.6, 'Correlation, $\\rho_{ij}$', transform=cb.transAxes, va='bottom', ha='center', fontsize=6.5)

for e in ('pdf', 'svg'): fig.savefig(f'{OUT}.{e}', dpi=600)
fig.savefig(f'{OUT}.png', dpi=450)
print('saved', OUT)
