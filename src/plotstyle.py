"""Shared figure style for the manuscript (replaces the rcParams block copied in every figure script)."""
import matplotlib as mpl
import matplotlib.pyplot as plt

INK, MUT, BAND, NOISE = '#1a1a1a', '#8c8c8c', '#ececec', '#ececec'
BCI, NEU = '#009E73', '#D55E00'                 # Okabe-Ito: real data / neutral model
BCI_L, NEU_L, NEU_LL = '#b3e2d1', '#f3c3a3', '#fae3d3'
DIV_COLORS = ['#24406e', '#5f7aa6', '#b4c1d6', '#efeeea', '#e3c2b2', '#c0765c', '#8a2f22']
YEAR = dict(zip(range(1, 9), (1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015)))   # BCI census years


def apply():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['TeX Gyre Heros', 'Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
        'axes.linewidth': 0.5, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
        'legend.frameon': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
        'mathtext.fontset': 'custom', 'mathtext.rm': 'TeX Gyre Heros',
        'mathtext.it': 'TeX Gyre Heros:italic', 'mathtext.bf': 'TeX Gyre Heros:bold'})


def diverging():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list('slate_terra', DIV_COLORS)


def mm_canvas(W, H):
    """Figure of W x H millimetres. Returns fig, ax_mm(x0, y0, w, h, **kw), plabel(x, y, text)."""
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax_mm = lambda x0, y0, w, h, **kw: fig.add_axes([x0 / W, y0 / H, w / W, h / H], **kw)
    plabel = lambda x, y, t: fig.text(x / W, y / H, t, fontsize=8, fontweight='bold', va='top', ha='left')
    return fig, ax_mm, plabel
