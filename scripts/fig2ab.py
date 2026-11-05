"""Fig. 2a-b cartoon: a) voter update with a square kernel, b) one patch of a toy SENM coarse-grained
at growing block sizes. Self-contained (simulates its own small lattice). Writes figures/fig2_ab.*"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.collections import PatchCollection

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
SEED      = 11
FX, FY    = 120, 60      # field the region is cut from
NU        = 2.5e-5
KERNEL    = 3            # square kernel half-width of the toy field
SWEEPS    = 900
N_SHOWN   = 5

# Panel b shows THE SAME square region five times over, side by side, each tile
# coarse-grained at its own block size. Same region every tile, so the only thing
# that changes across the strip is the resolution -- never the piece of forest.
REGION    = 40           # side of the region, in fine sites
BLOCKS    = [1, 2, 4, 5, 10]   # each must divide REGION
LN_AFTER  = 2            # the crossover falls at the seam after tile 2

# Okabe-Ito subset. Validated over ALL pairs on a white surface:
# deutan dE 11.0, tritan 8.6, normal 15.6 -- all above threshold.
SPECIES = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9"]
NEWSP   = "#7A3B9E"
OTHER   = "#DCDCDC"
INK     = "#141414"
RULE    = "#3C3C3C"
MUTED   = "#767676"
VACANT  = "#FFFFFF"

mpl.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["TeX Gyre Heros", "Helvetica", "Arial",
                         "Liberation Sans", "DejaVu Sans"],
    "font.size":        7,
    "text.color":       INK,
    "mathtext.fontset": "custom",
    "mathtext.rm":      "TeX Gyre Heros",
    "mathtext.it":      "TeX Gyre Heros:italic",
    "mathtext.bf":      "TeX Gyre Heros:bold",
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "savefig.bbox":     "tight",
    "savefig.pad_inches": 0.015,
})


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------

def neighbour_offsets(k):
    dx, dy = np.meshgrid(np.arange(-k, k + 1), np.arange(-k, k + 1))
    keep = ~((dx == 0) & (dy == 0))
    return dx[keep], dy[keep]

def simulate(lx, ly, nu, k, sweeps, rng, chunk_frac=50):
    """
    Voter model with speciation, one individual per site, periodic boundaries.
    Updated in small random chunks: within a chunk the sources are read from the
    current state, so it is asynchronous up to collisions that are rare at this
    chunk size and invisible in the resulting pattern.
    """
    n = lx * ly
    lat = np.arange(n).reshape(ly, lx)
    nxt = n
    odx, ody = neighbour_offsets(k)
    m = max(1, n // chunk_frac)

    for _ in range(sweeps):
        for _ in range(chunk_frac):
            xs = rng.integers(0, lx, m)
            ys = rng.integers(0, ly, m)
            spec = rng.random(m) < nu
            j = rng.integers(0, len(odx), m)
            src = lat[(ys + ody[j]) % ly, (xs + odx[j]) % lx]
            ns = spec.sum()
            if ns:
                src = src.copy()
                src[spec] = np.arange(nxt, nxt + ns)
                nxt += ns
            lat[ys, xs] = src
    return lat


def recolour(lat, n_shown):
    ids, counts = np.unique(lat, return_counts=True)
    top = ids[np.argsort(counts)[::-1][:n_shown]]
    idx = np.full(lat.shape, -1, dtype=int)
    for i, sp in enumerate(top):
        idx[lat == sp] = i
    return idx


def colour_of(i):
    return OTHER if i < 0 else SPECIES[i % len(SPECIES)]


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------
def cells(ax, idx, x0, y0, cell, ec="white", lw=0.0, override=None,
          alphas=None, zorder=2):
    """alphas, if given, is a per-cell array the same shape as idx."""
    ny, nx = idx.shape
    rects, cols = [], []
    for r in range(ny):
        for c in range(nx):
            col = override.get((r, c)) if override else None
            if col is None:
                col = colour_of(idx[r, c])
            if alphas is not None:
                col = mpl.colors.to_rgba(col, alphas[r, c])
            rects.append(Rectangle((x0 + c * cell, y0 + (ny - 1 - r) * cell),
                                   cell, cell))
            cols.append(col)
    ax.add_collection(PatchCollection(rects, facecolor=cols, edgecolor=ec,
                                      linewidth=lw, zorder=zorder))
    return nx * cell, ny * cell


def hairline(ax, x, y, w, h, color=RULE, lw=0.5, z=9):
    ax.add_patch(Rectangle((x, y), w, h, fill=False, ec=color, lw=lw, zorder=z))


# --------------------------------------------------------------------------
# Panel a  --  shared neighbourhood, then a branch
# --------------------------------------------------------------------------
def panel_a(ax, patch, focal, donor):
    ax.set_axis_off()
    ny, nx = patch.shape
    fr, fc = focal
    dr, dc = donor
    C = 1.0                                   # cell size of the main neighbourhood
    Z = 1.9                                   # cell size in the magnified outcomes

    # ---- the shared neighbourhood, drawn once; everything outside the kernel
    #      is dimmed so the eye goes straight to the eligible region
    ring = KERNEL + 0.5
    rr, cc = np.mgrid[0:ny, 0:nx]
    inside = np.maximum(abs(rr - fr), abs(cc - fc)) <= KERNEL     # square kernel
    alph = np.where(inside, 1.0, 0.38)
    w, h = cells(ax, patch, 0, 0, C, ec="white", lw=0.45,
                 override={(fr, fc): VACANT}, alphas=alph)
    hairline(ax, 0, 0, w, h)

    xf, yf = fc + 0.5, (ny - 1 - fr) + 0.5    # focal centre
    xd, yd = dc + 0.5, (ny - 1 - dr) + 0.5    # donor centre

    ax.add_patch(Rectangle((xf - ring, yf - ring), 2*ring, 2*ring, 
                           fill=False, ec=INK, lw=0.9, ls=(0, (2.4, 1.7)), zorder=6))

    # the vacancy, and the donor it will be filled from
    ax.add_patch(Rectangle((fc, ny - 1 - fr), C, C, facecolor=VACANT,
                           ec=INK, lw=0.9, zorder=7))
    ax.add_patch(Rectangle((dc, ny - 1 - dr), C, C,
                           facecolor=colour_of(patch[dr, dc]),
                           ec=INK, lw=0.9, zorder=7))
    ax.add_patch(FancyArrowPatch((xd, yd), (xf, yf), arrowstyle="-|>",
                                 mutation_scale=6.5, lw=0.85, color=INK,
                                 zorder=8, shrinkA=7, shrinkB=8,
                                 connectionstyle="arc3,rad=-0.22"))
    # radius tick
    ang = np.deg2rad(58)
    ax.plot([xf, xf + ring], [yf, yf],
            lw=0.7, color=INK, zorder=8, solid_capstyle="butt")
    ax.text(xf + ring * np.cos(ang) * 0.50 + 0.52,
            yf + ring * np.sin(ang) * 0.50 - 0.28, r"$k$",
            fontsize=7.5, ha="left", va="center", zorder=9)

    ax.text(w / 2, -1.15, "a tree dies; the site is refilled",
            ha="center", va="top", fontsize=6.6, color=MUTED, style="italic")

    # ---- the branch
    x_br = w + 2.2
    y_hi, y_lo = h * 0.93, h * 0.07
    ax.plot([x_br, x_br], [y_lo, y_hi], lw=0.7, color=RULE, zorder=6)
    ax.plot([w + 0.6, x_br], [yf, yf], lw=0.7, color=RULE, zorder=6)
    for yb, pr in ((y_hi, r"$1-\nu$"), (y_lo, r"$\nu$")):
        ax.add_patch(FancyArrowPatch((x_br, yb), (x_br + 2.1, yb),
                                     arrowstyle="-|>", mutation_scale=6.5,
                                     lw=0.7, color=RULE, zorder=6))
        ax.text(x_br + 1.05, yb + 0.45, pr, fontsize=7.2, ha="center",
                va="bottom", fontweight="bold")

    # 3x3 magnified crops of the outcome
    r0, c0 = fr - 1, fc - 1
    crop = patch[r0:r0 + 3, c0:c0 + 3]
    x_out = x_br + 2.9
    outs = [
        (y_hi, colour_of(patch[dr, dc]), INK, "local replacement",
         "a copy of a neighbour drawn\nfrom inside the kernel"),
        (y_lo, NEWSP, NEWSP, "speciation",
         "an individual of a species\nnot yet present"),
    ]
    for yb, fill, ec, name, blurb in outs:
        y0 = yb - 1.5 * Z
        cells(ax, crop, x_out, y0, Z, ec="white", lw=0.5,
              override={(1, 1): fill})
        hairline(ax, x_out, y0, 3 * Z, 3 * Z)
        ax.add_patch(Rectangle((x_out + Z, y0 + Z), Z, Z, fill=False,
                               ec=ec, lw=1.2, zorder=8))
        tx = x_out + 3 * Z + 1.1
        ax.text(tx, yb + 2.15, name, fontsize=7.4, ha="left", va="center",
                fontweight="bold")
        ax.text(tx, yb - 0.60, blurb, fontsize=6.3, ha="left", va="center",
                color=MUTED, linespacing=1.5)

    ax.set_xlim(-0.8, x_out + 3 * Z + 19.0)
    ax.set_ylim(-3.0, h + 2.6)
    ax.set_aspect("equal")
    ax.set_anchor("N")


# --------------------------------------------------------------------------
# Panel b  --  one field, resolution degrading left to right
# --------------------------------------------------------------------------
def coarse_slab(idx, b):
    """Majority species per block, and how dominant that majority is."""
    ny, nx = idx.shape
    ny2, nx2 = ny // b, nx // b
    maj = np.zeros((ny2, nx2), dtype=int)
    dom = np.zeros((ny2, nx2))
    for r in range(ny2):
        for c in range(nx2):
            blk = idx[r*b:(r+1)*b, c*b:(c+1)*b].ravel()
            v, ct = np.unique(blk, return_counts=True)
            j = int(np.argmax(ct))
            maj[r, c], dom[r, c] = v[j], ct[j] / blk.size
    return maj, dom


def panel_b(ax, region):
    """The same region, five times, each tile at its own resolution."""
    ax.set_axis_off()
    U = 0.235                                  # drawn size of one fine site
    S = REGION * U                             # tile side, identical for all tiles
    seam = 0.30
    x = 0.0
    lefts = []

    for b in BLOCKS:
        maj, _ = coarse_slab(region, b)
        cells(ax, maj, x, 0.0, U * b, ec="white",
              lw=0.0 if b == 1 else 0.30)
        hairline(ax, x, 0, S, S, lw=0.5)
        lefts.append(x)
        x += S + seam

    W = x - seam

    # ---- the l axis
    ya = -1.85
    ax.add_patch(FancyArrowPatch((0, ya), (W + 1.3, ya), arrowstyle="-|>",
                                 mutation_scale=7, lw=0.7, color=RULE))
    for xl in lefts:
        ax.plot([xl + S / 2, xl + S / 2], [ya - 0.25, ya + 0.25],
                lw=0.6, color=RULE)
    ax.text(W + 1.7, ya, r"$\ell$", fontsize=8, ha="left", va="center")

    # ---- the crossover, on the seam after tile LN_AFTER
    xn = lefts[LN_AFTER] - seam / 2
    ax.plot([xn, xn], [ya - 0.45, S + 1.25], lw=0.9, color=INK,
            ls=(0, (2.8, 1.7)), zorder=11)
    ax.text(xn, S + 1.6, r"$\ell_N$", fontsize=8.5, ha="center", va="bottom",
            fontweight="bold")

    # ---- what the two sides mean
    ax.text(xn / 2, ya - 1.15, "species-specific structure\nis resolved",
            fontsize=6.4, ha="center", va="top", color=MUTED, linespacing=1.5)
    ax.text((xn + W) / 2, ya - 1.15, "a neutral description\nis sufficient",
            fontsize=6.4, ha="center", va="top", color=MUTED, linespacing=1.5)

    ax.text(W / 2, S + 4.75, "the same patch of forest, seen through a coarser lens",
            fontsize=6.6, ha="center", va="bottom", color=MUTED, style="italic")

    ax.set_xlim(-0.9, W + 4.2)
    ax.set_ylim(ya - 5.2, S + 6.2)
    ax.set_aspect("equal")
    ax.set_anchor("N")


def main():
    rng = np.random.default_rng(SEED)
    print("simulating ...")
    lat = simulate(FX, FY, NU, KERNEL, SWEEPS, rng)
    idx = recolour(lat, N_SHOWN)
    print(f"  top {N_SHOWN} species cover {(idx >= 0).mean():.0%} of the field")

    # panel a: the most species-rich 9x9 window whose centre is a named species
    best, rc = -1.0, (0, 0)
    for r in range(FY - 9):
        for c in range(FX - 9):
            w = idx[r:r + 9, c:c + 9]
            if w[4, 4] < 0:
                continue
            _, ct = np.unique(w, return_counts=True)
            pr = ct / ct.sum()
            hh = -(pr * np.log(pr)).sum()
            if hh > best:
                best, rc = hh, (r, c)
    patch = idx[rc[0]:rc[0] + 9, rc[1]:rc[1] + 9].copy()
    focal = (4, 4)
    donor, far = None, -1.0
    for dr in range(9):
        for dc in range(9):
            d2 = (dr - 4)**2 + (dc - 4)**2
            if max(abs(dr - 4), abs(dc - 4)) <= KERNEL and patch[dr, dc] >= 0 \
               and patch[dr, dc] != patch[focal] and d2 > far:
                donor, far = (dr, dc), d2
    print(f"  panel a window {rc}, donor {donor}")

    # the single region panel b shows over and over: the most species-rich one
    bestr, r_rc = -1.0, (0, 0)
    for r in range(0, FY - REGION + 1, 4):
        for c in range(0, FX - REGION + 1, 4):
            wnd = idx[r:r + REGION, c:c + REGION]
            _, ct = np.unique(wnd, return_counts=True)
            pr = ct / ct.sum()
            hh = -(pr * np.log(pr)).sum()
            if hh > bestr:
                bestr, r_rc = hh, (r, c)
    region = idx[r_rc[0]:r_rc[0] + REGION, r_rc[1]:r_rc[1] + REGION].copy()
    print(f"  panel b region {r_rc}, coverage {(region >= 0).mean():.0%}")

    fig = plt.figure(figsize=(3.6, 2.95))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.30, 1.45], hspace=0.12)
    a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    panel_a(a1, patch, focal, donor)
    panel_b(a2, region)
    for ax, L in ((a1, "a"), (a2, "b")):
        ax.text(-0.035, 1.0, L, transform=ax.transAxes, fontsize=8.5,
                fontweight="bold", ha="right", va="top")
    for e in ("pdf", "svg", "png"):
        fig.savefig(f"figures/fig2_ab.{e}", dpi=600, transparent=True)
    plt.close(fig)
    print("wrote figures/fig2_ab (.pdf .svg .png)")


if __name__ == "__main__":
    main()
