"""Conspecific pair probability F(r): probability that two stems a distance r apart belong to the
same species (Methods, "Conspecific pair probability").

Stems are placed on a raster of nx x ny sites with spacing `a` metres (BCI: 1000 x 500 sites, a = 1;
neutral-model window: the central window of the lattice, a = lattice spacing). Pair counts at distance
r come from the spatial autocorrelation of each species' raster, by FFT with zero padding (the window
is not wrapped). F(r) = conspecific pairs / all pairs at the same r, which corrects for the edges.
Bins: 1 m up to 10 m, logarithmic above. This module is pure: loading stems is the caller's job.
"""
import numpy as np

EDGES = np.unique(np.round(np.concatenate([[0.5], np.arange(1.5, 10, 1), np.geomspace(10, 500, 34)]), 2))
RC = np.sqrt(EDGES[1:] * EDGES[:-1])


def _by_species(sp):
    """Index arrays, one per distinct label in sp."""
    _, codes = np.unique(sp, return_inverse=True)
    order = np.argsort(codes, kind='stable')
    return np.split(order, np.flatnonzero(np.diff(codes[order])) + 1)


class Domain:
    def __init__(self, nx=1000, ny=500, a=1.0):
        self.nx, self.ny, self.a = nx, ny, a
        self._s = (2 * nx, 2 * ny)
        X = np.minimum(np.arange(2 * nx), 2 * nx - np.arange(2 * nx))[:, None] * a
        Y = np.minimum(np.arange(2 * ny), 2 * ny - np.arange(2 * ny))[None, :] * a
        self._r = np.hypot(X, Y).ravel()

    def raster(self, x, y):
        """Counts per site of stems at positions x, y (metres); stems outside the window are dropped."""
        ix, iy = np.floor(np.asarray(x) / self.a).astype(int), np.floor(np.asarray(y) / self.a).astype(int)
        k = (ix >= 0) & (ix < self.nx) & (iy >= 0) & (iy < self.ny)
        return np.bincount(ix[k] * self.ny + iy[k], minlength=self.nx * self.ny).reshape(self.nx, self.ny).astype(float)

    def pairs(self, R):
        """Ordered pairs of distinct stems per distance bin."""
        A = np.fft.irfft2(np.abs(np.fft.rfft2(R, s=self._s)) ** 2, s=self._s)
        A[0, 0] -= R.sum()
        return np.histogram(self._r, EDGES, weights=A.ravel())[0]

    def F_points(self, x, y, sp):
        """F(r) for stems at (x, y) with species labels sp."""
        x, y, sp = np.asarray(x), np.asarray(y), np.asarray(sp)
        tot, num = np.zeros((self.nx, self.ny)), 0.0
        for idx in _by_species(sp):
            R = self.raster(x[idx], y[idx]); tot += R
            if len(idx) > 1:
                num = num + self.pairs(R)
        with np.errstate(invalid='ignore', divide='ignore'):    # empty distance bins give NaN
            return num / self.pairs(tot)

    def species_pairs(self, x, y, sp, min_n=200):
        """Pair counts of each species with >= min_n stems, and of all stems (for Fig. 3b / SI)."""
        x, y, sp = np.asarray(x), np.asarray(y), np.asarray(sp)
        out, tot = {'n': {}, 'pairs': {}}, np.zeros((self.nx, self.ny))
        for idx in _by_species(sp):
            R = self.raster(x[idx], y[idx]); tot += R
            if len(idx) >= min_n:
                out['pairs'][sp[idx[0]]] = self.pairs(R); out['n'][sp[idx[0]]] = len(idx)
        out['all'] = self.pairs(tot); out['rc'] = RC
        return out

    def F_labels(self, lab):
        """F(r) for a lattice map of species labels (nx, ny), one individual per site; also returns richness."""
        ids, cnt = np.unique(lab, return_counts=True)
        num = sum(self.pairs((lab == s).astype(float)) for s in ids[cnt > 1])
        with np.errstate(invalid='ignore', divide='ignore'):
            return num / self.pairs(np.ones(lab.shape)), len(ids)


def logfit(r, F, lo=25):
    """Logarithmic law a - b ln r fitted for r >= lo, evaluated at r."""
    fit = r >= lo
    b, a = np.polyfit(np.log(r[fit]), F[fit], 1)
    return a + b * np.log(r)
