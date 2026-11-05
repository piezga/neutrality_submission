"""BCI pair-probability data for Fig. 3 and the pairs SI figure (replaces compute_fofr.py).

    python scripts/compute_pairs.py all       -> fofr_bci.pkl        F(r), eight censuses, all stems
    python scripts/compute_pairs.py species   -> perspecies_c{4,8}.pkl   pair counts per species (>= 200 stems)
    python scripts/compute_pairs.py size      -> fofr_size.pkl       F(r) by dbh class, censuses 4 and 8 (needs a dbh column)

Stems are read from <raw>/bci.tree<c>.csv (columns x, y, name[, dbh]), the same files as the rest of the pipeline.
"""
import os, sys, pickle, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.pairs import Domain, RC

ap = argparse.ArgumentParser()
ap.add_argument('what', choices=['all', 'species', 'size'])
ap.add_argument('--raw', default='data/raw'); ap.add_argument('--out', default='data/processed/pairs')
args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
D = Domain(1000, 500, 1.0)


def load(c):
    d = pd.read_csv(f'{args.raw}/bci.tree{c}.csv')
    d = d[(d.x >= 0) & (d.x < 1000) & (d.y >= 0) & (d.y < 500)]
    if 'dbh' in d: d = d[d.dbh > 0]              # all stems with dbh > 0, as in the paper
    return d.reset_index(drop=True)


if args.what == 'all':
    out = {'rc': RC}
    for c in range(1, 9):
        d = load(c); out[c] = D.F_points(d.x.values, d.y.values, d.name.values)
        print('census', c, 'stems', len(d), flush=True)          # compare with the published stem counts
    pickle.dump(out, open(f'{args.out}/fofr_bci.pkl', 'wb'))

elif args.what == 'species':
    for c in (4, 8):
        d = load(c); out = D.species_pairs(d.x.values, d.y.values, d.name.values, min_n=200)
        pickle.dump(out, open(f'{args.out}/perspecies_c{c}.pkl', 'wb'))
        print('census', c, 'species with >= 200 stems:', len(out['n']), flush=True)

else:
    out = {'rc': RC}
    for c in (4, 8):
        d = load(c)
        if 'dbh' not in d: sys.exit('size classes need a dbh column in the raw files')
        for lab, lo, hi in [('10-49mm', 10, 50), ('50-99mm', 50, 100), ('>=100mm', 100, 1e9)]:
            e = d[(d.dbh >= lo) & (d.dbh < hi)]
            out[(c, lab)] = D.F_points(e.x.values, e.y.values, e.name.values); out[(c, lab, 'n')] = len(e)
            print('census', c, lab, len(e), flush=True)
    pickle.dump(out, open(f'{args.out}/fofr_size.pkl', 'wb'))
