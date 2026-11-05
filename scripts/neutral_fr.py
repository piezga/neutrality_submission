"""F(r) of the neutral model (Fig. 3b), with the SAME (K, nu) as the spectra of Fig. 2.

v2d samples the voter model at stationarity through its coalescent dual on a square NX x NX lattice; the
central window (1000/a) x (500/a) sites is read as the 1000 x 500 m plot, where a = metres per site.
a is a modelling choice and must be the one that defines the metre scale of the spectra
(500 m / lattice side for the SENM domain of Fig. 2), otherwise "k = 3 sites" means a different dispersal
distance in Fig. 2 and Fig. 3.

    python scripts/neutral_fr.py --side 1500 --spacing 1 [--K 14] [--nu 2e-5] [--nrel 3]

Needs neutral_model/v2d compiled with -DDUMP (printconf writes AConf_* = "kx ky species" lines).
Writes data/processed/pairs/cal_K<K>_nu<nu>.pkl  {'rc', 'F': [F(r) per realization], 'S': [species in window]}
"""
import os, sys, glob, pickle, argparse, subprocess, tempfile, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.pairs import Domain, RC, logfit

ap = argparse.ArgumentParser()
ap.add_argument('--side', type=int, required=True, help='lattice side NX = NY (sites)')
ap.add_argument('--spacing', type=float, required=True, help='metres per site')
ap.add_argument('--K', type=int, default=14)
ap.add_argument('--nu', type=float, default=2e-5)
ap.add_argument('--nrel', type=int, default=3)
ap.add_argument('--v2d', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'neutral_model', 'v2d'))
ap.add_argument('--out', default='data/processed/pairs')
args = ap.parse_args()

NX = NY = args.side; a = args.spacing
WX, WY = int(round(1000 / a)), int(round(500 / a))
if WX > NX or WY > NY:
    sys.exit(f'window {WX}x{WY} sites does not fit in the {NX}x{NY} lattice: raise --side or --spacing')
print(f'window {WX}x{WY} sites, spacing {a:.4f} m, kernel half-width {args.K} sites = {args.K * a:.2f} m')

work = tempfile.TemporaryDirectory(); run = work.name
t = time.time()
# v2d exits with status 1 even on success: check the output files instead
subprocess.run([os.path.abspath(args.v2d), str(NX), str(NY), f'{args.nu:e}', str(args.nrel), str(args.K)],
               cwd=run, capture_output=True)
files = sorted(glob.glob(f'{run}/AConf*'))
if len(files) < args.nrel:
    sys.exit('v2d did not write the expected AConf files (compiled without -DDUMP?)')

D = Domain(WX, WY, a); out = {'rc': RC, 'F': [], 'S': [], 'K': args.K, 'nu': args.nu, 'spacing_m': a}
for f in files:
    m = np.fromfile(f, sep=' ', dtype=np.int32).reshape(-1, 3)
    lab = np.empty((NX, NY), np.int32); lab[m[:, 0], m[:, 1]] = m[:, 2]
    lab = lab[(NX - WX) // 2:(NX + WX) // 2, (NY - WY) // 2:(NY + WY) // 2]
    F, S = D.F_labels(lab); out['F'].append(F); out['S'].append(S)

os.makedirs(args.out, exist_ok=True)
pickle.dump(out, open(f'{args.out}/cal_K{args.K}_nu{args.nu:g}.pkl', 'wb')); work.cleanup()
Fm = np.nanmean(out['F'], 0); fit = RC >= 25; b = -np.polyfit(np.log(RC[fit]), Fm[fit], 1)[0]
print(f'species in window {out["S"]}, slope b (r >= 25 m) = {b:.4f}, F(25 m) = {np.interp(25, RC, logfit(RC, Fm)):.4f}, {time.time() - t:.0f} s')
