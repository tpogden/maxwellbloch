""" Plot script. """

import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from maxwellbloch import mb_solve

# Args ------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Plotting script.")
parser.add_argument('--pdf', default=False, action="store_true",
    help="Output PDF file.")
parser.add_argument('--png', default=False, action="store_true",
    help="Output PNG file.")
parser.add_argument("--show", default=False, action="store_true",
    help="Show plot in GUI.")
parser.add_argument("--tex", default=False, action="store_true",
    help="Render with TeX.")
opts, unknown = parser.parse_known_args()
print('opts:', opts)

# Data ------------------------------------------------------------------------

MBS_JSON = 'src/mbs-vee-sech-sqrt1-sqrt3.json' 

mbs = mb_solve.MBSolve().from_json(MBS_JSON)
mbs.mbsolve(recalc=False)

# Plot ------------------------------------------------------------------------

sns.set_style('darkgrid')

fig = plt.figure(1, figsize=(10, 5))
ax = fig.add_subplot(211)
cmap_range = np.linspace(0.0, 1.0, 11)
cf = ax.contourf(mbs.tlist, mbs.zlist,
    np.abs(mbs.Omegas_zt[0]/(2*np.pi)),
     cmap_range, 
    cmap=plt.cm.Blues)
# ax.set_title('Rabi Frequency ($\Gamma / 2\pi $)')
ax.set_ylabel('Distance ($L$)')
for y in [0.0, 1.0]:
    ax.axhline(y, c='grey', lw=1.0, ls='dotted')
plt.colorbar(cf)
ax = fig.add_subplot(212)
cmap_range = np.linspace(0.0, 1.0, 11)
cf = ax.contourf(mbs.tlist, mbs.zlist,
    np.abs(mbs.Omegas_zt[1]/(2*np.pi)),
    cmap_range, 
    cmap=plt.cm.Greens)
# ax.set_title('Rabi Frequency ($\Gamma / 2\pi $)')
ax.set_xlabel('Time ($1/\Gamma$)')
ax.set_ylabel('Distance ($L$)')
for y in [0.0, 1.0]:
    ax.axhline(y, c='grey', lw=1.0, ls='dotted')
plt.colorbar(cf)

plt.tight_layout()

# Output ----------------------------------------------------------------------

path_no_ext = os.path.splitext(sys.argv[0])[0]
if opts.pdf:
    plt.savefig(f'{path_no_ext}.pdf')
    print(f'Saved {path_no_ext}.png.')
if opts.png:
    plt.savefig(f'{path_no_ext}.png')
    print(f'Saved {path_no_ext}.png.')
if opts.show:
    plt.show()

