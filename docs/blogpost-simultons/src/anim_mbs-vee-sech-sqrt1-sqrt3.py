""" Animation script. """

import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn as sns

from scipy.ndimage.interpolation import zoom

from maxwellbloch import mb_solve, fixed

# Args ------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Animation script.")
parser.add_argument('--mp4', default=False, action="store_true",
    help="Output MP4 file.")
parser.add_argument("--show", default=False, action="store_true",
    help="Show plot in GUI.")
parser.add_argument("--tex", default=False, action="store_true",
    help="Render with TeX.")
parser.add_argument('-p', '--fps', help='Frames per second', default=30,
    required=False)
parser.add_argument('-c', '--speed-of-light',
    help='Speed of Light in the system units.', default=0.1, required=False)
parser.add_argument('-z', '--zoom', help="To use interpolation on the output \
    data, select the order of interpolation. (e.g. 2, 4). Note this may \
    introduce numerical artefacts.", default=1,
    required=False)
opts, unknown = parser.parse_known_args()
print('opts:', opts)

# Data ------------------------------------------------------------------------

MBS_JSON = 'src/mbs-vee-sech-sqrt1-sqrt3.json' 

mbs = mb_solve.MBSolve().from_json(MBS_JSON)
mbs.mbsolve(recalc=False)

# Shift to fixed frame-f-reference
tlist_fixed = fixed.t_list(mbs, float(opts.speed_of_light))
rabi_freq_fixed_0 = fixed.rabi_freq(mbs, 0, float(opts.speed_of_light),
    part='real', interp_kind='linear')
rabi_freq_fixed_1 = fixed.rabi_freq(mbs, 1, float(opts.speed_of_light),
    part='real', interp_kind='linear')
t_min = np.min(tlist_fixed)
t_max = np.max(tlist_fixed)
zlist = mbs.zlist

# Interpolate (Zoom)
if zoom != 1:
    rabi_freq_fixed_0 = zoom(rabi_freq_fixed_0, int(opts.zoom), order=1)
    rabi_freq_fixed_1 = zoom(rabi_freq_fixed_1, int(opts.zoom), order=1)
    tlist_fixed = zoom(tlist_fixed, int(opts.zoom), order=1)
    zlist = zoom(zlist, int(opts.zoom), order=1)

# Anim ------------------------------------------------------------------------

Y_MIN = 0.0
Y_MAX = 1.0
Y_MAX_1 = 1.0
ATOMS_ALPHA = 0.2

sns.set_style("darkgrid")
pal = sns.color_palette("deep", 10)

fig = plt.figure(1, figsize=(12, 4))
ax = fig.add_subplot(111)
ax1 = ax.twinx()
ax1.grid(False)

line, = ax.plot([], [], lw=2, color=pal[2], clip_on=False)
line_1, = ax1.plot([], [], lw=2, color=pal[3], clip_on=False)

t_text = ax.text(0.90, 0.90, '', transform=ax.transAxes)

# Atom number density indicator area
ax.axvspan(0.0, 1.0, color=pal[0], alpha=ATOMS_ALPHA)

ax.set_xlim((mbs.z_min, mbs.z_max))
ax.set_ylim((Y_MIN, Y_MAX))
ax1.set_ylim((Y_MIN, Y_MAX_1))
ax.set_xlabel('Distance ($L$)')
ax.set_ylabel(
    '{0} Rabi Freq ($\Gamma / 2\pi $)'.format(mbs.atom.fields[0].label), 
        color=pal[2])
ax1.set_ylabel(
    '{0} Rabi Freq ($\Gamma / 2\pi $)'.format(mbs.atom.fields[1].label),
        color=pal[3])

plt.tight_layout()

def init():
    line.set_data([], [])
    line_1.set_data([], [])
    t_text.set_text('')
    # c_line.set_data([], [])
    # peak_line.set_data([], [])
    return line, line_1, t_text#, peak_line

def animate(i):
    
    x = zlist
    y_0 = rabi_freq_fixed_0[:, i]/(2.0*np.pi)
    y_1 = rabi_freq_fixed_1[:, i]/(2.0*np.pi)

    line.set_data(x, y_0)
    line_1.set_data(x, y_1)

    for coll in (ax.collections):
        ax.collections.remove(coll)

    for coll in (ax1.collections):
        ax1.collections.remove(coll)

    ax.fill_between(x, 0., y_0, alpha=0.5, color=pal[2], clip_on=False,
                    interpolate=True)
    ax1.fill_between(x, 0., y_1, alpha=0.2, color=pal[3], clip_on=False,
                    interpolate=True)

    # c_x = speed_of_light*tlist_fixed_frame[i]
    # c_line.set_data([c_x, c_x], c_y)

    # peak = np.argmax(y_0)
    # peak_line.set_data([x[peak], x[peak]], peak_y)

    t_text.set_text('t = {:.1f} $/\Gamma$'.format(tlist_fixed[i]))

    return line, line_1, t_text#, c_line, peak_line

# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(tlist_fixed), interval=25,
                               blit=False)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=opts.fps, metadata=dict(artist='Me'), bitrate=1800)

# Output ----------------------------------------------------------------------

path_no_ext = os.path.splitext(sys.argv[0])[0]
if opts.mp4:
    anim.save(f'{path_no_ext}.mp4', writer=writer)
    print(f'Saved {path_no_ext}.mp4.')
if opts.show:
    plt.show()
