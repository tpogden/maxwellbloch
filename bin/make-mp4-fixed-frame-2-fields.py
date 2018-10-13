#!/usr/bin/env python
# coding: utf-8

""" make-mp4-fixed-frame.py

    Make an MP4 video of a field solution in the fixed frame of reference.
"""

import argparse
import os

# Plot
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn as sns
import numpy as np

from scipy.ndimage.interpolation import zoom

from maxwellbloch import mb_solve, fixed

# Parse filename
parser = argparse.ArgumentParser(description="Takes an MBSolve problem \
    defined in a JSON file and outputs an MP4 video showing the propagation.")
parser.add_argument('-f', '--file', help='Path of input file.', required=True)
parser.add_argument('-s', '--save-path', help='Path to save output MP4.', 
    required=False, default='./')
parser.add_argument('-c', '--speed-of-light',
    help='Speed of Light in the system units.', default=0.1, required=False)
parser.add_argument('-m', '--y-min', help='Minimum of the y-axis.',
    default=0.0, required=False)
parser.add_argument('-y', '--y-max', help='Maximum of the y-axis maximum',
    default=1.0, required=False)
parser.add_argument('--y-max-1', help='Maximum of the second y-axis maximum',
                    default=1.0, required=False)
parser.add_argument('-z', '--zoom', help="To use interpolation on the output \
    data, select the order of interpolation. (e.g. 2, 4). Note this may \
    introduce numerical artefacts.", default=1,
    required=False)
parser.add_argument('-p', '--fps', help='Frames per second', default=30,
    required=False)
parser.add_argument('-a', '--atoms-alpha', help='Atoms alpha', default=0.2,
    required=False)
parser.add_argument('--c-line', action='store_true',
    help='Indicate the speed of light with vertical line', default=False)
parser.add_argument('--peak-line', action='store_true',
    help='Indicate the pulse peak with vertical line', default=False)

opts = parser.parse_args()
print(opts)

json_file = opts.file
save_path = opts.save_path
# TODO: CHECK FILE EXISTS
# TODO: CHECK SAVE PATH EXISTS

speed_of_light = float(opts.speed_of_light)
y_min = float(opts.y_min)
y_max = float(opts.y_max)
y_max_1 = float(opts.y_max_1)
z = int(opts.zoom)
fps = float(opts.fps)
atoms_alpha = float(opts.atoms_alpha)
show_c_line = opts.c_line
show_peak_line = opts.peak_line



mb_solve_00 = mb_solve.MBSolve().\
                from_json(json_file)

mb_solve_00.mbsolve(recalc=False, pbar_chunk_size=2)

### PLOT FIRST

fig = plt.figure(1) #, figsize=(16, 12))

# Probe
ax = fig.add_subplot(211)
cmap_range = np.linspace(0.0, 1.0e-3, 11)
cf = ax.contourf(mb_solve_00.tlist, mb_solve_00.zlist,
                 np.abs(mb_solve_00.Omegas_zt[0] / (2 * np.pi)),
                 cmap_range, cmap=plt.cm.Blues)
plt.colorbar(cf)

# Coupling
ax = fig.add_subplot(212)
cmap_range = np.linspace(0.0, 8.0, 11)
cf = ax.contourf(mb_solve_00.tlist, mb_solve_00.zlist,
                 np.abs(mb_solve_00.Omegas_zt[1] / (2 * np.pi)),
                 cmap_range, cmap=plt.cm.Greens)
plt.colorbar(cf)

plt.tight_layout()


### ANIMATION

tlist_fixed_frame = fixed.t_list(mb_solve_00, speed_of_light)
field_fixed_frame_0 = fixed.rabi_freq(mb_solve_00, 0, speed_of_light,
    part='real', interp_kind='linear')
field_fixed_frame_1 = fixed.rabi_freq(mb_solve_00, 1, speed_of_light,
    part='real', interp_kind='linear')

t_min = np.min(tlist_fixed_frame)
t_max = np.max(tlist_fixed_frame)

zlist = mb_solve_00.zlist

# Zoom
field_fixed_frame_0 = zoom(field_fixed_frame_0, z, order=1)
field_fixed_frame_1 = zoom(field_fixed_frame_1, z, order=1)
tlist_fixed_frame = zoom(tlist_fixed_frame, z, order=1)
zlist = zoom(zlist, z, order=1)

sns.set_style("darkgrid")
pal = sns.color_palette("deep", 10)

fig = plt.figure(2, figsize=(12, 4))
ax = fig.add_subplot(111)
ax1 = ax.twinx()
ax1.grid(False)

line, = ax.plot([], [], lw=2, color=pal[2], clip_on=False)
line_1, = ax1.plot([], [], lw=2, color=pal[3], clip_on=False)

t_text = ax.text(0.90, 0.90, '', transform=ax.transAxes)

c_line, = ax.plot([], [], lw=2, color=pal[1])
c_y = [y_min, y_max]

if show_c_line == False:  #  Hide speed of light indicator line
    c_line.set_visible(False)

peak_line, = ax.plot([], [], lw=2, color=pal[2])
peak_y = [y_min, y_max]

if show_peak_line == False:  # Hide pulse peak indicator line
    peak_line.set_visible(False)

# Atom number density indicator area
ax.axvspan(0.0, 1.0, color=pal[0], alpha=atoms_alpha)

ax.set_xlim((mb_solve_00.z_min, mb_solve_00.z_max))
ax.set_ylim((y_min, y_max))

ax1.set_ylim((y_min, y_max_1))

ax.set_xlabel('Distance ($L$)')
ax.set_ylabel(
    '{0} Rabi Freq ($\Gamma / 2\pi $)'.format(mb_solve_00.atom.fields[0].label), color=pal[2])
ax1.set_ylabel(
    '{0} Rabi Freq ($\Gamma / 2\pi $)'.format(mb_solve_00.atom.fields[1].label), color=pal[3])

plt.tight_layout()

def init():
    line.set_data([], [])
    line_1.set_data([], [])
    t_text.set_text('')
    c_line.set_data([], [])
    peak_line.set_data([], [])
    return line, line_1, t_text, peak_line

def animate(i):
    
    x = zlist
    y_0 = field_fixed_frame_0[:, i]/(2.0*np.pi)
    y_1 = field_fixed_frame_1[:, i]/(2.0*np.pi)

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

    c_x = speed_of_light*tlist_fixed_frame[i]
    c_line.set_data([c_x, c_x], c_y)

    peak = np.argmax(y_0)
    peak_line.set_data([x[peak], x[peak]], peak_y)

    t_text.set_text('t = {:.1f} $/\Gamma$'.format(tlist_fixed_frame[i]))

    return line, line_1, t_text, c_line, peak_line

### Simple plot for debugging

# def init():
#     line.set_data([], [])
#     line_1.set_data([], [])
#     return line, line_1

# def animate(i):
#     x = zlist
#     y_0 = np.abs(field_fixed_frame_0[:, i]) / (2.0 * np.pi)
#     y_1 = np.abs(field_fixed_frame_1[:, i]) / (2.0 * np.pi)
#     line.set_data(x, y_0)
#     line_1.set_data(x, y_1)

#     return line, line_1

# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(tlist_fixed_frame), interval=25,
                               blit=False)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)

print('Saving MP4')
file_stub = os.path.splitext(json_file)[0]

save_path = os.path.join(save_path, file_stub + '.mp4')
print(save_path)
anim.save(save_path, writer=writer)

plt.show()
