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
z = int(opts.zoom)
fps = float(opts.fps)
atoms_alpha = float(opts.atoms_alpha)
show_c_line = opts.c_line
show_peak_line = opts.peak_line



mb_solve_00 = mb_solve.MBSolve().\
                from_json(json_file)

mb_solve_00.mbsolve(recalc=False)

### ANIMATION

tlist_fixed_frame = fixed.t_list(mb_solve_00, speed_of_light)
field_fixed_frame = fixed.rabi_freq(mb_solve_00, 0, speed_of_light,
    part='real', interp_kind='cubic')

t_min = np.min(tlist_fixed_frame)
t_max = np.max(tlist_fixed_frame)

zlist = mb_solve_00.zlist

# Zoom
field_fixed_frame = zoom(field_fixed_frame, z)
tlist_fixed_frame = zoom(tlist_fixed_frame, z)
zlist = zoom(zlist, z)

sns.set_style("darkgrid")
pal = sns.color_palette("deep", 10)


fig = plt.figure(2, figsize=(12, 4))
ax = fig.add_subplot(111)

line, = ax.plot([], [], lw=2, color=pal[2], clip_on=False)
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

ax.set_xlabel('Distance ($L$)')
ax.set_ylabel('Rabi Frequency ($\Gamma / 2\pi $)')

plt.tight_layout()

def init():
    line.set_data([], [])
    t_text.set_text('')
    c_line.set_data([], [])
    peak_line.set_data([], [])
    return line, t_text, peak_line

def animate(i):
    
    x = zlist
    y = field_fixed_frame[:, i]/(2.0*np.pi)

    line.set_data(x, y)

    for coll in (ax.collections):
        ax.collections.remove(coll)

    ax.fill_between(x, 0., y, alpha=0.5, color=pal[2], clip_on=False,
                    interpolate=True)

    c_x = speed_of_light*tlist_fixed_frame[i]
    c_line.set_data([c_x, c_x], c_y)

    peak = np.argmax(y)
    peak_line.set_data([x[peak], x[peak]], peak_y)

    t_text.set_text('t = {:.1f} $/\Gamma$'.format(tlist_fixed_frame[i]))

    return line, t_text, c_line, peak_line

### Simple plot for debugging

# def init():
#     line.set_data([], [])
#     return line

# def animate(i):
#     x = zlist
#     y = np.abs(field_fixed_frame[:, i]) / (2.0 * np.pi)
#     line.set_data(x, y)

#     return line

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
