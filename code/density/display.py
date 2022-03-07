import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from multiprocessing import Pool

#plt.style.use("jcap")

import matplotlib as mpl
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.monospace"] = "Roboto mono"
mpl.rcParams["text.usetex"] = "true"
mpl.rcParams["figure.figsize"] = (6.5, 4.00)
mpl.rcParams["legend.framealpha"] = 0.5
mpl.rcParams["lines.linewidth"] = 2
mpl.rcParams["lines.markersize"] = 4
mpl.rcParams["font.size"] = 16
mpl.rcParams["legend.fontsize"] = 10

VERY_SMALL = 1
NUM_SLICES = 6
EXPAND_X=1.15 # Scale along x and y so that spheres look circular
AXIS_LIMIT = 1000


def make_gif(densities, pos_array, axis_name, cmap, fname, duration, percentile=99, balance=False):
    imgs = []
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)
    num_frames = maxes[2] - mins[2]
    max_frames = duration * 30
    skip_between = 1
    while len(range(0, num_frames, skip_between)) > max_frames:
        skip_between += 1
    print(f"Num frames: {num_frames}, skipping {skip_between}")
    for i in range(0, num_frames, skip_between):
        fig = make_frame(densities, pos_array, axis_name, cmap, percentile, balance, i + mins[2])
        fig.canvas.draw()
        imgs.append(Image.frombytes('RGB',
            fig.canvas.get_width_height(),fig.canvas.tostring_rgb()))
        plt.close()
    imgs[0].save(fp=fname, format='GIF', append_images=imgs,
                save_all=True, duration= int(duration * 1000 / num_frames), loop=0)

def make_frame(densities, pos_array, axis_name, cmap, percentile, balance, z_index):
    # Save imports until here for supercomputer run
    fig = plt.figure()
    csection = densities[:,:,z_index]
    want_min, want_max = np.nanpercentile(densities, 1), np.nanpercentile(densities, percentile)
    if balance:
        want_max = max(want_max, -want_min)
        want_min = -want_max
    c = plt.pcolormesh(pos_array, pos_array, csection.transpose()   , shading='auto',
        vmin=want_min,
        vmax=want_max,
        cmap=cmap)
    plt.colorbar(c, label=axis_name)
    plt.axis('equal')
    plt.xlabel("$y$ (m)")
    plt.ylabel("$x$ (m)")
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)
    plt.xlim(pos_array[mins[0]], pos_array[maxes[0]])
    plt.ylim(pos_array[mins[1]], pos_array[maxes[1]])
    plt.title("$z$ = {} m".format(str(pos_array[z_index])[:5]))
    fig.tight_layout()
    return fig

def show_cross_section(densities, pos_array, axis_name, cmap, fname, percentile=99, balance=False):
    plt.figure()
    csection = densities[:,:,densities.shape[-1]//2]
    plt.pcolormesh(pos_array, pos_array, csection.transpose(), shading='auto',
        vmin= np.nanpercentile(densities, 1),
        vmax= np.nanpercentile(densities, percentile),
        cmap=cmap)
    c = plt.colorbar()
    plt.axis('equal')
    plt.xlabel("$x$ (m)")
    plt.ylabel("$y$ (m)")
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)
    plt.xlim(pos_array[mins[0]], pos_array[maxes[0]])
    plt.ylim(pos_array[mins[1]], pos_array[maxes[1]])
    c.set_label(axis_name)
    plt.tight_layout()
    plt.savefig(fname)

def make_slices(densities, pos_array, axis_name, cmap, name, percentile=99, balance=False):
    fig = plt.figure(figsize=(8,5))
    ax = fig.gca(projection='3d')
    ax.set_axis_off()
    ax.grid(False)

    want_min, want_max = np.nanpercentile(densities, 1), np.nanpercentile(densities, percentile)
    if balance:
        want_max = max(want_max, -want_min)
        want_min = -want_max

    if np.nanmin(densities) != np.nanmax(densities):
        levels = np.linspace(want_min, want_max, 40)
    else:
        levels = np.linspace(0.99, 1.01, 40)
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)


    for i in np.linspace(mins[2]+1, maxes[2]-1, NUM_SLICES):
        i = int(i)
        z = pos_array[i]
        ax.contourf(pos_array, pos_array, z+densities[:,:,i]*VERY_SMALL,
            zdir='z', levels=z+VERY_SMALL*levels, cmap=cmap)

    ax.view_init(elev=10, azim=45)

    fig2 = plt.figure()
    ax2 = fig2.gca()
    contour_handle = ax2.contourf(pos_array, pos_array, densities[:,:,0], levels=levels, cmap=cmap)

    '''max_radius = max([max(abs(pos_array[mins[i]]), abs(pos_array[maxes[i]])) for i in range(len(densities.shape))])
    ax.set_xlim3d(-max_radius * EXPAND_X / 2, max_radius * EXPAND_X / 2)
    ax.set_ylim3d(-max_radius * EXPAND_X / 2, max_radius * EXPAND_X / 2)
    ax.set_zlim3d(-max_radius / 2, max_radius / 2)'''
    ax.set_xlim3d(-AXIS_LIMIT * EXPAND_X, AXIS_LIMIT * EXPAND_X)
    ax.set_ylim3d(-AXIS_LIMIT * EXPAND_X, AXIS_LIMIT * EXPAND_X)
    ax.set_zlim3d(-AXIS_LIMIT, AXIS_LIMIT)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$z$")

    c = fig.colorbar(contour_handle, ax=ax)
    c.set_label(axis_name)

    fig.tight_layout()
    fig.savefig(name+".pdf")
    fig.savefig(name+".png")

def make_figs(barename):
    print(barename.upper())
    with open("data/"+barename+".dat", 'rb') as f:
        densities = np.load(f)
    make_gif(densities, "figs/"+barename+".gif", 2.0)
    print("Gif done")
    show_cross_section(densities, "figs/"+barename+".png")
    show_cross_section(densities, "figs/"+barename+".pdf")
    print("Cross section done")
    make_slices(densities, "figs/"+barename+"-slice")