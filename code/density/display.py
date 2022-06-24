import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#plt.style.use("jcap")

FONT_SIZE = 16

import matplotlib as mpl
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.monospace"] = "Roboto mono"
mpl.rcParams["text.usetex"] = "true"
mpl.rcParams["figure.figsize"] = (4.25, 4)
mpl.rcParams["legend.framealpha"] = 0.5
mpl.rcParams["lines.linewidth"] = 2
mpl.rcParams["lines.markersize"] = 4
mpl.rcParams["font.size"] = FONT_SIZE
mpl.rcParams["legend.fontsize"] = 10

VERY_SMALL = 1
NUM_SLICES = 6
EXPAND_X=1.15 # Scale along x and y so that spheres look circular
AXIS_LIMIT = 1000
EPSILON = 0.00000001
NUM_CONTOURS = 6


def latexify(number):
    string = f"{number:.1e}"
    if 'e' in string:
        index = string.find('e')
        return string[:index] + "\\times 10^{" + str(int(string[index+1:])) + "}"
    return string


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
    if cmap == 'Greys_r': # Uncertainty
        want_min = np.nanmin(densities)
        
    if balance:
        want_max = max(want_max, -want_min)
        want_min = -want_max

    densities = np.clip(densities, want_min+(abs(want_min) * EPSILON), want_max-(abs(want_max) * EPSILON))


    # Does the colorbar have an arrow?
    extend = 'neither'
    if percentile < 98:
        if balance:
            extend = 'both'
        else:
            extend = 'max'

    c = plt.pcolormesh(pos_array, pos_array, csection.transpose(), shading='auto', vmin=want_min, vmax=want_max, cmap=cmap)
    plt.colorbar(c, extend=extend).set_label(axis_name)
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

def make_slices(densities, pos_array, axis_name, cmap, name, klm_error, percentile=99, balance=False):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_axis_off()
    ax.grid(False)

    want_min, want_max = np.nanpercentile(densities, 1), np.nanpercentile(densities, percentile)

    if cmap == 'Greys_r': # Uncertainty
        want_min = np.nanmin(densities)
        want_min = 0

    # Does the colorbar have an arrow?
    extend = 'neither'
    if percentile < 98:
        if balance:
            extend = 'both'
        else:
            extend = 'max'

    if balance:
        want_max = min(want_max, -want_min)
        want_min = -want_max
        if want_max < want_min:
            want_max, want_min = want_min, want_max

    densities = np.clip(densities, want_min+(abs(want_min) * EPSILON), want_max-(abs(want_max) * EPSILON))


    if want_min != want_max:
        levels = np.linspace(want_min, want_max, NUM_CONTOURS + 1)
    else:
        levels = np.linspace(0.99, 1.01, NUM_CONTOURS + 1)

    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)

    for i in np.linspace(mins[2]+1, maxes[2]-1, NUM_SLICES):
        i = int(i)
        z = pos_array[i]
        ax.contourf(pos_array, pos_array, z+densities[:,:,i]*VERY_SMALL,
            zdir='z', levels=z+VERY_SMALL*levels, cmap=cmap)

    ax.view_init(elev=8, azim=45)

    fig2 = plt.figure()
    ax2 = fig2.gca()
    contour_handle = ax2.contourf(pos_array, pos_array, densities[:,:,0], levels=levels, cmap=cmap, extend=extend)

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

    ax.set_title(f"$\\chi^2_r={str(klm_error)[:4]}$", y=0.9, fontdict={'fontsize': FONT_SIZE})

    axins = inset_axes(ax,
                    width="100%",  
                    height="5%",
                    loc='lower center',
                    borderpad=-6
                    )

    from matplotlib.ticker import FuncFormatter
    fmt = lambda x, pos: round(x, 2)
    c = fig.colorbar(contour_handle, ax=axins, extend=extend, orientation="horizontal", format=FuncFormatter(fmt))
    c.set_label(axis_name)


    fig.tight_layout()
    try:
        fig.savefig(name+".pdf")
        fig.savefig(name+".png")
    except Exception:
        pass