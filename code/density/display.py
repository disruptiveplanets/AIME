import numpy as np
import matplotlib.pyplot as plt
from setup import *
from PIL import Image
from multiprocessing import Pool

plt.style.use("jcap")




def make_gif(densities, fname, duration):
    imgs = []
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)
    num_frames = maxes[2] - mins[2]
    for i in range(num_frames):
        fig = make_frame(densities, i + mins[2])
        fig.canvas.draw()
        imgs.append(Image.frombytes('RGB',
            fig.canvas.get_width_height(),fig.canvas.tostring_rgb()))
        plt.close()
    imgs[0].save(fp=fname, format='GIF', append_images=imgs,
                save_all=True, duration= int(duration * 1000 / num_frames), loop=0)

def make_frame(densities, z_index):
    # Save imports until here for supercomputer run
    fig = plt.figure()
    csection = densities[:,:,z_index] / np.nanmean(densities)
    c = plt.pcolormesh(pos_array, pos_array, csection.transpose()   , shading='auto',
        vmin=0,
        vmax=np.nanmax(densities / np.nanmean(densities)),
        cmap='plasma')
    plt.colorbar(c, label="$\\eta / \overline{\\eta}$")
    plt.axis('equal')
    plt.xlabel("$x$ (m)")
    plt.ylabel("$y$ (m)")
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)
    plt.xlim(pos_array[mins[0]], pos_array[maxes[0]])
    plt.ylim(pos_array[mins[1]], pos_array[maxes[1]])
    plt.title("$z$ = {} m".format(str(pos_array[z_index])[:5]))
    fig.tight_layout()
    return fig

def show_cross_section(densities, fname):
    plt.figure()
    csection = densities[:,:,densities.shape[-1]//2] / np.nanmean(densities)
    plt.pcolormesh(pos_array, pos_array, csection.transpose(), shading='auto',
        vmin=0,
        vmax=np.nanmax(densities / np.nanmean(densities)),
        cmap='plasma')
    c = plt.colorbar()
    plt.axis('equal')
    plt.xlabel("$x$ (m)")
    plt.ylabel("$y$ (m)")
    mins = np.min(np.where(~np.isnan(densities)), axis=1)
    maxes = np.max(np.where(~np.isnan(densities)), axis=1)
    plt.xlim(pos_array[mins[0]], pos_array[maxes[0]])
    plt.ylim(pos_array[mins[1]], pos_array[maxes[1]])
    c.set_label("$\\eta / \overline{\\eta}$")
    plt.tight_layout()
    plt.savefig(fname)

def make_figs(barename):
    print(barename.upper())
    with open("data/"+barename+".dat", 'rb') as f:
        densities = np.load(f)
    make_gif(densities, "figs/"+barename+".gif", 2.0)
    print("Gif done")
    show_cross_section(densities, "figs/"+barename+".png")
    show_cross_section(densities, "figs/"+barename+".pdf")
    print("Cross section done")

if __name__ == "__main__":
    #make_figs(TAG+"-harmonic")
    make_figs(TAG+"-surface")
    #make_figs(TAG+"-likelihood")
    #make_figs(TAG+"-ensemble")
    plt.show()
