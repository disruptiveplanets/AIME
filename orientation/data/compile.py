import numpy as np
import matplotlib.pyplot as plt
import corner
import os

PULL = False
ORI_NAME = "ori-e-5"
SIGMA_DISPLAY = 3
N_BINS = 30
TRUES = np.array([0.39269908169, 0.05200629, -0.2021978, 0, 0, 0, 0, 0, 0, 0])

plt.style.use('jcap')

param_names = np.array(["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"])

if PULL:
    os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/orientation/fit_resolved/{ORI_NAME}-0-samples.npy .")
    os.system("scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/fit_resolved/scan-cadence/cad-00-0-samples.npy .")

with open("cad-00-0-samples.npy", 'rb') as f:
    vel_moments = np.load(f).reshape(-1, 10)

with open(f"{ORI_NAME}-0-samples.npy", 'rb') as f:
    o_moments = np.load(f).reshape(-1, 10)

num_select = min(vel_moments.shape[0], o_moments.shape[0])
vel_moments = vel_moments[np.random.choice(np.arange(0, vel_moments.shape[0]), num_select, replace=False)]
o_moments = o_moments[np.random.choice(np.arange(0, o_moments.shape[0]), num_select, replace=False)]

vel_resids = vel_moments - np.mean(vel_moments, axis=0)
o_resids = o_moments - np.mean(o_moments, axis=0)
ratio = np.mean(np.std(vel_resids, axis=0)) / np.mean(np.std(o_resids, axis=0))
print(f"Ratio {ratio}")
vel_resids /= ratio
vel_moments = vel_resids + np.mean(vel_moments, axis=0)

## Save new moments

for name, moments in zip(["scale-vel", "scale-ori"], [vel_moments, o_moments]):
    new_moments = np.zeros_like(moments)
    # Rotate so that K22 -> -K22 
    new_moments[:,0] = moments[:,0] # gamma_0
    new_moments[:,1] = -moments[:,1] # K22
    new_moments[:,2] = moments[:,2] # K20
    (new_moments[:,3], new_moments[:,4]) = (moments[:,4], -moments[:,3]) # K33
    (new_moments[:,5], new_moments[:,6]) = (-moments[:,5], -moments[:,6]) # K32
    (new_moments[:,7], new_moments[:,8]) = (-moments[:,8], moments[:,7]) # K31
    new_moments[:,9] = moments[:,9] # K30
    print(np.mean(new_moments, axis=0))

    np.save(f"{name}-0-samples.npy", new_moments.reshape(-1, 32, 10))

def make_corner():
    param_ranges = []
    for param_index in range(vel_resids.shape[1]):
        width = max(np.std(vel_resids[:, param_index]), np.std(o_resids[:, param_index]))
        param_ranges.append((-SIGMA_DISPLAY * width, SIGMA_DISPLAY * width))

    fig = corner.corner(vel_resids, labels=[f"$\Delta {f}$" for f in param_names], color="C0", range=param_ranges, bins=N_BINS)
    corner.corner(o_resids, fig=fig, color="C1", range=param_ranges, bins=N_BINS)
    plt.savefig("corner-compare.png")
    plt.savefig("corner-compare.pdf")

def make_unc():
    OFFSET = 0
    WIDTH = 0.7
    fig = plt.figure(figsize=(9, 3.5))
    ax_2 = plt.subplot2grid((1, 10+OFFSET), (0, 0), colspan=3)
    ax_3 = plt.subplot2grid((1, 10+OFFSET), (0, 3+OFFSET), colspan=7)
    
    vel_means = np.mean(vel_moments, axis=0)
    o_means = np.mean(o_moments, axis=0)

    mask_2 = [0, 1, 2]
    mask_3 = [3, 4, 5, 6, 7, 8, 9]
    handles = {}
    for (mask, ax) in zip([mask_2, mask_3], [ax_2, ax_3]):
        xs = np.arange(1, len(mask) + 1)
        parts1 = ax.violinplot(vel_moments[:,mask] - TRUES[mask], showextrema=False, widths=WIDTH)
        parts2 = ax.violinplot(vel_moments[:,mask] - TRUES[mask], showextrema=False, widths=WIDTH)
        handles['vel'] = ax.scatter(xs, vel_means[mask] - TRUES[mask], marker='_', c="C0", s=48)
        for pc in parts1['bodies']:
            pc.set_facecolor('C0')
        for pc in parts2['bodies']:
            pc.set_facecolor('none')
            pc.set_edgecolor('C0')
            pc.set_alpha(1.0)

        parts1 = ax.violinplot(o_moments[:,mask] - TRUES[mask], showextrema=False, widths=WIDTH)
        parts2 = ax.violinplot(o_moments[:,mask] - TRUES[mask], showextrema=False, widths=WIDTH)
        handles['o'] = ax.scatter(xs, o_means[mask] - TRUES[mask], marker='_', c="C1", s=48)
        for pc in parts1['bodies']:
            pc.set_facecolor('C1')
        for pc in parts2['bodies']:
            pc.set_facecolor('none')
            pc.set_edgecolor('C1')
            pc.set_alpha(1.0)

        ax.set_xticks(xs)
        ax.set_xticklabels([f"${f}$" for f in param_names[mask]])

    ax_2.set_ylabel("PPDs (centred on true values)")
    ax_2.set_ylim(-2.5e-5, 2.5e-5)
    ax_3.set_ylim(-0.3, 0.3)

    fig.legend([handles['vel'], handles['o']], ["Angular velocity data", "Orientation data"], loc='upper center', ncol=2)
    fig.tight_layout()
    fig.savefig("unc-compare.png")
    fig.savefig("unc-compare.pdf")

if __name__ == "__main__":
    make_unc()
    # plt.show()