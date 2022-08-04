import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use('jcap')

mpl.rcParams["font.size"] = 12

PLOT_ORDER = [
    "scan-perigee",
    "scan-vex",
    "scan-am",
    "scan-period",
    "probe-s-rho",
    "probe-s-theta",
    "scan-cadence",
    "observation-gap",
]

INCREASING = {
    "observation-gap": True,
    "probe-s-rho": True,
    "probe-s-theta": True,
    "scan-am": False,
    "scan-cadence": True,
    "scan-perigee": True,
    "scan-period": False,
    "scan-vex": True,
}
LABELS = {
    "scan-perigee": "$r_p$ (Earth radii)",
    "probe-s-theta": "$\sigma_\\theta$",
    "probe-s-rho": "$\sigma_P / P_\omega$",
    "scan-cadence": "$\Delta t$ (min)",
    "scan-period": "$P_\omega$ (hr)",
    "scan-am": "$a_\\mathcal{A}$ (m)",
    "scan-vex": "$v_\infty$ (km s$^{-1}$)",
    "observation-gap": "$T_\mathrm{gap}$ (hr)",
}

LOGS = {
    "scan-perigee": False,
    "probe-s-theta": True,
    "probe-s-rho": True,
    "scan-cadence": False,
    "scan-period": False,
    "scan-am": True,
    "scan-vex": False,
    "observation-gap": False,
}
TRUES = {
    "scan-perigee": 5,
    "probe-s-theta": 1e-2,
    "probe-s-rho": 1e-7,
    "scan-cadence": 2,
    "scan-period": 9,
    "scan-am": 1000,
    "scan-vex": 6,
    "observation-gap": 0,
}
COLORS = {
    "lumpy": "darkcyan",
    "fe": "olivedrab",
}
EXCLUDE = ["scan-vex"]
PULL = False
DIRECTORIES = ["lumpy", "fe"]

THRESHOLDS = ((1e-4, "dashed", "weak"),(1e-3, "solid", "strong"),)

if PULL:
    for directory in DIRECTORIES:
        for name in INCREASING.keys():
            os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/{directory}/{name}.npy {directory}/")

print('Thresholds\t\t'+'\t\t\t'.join([str(t[0]) for t in THRESHOLDS]))

fig, axs = plt.subplots(ncols=len(INCREASING)//2, nrows=2, figsize=(12,5), sharey=True)
axs = axs.reshape(-1)
handles = {}
for directory in DIRECTORIES:
    for i, name in enumerate(PLOT_ORDER):
        with open(f"{directory}/{name}.npy", 'rb') as f:
            uncs = np.load(f)

        uncs = uncs[~np.any(np.isnan(uncs), axis=1),:]

        if len(name) < 8:
            print(name, end='\t\t')
        else:
            print(name, end='\t')

        use_uncs = uncs[:, 2] # means

        with open(f"../../data/all-fig/{name}-x.npy", 'rb') as f:
            xs = np.load(f)[-len(uncs):]

        handles[directory] = (
            axs[i].plot(xs, uncs[:,2], color=COLORS[directory])[0],
            axs[i].fill_between(xs, uncs[:,0], uncs[:,4], alpha=0.3, color=COLORS[directory])
        )
        axs[i].fill_between(xs, uncs[:,1], uncs[:,3], alpha=0.3, color=COLORS[directory], label='foo')

        if directory == "lumpy":
            for threshold, style, thresh_name in THRESHOLDS:
                handles[thresh_name] = axs[i].axhline(y=threshold, c='r', linewidth=1, linestyle=style)
                if name in EXCLUDE:
                    continue

                locs = np.where(use_uncs > threshold)[0]
                if len(locs) == 0:
                    print(f"\t-\t\t", end="")
                    continue
                if INCREASING[name]:
                    index = np.min(locs)
                else:
                    index = np.max(locs)
                fraction = (threshold - use_uncs[index]) / (use_uncs[index - 1] - use_uncs[index])
                thresh_x = xs[index] - fraction * (xs[index] - xs[index - 1])

                if INCREASING[name]:
                    axs[i].axvspan(xmin=thresh_x, xmax=xs[-1], alpha=0.2, color='k')
                else:
                    axs[i].axvspan(xmin=xs[0], xmax=thresh_x, alpha=0.2, color='k')
                print(f"\t{thresh_x}", end="")
            print()

            axs[i].set_yscale('log')
            if LOGS[name]:
                axs[i].set_xscale('log')
            axs[i].axvline(x=TRUES[name], linewidth=1, c='k', linestyle='dotted')
            axs[i].set_xlabel(LABELS[name])
            axs[i].set_xlim(np.min(xs), np.max(xs))
            if name == "scan-perigee":
                perigee_tick_poses = [5, 15, 30, 45]
                axs[i].set_xticks(perigee_tick_poses)
                axs[i].set_xticklabels([str(p) for p in perigee_tick_poses])
            if i % 4 == 0:
                axs[i].set_ylabel("$\sigma_\\rho / \\rho$")

fig.legend([handles["lumpy"], handles["fe"], handles["weak"], handles["strong"]], ["Lumpy", "Finite element", "Weak cut-off", "Strong cut-off"], loc="upper center", ncol=2)
fig.savefig("all.png", bbox_inches='tight')
fig.savefig("all.pdf", bbox_inches='tight')
plt.show()