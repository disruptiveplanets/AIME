import os
import numpy as np
import matplotlib.pyplot as plt

THRESHOLDS = ((1e-4, False),(1e-3, False),)
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
PULL = True
DIRECTORY = "lumpy"
USE_MEANS = True

if PULL:
    for name in INCREASING.keys():
        os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/{DIRECTORY}/{name}.npy {DIRECTORY}/")


print('Thresholds\t\t'+'\t\t'.join([str(t[0]) for t in THRESHOLDS]))

for fname in os.listdir(DIRECTORY):
    if fname.endswith(".npy"):
        if fname[:-4] not in INCREASING:
            continue
        with open(DIRECTORY + "/" + fname, 'rb') as f:
            uncs = np.load(f)
        if len(fname[:-4]) < 8:
            print(fname[:-4], end='\t\t')
        else:
            print(fname[:-4], end='\t')

        if USE_MEANS:
            use_uncs = uncs[:, 5] # means
        else:
            use_uncs = uncs[:, 2] # medians
        
        plt.figure()
        xs = np.arange(len(uncs))
        plt.plot(xs,  uncs[:,5], label="mean", color="k")
        plt.plot(xs,  uncs[:,5] + uncs[:,6], label="mean+stdev", color="k", linewidth=1, linestyle='dashed')
        #plt.plot(xs,  uncs[:,5] - uncs[:,6], color="k", linewidth=1, linestyle='dashed')
        plt.fill_between(xs, uncs[:,1], uncs[:,3], label="1 sigma", alpha=0.3, color="C0")
        plt.fill_between(xs, uncs[:,0], uncs[:,4], label="2 sigma", alpha=0.3, color="C0")
        plt.plot(xs, uncs[:,2], label="median", color="C0",)
        plt.plot(xs, uncs[:,7], label="max/min", color="C1", linestyle='dotted')
        plt.plot(xs, uncs[:,8], color="C1", linestyle='dotted')

        plt.legend()
        plt.yscale('log')
        plt.title(fname[:-4])

        for threshold, wait in THRESHOLDS:
            plt.axhline(y=threshold, c='r')

            if INCREASING[fname[:-4]]:
                where_more_than_one = np.where(use_uncs>threshold)[0]
                if len(where_more_than_one) > 0:
                    thresh_index = where_more_than_one[0]
                else:
                    print("None", end='\t\t\t\t')
                    continue
                fraction = (threshold - use_uncs[thresh_index-1]) / (use_uncs[thresh_index] - use_uncs[thresh_index-1])
                thresh = thresh_index - 1 + fraction
                
                plt.axvline(x=thresh, c='k')
                print(f"inc.\t{thresh}", end='\t')

            else:
                where_less_than_one = np.where(use_uncs<threshold)[0]
                if len(where_less_than_one) > 0:
                    if wait:
                        # Cut out all the indices that happen before the "last fall": the place where the uncertainty first comes down.
                        last_fall = np.where(use_uncs>threshold)[0]
                        if len(last_fall) > 0:
                            last_fall = last_fall[-1]
                            indices_sub = where_less_than_one - last_fall
                            indices_sub = indices_sub[indices_sub >= 0]
                            if len(indices_sub) > 0:
                                thresh_index = indices_sub[0] + last_fall
                            else:
                                print("None", end='\t\t\t\t')
                                continue
                        else:
                            thresh_index = where_less_than_one[0]
                    else:
                        thresh_index = where_less_than_one[0]
                else:
                    print("None", end='\t\t\t\t')
                    continue
                fraction = (threshold - use_uncs[thresh_index]) / (use_uncs[thresh_index - 1] - use_uncs[thresh_index])
                thresh = thresh_index - fraction
                plt.axvline(x=thresh, c='k')
                print(f"dec.\t{thresh}", end='\t')
        print()
        plt.xlabel("Index")
        plt.ylabel("Uncertainty")
        plt.savefig(f"{fname[:-4]}.png")
plt.show()