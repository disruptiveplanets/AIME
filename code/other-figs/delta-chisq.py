import numpy as np

TRUE_CENTER_OF_MASS = [0, 15.839935068775791, 0]

TRUE_MOMENTS = {
    "sph-3": [
        
    ],
    "move-1.5": [
        -0.05195945615840046, -0.20085472167788565,
        0, -0.0003845813086399978,
        0, 0,
        0, -0.0019368958923091067,
    ],
} # Uncorrected

UNIFORM_MOMENTS = {
    "sph-3": [
        -0.05200629, -0.2021978,
        0, 0,
        0, 0,
        0, 0,
        0
    ],
    "move-1.5": [
       
    ], # Smoothed
} # Uncorrected

SAMPLE_NAMES = {
    "sph-3": "../density/samples/den-core-sph-3-0-samples.npy",
    "move-1.5": "../density/samples/den-core-move-1.5-0-samples.npy",
}
BULK_AMS = {
    "sph-3",
    "move-1.5": ,
}
SURFACE_AMS = {
    "sph-3": 1000,
    "move-1.5": ,
}
N_MOMENTS = 9


def correct_moments(moments, bulk, surf):
    ls = [2, 2, 3, 3, 3, 3, 3, 3, 3]
    results = []
    for moment, l in zip(moments, ls):
        results.append(moment * (bulk / surf)** (l - 2))
    return results

def chisq(means, true, cov):
    pass

def calc_dsq(name):
    print(name)
    samples = np.load("../density/samples/den-core-move-1.5-0-samples.npy").reshape(-1, 10)[:, 1:]
    print(samples.shape)

    corrected_samples = np.array([correct_moments(s) for s in samples])

    means = np.mean(corrected_samples, axis=0)
    print(means)
    cov = np.cov()
    print(cov.shape)

    true_unif = correct_moments(UNIFORM_MOMENTS[name])
    true_true = correct_moments(TRUE_MOMENTS[name])

    unif_chisq = chisq(means, true_unif, cov)
    true_chisq = chisq(means, true_true, cov)

    print(f"True chisq\t {true_chisq} ({true_chisq / N_MOMENTS})")
    print(f"Unif chisq\t {unif_chisq} ({unif_chisq / N_MOMENTS})")
    print(f"Delta chisq\t {unif_chisq - true_chisq} ({(unif_chisq - true_chisq) / N_MOMENTS})")

if __name__ == "__main__":
    calc_dsq("move-1.5")
    calc_dsq("sph-3")