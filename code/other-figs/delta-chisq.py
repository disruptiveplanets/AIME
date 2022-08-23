import numpy as np

TRUE_CENTER_OF_MASS = [0, 15.839935068775791, 0]

TRUE_MOMENTS = {
    "sph-3": [
        -0.05040959164683728, -0.19599016007639866,
        0, 0,
        0, 0,
        0, 0,
        0
    ],
    "double": [
        -0.05337731641927746, -0.2022901292396581,
        0, 0,
        0, 0,
        0, 0,
        0
    ],
    "move-1.5": [
        -0.05195945615840046, -0.20085472167788565,
        0, -0.0003845813086399978,
        0, 0,
        0, -0.0019368958923091067,
        0
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
    "double": [
        -0.05200629, -0.2021978,
        0, 0,
        0, 0,
        0, 0,
        0
    ],
    "move-1.5": [
       -0.051987377852116975, -0.20218802922577933,
        0, -0.00011727668026358638,
        0, 0,
        0, -0.0005743325909623937,
        0
    ], # Smoothed
} # Uncorrected

SAMPLE_NAMES = {
    "sph-3": "../density/samples/den-core-sph-3-0-samples.npy",
    "double": "../density/samples/den-core-double-0-samples.npy",
    "move-1.5": "../density/samples/den-core-move-1.5-0-samples.npy",
}
BULK_AMS = {
    "sph-3": (922.9234884822591, 1000),
    "double": (970.4652599064898, 1000),
    "move-1.5": (980.8811439828254, 999.8771707704473),
}
SURFACE_AMS = {
    "sph-3": 1000,
    "double": 1000,
    "move-1.5": 1000.1281468600504,
}
N_MOMENTS = 9


def correct_moments(moments, bulk, surf):
    ls = [2, 2, 3, 3, 3, 3, 3, 3, 3]
    results = []
    for moment, l in zip(moments, ls):
        results.append(moment * (bulk / surf)** (l - 2))
    return results

def chisq(means, true, inv_cov):
    return (means - true) @ inv_cov @ (means - true)

def calc_dsq(name):
    print(name)
    samples = np.load(SAMPLE_NAMES[name]).reshape(-1, 10)[:, 1:]
    bulk_true, bulk_unif = BULK_AMS[name]
    surf = SURFACE_AMS[name]

    corrected_samples = np.array([correct_moments(s, bulk_true, surf) for s in samples])

    means = np.mean(corrected_samples, axis=0)
    cov = np.cov(corrected_samples.transpose())
    inv_cov = np.linalg.inv(cov)


    true_unif = correct_moments(np.array(UNIFORM_MOMENTS[name]), bulk_unif, surf)
    true_true = correct_moments(np.array(TRUE_MOMENTS[name]), bulk_true, surf)



    print("Means", means)
    print("Truths", true_true)

    unif_chisq = chisq(means, true_unif, inv_cov)
    true_chisq = chisq(means, true_true, inv_cov)

    print(f"True chisq\t {true_chisq} ({true_chisq / N_MOMENTS})")
    print(f"Unif chisq\t {unif_chisq} ({unif_chisq / N_MOMENTS})")
    print(f"Delta chisq\t {unif_chisq - true_chisq} ({(unif_chisq - true_chisq) / N_MOMENTS})")

if __name__ == "__main__":
    # calc_dsq("sph-3")
    calc_dsq("double")
    # calc_dsq("move-1.5")