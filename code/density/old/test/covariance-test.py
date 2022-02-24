import corner
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import ortho_group

N_PARAMS = 10
N_TRIALS = 10000
ms = np.array([0,1,0,-1,2,1,0,-1,-2])

ortho = ortho_group.rvs(dim=10)
sigmas = np.random.random(N_PARAMS)
means = np.random.random(N_PARAMS)

params = []
for item in np.random.randn(N_TRIALS*N_PARAMS).reshape(N_TRIALS, N_PARAMS):
    params.append(np.matmul(ortho, np.matmul(item * sigmas+means, ortho.transpose())))
params = np.array(params).transpose()
params[0] -= np.mean(params[0])

hlms_local = np.zeros((N_PARAMS-1, N_TRIALS), dtype=np.complex)
hlms_local[0] = params[1]
hlms_local[1] = params[2] + 1j * params[3]
hlms_local[2] = params[4]
hlms_local[3] = -hlms_local[1].conj()
hlms_local[4] = params[5] + 1j * params[6]
hlms_local[5] = params[7] + 1j * params[8]
hlms_local[6] = params[9]
hlms_local[7] = -hlms_local[5].conj()
hlms_local[8] = hlms_local[4].conj()

exp_mesh = np.exp(1j * np.outer(ms, params[0]))
exp_mesh = 1 - 1j * np.outer(ms, params[0]) # Linear order
#exp_mesh = 1j * np.outer(ms, params[0])
hlms_global = hlms_local * exp_mesh

real_cov = np.cov(hlms_global)

def cov(a,b):
    return np.cov(a,b)[0,1]

i = 4
j = 3
new_cov = np.mean(exp_mesh[i]) * np.mean(exp_mesh[j]).conj() * cov(hlms_local[i], hlms_local[j]) + \
    np.mean(exp_mesh[i]) * np.mean(hlms_local[j]).conj() * cov(hlms_local[i], exp_mesh[j]) + \
    np.mean(hlms_local[i]) * np.mean(exp_mesh[j]).conj() * cov(exp_mesh[i], hlms_local[j]) + \
    np.mean(hlms_local[i]) * np.mean(hlms_local[j]).conj() * cov(exp_mesh[i], exp_mesh[j]) + \
    cov(exp_mesh[i], exp_mesh[j]) * cov(hlms_local[i], hlms_local[j]) + \
    cov(exp_mesh[i], hlms_local[j]) * cov(hlms_local[i], exp_mesh[j])

new_cov = (cov(hlms_local[i], hlms_local[j]) + np.mean(hlms_local[i]) * np.mean(hlms_local[j]).conj()) * \
    (np.var(params[0]) * ms[i] * ms[j] + 1) - 2 * np.mean(hlms_local[i]) * np.mean(hlms_local[j]).conj()+\
    (cov(hlms_local[i], params[0]) * ms[j] - 1j * np.mean(hlms_local[i])) * \
    (cov(params[0], hlms_local[j]) * ms[i] + 1j * np.mean(hlms_local[j]).conj())
print(new_cov, real_cov[i,j])