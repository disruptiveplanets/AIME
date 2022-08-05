import numpy as np
import quaternion

UNIFORM_ANGLE = 0

def randomize(model, y, sigma):
    if model == UNIFORM_ANGLE:
        return randomize_uniform_angle(y, sigma)
    
    else:
        raise Exception(f"Model {model} not implemented")

def log_likelihood(model, y, true, length, sigma):
    if model == UNIFORM_ANGLE:
        return like_uniform_angle(y, true, length, sigma)
    else:
        raise Exception(f"Model {model} not implemented")


def randomize_uniform_angle(y, sigma):
    axes = np.random.randn(len(y) * 3).reshape(3, -1)
    axes /= np.linalg.norm(axes, axis=0)
    axes *= np.random.randn(len(y)) * sigma
    deltas = quaternion.from_rotation_vector(axes.transpose())
    return y * deltas

def like_uniform_angle(y, true, length, sigma):
    # Assume y and true have dtype quaternion.
    s = 0
    for (q, q_star) in zip(y[:length], true[:length]):
        s += np.arccos(min(1, abs((q * q_star.conj()).real)))**2
    return -2 / (sigma*sigma) * s

if __name__ == "__main__":
    y = quaternion.as_quat_array(np.random.randn(8))
    for i in range(len(y)):
        if y[i].real < 0:
            y[i] *= -1
        y[i] /= np.sqrt(y[i].norm())

    sigmas = np.linspace(0, 1, 51)
    means = np.zeros_like(sigmas)
    std = np.zeros_like(sigmas)
    for (i, sigma) in enumerate(sigmas):
        likes = np.zeros(1000)
        for trial in range(len(likes)): 
            new_y = randomize_uniform_angle(y, sigma)
            likes[trial] = like_uniform_angle(new_y, y, 0, sigma)
        means[i] = np.mean(likes)
        std[i] = np.std(likes)

    import matplotlib.pyplot as plt
    plt.plot(sigmas, means)
    plt.fill_between(sigmas, means + std, means - std, alpha=0.3)
    plt.xlabel("Sigma")
    plt.ylabel("Like")
    plt.show()