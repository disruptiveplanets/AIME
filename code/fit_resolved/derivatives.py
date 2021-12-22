import numpy as np

def Gradient(function, step=1.0e-5):
    class g:
        def __init__(self, f, step):
            self.f = f
            self.eps = step

        def __call__(self, x):
            result = np.zeros_like(x)
            for i in range(len(x)):
                epsilon = np.zeros(len(x), dtype=np.float64)
                epsilon[i] = self.eps
                result[i] = (self.f(x+epsilon) - self.f(x-epsilon)) / epsilon[i] / 2
                    # sum over epsilon is its length
            return result

    return g(function, step)

def Hessian(function, step=1.0e-5):
    class h:
        def __init__(self, f, step):
            self.f = f
            self.eps = step

        def __call__(self, x):
            result = np.zeros((len(x), len(x)))
            for i in range(len(x)):
                for j in range(len(x)):
                    epsilona = np.zeros(len(x), dtype=np.float64)
                    epsilonb = np.zeros(len(x), dtype=np.float64)
                    epsilona[i] = self.eps
                    epsilonb[j] = self.eps
                    result[i,j] = (self.f(x+epsilona + epsilonb) + self.f(x - epsilona - epsilonb) \
                        - self.f(x+epsilona - epsilonb) - self.f(x+epsilonb - epsilona)) / epsilona[i] / epsilonb[j] / 4
                    # sum over epsilon is its length
            return (result + result.transpose()) / 2

    return h(function, step)

def minimize(f, start, eps=1.0e-5, grad_squared_stop=1e-5):
    g = Gradient(f, step=eps)
    last_gradient = g(start)
    last_point = start
    point = np.array([s for s in start])
    point[0] += 0.1
    gradient = g(point)
    while True:
        scale = np.abs(np.dot(point - last_point, gradient - last_gradient)) / np.sum((gradient - last_gradient)**2)
        new_point = point - scale * gradient
        new_gradient = g(new_point)
        print(f(new_point))
        last_point = point
        last_gradient = gradient
        point = new_point
        gradient = new_gradient

        if np.sum(gradient**2) / f(new_point) < grad_squared_stop:
            break
    return point
    

if __name__ == "__main__":
    def func(x):
        return np.sum(x**3)
    print(Gradient(func)(
        np.array([1.0, 2.0], dtype=np.float64)
    ))

    print(Hessian(func)(
        np.array([1.0, 2.0], dtype=np.float64)
    ))
