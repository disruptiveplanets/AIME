import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv
from PIL import Image
import cmath
from math import factorial, atan2

POINTS_TRY = 50



class Density:
    def __init__(self, grid_size, radius, klms, fn):
        self.fn = fn
        self.maxl = int(np.sqrt(len(klms)))
        assert(len(klms) == self.maxl * self.maxl)
        self.maxl -= 1
        self.volume_element = (2 * radius / grid_size) ** 3
        self.klms = self._process_klms(klms)

        self.mask = np.zeros((grid_size, grid_size, grid_size))
        self.xs = list(np.linspace(-radius, radius, grid_size))
        self.ys = list(np.linspace(-radius, radius, grid_size))
        self.zs = list(np.linspace(-radius, radius, grid_size))

        self._populate_mask()
        self._cut_mask()
        self._make_rlm()

        self.densities = np.copy(self.mask)
        self.densities[self.densities == 0] = np.nan
        self.trials = 0

    def populate(self, trials):
        for trial in range(trials):
            points = self._make_points()
            klmn = []
            for l in range(self.maxl+1):
                for m in range(-l, l+1):
                    line = []
                    for n in range(len(self.klms)):
                        line.append(self._get_klmn(points, l, m, n))
                    klmn.append(np.asarray(line))

            try:
                klmn_inv = np.linalg.inv(klmn)
            except:
                print("Error while inverting matrix")
                continue
            densities = np.matmul(klmn_inv,self.klms)

            for n in range(len(points)):
                self.densities[self.chunk_map==n] += np.real(densities[n])

            self.trials += 1

    def save_gif(self, path, duration=2.0, frames=None):
        if frames is None:
            frames = len(self.zs)
        imgs = []
        for i in range(frames+1):
            fig = self._make_frame(i / frames)
            fig.canvas.draw()
            imgs.append(Image.frombytes('RGB',
                fig.canvas.get_width_height(),fig.canvas.tostring_rgb()))
            plt.close()
        imgs[0].save(fp=path, format='GIF', append_images=imgs,
                 save_all=True, duration= int(duration * 1000 / frames), loop=0)

    def _process_klms(self, klms):
        half_klms = []
        for l in range(self.maxl + 1):
            for m in range(0, 2 * l - 1, 2):
                half_klms.append(complex(klms[l*l + m], klms[l*l + m+1]))
            half_klms.append(complex(klms[(l+1)*(l+1)-1], 0))

        full_klms = []
        for l in range(self.maxl + 1):
            for m in range(-l, l+1, 1):
                index = int(l * (l + 1) / 2) + l - abs(m)
                if m >= 0:
                    full_klms.append(half_klms[index])
                else:
                    full_klms.append(np.conj(half_klms[index]))
        return full_klms

    def _make_rlm(self):
        pass

    def _populate_mask(self):
        for i, x in enumerate(self.xs):
            for j, y in enumerate(self.ys):
                for k, z in enumerate(self.zs):
                    if self.fn(x, y, z):
                        self.mask[i,j,k]=1

    def _cut_mask(self):
        while np.count_nonzero(self.mask[:,:,0]) == 0:
            self.mask = np.delete(self.mask, 0, 2)
            del self.zs[0]
        while np.count_nonzero(self.mask[:,0,:]) == 0:
            self.mask = np.delete(self.mask, 0, 1)
            del self.ys[0]
        while np.count_nonzero(self.mask[0,:,:]) == 0:
            self.mask = np.delete(self.mask, 0, 0)
            del self.xs[0]
        while np.count_nonzero(self.mask[:,:,-1]) == 0:
            self.mask = np.delete(self.mask, -1, 2)
            del self.zs[-1]
        while np.count_nonzero(self.mask[:,-1,:]) == 0:
            self.mask = np.delete(self.mask, -1, 1)
            del self.ys[-1]
        while np.count_nonzero(self.mask[-1,:,:]) == 0:
            self.mask = np.delete(self.mask, -1, 0)
            del self.xs[-1]

    def _make_rlm(self):
        self.rlms = []
        for l in range(0, self.maxl+1):
            for m in range(-l, l+1):
                rlm = np.zeros_like(self.mask, dtype=complex)
                for i, x in enumerate(self.xs):
                    for j, y in enumerate(self.ys):
                        for k, z in enumerate(self.zs):
                            if self.mask[i,j,k] == 0: continue
                            r = np.sqrt(x*x+y*y+z*z)
                            rlm[i, j, k] = lpmv(m, l, z / r) / factorial(l + m)\
                                * r**l * cmath.exp(i * m * atan2(y, x))
                self.rlms.append(rlm)

    def _make_points(self):
        best_points = []
        smallest_volume = 0
        for trial in range(POINTS_TRY):
            points = self._distribute_points()
            map = self._get_chunk_map(points)
            vols = [np.sum(map==i) for i in range(len(points))]
            if np.min(vols) > smallest_volume:
                smallest_volume = np.min(vols)
                best_points = points
        self.chunk_map = self._get_chunk_map(best_points)
        return best_points

    def _get_chunk_map(self, points):
        map = np.zeros_like(self.mask)
        map[self.mask==0] = np.nan
        for i in range(len(self.xs)):
            for j in range(len(self.ys)):
                for k in range(len(self.zs)):
                    if self.mask[i,j,k] == 0: continue
                    map[i, j, k] = self._get_chunk_index(points, i, j, k)
        return map

    def _distribute_points(self):
        points = []
        while len(points) < len(self.klms):
            p = [np.random.randint(0, len(self.xs)),
                 np.random.randint(0, len(self.ys)),
                 np.random.randint(0, len(self.zs))]
            if self.mask[p[0], p[1], p[2]] == 1:
                points.append(p)
        return points

    def _get_klmn(self, points, l, m, n):
        return np.nansum((self.chunk_map==n) * self.rlms[l*l+l+m])\
            * self.volume_element

    def _get_chunk_index(self, points, i, j, k):
        return np.argmin(np.sum((points - np.asarray([i, j, k]))**2, axis=1))

    def _make_frame(self, progress):
        fig = plt.figure()
        slice = self.densities[:,:,
            int(min(progress * len(self.zs), len(self.zs) -1))] / self.trials
        c = plt.pcolormesh(self.xs, self.ys, slice, shading='auto',
            vmin=max(0, np.nanmin(self.densities) / self.trials),
            vmax=np.nanmax(self.densities) / self.trials,
            cmap='plasma')
        plt.colorbar(c)
        plt.axis('equal')
        plt.axis('off')
        plt.title("{}%".format(str(progress*100)[:5]))
        return fig

    def add_densities(self, densities, trials):
        self.trials += trials
        self.densities += densities

    def get_klm(self, i):
        return np.nansum(self.densities * self.rlms[i])\
            / self.trials * self.volume_element, self.klms[i]
