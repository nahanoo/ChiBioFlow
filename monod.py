from scipy.integrate import odeint
import numpy as np
from math import exp
from matplotlib import pyplot as plt


class Chemostat():
    def __init__(self, name, r):
        # Model paramters are stored in chemostat class
        self.name = name
        self.r = 0.28
        self.S = 4
        self.Ks = 1
        self.q = 0.2
        if self.name == 0:
            self.N = 0.1
        else:
            self.N = 0
            self.r = r


class Chain():
    def __init__(self, chems, r):
        self.chain = [Chemostat(c, r) for c in range(chems)]
        self.D = 0.15
        self.S = np.array([c.S for c in self.chain])
        self.rs = np.array([c.r for c in self.chain])
        self.ys = np.array([[c.N for c in self.chain]])
        self.y0 = np.array([])
        self.xs = np.arange(0, 100, 1)
        self.Ks = np.array([c.Ks for c in self.chain])
        self.qs = np.array([c.q for c in self.chain])
        self.Rs = np.array([[c.S for c in self.chain]])
        for i, y in enumerate(self.ys[0]):
            if i == 0:
                self.y0 = np.append(self.y0, y)
            else:
                self.y0 = np.append(self.y0, self.y0[-1] * self.D * self.xs[1])
        self.ys = np.append(self.ys, np.array([self.y0]), axis=0)

    def growth(self, N, t):
        N_in = np.array([])
        for c in self.chain:
            if c.name == 0:
                N_in = np.append(N_in, 0)
            else:
                N_in = np.append(N_in, N[c.name - 1])

        R = odeint(self.resources, self.S, [0, t], args=(N,))[-1]
        return self.rs * R / (self.Ks + R) * N - self.D * N + self.D * N_in

    def resources(self, R, t, N):
        S_in = np.array([])
        for c in self.chain:
            if c.name == 0:
                S_in = np.append(S_in, c.S)
            else:
                S_in = np.append(S_in, R[c.name - 1])
        return self.D * (S_in - R) - self.rs * R / (R + self.Ks) * N / self.qs

    def experiment(self, t):
        ys = odeint(self.growth, self.y0, self.xs[:-2])
        self.ys = np.append(self.ys, ys, axis=0)
        """for x, y in zip(chain.xs, chain.ys):
            rs = np.array(
                [np.array(odeint(self.resources, self.S, [0, x], args=(y,))[-1])])
            self.Rs = np.append(self.Rs, rs, axis=0)"""

    def plot(self):
        plt.plot(self.ys)
        plt.show()



chain = Chain(4, r=0.32)
chain.experiment(100)
chain.plot()
