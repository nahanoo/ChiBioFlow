from scipy.integrate import odeint
import numpy as np
from math import exp
from matplotlib import pyplot as plt


class Chemostat():
    def __init__(self, name):
        # Model paramters are stored in chemostat class
        self.name = name
        self.K = 1.5
        self.r = 0.28
        if self.name == 0:
            self.N = 0.1
        else:
            self.N = 0


class Chain():
    def __init__(self, chems):
        self.chain = [Chemostat(c) for c in range(chems)]
        self.D = 0.15
        self.rs = np.array([c.r for c in self.chain])
        self.Ks = np.array([c.K for c in self.chain])
        self.ys = np.array([[c.N for c in self.chain]])
        self.y0 = np.array([])
        for i, y in enumerate(self.ys[0]):
            if i == 0:
                self.y0 = np.append(self.y0, y)
            else:
                self.y0 = np.append(self.y0, self.y0[-1] * self.D)
        self.ys = np.append(self.ys, np.array([self.y0]), axis=0)

    def model(self, y, t):
        N_in = np.array([])
        for c in self.chain:
            if c.name == 0:
                N_in = np.append(N_in, 0)
            else:
                N_in = np.append(N_in, y[c.name - 1])
        return self.rs * (self.Ks - y) / self.Ks * y + self.D * N_in - self.D * y

    def experiment(self, t):
        xs = np.arange(0, t, 1)
        ys = odeint(self.model, self.y0, xs)
        self.ys = np.append(self.ys, ys, axis=0)

    def plot(self):
        plt.plot(self.ys)
        plt.show()


chain = Chain(4)
c0 = chain.chain[0]
chain.experiment(72)
