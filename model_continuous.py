from scipy.integrate import odeint
import numpy as np
from math import exp
from matplotlib import pyplot as plt


class Chemostat():
    def __init__(self):
        # Model paramters are stored in chemostat class
        self.N0 = 0.45
        self.Ks = 7957.142857142855
        self.M = 14650
        self.r = 0.24
        self.D = 0.15
        self.t = 400
        self.q = 6.296296296296296e-05
        self.ys = None
        self.xs = np.arange(0, self.t, 0.5)
        self.Rs = None

    def model(self, y, t):
        N = y[0]
        R = y[1]
        dN = self.r * R / (self.Ks + R) * N - self.D * N
        dR = self.D * self.M - self.D * R - N / self.q * self.r * R / (self.Ks + R)
        return [dN, dR]

    def experiment(self):
        self.ys = odeint(self.model, [self.N0, self.M], self.xs)

    def plot(self):
        plt.plot(self.xs, self.ys[:, 0])
        plt.show()


c = Chemostat()
c.experiment()
c.plot()
