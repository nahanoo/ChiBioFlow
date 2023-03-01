from scipy.integrate import odeint
import numpy as np
from math import exp
from matplotlib import pyplot as plt
from os.path import join
import pandas as pd


class Chemostat():
    def __init__(self):
        # Model paramters are stored in chemostat class
        # At Glug
        """self.N0 = 0.1
        self.Ks = 9548.571428571428
        self.M = 14650
        self.r = 0.26
        self.D = 0.14
        self.t = 72
        self.q = 6.296296296296296e-05
        self.ys = None
        self.xs = np.arange(0, self.t, 0.5)
        self.Rs = None"""
        """# Ct citrate
        self.N0 = 0.1
        self.Ks = 13129.815784939075
        self.M = 11553.3127304115
        self.r = 0.24
        self.D = 0.1
        self.t = 72
        self.q = 7.387404930499518e-05
        self.ys = None
        self.xs = np.arange(0, self.t, 0.5)"""
        self.Rs = None
        # Oa citrate
        self.N0 = 0.1
        self.Ks = 21924
        self.M = 11553.3127304115
        self.r = 0.38
        self.D = 0.1
        self.t = 72
        self.q = 7.387404930499518e-05/1.3
        self.ys = None
        self.xs = np.arange(0, self.t, 0.5)
        self.Rs = None
    def model_steady_state(self):
        return self.q * (self.M - self.Ks * self.D / (self.r - self.D))

    def model(self, y, t):
        N = y[0]
        R = y[1]
        dN = self.r * R / (self.Ks + R) * N - self.D * N
        dR = self.D * self.M - self.D * R - N / self.q * self.r * R / (self.Ks + R)
        return [dN, dR]

    def experiment(self):
        self.ys = odeint(self.model, [self.N0, self.M], self.xs)

    def plot_model(self):
        plt.plot(self.xs, self.ys[:, 0])
        plt.show()

    def compare_ss(self):
        df = pd.read_csv(join('M2','2022-05-03 15_41_19_M1_data.csv'))[['exp_time','od_measured']]
        df['exp_time'] = df['exp_time'] / 3600
        ss = self.model_steady_state()
        sgraph = [ss for x in df['exp_time']]
        plt.plot(df['exp_time'],df['od_measured'],label='OD ChiBio')
        plt.plot(df['exp_time'],sgraph,label='Modelled steady state')
        plt.xlabel('Time [h]')
        plt.ylabel('OD')
        plt.ylim(0)
        plt.legend()
        return plt
c = Chemostat()
c.experiment()
#fig = c.compare_ss()
c.plot_model()
print(c.model_steady_state())