import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sympy import symbols, solve, Eq, exp
from scipy.integrate import odeint


class Chemostat():
    def __init__(self, params):
        self.K = params['K']
        self.r = params['r']
        self.N = params['N']

        self.xs = np.ndarray(0)
        self.ys = np.ndarray(0)

    def model(self, N, t):
        return self.r * N * (1 - N/self.K)


class Chain():
    def __init__(self):
        self.chain = [Chemostat({'r': 0.456,
                                 'K': 0.8,
                                 'N': 0.1}),
                      Chemostat({'r': 0.456,
                                 'K': 1,
                                 'N': 0}),
                      Chemostat({'r': 0.456,
                                 'K': 1,
                                 'N': 0}),
                      Chemostat({'r': 0.456,
                                 'K': 1,
                                 'N': 0})
        ]
        self.volume = 20
        self.dilution_rate = 0.313
        self.transfer_rate = 1

    

    def dilute(self):
        v_trans = self.dilution_rate * self.volume / self.transfer_rate
        for counter, c in enumerate(self.chain):
            if counter == 0:
                N_in = 0
            else:
                N_in = self.chain[counter - 1].N

            c.N = (N_in * v_trans + c.N * self.volume) / (v_trans + self.volume) 
        
    def experiment(self, exp_time):
        intervals = exp_time * self.transfer_rate
        interval = 1 / self.transfer_rate
        for i in range(intervals):
            for c in self.chain:
                xs = interval * i + np.arange(0,interval,1/3600)
                c.xs = np.concatenate([c.xs,xs])
                ys = [e[0] for e in odeint(c.model,c.N,xs)]
                c.ys = np.concatenate([c.ys, ys])
                c.N = c.ys[-1]
            self.dilute()


k = Chain()
k.experiment(48)
for c in k.chain:
    plt.plot(c.xs,c.ys)

plt.show()
#plt.show()
"""    def dilute(self):
        return 20 * (self.Nt(3600, self.N0) - self.N0) / self.N0"""
