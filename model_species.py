import numpy as np
from scipy.integrate import odeint
from math import exp


class Specie():
    # Species class with modelling parameters
    def __init__(self, r, K, N, v, q):
        self.K_init = K
        self.K = self.K_init
        self.r = r
        self.N = N
        self.v = v
        self.q = q
        self.ys = np.ndarray(0)


class Chemostat():
    # Chemostat class
    def __init__(self, name):
        self.species = {'at': Specie(3, 0.27, 0.025, 0.08, 1),
                        'ct': Specie(0.09, 0.01, 0.025, 1, 1),
                        'ms': Specie(0.38, 0.062, 0.025, 1, 1),
                        'oa': Specie(0.28, 0.65, 0.025, 1, 0.1)}

        self.name = name
        # Concentrations of species are summed in total
        self.total = np.ndarray(0)
        self.sum_N()
        # Carrying capacity of chemostat
        self.K = 1.5
        self.K_init = self.K
        self.Ks = [self.K]

        self.dilution_factors = []

    def sum_N(self):
        # Sums concentrations of species
        total = sum([specie.ys for specie in self.species.values()])
        self.total = np.concatenate([self.total, total])
        return


class Chain():
    def __init__(self, chems):
        self.chain = [Chemostat(chem)
                      for chem in range(chems)]
        self.volume = 20
        # Transfer rates are dilutions per hour default is 2
        self.transfer_rate = 2
        self.dilution_rate = 0.31395
        self.xs = np.ndarray(0)
        self.x_dilutions = []

    def model(self, N, t, r, K, v):
        # Simple logistic model
        if v is None:
            return r * N * (1 - N/K)

        if v is not None:
            return (v / (v + exp(-t))) * r * N * (1 - N/K)

    def lag(self, t, q, v):
        return q/(q + exp(-v*t))

    def curved_lag(self, N, t, k, K, v, q):
        a = self.lag(t, q, v)
        return a * k * N * (1 - (N/K)**v)

    def get_dilution_rate(self, specie):
        # Calculates dilution rate for keeping one species steady state
        c1 = self.chain[0]
        # Calculating OD after 1 hour to calculate necessary dilution
        interval = 1 / self.transfer_rate
        s = c1.species[specie]
        xs = [0, interval]
        N1 = [e[0]
              for e in odeint(self.model, s.N, xs, args=(s.r, s.K, s.v))][-1]
        return ((2 * self.volume * (N1 - s.N)) / s.N) / self.volume

    def dilute(self):
        # Function that simulates a dilution row
        self.v_trans = self.dilution_rate * self.volume / self.transfer_rate
        for counter, c in enumerate(self.chain):
            N0 = 0
            N1 = 0
            for name, specie in c.species.items():
                N0 += specie.N
            for name, specie in c.species.items():
                if counter == 0:
                    N_in = 0
                else:
                    N_in = self.chain[counter - 1].species[name].N
                specie.N = (N_in * self.v_trans + specie.N *
                            self.volume) / (self.v_trans + self.volume)
                N1 += specie.N
            c.dilution_factors.append(N0/N1)

    def carrying_capacity(self):
        self.v_trans = self.dilution_rate * self.volume / self.transfer_rate
        for counter, c in enumerate(self.chain):
            if counter == 0:
                K_in = c.K_init
            else:
                K_in = self.chain[counter - 1].K
            c.K = (K_in * self.v_trans + c.K * self.volume) / \
                (self.v_trans + self.volume)
            c.Ks.append(c.K)

    def experiment(self, exp_time):
        if self.transfer_rate != 0:
            # Stores how many transfers are done in one experiment
            intervals = int(exp_time * self.transfer_rate)
            # Time between two dilutions
            interval = 1 / self.transfer_rate
            for i in range(intervals):
                xs = interval * i + np.arange(0, interval, 1/5)
                xs = np.append(xs, interval+interval*i)
                self.xs = np.concatenate([self.xs, xs])
                if i != 0:
                    self.dilute()
                    self.x_dilutions.append(xs[0])
                    self.carrying_capacity()
                for counter, c in enumerate(self.chain):
                    # Simulated time scale
                    net_N = 0
                    for name, specie in c.species.items():
                        # Modelled OD values
                        K = c.K * specie.K
                        ys = [e[0] for e in odeint(
                            self.curved_lag, specie.N, xs, args=(specie.r, K, specie.v, specie.q))]
                        specie.ys = np.concatenate([specie.ys, ys])
                        # Storing latest OD
                        specie.N = specie.ys[-1]
                        net_N += ys[-1] - ys[0]
                    #c.K = c.K - net_N

        if self.transfer_rate == 0:
            self.xs = np.arange(0, exp_time, 1/30)
            for c in self.chain:
                for name, specie in c.species.items():
                    K = c.K * specie.K
                    ys = [e[0] for e in odeint(
                        self.model, specie.N, self.xs, args=(specie.r, K, specie.v))]
                    specie.ys = np.concatenate([specie.ys, ys])
                    # Storing latest OD
                    #specie.N = specie.ys[-1]

        for c in self.chain:
            c.sum_N()
