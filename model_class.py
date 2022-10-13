import numpy as np
from scipy.integrate import odeint


class Chemostat():
    def __init__(self, params, name):
        # Model paramters are stored in chemostat class
        self.K = params['K']
        self.K_init = self.K
        self.r = params['r']
        self.N = params['N']

        # Modelling values
        self.xs = np.ndarray(0)
        self.ys = np.ndarray(0)

        self.name = name

    def model(self, N, t):
        # Simple logistic model
        return self.r * N * (1 - N/self.K)


class Chain():
    def __init__(self, temps):
        # We can create chains by passing temperatures of chemostasts
        # cs = Chain([28.0,33.0,380.43.0]) for example
        params = {28.0: {'r': 0.3,
                         'K': 2,
                         'N': 0.08},
                  33.0: {'r': 0.32,
                         'K': 0.85,
                         'N': 0.06},
                  38.0: {'r': 0.24,
                         'K': 0.26,
                         'N': 0.05},
                  43.0: {'r': 0.1,
                         'K': 0.07,
                         'N': 0.07}}

        self.chain = [Chemostat(params[temp], temp) for temp in temps]
        # We calculate ideal dilution rate
        self.volume = 20
        # Transfer rates are dilutions per hour default is 2
        self.transfer_rate = 2
        self.dilution_rate = 0.31395
        

    def get_dilution(self):
        c1 = self.chain[0]
        # Calculating OD after 1 hour to calculate necessary dilution
        N1 = [e[0] for e in odeint(c1.model, c1.N, [0, 0.5])][-1]
        return ((2 * self.volume * (N1 - c1.N)) / c1.N) / self.volume

    def dilute(self):
        # Function that simulates a dilution row
        self.v_trans = self.dilution_rate * self.volume / self.transfer_rate
        for counter, c in enumerate(self.chain):
            if counter == 0:
                N_in = 0
            else:
                N_in = self.chain[counter - 1].N

            c.N = (N_in * self.v_trans + c.N * self.volume) / \
                (self.v_trans + self.volume)

    def carrying_capacity(self):
        self.v_trans = self.dilution_rate * self.volume / self.transfer_rate
        for counter, c in enumerate(self.chain):
            if counter == 0:
                K_in = c.K_init
            else:
                K_in = self.chain[counter - 1].K

            c.K = (K_in * self.v_trans + c.K * self.volume) / \
                (self.v_trans + self.volume)

    def experiment(self, exp_time):
        if self.transfer_rate != 0:
            # Stores how many transfers are done in one experiment
            intervals = exp_time * self.transfer_rate
            # Time between two dilutions
            interval = 1 / self.transfer_rate
            for i in range(intervals):
                if i != 0:
                    self.dilute()
                    self.carrying_capacity()
                for counter,c in enumerate(self.chain):
                    # Simulated time scale
                    xs = interval * i + np.arange(0, interval, 0.25)
                    xs = np.append(xs, interval+interval*i)
                    c.xs = np.concatenate([c.xs, xs])
                    # Modelled OD values
                    ys = [e[0] for e in odeint(c.model, c.N, xs)]
                    c.ys = np.concatenate([c.ys, ys])
                    # Storing latest OD
                    c.N = c.ys[-1]
                    net_N = ys[-1] - ys[0]
                    c.K = c.K - net_N
                    
        if self.transfer_rate == 0:
            for c in self.chain:
                c.xs = np.arange(0, exp_time, 0.5)
                c.ys = [e[0] for e in odeint(c.model, c.N, c.xs)]
                c.N = c.ys[-1]
