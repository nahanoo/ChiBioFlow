import numpy as np
from scipy.integrate import odeint


class Specie():
    def __init__(self, r, K, N):
        self.K_init = K
        self.K = self.K_init
        self.r = r
        self.N = N

        self.ys = np.ndarray(0)

    def model(self, N, t):
        # Simple logistic model
        return self.r * N * (1 - N/self.K)


class Chemostat():
    def __init__(self, name):
        # Model paramters are stored in chemostat class
        self.species = {'at': Specie(0.32, 0.7, 0.02),
                        'ct': Specie(0.26, 0.5, 0.02),
                        'ms': Specie(0.35, 0.6, 0.02),
                        'oa': Specie(0.3, 0.8, 0.02)}

        self.name = name
        self.total = np.ndarray(0)
        self.sum_N()
        

    def sum_N(self):
        total = sum([specie.ys for specie in self.species.values()])
        self.total = np.concatenate([self.total,total])
        return

class Chain():
    def __init__(self, chems):
        # We can create chains by passing temperatures of chemostasts
        # cs = Chain([28.0,33.0,380.43.0]) for example

        self.chain = [Chemostat(chem) for chem in range(chems)]
        # We calculate ideal dilution rate
        self.volume = 20
        # Transfer rates are dilutions per hour default is 2
        self.transfer_rate = 2
        self.dilution_rate = 0.31395
        self.xs = np.ndarray(0)


    def get_dilution(self):
        c1 = self.chain[0]
        # Calculating OD after 1 hour to calculate necessary dilution
        N1 = [e[0] for e in odeint(c1.model, c1.N, [0, 0.5])][-1]
        return ((2 * self.volume * (N1 - c1.N)) / c1.N) / self.volume

    def dilute(self):
        # Function that simulates a dilution row
        self.v_trans = self.dilution_rate * self.volume / self.transfer_rate
        for counter, c in enumerate(self.chain):
            for name, specie in c.species.items():
                if counter == 0:
                    N_in = 0
                else:
                    N_in = self.chain[counter - 1].species[name].N
                specie.N = (N_in * self.v_trans + specie.N *
                            self.volume) / (self.v_trans + self.volume)


    def experiment(self, exp_time):
        if self.transfer_rate != 0:
            # Stores how many transfers are done in one experiment
            intervals = exp_time * self.transfer_rate
            # Time between two dilutions
            interval = 1 / self.transfer_rate
            for i in range(intervals):
                xs = interval * i + np.arange(0, interval, 0.25)
                xs = np.append(xs, interval+interval*i)
                self.xs = np.concatenate([self.xs,xs])
                for counter,c in enumerate(self.chain):
                    # Simulated time scale
                    for name,specie in c.species.items():
                        # Modelled OD values
                        ys = [e[0] for e in odeint(specie.model, specie.N, xs)]
                        specie.ys = np.concatenate([specie.ys, ys])
                        # Storing latest OD
                        specie.N = specie.ys[-1]
                        if counter != 0:
                            for name,specie in self.chain[counter].species.items():
                                if False: #specie.K < max_K:
                                    specie.K = self.chain[counter -1].species[name].N + self.chain[counter -1].species[name].K_init
                self.dilute()
            for c in self.chain:
                c.sum_N()
                
                
            

        if self.transfer_rate == 0:
            for c in self.chain:
                c.xs = np.arange(0, exp_time, 0.5)
                c.ys = [e[0] for e in odeint(c.model, c.N, c.xs)]
                c.N = c.ys[-1]