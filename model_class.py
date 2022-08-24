import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from curveball import models
from curveball import baranyi_roberts_model
from sympy import symbols, solve, Eq, exp


class Chemostat():
    def __init__(self, f, N0):
        self.K = symbols('K')
        self.r = symbols('r')
        self.t = symbols('t')
        self.N0 = N0
        self.N = symbols('N')
        self.fit(f)
        self.N1 = None

        self.xs = np.ndarray(0)
        self.ys = np.ndarray(0)

    def fit(self, f):
        self.df = pd.read_csv(
            f)[['exp_time', 'od_measured']]
        self.df.columns = ['Time', 'OD']
        self.df['Time'] = self.df['Time'] / 60**2
        ms = models.get_models(baranyi_roberts_model)
        logistic = ms[1]()
        fit = logistic.guess(data=self.df['OD'], t=self.df['Time'],param_guess={'y0':self.N0})
        self.params = {}
        for p in fit.keys():
            self.params[p] = fit[p].value


    def Nt(self, t, N0):
        return self.params['K']/(1 - (1 - self.params['K']/N0) * exp(-self.params['r']*t))

    def plot_fit(self):
        xs = self.df['Time']
        ys = [self.Nt(x, self.params['y0']) for x in self.df['Time']]
        ts = self.df['Time']
        ods = self.df['OD']
        plt.plot(xs, ys)
        plt.plot(ts,ods)
        plt.show()

    def plot_curve(self, t0, t1):
        xs = np.linspace(t0, t1, 60)
        ys = [self.Nt(x, self.N0) for x in xs]
        plt.plot(xs, ys)    
        plt.show()

    def get_x(self, N):
        fN = Eq(N, self.model())
        ft = solve(fN, self.t)[0]
        return float(ft.subs({self.K: self.params['K'], self.r: self.params['r']}))

    def get_dilution(self):
        return 20 * (self.Nt(3600,self.N0) - self.N0) / self.N0

    def model(self):
        return self.K/(1 - (1 - self.K/self.N0) * exp(-self.r*self.t))