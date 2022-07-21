from ctypes.wintypes import tagRECT
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
from curveball import models
from curveball import baranyi_roberts_model
from sympy import symbols, solve, Eq, exp, log

df = pd.read_csv(
    'data/overnight_06_15/M0/2021-04-16_23_37_40_M0_data.csv')[['exp_time', 'od_measured']]
df.columns = ['Time', 'OD']
models = models.get_models(baranyi_roberts_model)
logistic = models[1]()
fit = logistic.guess(data=df['OD'], t=df['Time'])
params = {}
for p in fit.keys():
    params[p] = fit[p].value

K = symbols('K')
r = symbols('r')
t = symbols('t')
N0 = symbols('N0')
N = symbols('N')

def Nt(params,t):
    return params['K']/(1 - (1 - params['K']/params['y0']) * exp(-params['r']*t))


def model():
    return K/(1 - (1 - K/N0) * exp(-r*t))


def plot_fit():
    xs = np.linspace(0, 3*3600, 60)
    ys = [Nt(params, x) for x in xs]
    plt.plot(xs, ys)
    plt.show()

def get_dilution():
    target = Nt(params,3600)
    volume = 20 - params['y0'] * 20 /target
    print(volume)

x = []
y = []
x0 = 0

for n in range(10):
    xs = np.linspace(x0,x0 + 3600,60)
    ys = [Nt(params,x) for x in np.linspace(0,3600,60)]
    x0 = xs[-1] + 60
    x = np.concatenate([x,xs])
    y = np.concatenate([y,ys])
plt.plot(x,y)
plt.show()

"""def get_dilution(interval):
    target = Nt(params,3600)
    fN = Eq(N, model())
    fN0 = solve(fN, N0)[0]
    a = fN0.subs({K:params['K'],r:params['r'],N:target,t:interval})"""

