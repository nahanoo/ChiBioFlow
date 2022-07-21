import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from curveball import models
from curveball import baranyi_roberts_model
from sympy import symbols, solve, Eq, exp
from model_class import Chemostat

c1 = Chemostat('data/overnight_06_15/M0/2021-04-16_23_37_40_M0_data.csv',0.08)
c2 = Chemostat('data/overnight_06_15/M0/2021-04-16_23_37_40_M0_data.csv',0.01)
chain = [c1]

dilution = c1.get_dilution()
for counter,c in enumerate(chain):
    t0 = c.get_x(c.N0)
    t1 = t0 + 3600
    N1 = c.Nt(3600,c.N0)
    c.N0 = c.Nt(3600,c.N0) * 20 / (20 + 10)
    print(c.N0)
    #c.plot_curve(t0,t1)

