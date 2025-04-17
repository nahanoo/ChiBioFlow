from equations import *
import numpy as np
from matplotlib import pyplot as plt

J = eq_J1_1
#solve(J,v1_1)[0].subs({J1_1:0.25,K1_1:10,R1:15})
v_max = 0.416666666666667
M = 15
q = 1 / M
R_star = solve(J, R1)
Ks = np.linspace(1E-2, 10, 1000)
Ds = np.linspace(0.001, 0.25, 50, 100)
K_grid, D_grid = np.meshgrid(Ks, Ds)
r_stars = -D_grid * K_grid / (D_grid - v_max)
plt.contourf(K_grid, D_grid, r_stars, levels=100, cmap="RdYlBu")
plt.colorbar(label="R* [mM]")
plt.contour(K_grid, D_grid, r_stars, levels=[1])
plt.xscale('log')
plt.show()

N_star = solve(Eq(0, eq_dR1_1.rhs), N1)[0]
n_stars = D_grid * q * (K_grid * M + r_stars *
                        (-K_grid + M - r_stars)) / (r_stars * v_max)
plt.contourf(K_grid, D_grid, n_stars, levels=100)
plt.colorbar(label="Abundance")
plt.xscale('log')
plt.show()
"""Vs = np.linspace(0.24, 0.5, 10)
V_grid, D_grid = np.meshgrid(Vs, Ds)
K = 0.1
r_stars = -D_grid * K / (D_grid - V_grid)
plt.contourf(V_grid, D_grid, r_stars, levels=100, cmap="RdYlBu")
plt.colorbar(label="R* [mM]")
#plt.contour(V_grid, D_grid, r_stars, levels=[1])
plt.show()
n_stars = D_grid * q * (K * M + r_stars *
                        (-K + M - r_stars)) / (r_stars * V_grid)
plt.contourf(V_grid, D_grid, n_stars, levels=100)
plt.colorbar(label="Abundance")
#plt.xscale('log')
plt.show()
"""
