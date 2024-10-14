from sympy import symbols, Eq, solve
from sympy import init_printing
import numpy as np
import pandas as pd
import plotly.graph_objects as go

init_printing()


"""
This is a chemostat math library for two species including cross-feeding of essential resources.
It's written for two species with up to two resources nad one cross-fed metabolite.
The index corresponds to the species or resource index. Index 3 is used for the cross-fed metabolite
"""

# Symbols
equations = symbols("dN1,dN2,dR1,dR2,dR3,J1_1,J1_2,J2_1,J2_2,J2_3")
dN1, dN2, dR1, dR2, dR3, J1_1, J1_2, J2_1, J2_2, J2_3 = equations
species = symbols("N1,N2")
N1, N2 = species
resources = symbols("R1,R2,R3")
R1, R2, R3 = resources
growth_parameters = symbols(
    "v1_1,v1_2,v2_1,v2_2,v2_3,K1_1,K1_2,K2_1,K2_2,K2_3,q1_1,q1_2,q2_1,q2_2,q2_3"
)
(
    v1_1,
    v1_2,
    v2_1,
    v2_2,
    v2_3,
    K1_1,
    K1_2,
    K2_1,
    K2_2,
    K2_3,
    q1_1,
    q1_2,
    q2_1,
    q2_2,
    q2_3,
) = growth_parameters
metabolite_production = symbols("a")
a = metabolite_production
medium = symbols("M1,M2,M3")
M1, M2, M3 = medium
D = symbols("D")

# Per capita growth rates
J1_1 = Eq(J1_1, v1_1 * R1 / (R1 + K1_1))
J1_2 = Eq(J1_2, v1_2 * R2 / (R2 / K1_2))
J2_1 = Eq(J2_1, v2_1 * R1 / (R1 + K2_1))
J2_2 = Eq(J2_2, v2_2 * R2 / (R2 + K2_2))
J2_3 = Eq(J2_3, v2_2 * R2 / (R2 + K2_2) * R3 / (R3 + K2_3))

# Species
# Species 1 consumes resource 1
N1_1 = Eq(dN1, J1_1.rhs * N1 - D * N1)
N1_2 = Eq(dN1, J1_2.rhs * N1 - D * N1)
N2_1 = Eq(dN2, J2_1.rhs * N2 - D * N2)
N2_2 = Eq(dN2, J2_2.rhs * N2 - D * N2)

# Resources
R1_1 = Eq(dR1, D * M1 - J1_1.rhs * N1 / q1_1 - D * R1)
R1_2 = Eq(dR2, D * M2 - J1_2.rhs * N1 / q1_2 - D * R2)
R2_1 = Eq(dR1, D * M1 - J2_1.rhs * N2 / q2_1 - D * R1)
R2_2 = Eq(dR2, D * M2 - J2_2.rhs * N2 / q2_2 - D * R2)
R1_2_1 = Eq(dR1, D * M1 - J1_1.rhs * N1 / q1_1 - J2_1.rhs * N2 / q2_1 - D * R1)
R1_2_2 = Eq(dR2, D * M1 - J1_2.rhs * N1 / q1_2 - J2_2.rhs * N2 / q2_2 - D * R2)


# Flow two species two niches
def two_species_two_niches():
    R1_eq = solve(Eq(D, J1_1.rhs), R1)[0]
    R2_eq = solve(Eq(D, J2_2.rhs), R2)[0]

    return solve(Eq(0, R1_1.rhs), N1)[0].subs({R1: R1_eq}).simplify(), solve(
        Eq(0, R2_2.rhs), N2
    )[0].subs({R2: R2_eq}).simplify()


# Flow two species one niche
def two_species_niche_1():
    R1_1_eq = solve(Eq(D, J1_1.rhs), R1)[0]
    R2_1_eq = solve(Eq(D, J2_1.rhs), R1)[0]
    return solve(Eq(0, R1_2_1.rhs), N1)[0].subs({R1: R1_1_eq}).simplify(), solve(
        Eq(0, R1_2_1.rhs), N2
    )[0].subs({R1: R2_1_eq}).simplify()


def two_species_niche_2():
    R1_2_eq = solve(Eq(D, J1_2.rhs), R2)[0]
    R2_2_eq = solve(Eq(D, J2_2.rhs), R2)[0]
    return solve(Eq(0, R1_2_2.rhs), N1)[0].subs({R2: R1_2_eq}).simplify(), solve(
        Eq(0, R1_2_2.rhs), N2
    )[0].subs({R2: R2_2_eq}).simplify()

def two_species_cross_feeding():
    


df = pd.read_csv("parameters.csv")
params = {key: value for key, value in zip(df["parameter"], df["value"])}


def plot_rel_abund(N1, N2):
    Ds = np.linspace(0, 0.2, 100)
    N1_array = []
    N2_array = []
    for D in Ds:
        params["D"] = D
        N1_abs = N1.subs(params)
        if not N1_abs.is_real:
            N1_array.append(0)
        elif N1_abs < 0:
            N1_array.append(0)
        else:
            N1_array.append(float(N1_abs))

        N2_abs = N2.subs(params)
        if not N2_abs.is_real:
            N2_array.append(0)
        elif N2_abs < 0:
            N2_array.append(0)
        else:
            N2_array.append(float(N2_abs))
    N1_rel = np.array(N1_array) / (np.array(N1_array) + np.array(N2_array))
    N2_rel = np.array(N2_array) / (np.array(N1_array) + np.array(N2_array))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Ds, y=N1_rel, name="N1"))
    fig.add_trace(go.Scatter(x=Ds, y=N2_rel, name="N2"))
    fig.show()


N1, N2 = two_species_two_niches()
plot_rel_abund(N1, N2)
