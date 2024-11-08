from sympy import symbols, Eq, solve
from sympy import init_printing
import pandas as pd
from sympy import lambdify

init_printing()

# Load parameters
df = pd.read_csv("parameters.csv")
params = {key: value for key, value in zip(df["parameter"], df["value"])}
"""
This is a chemostat math library for two species including cross-feeding of essential resources.
It's written for two species with up to two resources nad one cross-fed metabolite.
The index corresponds to the species or resource index. Index 3 is used for the cross-fed metabolite
"""

# Symbols
equations = symbols(
    "dN1,dN2,dN12_1,dN12_2,dN12_13,dR1,dR2,dR3,J1_1,J1_2,J2_1,J2_2,J2_3")

dN1, dN2, dN12_1, dN12_2, dN12_13, dR1, dR2, dR3, J1_1, J1_2, J2_1, J2_2, J2_3 = equations
species = symbols("N1,N2")
N1, N2 = species
resources = symbols("R1,R2,R3")
R1, R2, R3 = resources
growth_parameters = symbols(
    "v1_1,v1_2,v2_1,v2_2,v2_3,K1_1,K1_2,K2_1,K2_2,K2_3,q1_1,q1_2,q2_1,q2_2,q2_3"
)

v1_1, v1_2, v2_1, v2_2, v2_3, K1_1, K1_2, K2_1, K2_2, K2_3, q1_1, q1_2, q2_1, q2_2, q2_3 = growth_parameters
metabolite_production = symbols("a")
a = metabolite_production
medium = symbols("M1,M2,M3")
M1, M2, M3 = medium
D = symbols("D")

# Per capita growth rates
eq_J1_1 = Eq(J1_1, v1_1 * R1 / (R1 + K1_1))
eq_J1_2 = Eq(J1_2, v1_2 * R2 / (R2 + K1_2))
eq_J2_1 = Eq(J2_1, v2_1 * R1 / (R1 + K2_1))
eq_J2_2 = Eq(J2_2, v2_2 * R2 / (R2 + K2_2))
eq_J2_3 = Eq(J2_3, v2_3 * R3 / (R3 + K2_3))

# Species dynamics
# Equations
eq_dN1_1 = Eq(dN1, eq_J1_1.rhs * N1 - D * N1)
eq_dN1_2 = Eq(dN1, eq_J1_2.rhs * N1 - D * N1)
eq_dN2_1 = Eq(dN2, eq_J2_1.rhs * N2 - D * N2)
eq_dN2_2 = Eq(dN2, eq_J2_2.rhs * N2 - D * N2)
eq_dN12_1 = Eq(dN12_1, eq_J1_1.rhs * N1 + eq_J2_1.rhs * N2 - D * N1 - D * N2)
eq_dN12_2 = Eq(dN12_2, eq_J1_2.rhs * N1 + eq_J2_2.rhs * N2 - D * N1 - D * N2)
eq_dN2_13 = Eq(dN2, eq_J2_3.rhs * eq_J2_1.rhs * N2 - D * N2)

# Functions
f_dN1_1 = lambdify((N1, R1, *params.keys()), eq_dN1_1.rhs, "numpy")
f_dN1_2 = lambdify((N1, R2, *params.keys()), eq_dN1_2.rhs, "numpy")
f_dN2_1 = lambdify((N2, R1, *params.keys()), eq_dN2_1.rhs, "numpy")
f_dN2_2 = lambdify((N2, R2, *params.keys()), eq_dN2_2.rhs, "numpy")
f_dN12_1 = lambdify((N1, N2, R1, *params.keys()), eq_dN12_1.rhs, "numpy")
f_dN12_2 = lambdify((N1, N2, R2, *params.keys()), eq_dN12_2.rhs, "numpy")
f_dN2_13 = lambdify((N2, R1, R3, *params.keys()), eq_dN2_13.rhs, "numpy")

# Resources
# Equations
eq_dR1_1 = Eq(dR1, D * M1 - eq_J1_1.rhs * N1 / q1_1 - D * R1)
eq_dR1_2 = Eq(dR2, D * M2 - eq_J1_2.rhs * N1 / q1_2 - D * R2)
eq_dR2_1 = Eq(dR1, D * M1 - eq_J2_1.rhs * N2 / q2_1 - D * R1)
eq_dR2_2 = Eq(dR2, D * M2 - eq_J2_2.rhs * N2 / q2_2 - D * R2)
eq_dR12_1 = Eq(
    dR1, D * M1 - eq_J1_1.rhs * N1 / q1_1 - eq_J2_1.rhs * N2 / q2_1 - D * R1)
eq_dR12_2 = Eq(
    dR2, D * M2 - eq_J1_2.rhs * N1 / q1_2 - eq_J2_2.rhs * N2 / q2_2 - D * R2)
eq_dR2_3 = Eq(dR3,
              a * eq_J1_1.rhs * N1 / q1_1 - eq_J2_3.rhs * N2 / q2_3 - D * R3)

# Functions
f_dR1_1 = lambdify((N1, R1, *params.keys()), eq_dR1_1.rhs, "numpy")
f_dR1_2 = lambdify((N1, R2, *params.keys()), eq_dR1_2.rhs, "numpy")
f_dR2_1 = lambdify((N2, R1, *params.keys()), eq_dR2_1.rhs, "numpy")
f_dR2_2 = lambdify((N2, R2, *params.keys()), eq_dR2_2.rhs, "numpy")
f_dR2_3 = lambdify((N1, N2, R1, R3, *params.keys()), eq_dR2_3.rhs, "numpy")
f_dR12_1 = lambdify((N1, N2, R1, *params.keys()), eq_dR12_1.rhs, "numpy")
f_dR12_2 = lambdify((N1, N2, R2, *params.keys()), eq_dR12_2.rhs, "numpy")
