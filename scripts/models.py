import numpy as np
from scipy.integrate import odeint

#xs = np.linspace(0, 10000, 10000)


def species_1_niche_1(t, f_dN1_1, f_dR1_1, params):

    def model(y, t, f_dN1_1, f_dR1_1, params):
        N1, R1 = y
        return [
            f_dN1_1(N1, R1, *params.values()),
            f_dR1_1(N1, R1, *params.values())
        ]

    return odeint(
        model,
        [params["N01"], params["M1"]],
        t,
        args=(f_dN1_1, f_dR1_1, params),
    )


def species_1_niche_2(t, f_dN1_2, f_dR1_2, params):

    def model(y, t, f_dN1_2, f_dR1_2, params):
        N1, R2 = y
        return [
            f_dN1_2(N1, R2, *params.values()),
            f_dR1_2(N1, R2, *params.values()),
        ]

    return odeint(
        model,
        [params["N01"], params["M2"]],
        t,
        args=(f_dN1_2, f_dR1_2, params),
    )


def species_2_niche_1(t, f_dN2_1, f_dR2_1, params):

    def model(y, t, f_dN2_1, f_dR2_1, params):
        N2, R1 = y
        return [
            f_dN2_1(N2, R1, *params.values()),
            f_dR2_1(N2, R1, *params.values()),
        ]

    return odeint(
        model,
        [params["N02"], params["M1"]],
        t,
        args=(f_dN2_1, f_dR2_1, params),
    )


def species_2_niche_2(t, f_dN2_2, f_dR2_2, params):

    def model(y, t, f_dN2_2, f_dR2_2, params):
        N2, R2 = y
        return [
            f_dN2_2(N2, R2, *params.values()),
            f_dR2_2(N2, R2, *params.values()),
        ]

    return odeint(
        model,
        [params["N02"], params["M2"]],
        t,
        args=(f_dN2_2, f_dR2_2, params),
    )


def species_12_niche_1(t, f_dN1_1, f_dN2_1, f_dR12_1, params):

    def model(y, t, f_dN1_1, f_dN2_1, f_dR12_1, params):
        N1, N2, R1 = y
        return [
            f_dN1_1(N1, R1, *params.values()),
            f_dN2_1(N2, R1, *params.values()),
            f_dR12_1(N1, N2, R1, *params.values()),
        ]

    return odeint(
        model,
        [params["N01"], params["N02"], params["M1"]],
        t,
        args=(f_dN1_1, f_dN2_1, f_dR12_1, params),
    )


def species_12_niche_2(t, f_dN1_2, f_dN2_2, f_dR12_2, params):

    def model(y, t, f_dN1_2, f_dN2_2, f_dR12_2, params):
        N1, N2, R2 = y
        return [
            f_dN1_2(N1, R2, *params.values()),
            f_dN2_2(N2, R2, *params.values()),
            f_dR12_2(N1, N2, R2, *params.values()),
        ]

    return odeint(
        model,
        [params["N01"], params["N02"], params["M2"]],
        t,
        args=(f_dN1_2, f_dN2_2, f_dR12_2, params),
    )


def species_12_niche_12(t, f_dN1_1, f_dN2_2, f_dR1_1, f_dR2_2, params):

    def model(y, t, f_dN1_1, f_dN2_2, f_dR1_1, f_dR2_2, params):
        N1, N2, R1, R2 = y
        return [
            f_dN1_1(N1, R1, *params.values()),
            f_dN2_2(N2, R2, *params.values()),
            f_dR1_1(N1, R1, *params.values()),
            f_dR2_2(N2, R2, *params.values()),
        ]

    return odeint(
        model,
        [params["N01"], params["N02"], params["M1"], params["M2"]],
        t,
        args=(f_dN1_1, f_dN2_2, f_dR1_1, f_dR2_2, params),
    )


def species_12_niche_21(t, f_dN1_2, f_dN2_1, f_dR1_2, f_dR2_1, params):

    def model(y, t, f_dN1_2, f_dN2_1, f_dR1_2, f_dR2_1, params):
        N1, N2, R1, R2 = y
        return [
            f_dN1_2(N1, R2, *params.values()),
            f_dN2_1(N2, R1, *params.values()),
            f_dR2_1(N2, R1, *params.values()),
            f_dR1_2(N1, R2, *params.values()),
        ]

    return odeint(
        model,
        [params["N01"], params["N02"], params["M1"], params["M2"]],
        t,
        args=(f_dN1_2, f_dN2_1, f_dR1_2, f_dR2_1, params),
    )


def species_12_niche_13(t, f_dN1_1, f_dN2_13, f_dR12_13, f_dR2_3, params):

    def model(y, t, f_dN1_1, f_dN2_13, f_dR12_13, f_dR2_3, params):
        N1, N2, R1, R3 = y
        return [
            f_dN1_1(N1, R1, *params.values()),
            f_dN2_13(N2, R1, R3, *params.values()),
            f_dR12_13(N1, N2, R1, *params.values()),
            f_dR2_3(N1, N2, R1, R3, *params.values()),
        ]

    return odeint(
        model,
        [params["N01"], params["N02"], params["M1"], params["M3"]],
        t,
        args=(f_dN1_1, f_dN2_13, f_dR12_13, f_dR2_3, params),
    )
