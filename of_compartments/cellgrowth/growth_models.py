import numpy as np

from .compartment_model import Case


def growthModel_xing_simplified(t: float, z: np.ndarray, p: dict, case: Case):
    """
    Growth model based on Xing, Bishop & Leister et al. (2010)
    Modeling Kinetics of a Large-Scale Fed-Batch CHO Cell Culture
    by Markov Chain Monte Carlo Method, Biotechnology Progress.

    Neglects aspects related to monoclonal antibody production as well
    as the concentration threshold on ammonia and lactate inhibition.
    """

    S = z[0, :]
    G = z[1, :]
    L = z[2, :]
    A = z[3, :]
    DO = z[4, :]
    X = z[5, :]
    mu = (
        p["mu_m"]
        * (S / (S + p["K_S"]))
        * (G / (G + p["K_G"]))
        * (p["KI_L"] / (L + p["KI_L"]))
        * (p["KI_A"] / (A + p["KI_A"]))
    )
    mu_d = p["mu_dm"] * (L / (L + p["KD_L"])) * (A / (A + p["KD_A"]))
    m_G = p["a_1"] * G / (p["a_2"] + G)
    dS = -X * (p["m_S"] + (mu) / p["Y_XS"])
    dG = -X * (m_G + (mu) / p["Y_XG"]) - p["d_G"] * G
    dL = X * (p["m_S"] + (mu) / p["Y_XS"]) * p["Y_LS"]
    dA = X * ((mu) / p["Y_XG"]) * p["Y_AG"] - p["r_A"] * X + p["d_G"] * G
    kLas = case.compartment_data.kLa

    dO = kLas * (p["DO_eq"] - DO) - X * p["OUR"]
    dS[S <= 0] = 0
    dL[S <= 0] = 0
    dG[G <= 0] = 0

    growth_deltas = np.zeros(z.shape)
    growth_deltas[0, :] = dS
    growth_deltas[1, :] = dG
    growth_deltas[2, :] = dL
    growth_deltas[3, :] = dA
    growth_deltas[4, :] = dO
    growth_deltas[5, :] = X * (mu - mu_d) - X * case.k_d
    return growth_deltas


params_xing_simplified = {
    "m_S": 6.92e-11,  # mmol / cell / h
    "a_1": 3.2e-12,  # mmol / cell / h
    "a_2": 2.1,  # mM
    "mu_m": 0.029,  # 1 / h
    "mu_dm": 0.016,  # 1 / h
    "K_S": 0.084,  # mM
    "K_G": 0.047,  # mM
    "KI_L": 43,  # mM
    "KI_A": 6.51,  # mM
    "KD_L": 45.8,  # mM
    "KD_A": 6.51,  # mM
    "d_G": 7.2e-3,  # 1 / h
    "r_A": 6.3e-12,  # mmol / cell / h
    "Y_XS": 1.69e8,  # cell / mmol
    "Y_XG": 9.74e8,  # cell / mmol
    "Y_LS": 1.23,  # mol / mol
    "Y_AG": 0.67,  # mol / mol
    "DO_eq": 1.07,  # mM
    "OUR": 3.2e-10,  # mmol / cell h
}
