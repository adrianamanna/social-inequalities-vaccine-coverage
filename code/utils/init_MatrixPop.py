import numpy as np
import pandas as pd
import pickle

from utils.matrix_tools_sy import (
    get_Pd1_Md1,
    get_P_d1d2,
    get_M_d1d2,
    assortative_mixing,
    random_mixing,
)


def init_MatrixAge(country, sim_type, nw=None):
    if sim_type == 'synthetic':
        data = pd.read_pickle("./data/data.pkl")
        P_d1, N_pop, M_d1, M_d1_tot = get_Pd1_Md1(data, country)
        N_d1 = {key: P_d1[key] * N_pop for key in P_d1.keys()}
    
    elif sim_type == "mazsk":
        with open("./data/init_population/mazsk/Ms_" + nw + ".pkl", "rb") as f:
            Ms = pickle.load(f)
        with open("./data/init_population/mazsk/Ns_" + nw + ".pkl", "rb") as f:
            Ns = pickle.load(f)

        # re-naming keys so that age ranges from 0-7 and ses from 0 to3
        N_d1 = {i: v for i, v in enumerate(Ns["d1"].values())}
        M_d1 = Ms["d1"]
        M_d1_tot = Ms["tot_d1"]

    lev_d1 = list(N_d1.keys())
    return N_d1, M_d1, M_d1_tot, lev_d1


def init_PopSy_v2(country, n_dim2, dist_Pd2):
    data = pd.read_pickle("./data/data.pkl")
    P_d1, N_pop, M_d1, M_d1_tot = get_Pd1_Md1(data, country)
    N_d1 = {key: P_d1[key] * N_pop for key in P_d1.keys()}

    N_d2 = {key: dist_Pd2[key] * N_pop for key in range(len(dist_Pd2))}

    P_d1d2 = get_P_d1d2(P_d1, n_dim2, dist_Pd2)
    N_d1d2 = {key: P_d1d2[key] * N_pop for key in P_d1d2.keys()}

    lev_d1 = np.unique([int(_[0]) for _ in list(N_d1d2.keys())])
    lev_d2 = np.unique([int(_[1]) for _ in list(N_d1d2.keys())])
    lev_d1d2 = list(N_d1d2.keys())

    return N_d1, N_d2, N_d1d2, lev_d1, lev_d2, lev_d1d2, M_d1_tot


def init_MatrixSy_v2(M_d1_tot, N_d1d2, n_dim2, mixing_type, args):
    if mixing_type == "assortative_mixing":
        r, activity, ds = args
        m_diag, m_off_diag = assortative_mixing(r, activity, ds)

    elif mixing_type == "random_mixing":
        dist_Pd2 = args[0]
        m_diag, m_off_diag = random_mixing(dist_Pd2)
    else:
        raise ValueError("- mixing_type - not known.")

    M_d1d2_tot, M_d1d2 = get_M_d1d2(M_d1_tot, n_dim2, N_d1d2, m_diag, m_off_diag)

    return M_d1d2


def init_MatrixData(nw):
    with open("./data/init_population/mazsk/Ms_" + nw + ".pkl", "rb") as f:
        Ms = pickle.load(f)
    with open("./data/init_population/mazsk/Ns_" + nw + ".pkl", "rb") as f:
        Ns = pickle.load(f)

    # re-naming keys so that age ranges from 0-7 and ses from 0 to3
    N_d1 = {i: v for i, v in enumerate(Ns["d1"].values())}
    N_d2 = {i: v for i, v in enumerate(Ns["d2"].values())}

    keys_d1d2 = [(i, j) for i in N_d1.keys() for j in N_d2.keys()]
    N_d1d2 = {(i, j): v for (i, j), v in zip(keys_d1d2, Ns["d1d2"].values())}
    M_d1d2 = {(i, j): v for (i, j), v in zip(keys_d1d2, Ms["d1d2"].values())}
    for ij in M_d1d2.keys():
        M_d1d2[ij] = {(i, j): v for (i, j), v in zip(keys_d1d2, M_d1d2[ij].values())}

    lev_d1 = np.unique([int(_[0]) for _ in list(N_d1d2.keys())])
    lev_d2 = np.unique([int(_[1]) for _ in list(N_d1d2.keys())])
    lev_d1d2 = list(N_d1d2.keys())
    return N_d1d2, N_d1, N_d2, M_d1d2, lev_d1, lev_d2, lev_d1d2
