import pandas as pd
import numpy as np
import math

from utils.init_MatrixPop import init_PopSy_v2, init_MatrixData


def load_RES_Gab(E, M, NPI_S, ifr, type_res, country):
    # - MODEL Gab index of result:
    # [A] p0 res
    # [B] [p2, res] = 1 means we look at res
    # [D] dict NSEIRD Inew Dnew
    # [E] [0: nonvax, 1:vax]

    RES_Gab_r0s = pd.read_pickle(
        f"./res/epi/prep/PREP_{type_res}_Gab_{E}_{M}_{NPI_S}_IFR.pkl"
    )

    return RES_Gab_r0s


def load_RES(E, M, M_Cij, type_res):
    NPI_S = "" 
    # - MODEL Gab index of result:
    # [A] p0 res
    # [B] [p2, res] = 1 means we look at res
    # [C] R0 indx
    # [D] dict NSEIRD Inew Dnew
    # [E] [0: nonvax, 1:vax]

    RES_Gab_r0s = pd.read_pickle(
        f"./res/PREP_{type_res}_Gab_{E}_{M}_{NPI_S}.pkl"
    )

    # - MODEL Cij
    # [A] R0 indx
    # [B] 0:res - 1:res age
    # [C] dict NSEIRD Inew Dnew
    # [D] [0: nonvax, 1:vax]
    res_Cij_r0s = pd.read_pickle(
        f"./res/PREP_all_Cij_{E}_{M_Cij}_{NPI_S}.pkl"
    )

    return RES_Gab_r0s, res_Cij_r0s


def load_N(config_mat_var, config_mat, country):
    if config_mat_var["sim_type"] == "synthetic":
        N_d1, N_d2, N_d1d2, lev_d1, lev_d2, lev_d1d2, M_d1_tot = init_PopSy_v2(
            country, config_mat["n_dim2"], config_mat_var["dist_Pd2"]
        )
    elif config_mat_var["sim_type"] == "mazsk":
        N_d1d2, N_d1, N_d2, M_d1d2, lev_d1, lev_d2, lev_d1d2 = init_MatrixData(
            config_mat_var["nw"]
        )

    return N_d1d2, N_d1, N_d2


# function for rates
def load_RATES(E, M, NPI_S, model):
    RES = pd.read_pickle(
        f"./res/RATES_{model}_{E}_{M}_{NPI_S}.pkl"
    )
    return RES


# --
def des_epi(epi_df, ses_labs, C):
    res = []
    grouped_values = epi_df[C].apply(
        lambda x: pd.DataFrame(x).groupby(level=1, axis=1).sum().to_dict()
    )

    for ses in ses_labs:
        c_values = np.vstack(
            grouped_values.apply(
                lambda x: np.array([_ for _ in x[ses].values()]).flatten()
            )
        )
        res.append(
            pd.DataFrame(
                c_values,
                columns=["{}_{}".format(C, i) for i in range(c_values.shape[1])],
            ).describe()
        )
    return res


def des_epi_age(epi_df_, labs, C):
    res = []
    grouped_values = epi_df_[C].apply(
        lambda x: pd.DataFrame(x).groupby(level=0, axis=1).sum().to_dict()
    )
    for age_group in labs:
        c_values = np.vstack(
            grouped_values.apply(
                lambda x: np.array([_ for _ in x[age_group].values()]).flatten()
            )
        )
        res.append(
            pd.DataFrame(
                c_values,
                columns=["{}_{}".format(C, i) for i in range(c_values.shape[1])],
            ).describe()
        )
    return res


def des_d1d2(epi_df, key, C):
    # epi_df = epi_df_.copy()
    keys = epi_df[C][0].keys()
    res = []
    for key in keys:
        epi_df[f"{C}_{key}"] = epi_df[C].apply(
            lambda x: np.array([_ for _ in x[key].values()]).flatten()
        )
        res.append(
            [key, pd.DataFrame(epi_df["{}_{}".format(C, key)].to_list()).describe()]
        )
    return res



def des_epi_all(epi_df,C): 
    #epi_df = epi_df_.copy() 
    epi_df[f'{C}_all'] = epi_df[C].apply(lambda x: pd.DataFrame(x).sum(axis = 1).to_dict()) 
    res_all = pd.DataFrame(epi_df[f'{C}_all'].to_list()).describe() 
    return res_all

def aggregate_d1(res, COMP):
    return np.sum(res[COMP], axis=1)


def aggregate_d2(res, COMP):
    return np.sum(res[COMP], axis=0)


def aggregate_all(res, COMP):
    # print(res[COMP])
    return np.sum(res[COMP])


def des_rates(RATES):
    p75 = np.percentile(RATES, 75, axis=0)
    p25 = np.percentile(RATES, 25, axis=0)
    p50 = np.percentile(RATES, 50, axis=0)
    return p25, p50, p75


def prep_RATES_all(RES, ip, rate_type):
    niter = range(len(RES[ip]))
    res = RES[ip]

    if rate_type == "AR":
        COMPv, COMPnv = "Ev_newT", "E_newT"
    else:
        COMPv, COMPnv = "DvT", "DT"

    RATE_v = [
        aggregate_all(res[i], COMPv)
        / (
            aggregate_all(res[i], "NvT")
            - aggregate_all(res[i], "E_vaxT")
            - aggregate_all(res[i], "R_vaxT")
        )
        for i in niter
    ]
    RATE_nv = [
        aggregate_all(res[i], COMPnv)
        / (
            aggregate_all(res[i], "NT")
            + aggregate_all(res[i], "E_vaxT")
            + aggregate_all(res[i], "R_vaxT")
        )
        for i in niter
    ]
    RATE = [
        (aggregate_all(res[i], COMPv) + aggregate_all(res[i], COMPnv))
        / (aggregate_all(res[i], "NvT") + aggregate_all(res[i], "NT"))
        for i in niter
    ]

    RATE_v = des_rates(RATE_v)
    RATE_nv = des_rates(RATE_nv)
    RATE = des_rates(RATE)

    return RATE_v, RATE_nv, RATE


def prep_RATES_d1(RES, ip, rate_type):
    niter = range(len(RES[ip]))
    res = RES[ip]

    if rate_type == "AR":
        COMPv, COMPnv = "Ev_newT", "E_newT"
    else:
        COMPv, COMPnv = "DvT", "DT"

    RATE_v = [
        aggregate_d1(res[i], COMPv)
        / (
            aggregate_d1(res[i], "NvT")
            - aggregate_d1(res[i], "E_vaxT")
            - aggregate_d1(res[i], "R_vaxT")
        )
        for i in niter
    ]
    RATE_nv = [
        aggregate_d1(res[i], COMPnv)
        / (
            aggregate_d1(res[i], "NT")
            + aggregate_d1(res[i], "E_vaxT")
            + aggregate_d1(res[i], "R_vaxT")
        )
        for i in niter
    ]
    RATE = [
        (aggregate_d1(res[i], COMPv) + aggregate_d1(res[i], COMPnv))
        / (aggregate_d1(res[i], "NvT") + aggregate_d1(res[i], "NT"))
        for i in niter
    ]

    RATE_v = des_rates(RATE_v)
    RATE_nv = des_rates(RATE_nv)
    RATE = des_rates(RATE)

    return RATE_v, RATE_nv, RATE


def prep_RATES_d2(RES, ip, rate_type):
    niter = range(len(RES[ip]))
    res = RES[ip]

    if rate_type == "AR":
        COMPv, COMPnv = "Ev_newT", "E_newT"
    else:
        COMPv, COMPnv = "DvT", "DT"

    RATE_v = [
        aggregate_d2(res[i], COMPv)
        / (
            aggregate_d2(res[i], "NvT")
            - aggregate_d2(res[i], "E_vaxT")
            - aggregate_d2(res[i], "R_vaxT")
        )
        for i in niter
    ]
    RATE_nv = [
        aggregate_d2(res[i], COMPnv)
        / (
            aggregate_d2(res[i], "NT")
            + aggregate_d2(res[i], "E_vaxT")
            + aggregate_d2(res[i], "R_vaxT")
        )
        for i in niter
    ]
    RATE = [
        (aggregate_d2(res[i], COMPv) + aggregate_d2(res[i], COMPnv))
        / (aggregate_d2(res[i], "NvT") + aggregate_d2(res[i], "NT"))
        for i in niter
    ]

    RATE_v = des_rates(RATE_v)
    RATE_nv = des_rates(RATE_nv)
    RATE = des_rates(RATE)

    return RATE_v, RATE_nv, RATE


###


def rates_d1d2_inferred(RES_Cij, RES, ip, rate_type, N_d2):
    # this new function accounts for tha age dimention
    lev_d2 = range(len(N_d2))
    niter = range(len(RES[ip]))

    res_Cij = RES_Cij[0]
    res_Gab = RES[ip]

    if rate_type == "AR":
        COMPv, COMPnv = "Ev_newT", "E_newT"
    else:
        COMPv, COMPnv = "DvT", "DT"

    RATECij_v_age, RATECij_nv_age = [], []
    NvT_Gab_age, NT_Gab_age = [], []

    for age in range(8):
        RATECij_v_i = [
            res_Cij[i][COMPv][age]
            / (
                res_Cij[i]["NvT"][age]
                - res_Cij[i]["E_vaxT"][age]
                - res_Cij[i]["R_vaxT"][age]
            )
            for i in niter
        ]
        RATECij_v_i = [0 if math.isnan(x) else x for x in RATECij_v_i]

        RATECij_nv_i = [
            res_Cij[i][COMPnv][age]
            / (
                res_Cij[i]["NT"][age]
                + res_Cij[i]["E_vaxT"][age]
                + res_Cij[i]["R_vaxT"][age]
            )
            for i in niter
        ]

        NvT_Gab_i = [res_Gab[i]["NvT"][age] for i in niter]
        NT_Gab_i = [res_Gab[i]["NT"][age] for i in niter]

        RATECij_v_age.append(RATECij_v_i)
        RATECij_nv_age.append(RATECij_nv_i)
        NvT_Gab_age.append(NvT_Gab_i)
        NT_Gab_age.append(NT_Gab_i)

    RATE_d2_infered = [
        [
            sum(
                [
                    (
                        RATECij_v_age[age][i] * NvT_Gab_age[age][i][d2]
                        + RATECij_nv_age[age][i] * NT_Gab_age[age][i][d2]
                    )
                    for age in range(8)
                ]
            )
            / N_d2[d2]
            for d2 in lev_d2
        ]
        for i in niter
    ]

    # RATE_d2_infered = [[(RATECij_v[i]*NvT_Gab[i][d2]+RATECij_nv[i]*NT_Gab[i][d2])/N_d2[d2] for d2 in lev_d2] for i in niter]
    RATE_d2_infered = des_rates(RATE_d2_infered)

    return RATE_d2_infered  # RATE_d2_infered RATECij_v, RATECij_nv, NvT_Gab,
