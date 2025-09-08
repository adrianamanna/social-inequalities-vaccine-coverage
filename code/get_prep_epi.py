import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from utils.data_tools import des_epi, des_epi_all

from utils.EpiModels import EpiModel  # note model
from utils.init_EpiParams import init_model_params
from utils.tools import upload_yaml
import os

os.makedirs("res", exist_ok=True)

# . init config
config_params_torun = upload_yaml("config_params_torun")
params_to_run = config_params_torun["params_to_run"]

config_mat = upload_yaml("config_matrix")
config_epi = upload_yaml("config_epi")
params = upload_yaml("parameters")
country = config_mat["country"]

model_type = config_params_torun["model_type"]  # Gab, Cij
SCENARIOS = config_params_torun[params_to_run]

# . init vax priorities
PRIORITY_TO_RUN = config_epi["PRIORITIES_D2"][model_type]
print(PRIORITY_TO_RUN)


def run_epi(E, M, NPI_S, pvax_d2):
    config_mat_var = params["matrix_params"][M]
    config_epi_var = params["epi_params"][E]
    r0 = config_epi_var["R0s"][0]
    priority_dim2 = PRIORITY_TO_RUN[pvax_d2]

    model_args_len, vax_args_len, args, IFR, NPI_bool, M_npi, T_npi = init_model_params(
        model_type,
        NPI_S,
        config_mat,
        config_epi,
        config_mat_var,
        config_epi_var,
        params,
        priority_dim2,
    )

    model = EpiModel(
        model_type,
        config_epi_var["vaccination_type"],
        config_epi_var["vax_scenario"],
        config_epi["vax_saturation"],
        IFR,
        config_epi["mu_val"],
        config_epi["eps_val"],
        r0,
        config_epi["seed"],
        config_epi["VE"],
        config_epi_var["g1"],
        config_epi["g2"],
        config_epi_var["stop"],
        config_epi_var["t0"],
        config_epi_var["current_immunity"],
        model_args_len,
        vax_args_len,
        *args,
        NPI=NPI_bool,
        M_npi=M_npi,
        T_npi=T_npi,
    )

    RES = np.array(model.run_stoch(config_epi["N_iter"]), dtype=object)
    RES = pd.DataFrame(RES)
    RES.columns = [
        "S",
        "Sv",
        "E",
        "Ev",
        "I",
        "Iv",
        "R",
        "Rv",
        "D",
        "Dv",
        "E_new",
        "Ev_new",
        "I_new",
        "Iv_new",
        "D_new",
        "Dv_new",
        "N",
        "Nv",
        "t",
        "exceeding_vax_key",
        "S_vax",
        "E_vax",
        "R_vax",
    ]
    return RES


def prep_RES_all(RES):
    CCv = {}
    for _ in ["N", "S", "E", "I", "R", "D"]:
        C = des_epi_all(RES, _)
        Cv = des_epi_all(RES, "{}v".format(_))
        CCv["{}".format(_)] = [C, Cv]

    for _ in ["E_new", "I_new", "D_new"]:
        C = des_epi_all(RES, _)
        Cv = des_epi_all(RES, _[0] + "v_new")
        CCv["{}".format(_)] = [C, Cv]

    for _ in ["S_vax", "E_vax", "R_vax"]:
        C = des_epi_all(RES, _)
        CCv["{}".format(_)] = "NONE"
        CCv["{}".format(_)] = [C, Cv]

    return CCv


def prep_RES_ses(RES):
    lev_d2 = [0, 1, 2]
    CCv = {}
    for _ in ["N", "S", "E", "I", "R", "D"]:
        C = des_epi(RES, lev_d2, _)
        Cv = des_epi(RES, lev_d2, "{}v".format(_))
        CCv["{}".format(_)] = [C, Cv]

    for _ in ["E_new", "I_new", "D_new"]:
        C = des_epi(RES, lev_d2, _)
        Cv = des_epi(RES, lev_d2, _[0] + "v_new")
        CCv["{}".format(_)] = [C, Cv]

    for _ in ["S_vax", "E_vax", "R_vax"]:
        C = des_epi(RES, lev_d2, _)
        CCv["{}".format(_)] = "NONE"
        CCv["{}".format(_)] = [C, Cv]
    return CCv


# def prep_RES_age(RES):
#    lev_d1 = [0, 1, 2, 3, 4, 5, 6, 7]
#    CCv = {}
#    for _ in ["N", "S", "E", "I", "R", "D"]:
#        C = des_epi_age(RES, lev_d1, _)
#        Cv = des_epi_age(RES, lev_d1, "{}v".format(_))
#        CCv["{}".format(_)] = [C, Cv]
#    for _ in ["E_new", "I_new", "D_new"]:
#        C = des_epi_age(RES, lev_d1, _)
#        Cv = des_epi_age(RES, lev_d1, _[0] + "v_new")
#        CCv["{}".format(_)] = [C, Cv]
#
#    for _ in ["S_vax", "E_vax", "R_vax"]:
#        C = des_epi_age(RES, lev_d1, _)
#        CCv["{}".format(_)] = "NONE"
#        CCv["{}".format(_)] = [C, Cv]
#    return CCv


def run_epi_prep(E, M, NPI_S, pvax_d2):
    RES = run_epi(E, M, NPI_S, pvax_d2)
    res = prep_RES_all(RES)
    # res_age = prep_RES_age(RES)
    if model_type == "Gab":
        res_Gab_ses = prep_RES_ses(RES)
        return res, res_Gab_ses  # , res_age
    elif model_type == "Cij":
        return res  # , res_age


for E, M, NPI_S in SCENARIOS[model_type]:
    print("Scenario:{}-{}-{}".format(E, M, NPI_S))

    RES = [run_epi_prep(E, M, NPI_S, pvax_d2) for pvax_d2 in tqdm(PRIORITY_TO_RUN)]
    if model_type == "Gab":
        res, res_Gab_ses = np.transpose(RES)  # , res_age
        with open(
            "./res/PREP_all_{}_{}_{}_{}.pkl".format(model_type, E, M, NPI_S), "wb"
        ) as f:
            pickle.dump(res, f)
        with open(
            "./res/PREP_ses_{}_{}_{}_{}.pkl".format(model_type, E, M, NPI_S), "wb"
        ) as f:
            pickle.dump(res_Gab_ses, f)
        # with open(
        #    "./res/PREP_age_{}_{}_{}_{}.pkl".format(model_type, E, M, NPI_S), "wb"
        # ) as f:
        #    pickle.dump(res_age, f)

    elif model_type == "Cij":
        with open(
            "./res/PREP_all_{}_{}_{}_{}.pkl".format(model_type, E, M, NPI_S), "wb"
        ) as f:
            pickle.dump(RES, f)
