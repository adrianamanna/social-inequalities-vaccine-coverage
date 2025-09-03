# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import numpy as np
import pickle
from tqdm import tqdm

from utils.EpiModels import EpiModel  # note model
from utils.init_EpiParams import init_model_params
from utils.tools import upload_yaml


# . init config
config_params_torun = upload_yaml("config_params_torun")
config_mat = upload_yaml("config_matrix")
config_epi = upload_yaml("config_epi")
params = upload_yaml("parameters")
country = config_mat["country"]

# R0s              = config_params_torun['R0s']
model_type = config_params_torun["model_type"]
SCENARIOS = config_params_torun["params1&3"]  # 'params1&3', params2&4_NPIs, 'params5'

# . init priorities
PRIORITIES_D2 = config_epi["PRIORITIES_D2"]
PRIORITY_TO_RUN = PRIORITIES_D2[model_type]
print(PRIORITIES_D2)


def run_epi_rates(E, M, NPI_S, pvax_d2):
    config_mat_var = params["matrix_params"][M]
    config_epi_var = params["epi_params"][E]
    r0 = config_epi_var["R0s"][0]  # only one R0 NOTA BENE
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
        mode="sensitivity",
    )

    RES = np.array(model.run_stoch(config_epi["N_iter"]), dtype=object)
    return RES


for E, M, NPI_S in SCENARIOS[model_type]:
    print("Scenario:{E}-{M}-{NPI_S}".format(E=E, M=M, NPI_S=NPI_S))
    if M in ["Mc", "MO"] and NPI_S != "":
        continue
    RES = [run_epi_rates(E, M, NPI_S, pvax_d2) for pvax_d2 in tqdm(PRIORITY_TO_RUN)]
    with open(f"./res/RATES_{model_type}_{E}_{M}_{NPI_S}.pkl", "wb") as f:
        pickle.dump(RES, f)
