from utils.init_MatrixPop import (
    init_MatrixData,
    init_MatrixAge,
    init_PopSy_v2,
    init_MatrixSy_v2,
)
from utils.init_EpiModel import init_IFR, init_Omega


def init_model_params(
    model_type,
    NPI_S,
    config_mat,
    config_epi,
    config_mat_var,
    config_epi_var,
    sens_params,
    priority_dim2,
):
    # . init NPI
    if NPI_S != "":  # True:
        config_npi = sens_params["npi_params"][NPI_S[-1]]
        config_delta_C = sens_params["npi_params"]["delta_C"]
        T_npi = config_npi["day_npi"]
        NPI_bool = True
    else:
        NPI_bool = False
        M_npi = None
        T_npi = None

    # . init PoP and Matrix
    if model_type == "Cij":
        # . matrix
        N_d1, M_d1, M_d1_tot, lev_d1 = init_MatrixAge(
            config_mat["country"], config_mat_var["sim_type"], config_mat_var["nw"]
        )
        # . IFR
        IFR = init_IFR(lev_d1, config_epi["ifr_type"], model_type)
        args_model = [N_d1, M_d1, lev_d1]
        # . NPI
        if NPI_S != "":
            M_npi = {
                i: {
                    j: M_d1_tot[i][j] * config_delta_C[NPI_S[:3]] / N_d1[i]
                    for j in lev_d1
                }
                for i in lev_d1
            }

    elif model_type == "Gab":
        # . matrix
        if config_mat_var["sim_type"] == "synthetic":
            dist_Pd2 = config_mat_var["dist_Pd2"]

            N_d1, N_d2, N_d1d2, lev_d1, lev_d2, lev_d1d2, M_d1_tot = init_PopSy_v2(
                config_mat["country"], config_mat["n_dim2"], dist_Pd2
            )

            if config_mat_var["mixing_type"] == "assortative_mixing":
                args_sy_matrix = [
                    config_mat_var["r"],
                    config_mat_var["activity"],
                    config_mat_var["ds"],
                ]
            elif config_mat_var["mixing_type"] == "random_mixing":
                args_sy_matrix = [dist_Pd2]

            M_d1d2 = init_MatrixSy_v2(
                M_d1_tot,
                N_d1d2,
                config_mat["n_dim2"],
                config_mat_var["mixing_type"],
                args_sy_matrix,
            )

            # . NPI
            if NPI_S != "":
                M_d1_tot_NPI = {
                    i: {j: M_d1_tot[i][j] * config_delta_C[NPI_S[:3]] for j in lev_d1}
                    for i in lev_d1
                }
                if config_npi["mixing_type"] == "assortative_mixing":
                    args_sy_matrix_npi = [
                        config_npi["r"],
                        config_npi["activity"],
                        config_npi["ds"],
                    ]
                elif config_npi["mixing_type"] == "random_mixing":
                    args_sy_matrix_npi = [config_mat_var["dist_Pd2"]]
                M_npi = init_MatrixSy_v2(
                    M_d1_tot_NPI,
                    N_d1d2,
                    config_mat["n_dim2"],
                    config_npi["mixing_type"],
                    args_sy_matrix_npi,
                )
        elif config_mat_var["sim_type"] == "mazsk":
            N_d1d2, N_d1, N_d2, M_d1d2, lev_d1, lev_d2, lev_d1d2 = init_MatrixData(
                config_mat_var["nw"]
            )
        else:
            raise ValueError("Invalid sim_type")
        # . IFR
        IFR = init_IFR(lev_d1d2, config_epi["ifr_type"], model_type)
        args_model = [N_d1d2, N_d1, M_d1d2, lev_d1, lev_d2, priority_dim2]

    # . Vaccination
    vaccination_type = config_epi_var["vaccination_type"]

    if vaccination_type == "VDE":
        Omega_t = init_Omega(
            config_epi_var["vax_scenario"], config_mat_var["nw"], config_mat["country"]
        )
        args_vax = [Omega_t]
    elif vaccination_type == "NONE":
        args_vax = []

    model_args_len = len(args_model)
    vax_args_len = len(args_vax)

    args = args_model + args_vax

    return model_args_len, vax_args_len, args, IFR, NPI_bool, M_npi, T_npi
