import numpy as np
import random
import pandas as pd
import scipy.linalg as la


class EpiModel:
    def __init__(
        self,
        model_type: str,
        vaccination_type: str,
        vax_scenario: str,
        vax_saturation: float,
        IFR,
        mu_val,
        eps_val,
        R0,
        seeds,
        VE,
        g1,
        g2,
        stop,
        t0_epi,
        current_immunity,
        model_args_len,
        vax_args_len,
        *args,
        NPI=False,
        M_npi=None,
        T_npi=None,
        mode=None,
    ):
        self.model_type = model_type
        self.vaccination_type = vaccination_type
        self.vax_scenario = vax_scenario
        self.vax_saturation = vax_saturation
        self.IFR = IFR
        self.mu_val = mu_val
        self.eps_val = eps_val
        self.R0 = R0
        self.seeds = seeds
        self.VE = VE
        self.g1 = g1
        self.g2 = g2
        self.stop = stop
        self.t0_epi = t0_epi
        self.current_immunity = current_immunity
        self.NPI = NPI
        self.M_npi = M_npi
        self.T_npi = T_npi
        self.tot_Nv = 0
        self.exceeding_vax_key = {}
        self.mode = mode

        model_args = args[:model_args_len]
        vax_args = args[model_args_len : model_args_len + vax_args_len]

        if self.model_type == "Cij":
            self.initialize_Cij_model(*model_args)
        elif self.model_type == "Gab":
            self.initialize_Gab_model(*model_args)


        if self.vaccination_type == "VDE":
            self.initialize_VDE_model(*vax_args)
        else:
            raise ValueError("vaccination_type not valid")

        self.keys = list(self.N_lowerdim.keys())
        if self.keys != list(self.M.keys()):
            raise ValueError(
                "Population vector and Contact matrix must have the same keys!"
            )
        self.tot_N = sum([self.N_lowerdim[key] for key in self.keys])

        self.get_beta()

        if self.NPI is True:
            self.T_npi = T_npi + self.t0_epi

    def initialize_Cij_model(self, N_d1, M, lev_d1):
        self.N_lowerdim = N_d1
        self.M = M
        self.lev_d1 = lev_d1
        self.vax_saturation_d1 = {
            key: N_d1[key] * self.vax_saturation for key in N_d1.keys()
        }

    def initialize_Gab_model(self, N_d1d2, N_d1, M, lev_d1, lev_d2, priority_dim2):
        self.N_lowerdim = N_d1d2
        self.N_d1 = N_d1
        self.M = M
        self.lev_d1 = lev_d1
        self.lev_d2 = lev_d2
        self.priority_dim2 = priority_dim2
        self.vax_saturation_d1d2 = {
            key: N_d1d2[key] * self.vax_saturation for key in N_d1d2.keys()
        }


    def initialize_VDE_model(self, Omega_t):
        self.Omega_t = Omega_t


    # NB IS THE RATE IS EQUAL BY AGE: JUST A NUMBER
    # IF NOT: IT IS A VECTOR.

    def get_beta(self):
        TAB = pd.DataFrame(self.M)
        eigvals, eigvecs = la.eig(TAB)
        leading = eigvals.real[0]
        self.beta_val = self.mu_val * self.R0 / (leading)
        return

    # Initialize compartments
    def init_compartments(self):
        key_seed = random.choice(self.keys)
        compartments = [
            "S",
            "E",
            "I",
            "R",
            "D",
            "N",
            "E_new",
            "I_new",
            "D_new",
            "Sv",
            "Ev",
            "Iv",
            "Rv",
            "Dv",
            "Nv",
            "Ev_new",
            "Iv_new",
            "Dv_new",
            "S_vax",
            "E_vax",
            "R_vax",
        ]

        # Initialize dictionaries
        for compartment in compartments:
            setattr(self, compartment, {})

        for key in self.keys:
            for compartment in compartments:
                getattr(self, compartment)[key] = {0: 0}

            if key == key_seed:
                Ii = self.seeds
            else:
                Ii = 0

            self.I[key][0] = Ii
            self.I_new[key][0] = Ii

            self.N[key][0] = self.N_lowerdim[key]

        if self.current_immunity == 0:
            for key in self.keys:
                self.S[key][0] = self.N_lowerdim[key] - Ii

        elif self.current_immunity > 0:
            tot_R0 = self.tot_N * self.current_immunity
            for key in self.keys:
                self.R[key][0] = tot_R0 * self.N_lowerdim[key] / self.tot_N
                self.S[key][0] = (
                    self.N_lowerdim[key]
                    - Ii
                    - (tot_R0 * self.N_lowerdim[key] / self.tot_N)
                )

        return

    def trasmission_dynamic(self, t, m_key):
        if t > self.t0_epi:
            C_key = self.M[m_key]
            cc = []
            for kc in C_key:
                cc.append(
                    C_key[kc]
                    * (
                        (self.I[kc][t - 1] + ((1 - self.VE) * self.Iv[kc][t - 1]))
                        / self.N_lowerdim[kc]
                    )
                )
            lambda_val = self.beta_val * sum(cc)

            new_E = np.random.binomial(self.S[m_key][t - 1], 1.0 - np.exp(-lambda_val))
            new_Ev = np.random.binomial(
                self.Sv[m_key][t - 1], 1.0 - np.exp(-(1 - self.g1) * lambda_val)
            )

            new_I = np.random.binomial(self.E[m_key][t - 1], self.eps_val)
            new_Iv = np.random.binomial(self.Ev[m_key][t - 1], self.eps_val)

            new_RD = np.random.binomial(self.I[m_key][t - 1], self.mu_val)
            new_RDv = np.random.binomial(self.Iv[m_key][t - 1], self.mu_val)

            new_R, new_D = np.random.multinomial(
                new_RD, [(1 - self.IFR[m_key]), self.IFR[m_key]]
            )
            new_Rv, new_Dv = np.random.multinomial(
                new_RDv,
                [
                    (1 - ((1 - self.g2) * self.IFR[m_key])),
                    (1 - self.g2) * self.IFR[m_key],
                ],
            )
        else:
            new_E, new_Ev = 0, 0
            new_I, new_Iv = 0, 0
            new_R, new_Rv = 0, 0
            new_D, new_Dv = 0, 0

        return new_E, new_Ev, new_I, new_Iv, new_R, new_Rv, new_D, new_Dv

    def compute_vaccination_DE(self, t):
        if t > 948:  # 965
            if self.vax_scenario == "age_structure":
                self.Omega_t[t - 1] = {d1: 0 for d1 in self.lev_d1}
            elif self.vax_scenario == "age_random":
                self.Omega_t[t - 1] = 0

        if self.model_type == "Cij":
            # 1. random wrt population [Cij and Gab]
            if self.vax_scenario == "age_random":
                omega_key_t = {
                    key: self.N_lowerdim[key]
                    / sum(self.N_lowerdim.values())
                    * self.Omega_t[t - 1]
                    for key in self.keys
                }

            # 2. Cij and given by age
            elif self.vax_scenario == "age_structure":
                omega_key_t = self.Omega_t[t - 1]
                omega_key_t[0] = 0

        elif self.model_type == "Gab":
            dist_d2_byage = np.tile(self.priority_dim2, (len(self.lev_d1), 1))

            # 3. Gab and given by age and ses
            if self.vax_scenario == "age_structure":
                omega_key_t = {(0, d2): 0 for d2 in self.lev_d2}
                for age in self.Omega_t[t - 1].keys():
                    for d2 in self.lev_d2:
                        omega_key_t[(age, d2)] = (
                            self.Omega_t[t - 1][age] * dist_d2_byage[age][d2]
                        )

            # 4. Gab and random by age and given by ses
            elif self.vax_scenario == "age_random":
                omega_key_t = {}
                pN_d1 = {
                    age: self.N_d1[age] / sum(self.N_d1.values()) for age in self.lev_d1
                }
                for age in self.lev_d1:
                    for d2 in self.lev_d2:
                        omega_key_t[(age, d2)] = (
                            self.Omega_t[t - 1] * pN_d1[age] * dist_d2_byage[age][d2]
                        )

        omega_key_t = {key: np.round(value) for key, value in omega_key_t.items()}
        return omega_key_t

    def compute_vaccination_BE(self):
        tot_vax = self.pN_vax * self.tot_N

        if self.model_type == "Cij":
            # 1. random wrt population [Cij and Gab]
            if self.vax_scenario == "age_random":
                omega_key_t = {
                    key: self.N_lowerdim[key] / sum(self.N_lowerdim.values()) * tot_vax
                    for key in self.keys
                }

            # 2. Cij and given by age
            elif self.vax_scenario == "age_structure":
                omega_key_t = {
                    key: self.priority_dim1[key] * tot_vax for key in self.keys
                }

        elif self.model_type == "Gab":
            dist_d2_byage = np.tile(self.priority_dim2, (len(self.lev_d1), 1))

            # 3. Gab and given by age and ses
            if self.vax_scenario == "age_structure":
                omega_key_t = {(0, d2): 0 for d2 in self.lev_d2}
                for age in self.priority_dim1.keys():
                    for d2 in self.lev_d2:
                        omega_key_t[(age, d2)] = (
                            tot_vax * self.priority_dim1[age] * dist_d2_byage[age][d2]
                        )

            # 4. Gab and random by age and given by ses
            elif self.vax_scenario == "age_random":
                omega_key_t = {}
                pN_d1 = {
                    age: self.N_d1[age] / sum(self.N_d1.values()) for age in self.lev_d1
                }
                for age in self.lev_d1:
                    for d2 in self.lev_d2:
                        omega_key_t[(age, d2)] = (
                            tot_vax * pN_d1[age] * dist_d2_byage[age][d2]
                        )

        omega_key_t = {key: np.round(value) for key, value in omega_key_t.items()}
        return omega_key_t

    def update_compartments_trasmission(
        self, t, m_key, new_E, new_Ev, new_I, new_Iv, new_R, new_Rv, new_D, new_Dv
    ):
        self.S[m_key][t] = self.S[m_key][t - 1] - new_E
        self.Sv[m_key][t] = self.Sv[m_key][t - 1] - new_Ev

        self.E[m_key][t] = self.E[m_key][t - 1] + new_E - new_I
        self.Ev[m_key][t] = self.Ev[m_key][t - 1] + new_Ev - new_Iv

        self.I[m_key][t] = self.I[m_key][t - 1] + new_I - new_R
        self.Iv[m_key][t] = self.Iv[m_key][t - 1] + new_Iv - new_Rv

        self.R[m_key][t] = self.R[m_key][t - 1] + new_R
        self.Rv[m_key][t] = self.Rv[m_key][t - 1] + new_Rv

        self.D[m_key][t] = self.D[m_key][t - 1] + new_D
        self.Dv[m_key][t] = self.Dv[m_key][t - 1] + new_Dv

        self.E_new[m_key][t] = self.E_new[m_key][t - 1] + new_E
        self.Ev_new[m_key][t] = self.Ev_new[m_key][t - 1] + new_Ev

        self.I_new[m_key][t] = self.I_new[m_key][t - 1] + new_I
        self.Iv_new[m_key][t] = self.Iv_new[m_key][t - 1] + new_Iv

        self.D_new[m_key][t] = (
            self.D_new[m_key][t - 1] + new_D
        )  # nota bene è uguale a D
        self.Dv_new[m_key][t] = self.Dv_new[m_key][t - 1] + new_Dv

        self.N[m_key][t] = self.N[m_key][t - 1]
        self.Nv[m_key][t] = self.Nv[m_key][t - 1]

        self.S_vax[m_key][t] = self.S_vax[m_key][t - 1]
        self.E_vax[m_key][t] = self.E_vax[m_key][t - 1]
        self.R_vax[m_key][t] = self.R_vax[m_key][t - 1]

        return

    def update_compartments_vax(
        self, t, m_key, vax_in_S, vax_in_E, vax_in_R, exceeding_vax
    ):
        self.S[m_key][t] = self.S[m_key][t] - vax_in_S
        self.Sv[m_key][t] = self.Sv[m_key][t] + vax_in_S

        self.E[m_key][t] = self.E[m_key][t] - vax_in_E
        self.Ev[m_key][t] = self.Ev[m_key][t] + vax_in_E

        self.R[m_key][t] = self.R[m_key][t] - vax_in_R
        self.Rv[m_key][t] = self.Rv[m_key][t] + vax_in_R

        self.N[m_key][t] = self.N[m_key][t] - vax_in_S - vax_in_E - vax_in_R
        self.Nv[m_key][t] = self.Nv[m_key][t] + vax_in_S + vax_in_E + vax_in_R

        self.S_vax[m_key][t] = self.S_vax[m_key][t] + vax_in_S
        self.E_vax[m_key][t] = self.E_vax[m_key][t] + vax_in_E
        self.R_vax[m_key][t] = self.R_vax[m_key][t] + vax_in_R

        self.exceeding_vax_key[t][m_key] = exceeding_vax
        # if t> 10:
        #    print('i')
        return

    def distribute_vax(self, t, m_key, omega_key_t):
        SER_key = self.S[m_key][t] + self.E[m_key][t] + self.R[m_key][t]
        max_pop_to_vax = (self.vax_saturation * self.N_lowerdim[m_key]) - self.Nv[
            m_key
        ][t]
        max_vax_to_give = min(max_pop_to_vax, SER_key)

        if (max_vax_to_give > 0) and (omega_key_t > 0):
            if omega_key_t <= max_vax_to_give:
                # print('A')
                exceeding_vax = 0
                vax_in_S = np.floor(omega_key_t * self.S[m_key][t] / SER_key)
                vax_in_E = np.floor(omega_key_t * self.E[m_key][t] / SER_key)
                vax_in_R = np.floor(omega_key_t * self.R[m_key][t] / SER_key)
            elif omega_key_t > max_vax_to_give:
                # print('B')
                exceeding_vax = omega_key_t - max_vax_to_give
                vax_in_S = np.floor(max_vax_to_give * self.S[m_key][t] / SER_key)
                vax_in_E = np.floor(max_vax_to_give * self.E[m_key][t] / SER_key)
                vax_in_R = np.floor(max_vax_to_give * self.R[m_key][t] / SER_key)
        else:
            exceeding_vax = omega_key_t
            vax_in_S = 0
            vax_in_E = 0
            vax_in_R = 0

        return vax_in_S, vax_in_E, vax_in_R, exceeding_vax

    def compute_exceed_vax(self, t, SER_nv):
        exceed_omega_d1d2_t = {}
        if self.model_type == "Cij":
            for key in self.keys:
                if self.exceeding_vax_key[t][key] > 0:
                    exceed_omega_d1d2_t[key] = (
                        self.exceeding_vax_key[t][key]
                        * SER_nv[key]
                        / sum(SER_nv.values())
                    )
                else:
                    exceed_omega_d1d2_t[key] = 0

        elif self.model_type == "Gab":
            exceeding_vax_age = {
                d1: sum(
                    [
                        self.exceeding_vax_key[t][d1, d2]
                        for d2 in self.lev_d2
                        if (d1, d2) in self.exceeding_vax_key[t]
                    ]
                )
                for d1 in self.lev_d1
            }
            SER_nv_d1 = {
                d1: sum([SER_nv[d1, d2] for d2 in self.lev_d2 if (d1, d2) in SER_nv])
                for d1 in exceeding_vax_age.keys()
            }
            for key in SER_nv.keys():
                d1, d2 = key
                if exceeding_vax_age[d1] > 0:
                    exceed_omega_d1d2_t[(d1, d2)] = (
                        exceeding_vax_age[d1] * SER_nv[d1, d2] / SER_nv_d1[d1]
                    )
                else:
                    exceed_omega_d1d2_t[(d1, d2)] = 0
        return exceed_omega_d1d2_t

    def vaccination_dynamic(self, t, omega_key_t):
        # . vaccination
        self.exceeding_vax_key.setdefault(t, {})

        for m_key in self.keys:
            vax_in_C = self.distribute_vax(t, m_key, omega_key_t[m_key])
            self.update_compartments_vax(t, m_key, *vax_in_C)

        tot_exiding_vax = sum(
            self.exceeding_vax_key[t][_] for _ in self.exceeding_vax_key[t]
        )
        while tot_exiding_vax > 0:
            SER_nv = {_: self.S[_][t] + self.E[_][t] + self.R[_][t] for _ in self.keys}
            exceed_omega_t = self.compute_exceed_vax(t, SER_nv)
            for m_key in self.keys:
                vax_in_C = self.distribute_vax(t, m_key, exceed_omega_t[m_key])
                self.update_compartments_vax(t, m_key, *vax_in_C)
            tot_exiding_vax = sum(
                self.exceeding_vax_key[t][_] for _ in self.exceeding_vax_key[t]
            )
        return

    def compute_Omega_t_costant(self, t):
        self.Omega_t = {}
        SER_nv = {_: self.S[_][t] + self.E[_][t] + self.R[_][t] for _ in self.keys}
        SER_nv_tot = sum(SER_nv.values())
        omega_t = self.vax_rate * SER_nv_tot

        if self.vax_scenario == "age_random":
            self.Omega_t = {t - 1: omega_t}
        elif self.vax_scenario == "age_structure":
            self.Omega_t = {
                t - 1: {d1: omega_t * self.priority_dim1[d1] for d1 in self.lev_d1}
            }
        return

    def run(self):
        # print(self.beta_val)
        self.init_compartments()

        # vaccination before epi
        if self.vaccination_type == "VBE":
            omega_key_t = self.compute_vaccination_BE()
            self.vaccination_dynamic(0, omega_key_t)

        for t in range(1, self.stop + 1):
            # print(t, self.M[(0,0)][0, 0],self.M[(1,0)][1, 0])
            # .NPI
            if self.NPI is True:
                # new condition for NPI based on infected cases
                if t >= self.T_npi:
                    self.M = self.M_npi.copy()
            # .trasmission
            for m_key in self.keys:
                new_C = self.trasmission_dynamic(t, m_key)
                self.update_compartments_trasmission(t, m_key, *new_C)
            # vaccination during epi
            if self.vaccination_type == "VDE":
                omega_key_t = self.compute_vaccination_DE(t)
                self.vaccination_dynamic(t, omega_key_t)

            # vaccination during epi constant
            elif self.vaccination_type == "VDE_constant":
                self.compute_Omega_t_costant(t)
                omega_key_t = self.compute_vaccination_DE(t)
                self.vaccination_dynamic(t, omega_key_t)

        if self.mode == "sensitivity":
            if self.model_type == "Gab":
                dims = (len(self.lev_d1), len(self.lev_d2))
            elif self.model_type == "Cij":
                dims = len(self.lev_d1)

            # arrays = ['DnvT', 'DvT', 'RnvT', 'RvT', 'I_newT', 'Iv_newT', 'D_newT', 'Dv_newT', 'NvT']
            arrays = [
                "D",
                "Dv",
                "R",
                "Rv",
                "E_new",
                "Ev_new",
                "I_new",
                "Iv_new",
                "N",
                "Nv",
                "S_vax",
                "E_vax",
                "R_vax",
            ]
            data_arrays = {name + "T": np.zeros(dims) for name in arrays}

            for key in self.keys:
                if self.model_type == "Gab":
                    idx = (key[0], key[1])
                elif self.model_type == "Cij":
                    idx = key
                for name in arrays:
                    data_arrays[name + "T"][idx] = getattr(self, name)[key][self.stop]
            return data_arrays

        else:
            RES = [
                self.S,
                self.Sv,
                self.E,
                self.Ev,
                self.I,
                self.Iv,
                self.R,
                self.Rv,
                self.D,
                self.Dv,
                self.E_new,
                self.Ev_new,
                self.I_new,
                self.Iv_new,
                self.D_new,
                self.Dv_new,
                self.N,
                self.Nv,
                t,
                self.exceeding_vax_key,
                self.S_vax,
                self.E_vax,
                self.R_vax,
            ]
            return RES

    def run_stoch(self, N_iter):
        return [self.run() for _ in range(N_iter)]
