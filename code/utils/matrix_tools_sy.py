import numpy as np
import pandas as pd

# SYNTETIC MATRICES
def get_Pd1_Md1(data, country):
    """
    ::returns::
    - P_d1 :: dict
        key : age groups
        values: pop fracion of age group i
    - N_pop :: int
        total number of individuals in country
    - M_d1 :: dict
        key (i): age groups
        values: dict {age group j: average number of  contacts of i with age group j}
    - M_d1_tot: as M_d1 but with total number of contacts"""

    data_c = data[data["country"] == country].reset_index(drop=True)

    # matrix by age
    M_age = data_c["M"][0]
    M_d1 = pd.DataFrame(M_age.T).to_dict()

    # population by age
    P = data_c["N"][0] / sum(data_c["N"][0])
    P_d1 = pd.DataFrame(P).to_dict()[0]

    # total population
    N_pop = data_c["tot_N_pop"][0]

    # matrix of total number of contacts
    M_d1_tot = {}
    for i in M_d1.keys():
        M_d1_tot[i] = {}
        for j in M_d1[i].keys():
            M_d1_tot[i][j] = M_d1[i][j] * P_d1[i] * N_pop

    return P_d1, N_pop, M_d1, M_d1_tot


def get_P_d1d2(P_d1, n_dim2, dist_Pd2):
    n_dim1 = len(P_d1)
    P_d1d2 = {}

    for i in range(n_dim1):
        for alpha in range(n_dim2):
            P_d1d2[i, alpha] = P_d1[i] * dist_Pd2[alpha]

    return P_d1d2


def random_mixing(dist_p2):
    pa, pb, pc = dist_p2
    m = [
        [pa * pa, pa * pb, pa * pc],
        [pb * pa, pb * pb, pb * pc],
        [pc * pa, pc * pb, pc * pc],
    ]

    m = np.array(m)
    m_diag, m_off_diag = m.copy(), m.copy()
    return m_diag, m_off_diag


# Split matrix as assortative
def phi(r, p):
    return (1 - r) * p


def assortative_mixing(r, activity, ds):
    pa, pb, pc = activity
    ra, rb, rc = r
    da, db, dc = ds

    dab = (phi(ra, pa) + phi(rb, pb) - phi(rc, pc)) / 2
    dac = (phi(ra, pa) + phi(rc, pc) - phi(rb, pb)) / 2
    dbc = (phi(rb, pb) + phi(rc, pc) - phi(ra, pa)) / 2

    m_diag = [[ra * pa, dab, dac], [dab, rb * pb, dbc], [dac, dbc, rc * pc]]

    m_off_diag = [
        [ra * pa, (1 - ra) * pa * da, (1 - ra) * pa * (1 - da)],
        [(1 - rb) * pb * (1 - db), rb * pb, (1 - rb) * pb * db],
        [(1 - rc) * pc * (1 - dc), (1 - rc) * pc * dc, rc * pc],
    ]

    m_diag, m_off_diag = np.array(m_diag), np.array(m_off_diag)
    if np.any(np.array(m_diag) < 0) or np.any(np.array(m_off_diag) < 0):
        raise ValueError(
            f"Non-physical solution:: negative values in the ASS matrices \n:{m_diag, m_off_diag}"
        )

    return m_diag, m_off_diag


def get_M_d1d2(M_d1_tot, n_dim2, N_d1d2, m_diag, m_off_diag):
    n_dim1 = len(M_d1_tot)
    M_d1d2 = {}
    M_d1d2_tot = {}

    for i in range(n_dim1):
        for alpha in range(n_dim2):
            M_d1d2.setdefault((i, alpha), {})
            M_d1d2_tot.setdefault((i, alpha), {})

            for j in range(n_dim1):
                for beta in range(n_dim2):
                    if i == j:
                        m = m_diag.copy()
                    elif i > j:
                        m = m_off_diag.copy()
                    elif i < j:
                        m = m_off_diag.copy().T

                    M_d1d2_tot[(i, alpha)][(j, beta)] = M_d1_tot[i][j] * m[alpha][beta]
                    M_d1d2[(i, alpha)][(j, beta)] = (
                        M_d1d2_tot[(i, alpha)][(j, beta)] / (N_d1d2[i, alpha])
                    )

    return M_d1d2_tot, M_d1d2

# -------------------------------------------------------



