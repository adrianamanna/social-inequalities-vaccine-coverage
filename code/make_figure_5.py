import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np

from utils.data_tools import (
    prep_RATES_all,
    load_N,
    load_RATES,
    prep_RATES_d2,
    rates_d1d2_inferred,
)
from utils.fig_tools import plot_combined_pvax_legend
from utils.LABS import *
from utils.tools import upload_yaml

plt.rcParams.update(
    {"font.size": 9, "font.style": "normal", "font.family": "sans-serif"}
)

pvax_colors = ["#0081a7", "#6dc4bc", "#DC6A41", "#fdc856"]

pvax_colors_shadows = {
    0: ["#1a5f7a", "#0081a7", "#4da6c7"],  # Dark blue shadows
    1: ["#3a8b85", "#6dc4bc", "#a1ddd6"],  # Mint/teal shadows
    2: ["#b54427", "#DC6A41", "#e8956f"],  # Orange/rust shadows
    3: ["#d9a52d", "#fdc856", "#fdd985"],  # Yellow/gold shadows
}

lev_d1 = [0, 1, 2, 3, 4, 5, 6, 7]
lev_d2 = [0, 1, 2]


config_mat = upload_yaml("config_matrix")
config_epi = upload_yaml("config_epi")
params = upload_yaml("parameters")

# .bias dist2
PRIORITIES_D2 = config_epi["PRIORITIES_D2"]
P_d2_l = PRIORITIES_D2["Gab"]

p_vax_to_plot = ["p5", "p0", "p3", "p8"]
p_vax_positions = [list(P_d2_l.keys()).index(key) for key in p_vax_to_plot]

pvax_d2_l = {p: P_d2_l[p] for p in p_vax_to_plot}


NPI_S = ""

mol = [1000, 10000]
lab = ["Cases per 1K", "Deaths per 10K"]


x_all = np.linspace(-0.20, 0.20, len(lev_d2))
x_ses = np.array([x_all + p for p in range(len(pvax_d2_l.keys()))])

pos = [["a1", "a2", "a3", "a4"], ["b1", "b2", "b3", "b4"]]
title_ax = ["1. Non vaccinated", "2. Vaccinated", "3. All", "4."]

fig, axs = plt.subplot_mosaic(
    [
        ["a1", "a2", "a3", "a", "a4"],
        ["b1", "b2", "b3", "b", "b4"],
        ["x", "x", "x", "x", "x"],
    ],
    figsize=(9, 5.5),
    height_ratios=[1, 1, 0.3],
    width_ratios=[1, 1, 1, 0.1, 1],
    gridspec_kw={"wspace": 0.35, "hspace": 0.8},
)

axs["a1"].set_title("a  Scenario 1", loc="left", pad=15, fontsize=11, weight="bold")
axs["b1"].set_title("b  Scenario 2", loc="left", pad=15, fontsize=11, weight="bold")

ax_span = [[20, 4], [13, 1]]

for c, params_figure5 in enumerate([("E2_2w", "Ma"), ("E2_4w", "Mb")]):
    ax_nv = axs[pos[c][0]]
    ax_v = axs[pos[c][1]]
    ax_all = axs[pos[c][2]]
    ax_diff = axs[pos[c][3]]

    E, M = params_figure5
    RES = load_RATES(E, M, NPI_S, "Gab")
    RES_Cij = load_RATES(E, M, "", "Cij")

    config_mat_var = params["matrix_params"][M]
    N_d1d2, N_d1, N_d2 = load_N(config_mat_var, config_mat, "Hungary")

    for r, RATE_TYPE in enumerate(["AR"]):
        ax_nv.set_ylabel(lab[r], fontsize=9)

        for ip, p in enumerate(
            p_vax_positions
        ):  # pvax_d2_l#[0,2,3] # NOTA BENE: nella list gli indici degli p nel dizionario che vogliamo plottare
            RATE_v, RATE_nv, RATE = prep_RATES_all(RES, p, RATE_TYPE)
            Cnv25, Cnv, Cnv75 = RATE_nv
            Cv25, Cv, Cv75 = RATE_v
            C25, C, C75 = RATE

            Cnv, Cnv25, Cnv75 = Cnv * mol[r], Cnv25 * mol[r], Cnv75 * mol[r]
            Cv, Cv25, Cv75 = Cv * mol[r], Cv25 * mol[r], Cv75 * mol[r]
            C, C25, C75 = C * mol[r], C25 * mol[r], C75 * mol[r]

            cap_line = Line2D(
                [ip, ip],
                [Cnv25, Cnv75],
                color=pvax_colors[ip],
                linewidth=7.1,
                alpha=0.4,
                solid_capstyle="round",
            )
            ax_nv.add_line(cap_line)
            ax_nv.scatter(
                ip, Cnv, color=pvax_colors[ip], s=50, marker="_", linewidths=2
            )

            cap_line = Line2D(
                [ip, ip],
                [Cv25, Cv75],
                color=pvax_colors[ip],
                linewidth=7.1,
                alpha=0.4,
                solid_capstyle="round",
            )
            ax_v.add_line(cap_line)
            ax_v.scatter(ip, Cv, color=pvax_colors[ip], s=50, marker="_", linewidths=2)

            cap_line = Line2D(
                [ip, ip],
                [C25, C75],
                color=pvax_colors[ip],
                linewidth=7.1,
                alpha=0.4,
                solid_capstyle="round",
            )
            ax_all.add_line(cap_line)
            ax_all.scatter(ip, C, color=pvax_colors[ip], s=50, marker="_", linewidths=2)

            # plotting difference on last ax
            _, _, RATE = prep_RATES_d2(RES, p, "AR")
            RATE_d2_25, RATE_d2, RATE_d2_75 = RATE
            RATE_d2_25_infered, RATE_d2_infered, RATE_d2_75_infered = (
                rates_d1d2_inferred(RES_Cij, RES, p, "AR", N_d2)
            )  # qui cambio metodo #rates_d2_inferred

            for i_d2, d2 in enumerate(lev_d2):
                # .diff between Gab and Cij
                true50, true75, true25 = RATE_d2[d2], RATE_d2_75[d2], RATE_d2_25[d2]
                inf50, inf75, inf25 = (
                    RATE_d2_infered[d2],
                    RATE_d2_75_infered[d2],
                    RATE_d2_25_infered[d2],
                )

                diff50 = RATE_d2[d2] - RATE_d2_infered[d2]

                diff75 = (true50 - inf50) + (true75 - true50) + (inf75 - inf50)  #
                diff25 = (true50 - inf50) - (true50 - true25) + (inf50 - inf25)

                cap_line = Line2D(
                    [x_ses[ip][d2], x_ses[ip][d2]],
                    [diff25, diff75],
                    color=pvax_colors_shadows[ip][i_d2],
                    linewidth=7.1,
                    alpha=0.4,
                    solid_capstyle="projecting",
                )
                ax_diff.add_line(cap_line)
                ax_diff.scatter(
                    x_ses[ip][d2],
                    diff50,
                    color=pvax_colors_shadows[ip][i_d2],
                    s=50,
                    marker="_",
                    linewidths=2,
                )

            # setting on diff ax
            ylims = ax_diff.get_ylim()
            ax_diff.vlines(
                [0.5, 1.5, 2.5],
                ylims[0],
                ylims[1],
                color="lightgrey",
                ls="-",
                lw=1.3,
                alpha=0.3,
            )

            ax_diff.set_ylabel(r"$G_{\mathbf{a}\mathbf{b}}-f(C_{ij})$", fontsize=9)
            ax_diff.set_xlabel("dim 2", labelpad=1)
            ax_diff.set_xticks(x_ses.flatten(), ["L", "M", "H"] * 4, fontsize=7)

            if E == "E2_2w":
                ax_diff.set_ylim(-0.55, 0.55)
            elif E == "E2_4w":
                ax_diff.set_ylim(-0.15, 0.15)

        span = ax_span[c][r]
        ax_nv.set_ylim(Cnv - span, Cnv + span)
        ax_v.set_ylim(Cv - span, Cv + span)
        ax_all.set_ylim(C - span, C + span)


        for ia, ax_ in enumerate([ax_nv, ax_v, ax_all, ax_diff]):
            ax_.tick_params(axis="y", labelsize=8, pad=2)
            ax_.tick_params(axis="x", labelsize=8, pad=2)
            ax_.spines["top"].set_visible(False)
            ax_.spines["right"].set_visible(False)
            ax_.text(0.05, 0.95, title_ax[ia], fontsize=9, transform=ax_.transAxes)

        for ia, ax_ in enumerate([ax_nv, ax_v, ax_all]):
            ax_.set_xlim(-0.5, 3.5)
            ax_.set_xticks(range(len(p_vax_positions)), ["VD1", "VD2", "VD3", "VD4"])
            ax_.yaxis.set_major_locator(ticker.MaxNLocator(4))


datasets = [
    {
        "data": pvax_d2_l,
        "letters": ["L", "M", "H"],
        "colors": [
            pvax_colors_shadows[0],  # VD1 gets blue shadows [dark, medium, light]
            pvax_colors_shadows[1],  # VD2 gets mint shadows [dark, medium, light]
            pvax_colors_shadows[2],  # VD3 gets orange shadows [dark, medium, light]
            pvax_colors_shadows[3],  # VD4 gets yellow shadows [dark, medium, light]
        ],
        "labels": [r"VD1$^{*}$", r"VD2$^{*}$", r"VD3$^{*}$", r"VD4$^{*}$"],
    }
]

plot_combined_pvax_legend(axs["x"], datasets)
proportion = 0.35
box = axs["x"].get_position()
new_width = box.width * proportion
adjustment = (box.width - new_width) / 2
new_box = [box.x0 + adjustment, box.y0, new_width, box.height]
axs["x"].set_position(new_box)

axs["a"].set_visible(False)
axs["b"].set_visible(False)


os.makedirs("figs", exist_ok=True)
plt.savefig(f"./figs/figure_{E}_{M}.pdf", bbox_inches="tight")
plt.close()