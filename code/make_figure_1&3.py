import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker

import numpy as np
import os
os.makedirs("figs", exist_ok=True)

plt.rcParams.update(
    {"font.size": 9, "font.style": "normal", "font.family": "sans-serif"}
)

from utils.data_tools import (
    load_RES,
    load_N,
    load_RATES,
    prep_RATES_d2,
    rates_d1d2_inferred,
)

from utils.fig_tools import plot_combined_pvax_legend, thousands_formatter
from utils.LABS import *
from utils.tools import upload_yaml

npi_str = ""
vaccination_type = "VDE"
ifr = "age"

lev_d1 = [0, 1, 2, 3, 4, 5, 6, 7]
lev_d2 = [0, 1, 2]

pvax_colors = ["#0081a7", "#6dc4bc", "#DC6A41", "#fdc856"]
pvax_colors_shadows = {
    0: ["#1a5f7a", "#0081a7", "#4da6c7"],  # Dark blue shadows
    1: ["#3a8b85", "#6dc4bc", "#a1ddd6"],  # Mint/teal shadows
    2: ["#b54427", "#DC6A41", "#e8956f"],  # Orange/rust shadows
    3: ["#d9a52d", "#fdc856", "#fdd985"],  # Yellow/gold shadows
}


config_mat = upload_yaml("config_matrix")
config_epi = upload_yaml("config_epi")
params = upload_yaml("parameters")

P_d2_l = config_epi["PRIORITIES_D2"]["Gab"]

p_vax_to_plot = ["p5", "p0", "p3", "p8"]
p_vax_positions = [list(P_d2_l.keys()).index(key) for key in p_vax_to_plot]
print(p_vax_positions)

pvax_d2_l = {p: P_d2_l[p] for p in p_vax_to_plot}

x_all = np.linspace(-0.20, 0.20, len(lev_d2))
x_ses = np.array([x_all + p for p in range(len(pvax_d2_l.keys()))])


# -----
E, M, M_Cij = ("E2_S1_b", "M7", "C0")
S = E.split("_")[1]
print(E, M, S)
# -----

rate_type = "AR"
country = "Hungary"

pos = ["a", "b", "c", "d"]

config_mat_var = params["matrix_params"][M]
N_d1d2, N_d1, N_d2 = load_N(config_mat_var, config_mat, country)


fig, axs = plt.subplot_mosaic([["a", "b", "c", "d"], ["x", "x", "x", "x"]],
        figsize=(12, 3.5),
        height_ratios=[1, 0.25],
        width_ratios=[1, 0.8, 0.8, 0.8],
        gridspec_kw={"wspace": 0.38, "hspace": 0.39},
    )
ax_curves = axs[pos[0]]
ax_vax = axs[pos[1]]
ax_rates = axs[pos[2]]
ax_diff = axs[pos[3]]

RES_Gab_r0s, res_Cij_r0s = load_RES(E, M, M_Cij, "all")
RES_Gab_ses_r0s, _ = load_RES(E, M, M_Cij, "ses")

RES = load_RATES(E, M, "", "Gab")
RES_Cij = load_RATES(E, M_Cij, "", "Cij")

Iv = res_Cij_r0s[0][0]["I"][1].T
Inv = res_Cij_r0s[0][0]["I"][0].T

I = Iv + Inv

ax_curves.plot(Iv.index, I["50%"], color="k", label="Cij", ls=":")
ax_curves.fill_between(Iv.index, I["25%"], I["75%"], color="k", alpha=0.2)
ax_curves.yaxis.set_major_formatter(ticker.FuncFormatter(thousands_formatter))

for ip, p in enumerate(p_vax_positions): 
    Iv = RES_Gab_r0s[p]["I"][1].T
    Inv = RES_Gab_r0s[p]["I"][0].T
    I = Iv + Inv

    ax_curves.plot(
        Iv.index, I["50%"], color=pvax_colors[ip], label=f"{p}", ls="-", lw=1.8
        )
    ax_curves.fill_between(
            Iv.index, I["25%"], I["75%"], color=pvax_colors[ip], alpha=0.2
        )

    _, _, RATE = prep_RATES_d2(RES, p, "AR")
    RATE_d2_25, RATE_d2, RATE_d2_75 = RATE
    RATE_d2_25_infered, RATE_d2_infered, RATE_d2_75_infered = rates_d1d2_inferred(
    RES_Cij, RES, p, "AR", N_d2
        ) 

    for i_d2, d2 in enumerate(lev_d2):
        # proportion vaccinated in dim2
        i_vax_status = 1
        tot_vax_d2 = RES_Gab_ses_r0s[p]["N"][i_vax_status][d2].T["mean"].values
        tot_vax_d2_75 = RES_Gab_ses_r0s[p]["N"][i_vax_status][d2].T["75%"].values
        tot_vax_d2_25 = RES_Gab_ses_r0s[p]["N"][i_vax_status][d2].T["25%"].values

        i_vax_status = 0
        tot_not_vax_d2 = RES_Gab_ses_r0s[p]["N"][i_vax_status][d2].T["mean"].values
        tot_not_vax_d2_75 = (
            RES_Gab_ses_r0s[p]["N"][i_vax_status][d2].T["75%"].values
            )
        tot_not_vax_d2_25 = (
            RES_Gab_ses_r0s[p]["N"][i_vax_status][d2].T["25%"].values
            )

        percent_vax_d2 = (
                tot_vax_d2[-1] / (tot_vax_d2[-1] + tot_not_vax_d2[-1]) * 1000
            )
        percent_vax_d2_75 = (
                tot_vax_d2_75[-1] / (tot_vax_d2_75[-1] + tot_not_vax_d2_75[-1]) * 1000
            )
        percent_vax_d2_25 = (
                tot_vax_d2_25[-1] / (tot_vax_d2_25[-1] + tot_not_vax_d2_25[-1]) * 1000
            )

        cap_line = Line2D(
                [x_ses[ip][d2], x_ses[ip][d2]],
                [percent_vax_d2_25, percent_vax_d2_25],
                color=pvax_colors_shadows[ip][i_d2],
                linewidth=7.1,
                alpha=0.4,
                solid_capstyle="round",
            )
        ax_vax.add_line(cap_line)
        ax_vax.scatter(x_ses[ip][d2], percent_vax_d2, color=pvax_colors_shadows[ip][i_d2], s=50, marker='_')
 

        # Gab
        cap_line = Line2D(
                [x_ses[ip][d2], x_ses[ip][d2]],
                [RATE_d2_25[d2] * 1000, RATE_d2_75[d2] * 1000],
                color=pvax_colors_shadows[ip][i_d2],
                linewidth=7.1,
                alpha=0.4,
                solid_capstyle="round",
            )
        ax_rates.add_line(cap_line)
        ax_rates.scatter(
                x_ses[ip][d2],
                RATE_d2[d2] * 1000,
                color=pvax_colors_shadows[ip][i_d2],
                s=50,
                marker="_",
                linewidths=2,
            )

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

ax_curves.set_ylabel("Prevalence", labelpad=0)
ax_vax.set_ylabel("Vaccinated per 1K", labelpad=0)
ax_rates.set_ylabel("Cases per 1K", labelpad=0)
ax_diff.set_ylabel(r"$G_{\mathbf{a}\mathbf{b}}-f(C_{ij})$", fontsize=9)

ax_curves.set_xlabel("Time", labelpad=1)
ax_vax.set_xlabel("dim 2", labelpad=1)
ax_rates.set_xlabel("dim 2", labelpad=1)
ax_diff.set_xlabel("dim 2", labelpad=1)

ax_vax.set_xticks(x_ses.flatten(), ["1", "2", "3"] * 4)
ax_rates.set_xticks(x_ses.flatten(), ["1", "2", "3"] * 4)
ax_diff.set_xticks(x_ses.flatten(), ["1", "2", "3"] * 4)

ax_diff.set_ylim(-0.5, 0.5)

for ax in [ax_rates, ax_vax, ax_curves, ax_diff]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=8, pad=2)
    ax.tick_params(axis="x", labelsize=8, pad=3)

for ax in [ax_rates, ax_vax, ax_diff]:
    ylims = ax.get_ylim()
    ax.vlines(
            [0.5, 1.5, 2.5], ylims[0], ylims[1], color="grey", ls="-", lw=1.3, alpha=0.3
        )

    for _, t_ in zip([0.155, 0.395, 0.645, 0.895], ["VD1", "VD2", "VD3", "VD4"]):
        ax.text(
                _,
                0.93,
                t_,
                transform=ax.transAxes,
                fontsize=8,
                va="top",
                ha="right",
                color="grey",
            )

ax_rates.yaxis.set_major_locator(ticker.MaxNLocator(4))
ax_curves.yaxis.set_major_locator(ticker.MaxNLocator(5))
ax_vax.yaxis.set_major_locator(ticker.MaxNLocator(5))
ax_diff.yaxis.set_major_locator(ticker.MaxNLocator(5))
ax_curves.xaxis.set_major_locator(ticker.MaxNLocator(4))

# LEGEND
curves_gab = [
        (Line2D([0], [0], marker="", color=c, markerfacecolor="", ls="-", label=VD))
        for c, VD in zip(pvax_colors, ["VD1", "VD2", "VD3", "VD4"])
    ]
curves_cij = [
        (
            Line2D(
                [0],
                [0],
                marker="",
                color="black",
                markerfacecolor="",
                ls=":",
                label="$C_{ij}$",
            )
        )
    ]

fig.legend(
        handles=curves_gab,
        frameon=False,
        title=r"$G_{\mathbf{a}\mathbf{b}}$",
        fontsize=8,
        columnspacing=0.7,
        handlelength=1,
        handletextpad=0.7,
        labelspacing=0.1,
        bbox_to_anchor=(0.3, 0.86),
    )  # (0.70, 0.17) (0.64, 0.165)
fig.legend(
        handles=curves_cij,
        frameon=False,
        fontsize=9,
        columnspacing=0.7,
        handlelength=1.2,
        handletextpad=0.7,
        labelspacing=0.1,
        bbox_to_anchor=(0.3, 0.68),
    )  # (0.70, 0.17) (0.64, 0.165)

population = {"Population": [0.35, 0.45, 0.20], "Activity": [0.2, 0.4, 0.4]}
labels = ["pop dim2", "activity"]
datasets = [
        {
            "data": pvax_d2_l,
            "letters": ["1", "2", "3"],
            "colors": [
                pvax_colors_shadows[0],  # VD1 gets blue shadows [dark, medium, light]
                pvax_colors_shadows[1],  # VD2 gets mint shadows [dark, medium, light]
                pvax_colors_shadows[2],  # VD3 gets orange shadows [dark, medium, light]
                pvax_colors_shadows[3],  # VD4 gets yellow shadows [dark, medium, light]
            ],
            "labels": ["VD1", "VD2", "VD3", "VD4"],
        },
        {
            "data": population,
            "letters": ["1", "2", "3"],
            "colors": ["lightgrey", "lightgrey"],
            "labels": labels,
        },
    ]

plot_combined_pvax_legend(axs["x"], datasets)
proportion = 0.45  # Example proportion for resizing width
box = axs["x"].get_position()
new_width = box.width * proportion
adjustment = (box.width - new_width) / 2
new_box = [box.x0 + adjustment, box.y0, new_width, box.height]
axs["x"].set_position(new_box)

for key, letter in zip(["a", "b", "c", "d"], ["a", "b", "c", "d"]):
        axs[key].text(
            -0.06,
            1.17,
            letter,
            transform=axs[key].transAxes,
            fontsize=13,
            fontweight="bold",
            va="top",
            ha="right",
        )
plt.savefig(f"./figs/figure_{E}_{M}_{S}.pdf", bbox_inches="tight")
plt.show()