import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np

from utils.fig_tools import plot_combined_pvax_legend, thousands_formatter
from utils.data_tools import load_RES_Gab, load_N, load_RATES, prep_RATES_all
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

x_all = np.linspace(-0.20, 0.20, len(lev_d2))
x_ses = np.array([x_all + p for p in range(len(pvax_d2_l.keys()))])


COMPS = ["I", "D_new"]
RATE_TYPES = ["AR", "MR"]

country = "Hungary"

lab = ["Prevalence", "Deaths"]
lab0 = ["Cases per 1K", "Deaths per 10K"]

mol = [1000, 10000]


def make_figure_2(parameters):
    fig, axs = plt.subplot_mosaic(
        [["a", "b"], ["x", "x"]],
        figsize=(9.7, 3),
        height_ratios=[1, 0.3],
        width_ratios=[1, 1],
        gridspec_kw={"wspace": 0.4, "hspace": 0.7},
    )

    for param, axs_p in zip(parameters, [("a"), ("b")]):
        E, M, NPI_S = param
        print(f"Processing {E}, {M}, {NPI_S}")
        RES_Gab = load_RES_Gab(E, M, NPI_S, "all")
        RES = load_RATES(E, M, NPI_S, "Gab")
        config_mat_var = params["matrix_params"][M]
        N_d1d2, N_d1, N_d2 = load_N(config_mat_var, config_mat, country)

        if NPI_S == "":
            activity = config_mat_var["activity"]
        else:
            activity = params["npi_params"][NPI_S[-1]]["activity"]

        for r in [0]:
            RATE_TYPE = RATE_TYPES[r]
            COMP = COMPS[r]

            ax_curves = axs[axs_p[r]].inset_axes([0.0, 0, 0.615, 1])
            ax_rates = axs[axs_p[r]].inset_axes([0.76, 0, 0.30, 1])
            axs[axs_p[r]].axis("off")

            for ip, p in enumerate(
                p_vax_positions
            ):  
                Cv = RES_Gab[p][COMP][1].T
                Cnv = RES_Gab[p][COMP][0].T
                C = Cv + Cnv

                ax_curves.plot(Cv.index, C["50%"], color=pvax_colors[ip], ls="-", lw=1.8)
                ax_curves.fill_between(
                    Cv.index, C["25%"], C["75%"], color=pvax_colors[ip], alpha=0.2
                )
                ax_curves.yaxis.set_major_formatter(
                    ticker.FuncFormatter(thousands_formatter)
                )
                ax_curves.xaxis.set_major_locator(ticker.MaxNLocator(4))

                if NPI_S != "":
                    if E.startswith("E2_S1"):
                        ax_curves.axvline(x=65, color="darkgray", linestyle="--", lw=0.5)
                        ax_curves.text(
                            0.22,
                            0.3,
                            "NPI introduced",
                            transform=ax_curves.transAxes,
                            fontsize=6.5,
                            color="darkgrey",
                            rotation=90,
                            ha="right",
                        )
                    else:
                        ax_curves.axvline(
                            x=65 + 440, color="darkgray", linestyle="--", lw=0.5
                        )
                        ax_curves.text(
                            0.1,
                            0.3,
                            "NPI released",
                            transform=ax_curves.transAxes,
                            fontsize=6.5,
                            color="darkgray",
                            rotation=90,
                            ha="right",
                        )

                # AR
                RATE_v, RATE_nv, RATE = prep_RATES_all(RES, p, RATE_TYPE)
                Rate_all_25, Rate_all, Rate_all_75 = RATE

                Rate_all_25 = Rate_all_25 * mol[r]
                Rate_all = Rate_all * mol[r]
                Rate_all_75 = Rate_all_75 * mol[r]

                cap_line = Line2D(
                    [ip, ip],
                    [Rate_all_25, Rate_all_75],
                    color=pvax_colors[ip],
                    linewidth=7.1,
                    alpha=0.4,
                    solid_capstyle="round",
                )
                ax_rates.add_line(cap_line)
                ax_rates.scatter(ip, Rate_all, color=pvax_colors[ip], s=50, marker="_")
                ax_rates.set_ylabel(lab0[r], fontsize=9, labelpad=0)

                ax_rates.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
                ax_rates.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

            ax_rates.set_xlim(-0.5, 3.5)
            ax_rates.set_xticks(range(len(pvax_d2_l.keys())))
            ax_rates.set_xticklabels(["VD1", "VD2", "VD3", "VD4"], fontsize=8, rotation=90)

            ax_curves.set_ylabel(lab[r], fontsize=9)
            ax_curves.set_xlabel("Time", fontsize=9, labelpad=1)

            if E.startswith("E2_S1"):
                ax_curves.set_xlim(0, 300)
            else:
                ax_curves.set_xlim(450, 1000)

            for ax in [ax_curves, ax_rates]:
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.tick_params(axis="y", labelsize=8.5, pad=0)
                ax.tick_params(axis="x", labelsize=8.5, pad=0)

        ax_activity = axs[axs_p[0]].inset_axes([0.45, 0.75, 0.12, 0.3])

        ax_activity.spines["top"].set_visible(False)
        ax_activity.spines["right"].set_visible(False)
        ax_activity.bar(range(len(activity)), activity, color="lightgray")

        ax_activity.tick_params(axis="y", labelsize=6, pad=0)
        ax_activity.tick_params(axis="x", labelsize=6, pad=0.5)
        ax_activity.set_xticks(range(len(activity)))
        ax_activity.set_ylabel("activity", fontsize=6, labelpad=0)
        ax_activity.set_xlabel("dim 2", fontsize=6, labelpad=0)


    for ax in axs:
        if ax == "a":
            axs[ax].text(
                0.13,
                1.37,
                "a. Baseline",
                transform=axs[ax].transAxes,
                fontsize=11,
                fontweight="bold",
                va="top",
                ha="right",
            )
        elif ax == "b":
            title = "b. Tightened NPI" if E.startswith("E2_S1") else "b. Relaxed NPI"
            axs[ax].text(
                0.29,
                1.37,
                title,
                transform=axs[ax].transAxes,
                fontsize=11,
                fontweight="bold",
                va="top",
                ha="right",
            )

        if ax != "x":
            axs[ax].text(
                -0.02,
                1.18,
                "1.",
                transform=axs[ax].transAxes,
                fontsize=9,
                va="top",
                ha="right",
            )
            axs[ax].text(
                0.73,
                1.18,
                "2.",
                transform=axs[ax].transAxes,
                fontsize=9,
                va="top",
                ha="right",
            )


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
        }
    ]

    plot_combined_pvax_legend(axs["x"], datasets)
    proportion = 0.35  # Example proportion for resizing width
    box = axs["x"].get_position()
    new_width = box.width * proportion
    adjustment = (box.width - new_width) / 2
    new_box = [box.x0 + adjustment, box.y0, new_width, box.height]
    axs["x"].set_position(new_box)

    os.makedirs("figs", exist_ok=True)
    plt.savefig(f"./figs/figure_{E}_{M}_{NPI_S}.pdf", bbox_inches="tight")
    plt.close()
    return


params_to_run = [
    [("E2_S1_b", "M7", ""), ("E2_S1_b", "M7", "NPI_4")],  # figure 2
    [('E2_S2_b','M7', ''),  ('E2_S2_b','M7', 'rel_4')]   # figure 4
]
for param in params_to_run:
    make_figure_2(param)