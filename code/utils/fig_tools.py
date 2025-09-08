import numpy as np

def plot_combined_pvax_legend(ax, datasets):
    """
    Plot multiple datasets as grouped bar charts in a single axis.

    Parameters:
    -----------
    ax : matplotlib axis
        The axis to plot on
    datasets : list of dict
        Each dict should contain:
        - 'data': dict with categories and values (like pvax_d2_l or population)
        - 'letters': list of letters for bar labels
        - 'colors': list of colors for bars
        - 'labels': list of x-axis labels (optional, uses data keys if not provided)
    """

    bar_width = 0.18
    bar_spacing = 0.015  # White space between bars
    current_x_pos = 0
    all_x_positions = []
    all_labels = []

    for dataset_idx, dataset in enumerate(datasets):
        data = dataset["data"]
        letters = dataset["letters"]
        colors = dataset["colors"]
        labels = dataset.get("labels", list(data.keys()))

        # Add extra space before second dataset (population)
        if dataset_idx > 0:
            current_x_pos += 0.5  # Extra spacing between VDs and population

        for i, (category, values) in enumerate(data.items()):
            # Position bars within each category
            n_bars = len(values)
            if n_bars > 1:
                # Add spacing between bars
                bar_positions = np.array(range(n_bars)) * (bar_width + bar_spacing)
                bar_positions = bar_positions - (bar_positions.mean())
            else:
                bar_positions = [0]

            # Offset by current group position
            bar_positions = np.array(bar_positions) + current_x_pos
            all_x_positions.append(current_x_pos)

            # Plot bars with shadow colors
            if isinstance(colors[i], list):
                # If colors[i] is a list of shadow colors, use them for each bar
                bar_colors = colors[i][
                    : len(values)
                ]  # Take only as many colors as bars
            else:
                # If colors[i] is a single color, use it for all bars
                bar_colors = [colors[i]] * len(values)

            bars = ax.bar(
                bar_positions, values, width=bar_width, color=bar_colors, alpha=1
            )

            # Add letters on each bar
            for j, (bar, value) in enumerate(zip(bars, values)):
                if j < len(letters):
                    letter = letters[j]
                else:
                    letter = letters[
                        j % len(letters)
                    ]  # Cycle through letters if not enough
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value / 2,
                    letter,
                    ha="center",
                    va="center",
                    color="w",
                    fontsize=7,
                    fontweight="bold",
                )

            current_x_pos += 1.0  # Regular spacing between categories

        all_labels.extend(labels)

    # Set x-axis properties
    ax.set_xticks(all_x_positions)
    ax.set_xticklabels(all_labels, fontsize=8)

    # Remove spines and y-axis
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])

    return ax


def thousands_formatter(x, pos):
    if x == 0:
        return "0"
    else:
        return "{}K".format(int(x / 1000))
