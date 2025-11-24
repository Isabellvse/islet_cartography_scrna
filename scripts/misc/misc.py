import re
import matplotlib            # Base matplotlib functionality
import seaborn as sns
import matplotlib.pyplot as plt

def to_snake_case(name):
    # Convert to lowercase, replace spaces and parentheses with underscores, remove multiple underscores
    name = name.lower()
    name = re.sub(r"[^\w]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")

def set_my_theme():
    sns.set_style("white")  # similar to theme_classic background

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    plt.rcParams.update({
        # text + labels
        "axes.titlesize": 8,
        "axes.titleweight": "normal",
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 4,
        "legend.title_fontsize": 4,

        # line & tick colors
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "text.color": "black",

        # spine & tick thickness
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,

        # IMPORTANT: make ticks long & visible
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.bottom": True,
        "ytick.left": True,

        # background
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })