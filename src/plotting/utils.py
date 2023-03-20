import os
import sys
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pathlib as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from cycler import cycler
import seaborn as sns

import pandas as pd
import numpy as np
import scipy as sp


CURRENT_PATH    = os.path.abspath(".") + "/"
DATA_PATH       = CURRENT_PATH + "../../data/"
INPUT_PATH      = DATA_PATH + "input/"
OUTPUT_PATH     = DATA_PATH + "output/"
# FIGS_PATH       = CURRENT_PATH + "../../figs/"


def SET_SUBDIR(sub_directory=None):
    if isinstance(sub_directory, str):
        subdir = sub_directory.replace("/", "") + "/"
    else: 
        subdir = ""
    global FIGS_PATH
    FIGS_PATH = CURRENT_PATH + "../../figs/" + subdir

SET_SUBDIR()



def read_ASCII(filename, skiprows=0):
    if not filename.endswith(".txt"):
        filename += ".txt"
    try: 
        data = np.loadtxt(OUTPUT_PATH + filename, skiprows=skiprows)
    except FileNotFoundError:
        data = np.loadtxt(INPUT_PATH + filename, skiprows=skiprows)
    return np.asarray(data)


def save(filename, temp=True):
    if not filename.endswith(".pdf"):
        filename += ".pdf"
    file = FIGS_PATH + filename
    plt.savefig(file)
    if temp:
        file = FIGS_PATH + filename.replace(".pdf", ".png")
        plt.savefig(file)




def reset_matplotlib():
    plt.rcParams.update(plt.rcParamsDefault)


darkMode = False 


# COLOURS = ["#5d00ff", "#ef495e", "#ffa500", "#9989f2", "#c0c066", "#133bfc", "#cd7e26"]

# subCOLOURS = ["", "", "#FBAB18"]

# COLOURS = ['#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4']
COLOURS = ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8']
subCOLOURS = ["#DA8F84", "#83B5D2", "#CDC8EC", "w", "#F3C77C", "#B9D092", "#F7DBDC"]


# print(plt.style.library["fivethirtyeight"])


def _set_params():
    # face = plt.style.library["ggplot"]["axes.facecolor"]
    face = "#D4CCCD"
    custom_params = {
        "figure.figsize":(10,5), "figure.autolayout":True,
        "savefig.dpi":300, "savefig.bbox":"tight",
        "lines.linewidth":2.8, 
        "text.usetex":True, 
        "mathtext.fontset":"cm", "font.family":"STIXGeneral",
        "axes.labelsize":22, "axes.titlesize":22, "axes.titlelocation":"left", "axes.titleweight":"bold", 
        "axes.prop_cycle":cycler('color', COLOURS),
        # "axes.facecolor":face,
        "xtick.labelsize":16, "ytick.labelsize":16,
        "legend.fontsize":20, "legend.fancybox":True, "legend.frameon":True, "legend.shadow":True, 
        "font.size":18
        }
    
    if darkMode:
        custom_params["axes.facecolor"] = "#313332"
        custom_params["grid.color"] = "#616667"#"#BFBFBF"
        custom_params["legend.labelcolor"] = "w"

    for param in custom_params.keys():
        mpl.rcParams[param] = custom_params[param]

    sns.set_style("darkgrid", rc=custom_params)
    # sns.set_palette("hls")

_set_params()




def plot_colours():
    
    fig, ax = plt.subplots()
    x = np.linspace(0, 5*np.pi, 100)
    for i in range(len(COLOURS)):
        c = COLOURS[i].strip("#")
        ax.plot(x, np.sin(x+i*np.pi*0.3)*x-i*1.5, label=f"{i}: {c}")
    ax.legend()
    plt.show()

# plot_colours()



class LaTeX:
    def __init__(self):

        self.Hp = r"$\mathcal{H}$"
        self.gt = r"$\tilde{g}$"
        self.x = r"$x$"
        self.eta = r"$\eta$"
        self.tau = r"$\tau$"

        self.X_e = r"X_e"
        self.Y_p = r"Y_p"

        self.diff = r"$\mathrm{d}$"

        self.ped = lambda s: r"_\mathrm{%s}"%s
        self.ap = lambda s: r"^\mathrm{%s}"%s
        self.unit = lambda u: r"~\mathrm{%s}"%u
        self.inv = lambda A: r"%s^{-1}"%A

        self.frac = lambda a, b: r"$\frac{%s}{%s}$"%(a, b)

    def __call__(self, str):
        str = str.replace("$", "")
        return r"$%s$"%str

    def dv(self, f=None, x=None, n=1):
        x = x or self.x
        f = f or self.Hp
        n = int(n)

        d = self.diff
        dx = self.diff + x
        if n>1:
            d += "^{%i}" %n 
            dx += "^{%i}" %n 
        dfdx = d + f + "/" + dx
        return self(dfdx)
    
    def normal(self, mu=0, sigma=1, splitline=False):
        s1 = "\mathcal{N}(\mu\!=\!%.2f," %mu
        s2 = "\,\sigma\!=\!%.2f)" %sigma
        if splitline:
            return self(s1) + "\n" + self(s2)
            # s += "\n"
        else:
        # s += "\,\sigma\!=\!%.2f)" %sigma
            return self(s1+s2)




pinpoint_kw = dict(color="slategrey", lw=.9, ls=":", alpha=.5)
overplot_kw = dict(color="darkslategrey", ls="-.", alpha=.5, lw=2.2)


mark_axis_kw = dict(reset=True, direction="out", length=17, width=.9, color="slategrey", labelcolor="slategrey")
mark_xaxis_kw = dict(**mark_axis_kw, top=False)
mark_yaxis_kw = dict(**mark_axis_kw, left=False, right=True, labelleft=False, labelright=True)



""" not ready to delete this """





# background = "#313332"
# face = background
# title_font = "MathJax_SansSerif"


# # reset_matplotlib()
#STIX MathJax SansSerif

# title_font = fm.FontProperties(["sans-serif"])

# title_props = {"fontweight":"bold", "fontsize":22, "fontfamily":"cursive"}

# sns.set_style("darkgrid", {"axes.facecolor": face})
# matplotlib.rcParams["figure.figsize"] = (10,6)

# plt.rcParams["mathtext.fontset"] = "stix"
# plt.rcParams["font.family"] = "sans-serif"#"STIXGeneral"
# plt.rcParams["font.serif"] = ["Times New Roman"]#"STIXGeneral"
# plt.rcParams["text.usetex"] = True
# plt.rcParams["savefig.bbox"] = "tight"
# plt.rcParams["figure.autolayout"] = True


# plt.rc("axes", **{"titlesize":24, "labelsize":16, "titlelocation":"left", "titleweight":"bold"})
# plt.rc("figure", autolayout=True)
# plt.rc("lines", linewidth=2)

# sns.set(rc={'figure.figsize':(12, 6), "lines.linewidth":2.2}, font_scale=1.1, font="sans-serif")
# #FIXME
# colors = [
#     sns.color_palette('husl')[-3],
#     sns.color_palette('husl')[-2],
#     sns.color_palette('husl')[-1],
#     'mediumorchid',
#     sns.color_palette('deep')[-1],
#     sns.color_palette('dark')[-1]
# ]

# plt.rc("axes", titlesize=18, labelsize=16, prop_cycle=cycler('color', colors), titlelocation="left", titleweight="bold")
# plt.rc("axes", titlesize=18, labelsize=16, titlelocation="left", titleweight="bold")
# plt.rc("legend", fontsize=14, shadow=True, fancybox=True)
# plt.rc("figure", figsize=(10,6))

















class ConstantsAndUnits:

    def __init__(self) -> None:

        #Basic units (here we use SI)
        self.m           = 1.0;                         # Length (in meters)
        self.s           = 1.0;                         # Time (in seconds)
        self.kg          = 1.0;                         # Kilo (in kilos)
        self.K           = 1.0;                         # Temperature (in Kelvins)

        m, s, kg, K = self.m, self.s, self.kg, self.K

        #Derived units
        self.km          = 1e3 * m;                     # Kilometers
        self.N           = kg*m/(s*s);                  # Newton
        self.J           = self.N*m;                    # Joule
        self.W           = self.J/s;                    # Watt
        self.Mpc         = 3.08567758e22 * m;           # Megaparsec
        self.eV          = 1.60217653e-19 * self.J;     # Electronvolt
        
        #Physical constants    
        self.k_b         = 1.38064852e-23 * self.J/K;        # Bolzmanns constant
        self.m_e         = 9.10938356e-31 * kg;              # Mass of electron
        self.m_H         = 1.6735575e-27 * kg;               # Mass of hydrogen atom
        self.c           = 2.99792458e8 * m/s;               # Speed of light
        self.G           = 6.67430e-11 * self.N*m*m/(kg*kg); # Gravitational constant
        self.hbar        = 1.054571817e-34 * self.J*s;       # Reduced Plancks constant
        self.sigma_T     = 6.6524587158e-29 * m*m;           # Thomas scattering cross-section
        self.Lambda_2s1s = 8.227 / s;                        # Transition time between 2s and 1s in Hydrogen
        self.H0_over_h   = 100 * self.km/s/self.Mpc;         # H0 / h
        self.epsilon_0   = 13.605693122994 * self.eV;        # Ionization energy for the ground state of hydrogen
        self.xhi0        = 24.587387 * self.eV;              # Ionization energy for neutral Helium
        self.xhi1        = 4.0 * self.epsilon_0;             # Ionization energy for singly ionized Helium
                
        


