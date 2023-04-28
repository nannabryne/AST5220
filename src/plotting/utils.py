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

        self.__Greeks()

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
        
    def TITLE(self, title):
        words = title.split()
        s = ""
        for word in words:
            s += "\mathsf{%s}~"%word
        return self(s)
    

    def __Greeks(self):
        self.Psi = self("\mathit{\Psi}")
        self.Phi = self("\mathit{\Phi}")
        self.Theta = self("\mathit{\Theta}")





overplot_kw = dict(color="darkslategrey", ls="-.", alpha=.5, lw=2.2)

f_of_x_kw =     dict(color="royalblue",   alpha=.9)
dfdx_of_x_kw =  dict(color="forestgreen", alpha=.7, lw=1.4)
ddfdxx_of_x_kw =dict(color="dodgerblue",  alpha=.7, lw=1.4)



def pinpoint_x(ax, x_list, x_label_list=None, colour="slategrey", style=":", width=.9, alpha=.5, top=False, tick=True, line=True):
    if x_label_list is None:
        tick = False

    if line:
        line_kw = dict(color=colour, lw=width, ls=style, alpha=alpha)
        for x in x_list:
            ax.axvline(x, **line_kw)

    if tick:
        if top:
            tick_kw = dict(reset=True, direction="out", length=9, width=width, color=colour, labelcolor=colour, top=True, bottom=False, labeltop=True, labelbottom=False)
        else:
            tick_kw = dict(reset=True, direction="out", length=17, width=width, color=colour, labelcolor=colour, top=False, bottom=True, labeltop=False, labelbottom=True)
        ax.xaxis.set_tick_params("minor", **tick_kw)
        # ax.set_xticks(x_list, labels=[tex(x_label) for x_label in x_label_list], minor=True)
        ax.set_xticks(x_list, labels=x_label_list, minor=True)



def pinpoint_y(ax, y_list, y_label_list, colour="slategrey", style=":", width=.9, alpha=.5, right=False):
    line_kw = dict(color=colour, lw=width, ls=style, alpha=alpha)
    if right:
        tick_kw = dict(reset=True, direction="out", length=17, width=width, color=colour, labelcolor=colour, right=True, left=False, labelright=True, labelleft=False)
    else:
        tick_kw = dict(reset=True, direction="out", length=17, width=width, color=colour, labelcolor=colour, right=False, left=True, labelright=False, labelleft=True)

    for y in y_list:
        ax.axhline(y, **line_kw)

    ax.yaxis.set_tick_params("minor", **tick_kw)
    ax.set_yticks(y_list, labels=y_label_list, minor=True)





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
                
        



class ColourCycles:
    def __init__(self, name="default", ax=None):
        
        if isinstance(name, str):
            try:
                eval(f"self.{name}()")
            except ValueError:
                print("provide valid arg")
                exit()

        elif isinstance(name, (list, tuple)):
            self.colours = name
        
        self.n_colours = len(self.colours)

        if ax is not None:
            self.set_cycle(ax)


    def default(self):
        self.colours = ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8']
        
    def derivatives(self):
        self.colours = ["royalblue", "forestgreen", "dodgerblue"]

    def set_cycle(self, ax):
        ax.set_prop_cycle(cycler("color", self.colours))

    def __getitem__(self, it):
        return self.colours[it]
    
    def __call__(self):
        return self.colours
    







def reset_matplotlib():
    plt.rcParams.update(plt.rcParamsDefault)


darkMode = False 


# COLOURS = ["#5d00ff", "#ef495e", "#ffa500", "#9989f2", "#c0c066", "#133bfc", "#cd7e26"]

# subCOLOURS = ["", "", "#FBAB18"]

# COLOURS = ['#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4']

COLOURS = ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8']
subCOLOURS = ["#DA8F84", "#83B5D2", "#CDC8EC", "w", "#F3C77C", "#B9D092", "#F7DBDC"]


__colours = ColourCycles()
# print(plt.style.library["fivethirtyeight"])


def _set_params():
    # face = plt.style.library["ggplot"]["axes.facecolor"]
    face = "#D4CCCD"
    custom_params = {
        "figure.figsize":(10,5), "figure.autolayout":True,
        "savefig.dpi":300, "savefig.bbox":"tight",
        "lines.linewidth":2.8, 
        "text.usetex":True, 
        "mathtext.fontset":"custom", "mathtext.fallback":"stix", "font.family":"STIXGeneral",
        "axes.labelsize":22, "axes.titlesize":22, "axes.titlelocation":"left", "axes.titleweight":"bold", 
        "axes.prop_cycle":cycler('color', __colours()),
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
    colours = __colours()
    for i in range(len(colours)):
        c = colours[i].strip("#")
        ax.plot(x, np.sin(x+i*np.pi*0.3)*x-i*1.5, label=f"{i}: {c}")
    ax.legend()
    plt.show()

# plot_colours()





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


