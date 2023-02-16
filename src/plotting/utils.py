import os
import sys
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pathlib as pl
import matplotlib
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import numpy as np




CURRENT_PATH    = os.path.abspath(".") + "/"
INPUT_PATH      = CURRENT_PATH + "../../input/"
OUTPUT_PATH     = CURRENT_PATH + "../../output/"
DATA_PATH       = OUTPUT_PATH + "data/"
FIGS_PATH       = OUTPUT_PATH + "figures/"


CURRENT_PATH    = os.path.abspath(".") + "/"
DATA_PATH       = CURRENT_PATH + "../../data/"
INPUT_PATH      = DATA_PATH + "input/"
OUTPUT_PATH     = DATA_PATH + "output/"
FIGS_PATH       = CURRENT_PATH + "../../figs/"

# CURRENT_PATH    = pl.Path(__file__).parent
# OUTPUT_PATH     = CURRENT_PATH / pl.Path("/output")
# DATA_PATH       = OUTPUT_PATH / pl.Path("/data")
# FIGS_PATH       = OUTPUT_PATH / pl.Path("/figures")

def read_ASCII(filename, skiprows=0):
    if not filename.endswith(".txt"):
        filename += ".txt"
    try: 
        data = np.loadtxt(OUTPUT_PATH + filename)
    except FileNotFoundError:
        data = np.loadtxt(INPUT_PATH + filename)
    return np.asarray(data)
    


def read_ASCII2(filename, skiprows=0):
    if not filename.endswith(".txt"):
        filename += ".txt"
    try: 
        data = np.loadtxt(OUTPUT_PATH + filename)
    except FileNotFoundError:
        data = np.loadtxt(INPUT_PATH + filename)
    return np.asarray(data)




def save(filename, pdf=True):
    file = FIGS_PATH + filename.strip(".pdf")+".pdf"
    plt.savefig(file)




sns.set_style("darkgrid")
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "STIXGeneral"
plt.rcParams["text.usetex"] = True

plt.rc("figure", autolayout=True)
plt.rc("lines", linewidth=2)


#FIXME
colors = [
    sns.color_palette('husl')[-3],
    sns.color_palette('husl')[-2],
    sns.color_palette('husl')[-1],
    'mediumorchid',
    sns.color_palette('deep')[-1],
    sns.color_palette('dark')[-1]
]

plt.rc("axes", titlesize=18, labelsize=16, prop_cycle=cycler('color', colors))
plt.rc("legend", fontsize=14, shadow=True)

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
        self.J           = self.N*m;                         # Joule
        self.W           = self.J/s;                         # Watt
        self.Mpc         = 3.08567758e22 * m;           # Megaparsec
        self.eV          = 1.60217653e-19 * self.J;          # Electronvolt
        
        #Physical constants    
        self.k_b         = 1.38064852e-23 * self.J/K;        # Bolzmanns constant
        self.m_e         = 9.10938356e-31 * kg;         # Mass of electron
        self.m_H         = 1.6735575e-27 * kg;          # Mass of hydrogen atom
        self.c           = 2.99792458e8 * m/s;          # Speed of light
        self.G           = 6.67430e-11 * self.N*m*m/(kg*kg); # Gravitational constant
        self.hbar        = 1.054571817e-34 * self.J*s;       # Reduced Plancks constant
        self.sigma_T     = 6.6524587158e-29 * m*m;      # Thomas scattering cross-section
        self.lambda_2s1s = 8.227 / s;                   # Transition time between 2s and 1s in Hydrogen
        self.H0_over_h   = 100 * self.km/s/self.Mpc;              # H0 / h
        self.epsilon_0   = 13.605693122994 * self.eV;        # Ionization energy for the ground state of hydrogen
        self.xhi0        = 24.587387 * self.eV;              # Ionization energy for neutral Helium
        self.xhi1        = 4.0 * self.epsilon_0;             # Ionization energy for singly ionized Helium
                
        


