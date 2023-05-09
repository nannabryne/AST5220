from utils import *
import plotsrc_milestone_four as PLOT

const = ConstantsAndUnits()

tex = LaTeX()



power_data = read_ASCII("cells.txt")

ell, Dell = power_data[:,0], power_data[:,1]


fig, ax = plt.subplots()

ax.plot(ell, Dell)
# ax.set_xlim(2, 1200)
ax.set_xscale("log")

plt.show()