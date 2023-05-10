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





data = read_ASCII("power_something_ell6.txt")
k = data[:,0]

P_matter = data[:,1]

fig, ax = plt.subplots()
ax.plot(k, P_matter)
ax.set_yscale("log")
ax.set_xscale("log")

plt.show()