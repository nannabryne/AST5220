from utils import *
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()

fig, ax = plt.subplots()

c = ColourCycles(["red", "blue", "k", "yellow"], ax)

plt.plot([0,1,2], [0,3, 8])

plt.plot([0,1,2], [1,4, 9])


plt.plot([0,1,2], [-1,0, 1])

plt.show()