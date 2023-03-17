from utils import *
SET_SUBDIR("milestone2")
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()
tex = LaTeX()

rec = read_ASCII("recombination")

x, X_e = rec[:,0], rec[:,1]

# fig, ax = plt.subplots()
# ax.plot(x, X_e)
# ax.set_yscale("log")

# save("electron_fraction")



tau = np.where((x>0)|(x<-12), np.nan, rec[:,3])
# print(x[np.nanargmin(np.abs(tau-1))])
dtaudx = np.where(np.isnan(tau), np.nan, rec[:,4])
ddtaudxx = np.where(np.isnan(tau), np.nan, rec[:,5])


# fig, ax = plt.subplots()
# ax.plot(x, tau, label=tex.tau)
# ax.plot(x, -dtaudx, label=tex("-"+tex.dv(tex.tau)))
# ax.plot(x, ddtaudxx, label=tex.dv(tex.tau, n=2))
# ax.set_yscale("log")
# ax.legend()
# save("optical_depth_misc")



gt = np.where(np.isnan(tau), np.nan, rec[:,6])
dgtdx = np.where(np.isnan(tau), np.nan, rec[:,7])
ddgtdxx = np.where(np.isnan(tau), np.nan, rec[:,8])


fig, ax = plt.subplots()
ax.plot(x, gt, label=tex.gt)
ax.legend()

fig, ax = plt.subplots()
ax.plot(x, dgtdx, label=tex.dv(tex.gt))
ax.plot(x, ddgtdxx, label=tex.dv(tex.gt, n=2))
ax.legend()


plt.show()

