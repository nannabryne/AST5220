from utils import *
SET_SUBDIR("milestone2")
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()
tex = LaTeX()

rec = read_ASCII("recombination")

x_star = -6.985

x, X_e = rec[:,0], rec[:,1]

fig, ax = plt.subplots()
ax.plot(x, X_e, label=tex("X_e"))
ax.set_yscale("log")
ax.legend()
ax.set_xlabel(tex.x)

ax.axvline(x_star, **pinpoint_kw)
ax.xaxis.set_tick_params("minor", **mark_xaxis_kw)
ax.set_xticks([x_star], labels=[tex("x_*")], minor=True)

ax.axhline(0.99, **pinpoint_kw)
ax.yaxis.set_tick_params("minor", **mark_yaxis_kw)
ax.set_yticks([0.99], labels=[tex("0.99")], minor=True)

save("electron_fraction")




tau = np.where((x>0)|(x<-12), np.nan, rec[:,3])
# # print(x[np.nanargmin(np.abs(tau-1))])
dtaudx = np.where(np.isnan(tau), np.nan, rec[:,4])
dtaudx[0] = np.nan
dtaudx[-1] = np.nan
ddtaudxx = np.where(np.isnan(tau), np.nan, rec[:,5])
ddtaudxx[0] = np.nan
ddtaudxx[1] = np.nan
ddtaudxx[-1] = np.nan
ddtaudxx[-2] = np.nan

# tau = rec[:,3]
# dtaudx = rec[:,4]
# ddtaudxx = rec[:,5]

fig, ax = plt.subplots()
ax.plot(x, tau, c="royalblue", alpha=.7, label=tex.tau)
ax.plot(x, -dtaudx, c="forestgreen", alpha=.7, label=tex("-"+tex.dv(tex.tau)))
ax.plot(x, ddtaudxx, c="dodgerblue", alpha=.7, label=tex.dv(tex.tau, n=2))
ax.set_xlabel(tex.x)
ax.set_yscale("log")
ax.legend()

ax.axvline(x_star, **pinpoint_kw)
ax.xaxis.set_tick_params("minor", **mark_xaxis_kw)
ax.set_xticks([x_star], labels=[tex("x_*")], minor=True)

# ax.axhline(1.00, **pinpoint_kw)
# ax.yaxis.set_tick_params("minor", **mark_yaxis_kw)
# ax.set_yticks([0.999], labels=[tex("1.0")], minor=True)


save("optical_depth_misc")



gt = np.where(np.isnan(tau), np.nan, rec[:,6])
dgtdx = np.where(np.isnan(tau), np.nan, rec[:,7])
ddgtdxx = np.where(np.isnan(tau), np.nan, rec[:,8])


fig, ax = plt.subplots()
# ax.plot(x, gt/np.nanmax(gt), label=tex.gt)
# ax.plot(x, dgtdx/np.nanmax(dgtdx), label=tex.dv(tex.gt))
# ax.plot(x, ddgtdxx/np.nanmax(ddgtdxx), label=tex.dv(tex.gt, n=2))

ax.plot(x, gt, c="royalblue", alpha=.7, label=tex.gt)
ax.plot(x, dgtdx/10, c="forestgreen", alpha=.7, label=tex(tex.frac(1,10) + tex.dv(tex.gt)))
ax.plot(x, ddgtdxx/300, c="dodgerblue", alpha=.7, label=tex(tex.frac(1,300)+tex.dv(tex.gt, n=2)))

ax.axvline(x_star, **pinpoint_kw)
ax.xaxis.set_tick_params("minor", **mark_xaxis_kw)
ax.set_xticks([x_star], labels=[tex("x_*")], minor=True)

ax.set_xlabel(tex.x)
ax.set_xlim(-7.7, -5.7)
ax.legend()
save("visibility_function_misc")

plt.show()



