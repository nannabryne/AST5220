from utils import *
SET_SUBDIR("milestone2")
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()
tex = LaTeX()

rec = read_ASCII("recombination")

rec_Saha = read_ASCII("recombination_Saha_only")

x_rec = -6.985

x, Xe = rec[:,0], rec[:,1]
x_Saha, Xe_Saha = rec_Saha[:,0], rec_Saha[:,1]

fig, ax = plt.subplots()
# ax.plot(x, Xe, color=COLOURS[0], label=tex("X_e")+" ("+ tex("X_e"+tex.ap("Saha")+"=0.99") + ")")
# ax.plot(x_Saha, Xe_Saha, ls="--", color=COLOURS[0], alpha=.7, label=tex("X_e")+" ("+ tex("X_e"+tex.ap("Saha")+"=0.00") + ")")
ax.plot(x, Xe, color=COLOURS[0], label=tex("X_e")+ " (Saha+Peebles)")
ax.plot(x_Saha, Xe_Saha, label=tex("X_e")+ " (Saha)", **overplot_kw)

ax.set_yscale("log")
ax.legend()
ax.set_xlabel(tex.x)

ylim = (1e-4,1.5)
ax.set_ylim(ylim)
ax.set_xlim(-8.4,-3.6)

ax.axvline(x_rec, **pinpoint_kw)
ax.xaxis.set_tick_params("minor", **mark_xaxis_kw)
ax.set_xticks([x_rec], labels=[tex("x_*")], minor=True)

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
ax.plot(x, tau, label=tex.tau, **f_of_x_kw)
ax.plot(x, -dtaudx, label=tex("-"+tex.dv(tex.tau)), **dfdx_of_x_kw)
ax.plot(x, ddtaudxx, label=tex.dv(tex.tau, n=2), **ddfdxx_of_x_kw)
ax.set_xlabel(tex.x)
ax.set_yscale("log")
ax.legend()

ax.axvline(x_rec, **pinpoint_kw)
ax.xaxis.set_tick_params("minor", **mark_xaxis_kw)
ax.set_xticks([x_rec], labels=[tex("x_*")], minor=True)

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

ax.plot(x, gt,          label=tex.gt, **f_of_x_kw)
ax.plot(x, dgtdx/10,    label=tex(tex.frac(1,10) + tex.dv(tex.gt)), **dfdx_of_x_kw)
ax.plot(x, ddgtdxx/300, label=tex(tex.frac(1,300)+tex.dv(tex.gt, n=2)), **ddfdxx_of_x_kw)

ax.axvline(x_rec, **pinpoint_kw)
ax.xaxis.set_tick_params("minor", **mark_xaxis_kw)
ax.set_xticks([x_rec], labels=[tex("x_*")], minor=True)

ax.set_xlabel(tex.x)
ax.set_xlim(-7.7, -5.7)
ax.legend()
save("visibility_function_misc")

plt.show()



