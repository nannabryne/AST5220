from utils import *
SET_SUBDIR("milestone1")


class ColourCyclesI(ColourCycles):
    def densityparams(self):
        self.colours = ["#ffa500","#5d00ff", "#ef495e"]

    def densityparams_shy(self):
        self.colours = ["#fffb00", "#ae7fff", "#e7aeb0"]


tex = LaTeX()



def __colour_eras(ax, info, xlims=None, vlines=False, text=False):

    xA, xD = xlims or (-20, 5)
    xB, xC = info["x_eq"], info["x_Lambda"]
    y0, y1 = ax.get_ylim()

    c = ColourCyclesI("densityparams_shy")

    if text: 
        box_kw = dict(facecolor="slategrey", alpha=.1)
        ax.annotate("radiation era", (xA, y1), color=c[0], bbox=box_kw)
        ax.annotate("matter era", (xB, y1), color=c[1], bbox=box_kw)
        ax.annotate("DE era", (xC, y1), color=c[2], bbox=box_kw)

    ax.fill_betweenx([y0,y1], xA, xB, color=c[0], alpha=.1)
    ax.fill_betweenx([y0,y1], xB, xC, color=c[1], alpha=.1)
    ax.fill_betweenx([y0,y1], xC, xD, color=c[2], alpha=.1)

    if vlines:
        x_acc, x_0 = info["x_acc"], info["x_0"]
        for xi in [xB, x_acc, xC, x_0]:
            ax.axvline(xi, lw=.9, ls=":", color="w", alpha=.4)


def __compare_with_analytical(ax, analytical, which, ls="-."):

    rad = analytical["radiation"]
    mat = analytical["matter"]
    vac = analytical["vacuum"]


    kwargs = overplot_kw.copy()
    kwargs["ls"] = ls
    ax.plot(rad["x"], rad[which], **kwargs)
    ax.plot(mat["x"], mat[which], **kwargs)
    ax.plot(vac["x"], vac[which], **kwargs)





def DensityParameters(df, info, savefig=True):
    fig, ax = plt.subplots(figsize=(10,6))
    x = df["x"]
    c = ColourCyclesI("densityparams", ax)
    ax.plot(x, df["OmegaR"],  label=tex("\Omega" + tex.ped("r")))
    ax.plot(x, df["OmegaM"],  label=tex("\Omega" + tex.ped("m")))
    ax.plot(x, df["OmegaL"],  label=tex("\Omega_\Lambda"))


    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper right', ncols=3, borderaxespad=0.)
    ax.set_xlabel(tex.x)

    __colour_eras(ax, info, text=True)

    if savefig:
        save("density_params")


def HubbleDerivatives(df, info, analytical, savefig=True):
    fig, ax = plt.subplots()
    c = ColourCycles("derivatives", ax)

    x, Hp = df["x"], df["Hp"]
    ax.plot(x, df["dHpdx"]/Hp,   c=c[1], label=tex(tex.inv(tex.Hp)+ tex.dv()))
    ax.plot(x, df["ddHpdxx"]/Hp, c=c[2],  label=tex(tex.inv(tex.Hp)+ tex.dv(n=2))) 

    __compare_with_analytical(ax, analytical, "Hubble_der", ":")
    __compare_with_analytical(ax, analytical, "Hubble_doubleder", "--")
    
    ax.plot(0,-2)
    __colour_eras(ax, info,  vlines=False)
    
    ax.set_xlim(-16, 4)
    ax.set_ylim(-1.2, 1.4)
    
    ax.legend()
    ax.set_xlabel(tex.x)

        
    if savefig:
        save("Hubble_derivatives")


def ConformalTime_HubbleParameter(df, info, analytical, savefig=True):
    fig, ax = plt.subplots()


    ax.plot(df["x"], df["eta_c_org"]*df["Hp_org"], c="orangered", label=tex(tex.Hp+"\eta/c"))
    __colour_eras(ax, info,  vlines=False)
    __compare_with_analytical(ax, analytical, "eta_Hubble", "-.") # dunno about the last one

    pinpoint_x(ax, [info["x_acc"]], [tex("x" + tex.ped("acc"))])
    
    ax.set_xlim(-16, 4)
    ax.legend()    
    ax.set_xlabel(tex.x)
    ax.set_ylim(0.5, 4.5)

    if savefig:
        save("eta_Hubble")


def HubbleParameter(df, info, analytical, savefig=True):

    fig, ax = plt.subplots()
    c = ColourCycles("derivatives", ax)

    ax.plot(df["x"], df["Hp"], c=c[0], label=tex.Hp)

    __colour_eras(ax, info, vlines=False)
    __compare_with_analytical(ax, analytical, "Hubble", "-.")

    ax.set_ylabel(r"100 km s$^{-1}$ Mpc$^{-1}$")
    ax.set_xlabel(tex.x)


    ax.set_ylim(0.1, 1e5)
    ax.set_xlim(-16, 4)
    ax.legend()

    ax.set_yscale("log")


    if savefig:
        save("hubble_param")


def TimeMeasures(df, info, todays, savefig=True):

    fig, ax = plt.subplots()

    x = df["x"]
    ax.plot(x, df["t"],     c="dodgerblue",label=r"$t$")
    ax.plot(x, df["eta_c"], c="orangered",label=r"$\eta/c$")

    ax.set_ylabel("Ga")
    ax.set_xlabel(tex.x)
    # self._mark_acc(ax)
    ax.set_yscale("log")

    ax.set_xlim(-12, 2)
    ax.set_ylim(3e-8, 4e2)

    ax.legend()
    pinpoint_y(ax, [todays["t"], todays["eta_c"]], [tex("t_0"), tex("\eta_0/c")], right=True)

    if savefig:
        save("time_measures") 


def LuminosityDistance(df_obs, df_pri, df_post, savefig=True):
    fig, ax = plt.subplots(figsize=(10,4))

    ax.set_xscale("log")

    z = df_obs["z"]
    ax.errorbar(z, df_obs["dL"]/z, df_obs["err"]/z, label=tex("d_L"+ tex.ap("obs") + "/z"), **obs_err_kw)
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # avoid division by zero
    z = df_pri["z"]
    idx0 = np.argmin(np.abs(z))
    z[idx0] = 1e-3
    df_pri["dL"][idx0] = 1
    df_post["dL"][idx0] = 1

    ax.plot(z, df_pri["dL"]/z, alpha=.7, c=COLOURS[4], label=r"$d_L/z$ (pre MCMC)")
    ax.plot(z, df_post["dL"]/z, alpha=.7, c=COLOURS[5], label=r"$d_L/z$ (post MCMC) ")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.errorbar(z, df_obs["dL"]/z, df_obs["err"]/z, **obs_err_kw) # overplot ...

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"Gpc")

    ax.legend()

    if savefig:
        save("lum_distance_vs_redshift")


def Omega_PhaseSpace(df, todays_pri, savefig=True):
    
    # tmp
    org_x = 0.05+0.267
    org_y = 1-(org_x+5.50896e-05)

    org_x = todays_pri["OmegaM"]
    org_y = todays_pri["OmegaL"]

    chi2 = df["chi2"]
    Omega_M0 = df["OmegaM0"]
    Omega_L0 = df["OmegaL0"]

    chi2_min = np.min(chi2)
    sigma1, sigma2 = 3.53, 8.02 # sigma levels
    thresh1 = (chi2 < chi2_min + sigma1)
    thresh2 = (chi2 < chi2_min + sigma2)

    xx, yy, zz = r"$\Omega_{\mathrm{m}0}$", r"$\Omega_{\Lambda 0}$", r"$\chi2$"
    df = pd.DataFrame({xx:Omega_M0, yy:Omega_L0, zz:chi2})

    from scipy.stats import norm


    g = sns.JointGrid(data=df, x=xx, y=yy)
    fig, ax = g.figure, g.ax_joint
    ax_M, ax_L = g.ax_marg_x, g.ax_marg_y
    g.fig.set_figwidth(10)
    g.fig.set_figheight(7)

    
    # Plot histograms:
    g.plot_marginals(sns.histplot, kde=False, color=COLOURS[0], stat="density", alpha=.66)
    

    # Scatter plot confidence interval(s):
    c2, c1 = COLOURS[1], COLOURS[2]
    kwargs = dict(s=17, alpha=.7)

    second = ax.scatter(Omega_M0[thresh2], Omega_L0[thresh2], color=c2, edgecolor=c2, label=r"$2\sigma$ level", **kwargs)
    first = ax.scatter(Omega_M0[thresh1], Omega_L0[thresh1], color=c1, edgecolor=c1, label=r"$1\sigma$ level", **kwargs)
    
    # mark flat space lim
    x = np.linspace(-0.1, 1.1, 100)
    y = 1 - x
    ax.plot(x, y, ls="--", c="slategrey")
    # mark acceleration lim
    y = .5*x
    ax.plot(x, y, ls=":", c="slategrey")


    # ax.an notate(r"$\Omega_{\mathrm{K}0} = 0$", (x[60], y[60]), rotation=-30)

    # limits
    M_lims = (0, 1.)# (-0.02, np.max(Omega_M0[thresh2])+0.08)
    L_lims = (0, 1.5) #(-0.02, np.max(Omega_Lambda0[thresh2])+0.12)
    ax.set_xlim(*M_lims)
    ax.set_ylim(*L_lims)
    
    
    # Gaussians
    leg_kw = dict(fontsize=15, handlelength=.6, handletextpad=.1)#, borderpad=.2)

    MM = np.linspace(*M_lims, 100)
    mu, sigma = norm.fit(Omega_M0)
    OmM_pdf = norm.pdf(MM, mu, sigma)
    pdf, = ax_M.plot(MM, OmM_pdf)
    ax_M.legend([pdf], [tex.normal(mu, sigma)], **leg_kw)
    # ax_M.legend([pdf], [r"$\mathcal{N}(%.2f, %.2f)$" %(sigma, mu)], **leg_kw)
    

    LL = np.linspace(*L_lims, 100)
    mu, sigma = norm.fit(Omega_L0)
    OmL_pdf = norm.pdf(LL, mu, sigma)
    pdf, = ax_L.plot(OmL_pdf, LL)
    ax_L.legend([pdf], [tex.normal(mu, sigma, splitline=True)], borderaxespad=.05, **leg_kw)
    
    # add our original point
    c = COLOURS[4]
    kw = dict(lw=5, ms=16)
    org_p, = ax.plot(org_x, org_y, "+", c=c, label="fiducial", **kw)
    # add new point
    new_p, = ax.plot(Omega_M0[np.nanargmin(chi2)], Omega_L0[np.nanargmin(chi2)], "+", c=COLOURS[5], label="new best-fit", **kw)

    ax.legend(handles=[first, second, org_p, new_p], markerscale=1.7)
    # ax.legend(reverse=True, markerscale=2.)

    fig.tight_layout = True   

    if savefig:
        save("omega_phase_space")


def Curvature_PDF(df, todays_pri, savefig=True):
    from scipy.stats import norm

    Omega_K0 = df["OmegaK0"]
    fig, ax = plt.subplots(figsize=(10,4))
    sns.histplot(Omega_K0, ax=ax, stat="density", alpha=.66)
    
    K_lims = (-1., 1.)
    # x0, x1 = ax.get_xlim()
    ax.set_xlabel(r"$\Omega_{\mathrm{K}0}$")

    mu, sigma = norm.fit(Omega_K0)
    xx = np.linspace(*K_lims, 100)
    ax.plot(xx, sp.stats.norm.pdf(xx, mu, sigma), label=tex.normal(mu, sigma))
    ax.axvline(todays_pri["OmegaK"], color=COLOURS[4], lw=1.4, label="fiducial")
   
    ax.legend()

    if savefig:
        save("curvature_pdf")


def Hubble_PDF(df, todays_pri, savefig=True):

    H0 = df["H0"]
    fig, ax = plt.subplots(figsize=(10,4))
    # ax.hist(H0[thresh], 60)
    from scipy.stats import norm
    
    sns.histplot(H0, ax=ax, stat="density", alpha=.76)
    # x0, x1 = ax.get_xlim()
    h_lims = (0.65, .75)
    ax.set_xlabel(r"$h$")

    mu, sigma = norm.fit(H0)
    xx = np.linspace(*h_lims, 100)
 
    ax.plot(xx, sp.stats.norm.pdf(xx, mu, sigma), label=tex.normal(mu, sigma))
    ax.axvline(todays_pri["H"], color=COLOURS[4], lw=1.4, label="fiducial")
    ax.legend()
    


    if savefig:
        save("hubble_pdf")







