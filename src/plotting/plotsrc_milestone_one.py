from utils import *

""" CRAZY code, will fix this ASAP """




# colours related to some vital quantity
mainColours = {
    "R": "#ffa500", # radiation
    "M": "#5d00ff", # matter
    "L": "#ef495e",  # cosmological constant
    "Hubble": ["g", "b", "r"]
}

# colours derived from the above (supposed to be similar, but not the same)
subColours = {
    "R": "#fffb00",#'#fee090'
    "M": "#ae7fff",
    "L": "#e7aeb0"
}



tex = LaTeX()

class MilestoneI:

    def __init__(self):
        self.SaveMode()
        self.ShowMode()

        self.subdir = "milestone1/"

    def DefineBasics(self, x_array, indices):
        #   Define basics
        self.x = x_array
        self.idx = indices

    def _colour_eras(self, ax,  xlims=None, vlines=False, text=False):
        xA, xD = xlims or (np.nanmin(self.x), np.nanmax(self.x))
        xB, xC = self.x[self.idx["RM"]], self.x[self.idx["MDE"]]
        y0, y1 = ax.get_ylim()

        if text: 
            box_kw = dict(facecolor="slategrey", alpha=.1)
            ax.annotate("radiation era", (xA, y1), color=subColours["R"], bbox=box_kw)
            ax.annotate("matter era", (xB, y1), color=subColours["M"], bbox=box_kw)
            ax.annotate("DE era", (xC, y1), color=subColours["L"], bbox=box_kw)

        ax.fill_betweenx([y0,y1], xA, xB, color=subColours["R"], alpha=.1)
        ax.fill_betweenx([y0,y1], xB, xC, color=subColours["M"], alpha=.1)
        ax.fill_betweenx([y0,y1], xC, xD, color=subColours["L"], alpha=.1)

        if vlines:
            x_acc, x_0 = self.x[self.idx["acc"]], self.x[self.idx["0"]]
            for xi in [xB, x_acc, xC, x_0]:
                ax.axvline(xi, lw=.9, ls=":", color="w", alpha=.4)

    def _mark_milestones(self, ax):
        
        xB, xC, x_a, x_0 = self.x[self.idx["RM"]], self.x[self.idx["MDE"]], self.x[self.idx["acc"]], self.x[self.idx["0"]]
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())

        ax2.spines[["left", "right", "top"]].set_visible(False) 
        # ax2.axis("off")
        # c = "#dd00ff"
        # c = "#133bfc"
        c = "slategrey"
        ax2.xaxis.set_tick_params("major", reset=True, direction="out", color=c, length=17, width=0.9, top=False, grid_linewidth=0.0, labelsize=18, labelcolor=c)
        # ax2.xaxis.set_tick_params("minor", reset=True, direction="out", color=c, length=36, width=0.9, top=False, grid_linewidth=0.0, labelsize=16, labelcolor=c)
        ax2.xaxis.set_tick_params("minor", reset=True, direction="in", color=c, length=8, width=0.9, top=False, grid_linewidth=0.0, labelsize=18, labelcolor=c, pad=-24)
        ax2.grid(False)

        for xi in [xB, x_a, xC, x_0]:
            ax.axvline(xi, lw=.9, ls=":", color=c, alpha=.5)

        ax2.set_xticks([xB,  xC], labels=[r"$x_\mathrm{eq}$",r"$x_\Lambda$"], minor=True)
        ax2.set_xticks([x_a, x_0], labels=[r"$x_\mathrm{a}$",r"$x_0$"], minor=False)

    def _mark_acc(self, ax):
        # c, lc = "#5CE1E6", "#5CE1E6" # tick colour, line colour
        # c, lc = "k", "w"
        c = "slategrey"
        ax.xaxis.set_tick_params("minor", reset=True, direction="out", length=17, width=.9, top=False, color=c, labelcolor=c)
        x_a = self.x[self.idx["acc"]]
        ax.set_xticks([x_a], labels=[r"$x_\mathrm{acc}$"], minor=True)
        ax.axvline(x_a, lw=.9, ls=":", color=c, alpha=.5)
        
    def _compare_with_analytical(self, ax, val_list, ls=":", label=False):
        xA, xB, xC, xD = np.nanmin(self.x), self.x[self.idx["RM"]], self.x[self.idx["MDE"]], np.nanmax(self.x)

        kwargs = dict(color="darkslategrey", ls=ls, alpha=.5, lw=2.2)
        yAB, yBC, yCD = val_list
        if isinstance(yAB, (tuple, list, np.ndarray)):
            rAB, rBC, rCD = (self.x>=xA) & (self.x<xB), (self.x>=xB) & (self.x<xC), (self.x>=xC) & (self.x<=xD)
            # rad dom:
            if label:
                ax.plot(self.x[rAB], yAB[rAB], label=" ", **kwargs)  
            else:
                ax.plot(self.x[rAB], yAB[rAB], **kwargs)
            # matter dom:
            ax.plot(self.x[rBC], yBC[rBC], **kwargs)
            # CC dom:
            ax.plot(self.x[rCD], yCD[rCD], **kwargs)

        
        else:
            # rad dom:
            if label:
                ax.plot((xA, xB), (yAB, yAB), label=" ", **kwargs)  
            else:
                ax.plot((xA, xB), (yAB, yAB), **kwargs)
            # matter dom:
            ax.plot((xB, xC), (yBC, yBC), **kwargs)
            # CC dom:
            ax.plot((xC, xD), (yCD, yCD), **kwargs)

        if label:
            ax.legend(ncols=2)

    def SaveMode(self, mode="on"):
        if mode == "on":
            self.save = True
        elif mode == "off":
            self.save = False

    def ShowMode(self, mode="off"):
        if mode == "on":
            self.show = True
        elif mode == "off":
            self.show = False


    def DensityParameters(self, Omega_M, Omega_Lambda, Omega_R):
    
        fig, ax = plt.subplots(figsize=(10,6))
        ax.plot(self.x, Omega_M,      c=mainColours["M"],  label=tex("\Omega" + tex.ped("m")))
        ax.plot(self.x, Omega_Lambda, c=mainColours["L"],  label=tex("\Omega_\Lambda"))
        ax.plot(self.x, Omega_R,      c=mainColours["R"],  label=tex("\Omega" + tex.ped("r")))


        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper right', ncols=3, borderaxespad=0.)
        ax.set_xlabel(tex.x)
        # ax.set_title("Density parameters", fontweight="bold", fontname="sans-serif")
        
        # self._mark_milestones(ax)
        self._mark_acc(ax)
        self._colour_eras(ax, text=True)
    

        if self.save:
            save(self.subdir + "density_params")
        if self.show:
            plt.show()

    
    def HubbleDerivatives(self, Hp, dHp_dx, ddHp_dxx):
        fig, ax = plt.subplots()

        # ax.plot(self.x, dHp_dx/Hp,   c="forestgreen", label=r"$\mathcal{H}^{-1}\mathrm{d}\mathcal{H}/\mathrm{d}x$")
        ax.plot(self.x, dHp_dx/Hp,   c="forestgreen", label=tex(tex.inv(tex.Hp)+ tex.dv()))
        ax.plot(self.x, ddHp_dxx/Hp, c="dodgerblue", label=tex(tex.inv(tex.Hp)+ tex.dv(n=2))) 
       
        self._compare_with_analytical(ax, (-1, -.5, 1), ":", label=True)    #dHpdx
        self._compare_with_analytical(ax, (1, .25, 1), "--", label=True)    #ddHpdxx
        ax.plot(0,-2)
        self._colour_eras(ax, vlines=False)
        self._mark_acc(ax)
        ax.set_xlim(-16, 4)
        ax.set_ylim(-1.2, 1.4)
        
        # ax.legend(ncols=2)
        ax.set_xlabel(tex.x)

        # self._mark_milestones(ax)

            
        if self.save:
            save(self.subdir + "Hubble_derivatives")
        if self.show:
            plt.show()

    
    def ConformalTime_HubbleParametter(self, eta_c, Hp):
        fig, ax = plt.subplots()
        ax.plot(self.x, eta_c*Hp, c="orangered",label=tex(tex.Hp+"\eta/c"))
        self._colour_eras(ax, vlines=False)
        self._compare_with_analytical(ax, (1, 2, 100), "-.", label=True) # dunno about the last one
        self._mark_acc(ax)
        ax.set_xlim(-16, 4)
        # ax.legend(ncols=2)    
        ax.set_xlabel(tex.x)
        ax.set_ylim(0.5, 4.5)

        if self.save:
            save(self.subdir + "eta_Hubble")
        if self.show:
            plt.show()


        

    
    def HubbleParameter(self, Hp, OmegaR0, OmegaM0, OmegaL0):

        fig, ax = plt.subplots()

        ax.plot(self.x, Hp, c="royalblue",label=tex.Hp)
        H0 = Hp[self.idx["0"]]

        
        self._colour_eras(ax, vlines=False)
        # self._mark_milestones(ax)

        H_AB = H0 * np.sqrt(OmegaR0) * np.exp(-self.x)
        H_BC = H0 * np.sqrt(OmegaM0) * np.exp(-0.5*self.x)
        H_CD = H0 * np.sqrt(OmegaL0) * np.exp(self.x)
        self._compare_with_analytical(ax, (H_AB, H_BC, H_CD), "-.", label=True)
        self._mark_acc(ax)
        ax.set_ylabel(r"100 km s$^{-1}$ Mpc$^{-1}$")
        ax.set_xlabel(tex.x)
        # ax.set_xlim(np.min(x), 0)

        ax.set_ylim(0.1, 1e5)
        ax.set_xlim(-16, 4)

        ax.set_yscale("log")


        # ax.legend(ncols=2)

        if self.save:
            save(self.subdir + "hubble_param")
        if self.show:
            plt.show()



    def TimeMeasures(self, t, eta_c):
        fig, ax = plt.subplots()

        ax.plot(self.x, t,     c="dodgerblue",label=r"$t$")
        ax.plot(self.x, eta_c, c="orangered",label=r"$\eta/c$")

        ax.set_ylabel("Ga")
        ax.set_xlabel(tex.x)
        # self._mark_acc(ax)
        ax.set_yscale("log")

        # c = "slategrey"
        # ax.yaxis.set_tick_params("minor", reset=True, direction="out", length=17, width=.9, top=False, color=c, labelcolor=c)
        # x_a = self.x[self.idx["acc"]]
        # ax.set_xticks([x_a], labels=[r"$x_\mathrm{acc}$"], minor=True)
        # ax.axvline(x_a, lw=.9, ls=":", color=c, alpha=.5)


        ax.set_xlim(-12, 2)
        ax.set_ylim(3e-8, 4e2)

        ax.legend()

        # ax2 = ax.twinx()
        # ax2.set_ylim(ax.get_ylim())
        # ax2.spines[["left", "bottom", "top"]].set_visible(False) 
        c = "slategrey"
        ax.yaxis.set_tick_params("minor", reset=True, direction="out", length=17, width=.9, left=False, right=True, color=c, labelcolor=c,  labelleft=False, labelright=True)
        t_0 = t[self.idx["0"]]
        eta_0 = eta_c[self.idx["0"]]
        ax.set_yticks([t_0, eta_0], labels=[r"$t_0$", r"$\eta_0/c$"], minor=True)
        ax.axhline(t_0, lw=.9, ls=":", color=c, alpha=.5)
        ax.axhline(eta_0, lw=.9, ls=":", color=c, alpha=.5)

        if self.save:
            save(self.subdir + "time_measures") 
        if self.show:
            plt.show()


    
    def LuminosityDistance(self, z_obs, dL_obs, err_obs, z_comp, dL_comp):
        fig, ax = plt.subplots(figsize=(10,4))

        ax.set_xscale("log")
        # ax.set_yscale("log")
    
        ax.errorbar(z_obs, dL_obs/z_obs, err_obs/z_obs, elinewidth=1.1, capsize=2, linestyle="", marker="o", ms=4, label=tex("d_L"+ tex.ap("obs") + "/z"))
        # ax.plot(z_obs, dL_obs, 'o')#,         c="orangered")
        
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # ax.plot(10000, 0)
        # print("\n", dL_comp)
        # print(z_comp)
        # print(dL_comp/z_comp)
        # reg = (z_comp>xlim[0]-0.1) & (z_comp<xlim[1]+1)
        # dL_comp = dL_comp[reg]
        # z_comp = z_comp[reg]
        # y = np.where(~np.isnan(dL_comp/z_comp), dL_comp/z_comp, np.nan) # help
        ax.plot(z_comp, dL_comp/z_comp, alpha=.7, label=r"$d_L/z$")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"Gpc")

        ax.legend()

        if self.save:
            save(self.subdir + "lum_distance_vs_redshift")
        if self.show:
            plt.show()

    def OmegaM_OmegaLambda(self, Omega_M0, Omega_Lambda0, chi2):
        
        # tmp
        org_x = 0.05+0.267
        org_y = 1-(org_x+5.50896e-05)

        chi2_min = np.min(chi2)
        sigma1, sigma2 = 3.53, 8.02 # sigma levels
        thresh1 = (chi2 < chi2_min + sigma1)
        thresh2 = (chi2 < chi2_min + sigma2)

        xx, yy, zz = r"$\Omega_{\mathrm{m}0}$", r"$\Omega_{\Lambda 0}$", r"$\chi2$"
        df = pd.DataFrame({xx:Omega_M0, yy:Omega_Lambda0, zz:chi2})

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

        second = ax.scatter(Omega_M0[thresh2], Omega_Lambda0[thresh2], color=c2, edgecolor=c2, label=r"$2\sigma$ level", **kwargs)
        first = ax.scatter(Omega_M0[thresh1], Omega_Lambda0[thresh1], color=c1, edgecolor=c1, label=r"$1\sigma$ level", **kwargs)
        
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
        mu, sigma = norm.fit(Omega_Lambda0)
        OmL_pdf = norm.pdf(LL, mu, sigma)
        pdf, = ax_L.plot(OmL_pdf, LL)
        # leg = ax_L.legend([pdf], [r"$\mathcal{N}(\sigma\!=\!%.2f, \mu\!=\!%.2f)$" %(sigma, mu)], fontsize=14, markerscale=0.7, markerfirst=False)
        # leg = ax_L.legend([pdf], [r"$\mathcal{N}(%.2f,$"%sigma+"\n"+ r"$\, %.2f)$" %mu], borderaxespad=.1, **leg_kw)
        ax_L.legend([pdf], [tex.normal(mu, sigma, splitline=True)], borderaxespad=.05, **leg_kw)
        # text = leg.get_texts()[0]
        # text.set_rotation(-90)
        # reg = leg.get_frame()
        # reg.set_angle(-90)
        # ax_L.legend()
       
       
       # add our original point
        c = COLOURS[4]
        kw = dict(lw=5, ms=16)
        org_p, = ax.plot(org_x, org_y, "+", c=c, label="fiducial", **kw)
        # add new point
        new_p, = ax.plot(Omega_M0[np.nanargmin(chi2)], Omega_Lambda0[np.nanargmin(chi2)], "+", c=COLOURS[5], label="new best-fit", **kw)

        # ax.annotate(r"$(\Omega_{\mathrm{m}0}, \Omega_{\Lambda 0}) = (%.3f, %.3f)$" %(org_x, org_y), (org_x, org_y), (org_x-0.1, org_y-0.35), arrowprops={"fc":c, "shrink":0.1}, color=c)
        
        ax.legend(handles=[first, second, org_p, new_p], markerscale=1.7)
        # ax.legend(reverse=True, markerscale=2.)



        fig.tight_layout = True   

        if self.save:
            save(self.subdir + "omega_phase_space")
        if self.show:
            plt.show()

    def OmegaK_PDF(self, Omega_K0):
        from scipy.stats import norm

        fig, ax = plt.subplots(figsize=(10,4))
        sns.histplot(Omega_K0, ax=ax, stat="density", alpha=.66)
        
        K_lims = (-1., 1.)
        # x0, x1 = ax.get_xlim()
        ax.set_xlabel(r"$\Omega_{\mathrm{K}0}$")

        mu, sigma = norm.fit(Omega_K0)
        xx = np.linspace(*K_lims, 100)
        ax.plot(xx, sp.stats.norm.pdf(xx, mu, sigma), label=tex.normal(mu, sigma))
        ax.axvline(0.0, color=COLOURS[4], lw=1.4, label="fiducial")
        # kw = dict(lw=8, ms=22)
        # org_p, = ax.plot(0, 0, "+", c=COLOURS[4], label="fiducial", **kw)
        ax.legend()

        if self.save:
            save(self.subdir + "curvature_pdf")
        if self.show:
            plt.show()


    def HubblePDF(self, H0):
        # thresh = chi2<np.min(chi2)+3.53
        fig, ax = plt.subplots(figsize=(10,4))
        # ax.hist(H0[thresh], 60)
        from scipy.stats import norm
        #
        sns.histplot(H0, ax=ax, stat="density", alpha=.76)
        # x0, x1 = ax.get_xlim()
        h_lims = (0.65, .75)
        ax.set_xlabel(r"$h$")
        # sigma = np.std(H0)
        # mu = np.mean(H0)
        mu, sigma = norm.fit(H0)
        xx = np.linspace(*h_lims, 100)
        # pdf = sp.stats.norm.pdf(H0)
        ax.plot(xx, sp.stats.norm.pdf(xx, mu, sigma), label=tex.normal(mu, sigma))
        ax.axvline(0.67, color=COLOURS[4], lw=1.4, label="fiducial")
        # kw = dict(lw=8, ms=22)
        # org_p, = ax.plot(0.67, 0, "+", c=COLOURS[4], label="fiducial", **kw)
        ax.legend()
        


        if self.save:
            save(self.subdir + "hubble_pdf")
        if self.show:
            plt.show()



    



