from utils import *




class MilestoneI:

    def __init__(self):
        self.SaveMode()
        self.ShowMode()

        self.subdir = "milestone1/"

    def DefineBasics(self, x_array, indices):
        #   Define basics
        self.x = x_array
        self.idx = indices

    def _mark_milestones(self, ax):
        # FIXME

        ax.axvline(self.x[self.idx["RM"]],  ls="--", lw=1, color="orangered")
        ax.axvline(self.x[self.idx["MDE"]], ls="--", lw=1, color="olive")
        ax.axvline(self.x[self.idx["acc"]], ls="--", lw=1, color="firebrick")

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
        fig, ax = plt.subplots()

        ax.plot(self.x, Omega_M, label=r"$\Omega_\mathrm{M}(x)$")
        ax.plot(self.x, Omega_Lambda, label=r"$\Omega_\Lambda(x)$")
        ax.plot(self.x, Omega_R, label=r"$\Omega_\mathrm{R}(x)$")

        ax.legend()
        ax.set_xlabel(r"$x$")

        self._mark_milestones(ax)

        if self.save:
            save(self.subdir + "density_params")
        if self.show:
            plt.show()

    
    def MiscellaneousFunctions(self, eta_c, Hp, dHp_dx, ddHp_dxx):
        fig, ax = plt.subplots()

        ax.plot(self.x, ddHp_dxx/Hp, label=r"$\mathcal{H}''(x)/\mathcal{H}(x)$")
        ax.plot(self.x, dHp_dx/Hp,   label=r"$\mathcal{H}'(x)/\mathcal{H}(x)$")
        ax.plot(self.x, eta_c*Hp,    label=r"$\eta(x)\mathcal{H}(x)/c$")

        ax.legend()
        ax.set_xlabel(r"$x$")

        self._mark_milestones(ax)

        if self.save:
            save(self.subdir + "misc_Hubble")
        if self.show:
            plt.show()

    
    def HubbleParameter(self, Hp):

        fig, ax = plt.subplots()

        ax.plot(self.x, Hp, label=r"$\mathcal{H}(x)$")

        ax.set_ylabel(r"100 km/s / Mpc")
        ax.set_xlabel(r"$x$")
        # ax.set_xlim(np.min(x), 0)

        ax.legend()

        if self.save:
            save(self.subdir + "hubble_param")
        if self.show:
            plt.show()



    def TimeMeasures(self, t, eta_c):
        fig, ax = plt.subplots()

        ax.plot(self.x, t,     label=r"$t(x)$")
        ax.plot(self.x, eta_c, label=r"$\eta(x)/c$")

        ax.set_ylabel("Gyrs")
        ax.set_xlabel(r"$x$")
        ax.set_title("hei", fontweight="bold")

        ax.legend()

        if self.save:
            save(self.subdir + "time_measures") 
        if self.show:
            plt.show()


    
    def LuminosityDistance(self, z_obs, dL_obs, err_obs, z_comp, dL_comp):
        fig, ax = plt.subplots()
    
        ax.errorbar(z_obs, dL_obs, err_obs, elinewidth=1.1, capsize=2, linestyle="", marker="o", ms=4, label=r"$d_L^\mathrm{obs}(z)$")
        # ax.plot(z_obs, dL_obs, 'o')#,         c="orangered")
        ax.set_xscale("log")
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # ax.plot(10000, 0)


        ax.plot(z_comp, dL_comp, alpha=.5, label=r"$d_L(z)$")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"$\mathrm{[Gpc]}$")

        ax.legend()

        if self.save:
            save(self.subdir + "lum_distance_vs_redshift")
        if self.show:
            plt.show()


    def OmegaM_OmegaLambda(self, Omega_M0, Omega_Lambda0, chi2):

        thresh = chi2<np.min(chi2)+3.53

        # fig, ax = plt.subplots(figsize=(15,6))
        # fig.set_facecolor("#313332")
        # ax.patch.set_facecolor("grey")
        # ax.scatter(Omega_M0[thresh], Omega_Lambda0[thresh], c=chi2[thresh], alpha=0.6, s=5)

        # ax.plot((0, 1), (1, 0), ls='--', c='k')
        # # ax.set_aspect("equal")
        # ax.set_xlim(0, 1)
        # ax.set_ylim(0,1.4)

        # ax.set_xlabel(r"$\Omega_{\mathrm{M}0}$")
        # ax.set_ylabel(r"$\Omega_{\Lambda0}$")

        xx, yy, zz = r"$\Omega_{\mathrm{M}0}$", r"$\Omega_{\Lambda 0}$", r"$\chi2$"

        df = pd.DataFrame({xx:Omega_M0[thresh], yy:Omega_Lambda0[thresh], zz:chi2[thresh]})
        # sns.warn_singular=False
        
        g = sns.JointGrid(data=df, x=xx, y=yy)
        fig = g.figure
        g.fig.set_figwidth(10)
        g.fig.set_figheight(6)

        g.plot_joint(sns.scatterplot, alpha=.5, legend=False)
        g.plot_joint(sns.kdeplot, levels=3)
        g.plot_marginals(sns.kdeplot)
        

        fig.tight_layout = True    


        if self.save:
            save(self.subdir + "omega_phase_space")
        if self.show:
            plt.show()


    def HubblePDF(self, H0, chi2):
        thresh = chi2<np.min(chi2)+3.53
        fig, ax = plt.subplots()
        # ax.hist(H0[thresh], 60)
        #
        sns.histplot(H0[thresh], kde=True, ax=ax)
        ax.set_xlabel(r"$H_0$ [100 km/s / Mpc]")
        ax.set_title("Posteriori", fontfamily="sans-serif")


        if self.save:
            save(self.subdir + "hubble_pdf")
        if self.show:
            plt.show()



    

