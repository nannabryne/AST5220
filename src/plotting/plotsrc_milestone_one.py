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

        ax.legend()

        if self.save:
            save(self.subdir + "time_measures") 
        if self.show:
            plt.show()


    
    def LuminosityDistance(self, z_obs, dL_obs, err_obs, z_comp, dL_comp):
        fig, ax = plt.subplots()
    
        ax.errorbar(z_obs, dL_obs, err_obs, c="orangered", label=r"$d_L^\mathrm{obs}(z)$")
        ax.plot(z_obs, dL_obs, 'o',         c="orangered", label=r"$d_L(z)$")
        ax.set_xscale("log")
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()


        ax.plot(z_comp, dL_comp)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"$\mathrm{[Gpc]}$")

        ax.legend()

        if self.save:
            save(self.subdir + "lum_distance_vs_redshift")
        if self.show:
            plt.show()


    def OmegaM_OmegaLambda(self, Omega_M0, Omega_Lambda0, Chi2):

        thresh = Chi2<np.min(Chi2)+3.53

        fig, ax = plt.subplots()
        ax.scatter(Omega_M0[thresh], Omega_Lambda0[thresh], c=Chi2[thresh], alpha=0.6, s=5)

        ax.plot((0, 1), (1, 0), ls='--', c='k')
        # ax.set_aspect("equal")
        ax.set_xlim(0, 1)
        ax.set_ylim(0,1.4)

        ax.set_xlabel(r"$\Omega_{\mathrm{M}0}$")
        ax.set_ylabel(r"$\Omega_{\Lambda0}$")

        if self.save:
            save(self.subdir + "omega_phase_space")
        if self.show:
            plt.show()


    def HubblePDF(self, H0, Chi2):
        thresh = Chi2<np.min(Chi2)+3.53
        fig, ax = plt.subplots()
        ax.hist(H0[thresh], 60)
        ax.set_xlabel(r"$H_0$ [100 km/s / Mpc]")

        if self.save:
            save(self.subdir + "hubble_pdf")
        if self.show:
            plt.show()



    

