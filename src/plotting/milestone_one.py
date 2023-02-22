from utils import *
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()



# var_labels = {
#     "x":            r"$x$",
#     "eta":          r"$\eta$",
#     "t":            r"$t$",
#     "Hp":           r"$\mathcal{H}$", 
#     "dHpdx":        r"$\mathcal{H}'$",
#     "ddHpdxx":      r"$\mathcal{H}''$",
#     "OmegaM":       r"$\Omega_\mathrm{M}$",
#     "OmegaLambda":  r"$\Omega_\Lambda$",
#     "OmegaR":       r"$\Omega_\mathrm{R}$"
# }

# SI_unit_labels = {
#     "eta":  "m",
#     "t":    r"s",
#     "Hp":   r"$\mathrm{s}^{-1}$"
# }
# SI_unit_labels["dHdx"] = SI_unit_labels["Hp"]
# SI_unit_labels["ddHdxx"] = SI_unit_labels["Hp"]

SI2sensible = {
    "eta":  1/const.c,
    "Hp":   const.Mpc/(100*const.km/const.s),
    "dL":   1/const.Mpc * 1e-3
}
SI2sensible["dHpdx"] = SI2sensible["Hp"]
SI2sensible["ddHpdxx"] = SI2sensible["Hp"]

Gyr2sec =  1e9 * 365*24*60*60
sec2Gyr = 1/Gyr2sec

Hubble_conv_fac = const.Mpc/(100*const.km/const.s)
inv_Hubble_conv_fac = 1/Hubble_conv_fac

# sensible_unit_labels = {
#     "eta":  "s",
#     "t":    "s",
#     "Hp":   "100 km/s / Mpc",
#     "dL":   "Gpc"
# }


class SimpleCosmo:

    def __init__(self, filename, use_sensible_units=True):
        data = read_ASCII(filename)

        self.x              = data[:,0]
        self.eta            = data[:,1]
        self.t              = data[:,2]
        self.Hp             = data[:,3]
        self.dHpdx          = data[:,4]
        self.ddHpdxx        = data[:,5]
        self.OmegaM         = data[:,6]+data[:,7]    # total matter
        self.OmegaLambda    = data[:,8]
        self.OmegaR         = data[:,9]+data[:,10]   # total radiation

        
        self.z = self.compute_redshift()
        self.dL = self.compute_luminosity_distance()

        self.locate_milestones()

        self.plot = PLOT.MilestoneI()
        self.plot.DefineBasics(self.x, self.idx)

        self.eta_c = self.eta/const.c
        self.t_Gyr = self.t*sec2Gyr
        self.eta_c_Gyr = self.eta_c*sec2Gyr

        self.dL_Gpc = self.dL/const.Mpc*1e-3
        

    def locate_milestones(self):
        jump = 50   # avoid the inflation era
        self.idx = {
            "RM":   np.argmin(np.abs(self.OmegaM-self.OmegaR)[jump:])+jump,         # matter-radiation equality
            "MDE":  np.argmin(np.abs(self.OmegaM-self.OmegaLambda)[jump:])+jump,    # matter-dark energy transition
            "acc":  np.argmin(np.abs(self.dHpdx)[jump:])+jump,                      # universe starts accelerating
            "0":    np.argmin(np.abs(self.x)[jump:])+jump                           # today
        }
    
    def compute_luminosity_distance(self):
        Chi = self.eta[np.argmin(np.abs(self.x))] - self.eta
        dL = Chi*np.exp(-self.x) 
        return dL

    def compute_redshift(self):
        z = np.exp(-self.x)-1
        return z

    
    def make_table(self):
        print("                   x      z       t    ")
        print("                                [Gyr]")
        RM_idx = self.idx["RM"]
        for m in self.idx.keys():
            idx = self.idx[m]
            print(f"{m:10s}: {self.x[idx]:6.2f}  {self.z[idx]:6.2f}  {self.t_Gyr[idx]:10.6f}")


    def plot_density_params(self):
        self.plot.DensityParameters(self.OmegaM, self.OmegaLambda, self.OmegaR)

    def plot_misc_functions(self):
        self.plot.MiscellaneousFunctions(self.eta_c, self.Hp, self.dHpdx, self.ddHpdxx)

    def plot_time_measures(self):
        self.plot.TimeMeasures(self.t_Gyr, self.eta_c_Gyr)

    def plot_lum_dist(self, supernova_data):
        self.plot.LuminosityDistance(supernova_data.z_obs, supernova_data.dL_obs, supernova_data.err_obs, self.z, self.dL_Gpc)


class SupernovaFitting:

    def __init__(self, data_filename, result_filename):
        obs_data = read_ASCII(data_filename, 1)   
        self.z_obs = obs_data[:,0]   
        self.dL_obs = obs_data[:,1]  
        self.err_obs = obs_data[:,2]

        fit_data = read_ASCII(result_filename, 100)
        self.Chi2 = fit_data[:,0]
        self.H0 = fit_data[:,1]*const.H0_over_h
        self.OmegaM0 = fit_data[:,2]
        self.Omegak0 = fit_data[:,3]

        self.Omegagamma0 = self.compute_photon_density()
        self.OmegaLambda0 = self.compute_dark_energy_density()

        self.plot = PLOT.MilestoneI()

    def compute_photon_density(self):
        Omegagamma = 2*np.pi**2/30 * (const.k_b*2.755)**4 * 8*np.pi*const.G/(3*self.H0*self.H0)
        return Omegagamma
    
    def compute_dark_energy_density(self):
        OmegaLambda = 1 - self.OmegaM0 - self.Omegak0 - self.Omegagamma0
        return OmegaLambda

    def plot_Omega_phase_space(self):
        self.plot.OmegaM_OmegaLambda(self.OmegaM0, self.OmegaLambda0, self.Chi2)

    def plot_Hubble_pdf(self):
        self.plot.HubblePDF(self.H0*Hubble_conv_fac, self.Chi2)






test = SimpleCosmo("cosmology")
test.make_table()
test.plot_density_params()
test.plot_misc_functions()
test.plot_time_measures()
test2 = SupernovaFitting("supernovadata", "mcmc_fitting")
test.plot_lum_dist(test2)

test2.plot_Omega_phase_space()
test2.plot_Hubble_pdf()

