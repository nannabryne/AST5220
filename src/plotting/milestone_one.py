from utils import *
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()



Gyr2sec =  1e9 * 365*24*60*60
sec2Gyr = 1/Gyr2sec

Hubble_conv_fac = const.Mpc/(100*const.km/const.s)
inv_Hubble_conv_fac = 1/Hubble_conv_fac

Ga2sec =  1e9 * 365*24*60*60
sec2Ga = 1/Ga2sec

Gpc2m = const.Mpc*1e3
m2Gpc = 1/Gpc2m



INFO = dict(x_eq=-8.132, x_acc=-0.4869, x_Lambda=-0.2558, x_0=0)


'''
Read basic data
'''


data = read_ASCII("background_cosmology")


Hp_plus_derivatives = dict(
    Hp      = data[:,3]*Hubble_conv_fac,
    dHpdx   = data[:,4]*Hubble_conv_fac,
    ddHpdxx = data[:,5]*Hubble_conv_fac
)

Omegas = dict(
    OmegaM  = data[:,6]+data[:,7],   # total matter
    OmegaL  = data[:,8],             # vacuum
    OmegaR  = data[:,9]+data[:,10]   # total radiation
)

dL = dict(dL=data[:,12]*m2Gpc, z=np.exp(-data[:,0])-1)

cosmo_df = pd.DataFrame(dict(
    x       = data[:,0],
    eta_c_org = data[:,1]/const.c, 
    eta_c   = data[:,1]/const.c * sec2Ga,
    t       = data[:,2] * sec2Ga,
    Hp_org  = data[:,3],
    **Hp_plus_derivatives,
    **Omegas,
    **dL
))


xx = cosmo_df["x"]
idx0 = np.nanargmin(np.abs(xx))
todays = dict(
    OmegaR  = cosmo_df["OmegaR"][idx0],
    OmegaM  = cosmo_df["OmegaM"][idx0],
    OmegaL  = cosmo_df["OmegaL"][idx0],
    OmegaK  = 0.0,
    H       = cosmo_df["Hp"][idx0]*Hubble_conv_fac*const.H0_over_h,
    t       = cosmo_df["t"][idx0],
    eta_c   = cosmo_df["eta_c"][idx0]
)



x_rad = xx[(xx<=INFO["x_eq"])]
x_mat = xx[(xx>INFO["x_eq"]) & (xx<=INFO["x_Lambda"])]
x_vac = xx[(xx>INFO["x_Lambda"])]


ones_rad = np.ones(len(x_rad))
ones_mat = np.ones(len(x_mat))
ones_vac = np.ones(len(x_vac))

analytical_predictions = dict(
    x = xx,
    x_rad = x_rad,
    x_mat = x_mat,
    x_vac = x_vac,
    Hubble = np.append(np.append(np.sqrt(todays["OmegaR"])*np.exp(-x_rad), np.sqrt(todays["OmegaM"])*np.exp(-0.5*x_mat)),  np.sqrt(todays["OmegaL"])*np.exp(x_vac)),
    Hubble_der = np.append(np.append(-ones_rad, -0.5*ones_mat), ones_vac),
    Hubble_doubleder = np.append(np.append(ones_rad, -0.25*ones_mat), ones_vac),
    eta_Hubble = np.append(np.append(ones_rad, -2*ones_mat), np.nan*ones_vac)
)


H0 = todays["H"]
analytical_rad = dict(
    x = x_rad,
    Hubble = np.sqrt(todays["OmegaR"])*np.exp(-x_rad) * H0,
    Hubble_der = -np.ones(len(x_rad)),
    Hubble_doubleder = np.ones(len(x_rad)),
    eta_Hubble = np.ones(len(x_rad)) 
)

analytical_mat = dict(
    x = x_mat,
    Hubble = np.sqrt(todays["OmegaM"])*np.exp(-0.5*x_mat) * H0,
    Hubble_der = -0.5*np.ones(len(x_mat)),
    Hubble_doubleder = 0.25*np.ones(len(x_mat)),
    eta_Hubble = 2*np.ones(len(x_mat)) 
)

analytical_vac = dict(
    x = x_vac,
    Hubble = np.sqrt(todays["OmegaL"])*np.exp(x_vac) * H0,
    Hubble_der = np.ones(len(x_vac)),
    Hubble_doubleder = np.ones(len(x_vac)),
    eta_Hubble = np.nan*np.ones(len(x_vac)) 
)

analytical_predictions = pd.DataFrame(dict(radiation=analytical_rad, matter=analytical_mat, vacuum=analytical_vac))




'''
Do all the plotting prior to fitting
'''

PLOT.DensityParameters(cosmo_df, INFO)

PLOT.HubbleDerivatives(cosmo_df, INFO, analytical_predictions)

PLOT.ConformalTime_HubbleParameter(cosmo_df, INFO, analytical_predictions)

PLOT.HubbleParameter(cosmo_df, INFO, analytical_predictions)

PLOT.TimeMeasures(cosmo_df, INFO, todays)


'''
Read MCMC data
'''


obs_data = read_ASCII("supernovadata", 1)   
z_obs = obs_data[:,0]   
dL_obs = obs_data[:,1]  
err_obs = obs_data[:,2]

obs_df = pd.DataFrame(dict(
    z  =  obs_data[:,0],
    dL = obs_data[:,1],
    err = obs_data[:,2]
))

fit_data = read_ASCII("mcmc_fitting", 200)
result_mcmc_data = read_ASCII("revised_background_cosmology")

mcmc_df = pd.DataFrame(dict(
    chi2 = fit_data[:,0],
    H0 = fit_data[:,1]*const.H0_over_h*Hubble_conv_fac,
    OmegaM0 = fit_data[:,2],
    OmegaK0 = fit_data[:,3]
))
mcmc_df["OmegaL0"] = 1 - (mcmc_df["OmegaM0"] + mcmc_df["OmegaK0"])

newcosmo_df = pd.DataFrame(dict(
    dL = result_mcmc_data[:,12]*m2Gpc
))



'''
Do all the plotting posterior to fitting
'''


PLOT.LuminosityDistance(obs_df, cosmo_df, newcosmo_df)

PLOT.Omega_PhaseSpace(mcmc_df, todays)

PLOT.Curvature_PDF(mcmc_df, todays)

PLOT.Hubble_PDF(mcmc_df, todays)

plt.show()



