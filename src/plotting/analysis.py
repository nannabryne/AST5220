from utils import *
const = ConstantsAndUnits()
tex = LaTeX()



'''
*********************************
MILESTONE I: Background cosmology
*********************************
'''


def Milestone1():

    import plotsrc_milestone_one as PLOT1

    Hubble_conv_fac = const.Mpc/(100*const.km/const.s)
    inv_Hubble_conv_fac = 1/Hubble_conv_fac


    Ga2sec =  1e9 * 365*24*60*60
    sec2Ga = 1/Ga2sec

    Gpc2m = const.Mpc*1e3
    m2Gpc = 1/Gpc2m

    INFO = dict(x_eq=-8.132, x_acc=-0.4869, x_Lambda=-0.2558, x_0=0)



    #   Read basic data


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


    #   Do all the plotting prior to fitting

    PLOT1.DensityParameters(cosmo_df, INFO)
    PLOT1.HubbleDerivatives(cosmo_df, INFO, analytical_predictions)
    PLOT1.ConformalTime_HubbleParameter(cosmo_df, INFO, analytical_predictions)
    PLOT1.HubbleParameter(cosmo_df, INFO, analytical_predictions)
    PLOT1.TimeMeasures(cosmo_df, INFO, todays)



    #   Read MCMC data


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

    #   Do all the plotting posterior to fitting

    PLOT1.LuminosityDistance(obs_df, cosmo_df, newcosmo_df)
    PLOT1.Omega_PhaseSpace(mcmc_df, todays)
    PLOT1.Curvature_PDF(mcmc_df, todays)
    PLOT1.Hubble_PDF(mcmc_df, todays)




'''
***********************************
MILESTONE II: Recombination history
***********************************
'''


def Milestone2():

    import plotsrc_milestone_two as PLOT2


    #   Get data

    rec = read_ASCII("recombination")
    rec_Saha = read_ASCII("recombination_Saha_only")

    x_rec = -6.985
    Xe_fo = 2.026e-4
    INFO = dict(x_rec=-6.985, Xe_fo=2.026e-4)


    x, Xe = rec[:,0], rec[:,1]
    x_Saha, Xe_Saha = rec_Saha[:,0], rec_Saha[:,1]


    # dunno
    tau = np.where((x>0)|(x<-12), np.nan, rec[:,3])
    dtaudx = np.where(np.isnan(tau), np.nan, rec[:,4])
    dtaudx[0] = np.nan
    dtaudx[-1] = np.nan
    ddtaudxx = np.where(np.isnan(dtaudx), np.nan, rec[:,5])
    ddtaudxx[-1] = np.nan
    ddtaudxx[-2] = np.nan


    tau_plus_derivatives = dict(
        tau=tau,
        dtaudx=dtaudx,
        ddtaudxx=ddtaudxx
    )

    gt_plus_derivatives = dict(
        gt = np.where(np.isnan(tau), np.nan, rec[:,6]),
        dgtdx = np.where(np.isnan(tau), np.nan, rec[:,7]),
        ddgtdxx = np.where(np.isnan(tau), np.nan, rec[:,8])
    )

    df = pd.DataFrame(dict(
        x=rec[:,0], 
        Xe=rec[:,1], 
        x_Saha=rec_Saha[:,0], 
        Xe_Saha=rec_Saha[:,1],
        Tb=2.7255*np.exp(-rec[:,0]),
        **tau_plus_derivatives,
        **gt_plus_derivatives
    ))


    #   Do all the plotting

    PLOT2.ElectronFraction(df, INFO)
    # PLOT2.BaryonTemperature(df, INFO)
    PLOT2.OpticalDepth(df, INFO)
    PLOT2.VisibilityFunction(df, INFO)





'''
**********************************
MILESTONE III: Growth of structure
**********************************
'''

def Milestone3():

    import plotsrc_milestone_three as PLOT3

    x_eq = -8.132   # RM-equality
    x_rec = -6.985  # recombination onset

    PLOT3.INIT(x_eq, x_rec)

    def create_dataframe(data):
        dic = dict(
            x           = data[:,0],
            deltac      = data[:,1],
            deltab      = data[:,2],
            uc          = data[:,3],
            ub          = data[:,4],
            Theta0      = data[:,5],
            Theta1      = data[:,6],
            Theta2      = data[:,7],
            Phi         = data[:,8],
            Psi         = data[:,9]
        )
        dic["deltagamma"]   = 4*dic["Theta0"]
        dic["ugamma"]       = -3*dic["Theta1"]

        df = pd.DataFrame(dic)
        
        return df


    k_list = [0.001, 0.01, 0.1]
    df_list = []

    for k in k_list:
        pert = read_ASCII(f"perturbations_k{k}.txt")
        df_list.append(create_dataframe(pert))


    k_colour = ColourCycles()


    k_label = lambda k: tex("k=%.3f" %k + tex.unit(tex.inv("Mpc")))

    k_dict = dict(label     = [k_label(k) for k in k_list], 
                val_str   = [str(k) for k in k_list],
                colour    = k_colour(),
                value     = k_list,
                value_SI  = [k/const.Mpc for k in k_list])


    PLOT3.VelocityPerturbations(df_list, k_dict)
    PLOT3.DensityPerturbations(df_list, k_dict)
    PLOT3.PhotonQuadrupole(df_list, k_dict)
    PLOT3.GravitationalPotential(df_list, k_dict)
    # PLOT3.SanityChecks(df_list)




'''
****************************************
MILESTONE IV: CMB & MATTER POWER SPECTRA
****************************************
'''

def Milestone4():

    import plotsrc_milestone_four as PLOT4

    # INSERT CODE HERE WHEN FINISED







if __name__ == "__main__":


    # n_milestones = len(sys.argv[1:])
    

    milestoneI, milestoneII, milestoneIII, milestoneIV = False, False, False, False 

    milestoneIV = True

    if milestoneI:
        Milestone1()
    if milestoneII:
        Milestone2()
    if milestoneIII:
        Milestone3()
    if milestoneIV:
        Milestone4()
    
    plt.show()


