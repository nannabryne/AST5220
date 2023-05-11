from utils import *
import plotsrc_milestone_four as PLOT

const = ConstantsAndUnits()

tex = LaTeX()




C_ell_data = read_ASCII("cells.txt")
ell, Dell = C_ell_data[:,0], C_ell_data[:,1]

df_Dell = pd.DataFrame(dict(ell=C_ell_data[:,0], D_ell=C_ell_data[:,1]))





power_data = read_ASCII("power.txt")
k = power_data[:,0]
P = power_data[:,1]
Theta_ell = power_data[:,2:].T

ells = np.array([6, 100, 200, 500, 1000], dtype="int")

df_power = pd.DataFrame(dict(k=power_data[:,0], P=power_data[:,1]))


# Theta_dict = {f"\Theta_{ells[i]:d}": Theta_ell[i] for i in range(len(ells))}s
Theta_dict = {f"{ells[i]:d}": Theta_ell[i] for i in range(len(ells))}

df_transfer = pd.DataFrame(dict(k=power_data[:,0], **Theta_dict))


#   SORT OBSERVATIONAL DATA


Planck_data = read_ASCII("planck_low_ell_TT.txt")
galaxy_data = read_ASCII("galaxy_matter_power.txt")
WMAP_data   = read_ASCII("wmap_matter_power.txt")

df_obs_CMB = pd.DataFrame(dict(ell=Planck_data[:,0], 
                              D_ell=Planck_data[:,1],
                              err_up=Planck_data[:,2],
                              err_down=Planck_data[:,3]))

df_obs_matter1 = pd.DataFrame(dict(k=galaxy_data[:,0], 
                      P=galaxy_data[:,1],
                      err_up=galaxy_data[:,2]))

df_obs_matter2 = pd.DataFrame(dict(k=WMAP_data[:,0], 
                      P=WMAP_data[:,1],
                      err_up=WMAP_data[:,2]-WMAP_data[:,1]))


df_obs_matter = pd.concat([df_obs_matter1, df_obs_matter2], ignore_index=True)
df_obs_matter["err_down"] = df_obs_matter["err_up"] 



#   DO ALL THE PLOTTING

# PLOT.TransferFunction(df_transfer)
PLOT.CMBPowerSpectrum(df_Dell, df_obs_CMB)
PLOT.MatterPowerSpectrum(df_power, df_obs_matter)

plt.show()