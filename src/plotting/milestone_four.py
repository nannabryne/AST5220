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



PLOT.TransferFunction(df_transfer)
PLOT.CMBPowerSpectrum(df_Dell)
PLOT.MatterPowerSpectrum(df_power)

plt.show()