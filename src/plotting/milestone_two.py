from utils import *
import plotsrc_milestone_two as PLOT

const = ConstantsAndUnits()

#   get data:

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



PLOT.ElectronFraction(df, INFO)
# PLOT.BaryonTemperature(df, INFO)
PLOT.OpticalDepth(df, INFO)
PLOT.VisibilityFunction(df, INFO)

plt.show()

