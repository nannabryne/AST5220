from utils import *
import plotsrc_milestone_three as PLOT

const = ConstantsAndUnits()

tex = LaTeX()

x_eq = -8.132   # RM-equality
x_rec = -6.985  # recombination onset

PLOT.INIT(x_eq, x_rec)

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




PLOT.VelocityPerturbations(df_list, k_dict)

PLOT.DensityPerturbations(df_list, k_dict)

PLOT.PhotonQuadrupole(df_list, k_dict)

PLOT.GravitationalPotential(df_list, k_dict)


# PLOT.SanityChecks(df_list)


plt.show()











# Delete me:
data_cosmo = read_ASCII("background_cosmology")


x_cosmo, eta, Hp = data_cosmo[:,0], data_cosmo[:,1], data_cosmo[:,3]

Hp[x_cosmo>0] = np.nan

for i in range(len(k_list)):
    hor = np.abs(eta-1/k_dict["value_SI"][i]) 

    # hor = np.abs(Hp/const.c-k_dict["value_SI"][i]) 
    x_hor = x_cosmo[np.nanargmin(hor)]
    k = k_dict["val_str"][i]


    print(f"k = {k} Mpc^(-1) crosses horizon at x = {x_hor:.3f}")





# fig, ax = plt.subplots()


# for k in k_list:
#     ax.plot(x_cosmo, k/const.Mpc*eta)

# ax.axhline(1,color="k")
# ax.set_yscale("log")
# plt.show()



# Ga2sec =  1e9 * 365*24*60*60
# sec2Ga = 1/Ga2sec

# ck_inv_list = []
# for k in k_list:
#     ck_inv = 1/(k/const.Mpc * const.c) * sec2Ga
#     ck_inv_list.append(ck_inv)


# print(ck_inv_list)

# print(1/(1e-4/const.Mpc * const.c) * sec2Ga)
# print(1/(1e-5/const.Mpc * const.c) * sec2Ga)

# eta0  = 46.32*const.c * Ga2sec  # eta in sec
# print(1e-4/const.Mpc*eta0)




# print(1/(eta0/const.Mpc))

