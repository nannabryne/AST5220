from utils import *
import plotsrc_milestone_one as PLOT

const = ConstantsAndUnits()

tex = LaTeX()

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
    dic["ugamma"]   = -3*dic["Theta1"]

    df = pd.DataFrame(dic)
    
    return df

# testing ...


k_list = [0.001, 0.01, 0.1]
df_list = []

for k in k_list:
    pert = read_ASCII(f"perturbations_k{k}.txt")
    df_list.append(create_dataframe(pert))


# reprodusere farger ...
k_colour = ColourCycles(["dodgerblue", "orange", "forestgreen"])

k_label = lambda k: tex("k=%.3f" %k + tex.unit(tex.inv("Mpc")))

k_dict = dict(label     = [k_label(k) for k in k_list], 
              colour    = k_colour(),
              value     = k_list,
              df        = df_list)



fig1, ax1 = plt.subplots()

ax = ax1

for i, df in enumerate(df_list):
    ax.plot(df["x"], df["uc"], c=k_colour[i], alpha=.7, label=k_label(k_list[i]))
    ax.plot(df["x"], df["ub"], c=k_colour[i], alpha=.7, ls="--")

ax.set_yscale("log")
ax.set_title(tex("u" + tex.ped("{c,b}")))
ax.legend()
ax.set_ylim(1e-6,1e2)






fig2, ax2 = plt.subplots()

ax = ax2
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["deltac"], c=k_colour[i], alpha=.7, label=tex("k=%.3f" %k_list[i] + tex.unit(tex.inv("Mpc"))))
    ax.plot(df["x"], df["deltab"], c=k_colour[i], alpha=.7, ls="--")

ax.set_yscale("log")
ax.set_title(tex("\delta" + tex.ped("{c,b}")))
ax.legend()
ax.set_ylim(1e-2,1e6)




fig3, ax3 = plt.subplots()

ax = ax3
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["Theta0"], c=k_colour[i], alpha=.7, label=tex("k=%.3f" %k_list[i] + tex.unit(tex.inv("Mpc"))))

ax.set_title(tex(tex.Theta + "_0"))
ax.legend()
ax.set_ylim(-1,1)


fig4, ax4 = plt.subplots()

ax = ax4
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["Theta1"], c=k_colour[i], alpha=.7, label=tex("k=%.3f" %k_list[i] + tex.unit(tex.inv("Mpc"))))

ax.set_title(tex(tex.Theta + "_1"))
ax.legend()
ax.set_ylim(-1,1)


fig5, ax5 = plt.subplots()

ax = ax5
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["Phi"], c=k_colour[i], alpha=.7, label=tex("k=%.3f" %k_list[i] + tex.unit(tex.inv("Mpc"))))

ax.set_title(tex.Phi)
ax.legend()
ax.set_ylim(0,1)




ax_list = [ax1, ax2, ax3, ax4, ax5]

for ax in ax_list:
    # ax.set_xlim(-18, -8.3)
    ax.set_xlabel(tex.x)
    ax.axvline(-8.3, color="k", ls="-")
    pass




plt.show()


