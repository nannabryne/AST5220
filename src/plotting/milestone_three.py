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

    df = pd.DataFrame(dic)
    
    return df

# testing ...


k_list = [0.001, 0.01, 0.1]
df_list = []

for k in k_list:
    pert = read_ASCII(f"perturbations_k{k}.txt")
    df_list.append(create_dataframe(pert))

fig, ax = plt.subplots()

c = ColourCycles()
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["uc"], c=c[i], alpha=.7, label=tex("k=%f/\mathrm{Mpc}" %k_list[i]))
    ax.plot(df["x"], df["ub"], c=c[i], alpha=.7, ls="--")

ax.set_yscale("log")
ax.set_title(tex("u" + tex.ped("{b,c}")))
ax.legend()

fig, ax = plt.subplots()

c = ColourCycles()
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["deltac"], c=c[i], alpha=.7, label=tex("k=%.3f/\mathrm{Mpc}" %k_list[i]))
    ax.plot(df["x"], df["deltab"], c=c[i], alpha=.7, ls="--")

ax.set_yscale("log")
ax.set_title(tex("\delta" + tex.ped("{b,c}")))
ax.legend()


fig, ax = plt.subplots()

c = ColourCycles()
for i, df in enumerate(df_list):
    ax.plot(df["x"], df["Theta0"], c=c[i], alpha=.7, label=tex("k=%f/\mathrm{Mpc}" %k_list[i]))

ax.set_title(tex("\Theta_0"))
ax.legend()





plt.show()


