from utils import *

const = ConstantsAndUnits()



'''
Pt 1) Test code
'''

data = read_ASCII("cosmology")


x = data[:,0]
eta = data[:,1]
t = data[:,2]
Hp = data[:,3]
dHpdx = data[:,4]
ddHpdxx = data[:,5]
Omegam = data[:,6]+data[:,7]
OmegaLambda = data[:,8]
Omegar = data[:,9]



jump = 50
x_RM = x[np.argmin(np.abs(Omegam-Omegar)[jump:])+jump]
x_MDE = x[np.argmin(np.abs(Omegam-OmegaLambda)[jump:])+jump]
x_acc = x[np.argmin(np.abs(dHpdx)[jump:])+jump]
x_0 = x[np.argmin(np.abs(x)[jump:])+jump]





def mark_milestones(ax):
    ax.axvline(x_RM, ls="--", color="orangered")
    ax.axvline(x_MDE, ls="--", color="olive")
    ax.axvline(x_acc, ls="--", color="firebrick")


def various_Hubble():
    fig, ax = plt.subplots()

    ax.plot(x, ddHpdxx/Hp, label=r"$\mathcal{H}''(x)/\mathcal{H}(x)$")
    ax.plot(x, dHpdx/Hp, label=r"$\mathcal{H}'(x)/\mathcal{H}(x)$")
    ax.plot(x, eta*Hp/const.c, label=r"$\eta(x)\mathcal{H}(x)/c$")
    mark_milestones(ax)

    # ax.set_xlim(np.min(x), 0)

    ax.legend()

    plt.show()


various_Hubble()



# fig, ax = plt.subplots()
# ax.plot(x, eta*Hp/const.c)
# fig, ax = plt.subplots()
# ax.plot(x, eta)

# ax.set_yscale("log")
# fig, ax = plt.subplots()
# ax.plot(x, Hp)
# ax.set_yscale("log")

fig, ax = plt.subplots()
ax.plot(x, Omegam, label=r"$\Omega_\mathrm{m}$")
ax.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
ax.plot(x, Omegar, label=r"$\Omega_\mathrm{r}$")
mark_milestones(ax)
ax.legend()

# plt.show()
# save("test")
plt.show()






'''
Pt 2) Fit to supernova data
'''


data_sn = read_ASCII("mcmc_fitting", 200)

Chi2 = data_sn[:,0]
h = data_sn[:,1]
Omega_m = data_sn[:,2]
Omega_k = data_sn[:,3]

Chi2min = np.min(Chi2)
Omega_Lambda = 1 - Omega_m - Omega_k

thresh = Chi2<Chi2+3.53

Omega_Lambda0, Omega_m0 = np.meshgrid(Omega_m[thresh], Omega_Lambda[thresh])
# fig, ax = plt.subplots()
# ax.scatter(Omega_m[thresh], Omega_Lambda[thresh], c=Chi2[thresh])
# ax.set_aspect("equal")
# # plt.colorbar()
# # ax.contourf()

# plt.show()


# fig, ax = plt.subplots()

# ax.hist(h, 30)
# plt.show()
