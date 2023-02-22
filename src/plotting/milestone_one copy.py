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
OmegaM = data[:,6]+data[:,7]    # total matter
OmegaLambda = data[:,8]
OmegaR = data[:,9]+data[:,10]   # total radiation


z = np.exp(-x)-1
Chi = eta[np.argmin(np.abs(x))] - eta

dL = Chi*np.exp(-x) /const.Mpc *1e-3

jump = 50
x_RM = x[np.argmin(np.abs(OmegaM-OmegaR)[jump:])+jump]
x_MDE = x[np.argmin(np.abs(OmegaM-OmegaLambda)[jump:])+jump]
x_acc = x[np.argmin(np.abs(dHpdx)[jump:])+jump]
x_0 = x[np.argmin(np.abs(x)[jump:])+jump]





def mark_milestones(ax):
    ax.axvline(x_RM, ls="--", lw=1, color="orangered")
    ax.axvline(x_MDE, ls="--", lw=1, color="olive")
    ax.axvline(x_acc, ls="--", lw=1, color="firebrick")


def various_Hubble():
    fig, ax = plt.subplots()

    ax.plot(x, ddHpdxx/Hp, label=r"$\mathcal{H}''(x)/\mathcal{H}(x)$")
    ax.plot(x, dHpdx/Hp, label=r"$\mathcal{H}'(x)/\mathcal{H}(x)$")
    ax.plot(x, eta*Hp/const.c, label=r"$\eta(x)\mathcal{H}(x)/c$")
    mark_milestones(ax)

    # ax.set_xlim(np.min(x), 0)

    ax.legend()

    plt.show()

def cosmological_parameters():
    fig, ax = plt.subplots()
    ax.plot(x, OmegaM, label=r"$\Omega_\mathrm{M}$")
    ax.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
    ax.plot(x, OmegaR, label=r"$\Omega_\mathrm{R}$")
    mark_milestones(ax)
    ax.legend()

    plt.show()



def luminosity_distance_vs_redshift():
    obs_data = read_ASCII("supernovadata", 1)
    z_obs = obs_data[:,0]
    dL_obs = obs_data[:,1]
    err_obs = obs_data[:,2]

    fig, ax = plt.subplots()

    
    ax.errorbar(z_obs, dL_obs, err_obs, c="orangered")
    ax.plot(z_obs, dL_obs, 'o', c="orangered")
    ax.set_xscale("log")
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()


    ax.plot(z, dL)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$d_L\, \mathrm[Gpc]$")

    # ax.set_xlim(np.min(z_obs), np.max(z_obs))
    # ax.set_ylim(np.min(dL_obs), np.max(dL_obs))
    

    plt.show()


def scaled_Hubble():
    fig, ax = plt.subplots()

    ax.plot(x, Hp*const.Mpc/(100*const.km/const.s), label=r"$\mathcal{H}(x)$")
    mark_milestones(ax)

    ax.set_ylabel(r"100 km/s / Mpc")
    ax.set_xlabel(r"$x$")
    # ax.set_xlim(np.min(x), 0)

    ax.legend()

    plt.show()


def time_measures():
    fig, ax = plt.subplots()

    ax.plot(x, t, label=r"$t(x)$")
    ax.plot(x, eta/const.c, label=r"$\eta(x)/c$")
    mark_milestones(ax)

    ax.set_ylabel("s")
    ax.set_xlabel(r"$x$")

    # ax.set_xlim(np.min(x), 0)

    ax.legend()

    plt.show()


def mcmc_result():

    data_sn = read_ASCII("mcmc_fitting", 100)

    Chi2 = data_sn[:,0]
    h = data_sn[:,1]
    Omega_m = data_sn[:,2]
    Omega_k = data_sn[:,3]
    H_0 = h*const.H0_over_h

    Chi2min = np.min(Chi2)
    Omega_r = 2*np.pi**2/30 * (const.k_b*2.755)**4 * 8*np.pi*const.G/(3*H_0*H_0)
    Omega_Lambda = 1 - Omega_m - Omega_k - Omega_r

    thresh = Chi2<Chi2min+3.53

    fig, ax = plt.subplots()
    ax.scatter(Omega_m[thresh], Omega_Lambda[thresh], c=Chi2[thresh], alpha=0.6)

    ax.plot((0, 1), (1, 0), ls='--', c='k')
    # ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0,1.4)

    ax.set_xlabel(r"$\Omega_{\mathrm{M}0}$")
    ax.set_ylabel(r"$\Omega_{\Lambda0}$")
    # plt.colorbar()
    # ax.contourf()

    plt.show()


    fig, ax = plt.subplots()
    # ax.hist(Omega_Lambda[thresh], 50)
    ax.hist(H_0[thresh]*const.Mpc/(100*const.km/const.s), 60)

    # ax.set_xlabel(r"$\Omega_{\Lambda0}$")
    ax.set_xlabel(r"$H_0$ [100 km/s / Mpc]")

    plt.show()


mcmc_result()

# cosmological_parameters()

# various_Hubble()

# luminosity_distance_vs_redshift()

scaled_Hubble()

time_measures()



# fig, ax = plt.subplots()
# ax.plot(x, eta*Hp/const.c)
# fig, ax = plt.subplots()
# ax.plot(x, eta)

# ax.set_yscale("log")
# fig, ax = plt.subplots()
# ax.plot(x, Hp)
# ax.set_yscale("log")


# # plt.show()
# save("test")
# plt.show()






'''
Pt 2) Fit to supernova data
'''


# data_sn = read_ASCII("mcmc_fitting", 200)

# Chi2 = data_sn[:,0]
# h = data_sn[:,1]
# Omega_m = data_sn[:,2]
# Omega_k = data_sn[:,3]

# Chi2min = np.min(Chi2)
# Omega_Lambda = 1 - Omega_m - Omega_k

# thresh = Chi2<Chi2+3.53

# Omega_Lambda0, Omega_m0 = np.meshgrid(Omega_m[thresh], Omega_Lambda[thresh])
# fig, ax = plt.subplots()
# ax.scatter(Omega_m[thresh], Omega_Lambda[thresh], c=Chi2[thresh])
# ax.set_aspect("equal")
# # plt.colorbar()
# # ax.contourf()

# plt.show()


# fig, ax = plt.subplots()

# ax.hist(h, 30)
# plt.show()
