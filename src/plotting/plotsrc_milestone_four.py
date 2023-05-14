from utils import *

from plotsrc_milestone_one import ColourCyclesI
__maincolours = ColourCyclesI("densityparams")
__subcolours = ColourCyclesI("densityparams_shy")
ColourCMB = ColourCycles()
ColourMatter = ColourCycles()
ColourCMB = ColourCMB + __subcolours[0] + __maincolours[0]
ColourMatter = ColourMatter + __subcolours[1] + __maincolours[1]


SET_SUBDIR("milestone4")

tex = LaTeX()


k_axis_label = tex("k" + tex.unit("[h") + tex.unit(tex.inv("Mpc"))+"]" )

klims0 = (8e-4, 6e-1)
k_axis_label = tex("k/h" + tex.unit("[" + tex.inv("Mpc")+"]") )
klims0 = (0.001, 0.08)

k_eq = 0.0115/0.67





def __set_ell_label(fig, handles, ell_list, loc="center right"):

    ell_labels = [f"{ell:.0f}" for ell in ell_list]
    ell_title = tex("\ell = ")

    fig.legend(handles=handles, labels=ell_labels, title=ell_title, 
               loc=loc, ncols=len(ell_labels)+1, alignment="left", title_fontsize=16,
               handlelength=.8, fontsize=16, labelspacing=.3, columnspacing=.8, handletextpad=.4, 
               frameon=False)



def __set_D_ell_part_label(fig, handles, ap="[comp]", loc="center right"):

    comp_labels = [comp for comp in ["SW", "ISW", "Doppler", "pol"]]
    comp_title = tex("\mathrm{%s}=" %ap)

    fig.legend(handles=handles, labels=comp_labels, title=comp_title, 
               loc=loc, ncols=len(comp_labels)+1, alignment="left", title_fontsize=16,
               handlelength=.8, fontsize=16, labelspacing=.3, columnspacing=.8, handletextpad=.4, 
               frameon=False)



def TransferFunction_old(df, savefig=True):
    k = df["k"]
    df = df.drop("k", 1)

    fig1, ax1 = plt.subplots()

    '''
    labels: as in pert-plots!!
    '''

    for ell in list(df.columns):
        ax1.plot(k, df[ell], alpha=.6, label=tex(tex.Theta + tex.ped(ell)))

    ax1.set_ylim(-0.012, 0.012)


    fig2, ax2 = plt.subplots()

    for ell in list(df.columns):
        ax2.plot(k, df[ell]*df[ell]/k, alpha=.6, label=tex(tex.Theta + tex.ped(ell) + "^2" +  "/k"))

    ax2.set_ylim(0, 2e-3)

    
    for ax in [ax1, ax2]:
        # ax.set_xscale("log")
        ax.set_xlabel(k_axis_label)
        ax.set_xlim(klims0)
        ax.legend()
    


def TransferFunction(df, savefig=True):
    k = df["k"]
    df = df.drop("k", 1)

    fig, axes = plt.subplots(figsize=(10,8), nrows=2, sharex=False)

    '''
    labels: as in pert-plots!!
    '''

    c = ColourCycles()
    ell_handles = []
    ell_list = []

    ax1, ax2 = axes.flat

    Theta_kw = dict(alpha=.7, lw=1.8)
    for i, ell in enumerate(list(df.columns)):
        line, = ax1.plot(k, df[ell],  c=c[i],  **Theta_kw)
        ax2.plot(k, df[ell]*df[ell]/k, c=c[i],  **Theta_kw)
        ell_handles.append(line)
        ell_list.append(int(ell))

    ax1.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Theta + "_\ell"), **Theta_kw)
    ax1.legend(**legend_box_kw)

    ax2.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Theta + "_\ell" + "^2" +  "/k"), **Theta_kw)
    ax2.legend(**legend_box_kw)


    ax1.set_ylim(-0.009, 0.009)
    ax2.set_ylim(0, 2e-3)

    ax2.set_xlabel(k_axis_label)
    # ax2.set_xlim(klims0)

    ax1.set_xlim(klims0)
    ax2.set_xlim(klims0)#[0], klims0[1]/2.)

    
    for ax in [ax1, ax2]:
        # ax.set_xscale("log")
        # ax.set_xlabel(k_axis_label)
        # ax.set_xlim(klims0)
        ax.legend(**legend_box_kw)

    __set_ell_label(fig, ell_handles, ell_list)
    
    if savefig:
        save("Theta_ell")
    




    

def CMBPowerSpectrum(df, df_obs, savefig=True):
    fig, ax = plt.subplots(figsize=(10,6))

    main_kws = dict(c=ColourCMB[-1])
    part_kws = dict(lw=1.4, alpha=.7, ls="--")

    ell = df["ell"]


    #   actual function:

    ax.plot(ell, df["D_ell"], label=tex("\mathcal{D}"), **main_kws)

    #   cosmic variance:

    var = (np.sqrt( 2./(2.*ell+1.)) * df["D_ell"])
    ax.fill_between(ell, df["D_ell"]-var, df["D_ell"]+var, color=ColourCMB[-2], alpha=.1, ec="darkslategrey")
    
    #   components:

    c = ColourCycles("miscellaneous", ax)
    part_list = ["SW", "ISW", "Doppler", "pol"]
    handles = []
    for i, part in enumerate(part_list):
        line, = ax.plot(ell, df[f"D_ell_{part}"], **part_kws)
        handles.append(line)


    ax.plot(np.nan, np.nan, c="slategrey", label=tex("\mathcal{D}" + tex.ap("[comp]")), **part_kws)
    __set_D_ell_part_label(fig, handles, loc=(0.14, 0.82))


    #   observational data:

    err = np.zeros((2,len(df_obs["ell"])))
    err[0] = df_obs["err_down"]
    err[1] = df_obs["err_up"]
    ax.errorbar(df_obs["ell"], df_obs["D_ell"], err, label=tex("\mathcal{D}" + tex.ap("obs")), **obs_err_kw)

    #   actual function (again...):

    ax.plot(ell, df["D_ell"], **main_kws)

    
    # ----------


    ax.set_xlabel(tex(tex.ell))
    ax.set_ylabel(tex(tex.unit("\mu K" + "^2") ))
    
    ax.set_xscale("log")

    ax.legend()

    if savefig:
        save("CMB_power_spectrum")




def MatterPowerSpectrum(df, df_obs, savefig=True):
    fig, ax = plt.subplots(figsize=(10,4))


    ax.set_xscale("log")
    ax.set_yscale("log")

    # actual function:
    
    ax.plot(df["k"], df["P"], c=ColourMatter[-1], label=tex("P"))

    # observational data:

    err = np.zeros((2,len(df_obs["k"])))
    err[0] = df_obs["err_down"]
    err[1] = df_obs["err_up"]
    ax.errorbar(df_obs["k"], df_obs["P"], err, label=tex("P" + tex.ap("obs")), **obs_err_kw)

    # vertical lines:

    pinpoint_x(ax, [k_eq], [tex("k" + tex.ped("eq") + "/h")])

    # --------


    # ax.set_xlim(8e-5, 4e-1)

    ax.set_xlabel(tex("k" + tex.unit("[h") + tex.unit(tex.inv("Mpc"))+"]" ))
    ax.set_xlabel(k_axis_label)
    ax.set_ylabel(tex(tex.unit("h") + "^3" + tex.unit("Mpc") + "^{-3}"))
    

    ax.legend()

    if savefig:
        save("matter_power_spectrum")
    