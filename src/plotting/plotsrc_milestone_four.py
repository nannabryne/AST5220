from utils import *

import matplotlib.ticker as ticker

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

# klims0 = (8e-4, 6e-1)
# k_axis_label = tex("k" + tex.unit("[" + tex("h") + tex.inv("Mpc")+"]") )
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



def TransferFunction(df, savefig=True):
    k = df["k"]
    df = df.drop("k", 1)

    
    
    def plot_functions():

        cc = ColourCycles()
        c = ColourCycles([cc[0], cc[4], cc[1], cc[5], cc[2], cc[3]])

        ell_handles = []
        ell_list = []

        Theta_kw = dict(alpha=.7, lw=1.8)
        for i, ell in enumerate(list(df.columns)):
            line, = ax1.plot(k, df[ell],  c=c[i],  **Theta_kw)
            ax2.plot(k, df[ell]*df[ell]/k, c=c[i],  **Theta_kw)
            ell_handles.append(line)
            ell_list.append(int(ell))

        ax1.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Theta + "_\ell"), **Theta_kw)
        ax1.legend(**legend_box_kw)

        ax2.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Theta + "_\ell" + "^2" + "/k"), **Theta_kw)
        ax2.legend(**legend_box_kw)

        for ax in [ax1, ax2]:
            ax.legend(**legend_box_kw)

        
        __set_ell_label(fig, ell_handles, ell_list, "upper right")


    #   create sexy figure:

    subpl_kw = dict(top=0.903,bottom=0.068,left=0.103,right=0.974,hspace=0.324)
    grids_kw = dict(height_ratios=(7,5))

    fig, axes = plt.subplots(figsize=(10,7.2), nrows=2, sharex=False, gridspec_kw=grids_kw)
    fig.set_tight_layout(False)
    fig.subplots_adjust(**subpl_kw)
    ax1, ax2 = axes.flat


    

    ax1.set_ylim(-0.0066, 0.0066)
    ax2.set_ylim(0, 2.1e-3)
    ax2.set_ylabel(tex(tex.unit(tex.inv("h")) + tex.unit("Mpc")))

    ax1.set_xlabel(k_axis_label)
    # ax2.set_xlabel(k_axis_label)

    kmin1 = 0.001
    kmin2 = 0.001
    kmax1 = 0.14
    kmax2 = 0.069
    ax1.set_xlim(kmin1, kmax1)
    ax2.set_xlim(kmin2, kmax2)

    #   plot connection line:

    from matplotlib.patches import ConnectionPatch

    con_kw = dict(color=plt.rcParams["axes.facecolor"], linestyle="--", alpha=.9, linewidth=1.1)

    con_l = ConnectionPatch(xyA=(kmin1,1), coordsA=ax2.get_xaxis_transform(), 
                            xyB=(kmin2,0), coordsB=ax1.get_xaxis_transform(),
                            **con_kw)

    con_r = ConnectionPatch(xyA=(1,1), coordsA=ax2.transAxes, 
                            xyB=(kmax2,0), coordsB=ax1.get_xaxis_transform(),
                            **con_kw)
    fig.add_artist(con_l)
    fig.add_artist(con_r)

    line_kw = dict(**con_kw)
    line_kw["color"] = "w"
    line_kw["alpha"] = line_kw["alpha"]*0.6

    for ax in [ax1, ax2]:
        ax.axvline(kmin2, **line_kw)
        ax.axvline(kmax2, **line_kw)


    # plot actual functions:

    plot_functions()
    

    
    if savefig:
        save("Theta_ell")
    




    

def CMBPowerSpectrum(df, df_obs, savefig=True):

    


    fig, ax = plt.subplots(figsize=(10,6))

    main_kws = dict(c=ColourCMB[-1], alpha=.9)
    part_kws = dict(lw=1.4, alpha=.69, ls="--")

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
    __set_D_ell_part_label(fig, handles, loc=(0.14, 0.76))


    #   observational data:

    err = np.zeros((2,len(df_obs["ell"])))
    err[0] = df_obs["err_down"]
    err[1] = df_obs["err_up"]
    ax.errorbar(df_obs["ell"], df_obs["D_ell"], err, label=tex("\mathcal{D}" + tex.ap("obs")), **obs_err_kw)

    #   actual function (again...):

    # ax.plot(ell, df["D_ell"], **main_kws)

    
    # ----------


    ax.set_xlabel(tex(tex.ell))
    ax.set_ylabel(tex(tex.unit("\mu K" + "^2") ))
    
    ax.set_xscale("log")
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

    # minors = [5,50,500]
    # ax.get_xaxis().set_minor_locator(ticker.FixedLocator(minors))
    # ax.get_xaxis().set_minor_formatter(ticker.FixedFormatter([str(l) for l in minors]))

    def ang_scale(l):
        theta = np.array(l, float)
        near_zero = np.isclose(l, 0)
        theta[near_zero] = np.inf
        theta[~near_zero] = 180/l[~near_zero]
        return theta
    
    def inverse(theta):
        l = np.array(theta, float)
        near_zero = np.isclose(theta, 0)
        l[near_zero] = np.inf
        l[~near_zero] = np.pi/theta[~near_zero]
        return l

    secax = ax.secondary_xaxis("top", functions=(ang_scale, inverse))

    
    # secax.set_xlabel(tex(tex.theta + tex.unit("[^\circ]")))
    
    majors = [90,45,10,2,1,0.2]
    labels = [tex("%s ^\circ") %str(theta) for theta in majors]
    colour = "slategrey"
    width = .7
    tick_kw = dict(reset=True, direction="out", length=0, width=width, color=colour, labelcolor=colour, top=True, bottom=False, labeltop=True, labelbottom=False, grid_alpha=.4, grid_color=colour, grid_linewidth=width, grid_linestyle="-")
    # secax.set_xlabel(tex(tex.theta), fontdict={"color":colour})
    secax.xaxis.set_major_locator(ticker.FixedLocator(majors))
    secax.xaxis.set_minor_locator(ticker.NullLocator())
    secax.grid(True)
    secax.xaxis.set_tick_params("major", **tick_kw)
    secax.set_xticklabels(labels)


    
    

    

    ax.legend()

    if savefig:
        save("CMB_power_spectrum")




def MatterPowerSpectrum(df, df_obs, savefig=True):
    fig, ax = plt.subplots(figsize=(10,6))


    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

    # actual function:
    ax.plot(df["k"], df["P"], c=ColourMatter[-1], label=tex("P"))
    ylims = ax.get_ylim()
    # ax.set_ylim(ylims)

    comp_kw = dict(c="darkslategrey", alpha=.3)
    comp_kw = dict(c=ColourMatter[-2], alpha=.3)

    ax.plot(df["k"], df["P_c"], label=tex("P" + tex.ped("c")), **comp_kw)
    ax.plot(df["k"], df["P_b"], label=tex("P" + tex.ped("b")), ls="--", **comp_kw)

    # # primordial:
    # ax.plot(df["k"], df["P_R"], c="darkslategrey", label=tex("P_\mathcal{R}"), alpha=.3)
    

    # ks = 0.0144790*4
    # ks = np.pi*0.67/147 *4
    # ax.axvline(ks*1, ls="--", c="r")
    # ax.axvline(ks*2, ls="--", c="r")
    # ax.axvline(ks*3, ls="--", c="r")
    # ax.axvline(ks*4, ls="--", c="r")
    # ax.axvline(ks*5, ls="--", c="r")
    # ax.axvline(ks*6, ls="--", c="r")
    # ax.axvline(ks*7, ls="--", c="r")
    # ax.axvline(ks*8, ls="--", c="r")

    # observational data:

    err = np.zeros((2,len(df_obs["k"])))
    err[0] = df_obs["err_down"]
    err[1] = df_obs["err_up"]
    ax.errorbar(df_obs["k"], df_obs["P"], err, label=tex("P" + tex.ap("obs")), **obs_err_kw)

    # vertical lines:

    pinpoint_x(ax, [k_eq], [tex("k" + tex.ped("eq") )])

    # --------


    # ax.set_xlim(8e-5, 4e-1)

    ax.set_xlabel(tex("k" + tex.unit("[h") + tex.unit(tex.inv("Mpc"))+"]" ))
    ax.set_xlabel(k_axis_label)
    ax.set_ylabel(tex(tex.unit("h") + "^3" + tex.unit("Mpc") + "^{-3}"))
    

    ax.legend()

    if savefig:
        save("matter_power_spectrum")
    