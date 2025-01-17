from utils import *
SET_SUBDIR("milestone3")

tex = LaTeX()


# this is dumb
def INIT(x_RM_equality, x_recombination, x_DE_domination):
    global x_eq, x_rec, x_Lambda#, x_eq_args0, x_eq_args1
    x_eq = x_RM_equality
    x_rec = x_recombination
    x_Lambda = x_DE_domination

    # global x_eq_args0, x_eq_args1, x_rec_args0, x_rec_args1, x_eq_rec_args0, x_eq_rec_args1

    # # for only eq:
    # x_eq_args0  = dict(x_list=[x_eq], x_label_list=[tex("x" + tex.ped("eq"))], style="-.", top=True)
    # x_eq_args1  = dict(x_list=[x_eq], style="-.")

    # # for only rec:
    # x_rec_args0 = dict(x_list=[x_rec], x_label_list=[tex("x" + tex.ped("*"))], style="-.", top=True)
    # x_rec_args1 = dict(x_list=[x_rec], style="-.")

    # # for both:
    # x_eq_rec_args0  = dict(x_list=[x_eq, x_rec], x_label_list=[tex("x" + tex.ped("eq")), tex("x" + tex.ped("*"))], style="-.", top=True)
    # x_eq_rec_args1  = dict(x_list=[x_eq, x_rec], style="-.")





# go-to kwargs:
alpha0  = .7
lw0     = 2.3
ls0     = "-"

xlims0  = (-15,0)
gs_kw0  = dict(height_ratios=(6,5))


def __set_k_label(fig, handles, k_dict, loc="center right"):

    k_labels = k_dict["val_str"]
    k_title = tex("k" + tex.unit("["+tex.inv("Mpc") +"]") + "=")

    fig.legend(handles=handles, labels=k_labels, title=k_title, 
               loc=loc, ncols=len(k_labels)+1, alignment="left", title_fontsize=16,
               handlelength=.8, fontsize=16, labelspacing=.3, columnspacing=.8, handletextpad=.4, 
               frameon=False)



def __vertical_lines(ax, k_dict, labels_on=True):

    pinpoint_args = dict(x_list=[(x_eq-.1, x_eq+.22), (x_rec-0.14, x_rec+0.26)], top=True, alpha=.12, style=None)

    if labels_on:
        pinpoint_args["x_label_list"] = [tex("\sim x" + tex.ped("eq")), tex("\sim x" + tex.ped("*"))]

    pinpoint_span_x(ax, **pinpoint_args)

    # line_kw = dict(alpha=.9, lw=1.9, dashes=(0.3, 1.2), gapcolor="w")
    # line_kw = dict(alpha=.9, marker="o", picker=0.1)       
    for i in range(len(k_dict["value"])):
        # ax.axvline(k_dict["x_enter"][i], color=k_dict["colour"][i], **line_kw)
        x_enter = k_dict["x_enter"][i]
        pad_x = 0.07
        hatch = "----"
        hatch = "/////"
        ax.axvspan(x_enter-pad_x, x_enter+pad_x, color=k_dict["colour"][i], ec=None, alpha=.5, hatch=hatch)






def __MatterPerturbations(which, df_list, k_dict):

    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(10,6), gridspec_kw=gs_kw0)
    ax1, ax2 = axes.flat

    k_handles = []

    c_kw = dict(alpha=alpha0, lw=lw0, ls=ls0)
    b_kw = dict(alpha=alpha0, lw=lw0, ls="--")
    g_kw = dict(alpha=alpha0, lw=.8, ls=":")
    gg_kw= dict(alpha=alpha0, lw=lw0, ls=ls0)

    which = which.strip().lower()
    assert which == "u" or which == "delta"

    for i, df in enumerate(df_list):
        line, = ax1.plot(df["x"], np.abs(df[f"{which}c"]), c=k_dict["colour"][i], **c_kw)
        ax1.plot(df["x"], np.abs(df[f"{which}b"]), c=k_dict["colour"][i], **b_kw)
        ax1.plot(df["x"], np.abs(df[f"{which}gamma"]), c=k_dict["colour"][i], **g_kw)
        
        ax2.plot(df["x"], df[f"{which}gamma"], c=k_dict["colour"][i], **gg_kw)
        k_handles.append(line)

    __vertical_lines(ax1, k_dict, False)
    __vertical_lines(ax2, k_dict, True)
        
    if which == "delta":
        which = "\delta"
    
    abs_label_dict = {s:tex("|" + which + tex.ped(s) + "|") for s in ["c", "b", "\gamma"]}
    ax1.plot(np.nan, np.nan, c="slategrey", label=abs_label_dict["c"], **c_kw)
    ax1.plot(np.nan, np.nan, c="slategrey", label=abs_label_dict["b"], **b_kw)
    ax1.plot(np.nan, np.nan, c="slategrey", label=abs_label_dict["\gamma"], **g_kw)
    ax2.plot(np.nan, np.nan, c="slategrey", label=tex(which + tex.ped("\gamma")), **gg_kw)


    u_legend_kw = dict(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper left', borderaxespad=0.)
    ax1.legend(ncols=3, **u_legend_kw)
    ax2.legend(**u_legend_kw)

    ax1.set_yscale("log")
    ax2.set_xlabel(tex.x)
    ax2.set_xlim(*xlims0)

    __set_k_label(fig, k_handles, k_dict, "center right")

    return fig, axes










def DensityPerturbations(df_list, k_dict, savefig=True):

    labels = dict(c     = tex("| \delta"+tex.ped("c") + "|"), 
                  b     = tex("| \delta"+tex.ped("b")+ "|"),
                  gamma = tex("\delta"+tex.ped("\gamma") + "=4" + tex.Theta + "_0"))
    
    fig, axes =__MatterPerturbations("delta", df_list, k_dict)
    ylims = (1.3e-2, axes[0].get_ylim()[1])
    axes[0].set_ylim(ylims)

    if savefig:
        save("density_perturbations")


def VelocityPerturbations(df_list, k_dict, savefig=True):

    labels = dict(c     = tex("| u"+tex.ped("c") + "|"), 
                  b     = tex("| u"+tex.ped("b")+ "|"),
                  gamma = tex("u"+tex.ped("\gamma") + "=-3" + tex.Theta + "_1"))

    fig, axes = __MatterPerturbations("u", df_list, k_dict)
    ylims = (0.6e-4, axes[0].get_ylim()[1])
    axes[0].set_ylim(ylims)

    if savefig:
        save("velocity_perturbations")



def PhotonQuadrupole(df_list, k_dict, savefig=True):

    fig, ax = plt.subplots()
    k_handles = []
    Theta2_kw = dict(alpha=alpha0, lw=lw0, ls=ls0)
    for i, df in enumerate(df_list):
        line, = ax.plot(df["x"], df["Theta2"], c=k_dict["colour"][i], **Theta2_kw)
        k_handles.append(line)

    ax.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Theta + "_2"), **Theta2_kw)
    __vertical_lines(ax, k_dict, True)
    legend_kw = dict(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper left', borderaxespad=0.)
    ax.legend(**legend_kw)

    __set_k_label(fig, k_handles, k_dict, "upper right")
    ax.set_xlabel(tex.x)
    ax.set_xlim(-11.2,0)
    


    if savefig:
        save("photon_quadrupole")



def GravitationalPotential(df_list, k_dict, savefig=True):

    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(10,6), gridspec_kw=gs_kw0)
    ax1, ax2 = axes.flat
    k_handles = []

    Phi_kw      = dict(alpha=alpha0, lw=lw0, ls=ls0)
    Phi_Psi_kw  = dict(alpha=alpha0, lw=lw0, ls=ls0)
    Psi_kw      = dict(alpha=.3,     lw=lw0, dashes=(6.3, 1.1))

    for i, df in enumerate(df_list):
        line, = ax1.plot(df["x"], df["Phi"], c=k_dict["colour"][i], **Phi_kw)
        ax1.plot(df["x"], -df["Psi"], c=k_dict["colour"][i], **Psi_kw)


        ax2.plot(df["x"], (df["Phi"] + df["Psi"]), c=k_dict["colour"][i], **Phi_Psi_kw)
        k_handles.append(line)
    

    __vertical_lines(ax1, k_dict, False)
    __vertical_lines(ax2, k_dict, True)

    ax1.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Phi), **Phi_kw)
    ax1.plot(np.nan, np.nan, c="slategrey", label=tex("-" + tex.Psi), **Psi_kw)
    ax2.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Psi + "+" + tex.Phi), **Phi_Psi_kw)
    legend_kw = dict(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper left', borderaxespad=0., ncols=1)
    ax1.legend(**legend_kw)
    ax2.legend(**legend_kw)
    ax2.set_xlim(*xlims0)

    __set_k_label(fig, k_handles, k_dict)
    ax2.set_xlabel(tex.x)
    # ax1.axvline(x_Lambda)

    # ylim = ax1.get_ylim()
    # ylims = (ylim[0], ylim[1]+0.08)
    # ax1.set_ylim(ylims)


    

    if savefig:
        save("gravitational_potential")


def GravitationalPotential_alternative(df_list, k_dict, savefig=False):

    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(10,6), gridspec_kw=gs_kw0)
    ax1, ax2 = axes.flat
    k_handles = []

    Phi_kw      = dict(alpha=alpha0, lw=lw0, ls=ls0)
    Phi_Psi_kw  = dict(alpha=alpha0, lw=lw0, ls=ls0)
    Psi_kw      = dict(alpha=.5, lw=lw0, ls="--")

    for i, df in enumerate(df_list):
        line, = ax1.plot(df["x"], df["Phi"], c=k_dict["colour"][i], **Phi_kw)
        kval = k_dict["value"][i]
        ax1.plot(df["x"], np.abs(df["Psi"]), c=k_dict["colour"][i], **Psi_kw)

        ax2.plot(df["x"], np.abs(df["Psi"]), c=k_dict["colour"][i], **Phi_Psi_kw)
        k_handles.append(line)
    

    __vertical_lines(ax1, k_dict, False)
    __vertical_lines(ax2, k_dict, True)

    ax1.plot(np.nan, np.nan, c="slategrey", label=tex(tex.Phi), **Phi_kw)
    ax2.plot(np.nan, np.nan, c="slategrey", label=tex("|" + tex.Psi + "|"), **Phi_Psi_kw)
    legend_kw = dict(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper left', borderaxespad=0.)
    ax1.legend(**legend_kw)
    ax2.legend(**legend_kw)
    ax2.set_xlim(*xlims0)

    __set_k_label(fig, k_handles, k_dict)
    ax2.set_xlabel(tex.x)



    if savefig:
        save("gravitational_potential_alt")




def SanityChecks(df_list):

    '''
    For by-eye comparing with Hans' results
    '''
    k_list = [0.001, 0.01, 0.1]
    k_colour = ColourCycles(["dodgerblue", "orange", "forestgreen"])
    k_label_func = lambda k: tex("k=%.3f" %k + tex.unit(tex.inv("Mpc")))
    k_label = [k_label_func(k) for k in k_list]

    kw0 = dict(alpha=alpha0, lw=lw0, ls=ls0)
    kw1 = dict(alpha=alpha0, lw=lw0, ls="--")

    #   DENSITY PERTURBATIONS

    fig1, ax1 = plt.subplots(figsize=(8,7))

    ax = ax1
    for i, df in enumerate(df_list):
        ax.plot(df["x"], np.abs(df["deltac"]), c=k_colour[i], label=k_label[i], **kw0)
        ax.plot(df["x"], np.abs(df["deltab"]), c=k_colour[i], **kw1)

    ax.set_title(tex("\delta"+ tex.ped("{c,b}")))
    ax.set_yscale("log")



    #   VELOCITY PERTURBATIONS

    fig2, ax2 = plt.subplots(figsize=(8,7))

    ax = ax2
    for i, df in enumerate(df_list):
        ax.plot(df["x"], np.abs(df["uc"]), c=k_colour[i], label=k_label[i], **kw0)
        ax.plot(df["x"], np.abs(df["ub"]), c=k_colour[i], **kw1)

    ax.set_title(tex("u"+ tex.ped("{c,b}")))
    ax.set_yscale("log")
    


    #   PHOTON MONOPOLE

    fig3, ax3 = plt.subplots(figsize=(8,7))

    ax = ax3
    for i, df in enumerate(df_list):
        ax.plot(df["x"], df["Theta0"], c=k_colour[i], label=k_label[i], **kw0)

    ax.set_title(tex(tex.Theta + "_0"))


    #   PHOTON DIPOLE

    fig4, ax4 = plt.subplots()

    ax = ax4
    for i, df in enumerate(df_list):
        ax.plot(df["x"], df["Theta1"], c=k_colour[i], label=k_label[i], **kw0)

    ax.set_title(tex(tex.Theta + "_1"))

    fig5, ax5 = plt.subplots(figsize=(8,7))


    #   GRAVITATIONAL POTENTIAL 

    ax = ax5
    for i, df in enumerate(df_list):
        ax.plot(df["x"], df["Phi"], c=k_colour[i], label=k_label[i], **kw0)

    ax.set_title(tex.Phi)


    #


    reproduce_lims = [(1e-1,1e5), (1e-6,1e2), (-0.5,1.0), (-0.4,0.4), (0.1,0.8) ]

    ax_list = [ax1, ax2, ax3, ax4, ax5]

    for i, ax in enumerate(ax_list):
        # ax.set_xlim(-18, -8.3)
        ax.legend()
        ax.set_xlabel(tex.x)
        # ax.axvline(-8.3, color="k", ls="-")s
        # ax.set_ylim(reproduce_lims[i])

    plt.show()
