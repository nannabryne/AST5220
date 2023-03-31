from utils import *
SET_SUBDIR("milestone2")

tex = LaTeX()


#   for drawing nice line at recombination:
rec_kw = dict(x_label_list=[tex("x_*")], colour="k", width=1.3, top=True, alpha=.9)




def ElectronFraction(df, info, savefig=True):
    
    fig, ax = plt.subplots()


    ax.plot(df["x"], df["Xe"], color=COLOURS[0], label=tex("X_e"))
    ax.plot(df["x_Saha"], df["Xe_Saha"], **overplot_kw)
    # ax.set_title(tex.TITLE("Electron fraction"))
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel(tex.x)

    ylim = (1e-4,1.5)
    ax.set_ylim(ylim)
    ax.set_xlim(-8.4,-3.6)

    pinpoint_x(ax, [info["x_rec"]], **rec_kw)
    pinpoint_y(ax, [info["Xe_fo"]], [tex("X_{e}" + tex.ap("(fo)"))], right=True)

    if savefig:
        save("electron_fraction")


def BaryonTemperature(df, info, savefig=False):
    fig, ax = plt.subplots()

    ax.plot(df["x"], df["Tb"], color=COLOURS[2], label=tex("T" + tex.ped("b")))
    ax.set_title(tex.TITLE("Baryon temperature"))
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel(tex.x)
    ax.set_ylabel(tex(tex.unit("K")))
    pinpoint_x(ax, [info["x_rec"]], **rec_kw)
    ax.set_xlim(-8.4,-3.6)

    if savefig:
        save("baryon_temperature")


def OpticalDepth(df, info, savefig=True):
    
    fig, ax = plt.subplots()
    x = df["x"]

    ax.plot(x, df["tau"], label=tex.tau, **f_of_x_kw)
    ax.plot(x, -df["dtaudx"], label=tex("-"+tex.dv(tex.tau)), **dfdx_of_x_kw)
    ax.plot(x, df["ddtaudxx"], label=tex.dv(tex.tau, n=2), **ddfdxx_of_x_kw)
    ax.set_xlabel(tex.x)
    ax.set_yscale("log")
    ax.set_xlim(-9.3,-3.2)
    ax.set_ylim(1e-5, 1e4)
    ax.legend()

    pinpoint_x(ax, [info["x_rec"]],  **rec_kw)



    if savefig:
        save("optical_depth_misc")


def VisibilityFunction(df, info, savefig=True):
    fig, ax = plt.subplots()
    x = df["x"]

    ax.plot(x, df["gt"],          label=tex.gt, **f_of_x_kw)
    ax.plot(x, df["dgtdx"]/10,    label=tex(tex.frac(1,10) + tex.dv(tex.gt)), **dfdx_of_x_kw)
    ax.plot(x, df["ddgtdxx"]/300, label=tex(tex.frac(1,300)+tex.dv(tex.gt, n=2)), **ddfdxx_of_x_kw)


    pinpoint_x(ax, [info["x_rec"]], **rec_kw)

    ax.set_xlabel(tex.x)
    ax.set_xlim(-7.7, -5.7)
    ax.legend()


    if savefig:
        save("visibility_function_misc")