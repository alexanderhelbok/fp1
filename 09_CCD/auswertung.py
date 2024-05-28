from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

# define plot params
# plt.style.use("Source.mplstyle")
# plt.rc('text', usetex=True)  # enable use of LaTeX in matplotlib
# plt.rc('font', family="serif", size=14)  # font settings

# mpl.rcParams["axes.labelpad"] = 10

h, w = 1024, 1536

def load_data(file):
    hdul = fits.open(file)
    data = hdul[0].data
    return data

def plot_data(ax, data, cmap="seismic", **cbar_kwargs):
    ax.imshow(data, cmap=cmap, origin="lower")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(ax.images[0], cax=cax, **cbar_kwargs)
    # plt.show()


# calculate mean of multiple images
def calc_mean(file, stop, start=3):
    arr = np.zeros((h, w, stop-start))
    for i in range(start, stop):
        arr[:, :, i-start] = load_data(f"/home/taco/Documents/fp1/09_CCD/data/{file}_{i}.fits")
    mean = np.mean(arr, axis=2)
    return mean

# plot histogram of data
def plot_hist(ax, data, **kwargs):
    binarr = data.flatten()
    ax.hist(binarr, bins=range(int(min(binarr)), int(max(binarr)), 1), color='white', edgecolor='black', histtype='stepfilled', **kwargs)
    # plt.show()

def plot_rect(ax, x, y, w, h):
    ax.plot([x, x], [y, y+h], color='red')
    ax.plot([x, x+w], [y, y], color='red')
    ax.plot([x, x+w], [y+h, y+h], color='red')
    ax.plot([x+w, x+w], [y, y+h], color='red')

# cut out a windwo of pixels
def cut_window(data, x, y, w, h):
    return data[y:y+h, x:x+w]


def parabola(x, a, b, c):
    return a*x**2 + x/b + c


def parabola2(x, a, b):
    return a*x**2 + x/b + std_ron


def line(x, a, b):
    return a*x + b
# %%
# load 20 bias images and calculate the mean for each pixel
bias = calc_mean("/EX1/Ex1_bias", 23)
# exclude 50 first and 50 last pixels in x and y direction
bias = bias[50:-50, 50:-50]
fig, ax = plt.subplots()
plot_data(ax, bias)

fit, ax = plt.subplots()
# plot histogram of bias values, with binsize 1
plot_hist(ax, bias)
# %%
# =============== Ex1 ===============
# calculate mean and std_dev of bias
mean_bias = np.mean(bias)
std_dev_bias = np.std(bias)
print(f"Mean of bias: {mean_bias}")
print(f"Standard deviation of bias: {std_dev_bias}")

# subtract the bias from one image
data = load_data("/home/taco/Documents/fp1/09_CCD/data/EX1/Ex1_bias_13.fits")
data = data[50:-50, 50:-50] - bias
fig, ax = plt.subplots()
plot_data(ax, data)
fit, ax = plt.subplots()
plot_hist(ax, data)

std_ron = np.std(data)
mean_ron = np.mean(data)
print(f"Mean of data: {mean_ron}")
print(f"Standard deviation of data: {std_ron}")


# %%
# ================== Ex2 ==================
g_all = np.zeros(2)
for (i, color) in enumerate(["Green", "Blue"]):
    # load bias frames
    bias = calc_mean(f"/EX2/{color}/ex2_bias", 8)
    # load flat frame
    flat = load_data(f"/home/taco/Documents/fp1/09_CCD/data/EX2/{color}/ex2_flat_0.fits")
    # correct flat frame
    flat = flat - bias
    # plot flat frame
    # plot_data(flat)

    # calculate g factor from mean and std_dev of flat frame cutting out a window
    cutflat = cut_window(flat, 600, 150, 200, 200)
    mean_flat = np.mean(cutflat)
    std_dev_flat = np.std(cutflat)
    g = mean_flat / std_dev_flat**2
    g_all[i] = g
    print(f"Mean of {color} flat: {mean_flat}")
    print(f"Standard deviation of {color} flat: {std_dev_flat}")
    print(f"g factor of {color} flat: {g}")

    fig = plt.figure()

    gs = fig.add_gridspec(2,2)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    plot_data(ax1, flat)
    # plot red rectangle
    plot_rect(ax1, 600, 150, 200, 200)

    plot_data(ax3, cutflat)
    plot_hist(ax2, cutflat, density=True)

    ax2.set_xlabel("Pixel Sättigung (ADU)")
    ax2.set_ylabel("Relative Häufigkeit")

    # change width of ax2
    if i != 0:
        ax2.set_aspect(20000)
    else:
        ax2.set_aspect(40000)
    plt.tight_layout()

    # offset ax3 tick labels
    ax3.set_xticklabels([0, 600, 650, 700, 750])
    ax3.set_yticklabels([0, 150, 200, 350, 400])

    ax3.spines['bottom'].set_color('r')
    ax3.spines['top'].set_color('r')
    ax3.spines['right'].set_color('r')
    ax3.spines['left'].set_color('r')
    # increase width of spines
    ax3.spines['bottom'].set_linewidth(2)
    ax3.spines['top'].set_linewidth(2)
    ax3.spines['right'].set_linewidth(2)
    ax3.spines['left'].set_linewidth(2)

    # plt.savefig(f"/home/taco/Documents/fp1/09_CCD/bilder/ex2_{color}.pdf", bbox_inches='tight')




# %%
# ================== Ex3 ==================
texp = [0.5, 1, 2, 3, 4, 5, 6, 7.5, 10, 15, 20, 22.5, 25, 26, 27, 28, 29, 30]
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
for (j, color) in enumerate(["Green", "Blue"]):
    means = unp.uarray(np.zeros(len(texp)), np.zeros(len(texp)))
    # load bias
    bias = calc_mean(f"/Ex3/{color}/ex3_bias", 8)
    for (i, t) in enumerate(texp):
        # load flat frame
        flat = load_data(f"/home/taco/Documents/fp1/09_CCD/data/Ex3/{color}/ex3_flat_{i}.fits")
        flat = flat - bias

        # calculate g factor from mean and std_dev of flat frame cutting out a window
        cutflat = cut_window(flat, 600, 150, 200, 200)
        mean_flat = np.mean(cutflat)
        std_dev_flat = np.std(cutflat)
        means[i] = unc.ufloat(mean_flat, std_dev_flat)
        # print(f"Mean of flat with texp {t}: {mean_flat}")
        # print(f"Standard deviation of flat with texp {t}: {std_dev_flat}")

    # norm data using mean of 10, 15, 20s exposure
    norm = np.mean(means[9:12]/texp[9:12])
    newmeans = means / norm.n

    ax[j].errorbar(texp, unp.nominal_values(newmeans/texp), yerr=unp.std_devs(newmeans/texp), fmt="o", capsize=3, ecolor="k", label=f"Datenpunkte")
    if j == 0:
        ax[j].set_title(f"Grüner Filter")
    else:
        ax[j].set_title(f"Blauer Filter")
    # dashed line for mean
    ax[j].axhline(1, color="gray", linestyle="--")
    ax[j].set_ylabel(r"$\frac{N^{\mathrm{ADU}}}{t_{\mathrm{exp}}N^{\mathrm{ADU}}_{\mathrm{mean}}}$")

plt.xlabel("Belichtungszeit (s)")
plt.tight_layout()
# plt.savefig("/home/taco/Documents/fp1/09_CCD/bilder/ex4.pdf", bbox_inches='tight')

# %%
# ================== Ex5 ==================
g2_all = unp.uarray(np.zeros(2), np.zeros(2))
# plot sigma^2 vs mean for green and blue filter
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
for (j, color) in enumerate(["Green", "Blue"]):
    means = unp.uarray(np.zeros(len(texp)), np.zeros(len(texp)))
    # load bias
    bias = calc_mean(f"/Ex3/{color}/ex3_bias", 8)
    for (i, t) in enumerate(texp):
        # load flat frame
        flat = load_data(f"/home/taco/Documents/fp1/09_CCD/data/Ex3/{color}/ex3_flat_{i}.fits")
        flat = flat - bias

        cutflat = cut_window(flat, 600, 150, 200, 200)
        mean_flat = np.mean(cutflat)
        std_dev_flat = np.std(cutflat)
        means[i] = unc.ufloat(mean_flat, std_dev_flat)
        # print(f"Mean of flat with texp {t}: {mean_flat}")
        # print(f"Standard deviati  on of flat with texp {t}: {std_dev_flat}")

    # fit parabola to data
    stop = 7

    popt, pcov = curve_fit(parabola, unp.nominal_values(means[:-stop]), unp.std_devs(means[:-stop]) ** 2, p0=[1e-6, 2, 0])
    perr = np.sqrt(np.diag(pcov))
    uncpopt = unp.uarray(popt, perr)
    # print(f"{color} filter: k = {uncpopt[0]:.1uS}, g = {uncpopt[1]:.1uS}, c = {uncpopt[2]:.1uS}")

    popt2, pcov2 = curve_fit(parabola2, unp.nominal_values(means[:-stop]), unp.std_devs(means[:-stop]) ** 2, p0=[1e-6, 0.5])
    perr = np.sqrt(np.diag(pcov2))
    uncpopt2 = unp.uarray(popt2, perr)
    print(f"{color} filter: k = {uncpopt2[0]:.2uS}, g = {uncpopt2[1]:.1uS}")

    g2_all[j] = uncpopt[1]

    ax[j].scatter(unp.nominal_values(means), unp.std_devs(means) ** 2, fc="white", ec="k", s=30, lw=1, label="Messpunkte")
    ax[j].scatter(unp.nominal_values(means[:-stop]), unp.std_devs(means[:-stop]) ** 2, fc="crimson", ec="k", s=30, lw=1, label="Messpunkte Fit")
    xlims, ylims = ax[j].get_xlim(), ax[j].get_ylim()
    ax[j].plot(unp.nominal_values(means), parabola(unp.nominal_values(means), *popt), label="Fit", c="C1")
    # ax[j].plot(unp.nominal_values(means), parabola2(unp.nominal_values(means), *popt2), label="Fit", c="C1")

    ax[j].text(0.05, 0.85, r"$f(x) = k_\lambda^2 x^2 + \frac{x}{g} + c$", transform=ax[j].transAxes, fontsize=14)
    if j == 0:
        ax[j].text(0.05, 0.72, r"$k_\lambda = 3.7(2) \cdot 10^{-6}$", transform=ax[j].transAxes, fontsize=14)
        ax[j].text(0.05, 0.59, r"$g = 2.11(4)$", transform=ax[j].transAxes, fontsize=14)
        ax[j].text(0.05, 0.46, r"$c = -30(90)$", transform=ax[j].transAxes, fontsize=14)
    else:
        ax[j].text(0.05, 0.72, r"$k_\lambda = 2.4(2) \cdot 10^{-6}$", transform=ax[j].transAxes, fontsize=14)
        ax[j].text(0.05, 0.59, r"$g = 2.15(5)$", transform=ax[j].transAxes, fontsize=14)
        ax[j].text(0.05, 0.46, r"$c = -30(90)$", transform=ax[j].transAxes, fontsize=14)
    rect = mpl.patches.FancyBboxPatch((0.06, 0.43), 0.26, 0.5, linewidth=1.5, edgecolor="C1", facecolor="none",
                                      transform=ax[j].transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))

    # ax[j].text(0.05, 0.85, r"$f(x) = k_\lambda^2 x^2 + \frac{x}{g} + 5.95$", transform=ax[j].transAxes, fontsize=14)
    # if j == 0:
    #     ax[j].text(0.05, 0.72, r"$k_\lambda = 3.72(12) \cdot 10^{-6}$", transform=ax[j].transAxes, fontsize=14)
    #     ax[j].text(0.05, 0.59, r"$g = 2.13(2)$", transform=ax[j].transAxes, fontsize=14)
    # else:
    #     ax[j].text(0.05, 0.72, r"$k_\lambda = 2.51(13) \cdot 10^{-6}$", transform=ax[j].transAxes, fontsize=14)
    #     ax[j].text(0.05, 0.59, r"$g = 2.16(2)$", transform=ax[j].transAxes, fontsize=14)
    # rect = mpl.patches.FancyBboxPatch((0.06, 0.53), 0.3, 0.4, linewidth=1.5, edgecolor="C1", facecolor="none",
    #                                   transform=ax[j].transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
    ax[j].add_patch(rect)

    ax[j].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax[j].ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    if j == 0:
        ax[j].set_title(f"Grüner Filter")
    else:
        ax[j].set_title(f"Blauer Filter")
    # dashed line for mean
    ax[j].set_ylabel(r"$\sigma_{\mathrm{ADU}}^2$")
    ax[j].legend(loc=(0.5, 0.07))
    ax[j].set_ylim(ylims)
    ax[j].set_xlim(xlims)
plt.xlabel(r"$N_{\mathrm{ADU}}$")
plt.tight_layout()
# plt.savefig("/home/taco/Documents/fp1/09_CCD/bilder/ex5_2.pdf", bbox_inches='tight')

# %%
# ================== Blooming ==================
# Set up figure and image grid
fig = plt.figure(figsize=(9.75, 3))

grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1, 2),
                 axes_pad=0.15,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=.15,
                 )
for (i, ax) in enumerate(grid):
    # load data
    data = load_data(f"/home/taco/Documents/fp1/09_CCD/data/blooming/blooming_{i}.fits")
    # cut out a window
    data = cut_window(data, 600, 0, 700, h)
    # plot data
    im = ax.imshow(data, cmap="viridis", origin="lower", vmin=550, vmax=2000)

thecb = ax.cax.colorbar(im)
thecb.set_label("Pixel Sättigung (ADU)", labelpad=10)

# plt.savefig("/home/taco/Documents/fp1/09_CCD/bilder/blooming.pdf", bbox_inches='tight')

# %%
# ================== Ex6 ==================
g_mean = np.mean(g2_all)
std_ron_e = g_mean * std_ron
print(f"Mean g factor: {g_mean:.1uS}")
print(f"Standard deviation of readout noise: {std_ron_e}")
# %%
from scipy.signal import find_peaks, savgol_filter
# ================== Ex7 ==================
Temps = np.array([-20, -12, -4, 4, 12, 20])
texp = [180, 180, 120, 120, 60, 60]
temp = unp.uarray(np.zeros(len(Temps)), np.zeros(len(Temps)))
dark_hist_all = np.zeros((len(Temps), 1700))
peakarr = []
fig, ax = plt.subplots(2, 3, figsize=(7, 4.5), sharex=True, sharey=True)
for (i, T) in enumerate(Temps[0:6]):
    # load bias
    bias = calc_mean(f"/Ex8/T{T}/bias", 8)
    # load dark frames
    dark = calc_mean(f"/Ex8/T{T}/dark", 3, start=0)
    dark = dark - bias

    # plot
    # if i == 2:
    #     fig, ax = plt.subplots(2, 1, figsize=(7, 6))
    #     plot_data(ax[0], dark)
    #     plot_hist(ax[1], dark)
    #     ax[0].set_title(f"Temperatur: {T}°C")
    #     # ax[1].set_yscale("log")

    # fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    if i < 3:
        myax = ax[0, i]
    else:
        myax = ax[1, 5-i]

    # find peaks in histogram
    histarr = np.histogram(dark.flatten(), bins=range(65535))[0][:1700]
    # smooth data
    if i < 3:
        histarr = savgol_filter(histarr, 5, 3)
        peaks, _ = find_peaks(histarr, height=5, threshold=1, prominence=5, distance=60)
    elif i < 5:
        histarr = savgol_filter(histarr, 5, 2)
        peaks, _ = find_peaks(histarr, height=10, threshold=1, prominence=1, distance=150)
    else:
        histarr = savgol_filter(histarr, 31, 3)
        peaks, _ = find_peaks(histarr, height=5, distance=350, prominence=3)
    print(peaks, histarr[peaks])
    peakarr.append(peaks)
    myax.scatter(peaks, histarr[peaks], color="red", s=10)
    # myax.plot(histarr, color="black")

    plot_hist(myax, dark)
    myax.set_title(f"T = {T}°C")
    myax.set_yscale("log")
    myax.set_xlim(-100, 1600)
    # myax.set_ylim(1, 1e5)

    dark_hist_all[i] = np.histogram(dark.flatten(), bins=range(65535), density=True)[0][:1700]

    # exclude pixels with values > 10000
    dark = dark[~(dark > 1000)]
    # calculate dark current
    mean_dark = unc.ufloat(np.mean(dark), np.std(dark))/(texp[i]*g_mean)
    # print(mean_dark)
    temp[i] = mean_dark

plt.tight_layout()
Temps = Temps + 273.15
# %%
# plt imshow of dark current histograms
fig, ax = plt.subplots()
ax.imshow(dark_hist_all.T, cmap="viridis", origin="lower", aspect="auto", norm=mpl.colors.LogNorm())
plt.colorbar(ax.images[0], ax=ax)
# %%
y = unp.log(temp/Temps**(3/2))
# y = temp/Temps**(3/2)
x = 1/Temps

# fit line to data
popt, pcov = curve_fit(line, x[2:], unp.nominal_values(y[2:]), p0=[-8000, 20])
perr = np.sqrt(np.diag(pcov))
uncpopt = unp.uarray(popt, perr)
print(popt)
E_g = -uncpopt[0]*(2*8.617333262145e-5) # eV
print(f"E_g = {E_g:.1uS}")
print(f"b = {uncpopt[1]:.1uS}")
popt2, pcov2 = curve_fit(line, x[:2], unp.nominal_values(y[:2]), p0=[-8000, 20])
perr = np.sqrt(np.diag(pcov2))
print(popt2)
E_g2 = -popt2[0]*(2*8.617333262145e-5) # eV
print(f"E_g2 = {E_g2}")
print(f"b2 = {popt2[1]}")

newx = np.linspace(0.0025, 0.005, 100)

fig, ax = plt.subplots(figsize=(7, 4))
ax.errorbar(1/Temps[2:], unp.nominal_values(y[2:]), yerr=unp.std_devs(y[2:]), fmt="o", capsize=3, ecolor="k", mfc="C1", mec="k")
ax.errorbar(1/Temps[:2], unp.nominal_values(y[:2]), yerr=unp.std_devs(y[:2]), fmt="o", capsize=3, ecolor="k", mfc="C0", mec="k")
xlims, ylims = ax.get_xlim(), ax.get_ylim()

ax.plot(newx, line(newx, *popt), c="C1")
ax.plot(newx[45:], line(newx[45:], *popt2), c="C0")

ax.text(0.69, 0.86, r"$f(T) = -\frac{E_g}{2k_BT} + b$", transform=ax.transAxes, fontsize=14)
ax.text(0.69, 0.76, r"$E_g = 1.1(2)\ eV$", transform=ax.transAxes, fontsize=14)
ax.text(0.69, 0.66, r"$b = 13(4)$", transform=ax.transAxes, fontsize=14)
rect = mpl.patches.FancyBboxPatch((0.69, 0.64), 0.26, 0.29, linewidth=1.5, edgecolor="C1", facecolor="none",
                                      transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

ax.text(0.05, 0.28, r"$f(T) = -\frac{E_g}{2k_BT} + b$", transform=ax.transAxes, fontsize=14)
ax.text(0.05, 0.18, r"$E_g = -0.34\ eV$", transform=ax.transAxes, fontsize=14)
ax.text(0.05, 0.08, r"$b = -19.37$", transform=ax.transAxes, fontsize=14)
rect = mpl.patches.FancyBboxPatch((0.05, 0.06), 0.26, 0.29, linewidth=1.5, edgecolor="C0", facecolor="none",
                                      transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

# add legend
ax.errorbar(0, 0, yerr=0, fmt="o", capsize=3, ecolor="k", mfc="gray", mec="k", label="Messpunkte")
ax.plot([-1, 0], [-1, 0], c="gray", label="Fit")

print(xlims)
ax.set_xlim(xlims)
ax.set_ylim(ylims[0]-0.2, ylims[1])
ax.set_xlabel(r"$T^{-1}\ (\mathrm{K}^{-1})$")
ax.set_ylabel(r"$\ln\left(\frac{I}{T^{3/2}}\right)$")

plt.legend(loc="upper center")
plt.tight_layout()
# plt.savefig("/home/taco/Documents/fp1/09_CCD/bilder/ex7.pdf", bbox_inches='tight')

# %%
# Temps = Temps + 273.15

# y = temp/Temps**(3/2)
x = 1/Temps
# print(peakarr)
# plot peaks
for (i, peaks) in enumerate(peakarr):
    print(peaks)
    for (j, peak) in enumerate(peaks):
        y = np.log(peak / Temps[i] ** (3 / 2))
        plt.scatter(x[i], y, c=f"C{j}")







