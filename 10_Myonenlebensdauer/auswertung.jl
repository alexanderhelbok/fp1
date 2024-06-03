using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings, SpecialFunctions

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

# increase axes label pad using rcParams
mpl.rcParams["axes.labelpad"] = 10 

@py import numpy as np

# load data as array from txt file
df = CSV.read(string(@__DIR__, "/data/Poisson.txt"), DataFrame, header = ["N"])

# calculate statistics
Ntot = sum(df.N)
Stdtot = sqrt(Ntot)

μ, σ = mean(df.N), std(df.N)
stderr = σ/sqrt(2*(length(df.N) - 1))
Nstd = measurement(σ, stderr)
Nmean = measurement(μ, σ/sqrt(length(df.N)))
Stdpoisson = sqrt(Nmean)
μ*Stdtot/Ntot
err(Nmean)

# create histogram data
h = np.histogram(df.N, bins = [minimum(df.N):1:maximum(df.N);], density = false)
# convert to julia array
hy = pyconvert(Array, h[0])
hx = pyconvert(Array, h[1])[1:end-1]
# fit gauss to histogram
gauss(x, p) = @. p[3]/(p[2]*sqrt(2*pi)) * exp(-0.5 * ((x - p[1])/p[2])^2)
# poisson(x, p) = @. 400 * exp(-p[1]) * p[1]^x / gamma(x + 1)

popt, ci = bootstrap(gauss, hx, hy, p0=[μ, σ, 400], unc=true, its=10000, redraw=true, p=0.68)
# popt2, ci2 = bootstrap(poisson, hx, hy, p0=[μ], unc=true, its=10000, redraw=false)

# calculate chi2
gausschi = chisq(hy, gauss(hx, nom.(popt)), dof = length(hy) - 2)
# poissonchi = chisq(hy, poisson(hx, nom.(popt2)), dof = length(hy) - 1)

# create histogram
begin
ax = subplots(figsize = (7, 4.2))[1]
ax.hist(df.N,
        bins=[minimum(df.N):1:maximum(df.N);],
        density=false,
        label="Messwerte",
        edgecolor="black",
        alpha=0.5,
        align="left",
        zorder=0) 
ax.scatter(hx, gauss(hx, nom.(popt)), c="crimson", edgecolor="black", label="Gaussverteilung", s=40)
ax.fill_between(ci.x, ci.c0, ci.c1, color="crimson", alpha=0.3, label=L"1\sigma\ \mathrm{Konfidenzband}")

ax.plot(ci.x, gauss(ci.x, nom.(popt)), c="black", zorder=0, ls="--")   
# ax.plot(ci2.x, poisson(ci2.x, nom.(popt2)), c="black", zorder=0, ls="--")
# ax.scatter(hx, poisson(hx, nom.(popt2)), c="gold", edgecolor="black", label="Poissonverteilung", s=40)

text(0.05, 0.88, L"A = 4.0(2) \cdot 10^{2}", fontsize=14, transform=ax.transAxes)
text(0.05, 0.8, L"\mu = 23.3(3)", fontsize=14, transform=ax.transAxes)
text(0.05, 0.72, L"\sigma = 4.6(3)", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.06, 0.7), 0.2, 0.24, linewidth=1.5, edgecolor="crimson", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

ax.tick_params(axis="x", which="minor", bottom=false, top=false)
ax.tick_params(axis="x", which="major", top=false, direction="out", length=4)
# ax.set_xticks(np.arange(0, 20, 2, dtype=Int))

# text(0.04, 0.6, L"\mu = 24.0(3)", fontsize=14, transform=ax.transAxes)
# rect = mpl.patches.FancyBboxPatch((0.05, 0.58), 0.17, 0.08, linewidth=1.5, edgecolor="gold", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
# ax.add_patch(rect)

xlabel(L"N")
ylabel(L"\mathrm{Absolute\ Haufigkeit}")

legend(labels=[L"\mathrm{Messwerte}", L"\mathrm{Gaussverteilung}", L"1\sigma\ \mathrm{Konfidenzband}"])
tight_layout()
# savefig(string(@__DIR__, "/bilder/gauss1.pdf"), bbox_inches="tight")
end

# do the same with errors
# drop 0 in data
hx = hx[hy .> 0]
hy = hy[hy .> 0]
hyerr = sqrt.(hy)

popt, ci = bootstrap(gauss, hx, hy, p0=[μ, σ, 400], yerr=hyerr, unc=true, its=10000, redraw=false, p=0.68)
1000/popt[1]
chisq(hy, gauss(hx, nom.(popt)), sigma=hyerr, dof = length(hy) - 3)
begin
ax = subplots(figsize = (7, 4.2))[1]
ax.hist(df.N,
        bins=[minimum(df.N):1:maximum(df.N);],
        density=false,
        label="Messwerte",
        edgecolor="black",
        alpha=0.5,
        align="left",
        zorder=0)
c = ax.scatter(hx, gauss(hx, nom.(popt)), c="crimson", edgecolor="black", s=40, zorder=10)
d = ax.fill_between(ci.x, ci.c0, ci.c1, color="crimson", alpha=0.3, zorder=5)

# add uncertainties as dashed histogram 
b, = ax.bar(hx, 2*hyerr, bottom=hy-hyerr, color="none", edgecolor="black", hatch="xxxx", lw=1, zorder=4, alpha=0.7)
# ax.bar(hx, 2*hystd, bottom=hy-hystd, color="none", edgecolor="black", hatch="\\\\", lw=0, label="Statistische Unsicherheit", zorder=0)
ax.plot(ci.x, gauss(ci.x, nom.(popt)), c="black", zorder=0, ls="--")
# ax.plot(ci2.x, poisson(ci2.x, nom.(popt2)), c="black", zorder=0, ls="--")
# ax.scatter(hx, poisson(hx, nom.(popt2)), c="gold", edgecolor="black", label="Poissonverteilung", s=40)

text(0.05, 0.88, L"A = 4.0(3) \cdot 10^{2}", fontsize=14, transform=ax.transAxes)
text(0.05, 0.8, L"\mu = 23.3(3)", fontsize=14, transform=ax.transAxes)
text(0.05, 0.72, L"\sigma = 4.6(3)", fontsize=14, transform=ax.transAxes)
text(0.05, 0.64, L"\chi_\nu = 1.19", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.06, 0.63), 0.2, 0.31, linewidth=1.5, edgecolor="crimson", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

ax.tick_params(axis="x", which="minor", bottom=false, top=false)
ax.tick_params(axis="x", which="major", top=false, direction="out", length=4)

xlabel(L"N")
ylabel(L"\mathrm{Absolute\ Haufigkeit}")

handles, keys = ax.get_legend_handles_labels()

legend(handles = [handles[0], b, c, d], labels=[L"\mathrm{Messwerte}", L"\mathrm{Unsicherheit}", L"\mathrm{Gaussverteilung}", L"1\sigma\ \mathrm{Konfidenzband}"])
tight_layout()
# savefig(string(@__DIR__, "/bilder/gauss2.pdf"), bbox_inches="tight")
end


# load exp data
df = CSV.read(string(@__DIR__, "/data/ex2.txt"), DataFrame, header = ["N"])

# bin data
h = np.histogram(df.N, bins = [minimum(df.N):2:maximum(df.N);], density = false)
hx = pyconvert(Array, h[1])[1:end-1] .+ 0.5
hy = pyconvert(Array, h[0])
hx = hx[hy .> 0]
hy = hy[hy .> 0]
hyerr = sqrt.(hy)

htot = sum(hy)

# calculate values
tmin = minimum(df.N)
tmax = maximum(df.N)
tstd = std(df.N)
tmean = sum(hx .* hy) / htot ± tstd / sqrt(htot)
err(tmean)/nom(tmean)

# fit exp
model(x, p) = @. 2*p[2]/p[1] * exp(-x/p[1])
popt, ci = bootstrap(model, hx, hy, p0=[50., htot], yerr=hyerr, unc=true, its=10000, redraw=false, p=0.68)

1000/popt[1]

# cut ci where y = 1
ci = ci[model(ci.x, nom.(popt)) .> 1, :]
# plot data
begin
ax = subplots(figsize = (7, 4.2))[1]
ax.errorbar(hx, hy, yerr=hyerr, mfc="C0", fmt="o", capsize=3, label=L"\mathrm{Messwerte}", color="black", zorder=5)
xlims, ylims = ax.get_xlim(), ax.get_ylim()

ax.plot(ci.x, model(ci.x, nom.(popt)), color="crimson", label=L"\mathrm{Fit}", zorder=10, lw=2)
# ax.fill_between(ci.x, ci.c0, ci.c1, color="C1", alpha=0.5, label=L"1\sigma\ \mathrm{Konfidenzband}", zorder=8)

text(0.63, 0.86, L"N(\Delta t) = \frac{2N_0}{T_m}\exp\left( \frac{\Delta t}{T_m} \right)", fontsize=14, transform=ax.transAxes)
text(0.63, 0.72, L"N_{\text{ges}} = 2.44(6)\cdot 10^3", fontsize=14, transform=ax.transAxes)
text(0.63, 0.65, L"T_m = 42.3(1.6)\ \mathrm{ms}", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.64, 0.63), 0.32, 0.32, linewidth=1.5, edgecolor="crimson", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

yscale("log")
xlim(xlims)
ylim(ylims)

xlabel(L"\Delta t\ (\mathrm{ms})")
ylabel(L"\mathrm{Absolute\ Haufigkeit}")

legend(loc="lower left")
tight_layout()
# savefig(string(@__DIR__, "/bilder/exp1.pdf"))
end


# load exp data
df = CSV.read(string(@__DIR__, "/data/ex3.txt"), DataFrame, header = ["N"])
df.N *= 1e3

# bin data
h = np.histogram(df.N, bins = [minimum(df.N):0.1:maximum(df.N);], density = false)
hx = pyconvert(Array, h[1])[1:end-1] .+ 0.5
hy = pyconvert(Array, h[0])
hyerr = sqrt.(hy)

htot = sum(hy)

# fit exp
model(x, p) = @. p[3]/(10*p[1]) * exp(-x/p[1]) + p[2]
popt, ci = bootstrap(model, hx[2:end], hy[2:end], p0=[3., 18., 2000.], yerr=hyerr[2:end], unc=true, its=1000, redraw=false, p=0.68)

chisq(hy[2:end], model(hx[2:end], nom.(popt)), sigma=hyerr[2:end], dof = length(hy) - 3)

# interpolate exo to 0
model(0, popt)  
# integrate exp
-2*popt[3]*(exp(-1e3/popt[1]) - 1)

texp = 3.3 ± 0.1
(model(texp, popt) - popt[2])/model(texp, popt)

# cut ci where y = 1
ci = ci[model(ci.x, nom.(popt)) .> 1, :]
# plot data
begin
ax = subplots(figsize = (7, 4.2))[1]
ax.errorbar(hx, hy, yerr=hyerr, mfc="C0", fmt="o", capsize=3, label=L"\mathrm{Messwerte}", color="black", zorder=5)
ax.errorbar(hx[1], hy[1], yerr=hyerr[1], mfc="white", fmt="o", capsize=3, color="black", zorder=5)
yscale("log")
xlims, ylims = ax.get_xlim(), ax.get_ylim()

ax.plot(ci.x, model(ci.x, nom.(popt)), color="crimson", label=L"\mathrm{Fit}", zorder=10, lw=2)
# ax.fill_between(ci.x, ci.c0, ci.c1, color="C1", alpha=0.5, label=L"1\sigma\ \mathrm{Konfidenzband}", zorder=8)

text(0.43, 0.87, L"N(\Delta t) = \frac{N_{\mathrm{ges}}}{10\tau}\exp\left( \frac{\Delta t}{\tau} \right) + c ", fontsize=14, transform=ax.transAxes)
text(0.43, 0.75, L"\tau = 2.00(3)\  \mu\mathrm{s}", fontsize=14, transform=ax.transAxes)
text(0.43, 0.65, L"N_{\mathrm{ges}} = 2.41(2) \cdot 10^4", fontsize=14, transform=ax.transAxes)
text(0.43, 0.55, L"c = 17.9(5)", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.44, 0.53), 0.38, 0.42, linewidth=1.5, edgecolor="crimson", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlim(xlims)
ylim(ylims)

xlabel(L"\Delta t\ (\mu\mathrm{s})")
ylabel(L"\mathrm{Absolute\ Haufigkeit\ } N")

legend(loc="lower left")
tight_layout()
# savefig(string(@__DIR__, "/bilder/exp2.pdf"))   
end

x = [1:50;]
chis, taus, Ns, cs = zeros(x[end]), Array{Measurement}(undef, x[end]), Array{Measurement}(undef, x[end]), Array{Measurement}(undef, x[end])
for i in x
    # println(i)

    popt, ci = bootstrap(model, hx[i:end], hy[i:end], p0=[3., 18., 2000.], yerr=hyerr[i:end], unc=true, its=1000, redraw=false, p=0.68)

    chi = chisq(hy[i:end], model(hx[i:end], nom.(popt)), sigma=hyerr[i:end], dof = length(hy[i:end]) - 3)
    # println("χ² = ", chi)
    println("τ = ", popt[1], "N = ", popt[3], "c = ", popt[2])
    chis[i] = chi
    taus[i] = popt[1]
    Ns[i] = popt[3]
    cs[i] = popt[2]
end

[print("$(round(i, digits=2)),") for i in chis]
labels = [L"\tau", L"N", L"c", L"\chi_\nu^2"]
# plot data
begin
ax = subplots(4, 1, figsize = (7, 6), sharex=true)[1]
for (i, y) in enumerate([taus, Ns, cs, chis])
    ax[i-1].errorbar(x[2:end], nom.(y)[2:end], yerr=err.(y)[2:end], fmt="o", capsize=3, label=labels[i])
    # ax[i-1].legend()
end
xlabel(L"\textrm{Startwert\ für\ Fit}")
ax[3].set_ylabel(L"\chi_\nu^2")
ax[0].set_ylabel(L"\tau\ (\mu\mathrm{s})")
ax[1].set_ylabel(L"N_{\mathrm{ges}}")
ax[2].set_ylabel(L"c")

tight_layout()
# savefig(string(@__DIR__, "/bilder/test.pdf"), bbox_inches="tight")
end