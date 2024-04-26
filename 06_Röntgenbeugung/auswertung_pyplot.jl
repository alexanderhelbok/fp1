using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings
using PhysicalConstants.CODATA2018: c_0, h, e

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")  
@py import matplotlib as mpl

# define Constants
begin
d = 201.4e-12
U = 35 * u"kV"
end

theta_to_d(theta) = @. 2*d*sin(theta*π/180)

# load data
begin
df = CSV.read(string(@__DIR__, "/data/Ex1.csv"), DataFrame, header=["V", "N1", "N2", "N3"], skipto=2)
df.N = measurement.(mean([df.N1, df.N2, df.N3]), sqrt.(mean([df.N1, df.N2, df.N3])))
end

# calculate weighted mean of plateau region
start = 9

mean_N = mean(df.N[start:end])

model(x, p) = @. p[1]*exp(-(x - p[3])/p[2]) + p[4]


# popt, ci = bootstrap(model, df.V[start:end-8], nom.(df.N[start:end-8]), yerr=err.(df.N[start:end-8]), p0=[.1, .4, 330., nom.(mean_N)], redraw=false, unc=true, xlim=[300, nom.(df.V[end]) + 50])
popt, ci = bootstrap(model, df.V[start:end], nom.(df.N[start:end]), p0=[.1, 2.5, 340., nom.(mean_N)], redraw=false, unc=true, xlim=[300, nom.(df.V[end]) + 50])
nom.(popt[3]) + popt[2]
popt[4]

# plot data
begin
ax = subplots(figsize=(7, 4.5))[1]
ax.errorbar(df.V, nom.(df.N), err.(df.N), fmt=".", capsize=3, mfc="white", mec="k", ms=11, ecolor="k", label=L"\mathrm{Messwerte}")  
ax.errorbar(df.V[start:end], nom.(df.N)[start:end], err.(df.N)[start:end], fmt=".", capsize=3, mfc="crimson", mec="k", ms=11, ecolor="k")  

xlims, ylims = ax.get_xlim(), ax.get_ylim()
ax.plot(ci.x, model(ci.x, nom.(popt)), label=L"\mathrm{Fit}", c="C1", zorder=0)

text(0.6, 0.45, L"$N(U) = A\mathrm{e}^{-\frac{U - U_0}{\tau}} + N_0$", transform=ax.transAxes, fontsize=14)
text(0.6, 0.34, L"$U_0 + \tau = 340.3(4)\ (\mathrm{V})$", transform=ax.transAxes, fontsize=14)
text(0.6, 0.24, L"$N_0 = 1012(9) \ (\mathrm{1/s})$", transform=ax.transAxes, fontsize=14)
rect = mpl.patches.FancyBboxPatch((0.6, 0.22), 0.33, 0.3, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

axins = ax.inset_axes([0.225, 0.15, 0.3, 0.35])

axins.errorbar(df.V, nom.(df.N), err.(df.N), fmt=".", capsize=3, mfc="white", mec="k", ms=11, ecolor="k")  
axins.errorbar(df.V[start:end], nom.(df.N)[start:end], err.(df.N)[start:end], fmt=".", capsize=3, mfc="crimson", mec="k", ms=11, ecolor="k")  
axins.plot(ci.x, model(ci.x, nom.(popt)), label="Mittelwert", c="C1", zorder=0)

# change fontsize for axins
for item in (axins.get_xticklabels() + axins.get_yticklabels())
    item.set_fontsize(12)
end

axins.set_xlim(312, 338)
axins.set_ylim(900, 1490)

# change visibility of connector lines
rectpatch, connects = ax.indicate_inset_zoom(axins, edgecolor="gray", alpha=1, zorder=0)
connects[3].set_visible(true)
connects[2].set_visible(false)
connects[1].set_visible(false)
connects[0].set_visible(true)


ax.set_xlim(xlims)
ax.set_ylim(ylims)

xlabel(L"U\ (\mathrm{V})")
ylabel(L"N\ (\mathrm{1/s})")

legend(loc="upper right")
tight_layout()
savefig(string(@__DIR__, "/bilder/Ex1.pdf"), bbox_inches="tight")
end


# =============== Ex3 ===============
# load data
begin
ax = subplots(1, 2, figsize=(8, 4.5))[1]
Voltages = [15:2.5:30;]
λ0 = []
for (i, V) in enumerate(Voltages)
    println(i)
    df = CSV.read(string(@__DIR__, "/data/Ex3_$i.csv"), DataFrame, header=["ang", "N"], skipto=3)
    df.N = measurement.(df.N, sqrt.(df.N))
    df.lam = theta_to_d(df.ang)*1e12

    mean_init = mean(df.N[1:1])
    idx_init = findfirst(df.N .> mean_init + 2*nom(mean(df.N[1:4]))) - 1
    # calculate mean of first datapoints
    mean_N = mean(df.N[1:idx_init])
    # set error to std
    mean_N = measurement(nom(mean_N), nom.(std(df.N[1:idx_init])))

    # find first datapoint where N is above 3 sigma of the mean
    idx = findfirst(df.N .> mean_N + 3*err(mean_N))
    # # interpolate to find the exact value
    λ0 = [λ0; measurement((df.lam[idx-1] + df.lam[idx])/2, (df.lam[idx] - df.lam[idx-1])/2)]

    if i == 5
        ax[0].errorbar(df.lam, nom.(df.N), yerr=err.(df.N), label=L"\mathrm{Messwerte}", fmt=".", capsize=3, mfc="white", mec="k", ms=11, ecolor="k")  
        ax[0].errorbar(df.lam[1:idx_init], nom.(df.N[1:idx_init]), yerr=err.(df.N[1:idx_init]), fmt=".", capsize=3, mfc="crimson", mec="k", ms=11, ecolor="k")
        ax[0].plot(df.lam, nom.(mean_N) * ones(length(df.lam)), label=L"\mathrm{Mittelwert}", color="C1")
        ax[0].fill_between(df.lam, nom.(mean_N) - 3*err.(mean_N), nom.(mean_N) + 3*err.(mean_N), color="C1", alpha=0.3, label=L"$3\sigma\ \mathrm{Konfidenzband}$")

        ax[0].set_xlabel(L"\lambda\ (\mathrm{pm})")
        ax[0].set_ylabel(L"N\ (\mathrm{1/s})")

        ax[0].legend()
    end
end

λ0
model(x, p) = @. p[1]/x

popt, ci = bootstrap(model, Voltages, nom.(λ0), yerr=err.(λ0), p0=[1000.], redraw=false, unc=true, p=0.66)

χ = chisq(λ0, model(Voltages, nom.(popt)), sigma=err.(λ0), pcount=1)

h_exp = popt[1]*u"kV*pm"*e/c_0 |> u"J*s"

ax[1].errorbar(Voltages, nom.(λ0), yerr=err.(λ0), label=L"\mathrm{Messwerte}", fmt=".", capsize=3, mfc="crimson", mec="k", ms=11, ecolor="k")
xlims, ylims = ax[1].get_xlim(), ax[1].get_ylim()
ax[1].plot(ci.x, model(ci.x, nom.(popt)), label=L"\mathrm{Fit}", color="C1")
# fill_between(ci.x, ci.c0, ci.c1, color="C1", alpha=0.3, label="a")

ax[1].text(0.2, 0.88, L"$\lambda_{\mathrm{min}}(U) = \frac{h_{\mathrm{exp}}c}{eU}$", transform=ax[1].transAxes, fontsize=14)
ax[1].text(0.2, 0.78, L"$h_{\mathrm{exp}} = 6.626(15) \cdot 10^{-34}\ (\mathrm{Js})$", transform=ax[1].transAxes, fontsize=14)
# ax[1].text(0.2, 0.74, L"$\chi_{\nu} \approx 0.2$", transform=ax[1].transAxes, fontsize=14)
rect = mpl.patches.FancyBboxPatch((0.2, 0.76), 0.74, 0.19, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax[1].transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))

ax[1].add_patch(rect)

ax[1].set_xlim(xlims)
ax[1].set_ylim(ylims)

ax[1].set_xlabel(L"U\ (\mathrm{kV})")
ax[1].set_ylabel(L"\lambda_{\mathrm{min}}\ (\mathrm{pm})")

legend(loc="lower left")
tight_layout()
# savefig(string(@__DIR__, "/bilder/Ex3.pdf"), bbox_inches="tight")
end

h_exp = λ*u"pm"*e*U / (c_0)



# load new data
begin
df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
dfAl = CSV.read(string(@__DIR__, "/data/Ex4_1Al.csv"), DataFrame, header=["V", "ang", "N1", "N2", "N3"], skipto=2)
dfZn = CSV.read(string(@__DIR__, "/data/Ex4_1Zn.csv"), DataFrame, header=["V", "ang", "N1", "N2", "N3"], skipto=2)
for df in [dfAl, dfZn]
    df.N = measurement.(mean([df.N1, df.N2, df.N3]), sqrt.(mean([df.N1, df.N2, df.N3])))
    df.lam = theta_to_d(df.ang)
end
end

df
theta_to_d(unique(dfAl.ang))


df.N[findfirst(df.ang .== unique(dfAl.ang)[1])]
reference = [df.N[findfirst(df.ang .== unique(dfAl.ang)[1])], df.N[findfirst(df.ang .== unique(dfAl.ang)[2])]]
reference = [1, 1]

model(x, p) = @. p[1]*exp(-p[2]*(x - p[3]))

labels = ["Aluminium 104pm", "Aluminium 73pm", "Zink 104pm", "Zink 73pm"]
ρAl = 2.7

# plot data
begin
ax = subplots()[1]

for (j, df) in enumerate([dfAl, dfZn])
    for (i, lam) in enumerate(unique(df.lam))
        mask = df.lam .== lam
        # add reference to data
        # x, y = [0.; df.V[mask]], [measurement(1., 0.); df.N[mask]]
        x, y = df.V[mask], df.N[mask]
        myerrorbar(df.V[mask], df.N[mask] ./ reference[i], fmt=".", capsize=3, mfc="C$(2*j+i-3)", mec="k", ms=11, ecolor="k", label=labels[2*j+i-2])

        # if df == dfZn
            popt, ci = bootstrap(model, x, nom.(y), yerr=err.(y), p0=[50, 20., 0.25, .5], unc=true, redraw=false, its=1000)

            if j == 1
                println(popt[2]*10/densitydict["Al"])
            else
                println(popt[2]*10/densitydict["Zn"])
            end
            # popt = [0.5, 20, 0.1, .5]
            newx = [0:0.001:0.1;]
            plot(newx, model(newx, nom.(popt)) ./ reference[i], color="C$(2*j+i-3)", lw=2, alpha=0.7)
        # end
        # scatter(df.V[mask], df.N[mask], color="black", s=30, zorder=10)
    end
end
xlims, ylims = ax.get_xlim(), ax.get_ylim()

# legend

# myerrorbar(0, 1, label="Aluminium", fmt=".", capsize=3, mfc="white", mec="k", ms=11, ecolor="k")



xlabel(L"d\ (\mathrm{mm})")
ylabel(L"N\ (\mathrm{1/s})")

legend()
tight_layout()
# savefig(string(@__DIR__, "/bilder/Ex4_1.pdf"), bbox_inches="tight")
end

thicknessdict = Dict("Al" => 0.04, "Zn" => 0.05, "Cu" => 0.025, "Sn" => 0.025, "Ni" => 0.025)
densitydict = Dict("Al" => 2.7, "Zn" => 7.14, "Cu" => 8.96, "Sn" => 7.3, "Ni" => 8.9)
Zdict = Dict("Al" => 13, "Zn" => 30, "Cu" => 29, "Sn" => 50, "Ni" => 28)

# load data
begin
I0df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
I0df.N = measurement.(I0df.N, sqrt.(I0df.N))
I0df.lam = theta_to_d(I0df.ang)
df = CSV.read(string(@__DIR__, "/data/Ex4_Cu.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang)
I0df = filter(row -> row.ang in unique(df.ang), I0df)
end

μ = @. -log(df.N / I0df.N)/thicknessdict["Zn"]
h*c_0./(df.lam[2:end] .* u"m") .|> u"eV"

@. μ[2:end]*u"1/mm" / densitydict["Cu"] / u"g/cm^3" / u"cm^2/g" |> NoUnits

# plot data
begin
ax = subplots(figsize=(7.5, 3.5))[1]
# myerrorbar(df.lam[2:end], (μ./maximum(nom.(μ)))[2:end].^(1/3), fmt="o", label="Messwerte")
myerrorbar(h*c_0./(df.lam[2:end] .* u"m") ./ u"MeV" .|> NoUnits , μ[2:end].*u"1/mm" ./ densitydict["Cu"] ./ u"g/cm^3" ./ u"cm^2/g" .|> NoUnits, fmt="o", label="Messwerte")
# scatter(df.lam, nom.(df.N)./maximum(nom.(df.N)), color="black", s=30, zorder=10)
# scatter(I0df.lam, nom.(I0df.N)./maximum(nom.(I0df.N)), color="red", s=30, zorder=10)
# scatter(df.lam, nom.(df.N), color="black", s=30, zorder=10)
# scatter(I0df.lam, nom.(I0df.N), color="red", s=30, zorder=10)

# title("Zink 0.05mm")
# xlabel(L"\lambda\ (\mathrm{pm})")
xlabel(L"E \ (\mathrm{MeV})")
ylabel(L"\mu/\rho\ (\mathrm{cm^2/g})")

tight_layout()
# xscale("log")

# savefig(string(@__DIR__, "/bilder/Ex4_Zn.pdf"), bbox_inches="tight")
end