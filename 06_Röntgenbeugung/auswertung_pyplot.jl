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

model(x, p) = @. p[1]*exp(-p[2]*(x - p[3])) + p[4]

popt, ci = bootstrap(model, df.V[start:end], nom.(df.N[start:end]), yerr=err.(df.N[start:end]), p0=[1., 1., 30., 10.], redraw=false, unc=true)

# plot data
begin
ax = subplots(figsize=(7, 4.5))[1]
ax.errorbar(df.V, nom.(df.N), err.(df.N), fmt="o", ecolor="k", label="Messwerte", ms=4)   
scatter(df.V[start:end], df.N[start:end], color="black", s=30, zorder=10)

ax.plot(df.V[start:end], nom.(mean_N) .* ones(length(df.V[start:end])), label="Mittelwert", c="C1")
ax.plot(ci.x, model(ci.x, nom.(popt)), label="Mittelwert", c="C1")

text(0.65, 0.8, L"$\overline{N} = 1017(10)\ (\mathrm{1/s})$", transform=ax.transAxes, fontsize=14)
rect = mpl.patches.FancyBboxPatch((0.65, 0.78), 0.28, 0.09, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

axins = ax.inset_axes([0.3, 0.1, 0.3, 0.3])

axins.errorbar(df.V, nom.(df.N), err.(df.N), fmt="o", ecolor="k", label="Messwerte", ms=4)  

# change fontsize for axins
for item in (axins.get_xticklabels() + axins.get_yticklabels())
    item.set_fontsize(12)
end

axins.set_xlim(312, 330)
axins.set_ylim(900, 1500)

ax.indicate_inset_zoom(axins, edgecolor="gray", alpha=1, zorder=0)

xlabel(L"U\ (\mathrm{V})")
ylabel(L"N\ (\mathrm{1/s})")

legend()
tight_layout()
# savefig(string(@__DIR__, "bilder/Ex1.pdf"), bbox_inches="tight")
end

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
        ax[0].errorbar(df.lam, nom.(df.N), yerr=err.(df.N), fmt="o", label="Messwerte")
        ax[0].scatter(df.lam[1:idx_init], nom.(df.N[1:idx_init]), color="black", s=30, zorder=10)
        ax[0].plot(df.lam, nom.(mean_N) * ones(length(df.lam)), label="Mittelwert", color="red")
        ax[0].fill_between(df.lam, nom.(mean_N) - 3*err.(mean_N), nom.(mean_N) + 3*err.(mean_N), color="red", alpha=0.3)

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

# plot 
# begin
# ax = subplots()[1]
ax[1].errorbar(Voltages, nom.(λ0), yerr=err.(λ0), fmt="o", ecolor="k", ms=4, label=L"\mathrm{Messwerte}")
xlim, ylim = ax[1].get_xlim(), ax[1].get_ylim()
ax[1].plot(ci.x, model(ci.x, nom.(popt)), label=L"\mathrm{Fit}", color="C1")
# fill_between(ci.x, ci.c0, ci.c1, color="C1", alpha=0.3, label="a")

ax[1].text(0.2, 0.9, L"$\lambda_{\mathrm{min}}(U) = \frac{h_{\mathrm{exp}}c}{eU}$", transform=ax[1].transAxes, fontsize=14)
ax[1].text(0.2, 0.82, L"$h_{\mathrm{exp}} = 6.626(15) \cdot 10^{-34}\ (\mathrm{Js})$", transform=ax[1].transAxes, fontsize=14)
ax[1].text(0.2, 0.74, L"$\chi_{\nu} \approx 0.2$", transform=ax[1].transAxes, fontsize=14)
# rect = mpl.patches.FancyBboxPatch((0.5, 0.62), 0.4, 0.27, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax[1].transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
# ax[1].add_patch(rect)

ax[1].set_xlim(xlim)
ax[1].set_ylim(ylim)

ax[1].set_xlabel(L"U\ (\mathrm{kV})")
ax[1].set_ylabel(L"\lambda_{\mathrm{min}}\ (\mathrm{pm})")

legend(loc="lower left")
tight_layout()
savefig(string(@__DIR__, "/bilder/Ex3.pdf"), bbox_inches="tight")
end

h_exp = λ*u"pm"*e*U / (c_0)



# load new data
begin
df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
dfAl = CSV.read(string(@__DIR__, "/data/Ex4_1Al.csv"), DataFrame, header=["V", "ang", "N1", "N2", "N3"], skipto=2)
dfZn = CSV.read(string(@__DIR__, "/data/Ex4_1Zn.csv"), DataFrame, header=["V", "ang", "N1", "N2", "N3"], skipto=2)
for df in [dfAl, dfZn]
    df.N1 = measurement.(df.N1, sqrt.(df.N1))
    df.N2 = measurement.(df.N2, sqrt.(df.N2))
    df.N3 = measurement.(df.N3, sqrt.(df.N3))
    df.N = mean([df.N1, df.N2, df.N3])
    df.lam = theta_to_d(df.ang)
end
end

df

df.N[findfirst(df.ang .== unique(dfAl.ang)[1])]
reference = [df.N[findfirst(df.ang .== unique(dfAl.ang)[1])], df.N[findfirst(df.ang .== unique(dfAl.ang)[2])]]

model(x, p) = @. p[1]*exp(-p[2]*(x - p[3])) + p[4]

# plot data
begin
ax = subplots()[1]

for df in [dfAl, dfZn]
    for (i, lam) in enumerate(unique(df.lam))
        mask = df.lam .== lam
        # add reference to data
        x, y = [0.; df.V[mask]], [measurement(1., 0.); df.N[mask]]
        # x, y = df.V[mask], df.N[mask]

        # popt, ci = bootstrap(model, x, nom.(y), yerr=err.(y), p0=[0.0005, 1., 0.5, 10.], unc=true, redraw=false)

        println(popt)
        myerrorbar(df.V[mask], df.N[mask] ./ reference[i], fmt="o")   
        # plot(ci.x, model(ci.x, nom.(popt)) ./ reference[i], color="red")
        # scatter(df.V[mask], df.N[mask], color="black", s=30, zorder=10)
    end
end
myerrorbar(0, 1, fmt="o", label="Aluminium")

xlabel(L"U (\mathrm{V})")
ylabel(L"N (\mathrm{1/s})")

legend()
tight_layout()

end

thicknessdict = Dict("Al" => 0.04, "Zn" => 0.05, "Cu" => 0.025, "Sn" => 0.025, "Ni" => 0.025)
densitydict = Dict("Al" => 2.7, "Zn" => 7.14, "Cu" => 8.96, "Sn" => 7.3, "Ni" => 8.9)
Zdict = Dict("Al" => 13, "Zn" => 30, "Cu" => 29, "Sn" => 50, "Ni" => 28)

# load data
begin
I0df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
I0df.N = measurement.(I0df.N, sqrt.(I0df.N))
I0df.lam = theta_to_d(I0df.ang)
df = CSV.read(string(@__DIR__, "/data/Ex4_Al2.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang)
I0df = filter(row -> row.ang in unique(df.ang), I0df)
end

μ = @. -log(df.N / I0df.N)/thicknessdict["Cu"]

# plot data
begin
ax = subplots()[1]
myerrorbar(df.lam[1:end-10], (μ./maximum(nom.(μ)))[1:end-10].^(1/3), fmt="o", label="Messwerte")
scatter(df.lam, nom.(df.N)./maximum(nom.(df.N)), color="black", s=30, zorder=10)
scatter(I0df.lam, nom.(I0df.N)./maximum(nom.(I0df.N)), color="red", s=30, zorder=10)
# scatter(df.lam, nom.(df.N), color="black", s=30, zorder=10)
# scatter(I0df.lam, nom.(I0df.N), color="red", s=30, zorder=10)
end