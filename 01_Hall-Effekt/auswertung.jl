using CSV, Roots, LaTeXStrings, PythonPlot
using PhysicalConstants.CODATA2018: e, k_B

include(string(@__DIR__, "/../Source.jl"))
# include(string(@__DIR__, "/../SourceStatistics.jl"))
@py import matplotlib as mpl

mpl.use("pgf")
mpl.use("TkAgg")

# define constants
begin
    d = measurement(500, 1) * u"μm"
    B = measurement(252.6, 0.01) * u"mT"
    B0 = measurement(1.045, 0.01) * u"mT"
    I = measurement(500, 500*0.015+0.5) * u"μA"
    dB = B - B0
end

function T(R) 
    temp = @. 509/8 - 89/(2*R) + 43/30 * sqrt(R) + 161/113 * R^(11/10)
    return temp * u"K"
end

begin
# load data
tempdf = CSV.read(joinpath(@__DIR__, "data.csv"), DataFrame, header=["T1", "1a", "1b", "2a", "2b", "T2", "3a", "3b", "3aB", "3bB"], skipto=2)
# make measurement objects
for Temp in ["T1", "T2"]
    tempdf[!, Temp] = measurement.(tempdf[!, Temp], 0.05)
    tempdf[!, Temp] = T.(tempdf[!, Temp])
end
for i in ["1a", "1b", "2a", "2b", "3a", "3b", "3aB", "3bB"]
    tempdf[!, i] = measurement.(tempdf[!, i], 0.05) * u"mV"
end

# interpolate temperature
temp = tempdf[!, "T1"]
tempdf.T1 = (tempdf.T1 + tempdf.T2)/2
tempdf.T2[1:end-1] = (temp[2:end] + tempdf.T2[1:end-1])/2

begin
df = DataFrame()
df.T1 = tempdf.T1
df.T2 = tempdf.T2
df.U1 = (tempdf."1a" - tempdf."1b")/2
df.U2 = (tempdf."2a" - tempdf."2b")/2
df.U3 = (tempdf."3a" - tempdf."3b")/2
df.U3B = (tempdf."3aB" - tempdf."3bB")/2
df.newT = @. (df.T1 + df.T2)/2
end

Roomdf = df[1:4, :]

# drop last row
df = df[1:end-1, :]
# drop first 4 rows
df = df[5:end, :]
end

# print(maximum(RH), minimum.(RH))
# plot RH against T
# begin
# ax = subplots(figsize=(7, 4.5))[1]
# myerrorbar(df.T1, RH, fmt=".", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
# # myerrorbar(df.T1, ρ, fmt=".", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")

# xlabel(L"$T\ (\textrm{K})$")
# # ylabel(L"$\rho\ (\Omega\ \textrm{m})$")
# ylabel(L"$R_H\ (\Omega\ \textrm{T}^{-1})$")
# end

# room temp values
begin
UH = Roomdf.U3B .- Roomdf.U3
RH = @. UH * d / (dB * I) .|> u"m^3/C"

function f(x)
    x = nom.(x)
    fun(unbekannt) = @. cosh(log(2) / unbekannt * (x - 1)/(x+1)) - 0.5 * exp(log(2)/unbekannt)
    return find_zero.(fun, 1)
end

x = @. abs(Roomdf.U1/Roomdf.U2)

F = [f(X) for X in x]

ρ = @. pi * d/log(2) * F * (abs(Roomdf.U1) + abs(Roomdf.U2)) / (2 * I) .|> u"Ω*m"

μ = -RH ./ ρ .|> u"cm^2/(V*s)"

n = @. -1/(RH * e) .|> u"m^-3"

println("RH = ", mean(RH))
println("μ = ", mean(μ))
println("ρ = ", mean(ρ))
println("n = ", mean(n))
println("T = ", mean(Roomdf.newT))
end

# ========== RH ==========
begin
UH = df.U3B .- df.U3
RH = @. UH * d / (dB * I)

x = @. abs(df.U1/df.U2)

F = [f(X) for X in x]

ρ = @. pi * d/log(2) * F * (abs(df.U1) + abs(df.U2)) / (2 * I)

μ = -RH ./ ρ
μ .|> u"cm^2/(V*s)"

n = @. -1/(RH * e)
logn = @. log(n*u"m^3")
# end
df

func(x, p) = @. p[1]*x + p[2]

interval = 35:length(df.newT)-8
popt, ci = bootstrap(func, 
                    nom.(1 ./df.newT[interval]), 
                    nom.(logn[interval]), 
                    xerr=err.(1 ./df.newT[interval]), 
                    yerr=err.(logn[interval]) , 
                    p0=[1.0, 1.0], 
                    its=10000, 
                    unc=true,
                    xlim=[0, 0.015])

print(popt)
Ed = - 2*k_B * popt[1] * u"K" |> u"eV"

# mean of first 6 values
const_model(x, p) = @. p[1] *x/x
mean1, ccc = bootstrap(const_model, nom.(1 ./df.newT[1:6]), nom.(n[1:6]), xerr=err.(1 ./df.newT[1:6]), yerr=err.(n[1:6]), its=10000, unc=true, p0=[1e28])
println("Ed = ", Ed)
println("ND - NA = ", mean1[1] * u"1/m^3")
println(log(mean1[1]))
# end

# n
# ========== Plot1 ==========
# begin
ax = subplots(figsize=(7, 4.5))[1]
myerrorbar(1 ./ df.newT, logn, fmt=".", capsize=3,  mfc="white", mec="k", ms=7, ecolor="k")
myerrorbar(1 ./ df.newT[interval], logn[interval], fmt=".", capsize=3, mfc="deepskyblue", mec="k", ms=7, ecolor="k")
xlims, ylims = ax.get_xlim(), ax.get_ylim()

# legend
errorbar(1, 1, xerr=1, yerr=1, fmt=".", capsize=3,  mfc="silver", mec="k", ms=11, ecolor="k", label=L"\mathrm{Messwerte}")
plot(1, 1, c="gray", label=L"\mathrm{fit}")
fill_between((1, 2), 1, 1, alpha=0.3, color="gray", label=L"$2 \sigma\ \mathrm{Konfidenzband}$")

plot(ci.x, func(ci.x, nom.(popt)), "C0", zorder=10)
fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C0")

axhline(log(nom(mean1[1])), c="C1")

text(0.32, 0.43, L"$f(x) = a \cdot x + b$", transform=ax.transAxes, fontsize=14)
text(0.32, 0.33, L"$a = -53(3)\ \mathrm{K}$", transform=ax.transAxes, fontsize=14)
text(0.32, 0.23, L"$b = 46.73(2)$", transform=ax.transAxes, fontsize=14)

rect = mpl.patches.FancyBboxPatch((0.32, 0.21), 0.24, 0.27, linewidth=1.5, edgecolor="C0", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

text(0.65, 0.8, L"$\overline{\ln(n)} = 46.486(7)$", transform=ax.transAxes, fontsize=14)
rect = mpl.patches.FancyBboxPatch((0.65, 0.78), 0.25, 0.09, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

axins = ax.inset_axes([0.23, 0.6, 0.3, 0.3])

axins.errorbar(nom.(1 ./ df.newT[1:6]), nom.(logn[1:6]), xerr=err.(1 ./ df.newT[1:6]), yerr=err.(logn[1:6]), fmt=".", capsize=3, mfc="crimson", mec="k", ms=11, ecolor="k")
axins.errorbar(nom.(1 ./ df.newT[7:15]), nom.(logn[7:15]), xerr=err.(1 ./ df.newT[7:15]), yerr=err.(logn[7:15]), fmt=".", capsize=3, mfc="white", mec="k", ms=11, ecolor="k")
axins.axhline(log(nom(mean1[1])), c="C1", zorder=4)
axins.fill_between(xlims, log(nom(mean1[1]) + err(mean1[1])), log(nom(mean1[1]) - err(mean1[1])), alpha=0.3, color="C1")
axins.fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C0")

# change fontsize for axins
for item in (axins.get_xticklabels() + axins.get_yticklabels())
    item.set_fontsize(12)
end

axins.set_xlim(0.003, 0.004)
axins.set_ylim(46.46, 46.51)

ax.indicate_inset_zoom(axins, edgecolor="gray", alpha=1, zorder=0)

xlim(xlims[0]+0.001, xlims[1] + 0.006)
ylim(ylims)

xlabel(L"$T^{-1}\ (\mathrm{K}^{-1})$")
ylabel(L"$\ln(n\cdot \mathrm{m}^3)$")
legend(loc=(0.59, 0.4))
tight_layout()
# savefig(string(@__DIR__, "/bilder/plot1.pdf"), bbox_inches="tight")
end

# ========== Plot2 =========
begin
ax = subplots(figsize=(7, 4.5))[1]
myerrorbar(1 ./ df.newT, μ .|> u"cm^2/(V*s)", fmt=".", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
xlims, ylims = ax.get_xlim(), ax.get_ylim()
myerrorbar(1, 1, fmt=".", capsize=3, mfc="white", mec="k", ms=11, ecolor="k", label=L"$\mathrm{Messwerte}$")

xlim(xlims)
ylim(ylims)

xlabel(L"$T^{-1}\ (\textrm{K}^{-1})$")
ylabel(L"$\mu\ \left( \textrm{cm}^2\textrm{V}^{-1}\textrm{s}^{-1} \right)$")
legend()
tight_layout()
# savefig(string(@__DIR__, "/bilder/plot2.pdf"), bbox_inches="tight")
end

μ .|> u"cm^2/(V*s)"

# ========== Plot3 =========
logμ = @. log(μ/u"cm^2/(V*s)")
logT = @. log(df.newT/u"K")

interval1, interval2 = 1:29, 51:length(logT)
popt2, ci2 = bootstrap(func, nom.(logT[interval1]), nom.(logμ[interval1]), xerr=err.(logT[interval1]), yerr=err.(logμ[interval1]) , p0=[1.0, 1.0], its=10000, unc=true, xlim=[3.0, 6.5])
popt3, ci3 = bootstrap(func, nom.(logT[interval2]), nom.(logμ[interval2]), xerr=err.(logT[interval2]), yerr=err.(logμ[interval2]) , p0=[1.0, 1.0], its=10000, unc=true, xlim=[3.0, 6.5])

println(popt2, popt3)

begin
ax = subplots(figsize=(7, 4.5))[1]
myerrorbar(logT, logμ, fmt=".", capsize=3, mfc="white", mec="k", ms=7, ecolor="k", zorder=10)
myerrorbar(logT[interval1], logμ[interval1], fmt=".", capsize=3, mfc="deepskyblue", mec="k", ms=7, ecolor="k", zorder=10)
myerrorbar(logT[interval2], logμ[interval2], fmt=".", capsize=3, mfc="crimson", mec="k", ms=7, ecolor="k", zorder=10)
xlims, ylims = ax.get_xlim(), ax.get_ylim()

plot(ci2.x, func(ci2.x, popt2), zorder=5, c="C0")
plot(ci3.x, func(ci3.x, popt3), zorder=5, c="C1")
fill_between(ci2.x, ci2.c0, ci2.c1, alpha=0.3, color="C0", zorder=0)
fill_between(ci3.x, ci3.c0, ci3.c1, alpha=0.3, color="C1", zorder=0)

# legend
errorbar(1, 1, xerr=1, yerr=1, fmt=".", capsize=3,  mfc="silver", mec="k", ms=11, ecolor="k", label=L"\textrm{Messwerte}")
plot(1, 1, c="gray", label=L"\textrm{fit}")
fill_between((1, 2), 1, 1, alpha=0.3, color="gray", label=L"$2 \sigma\ \textrm{Konfidenzband}$")

text(0.04, 0.88, L"$f(x) = a_1 \cdot x + b_1$", transform=ax.transAxes, fontsize=14)
text(0.04, 0.78, L"$a_1 = 1.25(4)$", transform=ax.transAxes, fontsize=14)
text(0.04, 0.68, L"$b_1 = 3.75(17)$", transform=ax.transAxes, fontsize=14)

rect = mpl.patches.FancyBboxPatch((0.045, 0.66), 0.25, 0.27, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

text(0.705, 0.88, L"$f(x) = a_2 \cdot x + b_2$", transform=ax.transAxes, fontsize=14)
text(0.705, 0.78, L"$a_2 = -1.130(8)$", transform=ax.transAxes, fontsize=14)
text(0.705, 0.68, L"$b_2 = 15.06(5)$", transform=ax.transAxes, fontsize=14)

rect = mpl.patches.FancyBboxPatch((0.705, 0.66), 0.255, 0.27, linewidth=1.5, edgecolor="C0", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlim(xlims[0]-0.3, 6.6)
ylim(ylims)

xlabel(L"$\log(\textrm{K}/T)$")
ylabel(L"$\log\left( \mu/ \textrm{cm}^2\textrm{V}^{-1}\textrm{s}^{-1} \right)$")
legend()
tight_layout()
# savefig(string(@__DIR__, "/bilder/plot3.pdf"), bbox_inches="tight")
end


ρ .|> u"Ω*m"
RH .|> u"m^3/C"
n .|> u"m^-3"
μ .|> u"cm^2/(V*s)"
Ed

df.newT





using Plots

x = [1:10;]
y = [1:10;]

plot(x, y, fmt=:png)