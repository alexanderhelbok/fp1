using CSV, Unitful, Roots, LaTeXStrings, PythonPlot
using PhysicalConstants.CODATA2018: e, k_B

include(string(@__DIR__, "/../Source.jl"))
@py import matplotlib as mpl

plt.rc("text.latex", preamble="\\usepackage{mtpro2} \\usepackage{siunitx} \\usepackage{amsmath}")
# mpl.use("pgf")
# mpl.use("TkAgg")
plt.rcParams["pgf.texsystem"] = "lualatex"

# define constants
begin
    d = measurement(500, 1) * u"μm"
    B = measurement(252.6, 0.05) * u"mT"
    # B0 = measurement(1.045, 0.05) * u"mT"
    I = measurement(500, 1) * u"μA"
end

function T(R) 
    temp = @. 509/8 - 98/(2*R) + 43/30 * sqrt(R) + 161/113 * R^(11/10)
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

# drop last row
tempdf = tempdf[1:end-1, :]
# drop first 4 rows
tempdf = tempdf[5:end, :]

begin
df = DataFrame()
df.T1 = tempdf.T1
df.T2 = tempdf.T2
df.U1 = (tempdf."1a" - tempdf."1b")/2
df.U2 = (tempdf."2a" - tempdf."2b")/2
df.U3 = (tempdf."3a" - tempdf."3b")/2
df.U3B = (tempdf."3aB" - tempdf."3bB")/2
end
end

# ========== RH ==========
begin
UH = df.U3B .- df.U3
RH = @. UH * d / (B * I)

function f(x)
    x = nom.(x)
    fun(unbekannt) = @. cosh(log(2) / unbekannt * (x - 1)/(x+1)) - 0.5 * exp(log(2)/unbekannt)
    return find_zero.(fun, 1)
end

x = @. abs(df.U1/df.U2)

F = [f(X) for X in x]

ρ = @. pi * d/log(2) * F * (abs(df.U1) + abs(df.U2)) / (2 * I)

RH = -RH
μ = RH ./ ρ

n = @. 1/(RH * e)
logn = @. log(n*u"m^3")
# end
df


df.newT = @. (df.T1 + df.T2)/2

func(x, p) = @. p[1]*x + p[2]

interval = 35:length(df.newT)-8
popt, ci = bootstrap(func, 
                    nom.(1 ./df.newT[interval]), 
                    nom.(logn[interval]), 
                    xerr=myerr.(1 ./df.newT[interval]), 
                    yerr=myerr.(logn[interval]) , 
                    p0=[1.0, 1.0], 
                    its=1000, 
                    unc=true,
                    xlim=[0, 1])

Ed = - 2*k_B * popt[1] * u"K" |> u"eV"

# mean of first 6 values
const_model(x, p) = @. p[1] *x/x
mean1, ccc = bootstrap(const_model, nom.(1 ./df.newT[1:6]), nom.(n[1:6]), xerr=myerr.(1 ./df.newT[1:6]), yerr=myerr.(n[1:6]), its=1000, unc=true, p0=[1e28])
println("Ed = ", Ed)
println("NA - ND = ", mean1[1] * u"1/m^3")
# end

# n
# ========== Plot1 ==========
# begin
ax = subplots(figsize=(7, 4.5))[1]
myerrorbar(1 ./ df.newT, logn, fmt=".k", capsize=3, label=L"$\mathrm{data}$")
myerrorbar(1 ./ df.newT[interval], logn[interval], fmt=".r", capsize=3, label=L"$\mathrm{data}$")
xlims, ylims = ax.get_xlim(), ax.get_ylim()

plot(ci.x, func(ci.x, nom.(popt)), "C0", label=L"$\mathrm{fit}$")
fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C0", label=L"$2 \sigma\ \mathrm{confidence\ band}$")

axhline(log(nom(mean1[1])))
fill_between(xlims, log(nom(mean1[1]) + myerr(mean1[1])), log(nom(mean1[1]) - myerr(mean1[1])), alpha=0.3, color="C1", label=L"$\mathrm{mean\ of\ first\ 6\ values}$")

xlim(xlims)
ylim(ylims)

xlabel(L"$1/T\ (\textrm{1/K})$")
ylabel(L"$\log(n\cdot \textrm{m}^3)$")
legend()
tight_layout()
end

# ========== Plot2 =========
begin
ax = subplots(figsize=(7, 4.5))[1]
myerrorbar(1 ./ df.newT, μ, fmt=".k", capsize=3, label=L"$\mathrm{data}$")

xlabel(L"$1/T\ (\textrm{1/K})$")
ylabel(L"$\mu\ ()$")
legend()
tight_layout()
end


# ========== Plot3 =========
logμ = @. log(ustrip(μ))
logT = @. log(ustrip(df.newT))

interval1, interval2 = 1:29, 51:length(logT)
popt2, ci2 = bootstrap(func, nom.(logT[interval1]), nom.(logμ[interval1]), xerr=myerr.(logT[interval1]), yerr=myerr.(logμ[interval1]) , p0=[1.0, 1.0], its=1000, unc=true, xlim=[3.0, 6.5])
popt3, ci3 = bootstrap(func, nom.(logT[interval2]), nom.(logμ[interval2]), xerr=myerr.(logT[interval2]), yerr=myerr.(logμ[interval2]) , p0=[1.0, 1.0], its=10000, unc=true, xlim=[3.0, 6.5])

println(popt2, popt3)

begin
ax = subplots(figsize=(7, 4.5))[1]
myerrorbar(logT, logμ, fmt=".k", capsize=3, label=L"$\mathrm{data}$")
myerrorbar(logT[interval1], logμ[interval1], fmt=".b", capsize=3, label=L"$\mathrm{data}$")
myerrorbar(logT[interval2], logμ[interval2], fmt=".r", capsize=3, label=L"$\mathrm{data}$")
xlims, ylims = ax.get_xlim(), ax.get_ylim()

plot(ci2.x, func(ci2.x, popt2), label=L"$\mathrm{fit}$", zorder=10)
plot(ci3.x, func(ci3.x, popt3), label=L"$\mathrm{fit}$", zorder=10)
fill_between(ci2.x, ci2.c0, ci2.c1, alpha=0.3, color="C0", label=L"$2 \sigma\ \mathrm{confidence\ band}$")
fill_between(ci3.x, ci3.c0, ci3.c1, alpha=0.3, color="C1", label=L"$2 \sigma\ \mathrm{confidence\ band}$")

text(0.05, 0.88, L"$f(x) = a \cdot x + b$", transform=ax.transAxes, fontsize=14)
text(0.05, 0.78, L"$a = 1.12(3)$", transform=ax.transAxes, fontsize=14)
text(0.05, 0.68, L"$b = -11.77(14)$", transform=ax.transAxes, fontsize=14)

rect = mpl.patches.FancyBboxPatch((0.05, 0.66), 0.24, 0.27, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

text(0.72, 0.88, L"$f(x) = a \cdot x + b$", transform=ax.transAxes, fontsize=14)
text(0.72, 0.78, L"$a = -1.129(8)$", transform=ax.transAxes, fontsize=14)
text(0.72, 0.68, L"$b = -1.07(5)$", transform=ax.transAxes, fontsize=14)

rect = mpl.patches.FancyBboxPatch((0.715, 0.66), 0.24, 0.27, linewidth=1.5, edgecolor="C0", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlim(xlims[0]-0.1, 6.5)
ylim(ylims)

xlabel(L"$\log(\textrm{K}/T)$")
ylabel(L"$\log(\mu)$")
legend()
tight_layout()
end

