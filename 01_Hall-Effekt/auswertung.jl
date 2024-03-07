using CSV, Unitful, Roots, LaTeXStrings, PythonPlot
using PhysicalConstants.CODATA2018: e, k_B

include(string(@__DIR__, "/../Source.jl"))
mpl.use("TkAgg")

# @py import matplotlib.pyplot as plt
# @py import matplotlib as mpl
# plt.style.use("Source.mplstyle")
# include("Source.jl")

# mpl.use("pgf")
# plt.rcParams["pgf.texsystem"] = "xelatex"
# plt.rc("text", usetex=true)  # enable use of LaTeX in matplotlib
# plt.rc("font", family="sans-serif", serif="Times New Roman", size=14)  # font settings
# plt.rc("text.latex", preamble="\\usepackage{mtpro2} \\usepackage{siunitx}")

# define constants
begin
    d = measurement(500, 1) * u"μm"
    B = measurement(252.6, 0.05) * u"mT"
    B0 = measurement(1.045, 0.05) * u"mT"
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
    tempdf[!, Temp] = measurement.(tempdf[!, Temp], 0.5)
    tempdf[!, Temp] = T.(tempdf[!, Temp])
end
for i in ["1a", "1b", "2a", "2b", "3a", "3b", "3aB", "3bB"]
    tempdf[!, i] = measurement.(tempdf[!, i], 0.5) * u"mV"
end

# interpolate temperature
temp = tempdf[!, "T1"]
tempdf.T1 = (tempdf.T1 + tempdf.T2)/2
tempdf.T2[1:end-1] = (temp[2:end] + tempdf.T2[1:end-1])/2

# drop last row
tempdf = tempdf[1:end-1, :]

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
UH = df.U3B .- df.U3
RH = @. UH * d / (B0 * I)

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

n = @. log(1*u"m^3"/(RH * e))
df


df.newT = @. (df.T1 + df.T2)/2

func(x, p) = @. p[1]*x + p[2]
popt, ci = bootstrap(func, nom.(1 ./df.newT[1:end]), nom.(n[1:end]), xerr=errr.(1 ./df.newT[1:end]), yerr=errr.(n[1:end]) , p0=[1.0, 1.0], its=1000, unc=true)

Ed = - 2*k_B * popt[1] * u"K" |> u"eV"
# mean = mean(n[:8])^

begin
ax = subplots(figsize=(7, 4.5))[1]
# myerrorbar(1 ./df.newT, n, fmt=".k", capsize=3, label=L"$\mathrm{data}$")
# errorbar(1 ./df.newT, n, xerr=errr.(1 ./df.newT), yerr=errr.(n), fmt=".k", capsize=3, label=L"$\mathrm{data}$")
myerrorbar([1:0.01:1.5;], [1:0.01:1.5;])
tight_layout()
end

# ax = subplots(figsize=(7, 4.5))[1]

