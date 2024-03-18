using CSV, LaTeXStrings, GLMakie, Colors, CairoMakie
using CSV, LaTeXStrings, PythonPlot

include(string(@__DIR__, "/../Source.jl"))

# @py import matplotlib as mpl

# mpl.use("pgf")
# mpl.use("TkAgg")

# load data
calibrationdf = CSV.read(joinpath(@__DIR__, "data0.csv"), DataFrame, header=["gain", "mvpp", "mw"])

fun(x, p) = @. p[1] * x^2

calipopt, ci = bootstrap(fun, calibrationdf.mvpp, calibrationdf.mw, p0=[1.], unc=true)

chisq(calibrationdf.mw, fun(calibrationdf.mvpp, nom.(calipopt)), pcount=1)

println(calipopt)

# plot data
begin
fig = Figure()
ax = Axis(fig[1, 1]; kwargs..., xlabel = L"U_\mathrm{pp}\,\, (\mathrm{mVpp})", ylabel = L"P_\mathrm{RF}\,\, (\mathrm{W})")
scatter!(ax, calibrationdf.mvpp, calibrationdf.mw, color = :deepskyblue, label="data")
a = lines!(ax, calibrationdf.mvpp, fun(calibrationdf.mvpp, nom.(calipopt)), colormap=:tab10, label="fit")

legelems = [MarkerElement(color = :deepskyblue, marker = :circle, markersize = 15), 
            LineElement(color = :blue, colorrange = (2, 2), linewidth = 2)]

mylegend(fig, legelems, ["data", "fit"], 0.1, 0.6; legargs..., framevisible = false)

fig
end


Pout = measurement(0.831, 0.005)

# load data
tempdf = CSV.read(joinpath(@__DIR__, "data1.csv"), DataFrame, header=["U", "a", "b"], skipto=2)
tempdf.a = measurement.(tempdf.a, 0.01) ./ Pout
tempdf.b = measurement.(tempdf.b, 0.01) ./ Pout

mvpp_to_mW(U) = @. nom.(calipopt[1]) * U^2

tempdf.U2 = mvpp_to_mW(tempdf.U)

func(x, p) = @. p[1] * sin(pi/2 * x/p[2])^2 + p[3]
func(x, p) = @. p[1] * sin(p[2] * sqrt(abs(x)))^2 + p[3]

popt, ci = bootstrap(func, tempdf.U, nom.(tempdf.a), yerr=err.(tempdf.a), p0=[-0.9, 480., 1.], unc=true, xlim=(0, 700))
popt2, ci2 = bootstrap(func, tempdf.U, nom.(tempdf.b), yerr=err.(tempdf.b), p0=[1, 500, 0.0], unc=true, xlim=(0, 700))
# popt, ci = bootstrap(func, tempdf.U2, tempdf.b, p0=[0.8, 2., 0.03])

println(chisq(tempdf.a, func(tempdf.U, nom.(popt)), pcount=3), chisq(tempdf.b, func(tempdf.U, nom.(popt2)), pcount=3))

println((popt[2] + popt2[2]))

tempdf.a

ci

function inver(x)
    return sqrt(abs(x) ./ nom.(calipopt[1]))
end
function scali(x)
    if x > 0
        return mvpp_to_mW(x)
    else
        return x
    end
end
function Makie.inverse_transform(::typeof(scali))
    return x -> sqrt(abs(x) ./ nom.(calipopt[1]))
end
Makie.defaultlimits(::typeof(scali)) = (-Inf, Inf)
Makie.defined_interval(::typeof(scali)) = Makie.OpenInterval(-Inf, Inf)

tempdf

# plot data
begin
fig = Figure()
ax = Axis(fig[1, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, yminorticksvisible = false, xticksmirrored = false, xlabel = L"U_\mathrm{pp}\,\, (\mathrm{mVpp})", xaxisposition = :top, xticks = [0, 200, 300, 400, 500, 600], xminorticksvisible = false, xlabelsize=20, xticklabelsize=16)
ax3 = Axis(fig[1, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, yminorticksvisible = false, xticksmirrored = false, xaxisposition = :top, xticks = [200:20:640;], xminorticksvisible = false, xtickwidth = 1, xticksize = -8, xticklabelsvisible=false)
ax2 = Axis(fig[1, 1]; kwargs..., xticksmirrored = false, xlabel = L"P_\mathrm{RF}\,\, (\mathrm{W})", ylabel = L"\varepsilon", xticks = [0:0.2:1;], yticks = [0:0.2:1;], xlabelsize = 20, ylabelsize=20, xticklabelsize=16, yticklabelsize=16)
hidespines!(ax2)
hidespines!(ax3)

fill_between!(ax, ci2.x, ci2.c0, ci2.c1, color=colorant"#a4d0ef")
# fill_between!(ax, ci.x, ci.c0, ci.c1, color=colorant"#ffd4ae")

lines!(ax, ci2.x, func(ci2.x, nom.(popt2)), colormap=:tab10, label="fit")
# lines!(ax, tempdf.U, func(tempdf.U, nom.(popt)), colormap=:tab10, color=colorant"#e66c00", label="fit")

# errorbars!(ax2, tempdf.U2, nom.(tempdf.a), err.(tempdf.a), color=:black, label="a", whiskerwidth = 0)
scatter!(ax2, tempdf.U2, nom.(tempdf.b), colormap = :tab10, markersize=10)
# errorbars!(ax2, tempdf.U2, nom.(tempdf.a), err.(tempdf.b), color=:black, label="b", whiskerwidth = 0)
scatter!(ax2, tempdf.U2, nom.(tempdf.a), colormap = :tab10, markersize=10)
xlims!(ax, -0.04, inver(0.95))
xlims!(ax3, -0.04, inver(0.95))
xlims!(ax2, -0.04, 0.95)

ylims!(ax, -0.05, 1.05)
ylims!(ax2, -.05, 1.05)
ylims!(ax3, -.05, 1.05)

ax.xscale = scali
ax3.xscale = scali

legelems = [[LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#ff7f0e", marker = :circle, markersize = 10)], 
            [LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#1f77b4", marker = :circle, markersize = 10)],
            LineElement(color = colorant"#1f77b4", colorrange = (2, 2), linewidth = 2),
            PolyElement(color = colorant"#a4d0ef", points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])]

mylegend(fig, legelems, ["0. Order", "1. Order", "Fit", L"2\sigma\,\, \mathrm{Band}"], 0.5, 0.1; legargs..., framevisible = false)
# save(string(@__DIR__, "/bilder/fit.pdf"), fig)
fig
end
fig
# fit data

print(popt2[1], mvpp_to_mW(popt[2])*2)

GLMakie.activate!()
CairoMakie.activate!()

calipopt
