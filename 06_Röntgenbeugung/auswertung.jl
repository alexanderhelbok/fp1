using CSV, LaTeXStrings, GLMakie, CairoMakie, Colors
using PhysicalConstants.CODATA2018: c_0, h

include(string(@__DIR__, "/../Source.jl"))

CairoMakie.activate!()
GLMakie.activate!()

# define Constants
begin
d = 201.4

theta_to_d(theta) = @. 2*d*sin(theta*π/180)
lam_to_E(lam) = @. h*c_0/lam /u"keV"/u"pm" |> NoUnits
E_to_lam(E) = @. h*c_0/(E * u"keV") / u"pm" |> NoUnits
end

# load new data
begin
df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
backdf = CSV.read(string(@__DIR__, "/data/Ex2_hintergrund.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang) 
backdf.N = measurement.(backdf.N, sqrt.(backdf.N))

df.N .-= append!(backdf.N, zeros(length(df.N) - length(backdf.N)))
# new = df.N .- append!(backdf.N, zeros(length(df.N) - length(backdf.N)))
end

df.lam

kwargs = (yminorticksvisible = true,
    xminorticksvisible = true,
    spinewidth = 1.5,
    xminorticks = IntervalsBetween(5),
    yminorticks = IntervalsBetween(5),
    xtickwidth = 1.5,
    ytickwidth = 1.5,
    xticksize = -10,
    yticksize = -10,
    xminorticksize = -5,
    yminorticksize = -5,
    xticksmirrored = true,
    yticksmirrored = true,
    xgridvisible = false,
    ygridvisible = false,
    xgridwidth = 2,
    ygridwidth = 2,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelfont = "Times New Roman",
    yticklabelfont = "Times New Roman",
    xlabelfont = "Times New Roman",
    xlabelsize = 20,
    ylabelfont = "Times New Roman",
    ylabelsize = 20,
    xlabelpadding = 14,
    ylabelpadding = 14)



function mylegend(figure, elems, labels, pos; fig = figure[1, 1], rgs...)
	tempax = Axis(fig)
	hidedecorations!(tempax)
	hidespines!(tempax)
	leg_projetion = campixel(tempax.scene)
	@lift translate!(leg_projetion, Vec2f($(figure.scene.camera.resolution)[1]*pos[1], $(figure.scene.camera.resolution)[2]*pos[2]))
	Legend(leg_projetion, elems, labels; rgs...)
end

legargs = (labelfont = "Times New Roman", 
    labelsize = 20, 
    margin = ones(4).*18,
    patchlabelgap = 10,
    backgroundcolor = :transparent)

# plot data
begin
function inver(x)
    # return asind(x/(2*d))
    return E_to_lam(x)
end
function scali(x)
    # return theta_to_d(x)
    return lam_to_E(x)
end
function Makie.inverse_transform(::typeof(scali))
    return x -> inver(x)
end
Makie.defaultlimits(::typeof(scali)) = (-Inf, Inf)
Makie.defined_interval(::typeof(scali)) = Makie.OpenInterval(-Inf, Inf)
end

# lorentzian fit
multilorentzian(x, p) = @. p[1]/(1 + (x - p[2])^2/p[3]^2) + p[4]/(1 + (x - p[5])^2/p[6]^2) + p[7]/(1 + (x - p[8])^2/p[9]^2) + p[10] * x + p[11]

start, stop = 70, 220

popt, ci = bootstrap(multilorentzian, df.lam[start:stop], nom.(df.N)[start:stop], yerr=err.(df.N)[start:stop], p0=[500., 110., 2., 2000., 128., 2., 3500., 147., 1., -3.5, 670.], unc=true, redraw=false)

for i in [2, 5, 8]
    λ = measurement(nom(popt[i]), nom(popt[i + 1]))
    println("λ = $λ;  E_$i = $(lam_to_E(λ))")
end

begin
fig = Figure()
ax = Axis(fig[1, 1]; kwargs..., xticksmirrored = false, xlabel = L"\lambda\,\, (\mathrm{pm})", ylabel = L"N\,\, (\mathrm{1/s})", xticks = [50:25:250;], yticks = [0:1000:4000;])
ax2 = Axis(fig[1, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, xticksmirrored = false, xlabel = L"E\,\, (\mathrm{keV})", xaxisposition = :top)
# hidespines!(ax2)1800

ax.aspect = 7/2.75
ax2.aspect = 7/2.75

lines!(ax2, lam_to_E(df.lam), nom.(df.N), label="fit", linewidth=2, color=:transparent)
errorbars!(ax, df.lam, nom.(df.N), err.(df.N), label="Messwerte", whiskerwidth = 0, color=:black)
scatter!(ax, df.lam, nom.(df.N), colormap=:tab10, markersize=8)

# scatter!(ax, df.lam[start:stop], nom.(df.N)[start:stop], color=:blue, markersize=8)
lines!(ax, ci.x, multilorentzian(ci.x, nom.(popt)), color=2, colormap=:tab10, colorrange=(1, 10), linewidth=2)

# add text to plot
text!(ax, L"11.4(2)\,\, \mathrm{keV}", position=(65, 650), fontsize=18)
text!(ax, L"9.7(2)\,\, \mathrm{keV}", position=(85, 2000), fontsize=18)
text!(ax, L"8.41(7)\,\, \mathrm{keV}", position=(105, 3300), fontsize=18)

ax2.xscale = scali

xlims!(ax, df.lam[1] - 5, df.lam[end] + 5)
xlims!(ax2, inver.(ax.xaxis.attributes.limits[]))

# set lower limits
ylims!(ax, -150, nothing)

ylims!(ax2, ax.yaxis.attributes.limits[])

# legend
legelems = [[LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#1f77b4", marker = :circle, markersize = 10)], 
            LineElement(color = colorant"#ff7f0e", colorrange = (2, 2), linewidth = 2)]

mylegend(fig, legelems, ["Messdaten", "Fit"], (0.55, 0.3); legargs..., framevisible = false)
save(string(@__DIR__, "/bilder/Ex2.pdf"), fig)
fig
end
fig
 


begin
thicknessdict = Dict("Al" => 0.04, "Zn" => 0.05, "Cu" => 0.025, "Sn" => 0.025, "Ni" => 0.025)
densitydict = Dict("Al" => 2.7, "Zn" => 7.14, "Cu" => 8.96, "Sn" => 7.3, "Ni" => 8.9)
Zdict = Dict("Al" => 13, "Zn" => 30, "Cu" => 29, "Sn" => 50, "Ni" => 28)
end

# load data
begin
I0df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
I0df.N = measurement.(I0df.N, sqrt.(I0df.N))
I0df.lam = theta_to_d(I0df.ang)
end

stops, starts = [106, 60], [14, 7]
pos = [(0.7, 0.6), (0.7, 0.6)]

fig

fig = Figure(size = (800, 800))

for (i, material) in enumerate(["Cu", "Ni"])
df = CSV.read(string(@__DIR__, "/data/Ex4_$(material).csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang)
tempdf = filter(row -> row.ang in unique(df.ang), I0df)

μ = @. -log(df.N / tempdf.N)/thicknessdict[material]
# h*c_0./(df.lam[2:end] .* u"m") .|> u"eV"

μρ = @. μ*u"1/mm" / densitydict[material] / u"g/cm^3" / u"cm^2/g" |> NoUnits

model(x, p) = @. p[1] * (Zdict[material] * x)^3 + p[2]

println(argmax(μρ))
if i == 1
    interval = [1:stops[i]; 123:argmax(μρ)]
else
    interval = [1:stops[i]; 80:89]
end
stop = stops[i]

popt, ci = bootstrap(model, df.lam[interval], nom.(μρ)[interval], p0=[1e-9, 20.], yerr=err.(μρ)[interval], unc=true, redraw=false, xlim=(df.lam[argmax(μρ)] + 10, df.lam[1] - 10))
println(popt)
begin
function inver(x)
    # return asind(x/(2*d))
    return lam_to_E(x)
end
function scali(x)
    # return theta_to_d(x)
    return E_to_lam(x)
end
end

# plot data
if i == 1
    lamticks = [50:50:250;]
else
    lamticks = [50:25:200;]
end
ax = Axis(fig[i, 1]; kwargs..., xticksmirrored = false, ylabel = L"\mu/\rho\,\, (\mathrm{cm^2/g})", xlabel = L"E\,\, (\mathrm{keV})", yticks = [0:50:150;])
ax2 = Axis(fig[i, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, xlabel = L"\lambda\,\, (\mathrm{pm})", xaxisposition = :top, xticksmirrored = false, xticks=lamticks)
# hidespines!(ax2)

ax.aspect = 7/2.25
ax2.aspect = 7/2.25

Eion = lam_to_E((df.lam[argmax(μρ) + starts[i]] + df.lam[argmax(μρ)])/2 ± (df.lam[argmax(μρ) + starts[i]] - df.lam[argmax(μρ)])/2)
println(Eion)

vlines!(ax2, [df.lam[argmax(μρ)], df.lam[argmax(μρ) + starts[i]]], colormap=:tab10, linewidth=1, alpha=0.7)
fill_between!(ax, lam_to_E([df.lam[argmax(μρ) + starts[i]], df.lam[argmax(μρ)]]), [200, 200], [0, 0], color=colorant"#aadcee")

lines!(ax2, df.lam, nom.(df.N), label="fit", linewidth=2, color=:transparent)
errorbars!(ax, lam_to_E(df.lam), nom.(μρ), err.(μρ), label="Messwerte", whiskerwidth=0, color=:black)
scatter!(ax, lam_to_E(df.lam), nom.(μρ), colormap=:tab10, markersize=8)

scatter!(ax2, df.lam[interval], nom.(μρ)[interval], colormap=:tab10, colorrange=(1, 10), color=2, markersize=8)
lines!(ax2, ci.x, model(ci.x, nom.(popt)), color=4, colormap=:tab10, colorrange=(1, 10), linewidth=2)

# add text to plot
if i == 1
    text!(ax, L"E_K = 9.3(3)\,\, \mathrm{keV}", position=(10, 25), fontsize=18)
else
    text!(ax, L"E_K = 8.38(13)\,\, \mathrm{keV}", position=(8.75, 50), fontsize=18)
end

ax2.xscale = scali

xlims!(ax, lam_to_E(df.lam[end]) - 0.4, lam_to_E(df.lam[1]) + 0.4)
xlims!(ax2, inver.(ax.xaxis.attributes.limits[]))

# set lower limits
if i == 1
    ylims!(ax, 0, 165)
else
    ylims!(ax, 25, 165)
end
ylims!(ax2, ax.yaxis.attributes.limits[])

# legend
legelems = [LineElement(color = :transparent, colorrange = (2, 2), linewidth = 2),
            [LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#1f77b4", marker = :circle, markersize = 10)], 
            [LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#ff7f0e", marker = :circle, markersize = 10)], 
            LineElement(color = colorant"#d62728", colorrange = (2, 2), linewidth = 2)]

if i == 1
    axislegend(ax, legelems, [L"\textbf{Kupfer}", "Messdaten", "Messdaten Fit","Fit"], position=(0.95, 0.95); legargs..., framevisible = false)
else
    axislegend(ax, legelems, [L"\textbf{Nickel}", "Messdaten", "Messdaten Fit","Fit"], position=(0.95, 0.95); legargs..., framevisible = false)
end

save(string(@__DIR__, "/bilder/Ex4_3.pdf"), fig)
fig
end
fig = Figure(size = (800, 800))


CairoMakie.activate!()
GLMakie.activate!()



for (i, material) in enumerate(["Sn", "Zn", "Al"])
df = CSV.read(string(@__DIR__, "/data/Ex4_$(material).csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang)
tempdf = filter(row -> row.ang in unique(df.ang), I0df)

μ = @. -log(df.N / tempdf.N)/thicknessdict[material]
# h*c_0./(df.lam[2:end] .* u"m") .|> u"eV"

μρ = @. μ*u"1/mm" / densitydict[material] / u"g/cm^3" / u"cm^2/g" |> NoUnits

model(x, p) = @. p[1] * (Zdict[material] * x)^3 + p[2]

println(argmax(μρ))
# if i == 1
#     interval = [1:stops[i]; 123:argmax(μρ)]
# else
#     interval = [1:stops[i]; 80:89]
# end
# stop = stops[i]

# popt, ci = bootstrap(model, df.lam[interval], nom.(μρ)[interval], p0=[1e-9, 20.], yerr=err.(μρ)[interval], unc=true, redraw=false, xlim=(df.lam[argmax(μρ)] + 10, df.lam[1] - 10))
# println(popt)
function inver(x)
    # return asind(x/(2*d))
    return lam_to_E(x)
end
function scali(x)
    # return theta_to_d(x)
    return E_to_lam(x)
end

# plot data
if i == 1
    lamticks = [50:50:250;]
else
    lamticks = [50:25:200;]
end
ax = Axis(fig[i, 1]; kwargs..., xticksmirrored = false, ylabel = L"\mu/\rho\,\, (\mathrm{cm^2/g})", xlabel = L"E\,\, (\mathrm{keV})", yticks = [0:50:150;])
ax2 = Axis(fig[i, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, xlabel = L"\lambda\,\, (\mathrm{pm})", xaxisposition = :top, xticksmirrored = false, xticks=lamticks)
# hidespines!(ax2)

ax.aspect = 7/2.25
ax2.aspect = 7/2.25

# Eion = lam_to_E((df.lam[argmax(μρ) + starts[i]] + df.lam[argmax(μρ)])/2 ± (df.lam[argmax(μρ) + starts[i]] - df.lam[argmax(μρ)])/2)
# println(Eion)

# vlines!(ax2, [df.lam[argmax(μρ)], df.lam[argmax(μρ) + starts[i]]], colormap=:tab10, linewidth=1, alpha=0.7)
# fill_between!(ax, lam_to_E([df.lam[argmax(μρ) + starts[i]], df.lam[argmax(μρ)]]), [200, 200], [0, 0], color=colorant"#aadcee")

lines!(ax2, df.lam, nom.(df.N), label="fit", linewidth=2, color=:transparent)
errorbars!(ax, lam_to_E(df.lam), nom.(μρ).^(1/3), err.(μρ), label="Messwerte", whiskerwidth=0, color=:black)
scatter!(ax, lam_to_E(df.lam), nom.(μρ).^(1/3), colormap=:tab10, markersize=8)

# scatter!(ax2, df.lam[interval], nom.(μρ)[interval], colormap=:tab10, colorrange=(1, 10), color=2, markersize=8)
# lines!(ax2, ci.x, model(ci.x, nom.(popt)), color=4, colormap=:tab10, colorrange=(1, 10), linewidth=2)

# add text to plot
# if i == 1
#     text!(ax, L"E_K = 9.3(3)\,\, \mathrm{keV}", position=(10, 25), fontsize=18)
# else
#     text!(ax, L"E_K = 8.38(13)\,\, \mathrm{keV}", position=(8.75, 50), fontsize=18)
# end

ax2.xscale = scali

xlims!(ax, lam_to_E(df.lam[end]) - 0.4, lam_to_E(df.lam[1]) + 0.4)
xlims!(ax2, inver.(ax.xaxis.attributes.limits[]))

# set lower limits
if i == 1
    ylims!(ax, 0, 165)
elseif i == 2
    ylims!(ax, 25, 165)
else
    ylims!(ax, 0, 165)
end
ylims!(ax2, ax.yaxis.attributes.limits[])

# legend
legelems = [LineElement(color = :transparent, colorrange = (2, 2), linewidth = 2),
            [LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#1f77b4", marker = :circle, markersize = 10)], 
            [LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
                MarkerElement(color = colorant"#ff7f0e", marker = :circle, markersize = 10)], 
            LineElement(color = colorant"#d62728", colorrange = (2, 2), linewidth = 2)]

if i == 1
    axislegend(ax, legelems, [L"\textbf{Zinnober}", "Messdaten", "Messdaten Fit","Fit"], position=(0.95, 0.95); legargs..., framevisible = false)
else
    axislegend(ax, legelems, [L"\textbf{Zink}", "Messdaten", "Messdaten Fit","Fit"], position=(0.95, 0.95); legargs..., framevisible = false)
end

# save(string(@__DIR__, "/bilder/Ex4_3.pdf"), fig)
fig
end
fig=Figure()