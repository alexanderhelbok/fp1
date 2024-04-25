using CSV, LaTeXStrings, GLMakie, CairoMakie, Colors
using PhysicalConstants.CODATA2018: c_0, h

include(string(@__DIR__, "/../Source.jl"))

GLMakie.activate!()

# define Constants
begin
d = 201.4

end

theta_to_d(theta) = @. 2*d*sin(theta*π/180)

# load new data
begin
df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
backdf = CSV.read(string(@__DIR__, "/data/Ex2_hintergrund.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang) 
backdf.N = measurement.(backdf.N, sqrt.(backdf.N))

# df.N .-= append!(backdf.N, zeros(length(df.N) - length(backdf.N)))
new = df.N .- append!(backdf.N, zeros(length(df.N) - length(backdf.N)))
end

kwargs = (yminorticksvisible = true,
    xminorticksvisible = true,
    spinewidth = 2,
    xminorticks = IntervalsBetween(5),
    yminorticks = IntervalsBetween(5),
    xtickwidth = 2,
    ytickwidth = 2,
    xticksize = -14,
    yticksize = -14,
    xminorticksize = -7,
    yminorticksize = -7,
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

# plot data

function inver(x)
    return asind(x/(2*d))
    # return x/(2*d*(π/180) * 1e12)
end
function scali(x)
    return theta_to_d(x)
    # return @. 2*d*(x*π/180) * 1e12
end
function Makie.inverse_transform(::typeof(scali))
    return x -> asind(x/(2*d))
    # return x -> x/(2*d*(π/180) * 1e12)
end
Makie.defaultlimits(::typeof(scali)) = (-Inf, Inf)
Makie.defined_interval(::typeof(scali)) = Makie.OpenInterval(-Inf, Inf)


begin
fig = Figure()
ax = Axis(fig[1, 1]; kwargs..., xticksmirrored = false, xlabel = L"\lambda\,\, (\mathrm{pm})", ylabel = L"N (\mathrm{1/s})", xticks = [50:25:250;], yticks = [0:500:3500;])
ax2 = Axis(fig[1, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, xticksmirrored = false, xlabel = L"\theta\,\, (\mathrm{°})", xaxisposition = :top, xticks = [5:5:35;])
# hidespines!(ax2)

ax.aspect = 7/4.5
ax2.aspect = 7/4.5

errorbars!(ax, df.lam, nom.(df.N), err.(df.N), label="Messwerte")
scatter!(ax2, df.ang, nom.(df.N), color="black", markersize=4)
# scatter!(ax, [1:1:80;], [1:1:80;], color="black", s=30, zorder=10)
lines!(ax, df.lam, nom.(df.N), colormap=:tab10, label="fit")
lines!(ax, df.lam, nom.(new), colormap=:tab10, label="fit")
# fill_between!(ax, df.lam, nom.(df.N) - err.(df.N), nom.(df.N) + err.(df.N))

ax2.xscale = scali

xlims!(ax, df.lam[1] - 5, df.lam[end] + 5)
xlims!(ax2, inver.(ax.xaxis.attributes.limits[]))

# set lower limits
ylims!(ax2, -150, nothing)

ylims!(ax, ax2.yaxis.attributes.limits[])


fig
end


# load data
begin
Voltages = [15:2.5:30;]
for (i, V) in enumerate(Voltages[1])
    df = CSV.read(string(@__DIR__, "/data/Ex3_$i.csv"), DataFrame, header=["ang", "N"], skipto=3)
    df.N = measurement.(df.N, sqrt.(df.N))
    df.lam = theta_to_d(df.ang) 
end
end

plot()








model(x, p) = @. p[1]*exp(-p[2]*(x - p[3])) + p[4]




begin
fig = Figure()
ax = Axis(fig[1, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, yminorticksvisible = false, xticksmirrored = false, xlabel = L"U_\mathrm{pp}\,\, (\mathrm{mVpp})", xaxisposition = :top, xminorticksvisible = false, xlabelsize=20, xticklabelsize=16)
# ax3 = Axis(fig[1, 1]; kwargs..., yticklabelsvisible = false, yticksvisible = false, yminorticksvisible = false, xticksmirrored = false, xaxisposition = :top, xticks = [200:20:640;], xminorticksvisible = false, xtickwidth = 1, xticksize = -8, xticklabelsvisible=false)
# ax2 = Axis(fig[1, 1]; kwargs..., xticksmirrored = false, xlabel = L"P_\mathrm{RF}\,\, (\mathrm{W})", ylabel = L"\varepsilon", xticks = [0:0.2:1;], yticks = [0:0.2:1;], xlabelsize = 20, ylabelsize=20, xticklabelsize=16, yticklabelsize=16)
# hidespines!(ax2)
# hidespines!(ax3)

for df in [dfAl, dfZn]
    for lam in unique(df.lam)
        tempdf = filter(row -> row.lam == lam, df)
        # popt, ci = bootstrap(model, tempdf.V, nom.(tempdf.N), yerr=err.(tempdf.N), p0=[0.5, 0.1, 0.5, 0.1], unc=true, xlim=(0, 700))
        errorbars!(ax, tempdf.V, nom.(tempdf.N), err.(tempdf.N), label="Messwerte")
        # scatter!(ax, tempdf.V, tempdf.N, color="black", s=30, zorder=10)
        # lines!(ax, tempdf.V, model(tempdf.V, nom.(popt)), colormap=:tab10, label="fit")
    end
end

# fill_between!(ax, ci2.x, ci2.c0, ci2.c1, color=colorant"#a4d0ef")
# # fill_between!(ax, ci.x, ci.c0, ci.c1, color=colorant"#ffd4ae")

# lines!(ax, ci2.x, func(ci2.x, nom.(popt2)), colormap=:tab10, label="fit")
# # lines!(ax, tempdf.U, func(tempdf.U, nom.(popt)), colormap=:tab10, color=colorant"#e66c00", label="fit")

# # errorbars!(ax2, tempdf.U2, nom.(tempdf.a), err.(tempdf.a), color=:black, label="a", whiskerwidth = 0)
# scatter!(ax2, tempdf.U2, nom.(tempdf.b), colormap = :tab10, markersize=10)
# # errorbars!(ax2, tempdf.U2, nom.(tempdf.a), err.(tempdf.b), color=:black, label="b", whiskerwidth = 0)
# scatter!(ax2, tempdf.U2, nom.(tempdf.a), color=colorant"#ff7f0e", markersize=10)
# xlims!(ax, -0.04, inver(0.95))
# xlims!(ax3, -0.04, inver(0.95))
# xlims!(ax2, -0.04, 0.95)

# ylims!(ax, -0.05, 1.05)
# ylims!(ax2, -.05, 1.05)
# ylims!(ax3, -.05, 1.05)

# ax.xscale = scali
# ax3.xscale = scali

# legelems = [[LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
#                 MarkerElement(color = colorant"#ff7f0e", marker = :circle, markersize = 10)], 
#             [LineElement(color = :black, colorrange = (2, 2), linewidth = 2, points = Point2f[(0.5, 0), (0.5, 1)]),
#                 MarkerElement(color = colorant"#1f77b4", marker = :circle, markersize = 10)],
#             LineElement(color = colorant"#1f77b4", colorrange = (2, 2), linewidth = 2),
#             PolyElement(color = colorant"#a4d0ef", points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])]

# mylegend(fig, legelems, ["0. Order", "1. Order", "Fit", L"2\sigma\,\, \mathrm{Band}"], 0.5, 0.1; legargs..., framevisible = false)
# save(string(@__DIR__, "/bilder/fit.pdf"), fig)
fig
end
fig

df