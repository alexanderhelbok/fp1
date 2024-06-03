using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

# increase axes label pad using rcParams
mpl.rcParams["axes.labelpad"] = 10 

# load data delimeter is tab
begin
data = CSV.read(string(@__DIR__, "/data/data.csv"), DataFrame, skipto = 2, header = ["lam", "chi", "loggf", "fit_lam", "Wf"])
# filter out rows containing NaN
data = data[.!isnan.(data.fit_lam), :]
data = data[.!isnan.(data.Wf), :]
# filter data where Wf is equal fit_lam
data = data[data.Wf .!= data.fit_lam, :]
# filter data where fit_lam - lam is less than 0.15
data = data[abs.(data.fit_lam .- data.lam) .< 0.15, :]
# filter data where Wf is less than 0.12
data = data[data.Wf .< 0.12, :]
end

x = @. log(10, data.fit_lam) + data.loggf + log(10, Ni) - 5040/T*data.chi
# pop a single datapoint by index
data = data[1:end .!= findall(x .< -4.5)[2], :]

y = @. log(10, data.Wf/data.fit_lam^2) - data.loggf

model(x, p) = @. p[1] * x + p[2]

popt, ci = bootstrap(model, data.chi, y, p0=[5000., -15.5], unc=true, redraw=true, its=10000)

# calculate chisq
chisq(y, model(data.chi, popt))

T = -5040/popt[1]
Ni = 10^(popt[2])*100

# plot
begin
ax = subplots(figsize = (7, 4.5))[1]
# scatter([1:length(data.lam);], y, s=10, label = L"\mathrm{Data}")
scatter(data.chi, y, s=10, label = L"\mathrm{Data}", c = "k", zorder = 5)
# plot data where Wf > 0.12
xlims, ylims = ax.get_xlim(), ax.get_ylim()

plot(ci.x, model(ci.x, popt), color = "C1", label = L"\mathrm{Fit}", zorder = 10)
fill_between(ci.x, ci.c0, ci.c1, color = "C1", alpha = 0.2, label = L"2\sigma\ \mathrm{confidence\ band}", zorder = 8)

text(0.54, 0.86, L"f(\chi) = - \frac{5040}{T}\chi + \log\left( \frac{N_{\mathrm{Fe}}}{N_{\mathrm{H}}} \right)", fontsize=14, transform=ax.transAxes)
text(0.54, 0.74, L"T = 6100(300)\ K", fontsize=14, transform=ax.transAxes)
text(0.54, 0.62, L"\frac{N_{\mathrm{Fe}}}{N_{\mathrm{H}}} = 0.011(4)\ \%", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.55, 0.58), 0.4, 0.36, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlim(xlims)
ylim(ylims)

xlabel(L"\chi\ (\mathrm{eV})")
ylabel(L"\log \left( \frac{EW_{\lambda}}{\lambda^2} \right) - \log(\mathrm{gf})")

legend(loc="lower left")
tight_layout()
# savefig(string(@__DIR__, "/bilder/fit.pdf"))
end

x = @. log(10, data.fit_lam) + data.loggf + log(10, Ni) - 5040/T*data.chi
y = @. log(10, data.Wf/data.fit_lam)

# plot 
begin
ax = subplots(figsize = (7, 3.1))[1]
# 1:1 axis aspect ratio
ax.set_aspect("equal")

scatter(nom.(x), y, label = L"\mathrm{Data}", color = "black", s = 10)
# scatter data where x < -5
scatter(nom.(x[x .< -4.5]), y[x .< -4.5], label = L"\mathrm{Excluded\ data}", color = "crimson", s = 10)

xlims, ylims = ax.get_xlim(), ax.get_ylim()

# plot 45° line
# plot([minimum(nom.(-10)), maximum(nom.(0))], [minimum(nom.(-10)), maximum(nom.(0))], color = "black", linestyle = "--", label = "45°")

xlim(xlims)
ylim(ylims)

xlabel(L"\log \left( \frac{EW_{\lambda}}{\lambda} \right) + \log(\mathrm{gf}) + \log \left( \frac{N_{\mathrm{Fe}}}{N_{\mathrm{H}}} \right) - \frac{5040}{T} \chi")
ylabel(L"\log \left( \frac{EW_{\lambda}}{\lambda} \right)")

legend()
tight_layout()

# savefig(string(@__DIR__, "/bilder/exclude.pdf"))
end

# plot
begin
y2 = y[x .> -4.5]
x2 = x[x .> -4.5]
ax = subplots()[1]
# 1:1 axis aspect ratio
ax.set_aspect("equal")

scatter(nom.(x2), y2, label = L"\mathrm{Data}", color = "black", s = 10)
xlims, ylims = ax.get_xlim(), ax.get_ylim()
# scatter data where x > -3.3
# mark region using fill_between
fill_betweenx([-10, 0], -3.35, color = "crimson", alpha = 0.2)

text(0.36, .14, L"\mathrm{Linear\ regime}", ha = "right", color = "k", fontsize=14, transform=ax.transAxes)
text(0.36, .05, L"\propto N", color = "k", ha = "right", fontsize=14, transform=ax.transAxes)

text(0.43, .14, L"\mathrm{Saturation\ regime}", color = "crimson", fontsize=14, transform=ax.transAxes)
text(0.43, .05, L"\propto \sqrt{\ln N}", color = "crimson", fontsize=14, transform=ax.transAxes)
scatter(nom.(x2[x2 .> -3.35]), y2[x2 .> -3.35], label = L"\mathrm{Excluded\ data}", color = "crimson", s = 10)

# plot 45° line
plot([-10, 0], model([-10, 0], [1, -1.7]), color = "gray", linestyle = "--", label = L"45°\ \mathrm{line}", zorder = 0)

xlim(xlims)
ylim([ylims[0]-.1, ylims[1]]) 

xlabel(L"\log \left( \frac{EW_{\lambda}}{\lambda} \right) + \log(\mathrm{gf}) + \log \left( \frac{N_{\mathrm{Fe}}}{N_{\mathrm{H}}} \right) - \frac{5040}{T} \chi")
ylabel(L"\log \left( \frac{EW_{\lambda}}{\lambda} \right)")

legend(loc=(0.58, 0.3))
tight_layout()
# savefig(string(@__DIR__, "/bilder/exclude2.pdf"))
end

data2 = data[-4.5 .< x .< -3.35, :]

y = @. log(10, data2.Wf/data2.fit_lam^2) - data2.loggf

popt2, ci2 = bootstrap(model, data2.chi, y, p0=[5000., -15.5], unc=true, redraw=true)

# calculate chisq
chisq(y, model(data2.chi, popt))

T = -5040/popt2[1]
Ni = 10^(popt2[2])*100

# plot
begin
ax = subplots(figsize = (7, 4.5))[1]
# scatter([1:length(data2.lam);], y, s=10, label = L"\mathrm{Data2}")
scatter(data2.chi, y, s=10, label = L"\mathrm{Data}", c = "k", zorder = 0)
# plot data2 where Wf > 0.12
# scatter(data2.chi[data2.Wf .> 0.12], y[data2.Wf .> 0.12], s=10, label = "Data where W_f > 0.12")
xlims, ylims = ax.get_xlim(), ax.get_ylim()

plot(ci2.x, model(ci2.x, popt2), color = "C1", label = L"\mathrm{Fit}")
plot(ci.x, model(ci.x, popt), color = "gray", ls = "dashed", label = L"\mathrm{old\ Fit}", alpha = 0.5, zorder = 0)
fill_between(ci2.x, ci2.c0, ci2.c1, color = "C1", alpha = 0.2, label = L"2\sigma\ \mathrm{confidence\ band}")

text(0.54, 0.86, L"f(\chi) = - \frac{5040}{T}\chi + \log\left( \frac{N_{\mathrm{Fe}}}{N_{\mathrm{H}}} \right)", fontsize=14, transform=ax.transAxes)
text(0.54, 0.74, L"T = 5520(110)\ K", fontsize=14, transform=ax.transAxes)
text(0.54, 0.62, L"\frac{N_{\mathrm{Fe}}}{N_{\mathrm{H}}} = 0.043(7)\ \%", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.55, 0.58), 0.4, 0.36, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlim(xlims)
ylim(ylims)

xlabel(L"\chi\ (\mathrm{eV})")
ylabel(L"\log \left( \frac{EW_{\lambda}}{\lambda^2} \right) - \log(\mathrm{gf})")

legend(loc="lower left")
tight_layout()
# savefig(string(@__DIR__, "/bilder/fit2.pdf"))
end