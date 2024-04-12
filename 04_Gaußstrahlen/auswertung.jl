using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings
using PhysicalConstants.CODATA2018: c_0

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

# Constants
begin
λ = 633 * u"nm"
r = 150 * u"mm"
L = 150 * u"mm"
f = 100 * u"mm"
end

# load data
begin
data = CSV.read(string(@__DIR__, "/data/data.csv"), DataFrame, header=["z", "w1", "w2"], skipto=2)
data.z = measurement.(data.z, 0.1)
data.w1 = measurement.(data.w1/2, 3)
data.w2 = measurement.(data.w2/2, 3)
end

w(z, p) = @. p[1] * sqrt(1 + ((z - p[2])*1e4 *nom.(λ * 1e-3) / (π * p[1]^2))^2)

popt, ci = bootstrap(w, nom.(data.z), nom.(data.w1), xerr=err.(data.z), yerr=err.(data.w1), p0=[220., 20.], redraw=true, unc=true, p=0.68)
popt2, ci2 = bootstrap(w, nom.(data.z), nom.(data.w2), xerr=err.(data.z), yerr=err.(data.w2), p0=[81., 7.5], redraw=true, unc=true)

chi2 = chisq(nom.(data.w1), w(nom.(data.z), popt), sigma=err.(data.w1), pcount=2)
chi2 = chisq(nom.(data.w2), w(nom.(data.z), popt2), sigma=err.(data.w2), pcount=2)

println("w1 = $(popt[1]), z1 = $(popt[2])")
println("w2_theo = $(λ*f/(π*popt[1]*u"µm") |> u"µm"), w2_exp = $(popt2[1]), z2 = $(popt2[2]), $(popt2[2] + measurement(2.5, 0.2))")

# plot data
begin
ax = subplots(figsize=(7, 4))[1]
myerrorbar(data.z, data.w1, fmt=".", mfc="C0", mec="k", ecolor="k", ms=13, label=L"\mathrm{data}")

xlims, ylims = xlim(), ylim()

plot(ci.x, w(ci.x, nom.(popt)), label=L"\mathrm{fit}", color="C1")
fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, label=L"1 \sigma\ \mathrm{confidence\ band}", color="C0")

ax.text(0.54, 0.35, L"w(z) = w_0 \sqrt{1 + \left( \frac{(z - z_0) \lambda}{\pi w_0^2}  \right)^2}", fontsize=14, transform=ax.transAxes)
ax.text(0.54, 0.21, L"w_0 = 224(4)\ \mu\mathrm{m}", fontsize=14, transform=ax.transAxes)
ax.text(0.54, 0.09, L"z_0 = -18.45(17)\ \mathrm{cm}", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.55, 0.07), 0.4, 0.4, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlabel(L"z\ (\mathrm{cm})")
ylabel(L"w\ (\mu\mathrm{m})")

xlim(xlims)
ylim(ylims)

tight_layout()
legend()
# savefig(string(@__DIR__, "/bilder/w1.pdf"), bbox_inches="tight")
end

# plot data
begin
ax = subplots(figsize=(7, 4.4))[1]
myerrorbar(data.z, data.w2, fmt=".", mfc="C0", mec="k", ecolor="k", ms=13, label=L"\mathrm{data}")

xlims, ylims = xlim(), ylim()

plot(ci2.x, w(ci2.x, nom.(popt2)), label=L"\mathrm{fit}")
fill_between(ci2.x, ci2.c0, ci2.c1, alpha=0.3, label=L"2 \sigma\ \mathrm{confidence\ band}")

ax.text(0.29, 0.83, L"w(z) = w_0 \sqrt{1 + \left( \frac{(z - z_0) \lambda}{\pi w_0^2}  \right)^2}", fontsize=14, transform=ax.transAxes)
ax.text(0.29, 0.71, L"w_0 = 82(3)\ \mu\mathrm{m}", fontsize=14, transform=ax.transAxes)
ax.text(0.29, 0.63, L"z_0 = 7.51(15)\ \mathrm{cm}", fontsize=14, transform=ax.transAxes)
rect = mpl.patches.FancyBboxPatch((0.29, 0.62), 0.42, 0.32, linewidth=1.5, edgecolor="C1", facecolor="none", transform=ax.transAxes, boxstyle=mpl.patches.BoxStyle("Round", pad=0.02))
ax.add_patch(rect)

xlabel(L"z\ (\mathrm{cm})")
ylabel(L"w\ (\mu\mathrm{m})")

xlim(xlims)
ylim(ylims)

tight_layout()
legend(loc=(0.31, 0.32))
# savefig(string(@__DIR__, "/bilder/w2.pdf"), bbox_inches="tight")
end

df

# ========================================
# get time to frequency conversion factor
begin
df = CSV.read(joinpath(@__DIR__, "data/T0001ALL.CSV"), DataFrame, header=["t", "CH1", "peak1", "CH2", "peak2"], skipto=17)
df.CH2 /= maximum(df.CH2)
df = df[1200:3200, :]

FSR = c_0/(2*L) / u"GHz"|> NoUnits

# find 2 maxima in the data
maxima, height = ss.find_peaks(df.CH2, height=0.2, distance=100)
# convert to julia array
maxima = pyconvert(Array, maxima) .+ 1

dt = df.t[maxima[3]] - df.t[maxima[1]]

t_to_f(t) = @. FSR / dt * (t - t[1])
df.f = t_to_f(df.t)

didx = maxima[3] - maxima[1]

idx_to_f(idx) = @. FSR / didx * (idx - 1)
f_to_idx(f) = @. didx / FSR * f + 1 |> round |> Int

# popt, perr, ci = bootstrap(pseudo_voigt, t_to_f(df.t), df.CH2, p0=[1, 0.4, 0.5, 0.1, 0.1], redraw=true)

plot(df.f, df.CH2)
scatter(df.f[maxima], df.CH2[maxima])
end

lorentzian(x, p) = @. p[1] / ( 1 + ( ( x - p[2] )/ p[3] )^2 ) + p[4]

fourlorentzian(x, p) = @. p[1] / ( 1 + ( ( x - p[2] )/ p[3] )^2 ) + p[4] / ( 1 + ( ( x - p[5] )/ p[6] )^2 ) + p[7] / ( 1 + ( ( x - p[8] )/ p[9] )^2 ) + p[10]
multilorentzian(x, p) = @. p[1] / ( 1 + ( ( x - p[2] )/ p[3] )^2 ) + p[4] / ( 1 + ( ( x - p[5] )/ p[6] )^2 ) + p[7] / ( 1 + ( ( x - p[8] )/ p[9] )^2 ) + p[10]

begin
    ax = subplots(figsize=(7, 3.5))[1]
    # mid = maxima[2] - maxima[1]
    # threshhold = f_to_idx(0.5)
    start, stop = 1, 2000
    # mid
    popt, ci = bootstrap(multilorentzian, df.f[start:stop], df.CH2[start:stop], p0=[1., idx_to_f(maxima[1]), 0.006, .4, idx_to_f(maxima[2]), 0.006, 1., idx_to_f(maxima[3]), 0.006, -0.02], redraw=false, unc=true)

    println(popt)

    scatter(df.f[start:stop], df.CH2[start:stop], s=10, label=L"\mathrm{data}")

    xlims, ylims = xlim(), ylim()
    x = range(ci.x[1], ci.x[end], length=10000)
    plot(x, multilorentzian(x, nom.(popt)), label=L"\mathrm{multilorentzian\ fit}", color="C1", zorder=5)

    arrow(0.55, 0.5,  0.2, 0., length_includes_head=true, overhang=.1, head_width=0.03, head_length=0.03, transform=ax.transAxes, fc="gray")
    arrow(0.55, 0.5, -0.25, 0., length_includes_head=true, overhang=.1, head_width=0.03, head_length=0.03, transform=ax.transAxes, fc="white")
    ax.text(0.45, 0.53, L"\nu_{\mathrm{FSR}} \stackrel{!}{=} 1\ \mathrm{GHz}", fontsize=14, transform=ax.transAxes)

    arrow(0.4, 0.2,  0.1, 0., length_includes_head=true, overhang=.1, head_width=0.03, head_length=0.03, transform=ax.transAxes)
    arrow(0.4, 0.2, -0.1, 0., length_includes_head=true, overhang=.1, head_width=0.03, head_length=0.03, transform=ax.transAxes)
    ax.text(0.34, 0.23, L"0.51\ \mathrm{GHz}", fontsize=14, transform=ax.transAxes)

    arrow(0.65, 0.2,  0.1, 0., length_includes_head=true, overhang=.1, head_width=0.03, head_length=0.03, transform=ax.transAxes)
    arrow(0.65, 0.2, -0.11, 0., length_includes_head=true, overhang=.1, head_width=0.03, head_length=0.03, transform=ax.transAxes)
    ax.text(0.59, 0.23, L"0.49\ \mathrm{GHz}", fontsize=14, transform=ax.transAxes)

    xlabel(L"\nu\ (\mathrm{GHz})")
    ylabel(L"V\ \mathrm{(arb.u.)}")

    xlim(xlims)
    ylim(ylims)

    legend()
    tight_layout()
    # savefig(string(@__DIR__, "/bilder/multilorentzian.pdf"), bbox_inches="tight")
end

println("Δν1 = $(popt[5] - popt[2]), Δν2 = $(popt[8] - popt[5])")

begin
    df = CSV.read(joinpath(@__DIR__, "data/T0004.CSV"), DataFrame, header=["t", "CH1", "peak1", "CH2", "peak2"], skipto=17)
    df.CH2 /= maximum(df.CH2)
    # df = df[1:end, :]
    
    FSR = c_0/(2*L) / u"GHz"|> NoUnits
    
    # find 2 maxima in the data
    maxima, height = ss.find_peaks(df.CH2, height=0.2, distance=100)
    # convert to julia array
    maxima = pyconvert(Array, maxima) .+ 1
    
    dt = df.t[maxima[3]] - df.t[maxima[2]]
    
    t_to_f(t) = @. FSR / dt * (t - t[1])
    df.f = t_to_f(df.t)
    
    didx = maxima[3] - maxima[2]
    
    idx_to_f(idx) = @. FSR / didx * (idx - 1)
    f_to_idx(f) = @. didx / FSR * f + 1 |> round |> Int

    popt, perr, ci = bootstrap(fourlorentzian, df.f, df.CH2, p0=[1., idx_to_f(maxima[1]), 0.006, 1., idx_to_f(maxima[2]), 0.006, 1., idx_to_f(maxima[3]), 0.006, 1., idx_to_f(maxima[4]), 0.006, 0.], redraw=true)
    
    # calculate finesse
    R = 0.98
    F1 = FSR / (2 * popt[3])
    F2 = FSR / (2 * popt[6])
    F3 = FSR / (2 * popt[9])
    F4 = FSR / (2 * popt[12])
    Ftheo = π * sqrt(R) / (1 - R)
    println("F1 = $F1, F2 = $F2, F3 = $F3, F4 = $F4, Ftheo = $Ftheo")

    scatter(df.f, df.CH2)
    plot(ci.x, fourlorentzian(ci.x, nom.(popt)), color="C1")
    # scatter(df.f[maxima], df.CH2[maxima])
end

popt

for mid in maxima
    threshhold = f_to_idx(0.3)
    start, stop = mid - threshhold, mid + threshhold
    # mid
    popt, perr, ci = bootstrap(lorentzian, df.f[start:stop], df.CH2[start:stop], p0=[1, idx_to_f(mid), 0.006, -0.02], redraw=false)

    # begin
    scatter(df.f[start:stop], df.CH2[start:stop])
    # plot(x, lorentzian(x, [1, 0.43, 0.006, -0.02]))
    plot(ci.x, lorentzian(ci.x, nom.(popt)), label="lorentzian")
    # legend()
end



tempdf = CSV.read(joinpath(@__DIR__, "data/T0003ALL.CSV"), DataFrame, header=["t", "CH1", "peak1", "CH2", "peak2"], skipto=17)
df = tempdf
arr = [df.CH2[end-1299:end]; df.CH2[1:end-1300]]
df[1:1000, :] = tempdf[end-999:end, :]
df[2001:end, :] = tempdf[1:end-2000, :]
df.CH2 /= maximum(df.CH2)

# apply modulo to data
df = [df[1200:end, :], df[1:1200, :]]

x = [1:1:length(df.t);]

# df = df[1200:3200, :]
df.f = t_to_f(df.t)
plot(df.f, df.CH1)
plot(x, df.CH2)
plot(x/4000, arr)