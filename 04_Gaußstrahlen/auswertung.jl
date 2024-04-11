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
data.w1 = measurement.(data.w1/2, 2)
data.w2 = measurement.(data.w2/2, 2)
end

function w(z, p) 
    return @. p[1] * sqrt(1 + ((z - p[2])*1e4 *nom.(λ * 1e-3) / (π * p[1]^2))^2)
end
popt, ci = bootstrap(w, nom.(data.z), nom.(data.w1), xerr=err.(data.z), yerr=err.(data.w1), p0=[220., 20.], redraw=false, unc=true)
popt2, ci2 = bootstrap(w, nom.(data.z), nom.(data.w2), xerr=err.(data.z), yerr=err.(data.w2), p0=[81., 7.5], redraw=false, unc=true)

println("w1 = $(popt[1]), z1 = $(popt[2])")
println("w2_theo = $(λ*f/(π*popt[1]*u"µm") |> u"µm"), w2_exp = $(popt2[1]), z2 = $(popt2[2])")

# plot data
begin
# ax = subplots()[1]
myerrorbar(data.z, data.w1, fmt=".k", capsize=3, label="Messung 1")

xlims, ylims = xlim(), ylim()

plot(ci.x, w(ci.x, nom.(popt)), label="Fit")
fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, label=L"2 \sigma \mathrm{confidence band}")

xlabel(L"z (\mathrm{cm})")
ylabel(L"w (\mu\mathrm{m})")

xlim(xlims)
ylim(ylims)

tight_layout()
legend()
# savefig(string(@__DIR__, "/bilder/task1.pdf"), bbox_inches="tight")
end

# plot data
begin
# ax = subplots()[1]
myerrorbar(data.z, data.w2, fmt="o", label=L"\mathrm{Data}")

xlims, ylims = xlim(), ylim()

plot(ci2.x, w(ci2.x, nom.(popt2)), label=L"\mathrm{Fit}")
# plot(ci2.x, w(ci2.x, [82, 7.6]), label="Fit 1")
fill_between(ci2.x, ci2.c0, ci2.c1, alpha=0.3, label=L"2 \sigma \mathrm{confidence band}")

xlabel(L"z\ (\mathrm{cm})")
ylabel(L"w\ (\mu\mathrm{m})")

xlim(xlims)
ylim(ylims)

tight_layout()
legend()
# savefig(string(@__DIR__, "/bilder/task1.pdf"), bbox_inches="tight")
end

df

# ========================================
# get time to frequency conversion factor
begin
df = CSV.read(joinpath(@__DIR__, "data/T0001ALL.CSV"), DataFrame, header=["t", "CH1", "peak1", "CH2", "peak2"], skipto=17)
df.CH2 /= maximum(df.CH2)
df = df[1200:3200, :]
df.f = t_to_f(df.t)

FSR = c_0/(2*L) / u"GHz"|> NoUnits

# find 2 maxima in the data
maxima, height = ss.find_peaks(df.CH2, height=0.2, distance=100)
# convert to julia array
maxima = pyconvert(Array, maxima) .+ 1

dt = df.t[maxima[3]] - df.t[maxima[1]]

t_to_f(t) = @. FSR / dt * (t - t[1])
    
didx = maxima[3] - maxima[1]

idx_to_f(idx) = @. FSR / didx * (idx - 1)
f_to_idx(f) = @. didx / FSR * f + 1 |> round |> Int

# popt, perr, ci = bootstrap(pseudo_voigt, t_to_f(df.t), df.CH2, p0=[1, 0.4, 0.5, 0.1, 0.1], redraw=true)

plot(df.f, df.CH2)
scatter(df.f[maxima], df.CH2[maxima])
end

lorentzian(x, p) = @. p[1] / ( 1 + ( ( x - p[2] )/ p[3] )^2 ) + p[4]

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



df = CSV.read(joinpath(@__DIR__, "data/T0003ALL.CSV"), DataFrame, header=["t", "CH1", "peak1", "CH2", "peak2"], skipto=17)
df.CH2 /= maximum(df.CH2)
# df = df[1200:3200, :]
df.f = t_to_f(df.t)
plot(df.f, df.CH1)
plot(df.f, df.CH2)