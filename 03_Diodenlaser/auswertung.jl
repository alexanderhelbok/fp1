using CSV, PythonPlot, SavitzkyGolay, Intervals
using PhysicalConstants.CODATA2018: h, e, c_0

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

begin
# load data
df = CSV.read(joinpath(@__DIR__, "data/data_task1.csv"), DataFrame, header=["I", "P", "T"], skipto=2)
df.I = measurement.(df.I, 0.01)
df.P = measurement.(df.P, 0.01)

Temps = [measurement(i, 0.5) + 273.15 for i in [20, 25, 30]]

# parametrice line function with slope and x-intercept
lin(x, p) = @. p[1] * (x - p[2])
end

begin
dotcolors = ["deepskyblue", "crimson", "green"]
linecolors = ["C0", "C1", "C2"]
starts = [6, 11, 8]
λ = 670 * u"nm"
end

# plot data
begin
η = []

ax = subplots(figsize=(7, 4.5))[1]
axins = ax.inset_axes([0.1, 0.45, 0.45, 0.45])
ax2 = ax.twinx()
for Temp in 1:3
    tempdf = df[df.T .== Temp, :]
    # fit line to data
    popt, ci = bootstrap(lin, nom.(tempdf.I)[starts[Temp]:end], nom.(tempdf.P)[starts[Temp]:end], p0=[0.1, 0.1], unc=true)
    println("η = $(popt[1]), I_s = $(popt[2])")
    η = [η; popt[1]]

    ax.scatter(nom.(tempdf.I[starts[Temp]:end]), nom.(tempdf.P[starts[Temp]:end]), c=dotcolors[Temp], s=35, edgecolor="k", zorder=10)
    ax.scatter(nom.(tempdf.I[1:starts[Temp]]), nom.(tempdf.P[1:starts[Temp]]), c=dotcolors[Temp], s=35, edgecolor="k", zorder=10, marker="D")
    axins.scatter(nom.(tempdf.I[starts[Temp]:end]), nom.(tempdf.P[starts[Temp]:end]), c=dotcolors[Temp], s=45, edgecolor="k", zorder=10)
    axins.scatter(nom.(tempdf.I[1:starts[Temp]]), nom.(tempdf.P[1:starts[Temp]]), c=dotcolors[Temp], s=45, edgecolor="k", zorder=10, marker="D")

    ax.plot(ci.x, lin(ci.x, nom.(popt)), c=linecolors[Temp])
    axins.plot(ci.x, lin(ci.x, nom.(popt)), c=linecolors[Temp])
end

# hide ax2 ytick labels
ax2.set_yticklabels([])

# legend
ax.errorbar(-10, -10, fmt=".", capsize=3,  mfc="silver", mec="k", ms=18, ecolor="k", label=L"\mathrm{fit\ datapoints}")
ax.errorbar(-10, -10, fmt="D", capsize=3,  mfc="silver", mec="k", ms=10, ecolor="k", label=L"\mathrm{datapoints}")
ax.plot(-10, -10, c="gray", label=L"\mathrm{fit}")
# fill_between((1, 2), 1, 1, alpha=0.3, color="gray", label=L"$2 \sigma\ \mathrm{Konfidenzband}$")  

ax2.fill_between((-10, -20), 1, 1, alpha=0.7, color=linecolors[1], label=L"20.0(2)\ \mathrm{^\circ C}")
ax2.fill_between((-10, -20), 1, 1, alpha=0.7, color=linecolors[2], label=L"25.0(2)\ \mathrm{^\circ C}")
ax2.fill_between((-10, -20), 1, 1, alpha=0.7, color=linecolors[3], label=L"30.0(2)\ \mathrm{^\circ C}")

# axins.scatter(nom.(df.I), nom.(df.P), c="white", s=15, edgecolor="k")

# change fontsize for axins
for item in (axins.get_xticklabels() + axins.get_yticklabels())
    item.set_fontsize(12)
end

axins.set_xlim(12.1, 15.9)
axins.set_ylim(-0.5, 2.5)

# change visibility of connector lines
rectpatch, connects = ax.indicate_inset_zoom(axins, edgecolor="gray", alpha=1, zorder=0)
connects[3].set_visible(false)
connects[2].set_visible(true)

ax.set_xlim(-1, 40.4)
ax.set_ylim(-1, 24)
ax2.set_ylim(-1, 24)

ax.set_xlabel(L"$I\ (\mathrm{mA})$")
ax.set_ylabel(L"$P\ (\mathrm{mW})$")

ax.legend(loc=(0.687, 0.315))
ax2.legend(loc=(0.715, 0.04))
tight_layout()
# savefig(string(@__DIR__, "/bilder/task1.pdf"), bbox_inches="tight")
end

η 
ηd = @. η * e*λ / (h*c_0) * u"mW/mA" |> NoUnits


# ========================================
# get time to frequency conversion factor
begin
df = CSV.read(joinpath(@__DIR__, "data/T25/Current15.CSV"), DataFrame, header=["t", "V"], skipto=2)

FRS2 = 630

# find 2 maxima in the data
maxima, height = ss.find_peaks(df.V, height=0.5, distance=100)
# convert to julia array
maxima = pyconvert(Array, maxima)

dt = df.t[maxima[2]] - df.t[maxima[1]]

t_to_f(t) = @. FRS2 / dt * (t - t[1])

didx = maxima[2] - maxima[1]

idx_to_f(idx) = @. FRS2 / didx * (idx - 1)
f_to_idx(f) = @. didx / FRS2 * f + 1 |> round |> Int
# df.f = @. FRS2 / dt * (df.t - df.t[1])
end
FRS2/dt

df


currentdf = CSV.read(joinpath(@__DIR__, "data/Istepsize.csv"), DataFrame, header=["I1", "I2", "I3"], skipto=2)


function get_maximum(df; peaks=1, plt=false)
    if peaks != 1
        # get two maxima
        maxima, height = ss.find_peaks(df.V, height=0.5, distance=100)
        maxima = pyconvert(Array, maxima)
    
        mid = df.f[maxima[2]]
        start, stop = f_to_idx(mid - 50), f_to_idx(mid + 50)
        popt, ci = bootstrap(gauss, df.f[start:stop], df.V[start:stop], p0=[1.35, mid, 30, 0.06], unc=true)
    else
        maxima = argmax(df.V)
    
        mid = df.f[maxima]
        
        if f_to_idx(mid - 25) < 1
            start, stop = 1, f_to_idx(100 - mid)
        else
            start, stop = f_to_idx(mid - 25), f_to_idx(mid + 25)
        end
    
        popt, perr, ci = bootstrap(gauss, df.f[start:stop], df.V[start:stop], p0=[1.35, mid, 30, 0.06])
    end
    if plt == true
        ax = subplots(figsize=(7, 4.5))[1]
        ax.plot(df.f, df.V, label="data")
        ax.scatter(df.f[maxima], df.V[maxima], label="maxima", color="C1")
        ax.plot(ci.x, gauss(ci.x, nom.(popt)), label="fit")
    end
    return measurement(popt[2], popt[3])
end
    


intervals = [[Interval(3, 12), Interval(13,19), Interval(20, 24)],
            [Interval(3, 5), Interval(6, 12), Interval(13,21), Interval(22, 24)],
            [Interval(3, 5), Interval(6, 12), Interval(13,21), Interval(22, 24)]]

begin
lin(x, p) = @. p[1] * x + p[2]
gauss(x, p) = @. p[1] * exp(-4 * log(2) * (x - p[2])^2 / p[3]^2)

df, tempdf = DataFrame(), DataFrame()
maxima = measurement.(zeros(26), zeros(26))
Z = zeros(6000, length(currentdf.I1))
modifier = 0
ax = subplots(3, 1, figsize=(7, 4.5), sharex=true)[1]
for (j, Temp) in enumerate([25])
    for i in 1:26
        println(i)
        tempdf = CSV.read(joinpath(@__DIR__, "data/T$Temp/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
        tempdf.I = currentdf.I1[i] * ones(length(tempdf.t))
        tempdf.f = t_to_f(tempdf.t)

        Z[:, i] = tempdf.V

        if i == 7
            tempmax = get_maximum(tempdf, plt=true)
        else
            tempmax = get_maximum(tempdf)
        end
        # if i != 1
        #     if tempmax - modifier > maxima[i-1] + 100
        #         modifier = modifier + FRS2
        #     end
        # end
        maxima[i] = tempmax - modifier


        # if i == 23 || i == 22
        #     ax2 = subplots(figsize=(7, 4.5))[1]
        #     ax2.plot(tempdf.f, tempdf.V, label="data")
        #     ax2.scatter(tempdf.f[ma], tempdf.V[ma], label="maxima", color="C1", s=50, zorder=10)
        # end
        # if i > 22
        #     # ax.scatter(tempdf.I[ma], tempdf.f[ma] - 690, label="maxima", color="C1")
        #     ax.scatter(tempdf.I[1], nom.(ma) - FRS2)
        # else
        #     ax.scatter(tempdf.I[1], nom.(ma))
        #     # ax.scatter(tempdf.I[ma], tempdf.f[ma], label="maxima", color="C1")
        # end

        df = vcat(df, tempdf)
    end
    if j == 3
        maxima[1] = maxima[2]
    end

    # shift maxima
    # maxima = maxima .- nom(maximum(maxima))

    Z = Z[:, end:-1:1]

    ax[j-1].imshow(Z, aspect="auto", extent=(currentdf.I1[1], currentdf.I1[end], tempdf.f[1], tempdf.f[end]), cmap="Blues")
    ax[j-1].errorbar(currentdf.I1, nom.(maxima), yerr=err.(maxima), fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
    for (i, interval) in enumerate(intervals[j])
        # popt, ci = bootstrap(lin, currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), p0=[1., 1.], unc=true, its=100000)
        # ax[j-1].plot(ci.x, lin(ci.x, nom.(popt)), c="C$i")
        ax[j-1].errorbar(currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), fmt="o", capsize=3, mfc="C$i", mec="k", ms=7, ecolor="k")
        # ax[j-1].fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C$i")

        # println(popt)
    end
end
end

# t_to_f(df.t[end])

df
myerrorbar(currentdf.I1, maxima, fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
for (i, interval) in enumerate(intervals[1])
    imshow(Z, aspect="auto", extent=(currentdf.I1[end], currentdf.I1[1], -df.f[1], -df.f[end]), cmap="Blues")
    popt, ci = bootstrap(lin, currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), p0=[1., 1.], unc=true)
    plot(ci.x, lin(ci.x, nom.(popt)), c="C$i")
    myerrorbar(currentdf.I1[interval.first:interval.last], maxima[interval.first:interval.last], fmt="o", capsize=3, mfc="C$i", mec="k", ms=7, ecolor="k")
    fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C$i")

    println(popt)
end
maxima



df = CSV.read(joinpath(@__DIR__, "data/T30/Current1.CSV"), DataFrame, header=["t", "V"], skipto=2)
df.f = t_to_f(df.t)

mid = argmax(df.V)

start, stop = 1, f_to_idx(100 - df.f[mid])
popt, ci = bootstrap(voigt, df.f[start:stop], df.V[start:stop], p0=[1.35, mid, 30, 0.06], unc=true)

# plot data
begin   
# plot(df.t[start:stop], df.V[start:stop], label="data")
# plot(df.t[start:stop], smooth.y[start:stop], label="smoothed")
plot(df.f, df.V, label="data")
# plot(df.f[f_to_idx(600):f_to_idx(800)], df.V[f_to_idx(670):f_to_idx(720)])
# scatter(df.t[maxima], df.V[maxima], label="maxima", color="C1")
# plot(df.t, smooth.y, label="smoothed")
# plot(ci.x, voigt(ci.x, nom.(popt)), label="fit")
# plot(ci2.x, voigt(ci2.x, nom.(popt2)), label="fit smooth")
# plot(idx_to_f([1:1:6000;]), voigt(idx_to_f([1:1:6000;]), nom.([1.35, 600, 20, 0.06])), label="fit")
# plot(df.f[start:stop], df.V[start:stop])

# xlim(mid-50, mid+50)
end

# create 2d heatmap of the data
# fill matrix with data
Z = zeros(6000, length(currentdf.I1))
for i in 1:26
    df = CSV.read(joinpath(@__DIR__, "data/T25/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
    df.f = t_to_f(df.t)
    # df.V = savitzky_golay(df.V, 101, 2).y

    Z[:, i] = df.V
end
Z = Z[:, end:-1:1]
# flip Z


# plot heatmap
begin
ax = subplots(figsize=(7, 4.5))[1]
imshow(Z, aspect="auto", extent=(currentdf.I1[end], currentdf.I1[1], df.f[end], -df.f[1]), cmap="Blues")

xlabel(L"$I\ (\mathrm{mA})$")
ylabel(L"$f\ (\mathrm{GHz})$")

colorbar(label=L"$V\ (\mathrm{V})$")

tight_layout()
end



# use savgol to smooth data
smooth = savitzky_golay(df.V, 101, 2)

# fit a voigt profile to the data
voigt(x, p) = @. p[1] * exp(-2 * log(2) * (x - p[2])^2 / p[3]^2) + p[4]

popt, ci = bootstrap(voigt, df.t[start:stop], df.V[start:stop], p0=[1.35, -0.01875, 0.0003, 0.06], unc=true)
popt2, ci2 = bootstrap(voigt, df.t[start:stop], smooth.y[start:stop], p0=[1.35, -0.01875, 0.0003, 0.06], unc=true)




# voigt profile as convolution of gaussian and lorentzian
gauss(x, p) = @. p[1] * exp(-2 * log(2) * (x - p[2])^2 / p[3]^2)