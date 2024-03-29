using CSV, PythonPlot, SavitzkyGolay, Intervals, LaTeXStrings
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
a = [1, 2, 3, 2, 3, 4]

(max, max_idx) = findmax(df.V)


currentdf = CSV.read(joinpath(@__DIR__, "data/Istepsize.csv"), DataFrame, header=["I1", "I2", "I3"], skipto=2)


function get_maximum(df, height; plt=false)
    # df.V = savitzky_golay(df.V, 201, 2).y
    temparr = []
    # normalize data
    offset = 0.01 / maximum(df.V)
    gauss(x, p) = @. p[1] * exp(-2 * log(2) * (x - p[2])^2 / p[3]^2) - offset

    df.V /= maximum(df.V)
    # println(mean(df.V))
    if mean(df.V) < 0.5
        # get two maxima
        max, height = ss.find_peaks(df.V, height=height, distance=f_to_idx(40))
        max = pyconvert(Array, max) .+ 1

        if plt
            ax = subplots(figsize=(7, 4.5))[1]
            ax.plot(df.f, df.V, label="data")
            ax.scatter(df.f[max], df.V[max], label="maxima", color="C1")
        end
        
        threshhold = f_to_idx(25)
        for mid in max
            # start, stop = f_to_idx(mid - 50), f_to_idx(mid + 50)
                
            if mid - threshhold < 1
                start, stop = 1, 2*threshhold
            elseif mid + threshhold > length(df.f)
                start, stop = length(df.f) - 2*threshhold, length(df.f) 
            else
                start, stop = mid - threshhold, mid + threshhold
            end

            popt, perr, ci = bootstrap(gauss, df.f[start:stop], df.V[start:stop], p0=[1.35, idx_to_f(mid), 30, 0.06])
            temparr = [temparr; measurement(popt[2], popt[3])]
            if plt
                ax.plot(ci.x, gauss(ci.x, nom.(popt)), label="fit")
            end
        end
    end
    return temparr
end
    


intervals = [[Interval(3, 12), Interval(13,19), Interval(20, 24)],
            [Interval(3, 5), Interval(6, 12), Interval(13,21), Interval(22, 24)],
            [Interval(3, 5), Interval(6, 12), Interval(13,21), Interval(22, 24)]]

begin
lin(x, p) = @. p[1] * x + p[2]
function l(x, y1, y2)
    a = (y2 - y1) / (currentdf.I1[end] - currentdf.I1[1])
    b = y1 - a * currentdf.I1[1]
    return @. a * x + b
end

df, tempdf = DataFrame(), DataFrame()
Z = zeros(6000, length(currentdf.I1))
modifier = 0
y1, y2, shift = 0, 0, 0
ax = subplots(3, 1, figsize=(7, 4.5), sharex=true)[1]
for (j, Temp) in enumerate([20, 25, 30])
# for (j, Temp) in enumerate([20])
    maxima, maxima_I = [], []
    for i in 1:26
        println(i)
        tempdf = CSV.read(joinpath(@__DIR__, "data/T$Temp/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
        tempdf.I = currentdf.I1[i] * ones(length(tempdf.t))
        tempdf.f = t_to_f(tempdf.t)

        Z[:, i] = tempdf.V[end:-1:1] 

        # clip tempdf according to linear function
        if j == 1
            if i < 10
                y1, y2, shift = 600, -300, 500
            else
                y1, y2, shift = 1200, 300, 500
            end
        elseif j == 2
            y1, y2, shift = 800, 0, 800
        elseif j == 3
            y1, y2, shift = 300, 0, 600
        end
        tempdf.V[l(currentdf.I1[i], y1, y2) .> tempdf.f] *= 0 
        tempdf.V[l(currentdf.I1[i], y1, y2) .+ shift .< tempdf.f] *= 0
        # println(tempdf.f[1])
        # println(tempdf.f[end])
        if j == 1
            tempmax = get_maximum(tempdf, 0.925)
        elseif j == 2
            tempmax = get_maximum(tempdf, 0.91)
        elseif j == 3
            tempmax = get_maximum(tempdf, 0.85)
        end

        # if i != 1
        #     if tempmax - modifier > maxima[i-1] + 100
        #         modifier = modifier + FRS2
        #     end
        # end

        # append tempmax to maxima
        maxima = [maxima; tempmax .- modifier]
        maxima_I = [maxima_I; currentdf.I1[i] * ones(length(tempmax))]

        df = vcat(df, tempdf)
    end

    # shift maxima
    # maxima = maxima .- nom(maximum(maxima))

    # Z = Z[:, end:-1:1]

    ax[j-1].imshow(Z, aspect="auto", extent=(currentdf.I1[1], currentdf.I1[end], 0, 1270), cmap="seismic")
    ax[j-1].errorbar(maxima_I, nom.(maxima), yerr=err.(maxima), fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
    ax[j-1].plot(currentdf.I1, l(currentdf.I1, y1, y2), c="C0")
    ax[j-1].plot(currentdf.I1, l(currentdf.I1, y1, y2) .+ shift, c="C0")
    # for (i, interval) in enumerate(intervals[j])
    #     # popt, ci = bootstrap(lin, currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), p0=[1., 1.], unc=true, its=100000)
    #     # ax[j-1].plot(ci.x, lin(ci.x, nom.(popt)), c="C$i")
    #     ax[j-1].errorbar(currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), fmt="o", capsize=3, mfc="C$i", mec="k", ms=7, ecolor="k")
    #     # ax[j-1].fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C$i")

    #     # println(popt)
    # end
    ax[j-1].set_ylim(0, 1200)
end
end
imshow(Z, aspect="auto", extent=(currentdf.I1[1], currentdf.I1[end], tempdf.f[1], tempdf.f[end]), cmap="seismic")
# t_to_f(df.t[end])


# 2x2 grid of heatmaps
begin
df, tempdf = DataFrame(), DataFrame()
Z = zeros(6000, 26)
ax = subplots(2, 2, figsize=(10, 7.5), sharey=true)[1]
for (j, Temp) in enumerate([20, 25, 30])
    for i in 1:26
        tempdf = CSV.read(joinpath(@__DIR__, "data/T$Temp/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
        tempdf.f = t_to_f(tempdf.t)
        Z[:, i] = tempdf.V[end:-1:1] 
    end
    ax[round(Int, j/4), (j-1)%2].imshow(Z, aspect="auto", extent=(currentdf.I1[1], currentdf.I1[end], tempdf.f[1], tempdf.f[end]), cmap="seismic")
    ax[round(Int, j/4), (j-1)%2].set_title(latexstring("T = $Temp.0(5)\\ \\mathrm{^\\circ C}"))
    ax[round(Int, j/4), (j-1)%2].set_xlabel(L"I\ (\mathrm{mA})")

    ax[round(Int, j/4), (j-1)%2].tick_params(axis="both", direction="out", which="both", top=false, right=false)
end
Z = zeros(6000, 11)
for (i, Temp) in enumerate(20:30)
    println(i)
    tempdf = CSV.read(joinpath(@__DIR__, "data/Tempsweep/T$Temp.CSV"), DataFrame, header=["t", "V"], skipto=2)
    tempdf.f = t_to_f(tempdf.t)

    Z[:, i] = tempdf.V[end:-1:1] 
end
ax[1, 1].imshow(Z, aspect="auto", extent=(20, 30, tempdf.f[1], tempdf.f[end]), cmap="seismic")
ax[1, 1].set_title(L"I =  \mathrm{mA}")
ax[1, 1].set_xlabel(L"T\ (\mathrm{^\circ C})")
# change tick direction
ax[1, 1].tick_params(axis="both", direction="out", which="both", top=false, right=false)

ax[1, 0].set_ylabel(L"\Delta f\ (\mathrm{GHz})")
ax[0, 0].set_ylabel(L"\Delta f\ (\mathrm{GHz})")

tight_layout()
# savefig(string(@__DIR__, "/bilder/heatmaps.pdf"), bbox_inches="tight")
end


intervals = [[Interval(1, 6), Interval(7,14), Interval(15, 28), Interval(25, 29), Interval(30, 31)],
            [Interval(3, 6), Interval(7, 13), Interval(14,22), Interval(23, 25), Interval(26, 27)],
            []]

rest = [[([1, 4], "white"), ([19, 21, 23], "C4")],
        [([5], "white")],
        [([4, 13], "white")]]

begin
lin(x, p) = @. p[1] * x + p[2]
function l(x, y1, y2)
    a = (y2 - y1) / (currentdf.I1[end] - currentdf.I1[1])
    b = y1 - a * currentdf.I1[1]
    return @. a * x + b
end

df, tempdf = DataFrame(), DataFrame()
y1, y2, shift = 0, 0, 0
ax = subplots(1, 1, figsize=(7, 4.5), sharex=true)[1]
# for (j, Temp) in enumerate([20, 25, 30])
for (j, Temp) in enumerate([20])
    j = 1
    maxima, maxima_I = [], []
    for i in 1:26
        println(i)
        tempdf = CSV.read(joinpath(@__DIR__, "data/T$Temp/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
        tempdf.I = currentdf.I1[i] * ones(length(tempdf.t))
        tempdf.f = t_to_f(tempdf.t)

        # clip tempdf according to linear function
        if j == 1
            if i < 10
                y1, y2, shift = 600, -300, 500
            else
                y1, y2, shift = 1200, 300, 500
            end
        elseif j == 2
            y1, y2, shift = 800, 0, 800
        elseif j == 3
            y1, y2, shift = 300, 0, 600
        end
        tempdf.V[l(currentdf.I1[i], y1, y2) .> tempdf.f] *= 0 
        tempdf.V[l(currentdf.I1[i], y1, y2) .+ shift .< tempdf.f] *= 0

        if j == 1
            tempmax = get_maximum(tempdf, 0.925)
        elseif j == 2
            tempmax = get_maximum(tempdf, 0.91)
        elseif j == 3
            tempmax = get_maximum(tempdf, 0.85)
        end

        # append tempmax to maxima
        maxima = [maxima; tempmax]
        maxima_I = [maxima_I; currentdf.I1[i] * ones(length(tempmax))]

        df = vcat(df, tempdf)
    end
    if j == 1
        maxima[12:end] .-= FRS2
    end
    # shift maxima
    maxima = maxima .- nom(maxima[1])

    poptall, ciall = bootstrap(lin, maxima_I, nom.(maxima), yerr=err.(maxima), p0=[1., 1.], unc=true, its=10000)

    ax.errorbar(maxima_I, nom.(maxima), yerr=err.(maxima), fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
    ax.plot(maxima_I, lin(maxima_I, nom.(poptall)), c="gray", zorder=0, lw=2, alpha=0.5)
    if j != 3
        for (i, interval) in enumerate(intervals[j])

            popt, ci = bootstrap(lin, maxima_I[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), p0=[1., 1.], unc=true, its=100000)
            ax.plot(ci.x, lin(ci.x, nom.(popt)), c="C$i")
            ax.errorbar(maxima_I[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), fmt="o", capsize=3, mfc="C$i", mec="k", ms=7, ecolor="k")
            # ax[j-1].fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C$i")

            # println(popt)
        end
    elseif j == 3
        for (k, cut) in enumerate(cuts[2:end])
            tempmax, tempmax_I = maxima[maxima .> cut], maxima_I[maxima .> cut]
            tempmax, tempmax_I = tempmax[tempmax .< cuts[k]], tempmax_I[tempmax .< cuts[k]]
            
            popt, ci = bootstrap(lin, tempmax_I, nom.(tempmax), yerr=err.(tempmax), p0=[1., 1.], unc=true, its=10000)
            
            ax.plot(tempmax_I, lin(tempmax_I, nom.(popt)), c="C$k", zorder=0, lw=2, alpha=0.5)
            ax.errorbar(tempmax_I, nom.(tempmax), yerr=err.(tempmax), fmt="o", capsize=3, mfc="C$k", mec="k", ms=7, ecolor="k")
        end
    end
    for k in rest[j]
        ax.errorbar(maxima_I[k[1]], nom.(maxima[k[1]]), yerr=err.(maxima[k[1]]), fmt="o", capsize=3, mfc=k[2], mec="k", ms=7, ecolor="k")
    end
    # ax.set_ylim(0, 1200)
end
end

maxima[maxima .< 0]

cuts = [50, -20,  -100, -150, -200, -280, -400]

rest[1][1]
tempmax = maxima[maxima .< -280]
tempmax = tempmax[tempmax .> -400]

begin    
df, tempdf = DataFrame(), DataFrame()
Z = zeros(6000, 11)
# modifier = 0
# y1, y2, shift = 0, 0, 0
ax = subplots(figsize=(7, 4.5), sharex=true)[1]
# maxima, maxima_I = [], []
for (i, Temp) in enumerate(20:30)
    println(i)
    tempdf = CSV.read(joinpath(@__DIR__, "data/Tempsweep/T$Temp.CSV"), DataFrame, header=["t", "V"], skipto=2)
    tempdf.f = t_to_f(tempdf.t)

    Z[:, i] = tempdf.V[end:-1:1] 

    df = vcat(df, tempdf)
end

# shift maxima
# maxima = maxima .- nom(maximum(maxima))

# Z = Z[:, end:-1:1]

ax.imshow(Z, aspect="auto", extent=(20, 30, 0, 1270), cmap="seismic")
# ax.errorbar(maxima_I, nom.(maxima), yerr=err.(maxima), fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
# ax.plot(currentdf.I1, l(currentdf.I1, y1, y2), c="C0")
# ax.plot(currentdf.I1, l(currentdf.I1, y1, y2) .+ shift, c="C0")
# for (i, interval) in enumerate(intervals[j])
#     # popt, ci = bootstrap(lin, currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), p0=[1., 1.], unc=true, its=100000)
#     # ax[j-1].plot(ci.x, lin(ci.x, nom.(popt)), c="C$i")
#     ax[j-1].errorbar(currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), fmt="o", capsize=3, mfc="C$i", mec="k", ms=7, ecolor="k")
#     # ax[j-1].fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color="C$i")

#     # println(popt)
# end
# ax[j-1].set_ylim(0, 1200)
end



df


for i in 20:25
    println(i)
    df = CSV.read(joinpath(@__DIR__, "data/T30/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
    df.f = t_to_f(df.t)

    # df.V = savitzky_golay(df.V, 201, 2).y

    # normalize data
    offset = 0.01 / maximum(df.V)
    gauss(x, p) = @. p[1] * exp(-2 * log(2) * (x - p[2])^2 / p[3]^2) - offset

    df.V /= maximum(df.V)
    # println(mean(df.V))
    if mean(df.V) < 0.5
        # get two maxima
        max, height = ss.find_peaks(df.V, height=0.5, distance=f_to_idx(40))
        max = pyconvert(Array, max) .+ 1

        ax = subplots(figsize=(7, 4.5))[1]
        ax.plot(df.f, df.V, label="data")
        
        threshhold = f_to_idx(25)
        for mid in max
            # start, stop = f_to_idx(mid - 50), f_to_idx(mid + 50)
                
            if mid - threshhold < 1
                start, stop = 1, 2*threshhold
            elseif mid + threshhold > length(df.f)
                start, stop = length(df.f) - 2*threshhold, length(df.f) 
            else
                start, stop = mid - threshhold, mid + threshhold
            end

            popt, perr, ci = bootstrap(gauss, df.f[start:stop], df.V[start:stop], p0=[1.35, idx_to_f(mid), 30, 0.06])
            ax.plot(ci.x, gauss(ci.x, nom.(popt)), label="fit")
        end
        ax.scatter(df.f[max], df.V[max], label="maxima", color="C1")
    end
end

begin
df = CSV.read(joinpath(@__DIR__, "data/T25/Current8.CSV"), DataFrame, header=["t", "V"], skipto=2)
df.f = t_to_f(df.t)
df.V /= maximum(df.V)
get_maximum(df, 0.91, plt=true)
end
max, height = ss.find_peaks(df.V, height=0.6, distance=f_to_idx(40))
max = pyconvert(Array, max) .+ 1

begin
plot(df.f, df.V)
scatter(df.f[max], df.V[max])
hlines(0.6, df.f[1], df.f[end])
end

df = df[df.f .> l(df.f, 800, 400), :]

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

# voigt profile as convolution of gaussian and lorentzian
gauss(x, p) = @. p[1] * exp(-2 * log(2) * (x - p[2])^2 / p[3]^2) - 0.01


tempdf = CSV.read(joinpath(@__DIR__, "data/T30/Current24.CSV"), DataFrame, header=["t", "V"], skipto=2)
tempdf.f = t_to_f(tempdf.t)

tempdf.V[l(currentdf.I1[1], y1, y2) .> tempdf.f] *= 0

tempdf.V

# tempdf.I = currentdf.I1[i] * ones(length(tempdf.t))
tempdf.f = t_to_f(tempdf.t)

get_maximum(tempdf, plt=true)
tempdf.V = savitzky_golay(tempdf.V, 201, 2).y
