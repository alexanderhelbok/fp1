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

colors = ["C0", "C1", "C2", "C3", "C4", "C9"]
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
    # println("chisq =  $(chisq(nom.(tempdf.P)[starts[Temp]:end], lin(nom.(tempdf.I)[starts[Temp]:end], popt), pcount=2))")
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

FSR2 = 630

# find 2 maxima in the data
maxima, height = ss.find_peaks(df.V, height=0.5, distance=100)
# convert to julia array
maxima = pyconvert(Array, maxima)

dt = df.t[maxima[2]] - df.t[maxima[1]]

t_to_f(t) = @. FSR2 / dt * (t - t[1])

didx = maxima[2] - maxima[1]

idx_to_f(idx) = @. FSR2 / didx * (idx - 1)
f_to_idx(f) = @. didx / FSR2 * f + 1 |> round |> Int
# df.f = @. FSR2 / dt * (df.t - df.t[1])
end
FSR2/dt



currentdf = CSV.read(joinpath(@__DIR__, "data/Istepsize.csv"), DataFrame, header=["I1", "I2", "I3"], skipto=2)

cauchy(x, x0, g) = @. 1 / ( pi * g * ( 1 + ( ( x - x0 )/ g )^2 ) )
gauss( x, x0, s) = @. 1/ sqrt(2 * pi * s^2 ) * exp( - (x-x0)^2 / ( 2 * s^2 ) )

function pseudo_voigt(x, p)
    # fg = 2 * s * sqrt( 2 * log(2) )
    # fl = 2 * g
    fg = p[3]
    fl = p[4]
    f = abs( fg^5 +  2.69269 * fg^4 * fl + 2.42843 * fg^3 * fl^2 + 4.47163 * fg^2 * fl^3 + 0.07842 * fg * fl^4+ fl^5)^(1/5)
    eta = 1.36603 * ( fl / f ) - 0.47719 * ( fl / f )^2 + 0.11116 * ( f / fl )^3
    return @. p[1] * ( eta * cauchy( x, p[2], f) + ( 1 - eta ) * gauss( x, p[2], f ) )
end

function get_maximum(df, height; plt=false)
    # df.V = savitzky_golay(df.V, 201, 2).y
    temparr = []
    # normalize data
    offset = 0.01 / maximum(df.V)
    # gauss(x, p) = @. p[1] * exp(-4 * log(2) * (x - p[2])^2 / p[3]^2) - offset

    function pseudo_voigt(x, p)
        # fg = 2 * s * sqrt( 2 * log(2) )
        # fl = 2 * g
        fg = p[3]
        fl = p[4]
        f = abs( fg^5 +  2.69269 * fg^4 * fl + 2.42843 * fg^3 * fl^2 + 4.47163 * fg^2 * fl^3 + 0.07842 * fg * fl^4+ fl^5)^(1/5)
        eta = 1.36603 * ( fl / f ) - 0.47719 * ( fl / f )^2 + 0.11116 * ( f / fl )^3
        return @. p[1] * ( eta * cauchy( x, p[2], f) + ( 1 - eta ) * gauss( x, p[2], f ) ) - offset
    end

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

            # popt, perr, ci = bootstrap(gauss, df.f[start:stop], df.V[start:stop], p0=[1.35, idx_to_f(mid), 30, 0.06])
            popt, perr, ci = bootstrap(pseudo_voigt, df.f[start:stop], df.V[start:stop], p0=[10., idx_to_f(mid), 1., 2.], redraw=true)
            fV =  0.5346*popt[4] + sqrt(0.2166*popt[4]^2 + popt[3]^2)
            # temparr = [temparr; measurement(popt[2], popt[3])]
            temparr = [temparr; measurement(popt[2], fV)]
            if plt
                # ax.plot(ci.x, gauss(ci.x, nom.(popt)), label="fit")
                ax.plot(ci.x, pseudo_voigt(ci.x, nom.(popt)), label="fit")
            end
        end
    end
    return temparr
end
    
for i in 11:20
    tempdf = CSV.read(joinpath(@__DIR__, "data/T30/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
    tempdf.f = t_to_f(tempdf.t)

    get_maximum(tempdf, 0.85, plt=true)
    # title(i)
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
        #         modifier = modifier + FSR2
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
    #     # ax[j-1].plot(ci.x, lin(ci.x, nom.(popt)), c=colors[i])
    #     ax[j-1].errorbar(currentdf.I1[interval.first:interval.last], nom.(maxima[interval.first:interval.last]), yerr=err.(maxima[interval.first:interval.last]), fmt="o", capsize=3, mfc=colors[i], mec="k", ms=7, ecolor="k")
    #     # ax[j-1].fill_between(ci.x, ci.c0, ci.c1, alpha=0.3, color=colors[i])

    #     # println(popt)
    # end
    ax[j-1].set_ylim(0, 1200)
end
end


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
    ax[round(Int, j/4), (j-1)%2].set_title(latexstring("T = $Temp.0(2)\\ \\mathrm{^\\circ C}"))
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
ax[1, 1].set_title(L"I = 35.0(1)\ \mathrm{mA}")
ax[1, 1].set_xlabel(L"T\ (\mathrm{^\circ C})")
# change tick direction
ax[1, 1].tick_params(axis="both", direction="out", which="both", top=false, right=false)

ax[1, 0].set_ylabel(L"\Delta f\ (\mathrm{GHz})")
ax[0, 0].set_ylabel(L"\Delta f\ (\mathrm{GHz})")

tight_layout()
# savefig(string(@__DIR__, "/bilder/heatmaps.pdf"), bbox_inches="tight")
end

begin
modes = [[[2, 3, 5, 6], 
        [7, 8, 9, 10, 11, 12, 13, 14], 
        [15, 16, 17, 18, 20, 22, 24],
        [19, 21, 23, 25, 26, 27, 28, 29],
        [30, 31]],
        [[3, 4, 6],
        [7, 8, 9, 10, 11, 12, 13],
        [14, 15, 16, 17, 18, 19, 20, 21, 22],
        [23, 24, 25],
        [26, 27]],
        [[1, 2, 3, 5, 7, 9, 11, 14, 16],
        [6, 8, 10, 12, 15, 17, 18, 19, 20, 22, 24, 27, 30, 33, 36, 39],
        [21, 23, 26, 29, 32, 35, 38, 42, 45],
        [25, 28, 31, 34, 37, 41, 44, 47, 50, 53],
        [40, 43, 46, 49, 52],
        [48, 51]]]

lin(x, p) = @. p[1] * x + p[2]
function l(x, y1, y2)
    a = (y2 - y1) / (currentdf.I1[end] - currentdf.I1[1])
    b = y1 - a * currentdf.I1[1]
    return @. a * x + b
end

df, tempdf = DataFrame(), DataFrame()
y1, y2, shift = 0, 0, 0
ax = subplots(3, 1, figsize=(7, 9), sharex=true)[1]
for (j, Temp) in enumerate([20, 25, 30])
# for (j, Temp) in enumerate([20])
    # j = 3
    maxima, maxima_I, diffs = [], [],[]
    for i in 1:26
        println(i)
        tempdf = CSV.read(joinpath(@__DIR__, "data/T$Temp/Current$i.CSV"), DataFrame, header=["t", "V"], skipto=2)
        tempdf.I = currentdf.I1[i] * ones(length(tempdf.t))
        tempdf.f = t_to_f(tempdf.t)

        # clip tempdf according to linear function
        if j == 1
            if i < 6
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

        if j == 3
            diffs = [diffs; diff(tempmax)]
            println(diff(tempmax))
        end

        df = vcat(df, tempdf)
    end
    if j == 3
        println(diffs)
        println(mean(diffs))
    end
    if j == 1
        maxima[8:end] .-= FSR2
    end
    # shift maxima
    maxima = maxima .- nom(maxima[1])

    # calculate finesse
    Finesse = @. FSR2 / err(maxima)
    println(Finesse)

    poptall, ciall = bootstrap(lin, maxima_I, nom.(maxima), yerr=err.(maxima), p0=[1., 1.], unc=true, its=10000)
    println("all : $(poptall[1])")

    ax[j-1].errorbar(maxima_I, nom.(maxima), yerr=err.(maxima), fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
    xlims, ylims = ax[j-1].get_xlim(), ax[j-1].get_ylim()
    ax[j-1].plot(ciall.x, lin(ciall.x, nom.(poptall)), c="gray", zorder=0, lw=2, alpha=0.5, ls="--")
   
    for (k, mode) in enumerate(modes[j])
        # println(mode)
        tempmax, tempmax_I = maxima[mode], maxima_I[mode]
        popt, ci = bootstrap(lin, tempmax_I, nom.(tempmax), yerr=err.(tempmax), p0=[1., 1.], unc=true, its=10000, xlimconst=true, xlim=1)

        println("mode $k : $(popt[1])")

        ax[j-1].plot(ci.x, lin(ci.x, nom.(popt)), c=colors[k], zorder=0, lw=2, alpha=0.5)
        ax[j-1].errorbar(tempmax_I, nom.(tempmax), yerr=err.(tempmax), fmt="o", capsize=3, mfc=colors[k], mec="k", ms=7, ecolor="k")
    end
    ax[j-1].set_title(latexstring("T = $Temp.0(2)\\ \\mathrm{^\\circ C}"))
    ax[j-1].set_ylabel(L"\Delta f\ (\mathrm{GHz})")
    
    ax[j-1].set_xlim(xlims)
    ax[j-1].set_ylim(ylims)

    # ax.set_ylim(0, 1200)
end

# legend
errorbar(1, 1, xerr=1, yerr=1, fmt=".", capsize=3,  mfc="silver", mec="k", ms=11, ecolor="k", label=L"\textrm{Data}")
plot(1, 1, c="gray", label=L"\textrm{fit mode}")
plot(1, 1, c="gray", ls="--", label=L"\textrm{fit all}")

legax = ax[2].twinx()
legax.fill_between((1, 2), 1, 1, alpha=1, color="white", ec="k", label=latexstring(L"\textrm{other}"))
for (i, c) in enumerate(colors)
    legax.fill_between((1, 2), 1, 1, alpha=1, color=c, ec="k", label=latexstring("\\textrm{mode}\\ $i"))
end

legax.set_axis_off()
ax[2].tick_params(right=true, which="both")

# ax.legend()
legax.legend(loc="lower center", ncols=2, bbox_to_anchor =(0.63, -1.3))
ax[2].legend(bbox_to_anchor =(0.23, -1.15), loc="lower center")

ax[2].set_xlabel(L"I\ (\mathrm{mA})")
tight_layout()
subplots_adjust(hspace=0.3)
# savefig(string(@__DIR__, "/bilder/fitI.pdf"), bbox_inches="tight")
end



# voigt profile as convolution of gaussian and lorentzian
gauss(x, p) = @. p[1] * exp(-2 * log(2) * (x - p[2])^2 / p[3]^2) - 0.01


tempdf = CSV.read(joinpath(@__DIR__, "data/Tempsweep/T20.CSV"), DataFrame, header=["t", "V"], skipto=2)
tempdf.f = t_to_f(tempdf.t)

tempdf.V

# tempdf.I = currentdf.I1[i] * ones(length(tempdf.t))
tempdf.f = t_to_f(tempdf.t)

get_maximum(tempdf, 0.8, plt=true)
tempdf.V = savitzky_golay(tempdf.V, 201, 2).y

begin
modes = [[4, 5, 6], [7, 8], [9, 12], [10, 13], [11, 14]]
function l(x, y1, y2)
    a = (y2 - y1) / (Temps[end] - Temps[1])
    b = y1 - a * Temps[1]
    return @. a * x + b
end
    
df, tempdf = DataFrame(), DataFrame()
Temps = [20:1:30;]
Z = zeros(6000, length(Temps))
y1, y2, shift = 0, 0, 0
maxima, maxima_I = [], []
for (i, Temp) in enumerate(Temps)
    println(i)
    tempdf = CSV.read(joinpath(@__DIR__, "data/Tempsweep/T$Temp.CSV"), DataFrame, header=["t", "V"], skipto=2)
    tempdf.f = t_to_f(tempdf.t)

    Z[:, i] = tempdf.V[end:-1:1] 

    # clip tempdf according to linear function
    y1, y2, shift = 600, -150, 650
    
    tempdf.V[l(Temps[i], y1, y2) .> tempdf.f] *= 0 
    tempdf.V[l(Temps[i], y1, y2) .+ shift .< tempdf.f] *= 0

    tempmax = get_maximum(tempdf, 0.7)

    # if i != 1
    #     if tempmax - modifier > maxima[i-1] + 100
    #         modifier = modifier + FSR2
    #     end
    # end

    # append tempmax to maxima
    maxima = [maxima; tempmax]
    maxima_I = [maxima_I; Temps[i] * ones(length(tempmax))]

    df = vcat(df, tempdf)
end
# shift maxima
maxima .-= nom(maxima[1])
end

begin
ax = subplots(figsize=(7, 3.5), sharex=true)[1]

poptall, ciall = bootstrap(lin, maxima_I, nom.(maxima), yerr=err.(maxima), p0=[1., 1.], unc=true, its=10000)
println("all : $(poptall[1])")

ax.plot(ciall.x, lin(ciall.x, nom.(poptall)), c="gray", zorder=0, lw=2, alpha=0.5, ls="--")

# ax.imshow(Z, aspect="auto", extent=(x[1], x[end], 0, 1270), cmap="seismic")
ax.errorbar(maxima_I, nom.(maxima), yerr=err.(maxima), fmt="o", capsize=3, mfc="white", mec="k", ms=7, ecolor="k")
# ax.plot(Temps, l(Temps, y1, y2), c="C0")
# ax.plot(Temps, l(Temps, y1, y2) .+ shift, c="C0")

for (k, mode) in enumerate(modes)
    # println(mode)
    tempmax, tempmax_I = maxima[mode], maxima_I[mode]
    popt, ci = bootstrap(lin, tempmax_I, nom.(tempmax), yerr=err.(tempmax), p0=[1., -300.], unc=true, xlim=0.5, xlimconst=true)
    println("mode $k : $((popt[1]))")
    
    ax.plot(ci.x, lin(ci.x, nom.(popt)), c=colors[k], zorder=0, lw=2, alpha=0.5)
    ax.errorbar(tempmax_I, nom.(tempmax), yerr=err.(tempmax), fmt="o", capsize=3, mfc=colors[k], mec="k", ms=7, ecolor="k")
end

xlims, ylims = ax.get_xlim(), ax.get_ylim()

# legend
errorbar(1, 1, xerr=1, yerr=1, fmt=".", capsize=3,  mfc="silver", mec="k", ms=11, ecolor="k", label=L"\textrm{Data}")
plot(1, 1, c="gray", label=L"\textrm{fit mode}")
plot(1, 1, c="gray", ls="--", label=L"\textrm{fit all}")

legax = ax.twinx()
legax.fill_between((1, 2), 1, 1, alpha=1, color="white", ec="k", label=latexstring(L"\textrm{other}"))
for (i, c) in enumerate(colors[1:end-1])
    legax.fill_between((1, 2), 1, 1, alpha=1, color=c, ec="k", label=latexstring("\\textrm{mode}\\ $i"))
end

legax.set_axis_off()
ax.tick_params(right=true, which="both")

ax.set_xlabel(L"T\ (\mathrm{^\circ C})")
ax.set_ylabel(L"\Delta f\ (\mathrm{GHz})") 

ax.set_xlim(xlims)
ax.set_ylim(ylims)

ax.legend()
legax.legend(loc=(0.02, 0.05))

tight_layout()
# savefig(string(@__DIR__, "/bilder/fitTemp.pdf"), bbox_inches="tight")
end

mean(diffs[5:end])
