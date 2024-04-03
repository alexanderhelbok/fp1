using CSV, PythonPlot, SavitzkyGolay, Intervals, LaTeXStrings

include(string(@__DIR__, "/../Source.jl"))


cauchy(x, x0, g) = 1 / ( pi * g * ( 1 + ( ( x - x0 )/ g )^2 ) )
gauss( x, x0, s) = 1/ sqrt(2 * pi * s^2 ) * exp( - (x-x0)^2 / ( 2 * s^2 ) )

function pseudo_voigt(x, p)
    # fg = 2 * s * sqrt( 2 * log(2) )
    # fl = 2 * g
    fg = p[3]
    fl = p[4]
    f = ( fg^5 +  2.69269 * fg^4 * fl + 2.42843 * fg^3 * fl^2 + 4.47163 * fg^2 * fl^3 + 0.07842 * fg * fl^4+ fl^5)^(1/5)
    eta = 1.36603 * ( fl / f ) - 0.47719 * ( fl / f )^2 + 0.11116 * ( f / fl )^3
    return p[1] * ( eta * cauchy( x, p[2], f) + ( 1 - eta ) * gauss( x, p[2], f ) )
end

function bootstrap(fobj, xdata, ydata; xerr=zeros(length(xdata)), yerr=zeros(length(ydata)), 
    p=0.95, its=1000, samples=nothing, p0=nothing, smoothing=false, unc=false, xlim=0.1, xlimconst=false, redraw=true)
    """
    Bootstrap fit including confidence bands to a function fobj(xdata, p...) with length(p) parameters.
    Parameters:
    fobj: function to be fitted of the form fobj(xdata, p...)
    xdata: x data to be fitted
    ydata: y data to be fitted
    xerr: optional x error (default zeros(length(xdata)))
    yerr: optional y error (default zeros(length(ydata)))
    p: confidence interval (default 0.95)
    its: number of iterations (default 1000)
    p0: initial parameters; if none given, using 1 as starting value (default nothing)
    smoothing: whether to smooth the confidence band using Savitzky-Golay (default false)
    Returns:
    optimal parameters, 1 sigma error of parameters, x and y values for confidence band
    """
    if length(xdata) != length(ydata)
        throw(ArgumentError("x and y must be of the same size"))
    end
    if p < 0
        throw(ArgumentError("p must be positive"))
    elseif p > 1
        warn("p > 1, assuming p is a percentage")
        p = p / 100
    end
    if p0 === nothing
        throw(ArgumentError("p0 must be given (either as a vector or as a number)"))
    elseif typeof(p0) == Int
        p0 = ones(p0)
    end
    if samples === nothing
        samples = length(xdata)
    elseif typeof(samples) <: Number
        samples = round(Int, samples)
    end

    # initialize array for parameters and interpolated values for each iteration
    arr2 = Matrix{Float64}(undef, 1000, its)
    var = zeros(length(p0))
    sum = zeros(length(p0))
    # initialize DataFrame for confidence band
    if 0 <= xlim < 2 && xlimconst == false 
        span = maximum(xdata) - minimum(xdata)
        ci = DataFrame(x=range(minimum(xdata) - span*xlim, maximum(xdata) + span*xlim, length=1000), c0=Vector(undef, 1000), c1=Vector(undef, 1000), mean=Vector(undef, 1000))
    elseif xlimconst == true
        ci = DataFrame(x=range(minimum(xdata) - xlim, maximum(xdata) + xlim, length=1000), c0=Vector(undef, 1000), c1=Vector(undef, 1000), mean=Vector(undef, 1000))
    else
        ci = DataFrame(x=range(xlim[1], xlim[2], length=1000), c0=Vector(undef, 1000), c1=Vector(undef, 1000), mean=Vector(undef, 1000))
    end
    for i in 1:its
        if redraw 
            ind = rand(1:length(xdata), samples)
        else
            ind = [1:length(xdata);]
        end
        newx, newy = rand.(Normal.(xdata[ind], xerr[ind])), rand.(Normal.(ydata[ind], yerr[ind]))
        ax.scatter(newx, newy, alpha=0.1, c="gray", label=i)
        println(newx, newy)
        popt = curve_fit(fobj, newx, newy, p0).param
        arr2[:, i] = fobj(ci.x, popt)

        println(popt)

        sum += popt
        var += popt.^2
        plot(newx, fobj(newx, popt), alpha=0.1, c="gray", label=nothing)
        # scatter!(newx, newy, alpha=0.1, c=:gray, label=nothing)
    end

    pmean = sum/its
    popt = pmean
    perr = sqrt.(abs.(its/(its - 1)*(var/its - pmean.^2)))

    ci.c0 = quantile.(eachrow(arr2), 0.5 * (1 - p))
    ci.c1 = quantile.(eachrow(arr2), 1 - 0.5 * (1 - p))
    ci.mean = mean(arr2, dims=2)[:, 1]

    popt = curve_fit(fobj, ci.x, ci.mean, pmean).param
    # println(popt, pmean)
    # scatter(ci.x, ci.c0)
    # scatter(ci.x, ci.c1)
    # plot(ci.x, fobj(ci.x, popt), c="black", label="Fit")
    if smoothing
        # smooth the confidence band
        ci.c0 = ss.savgol_filter(ci.c0, 75, 1)
        ci.c1 = ss.savgol_filter(ci.c1, 75, 1)
    end
    if unc
        popt = measurement.(popt, perr)
        return popt, ci
    else
        return popt, perr, ci
    end
end

lin(x, p) = @. p[1] * x + p[2]

for i in 1:10
    popt = curve_fit(lin, x, y, [1., 1.]).param
    println(popt)
end

x = [1., 2.]
y = [1., 6.]

ax = subplots()[1]

popt, perr, ci = bootstrap(lin, x, y, p0=[5., 0.], its=10, redraw=false)

legend()
fV =  0.5346*popt[4] + sqrt(0.2166*popt[4]^2 + popt[3]^2)

scatter(x, y, label="Data")
plot(x, pseudo_voigt(x, [50, 1, 1, 0])) 
plot(ci.x, lin(ci.x, popt), label="Fit")
# draw  fwhm
axhline(maximum(y)/2, color="black", linestyle="--", label="FWHM")

pseudo_voigt(x, [1, 1, 1, 0])

df = CSV.read(joinpath(@__DIR__, "data/T25/Current15.CSV"), DataFrame, header=["t", "V"], skipto=2)




a = [1]

diff(a)