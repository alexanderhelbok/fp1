using DataFrames, Random, LsqFit, Statistics, Measurements, PythonCall, Distributions
@py import scipy.signal as ss


function chisq(obs, exp; sigma=nothing, dof=nothing, pcount=nothing)
    """
    Calculate chi-squared value for a given set of observed and expected values.
    Parameters:
    obs: observed values
    exp: expected values
    sigma: optional error on observed values (default nothing)
    dof: degrees of freedom (default 0)
    pcount: number of parameters (default 0)
    Returns:
    chi-squared value
    """
    if dof === nothing && pcount === nothing
        dof = 1
    elseif dof === nothing
        dof = length(obs) - pcount
    end

    if sigma === nothing
        return sum((obs - exp).^2)/dof
    else
        return sum(((obs - exp)./sigma).^2)/dof
    end
end


function de(fit, xdata, ydata, bounds; mut=0.8, crossp=0.7, popsize=20, its=1000, fobj=chisq, seed=nothing, sigma=nothing)
    """
    Differential evolution algorithm to fit a function fobj(xdata, p...) with length(p) parameters.
    Parameters:
    fit: LsqFit object
    xdata: x data to be fitted
    ydata: y data to be fitted
    bounds: bounds for each parameter
    mut: mutation factor (default 0.8)
    crossp: crossover probability (default 0.7)
    popsize: population size (default 20)
    its: number of iterations (default 1000)
    fobj: objective function (default chisq)
    seed: random seed (default None)
    sigma: optional error on observed values (default None)
    Returns:
    optimal parameters, 1 sigma error of parameters
    """
    if seed !== nothing
        Random.seed!(seed)
    end

    dimensions = length(bounds)
    # create population with random parameters (between 0 and 1)
    pop = rand(popsize, dimensions)
    # scale parameters to the given bounds
    bounds = permutedims(hcat(bounds...))
    min_b, max_b = bounds[:, 1], bounds[:, 2]
    diff = diff = abs.(max_b - min_b)
    pop_denorm = min_b' .+ pop .* diff'


    fitness = [fobj(ydata, fit(xdata, p), sigma=sigma) for p in eachrow(pop_denorm)]
    best_idx = argmin(fitness) 
    best = pop_denorm[best_idx, :]

    for i in 1:its
        for j in 1:popsize
            idxs = filter(x -> x != j, 1:dimensions)
            mu = @view pop_denorm[idxs[sample(1:dimensions-1, 3, replace=false)], :]
            mutant = mu[1, :] .+ mut .* (mu[2, :] .- mu[3, :])
            cross_points = rand(dimensions) .< crossp
            
            trial = [if cross_points[k] 
                        mutant[k] 
                    else 
                        pop_denorm[j, k] 
                    end for k in 1:dimensions]

            trial = clamp.(trial, bounds[:, 1], bounds[:, 2])
            
            f = fobj(ydata, fit(xdata, trial), sigma=sigma)
            if f < fitness[j]
                fitness[j] = f
                pop[j, :] = trial
                if f < fitness[best_idx]
                    best_idx = j
                    best = trial
                end
            end
        end
        yield(best, fitness[best_idx])
    end
end

function bootstrap(fobj, xdata, ydata; xerr=zeros(length(xdata)), yerr=zeros(length(ydata)), 
    p=0.95, its=1000, p0=nothing, smoothing=false, unc=false)
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

    # initialize array for parameters and interpolated values for each iteration
    arr2 = Matrix{Float64}(undef, 1000, its)
    var = zeros(length(p0))
    sum = zeros(length(p0))
    # initialize DataFrame for confidence band
    ci = DataFrame(x=range(minimum(xdata), maximum(xdata), length=1000), c0=Vector(undef, 1000), c1=Vector(undef, 1000), mean=Vector(undef, 1000))
    for i in 1:its
        ind = rand(1:length(xdata), length(xdata))
        newx, newy = rand.(Normal.(xdata[ind], xerr[ind])), rand.(Normal.(ydata[ind], yerr[ind]))
        popt = curve_fit(fobj, newx, newy, p0).param
        arr2[:, i] = fobj(ci.x, popt)

        sum += popt
        var += popt.^2
        # scatter!(newx, newy, alpha=0.1, c=:gray, label=nothing)
    end

    pmean = sum/its
    perr = sqrt.(its/(its - 1)*(var/its - pmean.^2))

    ci.c0 = quantile.(eachrow(arr2), 0.5 * (1 - p))
    ci.c1 = quantile.(eachrow(arr2), 1 - 0.5 * (1 - p))
    ci.mean = mean(arr2, dims=2)[:, 1]

    popt = curve_fit(fobj, ci.x, ci.mean, pmean).param
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

