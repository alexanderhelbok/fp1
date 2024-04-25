using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings
using PhysicalConstants.CODATA2018: c_0, h, e

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

# define Constants
begin
d = 201.4e-12
U = 35 * u"kV"
end

theta_to_d(theta) = @. 2*d*sin(theta*π/180)

# load data
begin
df = CSV.read(string(@__DIR__, "/data/Ex1.csv"), DataFrame, header=["V", "N1", "N2", "N3"], skipto=2)
df.N1 = measurement.(df.N1, sqrt.(df.N1))
df.N2 = measurement.(df.N2, sqrt.(df.N2))
df.N3 = measurement.(df.N3, sqrt.(df.N3))
df.N = mean([df.N1, df.N2, df.N3])
end

# calculate weighted mean of plateau region
start = 14

mean_N = mean(df.N[start:end])

# plot data
begin
ax = subplots()[1]
myerrorbar(df.V, df.N, fmt="o", label="Messwerte")   
scatter(df.V[16:end], df.N[16:end], color="black", s=30, zorder=10)

plot(df.V[start:end], mean_N * ones(length(df.V[start:end])), label="Mittelwert", color="red")

xlabel(L"U (\mathrm{V})")
ylabel(L"N (\mathrm{1/s})")
end

# load data
begin
Voltages = [15:2.5:30;]
λ0 = []
for (i, V) in enumerate(Voltages)
    println(i)
    df = CSV.read(string(@__DIR__, "/data/Ex3_$i.csv"), DataFrame, header=["ang", "N"], skipto=3)
    df.N = measurement.(df.N, sqrt.(df.N))
    df.lam = theta_to_d(df.ang)*1e12

    mean_init = mean(df.N[1:1])
    idx_init = findfirst(df.N .> mean_init + 2*nom(mean(df.N[1:4]))) - 1
    # calculate mean of first datapoints
    mean_N = mean(df.N[1:idx_init])
    # set error to std
    mean_N = measurement(nom(mean_N), nom.(std(df.N[1:idx_init])))

    # find first datapoint where N is above 3 sigma of the mean
    idx = findfirst(df.N .> mean_N + 3*err(mean_N))
    # # interpolate to find the exact value
    λ0 = [λ0; measurement((df.lam[idx-1] + df.lam[idx])/2, (df.lam[idx] - df.lam[idx-1])/2)]

    myerrorbar(df.lam, df.N, fmt="o", label="Messwerte")
    # scatter(df.lam[1:idx_init], nom.(df.N[1:idx_init]), color="black", s=30, zorder=10)
    # plot(df.lam, nom.(mean_N) * ones(length(df.lam)), label="Mittelwert", color="red")
    # fill_between(df.lam, nom.(mean_N) - 3*err.(mean_N), nom.(mean_N) + 3*err.(mean_N), color="red", alpha=0.3)
end
end
λ0
model(x, p) = @. p[1]/x + p[2]

popt, ci = bootstrap(model, Voltages, nom.(λ0), yerr=err.(λ0), p0=[1000., 0.1], redraw=false, unc=true)

h_exp = popt[1]*u"kV*pm"*e/c_0 |> u"J*s"

# plot 
begin
ax = subplots()[1]
myerrorbar(Voltages, λ0, fmt="o", label="Messwerte")
plot(ci.x, model(ci.x, nom.(popt)), label="Fit", color="red")

xlabel(L"U (\mathrm{kV})")
ylabel(L"\lambda_0 (\mathrm{pm})")
end

h_exp = λ*u"pm"*e*U / (c_0)



# load new data
begin
df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
dfAl = CSV.read(string(@__DIR__, "/data/Ex4_1Al.csv"), DataFrame, header=["V", "ang", "N1", "N2", "N3"], skipto=2)
dfZn = CSV.read(string(@__DIR__, "/data/Ex4_1Zn.csv"), DataFrame, header=["V", "ang", "N1", "N2", "N3"], skipto=2)
for df in [dfAl, dfZn]
    df.N1 = measurement.(df.N1, sqrt.(df.N1))
    df.N2 = measurement.(df.N2, sqrt.(df.N2))
    df.N3 = measurement.(df.N3, sqrt.(df.N3))
    df.N = mean([df.N1, df.N2, df.N3])
    df.lam = theta_to_d(df.ang)
end
end

df

df.N[findfirst(df.ang .== unique(dfAl.ang)[1])]
reference = [df.N[findfirst(df.ang .== unique(dfAl.ang)[1])], df.N[findfirst(df.ang .== unique(dfAl.ang)[2])]]

model(x, p) = @. p[1]*exp(-p[2]*(x - p[3])) + p[4]

# plot data
begin
ax = subplots()[1]

for df in [dfAl, dfZn]
    for (i, lam) in enumerate(unique(df.lam))
        mask = df.lam .== lam
        # add reference to data
        x, y = [0.; df.V[mask]], [measurement(1., 0.); df.N[mask]]
        # x, y = df.V[mask], df.N[mask]

        # popt, ci = bootstrap(model, x, nom.(y), yerr=err.(y), p0=[0.0005, 1., 0.5, 10.], unc=true, redraw=false)

        println(popt)
        myerrorbar(df.V[mask], df.N[mask] ./ reference[i], fmt="o")   
        # plot(ci.x, model(ci.x, nom.(popt)) ./ reference[i], color="red")
        # scatter(df.V[mask], df.N[mask], color="black", s=30, zorder=10)
    end
end
myerrorbar(0, 1, fmt="o", label="Aluminium")

xlabel(L"U (\mathrm{V})")
ylabel(L"N (\mathrm{1/s})")

legend()
tight_layout()

end

thicknessdict = Dict("Al" => 0.04, "Zn" => 0.05, "Cu" => 0.025, "Sn" => 0.025, "Ni" => 0.025)
densitydict = Dict("Al" => 2.7, "Zn" => 7.14, "Cu" => 8.96, "Sn" => 7.3, "Ni" => 8.9)
Zdict = Dict("Al" => 13, "Zn" => 30, "Cu" => 29, "Sn" => 50, "Ni" => 28)

# load data
begin
I0df = CSV.read(string(@__DIR__, "/data/Ex2.csv"), DataFrame, header=["ang", "N"], skipto=3)
I0df.N = measurement.(I0df.N, sqrt.(I0df.N))
I0df.lam = theta_to_d(I0df.ang)
df = CSV.read(string(@__DIR__, "/data/Ex4_Al2.csv"), DataFrame, header=["ang", "N"], skipto=3)
df.N = measurement.(df.N, sqrt.(df.N))
df.lam = theta_to_d(df.ang)
I0df = filter(row -> row.ang in unique(df.ang), I0df)
end

μ = @. -log(df.N / I0df.N)/thicknessdict["Cu"]

# plot data
begin
ax = subplots()[1]
myerrorbar(df.lam[1:end-10], (μ./maximum(nom.(μ)))[1:end-10].^(1/3), fmt="o", label="Messwerte")
scatter(df.lam, nom.(df.N)./maximum(nom.(df.N)), color="black", s=30, zorder=10)
scatter(I0df.lam, nom.(I0df.N)./maximum(nom.(I0df.N)), color="red", s=30, zorder=10)
# scatter(df.lam, nom.(df.N), color="black", s=30, zorder=10)
# scatter(I0df.lam, nom.(I0df.N), color="red", s=30, zorder=10)
end