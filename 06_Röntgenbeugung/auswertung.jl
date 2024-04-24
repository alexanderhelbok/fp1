using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings
using PhysicalConstants.CODATA2018: c_0, h

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

# define Constants
begin
d = 201.4e-12

end


# load data
begin
df = CSV.read(string(@__DIR__, "/data/Ex1.csv"), DataFrame, header=["V", "N1", "N2", "N3"], skipto=2)
df.N1 = measurement.(df.N1, sqrt.(df.N1))
df.N2 = measurement.(df.N2, sqrt.(df.N2))
df.N3 = measurement.(df.N3, sqrt.(df.N3))
df.N = mean([df.N1, df.N2, df.N3])
end

# plot data
begin
ax = subplots()[1]
myerrorbar(df.V, df.N, fmt="o", label="Messwerte")   

xlabel(L"U (\mathrm{V})")
ylabel(L"N (\mathrm{1/s})")
end