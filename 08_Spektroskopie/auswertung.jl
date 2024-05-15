using CSV, PythonPlot, SavitzkyGolay, LaTeXStrings
using PhysicalConstants.CODATA2018: c_0

include(string(@__DIR__, "/../Source.jl"))
mpl.use("pgf")
mpl.use("TkAgg")

# load data delimeter is tab
data = CSV.read(string(@__DIR__, "/data/data.csv"), DataFrame, skipto = 2, header = ["lam", "chi", "gf", "fit_lam", "Wf"])


data
# filter out rows containing NaN
data = data[.!isnan.(data.fit_lam), :]


isnan(data.fit_lam[2])