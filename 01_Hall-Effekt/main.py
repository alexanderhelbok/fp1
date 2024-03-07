import matplotlib.pyplot as plt

from Source import *
from scipy.optimize import fsolve
import astropy.constants as const


# load data
tempdf = pd.read_csv("01_Hall-Effekt/data.csv", delimiter=",")
# print(os.listdir())
# print(df)

tempdf.columns = ["T1", "1a", "1b", "2a", "2b", "T2", "3a", "3b", "3aB", "3bB"]

# make all data to ufloat with 0.1% error
# for col in tempdf.columns[:]:
#     tempdf[col] = tempdf[col].apply(lambda x: unc.ufloat(x, 0.1))
    # tempdf[col] = unp.uarray(tempdf[col], 0.1*np.ones(len(tempdf[col])))

def T(R):
    return 509/8 - 98/(2*R) + 43/30 * np.sqrt(R) + 161/113 * R**(11/10)

# create empty dataframe with dtype ufloat
df = pd.DataFrame()

# interpolate Temperature
df["T1"] = (tempdf.loc[:62, "T1"].values + tempdf.loc[:62, "T2"].values) / 2
df["T2"] = (tempdf.loc[1:, "T1"].values + tempdf.loc[:62, "T2"].values) / 2
# drop last 5 rows
df["U1"] = (tempdf.loc[:, "1a"] - tempdf.loc[:, "1b"]) / 2
df["U2"] = (tempdf.loc[:, "2a"] - tempdf.loc[:, "2b"]) / 2
df["U3"] = (tempdf.loc[:, "3a"] - tempdf.loc[:, "3b"]) / 2
df["U3B"] = (tempdf.loc[:, "3aB"] - tempdf.loc[:, "3bB"]) / 2
df = df[4:-1]
print(df)
# %%
# U1, U2, U3, U3B = unp.uarray(df["U1"], np.ones(62)/np.sqrt(2)), unp.uarray(df["U2"], np.ones(62)/np.sqrt(2)), unp.uarray(df["U3"], np.ones(62)/np.sqrt(2)), unp.uarray(df["U3B"], np.ones(62/np.sqrt(2)))
U1 = unp.uarray(df["U1"], np.ones(len(df))*0.035)*1e-3
U2 = unp.uarray(df["U2"], np.ones(len(df))*0.035)*1e-3
U3 = unp.uarray(df["U3"], np.ones(len(df))*0.035)*1e-3
U3B = unp.uarray(df["U3B"], np.ones(len(df))*0.035)*1e-3
T1 = unp.uarray(df["T1"], np.ones(len(df))*0.035)
T2 = unp.uarray(df["T2"], np.ones(len(df))*0.035)
# define constants
d = unc.ufloat(500*1e-6, 1e-6)
B = unc.ufloat(252.6 * 1e-3, 5*1e-5)
B0 = unc.ufloat(1.045*1e-3, 5*1e-5)
I = unc.ufloat(500 * 1e-6, 1e-8)

# calculate Hall coefficient
UH = U3B - U3
RH = UH * d / (B * I)

print(RH)
# %%
# plot data against T1
# %%
# solve transcendental equation for variable x using fsolve
def f(x):
    return fsolve(lambda unbekannt: np.cosh(np.log(2) / unbekannt * (x - 1)/(x+1)) - 0.5 * np.exp(np.log(2)/unbekannt), 1)


x = unp.nominal_values(1/unp.fabs(U1/U2))
temp = [f(X) for X in x]
F = np.ones(len(temp))
for i in range(len(temp)):
    F[i] = temp[i][0]

# calculate resistivity
rho = np.pi*d/np.log(2) * F * (unp.fabs(U1) + unp.fabs(U2)) / (2 * I)
print(rho)


# %%
RH = -RH
mu = RH/rho
# print(mu)

n = unp.log(1/(const.e.value * RH))

newT = (T1 + T2)/2

def func(x, a, b):
    return a * x + b


mean = np.mean(unp.nominal_values(n[:7]))

popt, pcov = curve_fit(func, 1/T(unp.nominal_values(newT[:-4])), unp.nominal_values(n[:-4]), sigma=unp.std_devs(n[:-4]))

Ed = -2*const.k_B.si.value * popt[0]
print(Ed, popt[0])

# print(n)
fig, ax = plt.subplots()
ax.errorbar(1/T(unp.nominal_values(newT[:-4])), unp.nominal_values(n[:-4]), yerr=unp.std_devs(n[:-4]), label="mu", fmt=".k", capsize=5)
ax.errorbar(1/T(unp.nominal_values(newT[-4:])), unp.nominal_values(n[-4:]), yerr=unp.std_devs(n[-4:]), label="mu", fmt=".r", capsize=5)
ax.plot(1/T(unp.nominal_values(newT[:-4])), func(1/T(unp.nominal_values(newT[:-4])), *popt), label="fit", zorder=5)
ax.axhline(mean, color="r", linestyle="--")
# plt.yscale("log")
plt.ylabel("log n [1/m**3]")
plt.xlabel("1/T [1/K]")
plt.show()
# %%


# plot data against newT
fig, ax = plt.subplots()
ax.errorbar(T(unp.nominal_values(newT)), unp.nominal_values(mu), yerr=unp.std_devs(mu), label="mu", fmt=".k", capsize=5)

# plt.yscale("log")
plt.ylabel("mu [m**2/Vs]")
plt.xlabel("T [K]")
plt.show()
# %%
newT = newT[:-2]
mu = mu[:-2]
mu = unp.log(mu)
newT = np.log(T(unp.nominal_values(newT)))
# %%


popt2, pcov2 = curve_fit(func, unp.nominal_values(newT[-5:]), unp.nominal_values(mu[-5:]), sigma=unp.std_devs(mu[-5:]))
popt3, pcov3 = curve_fit(func, unp.nominal_values(newT[:25]), unp.nominal_values(mu[:25]), sigma=unp.std_devs(mu[:25]))

# plot log mu against log T
fig, ax = plt.subplots()
ax.errorbar(unp.nominal_values(newT[10:-10]), unp.nominal_values(mu[10:-10]), yerr=unp.std_devs(mu[10:-10]), label="mu", fmt=".k", capsize=5)
ax.errorbar(unp.nominal_values(newT[:35]), unp.nominal_values(mu[:35]), yerr=unp.std_devs(mu[:35]), label="mu", fmt=".g", capsize=5)
ax.errorbar(unp.nominal_values(newT[-10:]), unp.nominal_values(mu[-10:]), yerr=unp.std_devs(mu[-10:]), label="mu", fmt=".r", capsize=5)

ax.plot(unp.nominal_values(newT[-10:]), func(unp.nominal_values(newT[-10:]), *popt2), label="fit", zorder=5)
ax.plot(unp.nominal_values(newT[:35]), func(unp.nominal_values(newT[:35]), *popt3), label="fit", zorder=5)
# ax.set_yscale("log")
# ax.set_xscale("log")

print(popt2, popt3)

plt.ylabel("log mu [m**2/Vs]")
plt.xlabel("log T [K]")
plt.show()
# %%
# plot data against T1
plt.plot(T(df["T1"]), df["1a"], label="1a")
plt.plot(T(df["T1"]), df["1b"], label="1b")
plt.plot(T(df["T1"]), df["2a"], label="2a")
plt.plot(T(df["T1"]), df["2b"], label="2b")
plt.plot(T(df["T1"]), df["3a"], label="3a")
plt.plot(T(df["T1"]), df["3b"], label="3b")
plt.plot(T(df["T1"]), df["3aB"], label="3aB")
plt.plot(T(df["T1"]), df["3bB"], label="3bB")
plt.plot(T(df["T1"]), df["T2"], label="T2")


# plt.xticks(np.arange(0, 100, 5))
# plt.yticks(np.arange(-50, 50, 5))
plt.xlabel("T1 [K]")
plt.ylabel("U [V]")
plt.legend()
plt.show()






# %%
print(np.log(T(unp.nominal_values(newT))))