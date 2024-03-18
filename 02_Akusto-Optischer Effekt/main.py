from Source import *
# %%

# load data
tempdf = pd.read_csv("/home/taco/Documents/fp1/02_Akusto-Optischer Effekt/data/data3.csv", delimiter=",", skiprows=1)

# tempdf.columns = ["mVPP", "0.Ordnung", "1.Ordnung"]
tempdf.columns = ["mVPP", "0.Ordnung"]

# plot
# plt.plot(tempdf["mVPP"]/0.006, tempdf["0.Ordnung"]/0.006, "x")
plt.plot(tempdf["mVPP"], tempdf["0.Ordnung"], "x")
plt.show()

# %%
from scipy.special import sinc
# fit sinc

def func(x, *p):
    # return p[0] * sinc(p[1]/4 * (x - p[2])/0.006)**2
    # return p[0] * np.sin(p[1] * np.sqrt(x))**2
    return p[0] * sinc(p[1]/4 * x * (1 - x))**2


x = np.linspace(min(tempdf["mVPP"]), max(tempdf["mVPP"]), 100)
y = tempdf["0.Ordnung"] / 0.791

popt, pcov = curve_fit(func, tempdf["mVPP"]/0.006, tempdf["0.Ordnung"], p0=[0.6, 0.0002, 1300])

# plot
fig, ax = plt.subplots()
ax.plot(tempdf["mVPP"]/0.006, tempdf["0.Ordnung"], "x")
ax.plot(x, func(x, *popt))
# ax.plot(x, func(x, *[0.6, 0.00002, 1300]))
plt.show()

# %%
def func(x, *p):
    return p[0] * sinc(p[1]/4 * (x - p[2]))**2
    # return p[0] * np.sin(p[1] * np.sqrt(x))**2
    # return p[0] * sinc(p[1]/4 * x/p[2] * (1 - x/p[2]))**2


x = np.linspace(min(tempdf["mVPP"]), max(tempdf["mVPP"]), 100)
y = tempdf["0.Ordnung"] / 0.791

popt, pcov = curve_fit(func, tempdf["mVPP"], y, p0=[0.8, 5, 7])
print(popt, np.sqrt(np.diag(pcov)))
# plot
fig, ax = plt.subplots()
ax.plot(tempdf["mVPP"], y, "x")
ax.plot(x, func(x, *popt))
# ax.plot(x, func(x, *[0.6, 5]))
plt.show()
# %%

def func(x, *p):
    return p[0] * np.sin(p[1] * x)**2 + p[2]

x = np.linspace(min(tempdf["mVPP"]), max(tempdf["mVPP"]), 1000) /100
y = tempdf["0.Ordnung"] / 0.791

# popt, pcov = curve_fit(func, tempdf["mVPP"], tempdf["0.Ordnung"], p0=[0.6, 0.001, 0.5])
popt, perr, ci = bootstrap(func, tempdf["mVPP"], tempdf["0.Ordnung"], p0=[0.6, 0.001, 0.5])
# popt2, pcov2 = curve_fit(func, tempdf["mVPP"], tempdf["1.Ordnung"], p0=[0.6, 0.001, 0.5])
popt2, perr2, ci2 = bootstrap(func, tempdf["mVPP"], tempdf["1.Ordnung"], p0=[0.6, 0.001, 0.5])

# Function x**(1/2)
def inverse(x):
    # return np.sqrt(x)
    return 0.200803 * (173. + 1. * np.sqrt(29929. + 9.96*10**6 * np.abs(x)))
def forward(x):
    return 2.49*10**(-6) * x**2 - 1.73*10**(-4) * x


tempdf["U2"] = forward(tempdf["mVPP"])

# x = np.array([0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,6.,6.3])
# y = func(x, *popt)

# print(x, y)

print(tempdf["mVPP"].values/100, tempdf["0.Ordnung"].values)

# %%

# plot
fig, ax = plt.subplots()
twinax = ax.twinx()
twinax.scatter(tempdf["U2"], tempdf["0.Ordnung"])
twinax.scatter(tempdf["U2"], tempdf["1.Ordnung"])
ax.plot(tempdf["mVPP"], func(tempdf["mVPP"], *popt))
# ax.fill_between(ci.x, ci.c0, ci.c1, alpha=0.5, color="C0")
# ax.plot(x[200:], func(x[200:], *popt2))
# ax.plot(x, func(x, *[0.6, 0.00002, 1300]))
# ax.plot([100, 200], [100, 400])

ax.set_xlim(-0.1, 1)
ax.set_xscale("function", functions=[forward, inverse])

# remove xticks on the top
ax.tick_params(top=False, which='both')
twinax.tick_params(bottom=False, which='both')
# twinx.
plt.show()


