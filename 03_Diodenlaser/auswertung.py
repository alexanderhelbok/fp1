# from Source import *
import joypy
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import pandas as pd
import numpy as np
# %%
# load data
# fill df with 0
df = pd.DataFrame(columns=["t", "V", "I"])
currentdf = pd.read_csv('/home/taco/Documents/fp1/03_Diodenlaser/data/Istepsize.csv')

for i in range(1, 27):
    tempdf = pd.read_csv(f'/home/taco/Documents/fp1/03_Diodenlaser/data/T25/Current{i}.CSV')
    # tempdf['Current'] = currentdf.iloc[i-1, 1]
    tempdf['Current'] = currentdf.iloc[26 - i, 1]
    tempdf.columns = ['t', 'V', 'I']
    tempdf.V += 0.02
    # tempdf.V = savgol_filter(tempdf.V, 101, 2)
    # vcat
    df = pd.concat([df, tempdf])

# for i in range(20, 30):
#     tempdf = pd.read_csv(f'/home/taco/Documents/fp1/03_Diodenlaser/data/Tempsweep/T{i}.CSV')
#     tempdf['Current'] = currentdf.iloc[i-20, 1]
#     tempdf.columns = ['t', 'V', 'I']
#
#     df = pd.concat([df, tempdf])

df.t = df.t - df.t.min()
print(df.head())
# %%
# plot data
fig, axes = joypy.joyplot(df, by="I", figsize=(10, 10), legend=False, kind="values", column="V", x_range=[0, 6000], overlap=3)
axes[-1].set_xticks(range(0, 6000, 1000), ["0", "1", "630", "3", "4", "1260"])

# axes[12].set_ymajor_locator(plt.MaxNLocator(3))
# for (i, ax) in enumerate(axes):
#     ax.set_ylabel(f"Current: {currentdf.iloc[i, 1]} mA")

plt.xlabel("frequency (GHz)")
# plt.ylabel("current (mA)")

# %%
# axes[12].yaxis.majorTicks[0].label1 = 2

plt.show()
