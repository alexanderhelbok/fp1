from matplotlib.colors import LinearSegmentedColormap

from Source import *
import os
from scipy.signal import argrelextrema, medfilt
from PIL import Image
# Define a Gaussian function to fit the data

#%%
def gauss(x, A, mu, sigma):
    return A * np.exp(-(x-mu)**2 / (2*sigma**2))

# Load the image
# image_path = '/home/taco/Documents/fp1/04_Gaußstrahlen/data/TM00_secondone.jpg'  # Replace with your image path
image_path = '/home/taco/Documents/fp1/04_Gaußstrahlen/data/WIN_20240410_11_25_21_Pro.jpg'  # Replace with your image path

image = Image.open(image_path)

original_image = np.asarray(image)
# Convert image to grayscale
gray_image = image.convert('L')
 
# Convert the grayscale image to a numpy array
image_array = np.array(gray_image)

# Sum the intensity of each column change axis to 1 if row is needed
column_intensity = image_array.sum(axis=0)


min_intensity_column = np.min(column_intensity)
max_intensity_column = np.max(column_intensity)
column_intensity = (column_intensity - min_intensity_column) / (max_intensity_column - min_intensity_column)

# Get the number of the column for each sum
columns = np.arange(len(column_intensity))

# Fit the summed column intensity data to a Gaussian
popt, pcov = curve_fit(gauss, columns, column_intensity, p0=[max(column_intensity), len(column_intensity)/2, len(column_intensity)/4])

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(columns, column_intensity, label='normalized Data', color='blue')
plt.plot(columns, gauss(columns, *popt), label='Gauss Fit', color='red')
plt.title('Summed Intensity of Each Column Fitted with Gaussian')
plt.xlabel('Column number')
plt.ylabel('Summed Intensity')
plt.legend()
plt.show()  # This will display the plot

# %%
plt.imshow(original_image)
# %%#
x = []
y = []
# redraw data from image_array
for i in range(len(image_array)):
    for j in range(len(image_array[i])):
        x += [j] * image_array[i][j]
        # x.append(j)
        # y.append(i)
        y += [i] * image_array[i][j]


column_intensity_x = image_array.sum(axis=0) / max(image_array.sum(axis=0))
column_intensity_y = image_array.sum(axis=1) / max(image_array.sum(axis=1))

poptx, pcovx = curve_fit(gauss, np.arange(len(column_intensity_x)), column_intensity_x)
popty, pcovy = curve_fit(gauss, np.arange(len(column_intensity_y)), column_intensity_y)
# %%
# Define the colormap colors
colors = [(0, 0, 0), (1, 0, 0)]  # black to red

# Create the colormap
cm = LinearSegmentedColormap.from_list('black_red', colors)

# use seabron to create scatterplot of endpoint with marginal histograms
df = pd.DataFrame({"x": x, "y": y})
g = sns.jointplot(data=df, x="x", y="y", kind="hist", marginal_kws=dict(bins=700, stat='density'), color="white", alpha=0)
g.ax_joint.imshow(original_image)

# remove grid
g.ax_joint.grid(False)

ax_marg_x = g.ax_marg_x
ax_marg_y = g.ax_marg_y
# ax_marg_x.cla()
# ax_marg_y.cla()

# ax_marg_x.scatter(np.arange(len(column_intensity_x)), column_intensity_x, color='blue')
# ax_marg_y.scatter(column_intensity_y, np.arange(len(column_intensity_y)), color='blue')
# ax_marg_x.plot(np.arange(len(column_intensity_x)), gauss(np.arange(len(column_intensity_x)), *poptx), color='red')
# ax_marg_y.plot(gauss(np.arange(len(column_intensity_y)), *popty), np.arange(len(column_intensity_y)), color='red')

# remove tick labels on marginal axes (and main)set_square_lim
ax_marg_x.tick_params(bottom=False, which='both', top=False, left=False, right=False)
ax_marg_y.tick_params(left=False, which='both', right=False, top=False, bottom=False)
# ax_marg_x.set_xticks([])
# ax_marg_y.set_yticks([])
g.ax_joint.tick_params(top=False, right=False, which='both')
# %%
# set square limits for the joint plot
# max_lim = set_square_lim(g.ax_joint)

# X = np.linspace(-max_lim, max_lim, 1000)

ax_marg_x = g.ax_marg_x
ax_marg_y = g.ax_marg_y

# plot the marginal histograms and the analytical normal distribution
ax_marg_x.plot(X, norm.pdf(X, 0, np.sqrt(n_runs/2)), color='red')
ax_marg_y.plot(norm.pdf(X, 0, np.sqrt(n_runs/2)), X, color='red')

# remove tick labels on marginal axes (and main)
ax_marg_x.tick_params(bottom=False, which='both', top=False)
ax_marg_y.tick_params(left=False, which='both', right=False)
g.ax_joint.tick_params(top=False, right=False, which='both')


print(np.mean(df.x), np.std(df.x), np.mean(df.y), np.std(df.y))
plt.tight_layout()
# plt.savefig("bilder/random_walk_2d_joint.pdf", bbox_inches="tight")
plt.show()