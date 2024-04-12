from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import hermite

from Source import *
import os
from scipy.signal import argrelextrema, medfilt
from PIL import Image
# Define a Gaussian function to fit the data

#%%
def gauss(x, A, mu, sigma):
    return A * np.exp(-(x-mu)**2 / (2*sigma**2))

# Define a function for the Gauss-Hermite polynomial
def hermite_poly(n, x):
    Hn = hermite(n, monic=True)
    return Hn(x)

# Define the Gauss-Hermite function for fitting
def gauss_hermite1(x, A, x0, sigma):
    # Use the hermite_poly function for the Hermite polynomial part
    return A * np.exp(-(x-x0)**2 / (2*sigma**2)) * (hermite_poly(1, (x-x0)/sigma))**2


def gauss_hermite2(x, A, x0, sigma):
    # Use the hermite_poly function for the Hermite polynomial part
    return A * np.exp(-(x-x0)**2 / (2*sigma**2)) * (hermite_poly(2, (x-x0)/sigma))**2


def gauss_hermiten(x, A, x0, sigma, n):
    # Use the hermite_poly function for the Hermite polynomial part
    return A * np.exp(-(x-x0)**2 / (2*sigma**2)) * (hermite_poly(n, (x-x0)/sigma))**2

# Load the image
image_names =  ['TM00_secondone.jpg', 'TM01.jpg', 'WIN_20240410_11_25_21_Pro.jpg', 'WIN_20240410_11_16_59_Pro.jpg']  # Replace with your image path
fig, ax = plt.subplots(2, 2, figsize=(10, 8), gridspec_kw={'hspace': 0.}, sharex=True, sharey=True)

# increase label font size
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

for (i, name) in enumerate(image_names):
    image = Image.open(f"/home/taco/Documents/fp1/04_Gaußstrahlen/data/{name}")

    original_image = np.asarray(image)
    # Convert image to grayscale
    gray_image = image.convert('L')

    # Convert the grayscale image to a numpy array
    image_array = np.array(gray_image) - np.min(np.array(gray_image))

    # Sum the intensity of each column change axis to 1 if row is needed
    column_intensity = image_array.sum(axis=0)


    min_intensity_column = np.min(column_intensity)
    max_intensity_column = np.max(column_intensity)
    column_intensity = (column_intensity - min_intensity_column) / (max_intensity_column - min_intensity_column)

    # Get the number of the column for each sum
    columns = np.arange(len(column_intensity))

    column_intensity_x = image_array.sum(axis=0) / max(image_array.sum(axis=0))
    column_intensity_y = image_array.sum(axis=1) / max(image_array.sum(axis=1))

    # fig, ax = plt.subplots(figsize=(10, 6))
    ax[i%2, i//2].set_aspect(1.)

    # ====== prep axes
    divider = make_axes_locatable(ax[i%2, i//2])
    # below height and pad are in inches
    ax_histx = divider.append_axes("top", 0.5, pad=0.1, sharex=ax[i%2, i//2])
    ax_histy = divider.append_axes("right", 0.5, pad=0.1, sharey=ax[i%2, i//2])

    ax[i%2, i//2].tick_params(which='both', top=False, right=False)
    ax[i%2, i//2].tick_params(which='both', direction="out")

    for axis in [ax_histx, ax_histy]:
        # remove spines
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        # remove ticks
        axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    ax_histy.spines['bottom'].set_visible(False)
    ax_histx.spines['left'].set_visible(False)


    # ====== plot
    if i in [0, 1, 3]:
        poptx, pcovx = curve_fit(gauss, np.arange(len(column_intensity_x)), column_intensity_x)
        ax_histx.plot(gauss(np.arange(len(column_intensity_x)), *poptx), color='crimson', lw=2)
    else:
        poptx, pcovx = curve_fit(gauss_hermite1, np.arange(len(column_intensity_x)), column_intensity_x, p0=[1.5, 800, 150])
        print(poptx)
        # ax_histx.plot(gauss_hermite(np.arange(len(column_intensity_x)), *poptx), color='crimson', lw=2)
        ax_histx.plot(gauss_hermite1(np.arange(len(column_intensity_x)), *poptx), color='crimson', lw=2)
    if i in [0, 2]:
        popty, pcovy = curve_fit(gauss, np.arange(len(column_intensity_y)), column_intensity_y)
        ax_histy.plot(gauss(np.arange(len(column_intensity_y)), *popty), np.arange(len(column_intensity_y)), color='crimson', lw=2)
    elif i == 1:
        popty, pcovy = curve_fit(gauss_hermite1, np.arange(len(column_intensity_y)), column_intensity_y, p0=[1, 400, 150])
        print(popty)
        ax_histy.plot(gauss_hermite1(np.arange(len(column_intensity_y)), *popty), np.arange(len(column_intensity_y)), color='crimson', lw=2)
    else:
        popty, pcovy = curve_fit(gauss_hermite2, np.arange(len(column_intensity_y)), column_intensity_y, p0=[1, 400, 150])
        print(popty)
        ax_histy.plot(gauss_hermite2(np.arange(len(column_intensity_y)), *popty), np.arange(len(column_intensity_y)), color='crimson', lw=2)


    ax[i%2, i//2].imshow(original_image, extent=[0, original_image.shape[1], 0, original_image.shape[0]])

    ax_histx.scatter(np.arange(len(column_intensity_x)), column_intensity_x, color='black', s=15)
    ax_histy.scatter(column_intensity_y, np.arange(len(column_intensity_y)), color='black', s=15)

    if i == 0:
        ax_histx.set_title("$\mathrm{TEM}_{00}$ mode", pad=15, fontsize=20)
    elif i == 1:
        ax_histx.set_title("$\mathrm{TEM}_{01}$ mode", pad=15, fontsize=20)
    elif i == 2:
        ax_histx.set_title("$\mathrm{TEM}_{10}$ mode", pad=15, fontsize=20)
    elif i == 3:
        ax_histx.set_title("$\mathrm{TEM}_{02}$ mode", pad=15, fontsize=20)

ax[1, 0].set_xlabel("pixel")
ax[1, 1].set_xlabel("pixel")
ax[0, 0].set_ylabel("pixel", labelpad=15)
ax[1, 0].set_ylabel("pixel", labelpad=15)

plt.tight_layout()
# plt.savefig("/home/taco/Documents/fp1/04_Gaußstrahlen/bilder/gaussmodes.pdf", bbox_inches='tight')
plt.show()


# %%
plt.plot(gauss_hermiten(np.arange(len(column_intensity_x)), 1, 400, 100, 6))
plt.show()

