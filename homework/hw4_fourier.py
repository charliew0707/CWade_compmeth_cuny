"""
HW4 – Fourier Analysis of an Eclipsing Binary Light Curve

Loads a TESS light curve (FITS file) for an eclipsing binary star, selects a
single well-sampled observing epoch, and computes its Fourier transform and
power spectrum. The light curve is then reconstructed with the inverse
transform using as few Fourier coefficients as possible to find the minimum
needed to capture the eclipsing behavior. The analysis is repeated on a second
epoch to check that the results are consistent.

Because TESS cadence is nearly even but has gaps from missing observations,
the gaps are filled by linear interpolation onto a uniform time grid and the
analysis is redone to see how much the gaps affect the result.

Data source: https://tessebs.villanova.edu/
"""

# Star File Name: 0007849727.fits

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

here = Path(__file__).parent
hdul = fits.open(here / "0007849727.fits")

hdul.info()                           # file structure
print(hdul[1].columns.names)          # check the real column names

times = hdul[1].data['times']
fluxes = hdul[1].data['fluxes']
ferrs = hdul[1].data['ferrs']

# Plot the entire light curve, this shows dense epochs with data
plt.figure(figsize=(12, 5))
plt.plot(times, fluxes, '.', ms=2)
plt.xlabel("Time")
plt.ylabel("Flux")
plt.title("Full Light Curve: 0007849727")
plt.show()


# specify data to range of data used: 2710 to 2780


# FFT the above clean data:


# create and plot power spectrum 

