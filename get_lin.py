import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

# Define the instrument name:
instrument = 'niriss'
# Define the precentage of non-linearity we are searching for:
x = 0.1

# First load the reference files for a given instrument/detector
# (one can download all these files from https://jwst-crds.stsci.edu):
if instrument == 'niriss':
    rfile = 'data/jwst_niriss_linearity_0011.fits'

def convert(raw_signal, coeff):
    """
    Given an input raw signal, returns the converted signal using the 
    non-linearity polynomial correction
    """
    return np.polyval(coeff[::-1],raw_signal)

# Extract data and headers. Data "d" is a [norder,np1,np2] array containing 
# norder coefficients of a polynomial transformation for each of the np1xnp2 
# pixels --- these coefficients are used to correct raw signal impacted by 
# non-linearity on the pixels.
d,h = fits.getdata(rfile,header=True)

# Define a solution matrix, i.e., a matrix that will save where the x% non-linearity occurs 
# for the different pixles:
solution_matrix = np.zeros([d.shape[1],d.shape[2]])

# Now we generate a range of counts on which to look for the x% level of non-linearity. Here 
# we will define it as in this ISR: http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2008-39.pdf. 
# Idea is to be as consistent as with the HST/WFC3 definition, which is "(...) the saturation 
# level is defined as the signal level at which the detector response becones non-linear by 5%". The 
# interpretation here is, according to that same ISR, that this is the signal level at which the measured 
# signal deviates from its linear behaviour by more than 5%. 
#
# Note we want to know the *linearized* signal at which this occurs, in order to compare it with the ETC results.
signal = np.arange(1000,66000,1.)

# Now iterate through each pixel, saving resulting "saturation at the x% level" on each pixel:
for i in range(d.shape[1]):
    for j in range(d.shape[2]):
        # Generate corrected signal:
        linearized_signal = convert(signal, d[:,i,j])
        # Calculate the percentage offset from the linear signal of the original signal (as in, e.g., 
        # Figure 2 of the above ISR):
        p = (linearized_signal/signal) - 1.
        # Calculate where the x% level of non-linearity deviation happens:
        idx = np.where(p<x)[0]
        if len(idx>0):
            print(i,j,linearized_signal[idx[-1]])
            solution_matrix[i,j] = linearized_signal[idx[-1]]

# Save the "x% saturation levels" image for future post-processing:
hdu = fits.PrimaryHDU(solution_matrix)
hdu1 = fits.HDUList([hdu])
hdu1.writeto('sm_'+str(x)+'_'+instrument+'.fits')
