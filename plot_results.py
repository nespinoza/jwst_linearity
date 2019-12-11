import matplotlib.pyplot as plt
plt.style.use('ggplot')
from astropy.io import fits
import numpy as np
import glob

instrument = 'nirspec-ns2'
plot_median = True

def tick_function(X, gain = 1.):
    """
    Convert tick labels from counts to electrons. 
    Taken/modified from: https://stackoverflow.com/questions/10514315/how-to-add-a-second-x-axis-in-matplotlib
    """
    V = gain*X
    return ["%.0f" % z for z in V]

if instrument == 'niriss':
    title = 'NIRISS'
    gain = 1.6 # Approximate gain of NIRISS in e/counts (https://jwst-docs.stsci.edu/near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-detector-overview)
elif instrument == 'nirspec-ns1':
    title = 'NIRSpec (NS1)'
    gain = 0.996 # From https://jwst-docs.stsci.edu/near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-performance#NIRSpecDetectorPerformance-Detectorgain
elif instrument == 'nirspec-ns2':
    title = 'NIRSpec (NS2)'
    # FULL mode gain:
    gain = 1.137 # From https://jwst-docs.stsci.edu/near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-performance#NIRSpecDetectorPerformance-Detectorgain
    # GENERIC mode gain:
    gain = gain*1.43
elif instrument == 'nircam-lwA' or instrument == 'nircam-lwB':
    detector = instrument[-1:]
    title = 'NIRCam (long-wavelength, detector '+detector+')'
    gain = 1.82 # From https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-performance
else:
    print('Warning: instrument not recognized/gain not defined. Electron values might be off as this assumes gain = 1.')
# Load result matrices:
files = glob.glob('*'+instrument+'.fits')
files.sort()
# Define plotting arguments:
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
# Iterate through the solution matrices:
for f in files:
    m = fits.getdata(f)
    # Flatten data matrix:
    mf = m.flatten()
    # Ommit pixels with corrections less than 1000, and larger than 65000:
    idx = np.where((mf>1000)&(mf<65000.))[0]
    mf = mf[idx]
    # Now plot distribution of pixels at which "x% saturation" happens:
    percent = np.double(f.split('_')[-2])*100.
    label_string = '{0:.0f} percent deviation'.format(percent)
    s = ax1.hist(mf,bins=100,label = label_string,alpha=0.8)
    if plot_median:
        # Plot line indicating median counts:
        max_bin = np.max(s[0])
        median_counts = np.median(mf)
        ax1.plot([median_counts,median_counts],[0,max_bin],'r--')
        ax1.text(median_counts + 1000,max_bin,'Median: {0:.0f} counts'.format(median_counts))
        ax1.text(median_counts + 1000,max_bin - 3000.,'(~{0:.0f} e$^-$)'.format(median_counts*gain))
ax1.set_title(title)
ax1.set_xlim([0.,70000.])
ax1.set_xlabel('Number of counts')
ax1.set_ylabel('Number of pixels with counts')
ax2.set_xlim(ax1.get_xlim())
new_tick_locations = np.arange(0,70000,10000)
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations,gain))
ax2.set_xlabel(r"Number of electrons")
ax1.legend()
plt.tight_layout()
plt.savefig(instrument+'.png')
