jwst linearity corrections
--------------------------

Some details regarding reference files:

- For NIRSpec, there are two files because there are two detectors: `jwst_nirspec_linearity_0024.fits` corresponds to NS1 
  and `jwst_nirspec_linearity_0023.fits` corresponds to NS2.

- For NIRCam there are several files, because there are 8 detectors on the short wavelength channel and 2 in the long-wavelength 
  one (see https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview) for an overview on 
  the detector placing/distribution/configuration. There are two files for the long-wavelenght channel 
  (`jwst_nircam_linearity_0052.fits` for detector A, `jwst_nircam_linearity_0049.fits` for detector B) and eight files for 
  the short-wavelength channel ()
