# X-rayAnalysis 
These programs execute the $Swift$ X-ray analysis for a single object in the following manner: 
  - Conduct a query in HEASARC for a given RA, DEC for a given named object within 8 arc sec offset.
  - Download the complete Window Timing (WT) and Photon Counting (PC) data for these sources.
  - Run **swift_analysis.py** to do analysis for WT mode, which cleans the files, create images, light cuves spectra and group them for further analysis.
  - Run **spec_fit.py** to conduct the spectral fit using XSPEC software with predefined model as powerlaw. One can modify it to add any models. This uses a tcl script file called tclex.xcm
  - This program generates various outputs, including the fitting results for each observation as well as a combined file which lists all the parameters for the fitted model as well as the calculated flux. 
  - Repeat all the above steps for the PC mode by using similar .py files, but with an addition of "_pc"
  
