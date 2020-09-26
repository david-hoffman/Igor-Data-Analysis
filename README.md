# Igor Data Analysis

I have written these [IGOR Pro](http://www.wavemetrics.com/index.html) procedures to analyze data collected during my doctoral research, as such the user interface is primarily through the command window and the documentation is incomplete. Nevertheless, I have found these procedures to be extraordinarily useful and I'd like to share them with the rest of the research community.

My research primarily uses ultrafast lasers to perform two main types of vibrational spectroscopy 1.) Femtosecond Stimulated Raman Spectroscopy, FSRS, 2.) and Impulsive Stimulated Raman Spectroscopy, ISRS (in preparation). The relevant publications can be found [here](http://scholar.google.com/citations?user=HGG__poAAAAJ&hl=en).

## Utilities.ipf
A set of utilities for the processing of time resolved vibrational spectroscopic data.

Some notes on specific functions

- `GroundSubtract2(timepoints,ground,startpt,endpt,type,[q])`: This function loops through all the time delays subtracting the ground state spectrum from the raw excited state spectra by scaling the excited state spectra to the ground state using a solvent feature. In this way fluctuations in Raman pump power are controlled for using the solvent as an internal standard. This is particularly important when the Raman pump is near a transient absorption feature as this algorithm will correct for any transient attenuation of the Raman pump.

- `SolventSubtract2(spectrum, solvent, sp, ep, type)`: This function scales the `solvent` _to_ the `spectrum`. This is necessary to be internally consistend with `GroundSubtract2`.

## DataImport.ipf
Functions designed to import data from a FSRS instrument.

This procedure file also adds items to the "Load Waves" menu:

	Menu "Load Waves"
		"-"
		"FSRS Data ... /1", LoadFSRSData()
		"Raw FSRS Data ... /S 1", LoadRawFSRSData()
		"-"
		"DAQ Scan Data ... /2", LoadDAQData()
		"Detector XC Data ... /S 2", LoadDetXCData()
		"-"
		"CW Data ... ", LoadCWData()
		"UV-vis Data ...", LoadUVvisData()
	End

The `LoadUVVisData()` procedure can be used to load any file in the JCAMP (.DX) format.

`LoadCWData()` can be used to load any generic txt file.

The LabVIEW programs which generate the FSRS data are available [here](https://www.github.com/david-hoffman/FSRS-LabVIEW).

## FitFuncs.ipf
A set of useful fitting functions for the ultrafast spectroscopist. Many different kinds of exponentials convoluted with a gaussian instrument response. 

_**NOTE:** the width of the gaussian IRF is defined to be the FWHM/(2*sqrt(ln(2))) which is consistent with Igor's built-in Gaussian function within the CurveFit dialog meaning that this parameter can be used directly._

## FitTimeSeries.ipf
Contains `FitTimeSeries` a versatile batch fitting function. Fits a sum of either Gaussians or Lorentzians to a subset of time resolved FSRS spectra. `Coefs` is a wave of coefficiencts for the peaks, i.e. `Coefs={0,Amp1,Center1,FWHM1,...,Ampn,Centern,FWHMn}`. The `0` is necessary for backwards compatibility. `Wavenumber` is a *text* wave containing the peak labels, this will determine the names of the output waves. I usually just use the approximate center frequency. **Make sure that your `Coefs` wave has the correct number of coefficients relative to your `Wavenumber` wave**.

Options include:

- `PolyNum` the order of the polynomial baseline. To remove the baseline input `PolyNum=0`. The default is cubic.

- `Plot=1` will plot the results in a pleasing manner.

- `PeakType=1` is Gaussian (default behaviour). `PeakType=0` is Lorentzian.

## FitTimeSeriesDisp.ipf
Contains `FitTimeSeries2` which is almost identical to `FitTimeSeries` except that it uses a dispersive lorentzian lineshape function which has one parameter for the amplitude of the real part and one for the amplitude of the dispersive part. If the initial values for either of these components is 0 than the function assumes that those peaks are totally real or dispersive and holds that component to 0.

## PerkinElmerImport.ipf
This file contains a single procedure, `LoadPEData()`. This procedure will ask the user which files to load and then will load Perkin Elmer's proprietary binary format into waves named after the file. It will include the experimental info stored in the binary file in the created wave's note. The procedure will also display the loaded waves.

## SignalProcessing.ipf
This Igor Procedure File (.ipf) contains procedures useful for analyzing ISRS data. In fact, these procedures should be useful for analyzing any data which can be reasonably described as a sum of damped sinusoids. The two main procedures are:
- `LPSVD(signal,[M,LFactor,RemoveBias])` which fits the data, in a *linear* least squares sense using the **Linear Prediction with Singular Value Decomposition** (LPSVD) algorithm. Estimates of the variances for the parameters returned by LPSVD are calculated as the [Cramer-Rao bound](http://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound) and returned in a wave named `sigma_LPSVD_coefs`.
- `Cadzow(signal, M, iters,[lFactor,q])` which filters the data using Cadzow's Composite Property Mapping Algorithm.
