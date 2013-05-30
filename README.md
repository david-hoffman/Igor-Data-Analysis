#Igor Data Analysis

The license for the code contained in this repository is contained in included the [`license.txt`](https://github.com/david-hoffman/Igor-Data-Analysis/blob/master/License.txt) file.

I have written these [IGOR Pro](http://www.wavemetrics.com/index.html) procedures to analyze data collected during my doctoral research, as such the user interface is primarily through the command window and the documentation is incomplete. Nevertheless, I have found these procedures to be extraordinarily useful and I'd like to share them with the rest of the research community.

My research primarily uses ultrafast lasers to perform two main types of vibrational spectroscopy 1.) Femtosecond Stimulated Raman Spectroscopy, FSRS ([DPH, et al., *PCCP* **2012**](pubs.rsc.org/en/content/articlehtml/2012/cp/c2cp23468h) and [DPH, et al., *JPCC* **2013**](http://pubs.acs.org/doi/abs/10.1021/jp400369b)), 2.) and Impulsive Stimulated Raman Spectroscopy, ISRS (DPH, et al. *in preparation*).

##Utilities.ipf
A set of utilities for the processing of time resolved vibrational spectroscopic data.

##DataImport.ipf
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

The `LoadUVVisData()` procedure can be used to load and file in the JCAMP (.DX) format.

`LoadCWData()` can be used to load any generic txt file.

The LabVIEW programs which generate the FSRS data are available [here](https://www.github.com/david-hoffman/FSRS-LabVIEW).

##FitFuncs.ipf
A set of useful fitting functions for the ultrafast spectroscopist. Many different kinds of exponentials convoluted with a gaussian instrument response. 

_**NOTE:** the width of the gaussian IRF is defined to be the FWHM/(2*sqrt(ln(2))) which is consistent with Igor's built-in Gaussian function within the CurveFit dialog meaning that this parameter can be used directly._

##PerkinElmerImport.ipf
This file contains a single procedure, `LoadPEData()`. This procedure will ask the user which files to load and then will load Perkin Elmer's proprietary binary format into waves named after the file. It will include the experimental info stored in the binary file in the created wave's note. The procedure will also display the loaded waves.

##SignalProcessing.ipf
This Igor Procedure File (.ipf) contains procedures useful for analyzing ISRS data. In fact, these procedures should be useful for analyzing any data which can be reasonably described as a sum of damped sinusoids. The two main procedures are:
- `LPSVD(signal,[M,LFactor,RemoveBias])` which fits the data, in a *linear* least squares sense using the **Linear Prediction with Singular Value Decomposition** (LPSVD) algorithm.
- `Cadzow(signal, M, iters,[lFactor,q])` which filters the data using Cadzow's Composite Property Mapping Algorithm.

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/david-hoffman/Igor-Data-Analysis/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/067a677b64204640d9421177434c208e "githalytics.com")](http://githalytics.com/david-hoffman/Igor-Data-Analysis)
