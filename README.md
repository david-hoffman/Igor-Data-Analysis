#ReadMe file for Igor Data Analysis

The license for the code contained in this repository is contained in included the license.txt file.

I have written these IgorPro (WaveMetrics Inc.) procedures to analyze data collected during my doctoral research, as such the user interface is not very friendly and the documentation is not complete. My research primarily uses ultrafast lasers to perform two main types of vibrational spectroscopy 1.) Femtosecond Stimulated Raman Spectroscopy, FSRS (DPH, et al., PCCP 2012, 14, 6298 and DPH, et al., JPCC 2013 in press), 2.) and Impulsive Stimulated Raman Spectroscopy, ISRS (DPH, et al., JPCA 2013 in preparation).

- FitFuncs.ipf
A set of useful fitting functions for the ultrafast spectroscopist. Many different kinds of exponentials convoluted with a gaussian instrument response. NOTE: the width of the gaussian IRF is defined to be the FWHM/(2*sqrt(ln(2))) which is consistent with Igor's built-in Gaussian function within the CurveFit dialog meaning that this parameter can be used directly.

- PerkinElmerImport.ipf
This file contains a single procedure, LoadPEData(). This procedure will ask the user which files to load and then will load Perkin Elmer's proprietary binary format into waves named after the file. It will include the experimental info stored in the binary file in the created wave's note. The procedure will also display the loaded waves.

- SignalProcessing.ipf
This Igor Procedure File (.ipf) contains procedures useful for analyzing ISRS data. In fact, these procedures should be useful for analyzing any data which can be reasonably described as a sum of damped sinusoids.
