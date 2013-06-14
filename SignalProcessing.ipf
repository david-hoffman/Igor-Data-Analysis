#pragma rtGlobals=3		// Use modern global access method.

Constant kSpeedOfLight=2.99792458e-5 //(cm/fs)

//LPSVD was developed by Tufts and Kumaresan (Tufts, D.; Kumaresan, R. IEEE Transactions on Acoustics,
//Speech and Signal Processing 1982, 30, 671 Ð 675.) as a method of harmonic inversion, i.e. decomposing
//a time signal into a linear combination of (decaying) sinusoids.
//
//A great reference that is easy to read for the non-EECS user is:
// Barkhuijsen, H.; De Beer, R.; BovŽe, W. M. M. .; Van Ormondt, D. J. Magn. Reson. (1969) 1985, 61, 465Ð481.
//
//This particular implementation was adapted, in part, from matNMR by Jacco van Beek
//http://matnmr.sourceforge.net/
//and  Complex Exponential Analysis by Greg Reynolds
//http://www.mathworks.com/matlabcentral/fileexchange/12439-complex-exponential-analysis/
//
//Author: David Hoffman (david.hoffman@berkeley.edu)
//Date: Feb, 2013

Function LPSVD(signal,[M,LFactor,RemoveBias])
	//A function that performs the linear prediction-singular value decomposition
	//of a signal that is assumed to be a linear combination of damped sinusoids
	Wave signal						//The signal to be analyzed
	Variable M						//Model order, user selectable
	Variable LFactor				//How to size the Hankel matrix, Tufts and Kumaresan suggest 1/3-1/2
	
	Variable RemoveBias			//if this variable is anything but default than bias will be removed from
									//the singular values of A
									
	//We want this function to be data folder aware
	//We'll do all our calculations in a specific data folder and then kill that folder at the end
	String savDF= GetDataFolder(1)	// Save current DF for restore.
	if( DataFolderExists("root:LPSVD_Data") )
		SetDataFolder root:LPSVD_Data
	else 
		NewDataFolder/O/S root:LPSVD_Data	// Our stuff goes in here.
	endif
	
	//Default number of prediction coefficients is half the number of points in the input wave
	If(ParamIsDefault(LFactor))
		LFactor = 1/2
	EndIf
	
	If(Lfactor>3/4)
		Print "You attempted to use an Lfactor greater than 3/4, it has been set to 3/4"
		LFactor=3/4
	EndIf
	//Timing
	Variable timerRef=startMSTimer
	
	Variable N=numpnts(signal)	//length of signal
	Variable L=floor(N*LFactor)	//Sizing of the Hankel matrix, i.e. the backward prediction matrix
	
	Rotate -1,Signal				//Shift the signal forward by 1
	WAVE  A = $Hankel(Signal,N-L,L)	//Generate the Hankel matrix
	Rotate 1,Signal					//Shift the signal back
	
	MatrixOp/O/C  A = Conj(A)		//Take the conjugate of the Hankel Matrix to form the prediction matrix
	
	Duplicate/C/O/R=[0,N-L-1] Signal, h	//Set up the data vector, the vector to be "predicted"
	h = Conj(h)						//Take the conjugate
	
	MatrixSVD A			//Perform an SVD on the Hankel Matrix
	
	//Create local variable for the output of the SVD operation
	WAVE/C S = W_W
	WAVE/C U = M_U
	WAVE/C VT = M_VT
	
	If(ParamIsDefault(M))			//We can estimate the model order if the user hasn't selected one
		M  = estimate_model_order(S,N,L)+8
		If(M>numpnts(S))
			Print "Estimation failed: M is set to max"
			M = numpnts(S)
		EndIf
		Print "Estimated model order: "+num2str(M)
	EndIF
	
	If(M>numpnts(S))
		M = numpnts(S)
		Print "M too large, set to max = "+num2str(M)
	EndIf
	
	//Now we can generate the LP coefficients
	Make/C/D/O/N=(L) LPCoefs=0
	
	Variable k=0
	Make/C/D/O/N=(DimSize(U,0)) U_K=0
	Make/C/D/O/N=(DimSize(VT,1)) VT_K=0
	
	//Remove bias
	If(ParamIsDefault(RemoveBias))
		RemoveBias = 1
	EndIf
	
	If(RemoveBias)
		//Here we subtract the arithmatic mean of the singular values determined to be
		//noise from the rest of the singular values as described in Barkhuijsen
		S -= mean(S,M,inf)
	EndIf
	
	
	S = 1/S*(p<M)	//invert S
	//Redimension the matrices to speed up the matrix multiplication step
	redimension/N=(M,-1) VT	//Make VT the "right" size
	Redimension/N=(-1,M) U	//Make U the "right" size
	MatrixOp/C/O LPCoefs= cmplx(-1,0)*VT^h x DiagRC(S,M,M) x U^h x h
	
	//Error check: are there any NaNs or INFs in LPCoefs?
	Variable V_numNans, V_numINFs
	WaveStats/Q LPCoefs
	If(V_numNans!=0 || V_numINFs!=0)
		Print "There has been an error generating the prediction-error filter polynomial"
		SetDataFolder savDF	// Restore current DF.
		timerRef = stopMSTimer(timerRef)
		return -2
	EndIf
	
	//I need to add 1 to the beginning of LPCoefs before taking roots
	InsertPoints 0, 1, LPCoefs
	LPCoefs[0]=1
	
	//I can now find the roots of B (assuming B represents the coefficients of a polynomial)
	Wave/C myRoots = $roots(LPCoefs)
	
	//Remove the poles that lie within the unit circle on the complex plane as directed by Kurmaresan
	//Actually it seems the correct thing to do is to remove roots with positive damping constants
	Duplicate/C/O myRoots usedRoots
	usedRoots = conj(ln(usedRoots))
	For(k=0;k<numpnts(myRoots);k+=1)
		//Use Kumerasan's criterion to separate components due to noise
		If(cabs(myRoots[k])>=1)
			usedRoots[k]=NaN
		EndIf
	EndFor
	//Remove the spurious roots
	WaveTransform/O zapNANs, usedRoots
	
	//Error checking: see if we removed all roots!
	If(numpnts(usedRoots)==0)
		Print "There has been an error finding the real poles"
		SetDataFolder savDF	// Restore current DF.
		timerRef = stopMSTimer(timerRef)
		return -1
	EndIf
	
	//Before continuing lets sort everything by the frequencies
	Make/D/O/N=(numpnts(usedRoots)) sortKeys
	sortKeys = imag(usedRoots)
	Sort sortKeys usedRoots
	
	//Lets make a matrix with dimension labels to store all our parameters
	Make/D/O/N=(numpnts(usedRoots),4) LPSVD_coefs = 0
	SetDimLabel 1, 0, amps, LPSVD_coefs
	SetDimLabel 1, 1, freqs, LPSVD_coefs
	SetDimLabel 1, 2, damps, LPSVD_coefs
	SetDimLabel 1, 3, phase, LPSVD_coefs
	
	//We can directly convert our poles into estimated damping factors and frequencies
	LPSVD_coefs[][%damps] = real(usedRoots[p])
	LPSVD_coefs[][%freqs] = imag(usedRoots[p])/(2*pi)
	
	//But we need to do a little more work to get the predicted amplitudes and phases
	//Here we generate our basis matrix
	Make/C/D/O/N=(numpnts(signal),numpnts(usedRoots)) basis
	basis = exp(p*usedRoots[q])
	//Take the inverse
	Wave pinvBasis = $pinv(basis)
	//And apply it to our signal to recover our predicted amplitudes
	//Amps here are complex meaning it has amplitude and phase information
	MatrixOp/C/O cAmps = pinvBasis x signal
	LPSVD_coefs[][%amps]=Real(r2polar(cAmps[p]))
	LPSVD_coefs[][%phase]=Imag(r2polar(cAmps[p]))
	
	//Calculate the errors
	WAVE Errors = $CalculateLPSVDError(LPSVD_coefs,Signal)
	
	//Move the results to the original data folder
	Duplicate/O LPSVD_coefs, $(savDF+"LPSVD_Coefs")
	Duplicate/O Errors, $(savDF+"sigma_LPSVD_Coefs")
	
	SetDataFolder savDF	// Restore current DF.
	
	//Clean up
	KillDataFolder root:LPSVD_Data
	
	//Output the calculation time
	Print "LPSVD time: "+num2str(stopMSTimer(timerRef)*1e-6)+" seconds"
End

Function/S CalculateLPSVDError(LPSVD_coefs,Data)
	//A function that estimates the errors on the LPSVD parameters using the Cramer-Rao
	//lower bound (http://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound).
	//This implementation is based on the work of Barkhuijsen et al (http://dx.doi.org/10.1016/0022-2364(86)90446-4)
	WAVE LPSVD_coefs	//Coefficients calculated from the LPSVD algorithm
	WAVE Data			//The data from which the LPSVD coefficients were calculated
	
	//***The first thing to do is to calculated the RMS of the residuals***
	//We reconstruct the model from the parameters
	reconstructSignal(LPSVD_coefs,"recon_"+NameOfWave(Data),numpnts(Data),DimDelta(data, 0),dataReal=0)
	//Note that the spacing doesn't really mattter and we want complex data returned (i.e. dataReal = 0)
	WAVE/C recon = $("recon_"+NameOfWave(Data))
	//Make the residual wave
	Duplicate/C/O recon $("res_"+NameOfWave(Data))
	WAVE/C res = $("res_"+NameOfWave(Data))
	res = data - recon
	
	//Calculate the RMS
	WaveStats/Q/C=4 res
	WAVE M_WaveStats
	Variable RMS = M_WaveStats[5]

	Duplicate/O LPSVD_coefs $("sigma_"+NameOfWave(LPSVD_coefs))
	Wave Errors = $("sigma_"+NameOfWave(LPSVD_coefs))
	Variable i,j	//For loop iterators
	
//	//The code block below calculates the asymptotic Cramer-Rao bounds according to dx.doi.org/10.1109/78.80943
//	//The speed up is not significant enough to actually use this code.
//	If(0)//Change this to test if assumptions hold...
//		For(i=0;i<DimSize(Errors,0);i+=1)
//			//Right now I'm using the results of an analytical theory
//			//dx.doi.org/10.1109/78.80943
//			Errors[i][0] = Sqrt(RMS^2*2*(-LPSVD_coefs[i][2]))
//			Errors[i][1] = Sqrt(RMS^2*4*(-LPSVD_coefs[i][2])^3/(-LPSVD_coefs[i][0])^2)
//			Errors[i][2] = Sqrt(RMS^2*4*(-LPSVD_coefs[i][2])^3/(-LPSVD_coefs[i][0])^2)
//			Errors[i][3] = Sqrt(RMS^2*2*(-LPSVD_coefs[i][2])/(-LPSVD_coefs[i][0])^2)
//		EndFor
//		Return  GetWavesDataFolder(Errors,2)//exit
//	EndIf
	
	//Next we need to generate the Fisher matrix
	Variable size = DimSize(LPSVD_coefs,0)*4
	Make/D/O/N=(size,size) FisherMat = 0
	Variable Chi0, Chi1, Chi2
	Variable Zeta0, Zeta1, Zeta2
	Variable ampi,dampi,freqi,phasei
	Variable ampj,dampj,freqj,phasej
	Variable mySum
	//We'll reuse res for the intermediate calculations
	//This implementation is based on the work of Barkhuijsen et al (http://dx.doi.org/10.1016/0022-2364(86)90446-4)
	For(i=0;i<DimSize(LPSVD_coefs,0);i+=1)
		ampi = LPSVD_coefs[i][%amps]
		freqi = LPSVD_coefs[i][%freqs]
		dampi = LPSVD_coefs[i][%damps]
		phasei = LPSVD_coefs[i][%phase]
		For(j=0;j<DimSize(LPSVD_coefs,0);j+=1)
			ampj = LPSVD_coefs[j][%amps]
			freqj = LPSVD_coefs[j][%freqs]
			dampj = LPSVD_coefs[j][%damps]
			phasej = LPSVD_coefs[j][%phase]
			
			res = exp(p*cmplx(dampi+dampj,2*pi*(freqi-freqj))+cmplx(0,1)*(phasei-phasej))
			WaveStats/Q/C=3 res
			WAVE M_WaveStats
			Chi0 =  M_WaveStats[23][0]
			Zeta0 =  M_WaveStats[23][1]
			
			res = p*exp(p*cmplx(dampi+dampj,2*pi*(freqi-freqj))+cmplx(0,1)*(phasei-phasej))
			WaveStats/Q/C=3 res
			WAVE M_WaveStats
			Chi1 =  M_WaveStats[23][0]
			Zeta1 =  M_WaveStats[23][1]
			
			res = p^2*exp(p*cmplx(dampi+dampj,2*pi*(freqi-freqj))+cmplx(0,1)*(phasei-phasej))
			WaveStats/Q/C=3 res
			WAVE M_WaveStats
			Chi2 =  M_WaveStats[23][0]
			Zeta2 =  M_WaveStats[23][1]
			
			//First Row
			FisherMat[4*i+0][4*j+0] = ampi*ampj*Chi2
			FisherMat[4*i+0][4*j+1] = -ampi*zeta1
			FisherMat[4*i+0][4*j+2] = ampi*ampj*zeta2
			FisherMat[4*i+0][4*j+3] = ampi*ampj*chi1
			//Second Row
			FisherMat[4*i+1][4*j+0] = ampj*zeta1
			FisherMat[4*i+1][4*j+1] = chi0
			FisherMat[4*i+1][4*j+2] = -ampj*chi1
			FisherMat[4*i+1][4*j+3] = ampj*zeta0
			//Third Row
			FisherMat[4*i+2][4*j+0] = -ampi*ampj*zeta2
			FisherMat[4*i+2][4*j+1] = -ampi*chi1
			FisherMat[4*i+2][4*j+2] = ampi*ampj*chi2
			FisherMat[4*i+2][4*j+3] = -ampi*ampj*zeta1
			//Fourth Row
			FisherMat[4*i+3][4*j+0] = ampi*ampj*chi1
			FisherMat[4*i+3][4*j+1] = -ampi*zeta0
			FisherMat[4*i+3][4*j+2] = ampi*ampj*zeta1
			FisherMat[4*i+3][4*j+3] = ampi*ampj*chi0

		EndFor
	EndFor
	MatrixInverse/O FisherMat		//Replace the Fisher matrix with its inverse
	FisherMat*=(2*RMS^2)
	//Fill up the Error wave with the errors.
	
	For(i=0;i<DimSize(Errors,0);i+=1)
		Errors[i][0] = Sqrt((FisherMat[1+i*4][1+i*4]))
		Errors[i][1] = Sqrt((FisherMat[0+i*4][0+i*4]))
		Errors[i][2] = Sqrt((FisherMat[2+i*4][2+i*4]))
		Errors[i][3] = Sqrt((FisherMat[3+i*4][3+i*4]))
	EndFor
	
	Return  GetWavesDataFolder(Errors,2)
End

STATIC Function/S roots(polyCoefs)
	//Uses the method described at http://mathworld.wolfram.com/PolynomialRoots.html
	//The assumed order of PolyCoefs is:
	//a0+...+an-1*x^(n-1)+an*x^n
	Wave/C polyCoefs
	
	Variable N=numpnts(PolyCoefs)
	Make/D/C/O/N=(N-1) normPolyCoefs
	normPolyCoefs=-polyCoefs[p+1]/polyCoefs[0] //Normalize the coefficients to that a0=1
	
	//Construct our matrix
	MatrixOp/C/O A=Identity(N-1,N-2)
	
	Concatenate/O {normPolyCoefs,A}, A2
	
	//Find the eigenvalues, which are the inverse of the roots of the polynomial.
	MatrixEigenV A2
	WAVE eigenValues = W_eigenValues
	
	//Invert eigenvalues to obtain real roots
	Duplicate/C/O eigenValues myRoots
	myRoots=1/myRoots
	
	//Clean up
	KillWaves/Z A, A2, eigenValues
	
	Return GetWavesDataFolder(myRoots,2)// string is full path to wave
End

STATIC Function estimate_model_order(s,N,L)
	//Adapted from from Complex Exponential Analysis by Greg Reynolds
	//http://www.mathworks.com/matlabcentral/fileexchange/12439-complex-exponential-analysis/
	// Use the MDL method as in Lin (1997) to compute the model
	// order for the signal. You must pass the vector of
	// singular values, i.e. the result of svd(T) and 
	// N and L. This method is best explained by Scharf (1992).
	Wave s
	Variable N,L
	
	Make/D/O/N=(L) MDL=0
	Duplicate/FREE/O S lnS
	lnS=ln(S)
	
	Variable i=0
	for(i=0;i<L;i+=1)
		MDL[i] = -N*sum(lnS,i,L)
		MDL[i] += N*(L-i)*ln(sum(s,i,L)/(L-i)) 
		MDL[i] += i*(2*L-i)*ln(N)/2;
	EndFor
	
	WaveStats/Q MDL
	
	Return V_minLoc
End

STATIC Function/S pinv(matrix)
	//A function that calculates the Moore-Penrose inverse
	//Adapted from MathWorks, pinv function
	//For some reason IgorPro doesn't include this functionality...
	//http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
	WAVE matrix
	
	//MatrixSVD Matrix
	MatrixSVD Matrix
	
	WAVE/C S = W_W
	WAVE/C U = M_U
	WAVE/C VT = M_VT
	
	//variable lengthS=numpnts(S)
	
	//Set any singular value less than the tolerance to zero
	Variable epsilon=5e-16		//Machine epsilon for double precision type
	Variable tol = wavemax(S)*epsilon*max(DimSize(matrix,0),DimSize(matrix,1))
	Duplicate/O S S2
	S2 = (s>tol)
	Variable rank = sum(S2)
	S=1/S*S2
	
	redimension/N=(-1 , rank) U //Make U the "right" size
	redimension/N=(rank, -1) VT //Make VT the "right" size
	MatrixOp/C/O M_inverse=VT^h x DiagRC(S,Rank,Rank) x U^h
	Return  GetWavesDataFolder(M_Inverse,2)
End

Function reconstructSignal(LPSVD_coefs,name,length,timeStep,[dataReal,ampcutoff,freqcutoff,dampcutoff])
	//A function that reconstructs the original signal in the time domain and frequency domain
	//from the LPSVD algorithms coefficients, which are passed as LPSVD_coefs
	
	WAVE LPSVD_coefs		//coefficients from the LPSVD algorithm
	String name				//Name of the generated waves
	Variable length			//Length of the time domain signal
	Variable timeStep		//Sampling frequency with which the signal was recorded, in fs
	Variable dataReal		//Should the output time domain data be real?
	Variable ampcutoff		//Cutoff for the amplitudes of the components
	Variable freqcutoff		//Cutoff for the frequency of the components
	Variable dampcutoff		//Cutoff for the damping of the components
	
	If(ParamIsDefault(dataReal))
		dataReal = 1
	EndIf
	
	If(ParamIsDefault(ampcutoff))
		ampCutoff = 0
	EndIf
	
	If(ParamIsDefault(freqcutoff))
		freqCutoff = 0
	EndIf
	
	If(ParamIsDefault(dampcutoff))
		dampCutoff = 0
	EndIf
	
	//Initialize time domain signal
	Make/C/D/O/N=(length) $name = 0
	WAVE/C timeDomain = $name
	
	Make/D/O/N=(2^10+1) $(name+"_FFT") = 0
	WAVE freqDomain = $(name+"_FFT")
	
	SetScale/P x 0, timeStep, "", timeDomain
	setscale/p x, 0, 1/(timestep*2^11*kspeedoflight),"", freqDomain
	
	//Now we can loop through our coefficients and build our wave
	Variable i=0
	
	Variable amp,damp,freq,phase
	
	For(i=0;i<DimSize(LPSVD_coefs,0);i+=1)
		damp = -LPSVD_coefs[i][%damps]/kSpeedOfLight/timestep/pi
		if((LPSVD_coefs[i][%amps])^2 > ampcutoff && damp > dampcutoff)
			amp =  LPSVD_coefs[i][%amps]
			damp = -LPSVD_coefs[i][%damps]/kSpeedOfLight/timestep/pi
			freq = LPSVD_coefs[i][%freqs]/kSpeedOfLight/timestep
			If( abs(Freq) > freqcutoff)		
				freqDomain += amp^2/((x-freq)^2+(damp/2)^2)
			EndIf
			//Keep in mind that LPSVD_coefs were constructed agnostic to the actual sampling
			//frequency so we will reconstruct it in the same way
			damp = LPSVD_coefs[i][%damps]
			phase = LPSVD_coefs[i][%phase]
			freq = LPSVD_coefs[i][%freqs]
			timeDomain += amp*exp(p*cmplx(damp,2*pi*freq)+cmplx(0,1)*phase)
		EndIf
	EndFor
	
	//Scale the frequency domain signal to match up with the FFT of the same signal
	freqDomain=(freqDomain)/(2*pi*kspeedoflight*timestep)^2
	
	If(dataReal)
		Redimension/R timeDomain
	EndIf
	
	Return 0
End

Function/S Cadzow(signal, M, iters,[lFactor,q])
	//Remove noise using the Cadzow composite property mapping method.
	//See Cadzow, J. A. IEEE Transactions on Acoustics, Speech and Signal Processing 1988, 36, 49 Ð62.
	//Adapted from from Complex Exponential Analysis by Greg Reynolds
	//http://www.mathworks.com/matlabcentral/fileexchange/12439-complex-exponential-analysis/
	
	Wave signal		//The signal to be filtered
	Variable M		//The expected number of signals (2 times the number of damped sinusoids
	Variable iters	//Number of iterations to be performed
	
	Variable lFactor	//User selectable factorization of the Hankel Matrix
	Variable q		//Verbose or not
	
	If(ParamIsDefault(lFactor))
		lFactor = 1/2
	EndIf
	
	If(ParamIsDefault(q))
		q=0
	Else
		q=1
	EndIf
	
	//We want this function to be data folder aware
	//We'll do all our calculations in a specific data folder and then kill that folder at the end
	String savDF= GetDataFolder(1)	// Save current DF for restore.
	if( DataFolderExists("root:Cadzow_Data") )
		SetDataFolder root:Cadzow_Data
	else 
		NewDataFolder/O/S root:Cadzow_Data	// Our stuff goes in here.
	endif
	
	//Timing
	Variable timerRef=startMSTimer
	
	Variable N = numpnts(signal);
	Variable L = floor(N*lFactor);
	 
	// T is the prediction matrix before filtering.
	Wave/C T = $Hankel(signal, N-L, L+1)
	T = conj(T)
	
	If(M>(N-L))
		M = N-L
		Print "M too large M set to: " + num2str(M)
	EndIf
	
	Variable i = 0
	Variable tol = 0
	Variable r = 0
	
	Print "Beginning Cadzow filtration, press ESC to abort, press CMD to check status."
	
	for(i=0;i<iters;i+=1)
		 
		// decompose T
		//MatrixSVD Matrix
		MatrixSVD T
	
		WAVE/C S = W_W
		WAVE/C U = M_U
		WAVE/C VT = M_VT
		   
		// check current rank
		tol = L*5e-16
		Duplicate/O S, S2
		S2 = (s>tol)
		r = sum(S2)
		
		If(q || (GetKeyState(0) & 1))
			printf "Cadzow iteration %d (rank is %d, target is %d).\r", i, r,M
		EndIf
		
		If(r <= M)
			Printf "Successful completion: "
			break
		ElseIf( r > M )
			//Filter the hankel matrix
			S = S*(p < M)
			Redimension/N=(-1,M) U
			Redimension/N=(M,-1) VT
			MatrixOp/C/O T = U x DiagRC(S,M,M) x VT
			// average to restore Hankel structure
			Wave recon_signal = $unHankelAvg(T)
			WAVE/C T = $Hankel(recon_signal,N-L,L+1)
		EndIf
		if (GetKeyState(0) & 32)	// Is Escape key pressed now?
			Printf "User abort: "
			Break
		EndIf
	EndFor
	 
	// need to extract data from matrix Tr
	T = conj(T)
	
	//Move the results to the original data folder
	Duplicate/O $unHankelAvg(T), $(savDF+nameOfWave(signal)+"_cad")
	WAVE nSignal = $(savDF+nameOfWave(signal)+"_cad")
	SetDataFolder savDF	// Restore current DF.
	
	//Clean up
	KillDataFolder root:Cadzow_Data
	
	//Scale the new signal appropriately
	CopyScales/P signal, nSignal
	
	//If the original signal was real, make the new signal real as well
	If((WaveType(signal) & 2^0) == 0)
		Redimension/R nSignal
	EndIf
	
	//Finish up the timing
	Variable microseconds = stopMSTimer(timerRef)
	Variable minutes = floor(microseconds/(60e6))
	Variable seconds = microseconds/(1e6)-minutes*60 
	
	if(!q)
		printf "Final rank is %d, target is %d, ", r,M
	EndIf
	
	Printf "%d iterations took ", i 
	If(minutes > 1)
		Printf "%g minutes and ",minutes
	ElseIf(minutes > 0)
		Printf "1 minute and "
	EndIf
	Printf "%g seconds, for %g sec/iter.\r",seconds,microseconds/(1e6)/i
	
	Return  GetWavesDataFolder($(nameOfWave(signal)+"_cad"),2)
End

 STATIC Function/S Hankel(signal, numRows, numCols)
	//A function that forms the Hankel Matrix from a signal
	//http://en.wikipedia.org/wiki/Hankel_matrix
	Wave signal
	Variable numRows
	Variable numCols
	
	Make/C/D/O/N=(numRows,numCols) myHankel=Signal[p+q]
	
	Return  GetWavesDataFolder(myHankel,2)
End

STATIC Function/S unHankelAvg(Hankel)
	//A function that takes a Hankel matrix and returns the original signal
	//that it was formed from by averaging along the anti-diagonals
	Wave/C Hankel		//The matrix to be inverted
	
	Variable numRows = DimSize(Hankel,0)
	Variable numCols = DimSize(Hankel,1)
	
	//Make the signal to be returned, make sure to set to zero!
	Make/C/D/O/N=(numRows+numCols-1) mySignal=0
	
	variable i=0,j=0
	Duplicate/C/O mySignal myNorm //Make the normalizing wave
	For(i=0;i<numRows;i+=1)
		For(j=0;j<numCols;j+=1)
			//Build up the signal and the norm
			mySignal[i+j]+=Hankel[i][j]
			myNorm[i+j] += cmplx(1,0)
		EndFor
	EndFor
	mySignal=mySignal/myNorm
	Return  GetWavesDataFolder(mySignal,2)
End

Function DrawPeaks(LPSVD_coefs,timestep,[freqcutoff,ampcutoff,dampcutoff])
//A quick procedure to draw lines indicating the peak values for the reconstructed frequency domain signal.
	WAVE LPSVD_coefs				//The LPSVD_coefs that you'd like to draw
	Variable timestep				//The time step of the original signal (for scaling the LPSVD coefs)
	Variable freqcutoff,ampcutoff,dampcutoff	//Frequency and amplitude cutoffs
	
	If(ParamIsDefault(ampcutoff))
		ampcutoff=0
	Endif
	
	If(ParamIsDefault(freqcutoff))
		freqcutoff=0
	Endif
	
	If(ParamIsDefault(dampcutoff))
		dampCutoff = 0
	EndIf
	//make sure our draw environment is good to go
	SetDrawLayer/K ProgBack
	SetDrawEnv linefgc= (34952,34952,34952),dash= 0,xcoord= bottom,ycoord= prel
	SetDrawEnv textxjust= 1,textyjust= 2,textrot= 90, save
	Variable i=0,Freq = 0,damp=0
	
	//Drawing the lines and the tags
	For(i=0;i<DimSize(LPSVD_coefs,0);i+=1)
		damp = -LPSVD_coefs[i][%damps]/kSpeedOfLight/timestep/pi
		if((LPSVD_coefs[i][%amps])^2 > ampcutoff && damp >dampcutoff)
			freq=LPSVD_coefs[i][%freqs]/kSpeedOflight/timestep
			if(freq>freqcutoff)
				If(round(freq) < 100)
					DrawLine freq, 1, freq, 0.25
				ElseIf(round(freq) < 1000)
					DrawLine freq, 1, freq, 0.3
				Else
					DrawLine freq, 1, freq, 0.35
				EndIf
				DrawText freq, 0, num2istr(freq)+" cm\\S-1\\M"
			EndIf
		EndIf
	EndFor
	//Change the Draw layer to front
	SetDrawLayer UserFront
End

Function OptimizeLPSVDCoefs(Data,LPSVD_coefs,[ampcutoff,freqcutoff,dampcutoff,holdfreqphase])
	Wave Data									//The original data
	Wave LPSVD_coefs							//Parameters to optimize
	Variable ampcutoff, freqcutoff,dampcutoff	//Cutoff parameters to remove spurious values
	Variable holdfreqphase						//hold the phases and frequencies constant during the fit
	
	If(ParamIsDefault(ampcutoff))
		ampcutoff=0
	EndIf
	
	If(ParamIsDefault(freqcutoff))
		freqcutoff=0
	EndIf
	
	If(ParamIsDefault(dampcutoff))
		dampcutoff=0
	EndIf
	
	If(ParamIsDefault(holdfreqphase))
		holdfreqphase=0
	EndIf
	
	//Make a copy of the LPSVD_coefs, we'll use this wave later
	//to repack to optimized variables
	Duplicate/O LPSVD_coefs $("opt"+NameOfWave(LPSVD_coefs))
	WAVE newLPSVD_coefs = $("opt"+NameOfWave(LPSVD_coefs))
	
	//Make a copy of Data and remove the scaling from the copy.
	Duplicate/O Data $("fit_"+nameofwave(Data))
	WAVE newData = $("fit_"+nameofwave(Data))
	SetScale/P x,0,1,"", newData
	
	Variable numComponents = DimSize(LPSVD_coefs,0)
	variable i = 0
	String removedComponents = ""
	For(i=numComponents;i>0;i-=1)
		If((newLPSVD_coefs[i-1][%amps])^2<ampcutoff || (-LPSVD_coefs[i-1][%damps]/kSpeedOfLight/dimdelta(data,0)/pi) < dampcutoff || abs(newLPSVD_coefs[i-1][%freqs])<freqcutoff)
			removedComponents += num2istr(abs(newLPSVD_coefs[i-1][%freqs])/kSpeedOfLight/DimDelta(data,0)) +", "
			DeletePoints (i-1),1, newLPSVD_coefs
			numComponents-=1
		Endif
	EndFor
	
	If(strlen(removedComponents))
		Print "The following frequency components were removed: " + removedComponents
	EndIf
	
	//unpack LPSVD_coefs into a regular coefficient wave
	//Make use of the fact that only half of the coefficients are necessary
	//Also, set any frequency below some tolerance to zero and hold it there
	Variable numCoefs =  ceil(numComponents/2)
	Make/D/O/N=(numCoefs*4) myCoefs
	String HoldStr = ""
	For(i=0;i<numCoefs;i+=1)
		myCoefs[4*i] = 2*LPSVD_coefs[i][%amps]
		myCoefs[4*i+1] = LPSVD_coefs[i][%damps]
		If(abs(LPSVD_coefs[i][%freqs])<1e-14)
			myCoefs[4*i+2] = 0
			myCoefs[4*i+3] = 0
		Else
			myCoefs[4*i+2] = LPSVD_coefs[i][%freqs]
			myCoefs[4*i+3] = LPSVD_coefs[i][%phase]
		EndIf
		If(holdfreqphase)
			HoldStr+="0011"
		Else
			HoldStr+="0000"
		EndIf
	EndFor
	
	//If there are an odd number of components the middle one is zero frequency
	If(numCoefs-floor(DimSize(LPSVD_coefs,0)/2))
		myCoefs[4*(numCoefs-1)] /= 2
	EndIf
	Variable V_FitNumIters
	Variable V_FitMaxIters=200
	//do the optimization (we're using funcfit, so we're minimizing the chi^2)
	FuncFit/H=holdstr/N/W=2/Q decayingSinusoids, myCoefs, newData
	
	Print "Number of interations: "+num2str(V_FitNumIters)
	//Well use the newData wave to hold the fit, why not?
	newData = decayingSinusoids(myCoefs,p)
	
	//Return the scaling
	CopyScales/P Data newData
	
	//Repack
	For(i=0;i<numCoefs;i+=1)
		newLPSVD_coefs[i][%amps] = myCoefs[4*i]/2
		newLPSVD_coefs[i][%damps] = myCoefs[4*i+1]
		newLPSVD_coefs[i][%freqs] = myCoefs[4*i+2]
		newLPSVD_coefs[i][%phase] = myCoefs[4*i+3]
		
		newLPSVD_coefs[2*numCoefs-i-1][%amps] = myCoefs[4*i]/2
		newLPSVD_coefs[2*numCoefs-i-1][%damps] = myCoefs[4*i+1]
		newLPSVD_coefs[2*numCoefs-i-1][%freqs] = -myCoefs[4*i+2]
		newLPSVD_coefs[2*numCoefs-i-1][%phase] = -myCoefs[4*i+3]
	EndFor
End

Function decayingSinusoids(w,t)
	//w[i] = amp
	//w[i+1] = damp
	//w[i+2] = freq
	//w[i+3] = phase
	Wave w
	Variable t
	
	Variable val=0
	Variable i=0
	Variable npts = numpnts(w)
	For(i=0;i<npts;i+=4)
		val += w[i]*exp(t*w[i+1])*Cos(2*pi*w[i+2]*t+w[i+3])
	EndFor
	
	Return val
End

Function PrintLPSVDCoefs(LPSVD_coefs,sigma_LPSVD_coefs,timestep,offset)
	WAVE LPSVD_coefs
	WAVE sigma_LPSVD_coefs
	Variable Timestep
	Variable offset
	
	Variable i=0,corrPhase
	Printf "  Amp\t\t+/-\t\tWidth\t\t+/-\t\tDamp\t\t+/-\t\tFreq\t\t+/-\t\tPhase\t\t+/-\t\tCorrPhase\r"
	For(i=0;i<DimSize(LPSVD_coefs,0);i+=1)
		Printf "%6.3f\t", LPSVD_coefs[i][%amps]
		Printf "%6.3f\t\t", sigma_LPSVD_coefs[i][%amps]
		
		Printf "%6.2f\t\t", -LPSVD_coefs[i][%damps]/kSpeedOfLight/timestep/pi
		Printf "%6.2f\t\t", sigma_LPSVD_coefs[i][%damps]/kSpeedOfLight/timestep/pi
		
		Printf "%6.2f\t\t", -1/(LPSVD_coefs[i][%damps]/timestep)/1000
		Printf "%6.2f\t\t", sqrt((1/(LPSVD_coefs[i][%damps]/timestep)/1000)^2*sigma_LPSVD_coefs[i]^2)
		
		Printf "%6d\t\t", LPSVD_coefs[i][%freqs]/kSpeedOfLight/timestep
		Printf "%6d\t\t", sigma_LPSVD_coefs[i][%freqs]/kSpeedOfLight/timestep
		
		Printf "%6d\t\t", LPSVD_coefs[i][%phase]/pi*180
		Printf "%6d\t\t", sigma_LPSVD_coefs[i][%phase]/pi*180
		corrPhase=LPSVD_coefs[i][%phase]-2*pi*LPSVD_coefs[i][%freqs]*offset/timestep
		Printf "%6d\r", (corrPhase-2*pi*floor(corrPhase/2/pi))/pi*180
	EndFor
End