#pragma rtGlobals=3		// Use modern global access method.

Constant kSpeedOfLight=2.99792458e-5 //(cm/fs)

//LPSVD was developed by Tufts and Kumaresan (Tufts, D.; Kumaresan, R. IEEE Transactions on Acoustics,
//Speech and Signal Processing 1982, 30, 671 Ð 675.) as a method of harmonic inversion, i.e. decomposing
//a time signal into a linear combination of (decaying) sinusoids.
//
//A great reference that is easy to read for the non-EECS user is:
// Barkhuijsen, H.; De Beer, R.; BovŽe, W. M. M. .; Van Ormondt, D. J. Magn. Reson. (1969) 1985, 61, 465Ð481.
//
//This particular implementation was adapted, in part, from matNMR (http://matnmr.sourceforge.net/)
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
	MatrixOp/C/O LPCoefs= (-1)*VT^h x DiagRC(S,M,M) x U^h x h
	
	//Error check: are there any NaNs or INFs in LPCoefs?
	Variable V_numNans, V_numINFs
	WaveStats/Q LPCoefs
	If(V_numNans!=0 || V_numINFs!=0)
		Print "There has been an error generating the prediction-error filter polynomial"
		SetDataFolder savDF	// Restore current DF.
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
		return -1
	EndIf
	
	//Before continuing lets sort everything by the frequencies
	Make/D/O/N=(numpnts(usedRoots)) sortKeys
	sortKeys = imag(usedRoots)
	Sort sortKeys usedRoots
	
	//Lets make a matrix with dimension labels to store all our parameters
	Make/D/O/N=(numpnts(usedRoots),5) LPSVD_coefs = 0
	SetDimLabel 1, 0, amps, LPSVD_coefs
	SetDimLabel 1, 1, amps2, LPSVD_coefs
	SetDimLabel 1, 2, freqs, LPSVD_coefs
	SetDimLabel 1, 3, damps, LPSVD_coefs
	SetDimLabel 1, 4, phase, LPSVD_coefs
	
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
	LPSVD_coefs[][%amps2]=magsqr(cAmps[p])
	LPSVD_coefs[][%phase]=Imag(r2polar(cAmps[p]))
	
	//Move the results to the original data folder
	Duplicate/O LPSVD_coefs, $(savDF+"LPSVD_Coefs")
	
	SetDataFolder savDF	// Restore current DF.
	
	//Clean up
	KillDataFolder root:LPSVD_Data
	
	//Output the calculation time
	Print "LPSVD time: "+num2str(stopMSTimer(timerRef)*1e-6)+" seconds"
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

Function reconstructSignal(LPSVD_coefs,name,length,timeStep,[ampcutoff,freqcutoff])
	//A function that reconstructs the original signal in the time domain and frequency domain
	//from the LPSVD algorithms coefficients, which are passed as LPSVD_coefs
	
	WAVE LPSVD_coefs		//coefficients from the LPSVD algorithm
	String name				//Name of the generated waves
	Variable length			//Length of the time domain signal
	Variable timeStep		//Sampling frequency with which the signal was recorded, in fs
	Variable ampcutoff		//Cutoff for the amplitudes of the components
	Variable freqcutoff		//Cutoff for the frequency of the components
	
	If(ParamIsDefault(ampcutoff))
		ampCutoff = 0
	EndIf
	
	If(ParamIsDefault(freqcutoff))
		freqCutoff = 0
	EndIf
	
	//Initialize time domain signal
	Make/D/O/N=(length) $name = 0
	WAVE timeDomain = $name
	
	Make/D/O/N=(2^10+1) $(name+"_FFT") = 0
	WAVE freqDomain = $(name+"_FFT")
	
	SetScale/P x 0, timeStep, "", timeDomain
	setscale/p x, 0, 1/(timestep*2^11*kspeedoflight),"", freqDomain
	
	//Now we can loop through our coefficients and build our wave
	Variable i=0
	
	Variable amp2,amp,damp,freq,phase
	
	For(i=0;i<DimSize(LPSVD_coefs,0);i+=1)
		if(LPSVD_coefs[i][%amps2] > ampcutoff)
			amp =  LPSVD_coefs[i][%amps]
			amp2 = LPSVD_coefs[i][%amps2]
			damp = -LPSVD_coefs[i][%damps]/kSpeedOfLight/timestep/pi
			freq = LPSVD_coefs[i][%freqs]/kSpeedOfLight/timestep
			If( abs(Freq) > freqcutoff)		
				freqDomain += amp2/((x-freq)^2+(damp/2)^2)
			EndIf
			//Keep in mind that LPSVD_coefs were constructed agnostic to the actual sampling
			//frequency so we will reconstruct it in the same way
			damp = LPSVD_coefs[i][%damps]
			phase = LPSVD_coefs[i][%phase]
			freq = LPSVD_coefs[i][%freqs]
			timeDomain += amp*exp(p*damp)*Cos(2*pi*freq*p+phase)
		EndIf
	EndFor
	
	//Scale the frequency domain signal to match up with the FFT of the same signal
	freqDomain=(freqDomain)/(2*pi*kspeedoflight*timestep)^2
	
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

Function/S unHankelAvg(Hankel)
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

Function DrawPeaks(LPSVD_coefs,timestep,[freqcutoff,ampcutoff])
//A quick procedure to draw lines indicating the peak values for the reconstructed frequency domain signal.
	WAVE LPSVD_coefs				//The LPSVD_coefs that you'd like to draw
	Variable timestep				//The time step of the original signal (for scaling the LPSVD coefs)
	Variable freqcutoff,ampcutoff	//Frequency and amplitude cutoffs
	
	If(ParamIsDefault(ampcutoff))
		ampcutoff=0
	Endif
	
	If(ParamIsDefault(freqcutoff))
		freqcutoff=0
	Endif
	
	//make sure our draw environment is good to go
	SetDrawLayer/K ProgBack
	SetDrawEnv linefgc= (0,0,0),dash=1,xcoord= bottom,ycoord= prel
	SetDrawEnv textxjust= 1,textyjust= 2,textrot= 90, save
	Variable i=0,Freq = 0
	
	//Drawing the lines and the tags
	For(i=0;i<DimSize(LPSVD_coefs,0);i+=1)
		if(LPSVD_coefs[i][%amps2] > ampcutoff)
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