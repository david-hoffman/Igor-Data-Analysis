#pragma rtGlobals=3		// Use modern global access method.
#include ":FitFuncs"

STATIC Constant kSpeedOfLight=2.99792458e-5 //(cm/fs)

Function FancifyISRS() : GraphStyle
	//A quick macro to format a plot in a good style
	//box the figure, add minor ticks to the bottom, and scale the left axis to be mOD
	ModifyGraph mirror=2,minor(bottom)=1,tick(left)=3,noLabel(left)=1,lblMargin(left)=10
	ModifyGraph rgb=(0,0,0)
	//Label the axes
	Label left "Spectral Intensity (a.u.)"
	Label bottom "Raman Shift (cm\\S-1\\M)"
	//Autoscale only visible data
	SetAxis/A=2 left
	
	String tnl = TraceNameList("",";",1)
	String tn= StringFromList(0,tnl, ";")
	
	Variable  offset=-Wavemax($tn)/10
	
	ModifyGraph prescaleExp(left)=-floor(log(wavemax($tn)))
	
	Variable i=1
	Do
		tn= StringFromList(i,tnl, ";")
		if (strlen(tn) ==0)
			break
		endif
		ModifyGraph lstyle[i]=3,offset[i]={0,offset}
		i+=1
	While(1)
End

Function/S CorrectMatrix(matrix,t0)
	Wave Matrix
	Wave t0
	variable i=0,length = DimSize(Matrix,0)//Length is the number of Pixels (1340 for the PIXIS)
	
	Duplicate/O Matrix $("Corr_"+nameofwave(Matrix))//This wave is as long as there are timepoints
	
	Wave myCorr = $("Corr_"+nameofwave(Matrix))
	
	Make/D/O/N=(DimSize(Matrix,1)) tempWave//This wave is as long as there are timepoints
	
	SetScale/P x DimOffset(Matrix,1),DimDelta(Matrix,1),"", tempWave//Set the scale so
	//that we can shift it around and perform a linear interpolation easily.
	
	//PadWave(tempWave,1000)
	
	for(i=0;i<length;i+=1)
		If(numtype(t0[i])==0)//not inf or Nan
			tempWave =Matrix[i][p]
			myCorr[i][]= tempWave(y+t0[i])
		Else
			myCorr[i][]=0
		EndIf
	EndFor
	Return GetWavesDataFolder(myCorr,2)// string is full path to wave
End

Function GenFFT(res,num)
	WAVE res
	Variable num
	variable i=0
	
	Variable Length = nextpow2(numpnts(res))
	Duplicate/O res res_TEMP
	
	FFT/Pad={length}/OUT=4/DEST=Res_TEMP_FFT res_TEMP
	
	Duplicate/O Res_TEMP_FFT $(nameofwave(res)+"_FFTMat")
	
	For(i=1;i<num;i+=1)
	
		res_TEMP=(p>i)*Res
		
		FFT/Pad={length}/OUT=4/DEST=Res_TEMP_FFT res_TEMP
		
		Concatenate {Res_TEMP_FFT}, $(nameofwave(res)+"_FFTMat")
		
	EndFor
	
	SetScale/p x, 0, 1/(length*deltax(res)*kspeedoflight), $(nameofwave(res)+"_FFTMat")
	SetScale/p y, leftx(res), deltax(res), $(nameofwave(res)+"_FFTMat")
	
	PauseUpdate
	NewWaterfall/N=GenFFT/W=(0,0,288,288) $(nameofwave(res)+"_FFTMat")
	ModifyWaterfall angle= 60
	ModifyWaterfall axlen= 0.6
	TextBox/C/N=text0/X=45.00/Y=-23.00/O=60/F=0/A=MC "End of Mask (fs)"
	//ModifyGraph width=288,height=288
	Label left "FFT Power (a.u.)"
	Label bottom "Raman Shift (cm\\S-1\\M\\u#2)"
	duplicate/O $(nameofwave(res)+"_FFTMat") waterfallcolor;waterfallcolor=q
	ModifyGraph zColor($(nameofwave(res)+"_FFTMat"))={waterfallcolor,*,*,BlueRedGreen256,0}
	ModifyGraph mode=7,hbFill=5
	ModifyGraph minor(bottom)=1
	ResumeUpdate
	KillWaves Res_TEMP_FFT, res_TEMP
End	

Function Calibrate2(pixel,wavelength,ranwave)
	WAVE pixel, wavelength,ranwave
	variable length=numpnts(ranwave)
	// This will perform a linear fit to wavelength vs Pixels for a known standard light source
	// For use in calibrating the axis for transient absorption
	
	Display/N=Calibration wavelength vs pixel as "Wavelength Calibration"
	ModifyGraph mode=3,marker=8
	CurveFit/Q line wavelength /X=pixel /D
	Wave w_coef, w_sigma
 	TextBox/C/N=text0/X=-2.50/Y=36.94  "Wavelength Calibration\rwavelength = a*pixel+b\ra= "+num2str(W_coef[1])+"±"+num2str(W_sigma[1])+"\rb = "+num2str(W_coef[0])+"±"+num2str(W_sigma[0])
	TextBox/C/N=text1/X=27.81/Y=4.46 "R\\S2\\M = " +num2str(V_r2)
	Label left "Wavelength (nm)"
	Label bottom "Pixel"
	ModifyGraph rgb($("fit_"+NameOfWave(wavelength)))=(0,0,65535)
	ModifyGraph mirror=2
	Print "offset = "+num2str(W_coef[0])
	Print "deltax = "+num2str(W_coef[1])
	//Store these as global variables for easy access in the future
	Variable/G wlDeltax = W_coef[1]
	Variable/G wlLeftx = W_coef[0]
	make/D/O/N=(length) nmshiftx=W_coef[0]+W_coef[1]*x
	Make/D/O/N=(length+1) nmshiftxForImages=W_coef[0]+W_coef[1]*x//creates a shiftx and a nmshiftx
End

Function XPMMax(Matrix,[startpnt,endpnt])
	//Same as XCMax except it uses the FindPeak function to find the center and width at each pixel
	//and then fits the data to the XPM model function (cross phase modulation)
	WAVE Matrix
	
	Variable startpnt, endpnt
	Variable i=0, length=DimSize(Matrix,0)//Length is the number of Pixels (1340 for the PIXIS)
	Variable V_FitError=1, t0=0, A=0, C=0,w=0, triedOnce=0
	Variable V_min, V_minloc, V_max, V_maxloc
	
	If(ParamIsDefault(endpnt))
		endpnt = length
	EndIf
	
	//Make waves to hold the important fit parameters as a function of pixel
	If(!WaveExists($(NameOfWave(Matrix)+"_Max")))
		Make/D/O/N=(length) $(NameOfWave(Matrix)+"_Max")=NaN, $(NameOfWave(Matrix)+"_Width")=NaN
		Make/D/O/N=(length) $(NameOfWave(Matrix)+"_MaxS")=NaN, $(NameOfWave(Matrix)+"_WidthS")=NaN
	Endif
	
	//Make the waves available to the function
	WAVE tempMax = $(NameOfWave(Matrix)+"_Max")
	WAVE tempWidth = $(NameOfWave(Matrix)+"_Width")
	WAVE tempMaxS = $(NameOfWave(Matrix)+"_MaxS")
	WAVE tempWidthS = $(NameOfWave(Matrix)+"_WidthS")
	
	Make/D/O/N=(DimSize(Matrix,1)) tempFit//This wave is as long as there are timepoints
	
	SetScale/P x DimOffset(Matrix,1),DimDelta(Matrix,1),"", tempFit//Set the scale so the fit parameters
	//mean something useful
	
	//Make my fit coefficients, I don't expect to need these in the future so I'll use a free wave
	//to keep my workspace clean
	Make/D/O/N=6 XPMCoefs
	
	//Start my timer!
	Variable timerRefNum = startMSTimer
	
	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
	PauseUpdate
	for(i=startpnt;i<endpnt;i+=1)
		
		If(mod(i,20)==0)
			//Progress indicator
			PrintF "*"
		EndIf
		
		tempFit = Matrix[i][p]//Assign column i to tempFit
		
		If(V_FitError)
		//If the last fit was faulty, try estimating the parameters again
			Wavestats/Q/Z tempFit
		
			//Based on the xpm model the width can be approximated by the separation of the two troughs
			//for normal phase XPM or two peaks for inverse XPM, which is 2 times the separation of peak
			//trough
		
			//Honestly this parameter estimation section could be greatly improved
			w=2*abs(V_minloc-V_maxloc)
		
			If(V_minloc<V_maxloc)
				//this is the regular situation
				A = V_max
				C = E*V_min/V_max-1
			Else
				A = V_min
				C = E*V_max/V_min-1
			EndIf
			t0 = (V_maxloc+V_minloc)/2
			
			//set my guess coefficients
			XPMCoefs = {w,t0,0,A,0,C}
		EndIf
		//DoUpdate
		V_FitError=0
		FuncFit/Q/N=1/W=2  XPM XPMCoefs tempFit /D /R /A=0
		
		//These declarations insures that the two waves exist, really only important if this macro is run
		//before W_Coef and W_Sigma are created.
		WAVE w_sigma
		//print w_coef
		If(!V_fiterror && (w_sigma[1]<10) &&  (abs(XPMCoefs[0]) > 30))
			//Fill the results waves
			tempMax[i]=XPMCoefs[1]
			tempWidth[i]=abs(XPMCoefs[0])
			tempMaxS[i]=w_sigma[1]
			tempWidthS[i]=w_sigma[0]
			triedOnce=0
		Else
			//With NaN's if that makes the most sense
			tempMax[i]=NaN
			tempWidth[i]=NaN
			tempMaxS[i]=NaN
			tempWidthS[i]=NaN
			If(!triedOnce)
				//I haven't tried once yet, so I'll go back again
				i-=1
				v_fiterror=1
				//Now I've tried once so I'll not do this again
				triedOnce=1
			Else
				triedOnce=0
			EndIF
		EndIf
	EndFor
	ResumeUpdate
	Print "/-"
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	Return 0
	
End

Function XPM(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) = y0+(A+B*(t-t0)+C*(t-t0)^2)*Exp(-((t-t0)/w)^2)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = w
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = A
	//CurveFitDialog/ w[4] = B
	//CurveFitDialog/ w[5] = C

	return w[2]+w[3]*(1+w[4]*(t-w[1])/Abs(w[0])+w[5]*((t-w[1])/w[0])^2)*Exp(-((t-w[1])/w[0])^2)
End

Function MatrixFFT(matrix,[startpnt,endpnt])
	//A function, similar to XCintegrate, except that it FFT's each column before integrating
	Wave matrix
	Variable startpnt,endpnt
	
	Variable i=0, length=DimSize(Matrix,0)//Length is the number of Pixels (1340 for the PIXIS)
	
	//Caculate the appropriate padding in case the number of points isn't a power of 2
	Variable padding = nextPow2(DimSize(Matrix,1))
	
	//Deal with the optional parameters
	If(ParamIsDefault(startpnt))
		startpnt=0
	EndIf

	If(ParamIsDefault(endpnt))
		endpnt=length
	EndIf
	
	//Make my waves
	Make/D/O/N=(DimSize(Matrix,1)) tempFFT//Wave to FFT
	Make/D/O/N=(Padding/2+1) $(nameofwave(matrix)+"_FFT_Sum")=0//Wave for the sum
	WAVE FFT_Sum = $(nameofwave(matrix)+"_FFT_sum")
	
	Make/D/O/N=(length,padding/2+1) $(nameofwave(matrix)+"_FFT")=0//Matrix FFT wave
	WAVE matrix_FFT = $(nameofwave(matrix)+"_FFT")
	
	Variable timerRefNum = startMSTimer
	
	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
	
	for(i=startpnt;i<endpnt;i+=1)
		If(mod(i,20)==0)
			PrintF "*"
		EndIf
		DoUpdate
		tempFFT = Matrix[i][p]//Assign column i to tempFit
		//Take the FFT and return the mag squared into tempFFT_Dest
		FFT/OUT=4/Winf=v/PAD={padding}/DEST=tempFFT_dest tempFFT
		//Add this to the ongoing sum
		FFT_sum+=tempFFT_dest
		//And assign this to the output matrix
		Matrix_FFT[i][]=tempFFT_dest[q]
	EndFor
	Print "/-"
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	//Figure out my scaling
	Variable k = 1/(DimDelta(Matrix,1)*padding*2.99792458*10^-5)
	//Set the scales
	SetScale/P y 0,k,"", Matrix_FFT
	SetScale/P x 0,k,"", FFT_sum
	
	Return 0
End

Function MatrixFFTwMask(matrix,t0,cutoff,[startpnt,endpnt])
	//A function, similar to XCintegrate, except that it FFT's each column before integrating
	Wave matrix,t0
	Variable cutoff
	Variable startpnt,endpnt
	
	Variable i=0, length=DimSize(Matrix,0)//Length is the number of Pixels (1340 for the PIXIS)
	
	//Caculate the appropriate padding in case the number of points isn't a power of 2
	Variable padding = nextPow2(DimSize(Matrix,1))
	If(padding<2^10)
		padding = 2^10
	EndIf
	
	//Deal with the optional parameters
	If(ParamIsDefault(startpnt))
		startpnt=0
	EndIf

	If(ParamIsDefault(endpnt))
		endpnt=length
	EndIf
	
	//Make my waves
	Make/D/O/N=(DimSize(Matrix,1)) tempFFT//Wave to FFT
	SetScale/P x, DimOffset(Matrix,1), DimDelta(Matrix,1),"", TempFFT
	
	Make/D/O/N=(Padding/2+1) $(nameofwave(matrix)+"_FFT_Sum")=0//Wave for the sum
	WAVE FFT_Sum = $(nameofwave(matrix)+"_FFT_sum")
	
	Make/D/O/N=(length,padding/2+1) $(nameofwave(matrix)+"_FFT")=0//Matrix FFT wave
	WAVE matrix_FFT = $(nameofwave(matrix)+"_FFT")
	SetScale/P x, DimOffset(Matrix,0), DimDelta(Matrix,0),"", matrix_FFT
	
	Variable timerRefNum = startMSTimer
	
	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
	
	for(i=startpnt;i<endpnt;i+=1)
		If(mod(i,20)==0)
			PrintF "*"
		EndIf
		DoUpdate
		tempFFT = (x>(cutoff+t0[i]))*Matrix[i][p]//Assign column i to tempFit
		//Take the FFT and return the mag squared into tempFFT_Dest
		FFT/OUT=4/PAD={padding}/DEST=tempFFT_dest tempFFT
		//Add this to the ongoing sum
		FFT_sum+=tempFFT_dest
		//And assign this to the output matrix
		Matrix_FFT[i][]=tempFFT_dest[q]
	EndFor
	Print "/-"
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	//Figure out my scaling
	Variable k = 1/(DimDelta(Matrix,1)*padding*2.99792458*10^-5)
	//Set the scales
	SetScale/P y 0,k,"", Matrix_FFT
	SetScale/P x 0,k,"", FFT_sum
	
	Return 0
End

Function XPM2(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(t) = y0+A*Exp(-((t-t0)/w)^2)*Sin((t-t0)*f+d)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = w
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = d
	//CurveFitDialog/ w[4] = l
	//CurveFitDialog/ w[5] = A

	return w[2]+w[5]*Exp(-((t-w[1])/w[0])^2)*Sin(2*pi*((t-w[1])/w[4]+w[3]))
End

Function scaleSolvent(s) : FitFunc
	//A simple structure fit function that takes a spectrum
	//shifts it and scales it to fit the data
	//useful for removing the solvent artifact
	Struct scaledSolventStruct &s
	
	Return s.coefw[0]*s.solvent(s.x-s.coefw[1])
End

Structure scaledSolventStruct
	//the structure associated with the above
	//structure fit function
	WAVE coefw
	Variable x
	WAVE solvent
EndStructure

Function FitTAwSol(Data, solvent, decay_params,solScale,soloffset,[holdoffset,q])
	//A function for fitting a transient absorption decay to a sum of exponentials convoluted with
	//a gaussian instrument response function (see FitFuncs.ipf for more details) and a scaled,
	//shifted, solvent TA.
	WAVE Data					//The TA of the sample
	WAVE Solvent				//The TA of the solvent
	WAVE decay_params		//Your guesses for the exponential parameters
	Variable solScale,solOffset	//Your guess for the solvent scale, shift parameters
	Variable holdOffset			//1 means that the shift parameter will be held constant
	Variable q					//quiet, no printout

	//Default behaviour is to allow shift parameter to vary
	If(ParamIsDefault(holdoffset))
		holdoffset = 0
	EndIf
	
	If(ParamIsDefault(q))
		q = 0
	EndIf
	
	//Print out the initial decay paramters so that the user knows what they are
	If(!q)
		Print decay_params
	EndIf
	
	//Make a temp wave holding the solvent trace
	Duplicate/O Solvent tempSol
	
	//Pad the temp wave so that we won't get an index out of range error when shifting
	PadWave(tempSol,1000)
	
	//Declare my structure for fitting
	Struct scaledSolventStruct solStruct
	
	//Initialize my fitting parameters for the solvent
	Make/D/O/N=2 scaleFactors = {solScale,soloffset}
	
	//Set up the hold string. If the constant offset is zero, hold it there
	String conv_exp_HS="00",scaleSolvent_HS
	
	//If the user sets the initial guess for the y offset to 0
	//it will be held there
	If(decay_params[2]==0)
		conv_exp_HS+="1"
	Else
		conv_exp_HS+="0"
	EndIf
	
	//If the user sets the initial guess for the constant offset to 0, it will be held there
	If(decay_params[3]==0)
		conv_exp_HS+="1"
	Else
		conv_exp_HS+="0"
	EndIf
	
	//Set up the fit string with the conv_expMulti
	String F_String = "{conv_expMulti, "+nameofwave(decay_params)+",HOLD=\""+conv_exp_HS+"\"}"
	
	//Apparently the WAVE keyword is necessary to create a wave reference
	WAVE solStruct.solvent = tempSol
	
	If(holdoffset)
		scaleSolvent_HS="01"
	Else
		scaleSolvent_HS=""
	EndIF
	//Add the scaleSolvent function to the fit string
	F_String += "{scaleSolvent, scaleFactors,HOLD=\""+scaleSolvent_HS+"\",STRC=solStruct}"
		
	//FIT!
	Variable V_FitError=0,V_FitNumIters=0, V_FitMaxIters=2000
	FuncFit/ODR=1/Q=(q)/L=(numpnts(data))/N=1/W=2 {String = F_String} Data /D /R
	If(!q)
		Print "Number of fit iterations = " + num2str(V_FitNumIters)
	EndIf
	
	//Adjust the scaling and position of the solvent data, if it's on the top graph
	If((StrSearch(TraceNameList("",";",1),NameOfWave(Solvent),2)>=0) && !V_FitError)
		ModifyGraph offset($NameOfWave(Solvent))={scaleFactors[1],0}
		ModifyGraph muloffset($NameOfWave(Solvent))={0,scaleFactors[0]}
	EndIf
			
	Return V_FitError
End

Function FitTAMatrix(data,solvent,decay_Coefs,solScale,soloffset,[startpnt,endpnt,holdoffset])
	//A function for fitting a transient absorption decay to a sum of exponentials convoluted with
	//a gaussian instrument response function (see FitFuncs.ipf for more details) and a scaled,
	//shifted, solvent TA. In this case this is done for every column of the matrix.
	WAVE Data					//The TA of the sample
	WAVE Solvent				//The TA of the solvent
	WAVE decay_Coefs		//Your guesses for the exponential parameters
	Variable solScale,solOffset	//Your guess for the solvent scale, shift parameters
	Variable holdOffset			//1 means that the shift parameter will be held constant
	
	//Optional variables, startpnt is the start point and endpnt is the end point
	//whild holdOffset is a boolean which controls whether the soloffset should be held fixed or not 
	Variable startpnt, endpnt
	//Other variables to be used later are declared here
	Variable i=0, length=DimSize(data,0)//Length is the number of Pixels (1340 for the PIXIS)
	Variable V_FitError=0
	
	String conv_exp_HS,scaleSolvent_HS, F_String,pixelerror=""
	
	Print "These are my initial decay coefficients: "
	Print "W\t\t= " +num2str(decay_Coefs[0])
	Print "t0\t\t= " +num2str(decay_Coefs[1])
	Print "y0\t\t= " +num2str(decay_Coefs[2])
	Print "Offset\t= " +num2str(decay_Coefs[3])
	For(i=4;i<numpnts(decay_Coefs);i+=2)
		Print "A"+num2str(i/2-1)+"\t\t= " +num2str(decay_Coefs[i])
		Print "tau"+num2str(i/2-1)+"\t= " +num2str(decay_Coefs[i+1])
	EndFor
	
	Print "And for copying: "
	Print decay_Coefs
	
	Print ""
	Print "Beginning fitting: press ESC to abort."
	//Set up my temp coefficients, Decay_coefs1 hold the last good set of coefficients
	Duplicate/O decay_Coefs Decay_Coefs1
	//temp_coefs holds the current set
	Duplicate/O decay_Coefs temp_Coefs
	
	If(!WaveExists($("RES_"+nameofwave(data))))
		//Only make a new wave if it doesn't already exist
		//Makes a wave to hold the residuals
		Make/D/O/N=(DimSize(data,0),DimSize(data,1)) $("RES_"+nameofwave(data))
		//Set the scale
		SetScale/P y DimOffset(data,1),DimDelta(data,1),"", $("RES_"+nameofwave(data))
	EndIf
	
	WAVE res_matrix=$("RES_"+nameofwave(data))
	
	//Deal with default behaviour
	If(ParamIsDefault(holdoffset))
		holdoffset = 0
	EndIf
	
	If(ParamIsDefault(startpnt))
		startpnt = 0
	EndIf
	
	If(ParamIsDefault(endpnt))
		endpnt = length
	EndIf
	
	//Make my waves
	//need to rename these programically
	Make/D/O/N=(DimSize(data,0)) t0Wave, widthWave, scaleWave, offsetwave
	
	//TempFit holds a column of the matrix to be fit, and tempSol holds the corresponding column of the solvent
	Make/D/O/N=(DimSize(data,1)) dataVector
	Make/D/O/N=(DimSize(solvent,1)) solVector
	
	SetScale/P x DimOffset(data,1),DimDelta(data,1),"", dataVector
	SetScale/P x DimOffset(solvent,1),DimDelta(solvent,1),"", solVector//Set the scale so the fit parameters
	//mean something useful
	
	Variable timerRefNum = startMSTimer
	
	Printf "Processing "+NameOfWave(data)+":\t -/"
	
	for(i=startpnt;i<endpnt;i+=1)
		
		If(mod(i,20)==0)
			PrintF "*"
		EndIf
		DoUpdate //update the graphs, can really slow things down, but its nice to watch
		
		dataVector = data[i](x)//Assign column i to tempFit
		solVector = solvent[i](x)//Assign column i to tempFit use x scaling, which should be the same
		
		V_fiterror = FitTAwSol(dataVector, solVector, temp_Coefs,solScale,soloffset,holdoffset=holdoffset,q=1)
		
		//Make wave references to the just created waves
		WAVE w_sigma
		WAVE res_tempfit = res_DataVector
		WAVE scaleFactors
		
		If(V_fiterror)
			scaleWave[i] = NaN
			offsetwave[i] = NaN
			widthWave[i] = NaN
			t0Wave[i] = NaN
			temp_coefs = decay_coefs1 //Make the next set the last good set
			pixelerror += num2str(i)+"\r"
		Else
			scaleWave[i] = scaleFactors[0]
			offsetwave[i] = scaleFactors[1]
			widthWave[i] = temp_coefs[0]
			t0Wave[i] = temp_Coefs[1]
			res_matrix[i][] = res_tempfit[q]
			decay_coefs1 = temp_coefs //Update the last good set
		EndIf
		
		if (GetKeyState(0) & 32)	// Is Escape key pressed now?
			Printf "User abort: "
			Break
		EndIf
		
	EndFor
	
	Print "/-"
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	Print pixelerror
	
	Return 0
	
End

 Function PadWave(myWave,n)
	//A function that will pad a wave with n 0s at the front and back
	//For use in fitMatrixTA so that the solvent response may be shifted in time
	
	WAVE myWave //the wave to pad
	Variable n// number of points to add to the front and the back
	
	Variable dx = deltax(myWave)//deltax for the original wave
	Variable lx = leftx(myWave)//leftx for the original wave
	
	Variable newlx = lx-dx*n//new leftx after padding
	
	InsertPoints numpnts(myWave), n, myWave//insert n points at the end
	InsertPoints 0, n, myWave//insert n points at the beginning
	
	SetScale/P x, newlx, dx, myWave//scale the newly padded wave
	
	Return 0
End

Function ISRSDecay(w,t) : FitFunc
	Wave w
	Variable t
	
	//Analytical convoluted exponential function/damped cosine functions
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation = A*(1 + Erf((t-t0)/w) +Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2)) (1 + Erf((t-t0)/w - w/(2*t1))))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=Heaviside amp
	//|  w[4]=A1, the decay (1/e) time of the first exponential
	//|  w[5]=t1 scaling coefficient
	//|  etc...
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 10
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = H
	//CurveFitDialog/ w[i] = a1
	//CurveFitDialog/ w[i+1] = t1
	//CurveFitDialog/ w[i+2] = f (in cm-1)
	//CurveFitDialog/ w[i+3] = d (in rads)
	
	t-=w[1]
	Variable val=w[3]*(1 + Erf(t/w[0])),npts=numpnts(w),i=4
	do
		if( i>=npts )
			break
		endif
		val+= w[i]*Exp((-4*t*w[i+1] + w[0]^2)/(4*w[i+1]^2))*(1 + Erf(t/w[0] - w[0]/(2*w[i+1])))*Cos(2*pi*t*w[i+2]*kspeedoflight+w[i+3])
		i+=4
	while(1)
	val/=2
	val+=w[2]
	return val
End

Function FitISRSDecay(Data, solvent, decay_params,solScale,soloffset,[holdoffset])
	WAVE Data
	WAVE Solvent
	WAVE decay_params
	
	Variable solScale,solOffset,holdOffset
	//Deal with default behaviour
	If(ParamIsDefault(holdoffset))
		holdoffset = 0
	EndIf
	
	Print "Initial Parameters: "
	Print decay_params
	
	Duplicate/O Solvent tempSol
	
	PadWave(tempSol,1000)
	
	//Declare my structure for fitting
	Struct scaledSolventStruct solStruct
	
	//Initialize my fitting parameters for the solvent
	Make/D/O/N=2 scaleFactors = {solScale,soloffset}
	DoUpdate //update the graphs, can really slow things down, but its nice to watch
	
	//Set up the hold string. If the constant offset is zero, hold it there
	String conv_exp_HS="00",scaleSolvent_HS
	
	If(decay_params[2]==0)
		conv_exp_HS+="1"
	Else
		conv_exp_HS+="0"
	EndIf
	
	If(decay_params[3]==0)
		conv_exp_HS+="1"
	Else
		conv_exp_HS+="0"
	EndIf
	
	Variable i=0
	
	//Go through decay params, if the oscillation period is 0 then fix those parameters in the fit
	//i.e. that component should be a pure decay.
	For(i=6;i<numpnts(decay_params);i+=+4)
		If(decay_params[i]==0)
			conv_exp_HS+="0011"
		Else
			conv_exp_HS+="0000"
		EndIf
	EndFor
	
	//Set up the fit string with the conv_expMulti
	String F_String = "{ISRSDecay, "+nameofwave(decay_params)+",HOLD=\""+conv_exp_HS+"\"}"
	
	//Apparently the WAVE keyword is necessary to create a wave reference
	WAVE solStruct.solvent = tempSol
	
	If(holdoffset)
		scaleSolvent_HS="01"
	Else
		scaleSolvent_HS=""
	EndIF
	//Add the scaleSolvent function to the fit string
	If(solScale!=0)
		F_String += "{scaleSolvent, scaleFactors,HOLD=\""+scaleSolvent_HS+"\",STRC=solStruct}"
	EndIF
	
	//FIT!
	Variable V_FitError=0,V_FitNumIters=0
	FuncFit/L=(numpnts(data))/N=1/W=2 {String = F_String} Data /D /R
	Print "Number of fit iterations = " + num2str(V_FitNumIters)
	
	//Adjust the scaling and position of the solvent data, if it's on the top graph
	If((StrSearch(TraceNameList("",";",1),NameOfWave(Solvent),2)>=0) && !V_FitError)
		ModifyGraph offset($NameOfWave(Solvent))={scaleFactors[1],0}
		ModifyGraph muloffset($NameOfWave(Solvent))={0,scaleFactors[0]}
	EndIf
	
	Return V_FitError
End

Function recreateISRS(w,n,stepsize)
	Wave w
	Variable n, stepsize
	
	Variable SpeedOfLight=2.99792458e-5
	
	Make/O/D/N=(n) rec_Osc = 0, rec_decay =0 
	
	SetScale/p x, -500, stepsize, rec_Osc,rec_decay
	
	Variable npts=numpnts(w)
	Variable i=4
	rec_decay=w[3]*(1 + Erf((x-w[1])/w[0]))
	do
		if( i>=npts )
			break
		endif
		If(w[i+2]!=0)
			rec_osc+= w[i]*Exp((-4*(x-w[1])*w[i+1] + w[0]^2)/(4*w[i+1]^2))*(1 + Erf((x-w[1])/w[0] - w[0]/(2*w[i+1])))*Cos(2*pi*(x-w[1])*w[i+2]*speedoflight+w[i+3])
		Else
			rec_decay+= w[i]*Exp((-4*(x-w[1])*w[i+1] + w[0]^2)/(4*w[i+1]^2))*(1 + Erf((x-w[1])/w[0] - w[0]/(2*w[i+1])))
		Endif
		i+=4
	while(1)
	rec_decay/=2
	rec_decay+=w[2]
	rec_osc/=2
	
	FFT/OUT=4/DEST=rec_Osc_FFT rec_Osc
	SetScale/p x, 0, 1/(n*stepsize*speedoflight), rec_Osc_FFT
	
End

//Function MatrixSolventSubtract(data,solvent,[startpt,endpt])
//	//This may be obselete soon
//	Wave data, solvent
//	variable startpt, endpt
//	Variable i=0, length=DimSize(data,0)//Length is the number of Pixels (1340 for the PIXIS)
//	Variable V_FitError=0, t0=0, A=0, C=0,w=0
//	Variable V_min, V_minloc, V_max, V_maxloc
//	
//	If(ParamIsdefault(startpt))
//		startpt = 0
//	Endif
//	
//	If(ParamIsdefault(endpt))
//		endpt = length
//	Endif
//	
//	//Make my waves
//	Make/D/O/N=(DimSize(data,1)) tempData, tempSol, tempData_ns//This wave is as long as there are timepoints
//	Make/D/O/N=(DimSize(data,0)) scale=0//This wave is as long as there are pixels
//	SetScale/P x DimOffset(data,1),DimDelta(data,1),"", tempData, tempSol, tempData_ns//Set the scale
//	
//	Duplicate/O Data $(nameofwave(data)+"_ns")
//	Wave Data_ns=$(nameofwave(data)+"_ns")
//	Data_ns=0
//	
//	Variable timerRefNum = startMSTimer
//	
//	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
//	
//	Variable delta=0, POI=0//POI is point of interest
//	
//	for(i=startpt;i<endpt;i+=1)
//		
//		tempData = Data[i][p]//Assign column i to tempData
//		tempSol = solvent[i][p]//Assign column i to tempSol
//			
//		Wavestats/Q/Z tempSol
//		
//		delta = abs(V_maxloc-V_minloc)
//		
//		If(abs(V_min)>abs(V_max))//We have an inverted shape XPM
//			POI = V_minloc-delta
//		Else //we have a normal shape
//			POI = V_maxloc-delta
//		EndIf
//		
//		scale[i]=abs(tempData(POI)/tempSol(POI))
//		If(scale[i]>1)
//			scale[i]=1
//		EndIf
//		
//		tempData_ns=tempData-scale[i]*tempSol
//		
//		Data_ns[i][]=tempData_ns[q]
//		DoUpdate
//		If(mod(i,20)==0)
//			PrintF "*"
//		EndIf
//		
//	EndFor
//	Print "/-"
//	
//	Variable microSeconds = stopMSTimer(timerRefNum)
//	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
//	
//	Return 0
//End
//
//Function MatrixSolventSubtract2(data,solvent,scales)
//	//This may be obselete soon
//	Wave data, solvent,scales
//	Variable i=0, length=numpnts(scales)//Length is the number of Pixels (1340 for the PIXIS)
//	
//	//Make my waves
//	Make/D/O/N=(DimSize(data,1)) tempData, tempSol, tempData_ns//This wave is as long as there are timepoints
//	SetScale/P x DimOffset(data,1),DimDelta(data,1),"", tempData, tempSol, tempData_ns//Set the scale
//	
//	Duplicate/O Data $(nameofwave(data)+"_ns")
//	Wave Data_ns=$(nameofwave(data)+"_ns")
//	Data_ns=0
//	
//	Variable timerRefNum = startMSTimer
//	
//	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
//	
//	for(i=0;i<length;i+=1)
//		
//		tempData = Data[i][p]//Assign column i to tempData
//		tempSol = solvent[i][p]//Assign column i to tempSol
//		
//		tempData_ns=tempData-scales[i]*tempSol
//		
//		Data_ns[i][]=tempData_ns[q]
//		DoUpdate
//		If(mod(i,20)==0)
//			PrintF "*"
//		EndIf
//		
//	EndFor
//	Print "/-"
//	
//	Variable microSeconds = stopMSTimer(timerRefNum)
//	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
//	
//	Return 0
//End
//
//Function FitTAMatrix(Matrix,decay_Coefs,[startpnt,endpnt])
//	//Same as XCMax except it uses the FindPeak function to find the center and width at each pixel
//	
//	WAVE Matrix, decay_Coefs
//	
//	Variable startpnt, endpnt
//	Variable i=0, length=DimSize(Matrix,0)//Length is the number of Pixels (1340 for the PIXIS)
//	Variable V_FitError=0
//	String holdString
//	
//	Duplicate/O decay_Coefs Decay_Coefs1
//	Duplicate/O decay_Coefs temp_Coefs
//	
//	If(!WaveExists($("RES_"+nameofwave(matrix))))
//		//Only make a new wave if it doesn't already exist
//		Make/D/O/N=(DimSize(Matrix,0),DimSize(Matrix,1)) $("RES_"+nameofwave(matrix))
//		SetScale/P y DimOffset(Matrix,1),DimDelta(Matrix,1),"", $("RES_"+nameofwave(matrix))//Set the scale so the fit parameters
//	EndIf
//	
//	WAVE res_matrix=$("RES_"+nameofwave(matrix))
//	
//	If(ParamIsDefault(endpnt))
//		endpnt = length
//	EndIf
//	
//	//Make my waves
//	Make/D/O/N=(DimSize(Matrix,0)) t0Wave
//	
//	Make/D/O/N=(DimSize(Matrix,1)) tempFit//This wave is as long as there are timepoints
//	
//	SetScale/P x DimOffset(Matrix,1),DimDelta(Matrix,1),"", tempFit//Set the scale so the fit parameters
//	//mean something useful
//	
//	Variable timerRefNum = startMSTimer
//	
//	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
//	
//	for(i=startpnt;i<endpnt;i+=1)
//		
//		If(mod(i,20)==0)
//			PrintF "*"
//		EndIf
//		DoUpdate
//		
//		tempFit = Matrix[i][p]//Assign column i to tempFit
//		
//		V_FitError=0
//		If(temp_coefs[3]==0)
//			holdString="0001"
//		Else
//			holdString=""
//		EndIf
//		FuncFit/H=holdString/Q/N=1/W=2  conv_expMulti temp_Coefs tempFit /D /R /A=0//automatically make fit an residuals
//		
//		WAVE w_sigma
//		WAVE res_tempfit
//		
//		//These declarations insures that the two waves exist, really only important if this macro is run
//		WAVE w_sigma
//		//print w_coef
//		If(V_fiterror)
//			temp_Coefs=NaN
//			w_sigma=NaN
//			t0Wave[i]=temp_Coefs[1]
//			temp_coefs=decay_coefs1 //Make the next set the last good set
//			Print "Shit" + num2str(i)
//		Else
//			res_matrix[i][]=res_tempfit[q]
//			t0Wave[i]=temp_Coefs[1]
//			decay_coefs1=temp_coefs //Update the last good set
//		EndIf
//		
//	EndFor
//	
//	Print "/-"
//	
//	Variable microSeconds = stopMSTimer(timerRefNum)
//	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
//	
//	Return 0
//	
//End