#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function fDispLorFit(w,x)
	WAVE w; Variable x
	
	Variable r= w[0]
	variable npts= numpnts(w),i=1
	variable newx = 0
	do
		if( i>=npts )
			break
		EndIf
		newx = (x-w[i+2])/(w[i+3]/2)
		r += (w[i]+w[i+1]*newx)/(newx^2+1)
		i+=4
	while(1)
	return r
End

Function Fano(w,x)
	WAVE w; Variable x
	
	Variable r= w[0]
	variable npts= numpnts(w),i=1
	variable newx = 0
	do
		if( i>=npts )
			break
		EndIf
		newx = (x-w[i+2])/(w[i+3]/2)
		r += w[i]*(w[i+1]+newx)^2/(newx^2+1)
		i+=4
	while(1)
	return r
End

Function fitTimeSeries2(timepoints, pnt1, pnt2,wavenumber,Coefs,[gnd,Wiggle,Width,Plot,PolyNum,subrangeStart,subrangeEnd,suffix])
//modified David Hoffman
	WAVE timepoints
	Variable pnt1, pnt2
	WAVE/T wavenumber
	WAVE Coefs
	WAVE gnd
	WAVE wiggle
	WAVE Width
	String Suffix
	Variable Plot
	Variable PolyNum
	Variable subrangeStart,subrangeEnd
	
	Variable includeGND = !ParamIsDefault(gnd)		//The user has indicated that the ground state should be included in the fitting procedure
	
	if(ParamIsDefault(Plot))
		Plot=0
	EndIf
	
	if(ParamIsDefault(suffix))
		If(includeGND)
			Suffix="_subg"
		Else
			Suffix="_withg"
		EndIf
	EndIf
	
	If(ParamIsDefault(PolyNum))
		PolyNum = 4
	EndIf
	
	If(ParamIsDefault(Wiggle))
		Make/N=(numpnts(wavenumber))/O/D/FREE Vlimit = 30
	Else
		Duplicate/O/FREE Wiggle Vlimit
	EndIf
	
	If(ParamIsDefault(Width))
		Make/N=(numpnts(wavenumber))/O/D/FREE MaxWidth = 100
	Else
		Duplicate/O/FREE Width MaxWidth
	EndIf
	
	IF(WaveExists(root:shiftx))
		WAVE shiftx=root:shiftx
	Else
		DoAlert/T="fitTimeSeries Failed!" 0, "Why are you trying to run this macro so early!?"
		Return -1
	EndIf
	
	//***************************************************************//
	// time to do some error checking//
	if(numpnts(coefs)!=(numpnts(wavenumber)*4+1)) //Check to see if there are enough
		//coefficients for the number of peaks you will fit to.
		Print "The number of peak coefficients DID NOT match the number of peaks"
		Return 0 // Exit the function
	EndIf
	//***************************************************************//
	
	//This string holds the name of the current time point being fit
	String currenttime
	
	//These variables hold the lengths of the various waves so the don't need to be calculated again.
	Variable lengthT=numpnts(timepoints),lengthW=numpnts(wavenumber)
	//Some index variables for use later
	Variable i,j,k
	
	//Let's see if the user has opted to use a subrange, if not lets set the subrange to the full range
	If(ParamIsDefault(subrangeStart))
		subrangeStart = 0
	EndIf
	
	If(ParamIsDefault(subrangeEnd))
		subrangeEnd = lengthT
	EndIf
	
	//Now that we know we don't have any fundamental errors, let's print out the coefficients,
	// the wavenumbers and the points so that the user can find them later if need be
	
	//Print a row with the wavenumbers
	Print " "
	Print "Initial guesses:"
	PrintF "Wavenumber:"
	For(i=0;i<lengthW;i+=1)
		PrintF "\t%s", wavenumber[i]
	EndFor
	PrintF "\r"
	//Print a row with the amp guesses
	PrintF "Amp:\t\t"
	For(i=0;i<lengthW;i+=1)
		PrintF "\t%g", Coefs[4*i+1]
	EndFor
	PrintF "\r"
	//Print a row with the disp guesses
	PrintF "Disp:\t\t"
	For(i=0;i<lengthW;i+=1)
		PrintF "\t%g", Coefs[4*i+2]
	EndFor
	PrintF "\r"
	//Print a row with the center guesses
	PrintF "Center:\t"
	For(i=0;i<lengthW;i+=1)
		PrintF "\t%g", coefs[4*i+3]
	EndFor
	PrintF "\r"
	//Print a row with the width guesses
	PrintF "Width:\t"
	For(i=0;i<lengthW;i+=1)
		PrintF "\t\t%g", coefs[4*i+4]
	EndFor
	PrintF "\r"
	Print " "
	Print "Now in copiable form"
	Print "Make/T/O/N="+num2str(lengthW)+" "+nameofwave(wavenumber)+"; DelayUpdate"
	Print Wavenumber
	Print "Make/D/O/N="+num2str(numpnts(coefs))+" "+nameofwave(coefs)+"; DelayUpdate"
	Print coefs
	Print " "
	Print "Point A is "+num2str(pnt1)+" ("+num2str(shiftx[pnt1])+" cm-1) and Point B is "+num2str(pnt2)+" ("+num2str(shiftx[pnt2])+" cm-1)."
	Print " "
	Print "We'll be fitting the timepoints in between " + num2str(timepoints[subrangeStart]) + " and " + num2str(timepoints[subrangeEnd-1])
	Print " "
	
	//Making parameter waves to keep track of everything
	Make/T/O/N=(lengthW) fitAreal,FitAdisp,fitfreq,fitwidth
	Make/T/O/N=(lengthW) fitArealerror,FitAdisperror,fitfreqerror,fitwidtherror
	
	For(i=0;I<lengthW;i+=1)
		fitAreal[i]="fitAreal_"+wavenumber[i]
		FitAdisp[i]="FitAdisp_"+wavenumber[i]
		fitfreq[i]="fitfreq_"+wavenumber[i]
		fitwidth[i]="fitwidth_"+wavenumber[i]
		
		fitArealerror[i]="fitArealerror_"+wavenumber[i]
		FitAdisperror[i]="FitAdisperror_"+wavenumber[i]
		fitfreqerror[i]="fitfreqerror_"+wavenumber[i]
		fitwidtherror[i]="fitwidtherror_"+wavenumber[i]
		
		
		Make/D/O/N=(lengthT) $fitAreal[i],$fitAdisp[i],$fitfreq[i],$fitwidth[i]
		Make/D/O/N=(lengthT) $fitArealerror[i],$fitAdisperror[i],$fitfreqerror[i],$fitwidtherror[i]
	EndFor
	
	
	//WAVE coef,coeftemp
	string temp
	variable fitsuccess=0,fitfail=0,lengthC=numpnts(coefs)
	
	Variable V_fitOptions = 4	// suppress progress window
	Variable V_FitError = 0
	Variable V_FitMaxIters=2000
	//make the hold string to hold the baseline to zero we will be including one later
	String H_string="1"
	for(i=1;i<lengthC;i+=4)
		If(coefs[i])
			//The mode is dispersive, let all variables float
			H_string+="0"
		Else
			//The mode is NOT dispersive, hold q to 0
			H_string+="1"
		EndIf
		If(coefs[i+1])
			H_string+="0"
		Else
			H_string+="1"
		EndIf
		H_string+="00"
	EndFor
	
	//Make the Epsilon WAVE
	Make/D/O/N=(lengthC) epsilonWave = 1e-6
	for(i=1;i<lengthC;i+=4)
		epsilonWave[i] = 1e-8//Amplitude Epsilon
	EndFor
	
	//make the constrating WAVE
	Variable numCs = 4
	Make/D/O/T/N=((lengthC-1)*numCs/4) CTextWave //Because for the frequency and width we need
	//TWO constraints
	
	Variable Alimit=1e-8 //constrain the peaks to be greater than zero
	//Variable Vlimit=wiggle // constrain the possible frequencies (important for closely spaced peaks)
	
	k=0
	for(i=0;i<((lengthC-1)*numCs/4);i+=numCs)
		j=i*4/numCs //this counter will take care of indexing the coefficient WAVE.
		
		//frequency constraint, adjustable by the user, look above
		CTextWave[i]="K"+num2str(j+3)+">"+num2str(coefs[j+3]-Vlimit[k])
		CTextWave[i+1]="K"+num2str(j+3)+"<"+num2str(coefs[j+3]+Vlimit[k])
		
		//width constraint, these values seem to work for the red table
		CTextWave[i+2]="K"+num2str(j+4)+">3"
		CTextWave[i+3]="K"+num2str(j+4)+"<"+num2str(MaxWidth[k])
		k+=1
	EndFor
	
	For(i=1;i<numpnts(coefs);i+=4)
		//Amplitude constraint, it must be positive
		If(coefs[i])
			k=numpnts(CTextWave)
			InsertPoints k, 1, CTextWave
			CTextWave[k]="K"+num2str(i)+">"+num2str(Alimit)
		EndIf
		If(coefs[i+1])
			k=numpnts(CTextWave)
			InsertPoints k, 1, CTextWave
			CTextWave[k]="K"+num2str(i+1)+">"+num2str(Alimit)
		EndIf
	EndFor
	Print "This is my constraint WAVE:"
	Print CTextWave
	Print " "
	Print "H_string is "+H_string
	
	Duplicate/O coefs tempPeak_Coefs
	
	String F_String=""
	F_String="{fDispLorFit, tempPeak_Coefs, hold=\""+H_string+"\",EPSW=epsilonWave}"
	
	If(PolyNum>2)
		Print "Including an order", polynum-1, "polynomial for the baseline"
		Make/D/O/N=(PolyNum) tempBaseln_Coefs=0
		F_String+="{poly_XOffset "+num2str(polynum)+", tempBaseln_Coefs}"
	ElseIf(PolyNum!=0)
		Print "Including line for the baseline"
		Make/D/O/N=2 tempBaseln_Coefs=0
		F_String+="{line, tempBaseln_Coefs}"
		Polynum=2
	Else
		Print "No baseline!"
	EndIf
	
	If(includeGND)
		//Declare my structure for fitting
		Struct scaledGroundStruct gndStruct
		//Apparently the WAVE keyword is necessary to create a wave reference
		WAVE gndStruct.gnd = gnd
		WAVE gndStruct.shift = shiftx
		
		Make/D/O/N=1 scaleFactor = -0.5
		
		Make/D/O/N=(LengthT) ScaleFactors, sigma_ScaleFactors
		
		F_String += "{scaleGround, scaleFactor,STRC=gndStruct}"
	Else
		Print "No ground included."
	EndIf
	
	Print ""
	
	Printf " Start fitting: -/"
	
	//String to hold the timepoints which weren't fit properly
	Variable v_avg
	String badFits = "Bad Fits:\r"
	//PauseUpdate
	for(i=subrangeStart;i<subrangeEnd;i+=1)
		currenttime = myTime(timepoints[i])+suffix
		
		//tempBaseln_Coefs = 0
		Wavestats/Q/R=[pnt1,pnt2] $currenttime
		tempBaseln_Coefs[0]=v_avg
		
		FuncFit/L=1340/N/Q/M=0/W=2 {string = F_String} $currenttime[pnt1, pnt2] /X=root:shiftx /D /C=CTextWave
		
		WAVE W_sigma
		
		If(includeGnd)
			ScaleFactors[i] = ScaleFactor[0]
			sigma_scaleFactors[i] = w_sigma[numpnts(w_sigma)-1]
		EndIf
		
		for(k=0;k<lengthW;k+=1)
			//Keeping track of all our parameters and putting them in reasonably
			//named waves!
			WAVE  wfitAreal=$fitAreal[k]
			WAVE  wfitAdisp=$fitAdisp[k]
			WAVE  wfitfreq=$fitfreq[k]
			WAVE  wfitwidth=$fitwidth[k]
			
			WAVE  wfitArealerror=$fitArealerror[k]
			WAVE  wfitAdisperror=$fitAdisperror[k]
			WAVE  wfitfreqerror=$fitfreqerror[k]
			WAVE  wfitwidtherror=$fitwidtherror[k]
			
			if(V_FitError!=0)
				wfitAreal[i]=NaN
				wfitAdisp[i]=NaN
				wfitfreq[i] =NaN
				wfitwidth[i]=NaN
				
				wfitArealerror[i]=NaN
				wfitAdisperror[i]=NaN
				wfitfreqerror[i] =NaN
				wfitwidtherror[i]=NaN
			else
				wfitAreal[i]=tempPeak_Coefs[4*k+1]
				wfitArealerror[i]=w_sigma[4*k+1]
				
				wfitAdisp[i]=tempPeak_Coefs[4*k+2]
				wfitAdisperror[i]=w_sigma[4*k+2]
				
				wfitfreq[i] =tempPeak_Coefs[4*k+3]
				wfitfreqerror[i] =w_sigma[4*k+3]
				
				wfitwidth[i]=abs(tempPeak_Coefs[4*k+4])
				wfitwidtherror[i]=w_sigma[4*k+4]
			EndIf
		EndFor
		If(V_FitError!=0)
			//Trying to fit a WAVE that doesn't exist shouldn't count against you.
			If(WaveExists($currenttime))
				badFits += currenttime+"\r"
				fitfail+=1
				Printf "X"// X for fail!
			Else
				Printf " DNE "//DNE=does not exist!
			EndIf
			//print coef
			duplicate/o coefs tempPeak_Coefs
			tempBaseln_Coefs = 0
		Else
			//Print "Fitting successful!"
			fitsuccess+=1
			Printf "O"//O for OK!
			If(PolyNum!=0)
				//Now I will make the fits and the baseline so that these can be plotted by the user.
				Make/D/O/N=1340 $("fit_"+currenttime+"_nb") = fdispLorFit(tempPeak_Coefs,shiftx)
				
				Make/D/O/N=1340 $("fit_"+currenttime+"_bl")
				WAVE tempBaseLineWave = $("fit_"+currenttime+"_bl")
				Duplicate/O $("fit_"+currenttime) $("fit_"+currenttime+"_bl2")
				WAVE tempBaseLineWave2 = $("fit_"+currenttime+"_bl2")
				tempBaseLineWave2=0
				
				//Clear the WAVE in the region in which we're working
				//tempBaseLineWave=tempBaseLineWave*(x<pnt1 && x>pnt2)
				tempBaseLineWave=tempBaseLineWave*(x<pnt1 || x>pnt2)
				//Now fill in only that region
				If(PolyNum>2)
					For(j=0;j<(PolyNum);j+=1)
						tempBaseLineWave+=(tempBaseln_Coefs[j]*(shiftx-shiftx[pnt1])^j)*(x>=pnt1 && x<=pnt2)
						tempBaseLineWave2+=(tempBaseln_Coefs[j]*(x-shiftx[pnt1])^j)
					EndFor
				Else
					tempBaseLineWave2=tempBaseln_Coefs[0]+tempBaseln_Coefs[1]*x
				EndIf
				//We're doing this so that different regions can be fit separately.
			
				WAVE myWave = $currenttime
				Make/D/O/N=1340 $(currenttime+"_nb")=myWave-tempBaseLineWave
			EndIf
		EndIf
		
		DoUpdate
		
		V_FitError=0
		if (GetKeyState(0) & 32)	// Is Escape key pressed now?
			Printf "User abort: "
			Plot=0
			Break
		EndIf
	EndFor
	
	DoUpdate
	
	PrintF "/-\r"
	Print " "
	
	Print "There were "+num2istr(fitsuccess)+" successful fittings and "+num2istr(fitfail)+" fitting failures in this run"
	If(fitfail!=0)
		Print badFits
	EndIf
	
	If(Plot!=0)//Plotting the results!
		for(k=0;k<lengthW;k+=1)
			tabfit2(freq=str2num(wavenumber[k]))
		EndFor
		
		Display/N=PlotFit
		for(k=0;k<lengthW;k+=1)
			PlotFitSub2("areal",str2num(wavenumber[k]),0)
			PlotFitSub2("adisp",str2num(wavenumber[k]),0)
		EndFor
		for(k=0;k<lengthW;k+=1)
			PlotFitSub2("freq",str2num(wavenumber[k]),1)
		EndFor
	EndIf
End

Function tabFit2([freq])
//Displays the results of peak fitting
//The user can chose whether they want to plot the amplitude, the frequency, the width or the area
	Variable freq
	WAVE timepoints = root:timepoints
	If(!WaveExists(timepoints))
		Print "Timepoints WAVE does not exist"
		return -1
	EndIf
	
	if(ParamIsDefault(freq))
	//If the user doesn't specify a frequency ask them for one
		freq=0
		Prompt freq,"Mode frequency:  "
		DoPrompt "Input the frequency of the mode you'd like to be tabulated",freq
		if( V_Flag )
			Print "Cancelled"
			return 0	// user canceled
		EndIf
	EndIf
	String sFreq=num2str(freq)
	If(!WaveExists($("fitAreal_"+sFreq)))
		Print "Error the waves don't exist"
		return -1
	EndIf
	Edit/N=TabFit timepoints as ("Mode " + sFreq + " Tabulated Data")
	AppendToTable $("fitAreal_"+sFreq),$("fitArealerror_"+sFreq)
	AppendToTable $("fitAdisp_"+sFreq),$("fitAdisperror_"+sFreq)
	AppendToTable $("fitfreq_"+sFreq),$("fitfreqerror_"+sFreq)
	AppendToTable $("fitwidth_"+sFreq),$("fitwidtherror_"+sFreq)
	ModifyTable style(Timepoints)=1,alignment(Timepoints)=1,format(Timepoints)=2
End

Function PlotFit2()
//Displays the results of peak fitting
//The user can chose whether they want to plot the amplitude, the frequency, the width or the area
	Variable freq=0
	Variable choice=1
	Variable plotAppend=1
	Prompt freq,"Mode frequency:  "
	Prompt choice,"I want to plot the: ",popup,"A Real;A Disp;Frequency;Width"
	Prompt plotAppend,"Do you want to append the plot to the top graph? ",popup,"No;Yes"
	DoPrompt "Input the frequency of the mode you'd like to be plotted",freq,choice,plotAppend
	if( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	String sFreq
	switch(choice)
		case 1:
			sFreq="Areal"
			break
		case 2:
			sFreq="Adisp"
			break
		case 3:
			sFreq="freq"
			break
		case 4:
			sFreq="width"
			break
		default:
			Print "Error in Switch statement"
			Return 0
	endswitch
	PlotFitSub2(sFreq,freq,plotAppend)
	Print "PlotFitSub2(\""+sFreq+"\","+num2str(freq)+","+num2str(plotAppend)+")"
End

//Helper function for the above plotFit() but can also be used standalone through the command line
Function PlotFitSub2(type,freq,plotAppend)
	String type
	Variable freq, plotAppend
	
	WAVE timepoints = root:timepoints
	If(!WaveExists(timepoints))
		Print "Timepoints WAVE does not exist"
		return -1
	EndIf
	
	String toPlot="fit"+type+"_"+num2str(freq)
	String error="fit"+type+"error_"+num2str(freq)
	if(!WaveExists($toPlot)||!WaveExists($error))
		Print toPlot+" doesn't exist"
		Return 0
	EndIf
	If(plotAppend==1)
		display/N=PlotFit $toPlot vs timepoints as ("Mode "+num2str(freq)+" "+type)
	Else
		AppendToGraph $toPlot vs timepoints
	EndIf
	//ErrorBars $toPlot Y,WAVE=($error,$error)
	//ModifyGraph mode($toPlot)=3,marker($toPlot)=8
	ModifyGraph mirror=2
	ModifyGraph prescaleExp(bottom)=-3
	Label bottom "Delay (ps)"
	SetGraphSizeACS()
	ModifyGraph minor=1
	strswitch(type)
		case "Areal":
			Label left "A\BReal\M (a.u.)"
			break
		case "Adisp":
			Label left "A\BDisp\M (a.u.)"
			break
		case "Width":
			Label left "Width (cm\S-1\M)"
			break
		case "freq":
			Label left "Frequency (cm\S-1\M)"
			ModifyGraph tkLblRot(left)=90, nticks(left)=3
			break
	endswitch
End
