#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Menu "Graph"
	Submenu "Graph Style"
		"Fancify FSRS/O 1", /Q, FancifyFSRS()
		"Fancify TA/O 2", /Q, FancifyTA()
		"-"
		"Color Waves/O 4", /Q, ColorWaves()
		"Rainbow Waves"
	End
End

Function FancifyFSRS()
	//A quick macro to format a plot in a good style
	//box the figure, add minor ticks to the bottom, and scale the left axis to be mOD
	ModifyGraph minor(bottom)=1,prescaleExp(left)=3
	//Label the axes
	Label left "Raman Gain (mOD)"
	Label bottom "Raman Shift (cm\\S-1\\M)"
	SetGraphSizeACS()
End

Function FancifyTA()
	//Similar to FancifyFSRS()
	ModifyGraph minor(bottom)=1,prescaleExp(bottom)=-3
	Label left "Transient Absorption (a.u.)"
	Label bottom "Delay (ps)"
	SetGraphSizeACS()
End

STATIC Function SetGraphSizeACS()
	//A simple function to set the graph size to ACS standards
	//A 3.3 square
	ModifyGraph mirror=2
	//Autoscale only visible data
	SetAxis/A=2 left
	//Set the margins and the graph width such that the figure is a 3.3 in square
	//the correct dimensions for a single column figure for ACS publications
	ModifyGraph margin(left)=36,margin(bottom)=36,margin(top)=12,margin(right)=12
	ModifyGraph width=192,height=192
	//3.3 in is a little small for most screens so we'll expand the figure by 50%
	ModifyGraph expand=1.5
End

Function graphKineticData(freq)
	//The below function can be used to quickly plot extracted amplitudes and frequencies 
	//on one plot with the fits of those extracted quantities. Similar to plotFit().
	Variable freq
	String sFreq = num2str(freq)
	WAVE timepoints=root:timepoints
	//Error checking
	If(!WaveExists($("fit_fitfreq_"+sfreq))||!WaveExists($("fit_fitamp_"+sfreq)))
		beep
		Print "You haven't fit the kinetics yet!"
		Return 0
	ElseIf(!WaveExists($("fitfreq_"+sfreq))||!WaveExists($("fitamp_"+sfreq)))
		beep
		Print "The waves you tried to plot do not exist!"
		Return 0
	EndIf
	
	//Now we can begin to display
	Display/N=KineticData  $("fitamp_"+sfreq) vs TImepoints as ("Mode "+sfreq)
	AppendToGraph/L=Freq $("fitfreq_"+sfreq) vs TImepoints
	AppendToGraph $("fit_fitamp_"+sfreq)
	AppendToGraph/L=L1 $("fit_fitfreq_"+sfreq)
	//Make the graph prettier
	ModifyGraph mode($("fitamp_"+sfreq))=3,mode($("fitfreq_"+sfreq))=3
	ModifyGraph marker($("fitamp_"+sfreq))=5,marker($("fitfreq_"+sfreq))=8
	ModifyGraph lStyle($("fit_fitamp_"+sfreq))=3,lStyle($("fit_fitfreq_"+sfreq))=3
	ModifyGraph rgb($("fit_fitamp_"+sfreq))=(0,0,0),rgb($("fit_fitfreq_"+sfreq))=(0,0,0)
	ModifyGraph lblPosMode(left)=1,lblPosMode(Freq)=1
	ModifyGraph freePos(L1)=0
	ModifyGraph axisEnab(left)={0.5,1}
	ModifyGraph axisEnab(L1)={0,0.5}
	//Set up the labeling
	ModifyGraph prescaleExp(left)=3
	Label left "Raman Gain (mOD)"
	Label Freq "Raman Shift (cm\\S-1\\M)"
	ModifyGraph prescaleExp(bottom)=-3
	Label bottom "Time (ps)"
	//Scale the axes
	SetAxis/A=2 left
	SetAxis/A=2 Freq
	//Add error bars
	ErrorBars $("fitamp_"+sfreq) Y,WAVE=($("fitamperror_"+sfreq),$("fitamperror_"+sfreq))
	ErrorBars $("fitfreq_"+sfreq) Y,WAVE=($("fitfreqerror_"+sfreq),$("fitfreqerror_"+sfreq))
	//Add some annotations
	TextBox/C/N=textfreq/B=1/F=0/A=MC/X=20/Y=-25 "\\F'Symbol'n\\B0\\M\\F]0 = ? ± ? cm\\S-1\\M\rA\\B1\\M = ? ± ? cm\\S-1\\M"
	AppendText/N=textfreq "\\F'Symbol't\\B1\\M\\F]0 = ? ± ? fs\rA\\B2\\M = ? ± ? cm\\S-1\\M\r\\F'Symbol't\\B2\\M\\F]0 = ? ± ? ps"
	TextBox/C/N=textamp/B=1/F=0/A=MB/X=20/Y=75 "\\F'Symbol't\\B1\\M\\F]0 = ? ± ? fs\r\\F'Symbol't\\B2\\M\\F]0 = ? ± ? fs"
	AppendText/N=textamp"\\F'Symbol't\\B3\\M\\F]0 = ? ± ? ps"
	//Draw a line in between plots
	SetDrawLayer UserFront
	DrawLine 0,0.5,1,0.5
EndMacro

Function StackPlot(base,offset,[shift,pnt,invert])
//Base is a string with one or more wildcards so that a family of waves may be plotted e.g. "*spec"
//offset is the amount you want the waves to be offset
	String base		//Match string to find the waves to plot
	Variable offset	//The amount to offset each wave relative to the one before
	WAVE shift		//Optional: xaxis, if the waves don't have internal scaling
	Variable pnt	//Optional:	offset the waves relative to this point
	Variable invert	//Optional: inver the ordering of the waves
	
	//Internal variables
	Variable i,totoffset
	String wn = ""
	String wl = WaveList(base,";","")//gathering family of waves
	If(!ParamIsDefault(invert))//If anything but default
		//Invert the list
		wl=InvertList(wl)
	EndIf
	
	//Set up display window
	Display/N=StackPlot as "Stack plot "+base
	
	i=0
	
	Do
		wn = StringFromList(i,wl)
		If (strlen(wn) ==0)
			//No more waves
			Break
		EndIf
		//Increase iterator
		i+=1
		//Append the wave
		If(!ParamIsDefault(shift))
			AppendToGraph $wn vs shift
		Else
			AppendToGraph $wn
		EndIf
	While(1)
	
	//Do the offsetting
	If(ParamIsDefault(pnt))
		WaveOffset(offset)
	Else
		WaveOffset(offset,pnt=pnt)
	EndIf
	
	//Color the waves
	ColorWaves()
End

Function WaveOffset(offset,[base,pnt,rev])
//Sequentially offsets the waves matching base by an amount offset in the top graph
//***NOTE: this will offset the waves in the order they are in on the graph***
	Variable offset	//The amount to offset the traces
	String Base		//Used if you only want a subset of the traces to be offset
	Variable pnt	//If used it sets the pixel pnt in each wave to 0 first before offsetting
	Variable rev	//Reverse order?
	
	String tl=TraceNameList("", ";",1)	//Pull all of the traces off the top graph
	
	If(!ParamIsDefault(base))
		tl = ListMatch(tl,base)			//remove the desired subset
	EndIf
	
	If(!ParamIsDefault(rev))
		tl=ReverseList(tl)
	EndIf
	
	String tn
	Variable i=0
	//Loop through the traces offseting each one
	do
		tn= StringFromList(i,tl)
		If( strlen(tn) == 0 )
			break
		EndIf
		Wave myTrace = $tn
		If(ParamIsDefault(pnt))
			ModifyGraph offset($tn)= {0,offset*i}	// Go through the traces on the top graph and offset them
		Else
			ModifyGraph offset($tn)= {0,offset*i-myTrace[pnt]}
		EndIf
		i += 1
	while (1)	// exit is via break statement
End

Function colorTimeWaves(timepoints) : Graphstyle 
//Colors the waves in the top graph in a pleasing manner and adds a legend
	WAVE timepoints
	
	Variable timerange = abs(timepoints[numpnts(timepoints)]-timepoints[0])
	
	Make/D/O/N=6 Red={16385,0,26411,39321,26214,0,65535}
	Make/D/O/N=6 Green={28398,0,1,1,0,0}
	Make/D/O/N=6 Blue={65535,65535,52428,31457,10485,0}
	
	Variable k,km 
	
	String tnl = TraceNameList("", ";",1)
	
	String name = "" //to hold the name of the WAVE
	String pORm, num
	
	k = ItemsInList(tnl)
	If (k < 2)
		return -1
	EndIf
	km=k//number of items
	Interpolate2/T=1/N=(km)/I=3/Y=RedInterp Red
	Interpolate2/T=1/N=(km)/I=3/Y=GreenInterp Green
	Interpolate2/T=1/N=(km)/I=3/Y=BlueInterp Blue
	PauseUpdate
	for(k=0;k<km;k+=1)
		name = StringFromList(k,tnl)
		SplitString/E=("(p|m)([[:digit:]]+)") name, pORm, num
		ModifyGraph rgb($name)=(Round(RedInterp[str2num(num)/timerange*km]),Round(GreenInterp[str2num(num)/timerange*km]),Round(BlueInterp[str2num(num)/timerange*km]))
	EndFor
	KillWaves Red, Green, Blue, RedInterp, GreenInterp, BlueInterp
End

Function labelTimes([pnt])
//Label time resolved spectra, for instance, those displayed using stackplot
	Variable pnt	//Where you want the labels to be
	//Set up the list of relevant waves
	String 	wl=RemoveFromList(WaveList("fit*",";",""),TraceNameList("", ";",1))
	
	Variable length = ItemsInList(wl);
	If(length == 0)
		Return -1//Function failed, no waves in top window!
	EndIf
	
	Variable i = 0
	Variable myOffset = 0 //variable for storing the offset
	
	String name = ""
	String timepoint, pORm, num
	
	//Before doing anything draw the rectangle
	WAVE XWave = $XWaveName("",StringFromList(0,wl))
	
	SetDrawEnv fsize=12
	
	If(ParamIsDefault(pnt)||!WaveExists(XWave))
		SetDrawEnv xcoord= prel,ycoord= left, textxjust= 1, save
	Else
		SetDrawEnv xcoord= bottom,ycoord= prel, textxjust= 1, linefgc= (65535,65535,65535), save
		DrawRect XWave[pnt-60],0,XWave[pnt+60],1
		SetDrawEnv ycoord= left,save
	EndIf
	
	For(i=0;i<length;i+=1)
		name = StringFromList(i,wl)
		SplitString/E=("(p|m)([[:digit:]]+)") name, pORm, num//Extract number info
		
		If(CmpStr(num,"")!=0)
			If(abs(str2num(num)) >= 1000)
				timepoint = num2str(str2num(num)/1000)+ " ps"
			Else
				timepoint = num + " fs"
			EndIf
			
			If(CmpStr(pORm,"m")==0 && str2num(num)!=0)
				timepoint = "-"+timepoint
			EndIf
			//Retrieve offset
			myOffset = GetNumFromModifyStr(traceinfo("",name,0),"offset","{",1)
			
			SetDrawEnv textrgb= (GetNumFromModifyStr(traceinfo("",name,0),"rgb","(",0),GetNumFromModifyStr(traceinfo("",name,0),"rgb","(",1),GetNumFromModifyStr(traceinfo("",name,0),"rgb","(",2))
			
			If(ParamIsDefault(pnt)||!WaveExists(XWave))
				DrawText 0.9,myOffset,timepoint
			Else
				WAVE XWave = $XWaveName("",name)
				WAVE myWave =$name
				DrawText  XWave[pnt],myOffset+myWave[pnt],timepoint
			EndIf
		EndIf
	EndFor
			
	Return 0 //Function ran successfully
End

Function InverseLegend([invert])
//Adds a legend to the graph in the right order and with fancy labels
	Variable invert
	String wl = ""
	If(ParamIsDefault(invert))
		invert = 1
		wl = InvertList(TraceNameList("", ";",1))
	Else
		wl = TraceNameList("", ";",1)
	EndIf
	
	//Collect a list of the waves from the top window and invert it
	Variable length = ItemsInList(wl);
	If(length == 0)
		Return -1//Function failed, no waves in top window!
	EndIf
	Variable i = 0
	
	String name = ""
	String timepoint, pORm, num
	
	Legend/C/N=iLegend/J "" //Make legend
	For(i=0;i<length;i+=1)
		name = StringFromList(i,wl)
		SplitString/E=("(p|m)([[:digit:]]+)") name, pORm, num//Extract number info
		
		If(abs(str2num(num)) >= 1000)
			timepoint = num2str(str2num(num)/1000)+ " ps"
		Else
			timepoint = num + " fs"
		EndIf
		
		If(CmpStr(pORm,"m")==0 && str2num(num)!=0)
			timepoint = "-"+timepoint
		EndIf
		AppendText/N=iLegend "\\s("+ name + ") " + timepoint
		//add to legend
	EndFor
	
	Return 0 //Function ran successfully
End

Function addTags(spectrum, minamp)
	WAVE spectrum
	Variable minamp
	
	WAVE myPeaks = $FindPeakLocations(spectrum, minamp=minamp, name="addTagsPeaks")
	
	Variable length = numpnts(myPeaks)
	Variable i = 0
	for(i=0; i<length; i+=1)
		addtag(spectrum, mypeaks[i])
	EndFor
End

Function removeTags()
	String myTags = AnnotationList("")
	variable i = 0
	String tagname = ""
	do
		tagname = stringfromlist(i,mytags)
		If( strlen(tagname)==0)
			break
		EndIf
		If(stringmatch(tagname, "tag*"))
			Tag/K/N=$tagname
		EndIf
		i+=1
	while(1)
End

Function markPeaks(peaklocations)
	WAVE peaklocations // the WAVE that holds the peak locations
	
	Variable length = numpnts(peaklocations)
	
	//Error checking
	If(length==0)
		Beep
		Print "There are no peaks!"
		Return -1
	EndIf
	
	variable i=0;
	
	//Set my draw environment to add the tags
	SetDrawEnv xcoord= bottom,textxjust= 1,textyjust= 2,textrot= 90, linefgc= (56797,56797,56797), save
	SetDrawLayer UserBack
	
	For(i=0;i<length;i+=1)
		DrawText peaklocations[i],0,num2str(round(peaklocations[i]))+" cm\\S-1\\M"
		DrawLine peaklocations[i],0.08,peaklocations[i],1
	EndFor
	
	Return 0
End

Function addtag(trace, pnt)
// Function adds a tag to the specified trace at the specified point with the wavenumber at that
// pnt as specified by shiftx
	WAVE trace //the trace
	Variable pnt // the pnt
	WAVE shiftx = root:shiftx // shiftx
	
	//error checking
	If(!WaveExists(shiftx))
		DoAlert/T="AddTag Failed!" 0, "There is no shiftx!"
		Return -1
	EndIf
	
	String myLabel
	
	sprintf myLabel, "%d cm\\S-1\\M", shiftx[pnt]
	
	print nameofwave(trace)
	
	Tag/B=1/A=MB/Y=10.00/X=0.00/C/N=$("tag"+num2istr(pnt))/O=90/F=0/L=1/TL=0 $nameofwave(trace), pnt, myLabel
	
	Return 0 // Function excuted normally
End

Function ColorWaves([invert]) 
//Colors the waves in the top graph in a pleasing manner
	Variable invert //do we want to change the way the waves are colored?
	If(ParamIsDefault(invert))
		invert=0
	EndIf
	
	String base
	Make/D/O/FREE/N=6 Red={16385,0,26411,39321,26214,0,65535}
	Make/D/O/FREE/N=6 Green={28398,0,1,1,0,0}
	Make/D/O/FREE/N=6 Blue={65535,65535,52428,31457,10485,0}
	
	Variable k,km 
	String tnl = TraceNameList("", ";", 1)
	
	k = ItemsInList(tnl)
	If (k < 2)
		return -1
	EndIf
	km=k
	Interpolate2/T=1/N=(km)/I=3/Y=RedInterp Red
	Interpolate2/T=1/N=(km)/I=3/Y=GreenInterp Green
	Interpolate2/T=1/N=(km)/I=3/Y=BlueInterp Blue
	
	If(invert)
		Reverse RedInterp, GreenInterp, BlueInterp
	EndIf
	
	PauseUpdate
	for(k=0;k<km;k+=1)
		ModifyGraph rgb[k]=(Round(RedInterp[k]),Round(GreenInterp[k]),Round(BlueInterp[k]))
	EndFor
	KillWaves RedInterp, GreenInterp, BlueInterp
End

Function RainbowWaves()
//Colors the waves in the top graph in a rainbow pattern
	String base
	Make/D/O/FREE/N=13 Red={65535,65535,65535,56172,31207,6241,0,0,0,0,0,6553,32767}
	Make/D/O/FREE/N=13 Green={0,24965,53052,65535,65535,65535,65535,65535,65535,43690,15603,0,0}
	Make/D/O/FREE/N=13 Blue={0,0,0,0,0,0,17873,41704,65535,65535,65535,65535,65535}
	Variable k,km 
	String tnl = TraceNameList("", ";", 1)
	
	k = ItemsInList(tnl)
	If (k < 2)
		return -1
	EndIf
	km=k
	Interpolate2/T=1/N=(km)/I=3/Y=RedInterp Red
	Interpolate2/T=1/N=(km)/I=3/Y=GreenInterp Green
	Interpolate2/T=1/N=(km)/I=3/Y=BlueInterp Blue
	PauseUpdate
	for(k=0;k<km;k+=1)
		ModifyGraph rgb[k]=(Round(RedInterp[k]),Round(GreenInterp[k]),Round(BlueInterp[k]))
	EndFor
	KillWaves RedInterp, GreenInterp, BlueInterp
End

Function displayTimePnt(timepoints)
	//A function to display a the ground subtracted and ground added spectra For a specified delay
	WAVE timepoints //WAVE containing all timepoints
	
	WAVE shiftx = root:shiftx //Set up my wave reference to the raman shift
		
	If(!WaveExists(shiftx)) //Check to see that it exists
		Print "Error, the shiftx WAVE does not exist"
		Return -2
	EndIf
	
	String timelist=""
	Variable i=0
	//Build a list to show in the pop up menu
	For(i=0;i<numpnts(timepoints);i+=1)
		timelist+=num2istr(timepoints[i])+";"
	EndFor
	
	//set up the pop up menu
	String timeStr //Variable to place the user's choice
	//setting up the prompt
	Prompt timeStr, "Time Points", popup timelist
	//displaying the prompt
	DoPrompt "Choose a time point to plot", timeStr
	
	If(V_Flag)
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	//Switch the string to a number
	Variable timePnt = str2num(timeStr)
	
	//switch number back to string but with a prepended "p" or "m"
	//in accordance with the output formatting of the FSRS instrument
	timeStr = myTime(timePnt)
	
	//Setting up my WAVE references so that I can easily access the waves
	WAVE withg = $(timeStr+"_withg")
	WAVE subg = $(timeStr+"_subg")
	
	//Forming a title string For the graph I'm going to make
	String titleString=""
	If(abs(timePnt) >= 1000)
		titlestring += num2str(timePnt/1000)+ " ps"
	Else
		titlestring+=  num2str(timePnt) + " fs"
	EndIf
	
	//Error Checking
	If(!WaveExists(subg))
		Print "Error, the subG WAVE does not exist"
		Return -1
	Else
		//No errors so I can make my graph
		Display/N=$timeStr subG vs shiftx
		titlestring+=" subg"
	EndIf
	
	If(!WaveExists(withg))
		Print "There is no withG WAVE"
	Else
		//If the withg WAVE exists I'll append it here
		AppendToGraph/W=$timeStr withG vs shiftx
		ModifyGraph/W=$timeStr rgb($nameofwave(withG))=(0,0,65535)
		titlestring+="/withg"
	EndIf
	
	//Make the graph look good, and set up the tools to draw a WAVE monotonic, For baseline correction
	FancifyFSRS()
	ShowTools/A
	GraphWaveDraw/M
	//Change the title of the window
	DoWindow/T $timeStr, titleString
	
	//Some output For the history
	Print "Displayed the "+num2str(timePnt)+" fs time point"
	
	Return 0
End