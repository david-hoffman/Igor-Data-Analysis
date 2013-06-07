#pragma rtGlobals=3		// Use modern global access method.stack
#include <Readback ModifyStr> //used to parse strings returned by TraceInfo command
#include ":Graphing"

STATIC Constant kSpeedOfLight=2.99792458e-5 //(cm/fs)

//Setting up my menus and hot keys
Menu "Macros"
	"Plot Fit", PlotFit()
	"Remove point/6", RemovePoints()
	"Calculate Bandwidth\ Duration /4",TransformLimit()
	"-"
	"Auto Baseline/3", /Q, autoBaseline()
End

Function init()
	//A quick initialization macro
	//Makes pixel and shift waves, and fills the shift waves with the peaks For cyclohexane
	Make/D/O/N=10 shift=NaN
	
	//A list of the solvents available
	String solvents = "Cyclohexane;DMSO;Dioxane;Ethanol;Benzene"
	
	Variable solvent = 0 //The user will set this after interacting with the prompt
	
	//Prompt the user to choose a solvent
	Prompt solvent, "What solvent would you like to use For calibration?",popup,solvents
	DoPrompt "Choose the solvent.", solvent
	
	If(V_Flag)
		Print "Cancelled"
		Return 0 //user canceled
	EndIf
	
	Switch(solvent)
		case 1:
			shift[0]= {383.81,426.62,801.39,1028.13,1157.64,1265.91,1347.57,1443.73} //CHEX
			break
		case 2:
			shift[0]= {309,333,382,667,698,953,1026,1042,1307,1417} //DMSO
			break
		case 3:
			shift[0]= {435,486.9,834.8,1015.2,1109.1,1127.5,1217.1,1304.6,1443.9}//Dioxane
			break
		case 4:
			shift[0]={433.5,883.3,1051.6,1095.2,1275.6,1453.7}//EtOH
			break
		case 5:
			shift[0]= {605.6,991.6,1178,1584.6,1606.4}//Benzene
			break
	endswitch
	
	//Resize the wave and form the pixel wave
	WaveTransform/O zapNaNs  shift
	Duplicate/O shift pixel
	pixel=NaN
	
	//Displays the pixel and shift waves For editing
	Edit/N=PixelShift pixel shift as "Pixel and Shift"
	
	//sets the max number of interations For curvefitting
	Variable/G V_FitMaxIters=200
	
	//Makes the SpeedOfLight constant available globally so that it can be used on the command line
	Variable/G SpeedOfLight = kSpeedOfLight
End

Function Calibrate(pixel,shift,laser,ranwave)
	//This will perform a linear fit to wavelength and calculated the wavelength
	// and wavenumber axes when given a calibration set of the waves PIXEL and SHIFT,
	//which are the known pixels and Raman shifts of one's calibrated solvent.  Also,
	//you need to include the laser wavelength in nm, LASER.
	
	WAVE pixel		//The pixel locations of those shifts
	WAVE shift		//The raman shifts of the calibration standard
	Variable laser		//The Raman laser wavelength, in nanometers
	WAVE ranwave	//Any spectrum, we just need to know how many points long the shift wave should be
	
	Variable length = numpnts(ranwave)
	
	Make/D/O/N=(numpnts(pixel)) calib_wavelength,res_shift
	
	calib_wavelength=1E7/(1E7/laser -shift)
	
	//Set up the graph
	Display/N=Calibration calib_wavelength vs pixel as "Shift Calibration"
	Label left "Wavelength (nm)"
	Label bottom "Pixel"
	ModifyGraph mirror=2
	ModifyGraph mode=3,marker=8
	
	//Perform a linear least squares fit
	CurveFit/Q line calib_wavelength /X=pixel /D
	
	//Declare references to the just created waves
	WAVE w_coef, w_sigma
	
	//Store these as global variables for easy access in the future
	Variable/G wlDeltax = W_coef[1]
	Variable/G wlLeftx = W_coef[0]
	
	//Color the fit line
	ModifyGraph rgb(fit_calib_wavelength)=(0,0,65535)
 	
 	//Make a residual shift wave
 	res_shift = shift-1e7*(1/laser - 1/(W_coef[1]*pixel+W_coef[0]))
 	WaveStats/Q res_shift
 	
 	//Add annotations to the graph
 	TextBox/C/N=text0/A=LT/X=0/Y=0  "Wavelength Calibration\rwavelength = a*pixel+b\ra = "+num2str(W_coef[1])+"±"+num2str(W_sigma[1])+"\rb = "+num2str(W_coef[0])+"±"+num2str(W_sigma[0])
	TextBox/C/N=text1/A=RB/X=0/Y=0 "Std Dev = " +num2str(V_sdev)+" cm\\S-1\\M\rR\S2\M = "+num2str(V_r2)
	
	//Create the shiftx and nmshiftx waves
	Make/D/O/N=(length) shiftx=1e7*((1/laser)-(1/(W_coef[1]*x+W_coef[0])))
	Make/D/O/N=(length+1) shiftxForImages=1e7*((1/laser)-(1/(W_coef[1]*x+W_coef[0])))
	Make/D/O/N=(length) nmshiftx=W_coef[0]+W_coef[1]*x
	
	KillWaves res_shift	//Clean up
	Return 0
End

Function TransformLimit()
//A small gui that will calculate the transform limit of a gaussian shaped pulse
	Variable tl = 2*Ln(2)/Pi//=0.441...

	Variable choice=1
	Variable Num=10
	
	Prompt choice,"Is the number in wavenumbers or fs?",popup,"Wavenumbers;femtoseconds"
	Prompt Num, "Enter the number to convert: "
	
	DoPrompt "Compute the transform limit",Num, choice
	
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	If(choice==1) // the user chose to convert from cm  to fs
		Print "A transform limited pulse with a bandwidth of "+num2str(Num)+" cm-1 has a duration of "+num2str(tl/kSpeedOfLight/num)+" fs"
	Else  // the user chose to convert from fs to cm
		Print "A transform limited pulse with a duration of "+num2str(Num)+" fs has a bandwidth of "+num2str(tl/kSpeedOfLight/num)+" cm\\S-1\\M"
	EndIf
End 

Function CalcFWHM(sigma)
//Convert the gaussian width (defined by Igor's fit function) to the full width at half maximum
	Variable sigma
	return 2*sqrt(ln(2))*sigma
End

Function/S ReverseList(list)
//A simple utility that takes an Igor list string and reverses the order
	String list
	
	String newList = ""
	Variable numItems = ItemsInList(list)
	
	Variable i=0
	For(i=0;i<numitems;i+=1)
		newList += StringFromList((numitems-i-1),list)+";"
	EndFor
	
	Return newList
End

Function DisplayWL(wl,[XWave])
	//A quick macro to display a list of waves
	String wl //the list to be displayed
	WAVE XWave //optional abscissa WAVE
	
	//Pull the length
	Variable length=ItemsInList(wl)
	
	If(length==0)
		//Error checking, in case the user has made a mistake
		Print "There was nothing to display."
		Return -1
	EndIf
	
	String TN=""//Trace name
	
	Display //set up the display window
	
	Variable i=0 //my iterator
	
	For(i=0;i<length;i+=1)
		//loop through the list pulling out each WAVE
		TN = StringFromList(i,wl, ";")
		If(WaveExists($TN))//Make sure the WAVE exists!
			//append the waves to the top graph, with or without an xwave as desired
			If(ParamIsDefault(XWave))
				AppendToGraph $TN
			Else
				AppendToGraph $TN vs XWave
			EndIf
		Else
			Print TN+" does not exist!"
		EndIf
	EndFor
	
	Return 0 //successful exit
End

Function autoBaseline()
	//A function to speed up baseline subtraction.
	//It will automatically find the withg WAVE and the polynomial and run
	//baseline_sub
	WAVE shiftx=root:shiftx
	
	If(!WaveExists(shiftx))
		Beep
		Print "You Raman Shift WAVE doesn't even exist yet!"
		Return 1
	EndIf
	
	String namesOfPoly = WaveList("w_*poly*",";","WIN:")//Get the names of the poly's in the top window
	String nameOfWithg=  WaveList("*withg",";","WIN:")//withg WAVE ***NOTE: CHANGED TO SUBG***
	String waveX, waveY,withg
	
	withg=StringFromList(0,nameOfWithg) //pull out withg
	waveY=StringFromList(0,namesOfPoly) //pull out w_ypoly
	wavex=StringFromList(1,namesOfPoly) //pull out w_xpoly
	
	Variable polyNum=0//Poly number
	
	sscanf LowerStr(waveY), "w_ypoly%u", polyNum//pull out polynum For use in baseline_sub
	
	If(WaveExists($withg) && WaveExists($waveY) && WaveExists($waveX))
		//Run baselin_sub
		baseline_sub($withg,shiftx,polynum)
	Else
		//Tell the user they made a mistake
		Beep
		Print "The proper waves do not exist in the top graph!"
		Return -1
	EndIf
	
	Return 0
End

function XCMaxG(Matrix)
	//This function takes a Matrix (presumably the DetXC) and fits each column to a gaussian
	//which is an approximation of what the cross correlation should look like, and then saves the center
	//frequency and width to two separate appropriately named waves.
	
	WAVE Matrix
	
	Variable V_FitError=0
	Variable i=0,length = DimSize(Matrix,0)//Length is the number of Pixels (1340 For the PIXIS)
	
	//Make my waves
	Make/D/O/N=(length) $(NameOfWave(Matrix)+"_Max")=NaN, $(NameOfWave(Matrix)+"_Width")=NaN
	Make/D/O/N=(length) $(NameOfWave(Matrix)+"_MaxS")=NaN, $(NameOfWave(Matrix)+"_WidthS")=NaN
	
	//Make the waves available to the function
	WAVE tempMax = $(NameOfWave(Matrix)+"_Max")
	WAVE tempWidth = $(NameOfWave(Matrix)+"_Width")
	WAVE tempMaxS = $(NameOfWave(Matrix)+"_MaxS")
	WAVE tempWidthS = $(NameOfWave(Matrix)+"_WidthS")
	
	Variable timerRefNum = startMSTimer
	
	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
	
	For(i=0;i<length;i+=1)
		
		If(mod(i,20)==0)
			PrintF "*"
		EndIf
		
		V_fitError = 0
		//Fit the ith column
		CurveFit/M=0/Q/N=1/W=2 gauss Matrix[i][*]
		
		//These declarations insures that the two waves exist, really only important If this macro is run
		//before W_Coef and W_Sigma are created.
		WAVE w_coef
		WAVE w_sigma
		//print w_coef
		If(!V_fiterror && (abs(w_coef[3]) < 400) &&  (abs(w_coef[3]) > 10) && w_sigma[2] < 10 && w_sigma[3] < 10)
			//Fill the results waves
			tempMax[i]=w_coef[2]
			tempWidth[i]=abs(w_coef[3])
			tempMaxS[i]=w_sigma[2]
			tempWidthS[i]=w_sigma[3]
		Else
			//With NaN's If that makes the most sense
			tempMax[i]=NaN
			tempWidth[i]=NaN
			tempMaxS[i]=NaN
			tempWidthS[i]=NaN
		EndIf
	EndFor
	
	Print "/-"
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	Return 0
	
End

function XCMaxK(Matrix)
	//This function takes a Matrix (presumably the DetXC) and fits each column to a conv_exp
	//which is an approximation of what the xc should look like, and then saves the center
	//frequency and width to two separate appropriately named waves.
	
	WAVE Matrix
	
	Variable V_FitError=0
	Variable i=0,length = DimSize(Matrix,0)//Length is the number of Pixels (1340 For the PIXIS)
	
	//Make my waves
	Make/D/O/N=(length) $(NameOfWave(Matrix)+"_Max")=NaN, $(NameOfWave(Matrix)+"_Width")=NaN
	Make/D/O/N=(length) $(NameOfWave(Matrix)+"_MaxS")=NaN, $(NameOfWave(Matrix)+"_WidthS")=NaN
	
	//Make the waves available to the function
	WAVE tempMax = $(NameOfWave(Matrix)+"_Max")
	WAVE tempWidth = $(NameOfWave(Matrix)+"_Width")
	WAVE tempMaxS = $(NameOfWave(Matrix)+"_MaxS")
	WAVE tempWidthS = $(NameOfWave(Matrix)+"_WidthS")
	
	Make/FREE/D/O/N=6 kerr_coefs
	
	Variable timerRefNum = startMSTimer
	
	Printf "Processing "+NameOfWave(Matrix)+":\t -/"
	
	V_fitError = -1
	
	For(i=0;i<length;i+=1)
		
		If(mod(i,20)==0)
			PrintF "*"
		EndIf
		
		//tempFit = Matrix[i][p]//Assign column i to tempFit
		
		If(V_FitError)
			//Makes some guesses based off a gaussian fit,
			//Only make new guesses If the previous fit failed
			CurveFit/M=0/Q/N=1/W=2 gauss Matrix[i][*]
			WAVE w_coef
			kerr_coefs={w_coef[3],w_coef[2],w_coef[0],w_coef[1]*0.1,w_coef[3]*2,w_coef[1]}
		EndIf
		
		V_FitError=0
		FuncFit/M=0/Q/N=1/W=2 kerr kerr_coefs Matrix[i][*]
		
		
		//These declarations insures that the two waves exist, really only important If this macro is run
		//before W_Coef and W_Sigma are created.
		WAVE w_sigma
		//print w_coef
		If(!V_fiterror && (w_sigma[1]<10) &&  (abs(Kerr_coefs[0]) > 30))
			//Fill the results waves
			tempMax[i]=kerr_coefs[1]
			tempWidth[i]=abs(kerr_coefs[0])
			tempMaxS[i]=w_sigma[1]
			tempWidthS[i]=w_sigma[0]
		Else
			//With NaN's If that makes the most sense
			tempMax[i]=NaN
			tempWidth[i]=NaN
			tempMaxS[i]=NaN
			tempWidthS[i]=NaN
		EndIf
	EndFor
	
	Print "/-"
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Processing time: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	Return 0
	
End

Function XCMaxN(Matrix)
	//Same as XCMaxG except it uses the FindPeak function to find the center and width at each pixel
	
	WAVE Matrix
	
	Variable i=0,length = DimSize(Matrix,0)//Length is the number of Pixels (1340 For the PIXIS)
	
	//Make my waves
	Make/D/O/N=(length) $(NameOfWave(Matrix)+"_Max")=NaN, $(NameOfWave(Matrix)+"_Width")=NaN
	Make/D/O/N=(length) $(NameOfWave(Matrix)+"_MaxS")=NaN, $(NameOfWave(Matrix)+"_WidthS")=NaN
	
	//Make the waves available to the function
	WAVE tempMax = $(NameOfWave(Matrix)+"_Max")
	WAVE tempWidth = $(NameOfWave(Matrix)+"_Width")
	WAVE tempMaxS = $(NameOfWave(Matrix)+"_MaxS")
	WAVE tempWidthS = $(NameOfWave(Matrix)+"_WidthS")
	
	Make/FREE/D/O/N=(DimSize(Matrix,1)) tempFit//This WAVE is as long as there are timepoints
	
	SetScale/P x DimOffset(Matrix,1),DimDelta(Matrix,1),"", tempFit//Set the scale so the fit parameters
	//mean something useful
	Variable V_PeakLoc=NaN, V_PeakWidth=NaN,V_PeakVal=0,V_sdev,V_avg
	
	For(i=0;i<length;i+=1)
		tempFit = Matrix[i][p]//Assign column i to tempFit
		WaveStats/Q tempFit
		If(V_avg>V_sdev)
			FindPeak/Q/M=(V_avg) tempFit
		Else
			FindPeak/Q/M=(V_sdev) tempFit
		EndIf
		//Fill the results waves
		//If(V_PeakVal>mean(tempFit)
		tempMax[i]=V_PeakLoc
		tempWidth[i]=V_PeakWidth/sqrt(2)
	EndFor
End

Function/S AverageWaves(baseName,avgName,[q])
//this function averages multiple waves with the same base name
//Igor has a built in package that is much more general but also more
//difficult to use
	String baseName		//Base string, make sure it contains wild cards
	String avgName		//The name you want For the result
	Variable q				//Quiet! Do you want history output or not
	
	If(ParamIsDefault(q))
		q = 1
	EndIf
	
	String wl = WaveList(baseName,";","")	//a list of waves matching the base name is generated
	Variable length =ItemsInList(wl)		//Number of Items in the list
	
	Variable V_numNans=0
	
	String toPrint = ""	//We will build up a printout string as we go instead of printing to the history directly
	
	String tn = ""	//Variables For the function
	If (length !=0)
		
		If(q)
			toPrint+= "Beginning averaging waves...\r"
		EndIf
		
		//Creates a WAVE that is the avg output WAVE
		Duplicate/O $StringFromList(0,wl, ";") $avgName, $(avgName+"_SDW")
		
		WAVE avg =$avgName					//create a WAVE reference to the just created WAVE
		WAVE SDW =$(avgName+"_SDW")	//create a WAVE reference to the just created WAVE
		
		avg=0
		SDW=0
		
		Variable i = 0,j=0
		For(i=0;i<length;i+=1)
			tn = StringFromList(i,wl, ";")	//Pull the next string name
			
			WaveStats/Q $tn					//Check to see if there are any NaNs			
			If(V_numNaNs==0)
				If(q)
					toPrint += "  Adding : " + tn+"\r"
				EndIf
				WAVE source = $tn
				avg += source
				SDW +=source^2
			Else
				If(q)
					toPrint+= "  "+tn+" had some NaNs and will NOT be included in the average\r"
				EndIf
				j+=1
			EndIf
		EndFor
		
		i-=j	//The number of waves ACTUALLY included in the average
		
		If(q)
			toPrint+= "  There were "+num2str(j)+" waves not included in the average\r"
			toPrint +="  Dividing by number of waves: "+num2istr(i)+"\r"
		EndIf
		
		avg /= i										//divide by number of waves
		SDW=Sqrt(SDW/(i-1)-i/(i-1)*avg^2)	//calculated the standard deviation
		
		If(q)
			toPrint += "  Done averaging waves!\r"
			Print toPrint
		EndIf
		
		Return GetWavesDataFolder(avg,2)	// return the full path to the wave
	Else
		DoAlert 0, "No waves could be found that match \"" + baseName + "\""
		Return ""
	EndIf
End

Function/S FindPeakLocations(spectrum, [minamp, name,plot])
//A function For automatically finding peaks in a spectrum
//Useful For determining the peak locations in calibration spectra

	WAVE Spectrum	//The spectrum in which to find peaks
	Variable minamp	//User definable minimum amplitude For found peaks
	Variable plot		//Do you want to plot the results
	String name		//User definable name For the peak locations
	
	If(ParamIsDefault(minamp))
		//If the user doesn't choose a minimum amplitude, choose one
		minamp = mean(spectrum)
		Print "Threshold = "+num2str(minamp)
	EndIf
	
	If(ParamIsDefault(plot))
		//The default is to plot
		plot=1
	EndIf
	
	If(ParamIsDefault(name))
		name = ""
	EndIf
	
	Variable V_Flag,V_PeakWidth,V_PeakLoc	//Make sure we have our variables ready to use
	Make/D/O/N=100 PeakLocations=NaN		//A wave in which to store the peak locations (pixel)
	
	Variable i=0	//Iterator
	do
		//Try to find a peak between the last known location and the end of the spectrum
		FindPeak/M=(minamp)/Q/R=[V_PeakLoc,] spectrum
		//If there was an error finding a peak or we reached the end of the spectrum then exit
		If(V_Flag!=0||V_PeakLoc>=numpnts(spectrum))
			Break
		Elseif(V_PeakWidth>4&&V_PeakWidth<60)
			PeakLocations[i]=V_PeakLoc
		EndIf
		V_PeakLoc+=1
		i+=1
	While(1)
	
	//Remove the unused locations
	WaveTransform/O zapNaNs peaklocations
	
	//If the user has selected a name than use it
	If(strlen(name)!=0)
		Duplicate/O PeakLocations $name
		KillWaves PeakLocations
	Else
		If(WinType("PeakLocations0")==0)//The window does NOT already exist
			Edit/N=PeakLocations peaklocations as "Peak Locations"
		EndIf
		//Switch name to default For later use
		Name = "PeakLocations"
	EndIf
	
	//Plot the results if desired 
	If(plot)
		Display spectrum as NameOfWave(Spectrum)+"'s Peak Locations"
		
		//A cheap way to visualize the peak locations
		AppendToGraph/L=peaks $name vs $name
		ModifyGraph mode($name)=1,rgb($name)=(0,0,0)
		
		SetAxis peaks 0,1
		ModifyGraph tick(peaks)=3,noLabel(peaks)=2,freePos(peaks)=0
		ModifyGraph mirror=2,standoff=0
		Legend/C/N=text0/J "\\s("+NameOfWave(Spectrum)+") Spectrum\r\\s(PeakLocations) Peak Locations"
	EndIf
	
	Return name
End

Function/S makeTimeWL(timepoints, prefix, suffix)
//A helper function that takes a list of times and returns a String list that can be manipulated
//using the StringFromList built in function
	WAVE timepoints		//Time points wave
	String prefix, suffix	//The prefix and suffix to add to the formated time string

	String wl = ""	//Initialize the wavelist string
	Variable i, length = numpnts(timepoints)
	
	//Loop throught the timepoints wave and format it in the appropriate manner
	For(i=0;i<length;i+=1)
		//Append the prefix and suffix and a semi colon
		wl += prefix + myTime(timepoints[i]) + suffix +";"
	EndFor
	
	//Return the list
	Return wl
End

Function/S AvgTimeWaves(timepoints, mainBase,[base1,base2,q])
//Averages all waves For a set of time delays
	WAVE timepoints	//must be a WAVE of the time points {-100,0,200...}
	String mainBase	//The mainbase is a string with wild cards that will be used to match wave names
	String base1		//base1 is the suffix added after the timedelay
	String base2		//base2 is the suffix added to the output wave name
	Variable q			//q, how much output to print
	
	String testStr
	String wl = ""
	
	//The optional parameter lets the user specify what he would like to be averages
	//In the default case it will be the excited state Raman gain spectra
	If(ParamIsDefault(base1))
		base1="exc*_spec"
	EndIf
	
	//We'll just append "_avg" to distinguish the average from the input
	If(ParamIsDefault(base2))
		base2="_avg"
	EndIf
	
	If(ParamIsDefault(q))
		q = 1	//Print out
	EndIf
	
	Variable i, length = numpnts(timepoints)
	String currenttime
	
	For(i=0;i<length;i+=1)
		currenttime = myTime(timepoints[i])
		testStr = AverageWaves(mainBase+currenttime+base1,currenttime+base2,q=q) //averages the current time point
		If (strlen(testStr) ==0)
			DoAlert 0, "You saw the previous alert, you probably mucked something up!"
			Return ""
		EndIf
		
		If(q)
			Print testStr
		EndIf
		
		wl += testStr +";"
	EndFor
	
	Return wl
End

Function/S myTime(myTime)
	Variable myTime
	
	If(myTime>0)
		Return "p"+num2istr(abs(myTime))
	Else
		Return "m"+num2istr(abs(myTime))
	EndIf
End

Function/S SpectraSubtract(spectrum1, spectrum2, peak,base,name)
//Removes one spectrum from another by scaling spectrum2 by a scale factor
//that's calculated from a common feature (such as a solvent line)
	WAVE spectrum1	//The main spectrum
	WAVE spectrum2	//The spectrum to be removed from spectrum1
	Variable peak		//The point at which the common feature is maximum
	Variable base		//A point which is the baseline
	String name		//The name of the wave in which the result of the subtraction should be placed
	
	//Make the wave to be returned
	Duplicate/O spectrum1 $name
	WAVE sub = $name
	
	//Calculated the scale factors
	Variable/G SpectraSubtractScale =(spectrum1[peak]-spectrum1[base])/(spectrum2[peak]-spectrum2[base])
	Variable/G SpectraSubtractScale_SD=0
	
	//Do the subtraction
	sub -= Spectrum2*SpectraSubtractScale
	
	//Time to propogate the error
	WAVE Spectrum1_SDW = $(NameOfWave(Spectrum1)+"_SDW")
	WAVE Spectrum2_SDW = $(NameOfWave(Spectrum2)+"_SDW")
	
	//Make sure the error waves exist, to ensure backwards compatibility.
	If(WaveExists(Spectrum1_SDW) && WaveExists(Spectrum2_SDW))
		Duplicate/O Spectrum1_SDW $(name+"_SDW")
		WAVE sub_SDW = $(name+"_SDW")
		
		//Propogate error
		SpectraSubtractScale_SD = (spectrum1_SDW[peak]/(spectrum2[peak]-spectrum2[base]))^2
		SpectraSubtractScale_SD += (spectrum1_SDW[base]/(spectrum2[peak]-spectrum2[base]))^2
		SpectraSubtractScale_SD += (spectrum2_SDW[peak]*(spectrum1[peak]-spectrum1[base])/(spectrum2[peak]-spectrum2[base])^2)^2
		SpectraSubtractScale_SD += (spectrum2_SDW[base]*(spectrum1[peak]-spectrum1[base])/(spectrum2[peak]-spectrum2[base])^2)^2
		
		SpectraSubtractScale_SD = Sqrt(SpectraSubtractScale_SD)
		
		sub_SDW = Sqrt(Spectrum1_SDW^2+(Spectrum2_SDW*SpectraSubtractScale)^2+(Spectrum2*SpectraSubtractScale_SD)^2)
	EndIf
	
	Return GetWavesDataFolder(sub,2)	//Return the full path to the result
End

Function/S SpectraSubtract2(spectrum1, spectrum2, startpt,endpt,name,type)
//Remove the spectrum2 from spectrum1 by fitting the specified peak to a
//lineshape function with a cubic baseline to accurately extract the height this is 
//done For both spectra in order to generate the appropriate scale factor.
	WAVE spectrum1	//Spectrum from which to subtract
	WAVE spectrum2	//this spectrum
	Variable startpt	//Defines the start of the fit region
	Variable endpt		//Defines the end of the fit region
	String name		//Name of the subtracted wave
	Variable type		//1 for gaussian 0 for lorentzian
	
	WAVE Spectrum1_SDW = $(NameOfWave(Spectrum1)+"_SDW")
	WAVE Spectrum2_SDW = $(NameOfWave(Spectrum2)+"_SDW")
	
	If(!WaveExists(Spectrum1_SDW) || !WaveExists(Spectrum2_SDW))
		Print "There are no standard deviation waves, please run average waves again"
		Return ""
	EndIf
	
	String myResults = SinglePeakArea(spectrum1,startpt,endpt,type)
	
	Variable Spec1Amp = str2num(StringFromList(0,myResults))
	Variable Spec1Amp_SD = str2num(StringFromList(1,myResults))
	
	myResults = SinglePeakArea(spectrum2,startpt,endpt,type)
	
	Variable Spec2Amp = str2num(StringFromList(0,myResults))
	Variable Spec2Amp_SD = str2num(StringFromList(1,myResults))
	
	//now make the sub WAVE
	Duplicate/O spectrum1 $name
	WAVE sub = $name
	
	//Make the scale factor
	Variable Scale = Spec1Amp/Spec2Amp
	
	//Do the subtraction
	sub -= Spectrum2*Scale
	
	//Error propagation on the scale factor
	Variable SpectraSubtract2Scale_SD = Sqrt((Spec1Amp_SD/Spec2Amp)^2+(Spec1Amp*Spec2Amp_SD/(Spec2Amp^2))^2)
	
	//Time to propogate the error on the actual waves
	Duplicate/O Spectrum1_SDW $(name+"_SDW")
	WAVE sub_SDW = $(name+"_SDW")
		
	sub_SDW = Sqrt(Spectrum1_SDW^2+(Spectrum2_SDW*Scale)^2+(Spectrum2*SpectraSubtract2Scale_SD)^2)
	
	String toReturn
	sprintf toReturn, "%.16g;%.16g;",Scale,SpectraSubtract2Scale_SD
	
	//return the path to the solvent removed WAVE
	Return toReturn
End

Function/S SinglePeakArea(spectrum,sp,ep,type)
	WAVE Spectrum	//The spectrum to fit
	Variable sp		//The starting point of the fit
	Variable ep		//The end point of the fit
	Variable type		//1 for gaussian 0 for lorenztian
	
	Variable PolyNum = 3
	
	//Generate guesses
	If(type)
		CurveFit/O/W=2/Q/N gauss spectrum[sp,ep]
	Else
		CurveFit/O/W=2/Q/N lor spectrum[sp,ep]
	EndIf
	
	Wave W_coef
	
	W_coef[0]=0	//Make sure the offset is set to zero, we will be including a baseline later
	
	Duplicate/O w_coef tempPeak_Coefs
	Make/D/O/N=(PolyNum+1) tempBaseln_Coefs=0
	
	//Set up the fit with the lineshape and the polynomial baseline
	String F_String=""
	If(type)
		F_String="{gauss, tempPeak_Coefs, hold=\"1\"}"
	Else
		F_String="{lor, tempPeak_Coefs, hold=\"1\"}"
	EndIf
	
	F_String+="{poly_XOffset "+num2str(PolyNum)+", tempBaseln_Coefs}"
	
	Variable V_FitError = 0
	Variable V_FitMaxIters=200
	//Do the fit
	
	FuncFit/W=2/Q/N {string = F_String} Spectrum[sp,ep]
	WAVE W_sigma
	
	If(V_FitError)
		//There was a problem with the fitting
		Return "Error"
	EndIf
	
	String toReturn = ""
	
	If(type)
		sprintf toReturn, "%.16g;%.16g;",(tempPeak_Coefs[1]*sqrt(pi)*tempPeak_Coefs[3]),sqrt(pi)*sqrt(w_sigma[1]^2+w_sigma[3]^2)
	Else
		sprintf toReturn, "%.16g;%.16g;",(tempPeak_Coefs[1]*pi/sqrt(tempPeak_Coefs[3])),pi/2*sqrt((4*tempPeak_Coefs[3]^2*w_sigma[1]^2+tempPeak_Coefs[1]^2*w_sigma[3]^2)/tempPeak_Coefs[3]^3)
	EndIf
	
	//I'm using strings as a data structure
	Return toReturn //The first position is the area, and the second is the error on the area
End

Function/S GroundSubtract1to1(timepoints,ground)
//subtracts the ground spectrum from all of the timepoints (1-to-1 subtraction)
//times points and ground should be averaged before using this function
//*****NOTE: NO NORMALIZATION IS DONE HERE****

	WAVE timepoints	//Time delays
	WAVE ground		//Ground state spectrum
	
	Variable i, length = numpnts(timepoints)
	String currenttime, wl = ""
	
	String tracename, newname
	
	//Loop through the time delays
	For(i=0;i<length;i+=1)
		currenttime = myTime(timepoints[i])
		//Some error checking, ensure that the waves exist
		If(!WaveExists($(currenttime+"_avg")))
			DoAlert/T="GroundSubtract Failed!" 0, currenttime+"_avg does not exist!"
			Return ""
		EndIf
		
		//Make the ground subtracted wave
		Duplicate/o $(currenttime+"_avg") $(currenttime+"_subG")
		WAVE sub = $(currenttime+"_subG")
		sub -= ground	//Do the subtraction
		wl += GetWavesDataFolder(sub,2) + ";" //Append the full path the the output
	EndFor
	
	//Return a list of ground subtracted waves
	Return wl
End

STATIC Function GroundSubGrapher()
	//Checks to see If the groundsubtracted WAVE is already displayed
	//If not then it displays it.
	If(WinType("GNDSubGraph")==0)//The window does NOT already exist
		//Here we make some assumptions about what already exists in the pxp
		//and how its formated.
		WAVE groundsubtracted = root:groundsubtracted
		WAVE timepoints = root:timepoints
		WAVE groundsubtracted_SDW = groundsubtracted_SDW
		//Display
		Display/N=GNDSubGraph groundsubtracted vs timepoints
		ModifyGraph mirror=2,minor(bottom)=1,prescaleExp(bottom)=-3
		Label left "Scale Factor"
		Label bottom "Delay (ps)"
		ErrorBars groundsubtracted Y,WAVE=(groundsubtracted_SDW,groundsubtracted_SDW)
		ModifyGraph mode=4,marker=8
	EndIf
End

Function/S GroundSubtract(timepoints,ground,peak,base,[q])
//A function that loops through all the time delays removing the ground
//spectrum from each one using SpectraSubtract
	WAVE timepoints	//Time delays
	WAVE ground		//Ground state spectrum
	Variable peak		//Location of the peak
	Variable base		//Location of the baseline
	Variable q			//Verbose or not

	If(ParamIsDefault(q))
		q=0	//Not verbose
	EndIf
	
	Variable i, length,scale
	String currenttime, wl = ""
	
	//Get some global variables set up so that we can pass variables from
	//SpectraSubtract to this function without using structures.
	Variable/G SpectraSubtractScale
	Variable/G SpectraSubtractScale_SD
	
	//Make a wave to store the scale factors
	Make/D/O/N=(numpnts(timepoints)) groundsubtracted, groundsubtracted_SDW
	
	length = numpnts(timepoints)
	For(i=0;i<length;i+=1)
		
		currenttime = myTime(timepoints[i])
		
		If(!WaveExists($(currenttime+"_avg")))
			DoAlert/T="GroundSubtract Failed!" 0, currenttime+"_avg does not exist!"
			Return ""
		EndIf
		
		//Build up out list of ground subtracted spectra
		wl += SpectraSubtract($(currenttime+"_avg"),ground,peak, base,(currenttime+"_subg")) + ";"
		
		//Update our scale factor waves
		groundsubtracted[i]=SpectraSubtractScale
		groundsubtracted_SDW[i]=SpectraSubtractScale_SD
		
		//Print out If wanted...
		If(q)
			If(abs(timepoints[i]) < 1000)
				Print currenttime+"\t\t"+num2str(SpectraSubtractScale)
			Else
				Print currenttime+"\t"+num2str(SpectraSubtractScale)
			EndIf
		EndIf
	EndFor
	
	//Graph our scale waves
	GroundSubGrapher()
	
	Return wl
End

Function/S GroundSubtract2(timepoints,ground,startpt,endpt,type,[q])
//A function that loops through all the time delays removing the ground
//spectrum from each one using SpectraSubtract2
	WAVE timepoints	//Time delays
	WAVE ground		//Ground state spectrum
	Variable startpt	//The start of the subrange which to fit
	Variable endpt		//The end of the subrange which to fit
	Variable type		//1 for gaussian 0 for lorenztian
	Variable q			//Verbose?
	
	If(ParamIsDefault(q))
		q=0 //Not verbose
	EndIf
	
	Variable i, length, scale, scale_sd
	String name, currenttime, wl = ""
	
	WAVE ground_SDW = $(NameOfWave(ground)+"_SDW")
	
	If(!WaveExists(ground_SDW))
		Print "There are no standard deviation wave for the gnd, please run average waves again"
		Return ""
	EndIf
	
	String excResults, gndResults = SinglePeakArea(ground,startpt,endpt,type)
	
	Variable excAmp, gndAmp = str2num(StringFromList(0,gndResults))
	Variable excAmp_SD, gndAmp_SD = str2num(StringFromList(1,gndResults))
	
	Make/D/O/N=(numpnts(timepoints)) groundsubtracted, groundsubtracted_SDW
	
	length = numpnts(timepoints)
	For(i=0;i<length;i+=1)
		
		currenttime = myTime(timepoints[i])
		
		If(!WaveExists($(currenttime+"_avg")))
			DoAlert/T="GroundSubtract2 Failed!" 0, currenttime+"_avg does not exist!"
			Return ""
		EndIf
	
		WAVE exc = $(currenttime+"_avg")
		name = currenttime+"_subg"
		
		WAVE EXC_SDW = $(currenttime+"_avg_SDW")
		
		excResults = SinglePeakArea(exc,startpt,endpt,type)
		
		excAmp = str2num(StringFromList(0,excResults))
		excAmp_SD = str2num(StringFromList(1,excResults))
		
		//now make the sub WAVE
		Duplicate/O exc $(currenttime+"_subg")
		WAVE sub = $name
		
		//Make the scale factor
		//***NOTE: Here we are scaling the excited state to the ground state***//
		scale = gndAmp/excAmp
		
		//Do the subtraction
		//The ground state is acting as our "standard" for Raman pump power"
		sub = scale*exc-ground
		
		//Error propagation on the scale factor
		scale_SD = Sqrt((gndAmp_SD/excAmp)^2+(gndAmp*excAmp_SD/(excAmp^2))^2)
		
		//Time to propogate the error on the actual waves
		Duplicate/O sub $(name+"_SDW")
		WAVE sub_SDW = $(name+"_SDW")
			
		sub_SDW = Sqrt(ground_SDW^2+(exc_SDW*Scale)^2+(exc*scale_sd)^2)
		
		//This can be much faster if we don't refit the ground over and over again. Will have to split
		//up spectrasubtract2 into subroutines
		wl += currenttime+"_subg;"
		groundsubtracted[i]=Scale
		groundsubtracted_SDW[i]=Scale_SD
		
		//Print out If wanted...
		If(q)
			If(abs(timepoints[i]) < 1000)
				Print currenttime+"\t\t"+num2str(Scale)
			Else
				Print currenttime+"\t"+num2str(Scale)
			EndIf
		EndIf
	EndFor
	
	GroundSubGrapher()
	
	Return wl
End

Function/S SolventSubtract(spectrum, solvent,peak,base)
//removes the solvent spectrum from the spectrum using SpectraSubtract
	WAVE spectrum		//The spectrum
	WAVE solvent			//The solvent wave
	Variable peak, base	//See SpectraSubtract
	
	Variable/G SpectraSubtractScale
	
	String toReturn = SpectraSubtract(spectrum, solvent, peak,base,(NameOfWave(spectrum)+"_ns"))
	Print SpectraSubtractScale
	
	Return toReturn
End

Function/S SolventSubtract2(spectrum, solvent, sp,ep,type)
//removes the solvent spectrum from the spectrum using SpectraSubtract2
	WAVE spectrum
	WAVE solvent
	Variable sp
	Variable ep
	Variable type	//1 for gaussian 0 for lorenztian
	
	//***NOTE: we are scaling the solvent to the spectrum, this is necessary to be internally consistent***
	String data = SpectraSubtract2(spectrum, solvent, sp,ep,(NameOfWave(spectrum)+"_ns"),type)
	
	Print "Scale factor = "+StringFromList(0,data)
	
	Return NameOfWave(spectrum)+"_ns"
End

Function AddAllGround(timepoints,ground,groundadded,[extraBase])
//adds the specified amount of ground in groundadded to each _subg WAVE
//ground added must be the same length as timepoints and in the same order
	WAVE timepoints	//Time delays
	WAVE ground		//Noiseless ground spectrum
	WAVE groundadded	//Amount of ground to add back at each time points
	String extraBase	//In case the naming structure is different
	
	String currenttime
	
	If(ParamIsDefault(extraBase))
		extraBase = ""
	EndIf
	
	//Loop through the time delays and add the ground
	Variable i, length = numpnts(timepoints)
	For(i=0;i<length;i+=1)
		currenttime = myTime(timepoints[i])
		WAVE tracewave = $(currenttime+extraBase+"_subg")
		If(!WaveExists(tracewave))
			Print currenttime+extraBase+"_subg is not a valid WAVE!"
		Else
			Make/D/O/N=(numpnts(tracewave)) $(currenttime+extraBase+"_withg")=tracewave+groundadded[i]*ground
		EndIf
	EndFor
End

Function AddGround(timepoint,amount,timepoints,ground,groundadded,[extraBase])
//this adds a set amount of ground to particular timepoint and updates groundadded accordingly
	WAVE ground			//Noiseless ground state spectrum
	WAVE timepoints		//Time delays
	WAVE groundadded		//Wave of ground added values
	Variable timepoint	//which time delay?
	Variable amount		//How much ground to add (this is additive to what's in groundadded)
	String extraBase

	If(ParamIsDefault(extraBase))
		extraBase = ""
	EndIf
	
	//Find the index for the given timepoint
	FindValue/v=(timepoint) timepoints
	
	//Do the addition
	String currenttime = myTime(timepoint)+"_withg"
	If(WaveExists($currenttime))
		WAVE withg = $currenttime
		withg=withg+amount*ground
		groundadded[V_value]=groundadded[V_value]+amount
	Else
		Print currenttime + " does not exist, please check your parameters"
	EndIf
End

Function KillWavesMulti(base)
//Allows you to kill waves with a certain base name
	String base	//DON'T FORGET TO INCLUDE THE * and be careful
	
	String wl = WaveList(base,";","")//build list of traces that meet criteria
	
	If (strlen(wl) ==0)
		Print "No waves match that description"
		Return 0
	EndIf
	
	String WaveToKill
	Variable i=0
	//Loop through the wave list, killing each wave
	Do
		WaveToKill = StringFromList(i,wl, ";")
		If (strlen(WaveToKill) ==0)
			break
		EndIf
		i+=1
		KillWaves/Z $WaveToKill
	While(1)
	Print wl
	Print num2str(i)+" waves killed"
	Return 0
End

Function/S baseline_sub(Spectrum,shift,polynum,[type])
//subtracts the specified polynomial after interpolating it relative to shiftx from wave1 with a 3 order spline
	WAVE Spectrum	//The spectrum to be baseline corrected
	WAVE shift		//The xwave used when drawing the baseline
	Variable polynum	//The number of the polynomial that is the baseline
	Variable type		//The type of interpolation to do, see Interpolate2's documentation
	
	If(ParamIsDefault(type))
		//Cubic spline is default
		type=2
	EndIf
	
	//Pull out some useful info
	String num=num2istr(polynum)
	String wn=NameOfWave(Spectrum)
	
	//Make the baseline subtracted wave
	Duplicate/O Spectrum $(wn+"_nb")
	WAVE nobase = $(wn+"_nb")
	
	//Make the baseline wave
	Duplicate/O Spectrum $(wn+"_bl")
	WAVE baseline = $(wn+"_bl")
	
	//Do the interpolation
	Interpolate2/I=3/T=(type)/X=shift/Y=baseline $("W_xpoly"+num),$("W_ypoly"+num)
	
	//Subtract the baseline
	nobase=Spectrum-baseline
	
	//Append the baseline trace to the top graph If Spectrum is present
	String wl=WaveList("*", ";", "WIN:" )
	If(stringmatch(wl,("*"+wn+";*"))&&!stringmatch(wl,("*"+wn+"_bl;*")))
		AppendToGraph baseline vs shift
		ModifyGraph rgb($NameOfWave(baseline))=(2,39321,1)
	EndIf
	
	Return GetWavesDataFolder(nobase,2)// string is full path to baseline subtracted data
End

Function CalcTA(timepoints, mainBase, [q])
//calculates the transient absorption at each time point
	WAVE timepoints	//The time delays
	String mainBase	//The base name for all the data ex. "Chex*"
	Variable q			//verbose ?
	
	If(ParamIsDefault(q))
		//Verbose = yes
		q = 1
	EndIf
	
	//Some timing
	Variable TimerRefNum = StartMSTimer
	
	//Check to see if calibrate has been run, if so we can set the wave scaling in nm
	NVAR/Z myLeftx = wlLeftx
	NVAR/Z myDeltax = wlDeltax
	
	Variable canScale = NVAR_Exists(myLeftx) && NVAR_Exists(myDeltax)
	
	If (!canScale)	// No such global numeric variable?
		Print "Calibrate has not been run yet, waves will not be scaled."
	Endif
	
	//First thing to do is to calculate the ground Raman pump on spectrum
	String pathToGndOn = averagewaves(mainBase+"*gr*pumpon","gnd_pumpon_avg",q=q)
	//Error checking
	If (strlen(pathToGndOn) ==0)
		DoAlert/T="TransA Failed!" 0, "You saw the previous alert, you probably mucked something up!"
		Return -1
	EndIf
	
	//Make my local references
	WAVE GndOn = $pathToGndOn
	WAVE GndOn_SDW = $(pathToGndOn+"_SDW")
	
	//Same thing for Raman pump off
	String pathToGndOff = averagewaves(mainBase+"*gr*pumpoff","gnd_pumpoff_avg",q=q)
	If (strlen(pathToGndOff) ==0)
		DoAlert/T="TransA Failed!" 0, "You saw the previous alert, you probably mucked something up!"
		Return -2
	EndIf
	WAVE GndOff = $pathToGndOff
	WAVE GndOff_SDW = $(pathToGndOff+"_SDW")
	
	//Average the excited state Raman pump on data
	String pumpOnWL = avgtimewaves(timepoints, mainBase,base1="exc*_pumpon",base2="_pumpon_avg",q=q)
	//Error checking
	If (strlen(pumpOnWL) ==0)
		Return -3
	EndIf
	
	//Average the excited state Raman pump off data
	String pumpOffWL = avgtimewaves(timepoints, mainBase,base1="exc*_pumpoff",base2="_pumpoff_avg",q=q)
	If (strlen(pumpOffWL) ==0)
		Return -4
	EndIf
	
	//Time to calculate the transient absorption
	string exWaveNameOn,exWaveNameOff,TAWaveNameOn,TAWaveNameOff,currenttime
	Variable i, length=numpnts(timepoints)
	For(i=0;i<length;i+=1)
		currenttime = myTime(timepoints[i])
		
		//Create wave names
		exWaveNameOn = currenttime+"_pumpon_avg"
		TAWaveNameOn = currenttime+"_TAROn"
		exWaveNameOff = currenttime+"_pumpoff_avg"
		TAWaveNameOff = currenttime+"_TAROff"
		
		//Make the waves
		Make/D/O/N=(numpnts(GndOn)) $TAWaveNameOn=0, $TAWaveNameOff=0
		Make/D/O/N=(numpnts(GndOn)) $(TAWaveNameOn+"_SDW")=0, $(TAWaveNameOff+"_SDW")=0
		
		//Make some local references
		WAVE exWaveOn = $exWaveNameOn
		WAVE exWaveOff = $exWaveNameOff
		
		WAVE exWaveOn_SDW = $(exWaveNameOn+"_SDW")
		WAVE exWaveOff_SDW = $(exWaveNameOff+"_SDW")
		
		WAVE TAROn = $TAWaveNameOn
		WAVE TAROff = $TAWaveNameOff
		
		WAVE TAROn_SDW = $(TAWaveNameOn+"_SDW")
		WAVE TAROff_SDW = $(TAWaveNameOff+"_SDW")
		
		If(q)
			Print "Calculating TA For "+currenttime
		EndIf
		
		//Do the calculation
		TAROn=-log(exWaveOn/GndOn)
		TAROff=-log(exWaveOff/GndOff)
		
		//Propogate the error
		TAROn_SDW=sqrt((exWaveOn_SDW/(exWaveOn*log(10)))^2+(GndOn_SDW/(GndOn*log(10)))^2)
		TAROff_SDW=sqrt((exWaveOff_SDW/(exWaveOff*log(10)))^2+(GndOff_SDW/(GndOff*log(10)))^2)
		
		If(canScale)
			SetScale/P x,myLeftx,myDeltax,"nm", TAROn, TAROff, TAROn_SDW, TAROff_SDW
		EndIF
		
	EndFor
	
	//Create TA matrices for easy contour viewing
	Concatenate/O  makeTimeWL(timepoints, "", "_TAROff"), TAROffMat
	Concatenate/O  makeTimeWL(timepoints, "", "_TAROn"), TAROnMat
	
	Concatenate/O  makeTimeWL(timepoints, "", "_TAROff_SDW"), TAROffMat_SDW
	Concatenate/O  makeTimeWL(timepoints, "", "_TAROn_SDW"), TAROnMat_SDW
	
	//Bin the TA by 100 pixels each
	For(i=0;i<13;i+=1)
		IntegrateMatrix(TAROffMat,"sumTAROff"+num2str(i),sp=100*i+20,ep=100*(i+1)+20)
		IntegrateMatrix(TAROnMat,"sumTAROn"+num2str(i),sp=100*i+20,ep=100*(i+1)+20)
	EndFor
	
	IntegrateMatrix(TAROffMat,"sumTAROff",sp=20,ep=1320)
	IntegrateMatrix(TAROnMat,"sumTAROn",sp=20,ep=1320)
	
	//Display the output
	If(WinType("TotalSumTAGraph")==0)//The window does NOT already exist
		Display/N=TotalSumTAGraph $"sumTAROn" $"sumTAROff" vs timepoints as "Total Sum TA"
		ErrorBars $"sumTAROn" Y,WAVE=($"sumTAROn_SDW",$"sumTAROn_SDW")
		ErrorBars $"sumTAROff" Y,WAVE=($"sumTAROff_SDW",$"sumTAROff_SDW")
		FancifyTA()
		ModifyGraph rgb(sumTAROn)=(0,0,0)
		Legend/C/N=text0/J "TA "+date()+"\r\\s(sumTAROn) Rpu On\r\\s(sumTAROff) Rpu Off"
	EndIf
	
	If(WinType("SumTAROffGraph")==0)//The window does NOT already exist
		Display/N=SumTAROffGraph as "Sum TA Rpu Off"
		For(i=0;i<13;i+=1)
			AppendToGraph $("sumTAROff"+num2str(i)) vs timepoints
			ErrorBars $("sumTAROff"+num2str(i)) Y,WAVE=($("sumTAROff"+num2str(i)+"_SDW"),$("sumTAROff"+num2str(i)+"_SDW"))
		EndFor
		FancifyTA()
		RainbowWaves()
	EndIf
	
	If(WinType("SumTAROnGraph")==0)//The window does NOT already exist
		Display/N=SumTAROnGraph as "Sum TA Rpu On"
		For(i=0;i<13;i+=1)
			AppendToGraph $("sumTAROn"+num2str(i)) vs timepoints
			ErrorBars $("sumTAROn"+num2str(i)) Y,WAVE=($("sumTAROn"+num2str(i)+"_SDW"),$("sumTAROn"+num2str(i)+"_SDW"))
		EndFor
		FancifyTA()
		RainbowWaves()
	EndIf
	
	Variable microSeconds = stopMSTimer(timerRefNum)
	Printf "Time to calculate: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	Return 0
End

Function FourierFilter(spectrum,cutoff_freq)
//applies a low pass filter with cut off cutoff to the specified spectrum
	WAVE spectrum
	Variable cutoff_freq//frequency of noise in your baseline in terms of pixels
	Variable length=numpnts(spectrum)
	Variable cutoff=length/cutoff_freq
	Variable FTlength
	FFT spectrum
	//now can store length of FT spectrum
	FTlength = numpnts(spectrum)
	DeletePoints floor(cutoff),(FTlength-floor(cutoff)), spectrum
	InsertPoints floor(cutoff),(FTlength-floor(cutoff)), spectrum
	IFFT/R spectrum
End

Function FFB(spectrum,cutoff_freq,n)
//applies a Butterworth filter of order n (http://en.wikipedia.org/wiki/Butterworth_filter)
//with cutoff to the specified spectrum
	WAVE/C spectrum//we need to let Igor know that spectrum may be complex
	Variable cutoff_freq,n//frequency of noise in your baseline in terms of pixels
	Variable length=numpnts(spectrum)
	Variable cutoff=length/cutoff_freq//calculating the cutoff in the fourier domain
	Variable FTlength
	//Take the FFT
	FFT spectrum
	spectrum=spectrum*1/(1+(p/(cutoff))^(2*(n)))//apply our filter
	IFFT/R spectrum//transform back to our original domain
End

Function/S FFWaves(base,cutoff_freq,n)
//Apply the same butterworth filter (as defined by n and cutoff_freq
//to all waves matching the specified base name
	String base
	Variable cutoff_freq,n
	Variable i=0
	
	String tracename, wl = WaveList(base,";","")
	
	do //applies the FourerFilter to each WAVE in list
		tracename = StringFromList(i,wl, ";")
		If (strlen(tracename) ==0)
			break
		EndIf
		FFB($tracename,cutoff_freq,n)
		i +=1
	while (1)
	
	Return wl
End

Function/S makeRaman(pos, intens, name, [type, width])
//For taking the frequencies and cross sections from a Gaussian output file and turning
// them into a spectrum with the specified width
	WAVE pos		//The frequencies of the vibrations
	WAVE intens	//The intensities of the vibrations
	String name	//The name of the output wave
	Variable type	//The line shape, 0 for lorentzian, anything else for gaussian
	Variable width	//The FWHM of the line shape convolved with the stick spectrum
	
	//Error checking!
	If(numpnts(pos) != numpnts(intens))
		Print "Error! " + NameOfWave(pos) + " is not the same length as " + NameOfWave(intens)
		Return ""
	EndIf
	
	If(ParamIsDefault(type))
		type = 0	//Lorentzian
	EndIf
	
	If(ParamIsDefault(width))
		//Width is FWHM For either line shape
		width = 7
	EndIf
	
	Variable i, length = numpnts(pos)
	
	//Make our spectrum, 0 to 4000 wavenumbers with a resolution of 1 wavenumber
	Make/D/O/N=4000 $name = 0
	WAVE dest = $name
	
	//Loop through our waves
	For(i=0;i<length;i+=1)
		If(type)
			//Type is anything but 0
			dest += intens[i]*16^(-((p-pos[i])/width)^2) //Gaussian
		Else
			//Type is 0
			dest += intens[i]/((p-pos[i])^2+(width/2)^2)*width^2/4 //Lorentzian
		EndIf
	EndFor
	
	//Make a little report in the WAVE's note
	Note dest, "Command that generated this WAVE: "
	Note/NOCR dest, "makeRaman("+NameOfWave(pos)+", "+NameOfWave(intens)+", \""
	Note/NOCR dest, name+"\", type="+num2str(type)+", width="+num2str(width)+")"
	Note dest, "Details: "
	Note dest, "Peak type = "
	If(type)
		Note/NOCR dest, "Gaussian"
	Else
		Note/NOCR dest, "Lorenzian"
	EndIf
	
	Return GetWavesDataFolder(dest,2) 
End

Function/S makeRINE(Location,Amp,width, shiftx, [scale])
//Takes the location, amplitude and gaussian width given by multipeak fit to simulate
//a RINE spectrum based on Dave McCamant's paper in 2005
	WAVE Location, Amp, width, shiftx
	Variable scale
	If(ParamIsDefault(scale))
		scale = 2 // how much to scale the gaussian widths
	EndIf
	
	Variable length=numpnts(width)
	
	Make/D/O/N=1340 RINE=0
	
	Variable i, g
	For(i=0;i<length;i+=1)
		g=width[i]/scale
		RINE+=Amp[i]*(shiftx-Location[i])*g/((shiftx-Location[i])^2+g^2)
	EndFor
	Return GetWavesDataFolder(RINE,2) 
End

Function/S InvertList(List)
//A simple function that inverts a string list
	String List		//The list to be inverted
	
	Variable length = ItemsInList(List);
	Variable i = 0
	
	String NewList = ""
	For(i=0;i<length;i+=1)
		NewList = StringFromList(i,List)+";"+NewList
	EndFor
	
	Return NewList
End

Function/S IntegrateMatrix(Matrix,name,[sp,ep])
	//This function takes a Matrix and integrates it along the frequency axis
	//between startpnt and endpnt
	
	WAVE Matrix	//The matrix to be integrated
	String name	//The name of the integration
	Variable sp	//Start point
	Variable ep	//End point (note that because we use a for loop this one isn't actually included in the integration)
	
	Variable i=0,length = DimSize(Matrix,0)//Length is the number of Pixels
	
	//Deal with the optional parameters
	If(ParamIsDefault(sp))
		sp=0
	EndIf

	If(ParamIsDefault(ep))
		ep=length
	EndIf
	
	Variable CalcError = 0
	//Does the matrix have an error wave associated with it?
	If(WaveExists($(NameOfWave(Matrix)+"_SDW")))
		CalcError = 1
		Duplicate/O/FREE $(NameOfWave(Matrix)+"_SDW") MatrixVar
		MatrixVar = MatrixVar^2
	EndIf
	
	//Make my waves
	Make/D/O/N=(DimSize(Matrix,1)) $name=0
	//Make the waves available to the function
	WAVE tempInt = $name
	//Transfer the wave scaling from the matrix to the integrated wave
	SetScale/P x DimOffset(Matrix,1),DimDelta(Matrix,1),"", tempInt	
	
	If(CalcError)
		Make/D/O/N=(DimSize(Matrix,1)) $(name+"_SDW")=0
		WAVE tempIntSDW = $(name+"_SDW")
		SetScale/P x DimOffset(Matrix,1),DimDelta(Matrix,1),"", tempIntSDW
	EndIf
	
	For(i=sp;i<ep;i+=1)
		tempInt += Matrix[i][p]//Assign column i to tempFit
		If(CalcError)
			TempIntSDW += MatrixVar[i][p] 
		EndIf
	EndFor
	
	If(calcError)
		TempIntSDW = sqrt(TempIntSDW)
	Endif
	
	//Add a note to the resulting wave
	Note tempInt, "The command that generated this wave is: "
	Note tempInt, "IntegrateMatrix("
	Note/NOCR tempInt, nameofwave(Matrix)+",\""+name+"\""
	Note/NOCR tempInt, ", sp="+num2str(sp)
	Note/NOCR tempInt, ", ep="+num2str(ep)
	Note/NOCR tempInt,")"
	
	Return GetWavesDataFolder(tempInt,2)// string is full path to wave
End

Function Reintegrate(base,sp,ep)
	String base	//Matching string
	Variable sp	//Lower integration limit
	Variable ep	//higher integration limit
	
	String wl = WaveList(base,";","")//a list of waves matching the base name is generated
	
	String tn = ""
	
	Variable i = 0
	Printf "Reintegrating, from %g to %g: ",sp,ep
	do //Loops through and integrates
		tn = StringFromList(i,wl, ";")
		If (strlen(tn) ==0)
			break
		EndIf
		Printf tn
		IntegrateMatrix($tn,"Int_"+tn,sp=sp,ep=ep)
		i+=1
	while (1)
	Printf ".\r"
End

Function/S SVDFilter(Matrix, perToKeep)
//A function that filters a matrix by performing a singular value decomposition
//and keeping only the first singular values whose sum is a certain percentage (perToKeep)
//of the sum of all the singular values.
	Wave Matrix			//The matrix to filter
	Variable perToKeep	//The percentage that defines the cutoff For the singular values
	
	//Set up our output matrix
	Duplicate/O Matrix $(nameofWave(Matrix)+"_smth")
	Wave output = $(nameofWave(Matrix)+"_smth")
	
	MatrixSVD Matrix	//Calculate the SVD of the matrix
	
	//Create local references to our generated waves
	WAVE/C S = W_W
	WAVE/C U = M_U
	WAVE/C VT = M_VT
	
	//Determine the index afterwhich all singular values should be set to zero
	Variable M = numToKeep(S,perToKeep)
	Print "Keeping "+num2str(M)+" singular values"
	
	//Redimension our U and VT waves so that the reconstruction is as fast as possible
	Redimension/N=(-1,M) U
	Redimension/N=(M,-1) VT
	//Generate the output
	MatrixOp/O output = U x DiagRC(S,M,M) x VT
	
	//Clean up
	KillWaves/Z S, U, VT
	
	//Return the path to the output wave
	Return GetWavesDataFolder(Output,2)
End

STATIC Function numToKeep(W,perToKeep)
//A helper function For the above filtering function
//you pass it the singular values as W and the percentage to keep
//as perToKeep, the function then returns the index afterwhich all
//singular values should be set to zero
	Wave W						//A wave containing the singular values from the SVD operation
	Variable perToKeep			//The percentage to keep (0 to 1)
	
	Variable mySum = sum(w)	//the sum of all the singular values.
	Variable x = 0					//A Variable to keep track of things
	
	//Some error checking to keep the output reasonable
	If(perToKeep>1 || perToKeep<0.01)
		return numpnts(w)
	EndIf
	
	Variable i = 0					//My iterator
	For(i=0;i<numpnts(w);i+=1)
		//Iterate through the singular values until x is greater than per to keep
		x+= w[i]/mySum
		If(x>=perToKeep)
			Break
		EndIf
	EndFor
	
	//Return the index
	Return i
End

Function WheresTheNAN(Matrix)
//Sometimes you just want to find the NaN in a wave
	Wave Matrix
	Variable imax = dimsize(matrix,0)
	Variable jmax = dimsize(matrix,1)
	
	Variable i,j
	
	For(i=0;i<imax;i+=1)
		For(j=0;j<jmax;j+=1)
			If(numtype(Matrix[i][j]))
				Print "i = " +num2str(i) + "; j = "+num2str(j)
			EndIf
		EndFor
	EndFor
End