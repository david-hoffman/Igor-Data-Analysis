#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include ":Utilities"		//Make sure that the utilities procedures are included
#include ":PerkinElmerImport"
#include ":Graphing"

//Set up my menus
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

Function LoadRawFSRSData()
//A function that loads raw FSRS data as formated by the FSRS instrumentation software
//for a single (Raman pump) chopper instrument and then calculates the Raman gain spectra as well
//see https://www.github.com/david-hoffman/FSRS-LabVIEW
	
	//The following open command does not actually open the selected file it just returns
	//the selected files paths in the string S_filename for later use
	Variable refNum
	Open/F="????"/D/R/MULT=1/M="Choose Raw FSRS Files ... " refNum
	
	//move the paths so a safer location
	String outputPaths = S_fileName
	
	//check to see If the user cancelled
	If (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	Else
		
		//If the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		//declare some dummy variables
		Variable i,j, k
		
		String path = ""
		String fileName=""
		
		//Loop through the selected files
		for(i=0; i<numFilesSelected; i+=1)
			//Pull a path
			path = StringFromList(i, outputPaths, "\r")
			//Pull the filename from the path
			filename = stringfromlist((ItemsInList(path,":")-1),path,":")
			print path
			
			LoadWave/N/D/L={0,0,1340,10,0}/O/G path//load the file
			
			Print "Done Loading!"
			//Here we process the raw data
			
			j=0		//naming iterator
			
			Printf "Processing - "
			
			For(k=0;k<V_flag;k+=2)
				If(mod(k,100)==0)
					printf "*"
				EndIf
				//Basically, we loop through all the waves loaded renaming them pump on pump off 
				//and creating gain spectra
				Duplicate/O $("WAVE"+num2istr(k)) $(filename+"_pumpoff"+num2istr(j))
				Duplicate/O $("WAVE"+num2istr(k+1)) $(filename+"_pumpon"+num2istr(j))
				
				//Make the WAVE assignment so we can use these later
				WAVE pumpoff = $(filename+"_pumpoff"+num2istr(j))
				WAVE pumpon = $(filename+"_pumpon"+num2istr(j))
				
				//Make the gain spectrum
				Make/D/N=1340/O $(filename+"_spec"+num2istr(j))=ln(pumpon/pumpoff)
				
				//Kill the waves as we go
				KillWaves $("WAVE"+num2istr(k)), $("WAVE"+num2istr(k+1))
			
				//Increment our iterators
				j+=1
				
			EndFor
			Printf " -\r Done Processing!\r"
		EndFor
	EndIf
End

Function LoadUVvisData()
//Author: David Hoffman
//This macro loads data from the lambda 2 spectrophotometer and puts it into a WAVE and then
//changes the scaling to the appropriate amount, i.e. there is no need to load the data as an X-Y pair
	Variable begin,step
	//The following open command does not actually open the selected file
	// it just returns the selected files paths in the string S_filename for later
	// use
	Variable refNum
	Open/F="PE Data Files (*.DX):.DX;"/D/R/MULT=1/M="UVvis Files ... " refNum
	
	//move the paths to a safer location
	String outputPaths = S_fileName

	//check to see If the user cancelled
	If (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	Else
		//If the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		Display/N=UVvis //Empty window to fill.
		
		Print "Loading JCAMP Data:"
		print OutPutPaths
		String filename, path
		Variable i	//Iterator
		for(i=0; i<numFilesSelected; i+=1)
			//Pull the path
			path = StringFromList(i, outputPaths, "\r")
			//Pull the filename from the path
			filename = StringFromList(0,StringFromList((ItemsInList(path,":")-1),path,":"),".")

			LoadWave/N/D/Q/O/G path//load the file
			//Make the waves accessible to the function
			WAVE xWave = wave0
			WAVE Data = wave1

			begin=xWave[0]
			step=xWave[1]-xWave[0]
			//Setting the WAVE scaling and assigning the units of the abscissa to nm
			SetScale/P x, begin, step,"nm",Data
			//Assigning the units of the ordinate to A (absorbance)
			SetScale d 0,0,"A", Data
			//Rename the wave
			Duplicate/O Data $(filename+"_spec")
			//Display it
			AppendToGraph $(filename+"_spec")
		EndFor
	EndIf
	ModifyGraph mirror=2,minor(bottom)=1
	Label bottom "Wavelength (\\U)"
	Print "\rFinished loading data!\r"
	Printf "%g spectra were loaded in this import\r", i
	
	//Clean up
	KillWaves xWave, data
	
	Return 0	//Successful execution
End

Function LoadCWData()
//A function that loads text files that have been converted from .SPE files
//in winspec (Princeton Instruments)
//***NOTE: this will load all the .txt files in a give folder***
	
	//Ask the user to choose the path to the folder where the data lives
	NewPath/M="Choose the folder with your data ..."/O path 
	
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	Variable i=0	//an iterator
	
	String fileList=IndexedFile(path,-1, ".txt")	//List of files in the folder
	Variable length = ItemsInList(fileList)			//How many files?
	String filename=""								//Holder for the filename
	
	Variable myMod = ceil(length/20)
	
	If (length==0)//check to see If folder is empy
		print "Break - - No files in folder"
		Return -1
	EndIf
	
	//Some print out
	printf "Loading CW spectra"
	PrintF ":  -/"
	
	//We want the load time =)
	Variable timerRefNum = startMSTimer
	
	String toPrint = ""		//A string to store the names of the files we're loading
	//Loop through all the files in the folder
	For(i=0;i<length;i+=1)
		If(!mod(i, myMod ))
			Printf "*"
		EndIf
		//Update the file name
		filename=StringFromList(i,fileList)
		
		LoadWave/N/D/Q/O/J/P=path filename//load the file
		//Make the waves accessible to the function
		WAVE wave0 = wave0
		WAVE wave1 = wave1
				
		If(mod(i,10)==0)
			toprint+= "\r"
		EndIf
		
		toPrint+= filename+", "

		//copy waves, and name them appropriately
		filename = StringFromList(0,filename,".")
		Duplicate/O wave1, $(filename+"_spec")
	EndFor
	
	Printf "/-\r"
	//Print out the names of the files loaded and some extra print out
	Print toPrint
	Printf "\r\rFinished loading data!\r\r"
	Printf "%g spectra were loaded in this import\r", i
	
	//Print out the time it took to load the files.
	Variable microSeconds = stopMSTimer(timerRefNum)
	
	Printf "Time to load: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	KillWaves/Z wave0,wave1//Removes temporary waves
	
	Return 0	//Successful execution
End

Function LoadFSRSData()
//A function that loads FSRS data as formated by the FSRS instrumentation software
//for a single (Raman pump) chopper instrument
//***NOTE: this will load ALL the files in a given folder***
//see https://www.github.com/david-hoffman/FSRS-LabVIEW
	
	Variable loadall		//For the prompt, 1 means only Spectra 2 means probe spectra too
	Prompt LoadAll,"What would you like to load?",popup,"Spectra;All"
	DoPrompt "Importing FSRS Data (1 chopper!)...",Loadall
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	//Ask the user to choose the path to the folder where the data lives
	NewPath/M="Choose the folder with your data ..."/O path 
	
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	Variable i=0	//an iterator
	
	String fileList=IndexedFile(path,-1, "????")	//List of files in the folder
	fileList = ListMatch(fileList,"!.DS_Store")		//Remove the .DS_Store File (made by macs)
	Variable length = ItemsInList(fileList)			//How many files?
	String filename=""								//Holder for the filename
	
	Variable myMod = ceil(length/20)
	
	If (length==0)//check to see If folder is empy
		print "Break - - No files in folder"
		Return -1
	EndIf
	
	//Some print out
	printf "Loading FSRS spectra"
	If(loadall==2)
		PrintF " and Pump/Probe"
	EndIf
	PrintF ":  -/"
	
	//We want the load time =)
	Variable timerRefNum = startMSTimer
	
	String toPrint = ""		//A string to store the names of the files we're loading
	//Loop through all the files in the folder
	For(i=0;i<length;i+=1)
		If(!mod(i, myMod ))
			Printf "*"
		EndIf
		//Update the file name
		filename=StringFromList(i,fileList)
		
		LoadWave/N/D/Q/O/J/P=path filename//load the file
		//Make the waves accessible to the function
		WAVE wave0 = wave0
		WAVE wave1 = wave1
		WAVE wave2 = wave2
		
		If(mod(i,10)==0)
			toprint+= "\r"
		EndIf
		
		toPrint+= filename+", "

		//copy waves, and name them appropriately
		Duplicate/O wave0, $(filename+"_spec")
		If(loadall==2)
			Duplicate/O wave1, $(filename+"_PumpOn")
			Duplicate/O wave2, $(filename+"_PumpOff")
		EndIf
	EndFor
	
	Printf "/-\r"
	//Print out the names of the files loaded and some extra print out
	Print toPrint
	Printf "\r\rFinished loading data!\r\r"
	Printf "%g spectra were loaded in this import\r", i
	
	//Print out the time it took to load the files.
	Variable microSeconds = stopMSTimer(timerRefNum)
	
	Printf "Time to load: %02d:%05.2f\r", microSeconds/60e6, mod(microSeconds,60e6)/1e6
	
	KillWaves/Z wave0,wave1,wave2//Removes temporary waves
	
	Return 0	//Successful execution
End

Function LoadFSRSData2()
//A function that loads FSRS data as formated by the FSRS instrumentation software
//for a two chopper instrument
//***NOTE: this will load ALL the files in a given folder***
//see https://www.github.com/david-hoffman/FSRS-LabVIEW
	Variable loadall
	
	Prompt LoadAll,"What would you like to load?",popup,"Spectra;TA;All"
	DoPrompt "Importing Red Table Data (2 chopper!)...",Loadall
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	NewPath/M="Choose the folder with your data ..."/O path 
	
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	Variable i = 0
	String filename=IndexedFile(path,i, "????")
	
	If (CmpStr(filename,".DS_Store")==0)//for Mac's, .DS_store is a hidden file in ever folder
		i+=1//this skips the .DS_store file
		filename=IndexedFile(path,i, "????")
	EndIf
	
	If (CmpStr(filename,"")==0)//check to see If folder is empy
		print "Break - - No files in folder"
		Return -1
	EndIf
	
	printf "Loading "
	String toPrint = ""
	Do
		LoadWave/N/D/Q/O/J/P=path filename//load the file
		WAVE wave0 = wave0
		WAVE wave1 = wave1
		WAVE wave2 = wave2
		WAVE wave3 = wave3
		WAVE wave4 = wave4
		WAVE wave5 = wave5
		WAVE wave6 = wave6
		WAVE wave7 = wave7
		If(mod(i,10)==0)
			toPrint+= "\r"
		EndIf
		toPrint+= filename+", "
		//copy waves
		If(loadall==1||loadall==3)
			Duplicate/O wave0, $(filename+"exc_spec")
			Duplicate/O wave1, $(filename+"gnd_spec")
		EndIf
		If(loadall==2||loadall==3)//choose 1 for this to get all the data
			Duplicate/O wave6, $(filename+"_TAROn")
			Duplicate/O wave7, $(filename+"_TAROff")
		EndIf
		If(loadall==3)
			Duplicate/O wave2, $(filename+"_PumpOnExc")
			Duplicate/O wave3, $(filename+"_PumpOffExc")
			Duplicate/O wave4, $(filename+"_PumpOnGnd")
			Duplicate/O wave5, $(filename+"_PumpOffGnd")
		EndIf
		i+=1	//increment
		filename=IndexedFile(path,i,"????")
	While (CmpStr(filename,"")!=0)
	Print toPrint
	Printf "\r\rFinished loading data!\r\r"
	Printf "%g spectra were loaded in this import\r", i
	
	//clean up	
	KillWaves/Z wave0,wave1,wave2,wave3,wave4,wave5,wave6,wave7//Removes temporary waves
	
	Return 0	//Successful load
End

Function LoadDAQData()
//A function that loads data aquired by the Lockin amplifer
//see https://www.github.com/david-hoffman/FSRS-LabVIEW

	//The following open command does not actually open the selected file
	//it just returns the selected files paths in the string S_filename for later
	//use
	Variable refNum
	Open/F="????"/D/R/MULT=1/M="Choose DAQ Files ... " refNum
	
	//move the paths so a safer location
	String outputPaths = S_fileName
	
	//check to see If the user cancelled
	If (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	Else
		//If the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		Variable i, j // declare some dummy variables
		
		String path, filename
		
		Display/N=DAQ
		for(i=0; i<numFilesSelected; i+=1)
			//pull the path
			path = StringFromList(i, outputPaths, "\r")
			//pull the filename
			filename = stringfromlist((ItemsInList(path,":")-1),path,":")
			
			LoadWave/N/D/Q/O/J path//load the file
			
			WAVE wave0 = wave0	//Delays
			WAVE wave1 = wave1	//Data
			WAVE wave2 = wave2	//St dev
			
			Print "Loading " + S_path+ " ... "
			If(V_flag>3||V_flag<2)//is the file a DAQ scan file?
				Print "File is not properly formatted, will cancel load of "+filename
			Else //It is a DAQ scan file
				Print "Load successful!"
				Duplicate/O wave0, $(filename+"D")
				Duplicate/O wave1, $(filename+"I")
				AppendToGraph $(filename+"I") vs $(filename+"D")
				ModifyGraph mode=3,marker=8
				If(V_flag==3) //It's a newer version of DAQ scan
					Duplicate/O wave2, $(filename+"S")
					ErrorBars $(filename+"I") Y,WAVE=($(filename+"S"),$(filename+"S"))
				EndIf
			EndIf
			
			//Make the formating of the graph pretty
			If(WaveMax(wave1)>2000)
				ModifyGraph prescaleExp(left)=-3
				Label left "Intensity (mV)"
			Else
				Label left "Intensity (µV)"
			EndIf
			
			If(WaveMax(wave0)>10000)
				ModifyGraph prescaleExp(bottom)=-3
				Label bottom "Delay (ps)"
			Else
				Label bottom "Delay (fs)"
			EndIf
			
			ModifyGraph mirror=2,minor(bottom)=1
			
			for(j=0;j<V_flag;j+=1)
				// Kill any left over waves
				Killwaves/Z $("WAVE"+num2str(j))
			EndFor
		EndFor
	EndIf
End

Function LoadDetXCData()
//A function that loads FSRS data as formated by the FSRS instrumentation software
//see https://www.github.com/david-hoffman/FSRS-LabVIEW
	String pathName=""
	String fileName="",TopWinName="",path=""
	Variable refNum,algoType
	
	//Build the prompt to ask the user how they would like to process this data (to be explained later)
	Prompt algoType,"How would you like to process this data?",popup,"Gauss;Kerr"
	//Display the prompt to the user
	DoPrompt "Detector Cross Correlation",algoType
	
	If( V_Flag )
		Print "Cancelled"
		return 0	// user canceled
	EndIf
	
	//The following open command does not actually open the selected file it just
	//returns the selected files paths in the string S_filename for later use
	Open/F="????"/D/R/MULT=1/M="Choose Det XC files" refNum
	
	//move the paths so a safer location
	String outputPaths = S_fileName
	
	//check to see If the user cancelled
	If (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	Else
		//If the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		Variable i,j // declare some iterators
		
		Variable V_avg, V_sdev,V_Min,V_max
		
		PauseUpdate //Things are about to get hairy
		
		for(i=0; i<numFilesSelected; i+=1)
			//Pull out the path from the list of paths
			path = StringFromList(i, outputPaths, "\r")
			//Pull out the filename from the path
			filename = stringfromlist((ItemsInList(path,":")-1),path,":")
			//Load the data, treating it as a matrix with the scaling WAVE on the side
			LoadWave/N/Q/O/J/M/U={0,1,0,0}/D/K=0/L={0,0,0,1,0} path
			
			Print "Loading " + S_path+ " ... "
			
			WAVE wave0 = wave0
			
			If(DimSize(wave0,1)!=1340)
				//Not the right size!
				Print "Load Failed, file not properly formatted"
				Return 1 //Exit Macro
			Else
				//Place the data in a properly named WAVE
				Duplicate/O wave0, $filename
				//Rotate the data so that when displayed the x axis is the pixel and the y is delay
				MatrixTranspose $filename
				
				If(algoType==1)
					XCMaxG($filename)//Call the subfunction that will do all the calculations, using a single
					//gaussian as the model for the instrument response
				Else
					XCMaxK($filename)//Call the subfunction that will do all the calculations, using a
					//gaussian and convoluted exponential for the instrument response
				EndIf
				
				//Has this file been loaded before?
				If(WinType(filename+"0")==0)
					//The window does NOT already exist, because we're doing this programmically,
					//the name will automatically have a 0 appended set up the host window
					Display/W=(0,0,350,400)/N=$filename
					
					//Just incase it doesn't simply have a 0 appended pull the name
					TopWinName = S_name
					
					DoWindow/T $TopWinName, "Detector Cross Correlation: "+filename
					//Display/HOST=$TopWinName/N=DetXCImage/W=(0,1,1,0.5)/FG=(FL,FT,FR,*)// The Image
					WAVE shiftxForImages = shiftxForImages
					WAVE shiftx = shiftx
					If(WaveExists(shiftxForImages)&&WaveExists(shiftx))
						AppendImage/L=detxc $filename vs {shiftxForImages,*}//Make an Image, not a contour
						AppendtoGraph/L=detxc $(filename+"_Max") vs shiftx
						ModifyImage $filename ctab= {*,*,Geo32,0}//Make the plot colorful
						Label bottom "Raman Shift (cm\\S-1\\M)"
					Else
						AppendImage/L=detxc $filename//Make an Image, not a contour
						AppendtoGraph/L=detxc $(filename+"_Max")
						ModifyImage $filename ctab= {*,*,Geo32,0}//Make the plot colorful
					EndIf
					Label detxc "Time Delay (fs)"
		
					ModifyGraph lsize($(filename+"_Max"))=2,rgb($(filename+"_Max"))=(0,0,0)
					
					If(WaveExists(shiftx))
						AppendToGraph/L=Center $(filename+"_Max") vs shiftx
						Appendtograph/R=width $(filename+"_Width") vs shiftx
						Label bottom "Raman Shift (cm\\S-1\\M)"
					Else
						AppendToGraph/L=Center $(filename+"_Max")
						Appendtograph/R=width $(filename+"_Width")
					EndIf
					
					ModifyGraph muloffset($(filename+"_width"))={0,2*sqrt(ln(2))}
					ModifyGraph axRGB(center)=(65535,0,0),tlblRGB(center)=(65535,0,0), alblRGB(center)=(65535,0,0)
					ModifyGraph axRGB(width)=(0,0,65535),tlblRGB(width)=(0,0,65535), alblRGB(width)=(0,0,65535)
					ModifyGraph rgb($(filename+"_Width"))=(0,0,65535)
					ModifyGraph mirror(detxc)=2,lblPos(Center)=50,lblPos(width)=50,axisEnab(detxc)={0.52,1}
					ModifyGraph axisEnab(Center)={0,0.48},axisEnab(width)={0,0.48},freePos(detxc)=0
					ModifyGraph freePos(Center)=0,freePos(width)=0, lblPos(detxc)=50
					ModifyGraph standoff(bottom)=0
					
					ModifyGraph margin(left)=36,margin(bottom)=36,margin(top)=12,margin(right)=36
					ModifyGraph height=312,width=168,expand=1.5
					
					Label Center "Time Delay (fs)"
					Label width "FWHM (fs)"
					ModifyGraph mirror(bottom)=2
					Legend/A=RC/x=5/Y=-10/C/N=Leg0/J "\\s("+filename+"_Max#1) t\\B0\\M\r\\s("+filename+"_Width) FWHM"
					
					SetDrawLayer UserFront
					DrawLine 0,0.52,1,0.52
					DrawLine 0,0.48,1,0.48
				Else
					TopWinName = filename+"0"
				EndIf
				
				//Processing the data and displaying values
				WaveStats/Q $(filename+"_Max")
				
				String chirp="\rÆ = "+ num2str(round(V_max-V_min))
				
				WaveStats/Q $(filename+"_width")
				TextBox/A=LC/X=5/Y=-10/W=$TopWinName/C/N=text2/X=64/Y=0 "Avg FWHM = "+num2str(round(2*sqrt(ln(2))*V_avg))+chirp
				
				//Now we integrate the signal over the whole window
				IntegrateMatrix($(filename),("Int_"+filename))
				//And process it
				WAVE w_coef, w_sigma
				CurveFit/N=1/W=2/Q gauss  $("Int_"+filename) /D
				
				Variable Center = w_coef[2]
				Variable sigma = w_coef[3]
				Variable s_Center = w_sigma[2]
				Variable s_sigma = w_sigma[3]
				
				If(algoType==2)//We need to do some extra processing
					Make/FREE/D/O/N=6 kerr_coefs={w_coef[3],w_coef[2],w_coef[0],w_coef[1]*0.1,w_coef[3]*2,w_coef[1]}
					FuncFit/Q/N=1/W=2 kerr kerr_coefs $("Int_"+filename) /D
					Center = kerr_coefs[1]
					sigma = kerr_coefs[0]
					s_Center = w_sigma[1]
					s_sigma = w_sigma[0]
				EndIf
				
				//Displaying the integrated XC
				If(WinType(filename+"1")==0)//The window does NOT already exist, because we're
				//Doing this programmically, the name will automatically have a 0 appended
					Display/N=$filename $("Int_"+filename) as ("Integrated XC: "+filename)
					TopWinName = S_name
					AppendToGraph/W=$TopWinName $("fit_Int_"+filename)
					ModifyGraph/W=$TopWinName mirror=2,minor(bottom)=1
					ModifyGraph/W=$TopWinName mode($("Int_"+filename))=3,marker($("Int_"+filename))=8
					ModifyGraph/W=$TopWinName rgb($("fit_Int_"+filename))=(0,0,65535)
					ModifyGraph/W=$TopWinName nticks(left)=0,lblMargin(left)=20
					Label/W=$TopWinName left "Integrated XC (a.u.)"
					Label/W=$TopWinName bottom "Time (fs)"
				Else
					TopWinName = filename+"1"
				EndIf
				Legend/W=$TopWinName/C/N=Leg0/J "\\JC"+Date()+"\r\\s(Int_"+filename+") XC"
				AppendText/W=$TopWinName/N=Leg0 "\\s(fit_Int_"+filename+") Fit"
				AppendText/W=$TopWinName/N=Leg0 "t\B0\M = "+num2str(round(Center))+" ± "+num2str(round(s_Center))+" fs"
				AppendText/W=$TopWinName/N=Leg0 "\F'Symbol'\f01s\]0 = "+num2str(round(sigma))+" ± "+num2str(round(s_sigma))+" fs"
				AppendText/W=$TopWinName/N=Leg0 "FWHM = "+num2str(round(2*sqrt(ln(2))*sigma))+" ± "+num2str(round(2*sqrt(ln(2))*s_sigma))+" fs"
				
			EndIf
			KillWaves/Z wave0
		EndFor
		
		ResumeUpdate
		
	EndIf
End

Function LoadISRSData()
//A function that loads ISRS data as formated by the ISRS instrumentation software
//for a one chopper instrument
	String pathName=""
	String fileName=""
	Variable refNum
	
	//The following open command does not actually open the selected file it just
	//returns the selected files paths in the string S_filename for later use
	Open/F="????"/D/R/MULT=1/M="Choose Det XC files" refNum
	
	//move the paths so a safer location
	String outputPaths = S_fileName
	
	//check to see if the user cancelled
	If (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	Else
		//if the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		Variable i,j // declare some dummy variables
		
		//Set up a window to later append traces to
		Display/N=IntegratedData as ("Integrated ISRS")
		
		//Get the actual window name (just in case)
		String win=WinName(0,1)
		
		For(i=0; i<numFilesSelected; i+=1)
			String path = StringFromList(i, outputPaths, "\r")
			filename = stringfromlist((ItemsInList(path,":")-1),path,":")
			LoadWave/N/Q/O/J/M/U={0,1,0,0}/D/K=0/L={0,0,0,1,0} path//load the file
			WAVE wave0=wave0
			Print "Loading " + S_path+ " ... "
			If(DimSize(wave0,1)!=1340)
				Print "Load Failed, file not properly formatted"
				Return 1 //Exit Macro
			Else
				//Give the data a better name
				Duplicate/O wave0, $(filename)
				//Transpose the data so that when displayed the x axis is pixel and the y is delay
				MatrixTranspose $(filename)
				
				//append the integrated ISRS signal
				AppendToGraph/W=$win $IntegrateMatrix($(filename),("Int_"+filename))
				DoUpdate
				
			EndIf
			KillWaves/Z wave0
		EndFor
		
		FancifyTA()
		ColorWaves()
		
	EndIf
	
	Return 0
End

Function LoadDAQISRSData()
	//The difference between this and the regular LoadDAQData()
	//is that this function assumes the data to be an equallly spaced
	//waveform and scales the data appropriately
	
	//Also display all data in one graph
	String pathName=""
	String fileName=""
	Variable refNum
	
	//The following open command does not actually open the selected file it just
	//returns the selected files paths in the string S_filename for later use
	Open/F="????"/D/R/MULT=1/M="Choose DAQ Files ... " refNum
	
	//move the paths so a safer location
	String outputPaths = S_fileName
	
	//check to see if the user cancelled
	if (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	else
		//if the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		Variable i,j, begin, step// declare some dummy variables
		
		Display/N=DAQ //set up the window in which to display the data
		
		for(i=0; i<numFilesSelected; i+=1)
			String path = StringFromList(i, outputPaths, "\r")
			filename = stringfromlist((ItemsInList(path,":")-1),path,":")
			LoadWave/N/D/Q/O/J path//load the file
			Wave wave0 = wave0
			Wave wave1 = wave1
			Wave wave2 = wave2
			Print "Loading " + S_path+ " ... "
			if(V_flag!=3)//is the file a DAQ scan file?
				Print "File is not properly formatted, will cancel load of "+filename
			Else //It is a DAQ scan file
				Print "Load successful!"
				
				WAVE wave0
				WAVE wave1
				
				begin = wave0[0]
				step = wave0[1]-wave0[0]
				
				SetScale/P x, begin, step,"fs",wave1 //Setting the wave scaling and assigning the units of the abscissa to nm
				SetScale d 0,0,"µV", wave1 //Assigning the units of the ordinate to A (absorbance)
				
				Duplicate/O wave1, $(filename+"I")
				AppendToGraph $(filename+"I")
				if(WaveMax(wave1)>2000)
					ModifyGraph prescaleExp(left)=-3
					Label left "Intensity (\U)"
				Else
					Label left "Intensity (\U)"
				EndIf
				If(WaveMax(wave0)>10000)
					ModifyGraph prescaleExp(bottom)=-3
					Label bottom "Delay (\U)"
				Else
					Label bottom "Delay (\U)"
				EndIf
				Duplicate/O wave2, $(filename+"S")
				ErrorBars $(filename+"I") Y,wave=($(filename+"S"),$(filename+"S"))
			EndIf
			ModifyGraph mode=3,marker=8
			ModifyGraph mirror=2,minor(bottom)=1
			for(j=0;j<V_flag;j+=1)
				// Kill any left over waves
				Killwaves/Z $("wave"+num2str(j))
			EndFor
		endfor
	EndIf
End