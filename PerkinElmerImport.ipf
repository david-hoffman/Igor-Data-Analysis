#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function loadPEData()
//A function that loads Perkin Elmer's propietary data format
//Initially designed to read data from their Lambda 2 and Lambda 6
//UV/vis spectrometers.
	String toReturn = ""
	String pathName=""
	Variable refNum

	//Header info to be stored in these variables
	String fileName=""
	Variable specType, dataType, numPoints
	Variable xMax,xMin, dataInterval, specManipFlag
	Variable ordType, scaleFactor, ordMax, ordMin, smoothFact
	Variable numAccum, expansionFactor, accessory
	Variable ordModeFormat, cycleNum, numCycles, cycleTime, scanSpeed, response,NIRSensitivity, slit
	Variable uvGain, detSwitch, cellNum, mySmooth, Lamp, accessory2, timeDriveWaveLength
	Variable instrumentID, methodName, dInterval, dOrder
	String ordUnits=PadString("",4,0), abUnits=PadString("",4,0)
	String myDate=PadString("",8,0), myTime=PadString("",8,0), idLine=PadString("",72,0)
	Variable globalHeader, localHeader, numBytes, scaleDef, viewOrdMax, viewOrdMin
	String region=PadString("",1,0)
	Variable recordStatus
	String autoSave=PadString("",1,0),autoPrint=PadString("",1,0)
	Variable numCells, concFactor
	String autoZero=PadString("",1,0),azeroRestore=PadString("",1,0)
	Variable sipperTime
	String sipperAutoPurge=PadString("",1,0)
	Variable sipper, analysisTimeCPRG ,sipperDelay, Temperature
	String aplot=PadString("",1,0), autoSampler=PadString("",1,0),trayDef=PadString("",9,0),flush=PadString("",1,0)
	Variable startloc, endloc, azeroLoc, concStartLoc,concEndLoc, concAzeroLoc, numReplicates

	//The following open command does not actually open the selected file
	// it just returns the selected files paths in the string S_filename for later
	// use
	Open/F="PE Data Files (*.SP):.SP;"/D/R/MULT=1/M="UVvis Files ... " refNum
	
	//move the paths to a safer location
	String outputPaths = S_fileName
	
	//check to see If the user cancelled
	If (strlen(outputPaths) == 0)
		Print "Cancelled"
		Return 0
	Else
		//If the user did select files calculate how many they chose
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")//The number of files the user selected
		
		Variable i,j // declare some dummy variables
		
		Display/N=UVvis //Empty window to fill.
		
		Print "Loading Perkin-Elmer PECSS Data:"
		
		for(i=0; i<numFilesSelected; i+=1)
			String path = StringFromList(i, outputPaths, "\r")
			Print path
			filename = stringfromlist((ItemsInList(path,":")-1),path,":")
			//Time to read in the important stuff!

			Open/R refNum as path
			
			FSetPos refNum, 8
			FBinRead/B=3/F=1/U refNum, specType
			
			FSetPos refNum, 9
			FBinRead/B=3/F=1/U refNum, dataType
			
			FSetPos refNum, 10
			FBinRead/B=3/F=3/U refNum, numPoints
			
			FSetPos refNum, 14
			FBinRead/B=3/F=3/U refNum, xMax
			xMax/=100
			
			FSetPos refNum, 18
			FBinRead/B=3/F=3/U refNum, xMin
			xMin/=100
			
			FSetPos refNum, 22
			FBinRead/B=3/F=2/U refNum, dataInterval
			dataInterval/=100
			
			FSetPos refNum, 24
			FBinRead/B=3/F=1/U refNum, specManipFlag
			
			FSetPos refNum, 25
			FBinRead/B=3/F=1/U refNum, ordType
			
			FSetPos refNum, 28
			FBinRead/B=3/F=3/U refNum, scaleFactor
			
			FSetPos refNum, 32
			FBinRead/B=3/F=3/U refNum, ordMin
			
			FSetPos refNum, 36
			FBinRead/B=3/F=3/U refNum, ordMax
			
			FSetPos refNum, 40
			FBinRead/B=3/F=2/U refNum, numAccum
			
			FSetPos refNum, 42
			FBinRead/B=3/F=1/U refNum, expansionFactor
			
			FSetPos refNum, 46
			FBinRead/B=3/F=1/U refNum, smoothFact
			
			FSetPos refNum, 48
			FBinRead/B=3/F=1/U refNum, accessory
			
			FSetPos refNum, 49
			FBinRead/B=3/F=1/U refNum, ordModeFormat
			
			FSetPos refNum, 50
			FBinRead/B=3/F=2/U refNum, cycleNum
			
			FSetPos refNum, 52
			FBinRead/B=3/F=2/U refNum, numCycles
			
			FSetPos refNum, 54
			FBinRead/B=3/F=3/U refNum, cycleTime
			cycleTime/=100
			
			FSetPos refNum, 58
			FBinRead/B=3/F=2/U refNum, scanSpeed
			scanSpeed/=10
			
			FSetPos refNum, 60
			FBinRead/B=3/F=1/U refNum, response
			response/=10
			
			FSetPos refNum, 61
			FBinRead/B=3/F=1/U refNum, NIRSensitivity
			
			FSetPos refNum, 62
			FBinRead/B=3/F=2/U refNum, slit
			slit/=10
			
			FSetPos refNum, 64
			FBinRead/B=3/F=2/U refNum, uvGain
			uvGain/=10
			
			FSetPos refNum, 66
			FBinRead/B=3/F=2/U refNum, detSwitch
			detSwitch/=10
			
			FSetPos refNum, 68
			FBinRead/B=3/F=1/U refNum, cellNum
			
			FSetPos refNum, 69
			FBinRead/B=3/F=2/U refNum, mySmooth
			
			FSetPos refNum, 71
			FBinRead/B=3/F=1/U refNum, Lamp
			
			FSetPos refNum, 72
			FBinRead/B=3/F=2/U refNum, accessory2
			
			FSetPos refNum, 76
			FBinRead/B=3/F=3/U refNum, timeDriveWaveLength
			timeDriveWaveLength/=100
			
			FSetPos refNum, 80
			FBinRead/B=3/F=2/U refNum, instrumentID
			
			FSetPos refNum, 82
			FBinRead/B=3 refNum, methodName
			
			FSetPos refNum, 136
			FBinRead/B=3/F=2/U refNum, dInterval
			
			FSetPos refNum, 138
			FBinRead/B=3/F=1/U refNum, dOrder
			
			FSetPos refNum, 146
			FBinRead/B=3 refNum, ordUnits
			
			FSetPos refNum, 150
			FBinRead/B=3 refNum, abUnits
			
			FSetPos refNum, 156
			FBinRead/B=3 refNum, myDate
			
			FSetPos refNum, 164
			FBinRead/B=3 refNum, myTime
			
			FSetPos refNum, 176
			FBinRead/B=3 refNum, idLine
			
			FSetPos refNum, 248
			FBinRead/B=3/F=1/U refNum, globalHeader
			
			FSetPos refNum, 249
			FBinRead/B=3/F=1/U refNum, localHeader
			
			FSetPos refNum, 250
			FBinRead/B=3/F=2/U refNum, numBytes
			
			FSetPos refNum, 252
			FBinRead/B=3/F=3/U refNum, scaleDef
			
			FSetPos refNum, 311
			FBinRead/B=3/F=3/U refNum, viewOrdMin
			
			FSetPos refNum, 315
			FBinRead/B=3/F=3/U refNum, viewOrdMax
			
			FSetPos refNum, 319
			FBinRead/B=3 refNum, region
			
			FSetPos refNum, 320
			FBinRead/B=3/F=1/U refNum, recordStatus
			
			FSetPos refNum, 321
			FBinRead/B=3 refNum, autoSave
			
			FSetPos refNum, 322
			FBinRead/B=3 refNum, autoPrint
			
			FSetPos refNum, 323
			FBinRead/B=3/F=1/U refNum, numCells
			
			FSetPos refNum, 324
			FBinRead/B=3/F=3/U refNum, concFactor
			
			FSetPos refNum, 328
			FBinRead/B=3 refNum, autoZero
			
			FSetPos refNum, 331
			FBinRead/B=3 refNum, azeroRestore
			
			FSetPos refNum, 332
			FBinRead/B=3/F=2/U refNum, sipperTime
			sipperTime/=10
			
			FSetPos refNum, 334
			FBinRead/B=3 refNum, sipperAutoPurge
			
			FSetPos refNum, 335
			FBinRead/B=3/F=1/U refNum, sipper
			
			FSetPos refNum, 336
			FBinRead/B=3/F=3/U refNum, analysisTimeCPRG
			analysisTimeCPRG/=100
			
			FSetPos refNum, 340
			FBinRead/B=3/F=2/U refNum, sipperDelay
			sipperDelay/=10
			
			FSetPos refNum, 342
			FBinRead/B=3/F=2/U refNum, Temperature
			Temperature/=10
			
			FSetPos refNum, 344
			FBinRead/B=3 refNum, aplot
			
			FSetPos refNum, 345
			FBinRead/B=3 refNum, autoSampler
			
			FSetPos refNum, 346
			FBinRead/B=3 refNum, trayDef
			
			FSetPos refNum, 355
			FBinRead/B=3 refNum, flush
			
			FSetPos refNum, 356
			FBinRead/B=3/F=2/U refNum, startloc
			
			FSetPos refNum, 358
			FBinRead/B=3/F=2/U refNum, endloc
			
			FSetPos refNum, 360
			FBinRead/B=3/F=2/U refNum, azeroLoc
			
			FSetPos refNum, 362
			FBinRead/B=3/F=2/U refNum, concStartLoc
			
			
			FSetPos refNum, 364
			FBinRead/B=3/F=2/U refNum, concEndLoc
			
			FSetPos refNum, 366
			FBinRead/B=3/F=2/U refNum, concAzeroLoc
			
			FSetPos refNum, 368
			FBinRead/B=3/F=2/U refNum, numReplicates
			
			//Load the WAVE from the binary file
			//GBLoadWave/Q/O/B/N=WAVE/T={32+64,4}/S=512/W=1/Y={0,ScaleDef/100/scaleFactor} path
			GBLoadWave/Q/O/B/N=WAVE/T={32,4}/S=512/W=1/Y={0,ScaleDef/100/scaleFactor} path
			WAVE wave0 = wave0
			//Make the wavenote
			Note wave0, "Spectrum path: "+path
			Note wave0, "Spectrum type:  "
			Switch(specType)
				case 0:
					Note/NOCR wave0, "IR"
					SetScale/I x, xMin, xMax, "cm\S-1\M", wave0
					break
				case 1:
					Note/NOCR wave0, "UV"
					SetScale/I x, xMin, xMax, "nm", wave0
					break
				case 2:
					Note/NOCR wave0, "LS"
					break
				case 11:
					Note/NOCR wave0, "TDRIVE"
					SetScale/I x, xMax, xMin, "s", wave0
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(specType)+")"
					break
			endswitch
			
			Note wave0, "Number of points: "+num2str(numPoints)
			Note wave0, "Abscissa maximum: " + num2str(xMax)
			If(specType == 11)
				Note/NOCR wave0, " s"
			Else
				Note/NOCR wave0, " nm"
			EndIf
			Note wave0, "Abscissa minimum: " + num2str(xMin)
			If(specType == 11)
				Note/NOCR wave0, " s"
			Else
				Note/NOCR wave0, " nm"
			EndIf
			Note wave0, "Data interval: "+num2str(dataInterval)
			If(specType == 11)
				Note/NOCR wave0, " s"
			Else
				Note/NOCR wave0, " nm"
			EndIf
			Note wave0, "Type of data manipulation: "
			Switch(specManipFlag)
				case 1:
					Note/NOCR wave0, "log"
					break
				case 2:
					Note/NOCR wave0, "diff"
					break
				case 3:
					Note/NOCR wave0, "arith"
					break
				case 16:
					Note/NOCR wave0, "Arith"
					break
				case 32:
					Note/NOCR wave0, "modified"
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(specManipFlag)+")"
					break
			endswitch
			Note wave0, "Ordinate type: "
			switch(ordType)
				case 1:
					Note/NOCR wave0, "A"
					SetScale d 0,0,"A", wave0
					break
				case 2:
					Note/NOCR wave0, "%T"
					SetScale d 0,0,"%T", wave0
					break
				case 3:
					Note/NOCR wave0, "I"
					SetScale d 0,0,"I", wave0
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(ordType)+")"
					break
			endswitch
			Note wave0, "Internal data scale factor: "+num2istr(scaleFactor)
			Note wave0, "Ordinate minimum: "+num2str(ordMin*ScaleDef/100/scaleFactor)
			Note wave0, "Ordinate maximum: "+num2str(ordMax*ScaleDef/100/scaleFactor)
			Note wave0, "Number of accumulations: "+num2str(numAccum)
			Note wave0, "Absorbance expansion factor: "+num2str(expansionFactor)
			Note wave0, "Smooth factor: " + num2str(smoothFact)
			Note wave0, "Accesory: "
			switch(accessory)
				case 1:
					Note/NOCR wave0, "manual"
					break
				case 2:
					Note/NOCR wave0, "Sipper"
					break
				case 3:
					Note/NOCR wave0, "Cell Changer"
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(accessory)+")"
					break
			endswitch
			Note wave0, "Ordinate: "
			switch(ordModeFormat)
				case 0:
					Note/NOCR wave0, "Transmission"
					break
				case 1:
					Note/NOCR wave0, "1st Derivative"
					break
				case 2:
					Note/NOCR wave0, "2nd Derivative"
					break
				case 3:
					Note/NOCR wave0, "3rd Derivative"
					break
				case 4:
					Note/NOCR wave0, "4th Derivative"
					break
				case 5:
					Note/NOCR wave0, "Absorbance"
					break
				case 6:
					Note/NOCR wave0, "Log epsilon"
					break
				case 7:
					Note/NOCR wave0, "FR"
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(ordModeFormat)+")"
			endswitch
			Note wave0, "Cycle number: "+num2str(cycleNum)
			Note wave0, "Number of cycles: "+ num2str(numCycles)
			Note wave0,  "Cycle time: "+num2str(cycleTime)+" s"
			Note wave0, "Scan speed: " + num2str(scanSpeed)
			Note wave0, "Response time: " +num2str(response)
			Note wave0, "NIR sensitivity: "+num2str(NIRSensitivity)
			Note wave0, "Slit width: "+num2str(slit)+" nm"
			Note wave0, "UV gain factor: " + num2str(uvGain*10)
			Note wave0, "Detector change wavelength: " + num2str(detSwitch)
			Note wave0, "Cell number: "+num2str(cellNum)
			Note wave0, "Smooth: "+num2str( mySmooth)
			Note wave0, "Lamp: "+num2str(Lamp)
			Note wave0, "Accessory: "
			switch(accessory2)
				case 1:
					Note/NOCR wave0, "5-cell Changer"
					break
				case 2:
					Note/NOCR wave0, "6-cell Changer"
					break
				case 3:
					Note/NOCR wave0, "13-cell Changer"
					break
				case 4:
					Note/NOCR wave0, "8-cell Changer"
					break
				case 8:
					Note/NOCR wave0, "Single Sipper"
					break
				case 9:
					Note/NOCR wave0, "Single Sipper + Multisampler"
					break
				case 10:
					Note/NOCR wave0, "Super Sipper"
					break
				case 11:
					Note/NOCR wave0, "Super Sipper + Multisampler"
					break
				case 12:
					Note/NOCR wave0, "AS-90"
					break
				case 13:
					Note/NOCR wave0, "Single Sipper + AS-90"
					break
				case 14:
					Note/NOCR wave0, "Super Sipper + AS-90"
					break
				case 31:
					Note/NOCR wave0, "Multisampler"
					break
				case 32:
					Note/NOCR wave0, "Temperature sensor"
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(accessory2)+")"
					break
			endswitch
			Note wave0, "Time drive wavelength: "+num2str(timeDriveWaveLength)+" nm"
			Note wave0, "Instrument Identification: "
			switch(instrumentID)
				case 505:
					Note/NOCR wave0, "Lambda 15,17"
					break
				case 509:
					Note/NOCR wave0, "Lambda 19 NIR"
					break
				case 510:
					Note/NOCR wave0, "Lambda 19 UV/VIS"
					break
				case 502:
					Note/NOCR wave0, "Lambda 2"
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(instrumentID)+")"
					break
			endswitch
			Note wave0, "STORE/RESTORE method name: "+num2str(methodName)
			Note wave0, "Derivative interval: "+num2str(dInterval)
			Note wave0, "Derivative order: "+num2str(dOrder)
			Note wave0, "Ordinate units: "+ordUnits
			Note wave0, "Abscissa unit: "+abUnits
			Note wave0, "Date: "+myDate
			Note wave0, "Time: " + myTime
			Note wave0, "Info: "+idLine
			Note wave0, "Global Header Format: "+num2str(globalHeader)
			Note wave0, "Local Header Format: "+ num2str(localHeader)
			Note wave0, "Number of bytes in extension: "+num2str(numBytes)
			Note wave0, "Scale definition * 100: "+num2str(scaleDef)
			Note wave0, "Ordinate minimum used in view: "+num2str(viewOrdMin)
			Note wave0, "Ordinate maximum used in view: "+num2str(viewOrdMax)
			Note wave0, "Memory region: "+region
			Note wave0, "Recorder Status: "
			switch(recordStatus)
				case 0:
					Note/NOCR wave0, "off"
					break
				case 1:
					Note/NOCR wave0, "on"
					break
				case 2:
					Note/NOCR wave0, "cont"
					break
				default:
					Note/NOCR wave0, "UNKNOWN ("+num2str(recordStatus)+")"
					break
			endswitch
			Note wave0, "Autosave: "+autoSave
			Note wave0, "Autoprint: "+autoPrint
			Note wave0, "Number of cells for CPRG: "+num2str(numCells)
			Note wave0, "Concentration factor for CONC: "+num2str(concFactor)
			Note wave0, "Setting for autozero: "+autoZero
			Note wave0, "Azero for RESTORE: "+azeroRestore
			Note wave0, "Sipper sampling time: "+num2str(sippertime)+" s"
			Note wave0, "Sipper Auto Purge: "+sipperAutoPurge
			Note wave0, "Sipper: "
			If(sipper)
				Note/NOCR wave0, "On"
			ElseIf(!sipper)
				Note/NOCR wave0, "Off"
			Else
				Note/NOCR wave0, "UKNOWN ("+num2str(sipper)+")"
			EndIf
			Note wave0, "Total analysis time for CPRG: "+num2str(analysisTimeCPRG)+" min"
			Note wave0, "Sipper delay time: "+num2str(sipperDelay)+" s"
			Note wave0, "Temperature: "+num2str(Temperature)
			Note wave0, "Autoplot: "+aplot
			Note wave0, "Autosampler active: "+autoSampler
			Note wave0, "Tray definition filenames: "+trayDef
			Note wave0, "Flush: "+flush
			Note wave0, "Start location: "+num2str(startloc)
			Note wave0, "End location: "+num2str(endloc)
			Note wave0, "Autozero location: "+num2str(azeroLoc)
			Note wave0, "CONC start location: "+num2str(concStartLoc)
			Note wave0, "CONC end location: "+num2str(concEndLoc)
			Note wave0, "CONC autozero location: "+num2str(concAzeroLoc)
			Note wave0, "Number of replicates: "+num2str(numReplicates)
			
			//Print note(wave0)
			
			Duplicate/O wave0 $(StringFromList(0,filename,".")+"_spec")
			
			AppendToGraph $(StringFromList(0,filename,".")+"_spec")
			
			//Print "/--------------------------------------/"
		EndFor
	EndIf
	ModifyGraph mirror=2,minor(bottom)=1
	Label left "Absorbance (\\U)"
	If(SpecType==11)
		Label bottom "Time (\\U)"
	Else
		Label bottom "Wavelength (\\U)"
	EndIf

	KillWaves/Z wave0
	Return 0
End