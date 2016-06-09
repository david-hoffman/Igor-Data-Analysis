#pragma rtGlobals=3		// Use modern global access method.
#include ":Utilities"

//At some point we should add foldering to keep the data in seperate places
//Make it so that all intermediate waves are kept in a seperate folder will the baseline subtracted data
//(Only the final one) is saved in the root folder.

//Interesting idea: use the fact that for the time series the baselines from one time point to the next will 
//be similar. Maybe I can use this to my advantage as a kind of adaptive fitting process...
//

Function BetterSDev(myWave,cutoff)
//A method for calculating a better estimate for the standard deviation by first removing
//outliers and then calculating the standard deviation.
	WAVE myWave
	Variable cutoff
	Variable v_npnts,v_avg,v_sdev
	WaveStats/Z/Q myWave
	Duplicate/O myWave myTempWave
	variable i=0
	For(i=0;i<v_npnts;i+=1)
		//Remove outliers
		If(abs(myTempWave[i]-v_avg)>(cutoff*V_sdev))
			myTempWave[i]=NaN
		EndIf
	EndFor
	//Recalculate wavestats
	WaveStats/Z/Q myTempWave
	KillWaves/Z myTempWave
	Return v_sdev
End

Function BetterSDev2(myWave)
//A method for calculating a better estimate for the standard deviation by iteratively removing
//outliers and then calculating the standard deviation.
	WAVE myWave //Wave to calculate the standard deviation
	
	Variable v_npnts,v_avg,v_sdev,old_npnts//Some variables
	
	Duplicate/FREE/O myWave myTempWave//make a temp wave to do the calculations on
	
	variable i=0
	
	WaveStats/Z/Q myTempWave
	
	Do
		old_npnts=v_npnts
	
		For(i=0;i<v_npnts;i+=1)
			//Remove outliers
			If(abs(myTempWave[i]-v_avg)>(2*V_sdev))
				myTempWave[i]=NaN
			EndIf
		EndFor
		//Recalculate wavestats
		WaveStats/Z/Q myTempWave
	While(old_npnts!=v_npnts)//Check to see if any outliers were removed, continue if yes

	Return v_sdev

End

Function ForcePnts(myWave, myPnts)
	WAVE myWave, myPnts
	If(0)//WaveExists($(NameOfWave(myWave)+"_pnts")))
		//Do something
		Concatenate/NP/O {$(NameOfWave(myWave)+"_pnts"),myPnts}, $(NameOfWave(myWave)+"_pnts")
		WAVE myPnts = $(NameOfWave(myWave)+"_pnts")
		Sort myPnts, myPnts
		//Remove duplicates, need to write a function to do this, I think ...
	Else
		Duplicate/o myPnts $(NameOfWave(myWave)+"_pnts")
	Endif
	Duplicate/O myPnts, $(NameOfWave(myWave)+"_Ys")
	WAVE myYs = $(NameOfWave(myWave)+"_Ys")
	
	myYs=myWave(myPnts)
End

Function RemoveBaselines(timepoints,[myPnts,InterpOrder])
// Goes through timepoints and removes the baselines
// If myPnts is not provided no points are forced
	WAVE timepoints, myPnts
	Variable InterpOrder //Used for the interpolation of the baseline
	
	If(ParamIsDefault(InterpOrder))
		InterpOrder=2
	EndIf
	
	Variable i, length = numpnts(timepoints)
	String currenttime
	
	for(i=0;i<length;i+=1)
		if(timepoints[i]<=0)
			currenttime="m"+num2istr(abs(timepoints[i]))
		else
			currenttime="p"+num2istr(timepoints[i])
		endif
		currenttime+="_withg"
		WAVE myWave = $currenttime
		If(!ParamIsDefault(myPnts))
			ForcePnts(myWave,myPnts)
		EndIf
		RemoveBaseline(myWave,interpOrder=InterpOrder)
	EndFor
End

STATIC Function devGrapher()
	//Function that graphs the differential waves and thresholds
	//Or updates said graph if it already exists
End

Function PickPoints(myWave,[DIF1Thresh,DIF2Thresh,DIF3Thresh,cutoff,extra])
//This function determines which points belong to the baseline by examining the first,
//second and third derivatives of the data
	WAVE myWave
	
	Variable DIF1Thresh, DIF2Thresh,DIF3Thresh
	
	Variable cutoff,extra//Cutoff is used for BetterSDev and extra determines whether to
	//strongly smooth the data first as an estimate for the baseline.
	
	Variable length = numpnts(myWave), i=0,j=0
	//My smoothing variables
	
	If(ParamIsDefault(extra))
		Extra=0
	EndIf
	
	Duplicate/O myWave myWave_smth
	Duplicate/O myWave myWave_tempBL
	//A first guess at what the baseline may be is a heavily smoothed version of the
	//original data
	If(Extra)
		FFB(myWave_tempBL,200,2)
		//This temporary baseline is removed from the original wave which is then smoothed
		myWave_smth=myWave-myWave_tempBL
	EndIf
	FFB(myWave_smth,10,3)//Comment this line out if the data
	//are already smoothed!
	
	//Change of plans!
	//1.) we heavily smooth the original data, extreme smoothing to remove any sharp features
	//2.) we subtract this heavily smoothed data from the original
	//3.) we pick points on this semi baseline subtracted data.
	//4.) we use these points on the original data!
	//Test seems to work, need algorithm to make point spacing more sparse..
	
	//Create the derivatives
	Differentiate myWave_smth/D=myWave_DIF1
	Differentiate myWave_DIF1/D=myWave_DIF2
	Differentiate myWave_DIF2/D=myWave_DIF3
	
	If(ParamIsDefault(cutoff))
		cutoff=0.1
	EndIf
	
	//Calculate the thresholds if the user hasn't specified them
	If(ParamIsDefault(DIF1Thresh))
		//Do something wave stats wise....
		//DIF1Thresh=1.5e-5
		DIF1Thresh=BetterSDev(myWave_DIF1,cutoff)
		//DIF1Thresh=BetterSDev2(myWave_DIF1)
	EndIf
	If(ParamIsDefault(DIF2Thresh))
		//Do something wave stats wise....
		//DIF2Thresh=2.5e-6
		DIF2Thresh=BetterSDev(myWave_DIF2,cutoff)
		//DIF2Thresh=BetterSDev2(myWave_DIF2)
	EndIf
	If(ParamIsDefault(DIF3Thresh))
		//Do something wave stats wise....
		//DIF3Thresh=5e-7
		DIF3Thresh=BetterSDev(myWave_DIF3,cutoff)
		//DIF3Thresh=BetterSDev2(myWave_DIF3)
	EndIf
	
	//Display the derivatives and thresholds for the user
	devGrapher()
	
	//Make a wave to hold the chosen points
	Make/D/O/N=(length) $(NameOfWave(myWave)+"_pnts")
	//Create a reference to this wave
	WAVE myPnts = $(NameOfWave(myWave)+"_pnts")
	
	//Add points to the wave
	For(i=0;i<length;i+=1)
		If(abs(myWave_DIF1[i])<DIF1Thresh && abs(myWave_DIF2[i])<DIF2Thresh && abs(myWave_DIF3[i])<DIF3Thresh)
			//Should I add this point to the baseline?
			myPnts[j]=i
			j+=1
		EndIf
	EndFor
	
	//Remove excess points
	Redimension/N=(j) myPnts
	
	//Refine the points
	refinePoints(myWave,myPnts)
	
	//Remove temporary waves
	KillWaves myWave_smth//,myWave_tempBL
	
	Return 0
End

STATIC Function refinePoints(myWave,myPnts)
	//Fix end points by setting them equal to myWave's end points
	//This is need to ensure that the interpolation, done later is
	//done correctly
	WAVE myWave, myPnts
	//First we need to check if the  previous algorithm already included the endpoints
	Variable V_Max, V_Min, length=numpnts(myWave)
	WaveStats/Z/Q myPnts
	
	If(V_Min!=0)//Did I include the first point?
		InsertPoints 0, 1, myPnts
	EndIf
	
	If(V_Max<(length-1))//Did I include the last one?
		InsertPoints numpnts(myPnts), 1, myPnts
		myPnts[numpnts(myPnts)]=length-1
	EndIf
	//The next part of this function will make sure that given sets of points are not too closely spaced
	//Steps:
	//1.) figure out cluster size and characteristics
	//	Beginning
	//	end
	//	number of points
	//2.) based on these characteristics
	//	we'll remove excess points so that the spline won't be interpolated through too many points
	
	Variable i=0,j=0,k=0,waveLength=numpnts(myPnts)
	
	//We'll store our cluster info in a multidimensional wave to keep track of everything
	//We'll try using a free wave as the cluster info won't be needed in the future
	Make/D/O/N=(100,2) ClusterInfo=0
	SetDimLabel 1,0,Begin, ClusterInfo
	SetDimLabel 1,1,End, ClusterInfo
	
	//Some variables to control the exit wave
	Variable clusterSpacing=15//Makes sure that the points are well separated
	
	For(i=0;i<waveLength;i+=1)
		//Points to removeare set to -1
		If((myPnts[i+1]-myPnts[i])<clusterSpacing&&(myPnts[i+1]-myPnts[i])>0)
			ClusterInfo[j][%Begin] = i
			Do
				i+=1
			While((myPnts[i+1]-myPnts[i])<clusterSpacing&&(myPnts[i+1]-myPnts[i])>0)
			ClusterInfo[j][%End] = i
			j+=1
		EndIf
	EndFor
	Redimension/N=(j,2) ClusterInfo
	
	Variable NumClusters=j,ClusterSize=0,newNumPnts=0
	Variable firstValue=0, incrementValue=0
	
	For(j=0;j<NumClusters;j+=1)
		
		ClusterSize=ClusterInfo[j][%End]-ClusterInfo[j][%Begin]+1
		newNumPnts=floor(ClusterSize/clusterSpacing)+1
		//Special (I think) case for only 1 point
		
		If(newNumPnts==1)
			
			myPnts[ClusterInfo[j][%Begin]]= (myPnts[ClusterInfo[j][%End]]+myPnts[ClusterInfo[j][%Begin]])/2
			
			For(i=(ClusterInfo[j][%Begin]+1);i<(ClusterInfo[j][%End]+1);i+=1)
				myPnts[i]=-1
			EndFor
	
		Else
			//Normal cases
			k=0
			firstValue=myPnts[ClusterInfo[j][%Begin]]
			incrementValue=(myPnts[ClusterInfo[j][%End]]-firstValue)/newNumPnts
		
			For(i=(ClusterInfo[j][%Begin]);i<(ClusterInfo[j][%End]+1);i+=1)
				If(k<newNumPnts+1)
					//Print k
					myPnts[i]=firstValue+k*incrementvalue
					k+=1
				Else
					myPnts[i]=-1
				EndIf
			EndFor
		
		EndIf
	
	EndFor
			
	Extract/O myPnts, myPnts, myPnts>=0//remove the unwanted points
End

Function/S RemoveBaseline(myWave,[interpOrder,Name])
//After the baseline points are chosen this will interpolate a spline between the points and
//remove the baseline from the data.
	Wave myWave
	Variable interpOrder
	String Name
	
	If(ParamIsDefault(interpOrder))
		interpOrder=2//Cubic Spline
	EndIf
	
	If(ParamIsDefault(Name))
		Name=NameOfWave(myWave)+"_nb"
	EndIf
	
	//My smoothing variables
	Variable cutoff=10, n=3
	
	Duplicate/FREE myWave myWave_smth
	FFB(myWave_smth,cutoff,n)//Comment this line out if the data
	//are already smoothed!
	
	Variable length = numpnts(myWave)
	
	Wave myPnts = $(NameOfWave(myWave)+"_pnts")
	If(!WaveExists(myPnts))
		Print "There was an error,"+NameOfWave(myPnts)+" doesn't exist"
		Return ""
	EndIf
	Duplicate/O myPnts, $(NameOfWave(myWave)+"_Ys")
	WAVE myYs = $(NameOfWave(myWave)+"_Ys")
	
	//Calculate points
	myYs = myWave_smth(myPnts)
	
	//Make my result waves
	Make/D/O/N=(length) $(NameOfWave(myWave)+"_bl")
	WAVE myWave_bl = $(NameOfWave(myWave)+"_bl")
	Duplicate/O myWave $Name
	WAVE myWave_nb = $Name
	
	Interpolate2/T=(interpOrder)/I=3/E=2/Y=myWave_bl myPnts, myYs//Cubic spline interpolation
	
	myWave_nb-=myWave_bl
	
	Return GetWavesDataFolder(myWave_nb, 2)
End

Function/S DoBaseline(myWave,[cutoff])
//Function that picks points and removes the baseline.
	Wave myWave
	Variable cutoff
	
	If(ParamIsDefault(cutoff))
		cutoff=0.1
	EndIf
	
	If(!WaveExists(myWave))
		DoAlert 0, "The wave you gave me does not exist"
		Return "-1"
	EndIf
	//pick the points
	pickpoints(myWave,cutoff=cutoff)
	Return RemoveBaseline(myWave)//Return a reference to the baseline subtracted data
End

Function/S DoBaseline2(myWave,num)
	//Need to clean up this function
	//Questions:
	//Can I make this generically adaptive? Can I make it so it can correctly choose the right number of passes?
	//Can I fix the naming scheme to make this less retarded?
	Wave myWave
	Variable num
	String wl=GetWavesDataFolder(myWave, 2)+";"//My wave list
	If(!WaveExists(myWave))
		DoAlert 0, "The wave you gave me does not exist"
		Return ""
	EndIf
	
	Variable i=0
	Variable myCutoff=0.1
	pickpoints(myWave,cutoff=myCutoff)
	String myName="",myName_nb=(NameOfWave(myWave)+"_nb"+num2str(i))
	RemoveBaseline(myWave,Name=myName_nb)
	wl+=GetWavesDataFolder($myName_nb, 2)+";"
	For(i=0;i<num;i+=1)
		myName=NameOfWave(myWave)+"_nb"+num2str(i)
		myName_nb=NameOfWave(myWave)+"_nb"+num2str(i+1)
		pickpoints($myName,cutoff=myCutoff+i*0.2/num)
		RemoveBaseline($myName,name=myName_nb)
		wl+=GetWavesDataFolder($(nameOfWave(myWave)+"_nb"+num2str(i+1)), 2)+";"
	EndFor
	Return wl
End

Function displayTimePnt2(nameStr)
	String nameStr
	//Error Checking
	WAVE shiftx = shiftx
	if(!WaveExists(shiftx))
		Print "Error, shiftx does not exist"
		Return -2
	EndIf
	if(!WaveExists($(nameStr+"_subG")))
		Print "Error, the subG wave does not exist"
		Return -1
	EndIf
	if(!WaveExists($(nameStr+"_withG")))
		Print "Error, the withG wave does not exist"
		Return 1
	EndIf
	PauseUpdate; Silent 1		// building window...
	Display /W=(28,44,428,252) $(nameStr+"_subG"),$(nameStr+"_withG") as nameStr
	ModifyGraph rgb($(nameStr+"_subG"))=(0,0,0),rgb($(nameStr+"_withG"))=(0,0,65535)
	//Check to see if the baseline waves exist
	if(WaveExists($(nameStr+"_withg_bl")))
		AppendToGraph $(nameStr+"_withg_bl")
		ModifyGraph rgb($(nameStr+"_withg_bl"))=(65535,0,0)
	EndIf
	if(WaveExists($(nameStr+"_withg_Ys"))&&WaveExists($(nameStr+"_withg_pnts")))
		AppendToGraph $(nameStr+"_withg_Ys") vs $(nameStr+"_withg_pnts")
		ModifyGraph rgb($(nameStr+"_withg_Ys"))=(0,0,0)
		ModifyGraph mode($(nameStr+"_withg_Ys"))=3,marker($(nameStr+"_withg_Ys"))=19
		ModifyGraph msize($(nameStr+"_withg_Ys"))=2
	EndIf
End
