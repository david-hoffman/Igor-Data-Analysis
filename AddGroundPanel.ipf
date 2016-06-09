#pragma rtGlobals=1		// Use modern global access method.
#include ":Utilities"

Menu "Macros"
	"Add Ground", AddGroundPanel()
	"-"
End

Function AddGroundPanel() : Panel
//Function that generates the panel with the attached graphs.
	
	String/G name_Timepoints="root:timepoints"
	String/G name_Groundadded="root:groundadded"
	String/G name_FitGnd="root:fit_OCPo2"
	
	Wave timepoints=$name_Timepoints
	Wave groundadded=$name_Groundadded
	Wave fitgnd=$name_FitGnd
	
	If(!WaveExists(timepoints)||!WaveExists(groundadded))
		DoAlert 0, "This panel will not work! The proper waves do not exist!"
		Return -1
	EndIf
	
	String timelist="\"_None_;"
	variable i=0
	//Build a list to show in the pop up menu
	for(i=0;i<numpnts(timepoints);i+=1)
		timelist+=num2istr(timepoints[i])+";"
	EndFor
	timelist+="\""
	
	PauseUpdate; Silent 1		// building window...
	Display/K=1/W=(0,0,600,800)/N=GAGraph
	Display/HOST=GAGraph/N=TimeTraces/W=(0,1,1,0.5)/FG=(FL,FT,FR,*)// Groundadded vs Timepoints
	Display/HOST=GAGraph/N=Groundadded/W=(0,0.5,1,1)/FG=(FL,*,FR,FB) Groundadded vs Timepoints
	ModifyGraph/W=GAGraph#Groundadded mirror=2,minor(bottom)=1,prescaleExp(left)=2,prescaleExp(bottom)=-3
	ModifyGraph/W=GAGraph#Groundadded  mode=4,marker=19
	Label/W=GAGraph#Groundadded left "Fraction of Ground Added (%)";Label/W=GAGraph#Groundadded bottom "Time (ps)"
	
	NewPanel/HOST=GAGraph/N=GAPanel/EXT=0/W=(0,0,300,200) as "Time to add some ground!"
	
	//PopupMenu whichTime noProc
	PopupMenu whichTime,pos={60,10},size={150,20},proc=chooseTimePoint,title="\f01Choose a time point"
	PopupMenu whichTime,mode=1,help={"Choose a timepoint to which you want to add ground"}
	Execute "PopupMenu whichTime,value="+timelist
	
	//PopupMenu for the GroundAdded
	PopupMenu whichGA,pos={30,40},size={180,20},proc=chooseGroundAdded,title="\f01Choose a GroundAdded Wave"
	PopupMenu whichGA,mode=1,help={"Choose a wave to store the groundadded"},value=WaveList("GroundAdded*",";","")
	
	//PopupMenu for the fit_gnd
	PopupMenu whichGND,pos={30,70},size={180,20},proc=chooseFitGnd,title="\f01Choose the fit ground"
	PopupMenu whichGND,mode=1,help={"Choose the fit ground"},value=WaveList("fit*",";","")
	
	//Slider
	Slider groundToAdd,pos={20,100},size={250,45},limits={0,1,0},vert= 0
	Slider groundToAdd, ticks=20, proc=AddGroundSlider, title="Ground to Add"
End

STATIC Function plotTimePnt(timeStr)
	String timeStr
	String currenttime
	
	Variable startUp = CmpStr(TraceNameList("GAGraph#TimeTraces",";",1),"")
	//Get axis data first
	If(startUp!=0)
		GetAxis/W=GAGraph#TimeTraces/Q bottom
	EndIf
	
	//clear the graph first
	KillWindow GAGraph#TimeTraces
	
	if(str2num(timeStr)<=0)
		currenttime="m"+num2istr(abs(str2num(timeStr)))
	else
		currenttime="p"+num2istr(str2num(timeStr))
	endif
	
	WAVE withg = $(currenttime+"_withg")
	WAVE subg = $(currenttime+"_subg")
	
	//Error Checking
	if(!WaveExists(subg))
		Print "Error, the subG wave does not exist"
		Return -1
	EndIf
	if(!WaveExists(withg))
		Print "Error, the withG wave does not exist"
		Return 1
	EndIf
	
	String titleString=""
	If(abs(str2num(timeStr)) >= 1000)
		titlestring += num2str(str2num(timeStr)/1000)+ " ps"
	Else
		titlestring+= timeStr + " fs"
	Endif
	titlestring+=" subg/withg"
	
	DoWindow/T GAGraph, titlestring
	Display/HOST=GAGraph/N=TimeTraces/W=(0,1,1,0.5)/FG=(FL,FT,FR,*) subg,withg vs shiftx
	ModifyGraph/W=GAGraph#TimeTraces rgb[0]=(0,0,0),rgb[1]=(0,0,65535)
	ModifyGraph/W=GAGraph#TimeTraces mirror=2,minor(bottom)=1,prescaleExp(left)=3
	Label/W=GAGraph#TimeTraces left "Raman Gain (mOD)"
	Label/W=GAGraph#TimeTraces bottom "Raman Shift (cm\\S-1\\M)"
	SetAxis/W=GAGraph#TimeTraces/A=2 left
	If(startUp!=0)
		SetAxis/W=GAGraph#TimeTraces bottom V_Min,V_Max
	EndIf
End

Function chooseTimePoint(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	
	String/G name_Timepoints
	String/G name_Groundadded
	
	Wave timepoints=$name_Timepoints
	Wave groundadded=$name_Groundadded
	
	switch( pa.eventCode )
		case 2: // mouse up
			If(CmpStr(pa.UserData,pa.popStr)==0)
				//The user has chosen the same timepoint, do nothing
				break
			EndIf
			
			//Remove customization for previous time point
			FindValue/v=(str2num(pa.UserData)) timepoints
			ModifyGraph/W=GAGraph#GroundAdded rgb($name_Groundadded[-V_Value-1])=(0,0,65535)
			ModifyGraph/W=GAGraph#GroundAdded msize($name_Groundadded[-V_Value-1])=10
			ModifyGraph/W=GAGraph#GroundAdded marker($name_Groundadded[-V_Value-1])=42
			
			If(CmpStr(pa.popStr,"_None_")==0)
				//If the user has changed to None
				Break
			EndIf
			
			plotTimePnt(pa.popStr)
			FindValue/v=(str2num(pa.popStr)) timepoints
			Slider groundToAdd, win=GAGraph#GAPanel, value= groundadded[V_value]
			//Cursor/W=GAGraph#GroundAdded/P A, $StringFromList(0,TraceNameList("GAGraph#GroundAdded",";",1)), V_Value
			ModifyGraph/W=GAGraph#GroundAdded rgb($name_Groundadded[V_Value])=(0,0,65535)
			ModifyGraph/W=GAGraph#GroundAdded msize($name_Groundadded[V_Value])=10
			ModifyGraph/W=GAGraph#GroundAdded marker($name_Groundadded[V_Value])=42
			
			//Store new timepoint in the userdata for this popup for use by other functions
			pa.UserData = pa.popStr
			break
		case -1: // control being killed
			break
	endswitch
	

	return 0
End

Function chooseGroundAdded(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	
	String/G name_Timepoints
	String/G name_Groundadded
	String/G name_FitGnd
	
	Wave timepoints=$name_Timepoints
	
	
	switch( pa.eventCode )
		case 2: // mouse up
			//Update all "withg" waves to reflect the new ground added
			AddAllGround($name_Timepoints,$name_FitGnd,$name_Groundadded)
			If(CmpStr(pa.UserData,pa.popStr)==0)
				//The user has chosen the same groundadded, do nothing
				break
			EndIf
			
			//Update user data and path to groundadded
			pa.UserData = pa.popStr
			name_Groundadded = pa.popStr
			
			Wave groundadded=$name_Groundadded//Have to do assignment here!
			If(!WaveExists($name_Timepoints))
				//There is no timepoints wave chosen yet
				//Do nothing and return control
				Return -1
			EndIf
			
			//Kill the subwindow so that we can recreate it with the new ground added
			KillWindow GAGraph#Groundadded
			
			//Create a new subwindow with the same name
			Display/HOST=GAGraph/N=Groundadded/W=(0,0.5,1,1)/FG=(FL,*,FR,FB) Groundadded vs Timepoints
			ModifyGraph/W=GAGraph#Groundadded mirror=2,minor(bottom)=1,prescaleExp(left)=2,prescaleExp(bottom)=-3
			ModifyGraph/W=GAGraph#Groundadded  mode=4,marker=19
			Label/W=GAGraph#Groundadded left "Fraction of Ground Added (%)";Label/W=GAGraph#Groundadded bottom "Time (ps)"
			
			//Make sure the slider is in the correct position
			FindValue/v=(str2num(GetUserData("GAGraph#GAPanel","whichTime",""))) timepoints
			Slider groundToAdd, win=GAGraph#GAPanel, value= groundadded[V_value]
			
			//Print "Cursor/W=GAGraph#GroundAdded/P A, "+name_Groundadded+", "+num2str(V_Value)
			//Print "A list of traces: "+TraceNameList("GAGraph#GroundAdded","\t",1)
			
			//Cursor/W=GAGraph#GroundAdded/P A, $StringFromList(0,TraceNameList("GAGraph#GroundAdded",";",1)), V_Value
			ModifyGraph/W=GAGraph#GroundAdded rgb($name_Groundadded[V_Value])=(0,0,65535),msize($name_Groundadded[V_Value])=10,marker($name_Groundadded[V_Value])=42
			break
		case -1: // control being killed
			break
	endswitch
	

	return 0
End

Function chooseFitGnd(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	
	String/G name_Timepoints
	String/G name_Groundadded
	String/G name_FitGnd
	
	switch( pa.eventCode )
		case 2: // mouse up
			
			If(CmpStr(pa.UserData,pa.popStr)==0)
				//The user has chosen the same groundadded, do nothing
				break
			EndIf
			
			pa.UserData = pa.popStr
			name_FitGnd = pa.popStr
			
			If(!WaveExists($name_Timepoints))
				//There is no timepoints wave chosen yet
				//Do nothing and return control
				Return -1
			EndIf
			//Update all "withg" waves to reflect the new ground added
			AddAllGround($name_Timepoints,$name_FitGnd,$name_Groundadded)
			
			break
		case -1: // control being killed
			break
	endswitch
	

	return 0
End

Function AddGroundSlider(sa) : SliderControl
	STRUCT WMSliderAction &sa
	
	String/G name_Timepoints
	String/G name_Groundadded
	String/G name_FitGnd
	
	Wave timepoints=$name_Timepoints
	Wave groundadded=$name_Groundadded
	
	If(!WaveExists($name_Timepoints))
		//The timepoints wave hasn't been chosen yet
		//Return control to calling function
		Return -1
	EndIf
	
	String S_myTime
	Variable myTime
	
	
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				S_myTime = GetUserData("GAGraph#GAPanel","whichTime","")
				If(CmpStr(S_myTime,"_None_")==0)
					//Do nothing
					Break
				EndIf
				myTime = str2num(S_myTime)
				Variable curval = sa.curval
				addground3(curval)
			endif
			break
	endswitch

	return 0
End

Static Function AddGround3(amount)
//this adds a set amount of ground to particular timepoint and updates groundadded accordingly
	Variable amount
	
	String/G name_Timepoints
	String/G name_Groundadded
	String/G name_FitGnd
	
	Wave timepoints=$name_Timepoints
	Wave groundadded=$name_Groundadded
	WAVE ground=$name_FitGnd
	
	Variable timepoint = str2num(GetUserData("GAGraph#GAPanel","whichTime",""))
	
	FindValue/v=(timepoint) timepoints
	
	String currenttime
	
	if(timepoint<=0)
		currenttime="m"+num2istr(abs(timepoint))
	else
		currenttime="p"+num2istr(timepoint)
	endif
	WAVE withg = $(currenttime+"_withg")
	WAVE subg = $(currenttime+"_subg")
	If(WaveExists(withg)&&WaveExists(subg))
		withg=subg+amount*ground
		groundadded[V_value]= amount
	Else
		DoAlert/T="Wave does not exist " 0, currenttime + " does not exist, please check your parameters"
	EndIf
End
