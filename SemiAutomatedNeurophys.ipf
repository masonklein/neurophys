 #pragma rtGlobals=1		// Use modern global access method.
//
// (updated 2017)
//
// These functions are designed for 3D confocal images that have 
// been pre-processed in NIS Elements before doing anything in Igor.  Specifically, the original 3D movie will have 
// been cropped around the region of interest in 
// order to make files smaller and the procedure run faster.  This is not strictly required
//
// Also assumes an XY projection (aka z-projection) was done and saved as a single TIFF sequence.  Igor can load these directly
// This means Igor won't need to keep them in memory, and we can safely
// use 14-bit images instead of the usual 8-bit ones with weird scaling issues.  
//
// Files for the raw 3D movie are assumed to have the Elements naming convention, which is NAME_t##_z##.  Presumably the
// number of digits for t and z can vary.  That is, we operate on 3D movies where each slice is a separate TIFF file.
// 
// Finally, there should be a text file with meta data from the movie recording (format "###_time"), and a stimulus
// text file (typically temperature, although this can me modified for anything time-varying), with format "###_temp".
// 
// To run the program, type SpinDiskND2("name") into the command line in Igor, where "name" is whatever prefix
// you want to apply to all the output waves (quote marks *should* be included).  Follow the on-screen prompts.  


Function SpinDiskND2(prefix)

	string prefix // the only input; all waves associated with movie will start with this string
	//
	// Initialize parameters if this is a new prefix (does nothing otherwise)
	SDND2_initParams(prefix)
	//
	// Display panel where user chooses the data path, points at image files, and
	// loads timing, z, and temperature data
	SDND2_getParams(prefix)
	PauseForUser $(winname(0,64))
	wave paramWave = $(prefix+"params")
	//
	// Declare some local variables to point at some parameters
	variable numObjects = paramWave[1]
	variable numStacks = paramWave[6]
	variable firstTime = paramWave[7]
	variable lastTime = paramWave[8]
	variable tempTimeStep = paramWave[12]
	variable numToAvg = paramWave[2]
	// 
	// Make graphs (empty at first) for temperature and neuron signals
	if(!waveexists($(prefix+"signal"+num2str(numToAvg))))
		Make/O/N=(numObjects,numStacks) $(prefix+"signal"+num2str(numToAvg)) = NaN
		print "HI -- it made the signal wave"
	endif
	wave signal = $(prefix+"signal"+num2str(numToAvg))
	wave temperature = $(prefix+"temp0")
	SetScale/I y,firstTime,lastTime,"s",signal
	SetScale/P x,0,tempTimeStep,"s",temperature
	//
	variable j
	Display temperature[][2]
	ModifyGraph axisEnab(left) = {0,0.4}, nticks(bottom)=20
	Label bottom "Time"
	Label left "Temp. (C)"
	for(j=0; j<numObjects; j+=1)
		AppendToGraph/C=(SDND2_TenColors("red",j),SDND2_TenColors("green",j),SDND2_TenColors("blue",j))/L=left2 signal[j][]
	endfor
	Label left2 "signal (a.u.)"
	ModifyGraph lsize=2,lblPos(left2)=49,axisEnab(left2)={0.5,1},freePos(left2)=0
	SetAxis left2 0,*
	Legend
	// 
	// Let user move the graph somewhere before starting the path/signal computation
	string graphName = winname(0,1)
	SDND2_OKpanel()
	PauseForUser $(winname(0,64)), $graphName
	//
	// User chooses the locations of objects of interest, and signal levels are computed
	SDND2_GetPathAndSignals(prefix)
	SetScale/I y,firstTime,lastTime,"s",signal
	//
	// New part added 09/14/2014
	// basically just display the temperature (actual time also, not scaled) on a separate graph
	// so the user can use it to choose the "zero point" value for the time scale
	// (typically this is when the SV first starts to change for a sine wave)
	// should also calculate the typical value and print it (perhaps on the graph itself)
	SDND2_OtherTempGraph(temperature)
	//
	// Even newer part added 09/17/2014
	// just need to compute the average for the *last* wave (which will typically be the background box)
	// and print it somewhere (either on the graph or elsewhere)
	make/N=(numStacks) $("bgTemporary") = signal[numObjects-1][p]
	wave bgTemp = $("bgTemporary")
	print faverage(bgTemp)
	killwaves bgTemp
	// 
End // of the main SpinDiskND2 function



// Assigns default values to parameters for the movie analysis -- if they have already
// been defined for this prefix, the function does nothing
Function SDND2_initParams(prefix)

	string prefix
	if(!waveexists($(prefix+"params")))
		//
		make/O/N=15 $(prefix+"params")
		wave paramWave = $(prefix+"params")
		//
		paramWave[0] = NaN		// data path # [nan]	
		paramWave[1] = 1		// number of objects
		paramWave[2] = 100 		// number of brightest pixels to average
		//
		paramWave[3] = 512		// # rows
		paramWave[4] = 512		// # columns
		paramWave[5] = 0		// # planes
		paramWave[6] = 0		// # stacks
		paramWave[7] = 0.3		// first time [0.3 seconds]
		paramWave[8] = 300		// last time [300 seconds]
		paramWave[9] = 5		// time step [5 seconds]
		paramWave[10] = 0		// first stack 
		paramWave[11] = 100		// last stack
		// 
		paramWave[12] = 0.2		// time between temperature points
		//	
		paramWave[13] = 0.222	// distance per pixel (XY) [222 nm]
		paramWave[14] = 0.400	// distance per pixel (Z) [400 nm]
	endif
End // of SDND2_initParams function 

// Panel where user sets parameters of the movie
Function SDND2_getParams(prefix) : Panel
	//
	// input prefix string, and a corresponding global variable to use within button functions
	string prefix
	string/G gPrefix = prefix
	//
	// parameter wave, which should already exist before running this function
	wave/Z paramWave = $(prefix+"params")
	if(!waveexists($(prefix+"params")))
		return 0
	endif
	//
	// Make the panel
	NewPanel/W=(760,100,1160,740)
	string panelTitle = prefix + " PARAMETERS"
	TitleBox tb0, fixedsize=1,fstyle=1,fsize=20,pos={15,8},size={370,32},title=panelTitle					// 40
	//
	// Data paths and objects
	TitleBox tb1,fsize=20,fixedsize=1,fstyle=1,pos={20,50},size={300,30},title="PATHS AND OBJECTS"	
	SetVariable sv0,fsize=20,pos={30,90},size={150,30},value=paramWave[0],title="data path",proc=SDND2_datapathPROC
	SetVariable sv1, fsize=20,pos={30,130},size={150,30},value=paramWave[1],title="# objects" 
	SetVariable sv2, fsize=20,pos={30,170},size={210,30},value=paramWave[2],title="# pixels to avg."		// 200
	//
	// Movie file information
	TitleBox tb2,fsize=20,fixedsize=1,fstyle=1,pos={20,210},size={300,30},title="MOVIE FILE INFO."	
	SetVariable sv3, fsize=20,pos={30,250},size={160,30},value=paramWave[5],title="# planes      "
	SetVariable sv4, fsize=20,pos={30,280},size={160,30},value=paramWave[6],title="# stacks      "
	SetVariable sv5, fsize=20,pos={30,310},size={190,30},value=paramWave[7],title="first time (s) "
	SetVariable sv6, fsize=20,pos={30,340},size={190,30},value=paramWave[8],title="last time (s) "		// 370
	Button button0, fsize=20,pos={230,280},size={120,30},proc=SDND2_loadmetaPROC,title="LOAD"
	//
	// Range of interest
	TitleBox tb3,fsize=20,fixedsize=1,fstyle=1,pos={20,380},size={300,30},title="RANGE OF INTEREST"
	SetVariable sv7, fsize=20,pos={30,420},size={150,30},value=paramWave[10],title="first stack"
	SetVariable sv8, fsize=20,pos={30,450},size={150,30},value=paramWave[11],title="last stack"			// 480
	//
	// Temperature
	TitleBox tb4,fsize=20,fixedsize=1,fstyle=1,pos={20,490},size={300,30},title="TEMPERATURE"
	SetVariable sv9, fsize=20,pos={30,530},size={180,30},value=paramWave[12],title="time step (s)" 		// 560
	Button button1, fsize=20,pos={230,530},size={120,30},proc=SDND2_loadtempPROC,title="LOAD"
	//
	// OK button 
	Button button2, fsize=20,pos={140,600},size={120,30},proc=SDND2_okbuttonPROC,title="OK"
	


End	// of SDND2_getParams function (panel)

Function SDND2_datapathPROC(ctrlName,varNum,varStr,varName) : SetVariableControl
	string ctrlName
	variable varNum	// the value in the SetVariable box, which is what we use below
	string varStr
	string varName
	//
	// Declare the data path names and get path information (PathInfo does not cause an error if the path doesn't exist)
	string dataPathName = "path" + num2str(varNum)
	string dataPathNameXY = dataPathName + "xy"
	//string dataPathNameYZ =  dataPathName + "yz"
	//string dataPathNameXZ = dataPathName + "xz"
	PathInfo $dataPathName
	//
	svar gPrefix
	wave paramWave = $(gPrefix+"params")
	//
	// Check if the data path already exists, and prompt the user to create it if it doesn't
	// (Also will create the three paths for the projection images)
	variable i
	if(V_flag==0)
		string/G gFileName, gFileNameXY//, gFileNameYZ, gFileNameXZ
		//
		print "select one of the 4D movie files"
		ImageLoad/Q/N=$"tempImage"/T=tiff
		wave tempImageWave = $("tempImage")
		paramWave[3] = dimsize(tempImageWave,0)
		paramWave[4] = dimsize(tempImageWave,1)
		KillWaves/Z tempImageWave
		variable splitPoint = strlen(S_fileName)-strsearch(S_fileName,"_t",0)
		for(i=0; i<splitPoint; i+=1)
			S_fileName = RemoveEnding(S_fileName)
		endfor
		gFileName = S_fileName
		NewPath/Q $dataPathName, S_path
		//
		print "select one of the XY projection image files"
		ImageLoad/Q/RTIO/P=$dataPathName/T=tiff
		splitPoint = strlen(S_fileName)-strsearch(S_fileName,"_t",0)
		for(i=0; i<splitPoint; i+=1)
			S_fileName = RemoveEnding(S_fileName)
		endfor
		gFileNameXY = S_fileName
		NewPath/Q $dataPathNameXY, S_path
		//
//		print "select one of the YZ projection image files"
//		ImageLoad/Q/RTIO/P=$dataPathName/T=tiff
//		splitPoint = strlen(S_fileName)-strsearch(S_fileName,"_t",0)
//		for(i=0; i<splitPoint; i+=1)
//			S_fileName = RemoveEnding(S_fileName)
//		endfor
//		gFileNameYZ = S_fileName
//		NewPath/Q $dataPathNameYZ, S_path
//		//
//		print "select one of the XZ projection image files"
//		ImageLoad/Q/RTIO/P=$dataPathName/T=tiff
//		splitPoint = strlen(S_fileName)-strsearch(S_fileName,"_t",0)
//		for(i=0; i<splitPoint; i+=1)
//			S_fileName = RemoveEnding(S_fileName)
//		endfor
//		gFileNameXZ = S_fileName
//		NewPath/Q $dataPathNameXZ, S_path
		//
	endif
End
//
// BUTTONS from the above parameter acquisition panel function...
//
// Loads meta data (i.e., times and z-positions for all images)
Function SDND2_loadmetaPROC(ctrlName) : ButtonControl
	string ctrlName
	//
	svar gPrefix
	//
	SDND2_LoadMeta(gPrefix)
	//
	DrawText 265,330,"LOADED!"	
End
//
// Loads temperature data (SV, PV, and time)
Function SDND2_loadtempPROC(ctrlName) : ButtonControl
	string ctrlName
	//
	svar gPrefix
	//
	SDND2_LoadTemp(gPrefix)
	//
	DrawText 265,580,"LOADED!"
End
//
// OK button that closes a panel
Function SDND2_okbuttonPROC(ctrlName) : ButtonControl
	string ctrlName
	DoWindow/K $(winname(0,64))
End
//
// FUNCTIONS from the above parameter acquisition panel function...
//
// Loads meta data into the parameter wave -- this comes from a text file
// 	generated in Elements and saved
Function SDND2_LoadMeta(prefix)
	//
	string prefix
	//
	wave paramWave = $(prefix+"params")
	//
	string pathStr = "path" + num2str(paramWave[0])
	//
	if(!waveexists($(prefix+"meta0")))
		LoadWave/O/A=$(prefix+"meta")/J/P=$pathStr/Q/L={0,1,0,1,3}/M
	endif
	wave meta = $(prefix+"meta0")
	//
	// NEW: 06/25/2014
	// Add a part that cleans up the meta file, which now saves in an annoying way that
	// stores x and y coordinates interleaved with the imaging stacks.  Practically speaking,
	// we want to delete rows in this wave where the 2nd and 3rd columns are NaN
	variable numPoints = dimsize(meta,0)
	variable k, index
	
	for(k=numPoints; k>0; k-=1)		
		index = k-1
		if(!(meta[index][1]>0))
			DeletePoints/M=0 index,1,meta
		endif	
	endfor
	//
	// Figure out the number of planes per stack, spacing between planes, 
	// and the number of stacks
	variable numPlanes
	variable i
	for(i=1; i<dimsize(meta,0); i+=1)
		if(meta[i][2]<meta[i-1][2])
			numPlanes = i
			break			
		endif	
	endfor
	paramWave[5] = numPlanes					// planes per stack
	paramWave[14] = meta[1][2] - meta[0][2]		// distance between planes
	paramWave[6] = dimsize(meta,0) / numPlanes	// number of stacks
	//
	paramWave[10] = 0							// first stack number
	paramWave[11] = paramWave[6] - 1			// last stack number
	//
	// Make a reduced wave that holds the average time for each stack instead of 
	// the time for each plane/slice 
	Make/O/N=(paramWave[6]) $(prefix+"time0")
	wave timeWave = $(prefix+"time0")
	variable j
	variable runningTotal
	for(i=0; i<paramWave[6]; i+=1)
		runningTotal=0
		for(j=i*paramWave[5]; j<(i+1)*paramWave[5]; j+=1)
			runningTotal += meta[j][0]
		endfor
		timeWave[i] = runningTotal / paramWave[5]	
	endfor
	paramWave[7] = timeWave[0]
	paramWave[8] = timewave[paramWave[6]-1]
	paramWave[9] = SDND2_MeanWaveStep(timeWave)
	//
End // of SDND2_LoadMeta function
//
Function SDND2_LoadTemp(prefix)
	string prefix
	//
	wave paramWave = $(prefix+"params")
	string pathStr = "path" + num2str(paramWave[0])
	//
	// Load the temperature values and times from the text file
	if(!waveexists($(prefix+"temp0")))
		LoadWave/O/A=$(prefix+"temp")/G/M/P=$pathStr/Q
		wave temp = $(prefix+"temp0")
		variable initialTemp = temp[0][0]
		temp[][0] -= initialTemp
	else
		wave temp = $(prefix+"temp0")
	endif
	//
	// Compute the average time spacing between measurements
	// (NOTE that the user should manually change it back if it 
	//	gives something nonsensical due to older versions of the 
	//	LabVIEW program allowing the time to reset in the text file)
	variable numIndices = dimsize(temp,0)
	variable i
	variable runningTotal = 0
	for(i=0;i<(numIndices-1);i+=1)
		runningTotal += temp[i+1][0] - temp[i][0]
	endfor
	paramWave[12] = runningTotal / (numIndices - 1)
	//	
End // of SDND2_LoadTemp function
//
Function SDND2_MeanWaveStep(aWave)
	wave aWave		// needs to be a numerical wave
	//
	variable numIndices = numpnts(aWave)
	//
	variable i
	variable runningTotal = 0
	for(i=0;i<(numIndices-1);i+=1)
		runningTotal += aWave[i+1] - aWave[i]
	endfor
	//
	variable output = runningTotal / (numIndices - 1)
	return output
End	// of MeanWaveStep function




//
// Updated/modified version of GetPath4 from the SpinDisk4.ipf procedure
Function SDND2_GetPath(prefix)
	string prefix
	//
	// Global variables
	variable/G gUpdate=1, gBack=0, gLines=0, gIgnore=0, gJump=0, gView=0
	svar gFileName, gFileNameXY, gFileNameYZ, gFileNameXZ
	variable defaultRad = 10
	//
	// Parameter wave
	wave paramWave = $(prefix+"params")
	variable numObjects = paramWave[1]
	variable numRows = paramWave[3]
	variable numCols = paramWave[4]
	variable numPlanes = paramWave[5]	
	variable numStacks = paramWave[6]
	variable startStack = paramWave[10]
	variable endStack = paramWave[11]
	variable numStackDigits = strlen(num2str(numStacks))
	//
	// Data paths
	string dataPathName = "path" + num2str(paramWave[0])
	string dataPathNameXY = dataPathName + "xy"
	string dataPathNameYZ =  dataPathName + "yz"
	string dataPathNameXZ = dataPathName + "xz"
	//
	// Make a new path with default values, if it doesn't already exist
	// Path has seven variables per object...
	// [0] = horizontal position
	// [1] = vertical position
	// [2] = horizontal radius
	// [3] = vertical radius
	// [4] = rotation angle
	// [5] = z-position
	// [6] = z "radius" 
	// Then the second object will take up indexes 7 through 13, and so on
	variable i
	if( !waveexists($(prefix+"path")) )
		Make/O/N=(7*numObjects,numStacks) $(prefix+"path") = 0
		wave path = $(prefix+"path")
		//
		for(i=0; i<numObjects; i+=1)
			path[7*i][startStack] = round(numRows/2)
			path[7*i+1][startStack] = round(numCols/2)
			path[7*i+2,7*i+3][startStack] = defaultRad
			path[7*i+4][startStack] = 0
			path[7*i+5][startStack] = round(numPlanes/2)
			path[7*i+6][startStack] = round(numPlanes/10)		
		endfor		
	else
		wave path = $(prefix+"path")
	endif
	//
	// Redimension the path wave if needed
	if(dimsize(path,0)/7 != numObjects)
		Redimension/N=(7*numObjects,-1) path
	endif
	//
	// Wave used to hold onto path values when stacks are skipped or jumped to
	Make/O/N=(7*numObjects,numStacks) $(prefix+"prevPath")
	wave prevPath = $(prefix+"prevPath")
	prevPath = path[p][startStack]
	//
	// Initialize projection images
	ImageLoad/N=$"xFrame"/O/P=$dataPathNameYZ/Q/T=tiff (gFileNameYZ+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	ImageLoad/N=$"yFrame"/O/P=$dataPathNameXZ/Q/T=tiff (gFileNameXZ+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	ImageLoad/N=$"zFrame"/O/P=$dataPathNameXY/Q/T=tiff (gFileNameXY+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	//
	// Build the graph to display / update as path is built
	variable graphSizeScale = 4
	//variable graphSizeScale = floor(400/numCols)
	NewImage zFrame; ShowInfo
	Appendimage/T=T2/L xFrame
	Appendimage/T/L=L2 yFrame
	ModifyGraph width=(graphSizeScale*numRows),height=(graphSizeScale*numCols)	
	ModifyGraph nticks=10,minor=1,freePos=0, fSize=8,btLen=3, mirror=0, lblPos=13, tlOffset=-2.00
	ModifyGraph tkLblRot(L2)=90
	SetAxis/A/R L2
	variable/C zProjFrac = SDND2_FindZProjFrac(prefix)
	ModifyGraph axisEnab(top)={0.00,real(zProjFrac)},axisEnab(left)={(1-imag(zProjFrac)),1.00}
	ModifyGraph axisEnab(L2)={0.00,(1-imag(zProjFrac)-0.05)},axisEnab(T2)={(real(zProjFrac)+0.05),1.00}
	string graphName = winName(0,1)
	//
	// Loop through the stacks (projection views) and let the user indicate the paths for the objects of interest
	variable j, k
	for(i=startStack; i<endStack; i+=1)
		//
		// Update the projections on the graph
		ImageLoad/N=$"xFrame"/O/P=$dataPathNameYZ/Q/T=tiff (gFileNameYZ+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		ImageLoad/N=$"yFrame"/O/P=$dataPathNameXZ/Q/T=tiff (gFileNameXZ+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		ImageLoad/N=$"zFrame"/O/P=$dataPathNameXY/Q/T=tiff (gFileNameXY+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		DoUpdate
		//
		// For each object, set the path values for this stack equal to the previous stack 
		// if the current stack's path values are blank (as would be the case for a previously
		// partially-completed path, or empty path, or part of a path with a new object added, etc.)
		// (Note that "previous stack" will mean the previous stack that doesn't have 0's or NaN's)
		for(j=0; j<numObjects; j+=1)
			if(path[7*j][i]==0 || path[7*j][i]==NaN)
				for(k=0; k<7; k+=1)
					path[7*j+k][i] = prevPath[7*j+k]
				endfor
			endif		
		endfor
		//
		// Do-while loop through the panel that draws outlines of the ROIs onto the projection graph
		// (this continues looping as long as gUpdate = 1)
		gUpdate = 1
		gBack = 0
		gLines = 0
		gView = 0
		variable viewBefore
		//
		do
			// Place cursors on the z-projection
			for(j=0; j<numObjects; j+=1)
				cursor/C=(65535,65535,0)/H=(gLines)/I/P/S=1 $(num2char(65+j)) $"zFrame" path[7*j][i], path[7*j+1][i]
			endfor
			//
			// Draw ovals on the z-projection (light blue hollow shapes)
			SetDrawLayer ProgFront
			SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=top,ycoord=left,save
			for(j=0; j<numObjects; j+=1)
				SetDrawEnv translate=path[7*j][i],path[7*j+1][i],rotate=path[7*j+4][i],rsabout
				DrawOval path[7*j][i]-path[7*j+2][i], path[7*j+1][i]-path[7*j+3][i], path[7*j][i]+path[7*j+2][i], path[7*j+1][i]+path[7*j+3][i]			
			endfor
			//
			// Draw rectangles onto y-projection
			SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=top,ycoord=L2,save
			for(j=0; j<numObjects; j+=1)
				DrawRect path[7*j][i]-path[7*j+2][i], path[7*j+5][i]-path[7*j+6][i], path[7*j][i]+path[7*j+2][i], path[7*j+5][i]+path[7*j+6][i]
			endfor
			//
			// Draw rectangles onto x-projection
			SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=T2,ycoord=left,save
			for(j=0; j<numObjects; j+=1)
				DrawRect path[7*j+5][i]-path[7*j+6][i], path[7*j+1][i]-path[7*j+3][i], path[7*j+5][i]+path[7*j+6][i], path[7*j+1][i]+path[7*j+3][i]
			endfor
			//
			// Display the panel, which allows the user to move cursors and change ROI 
			// size in all three dimensions
			// (Also, note the value of the "view" global variable before the panel is displayed)
			viewBefore = gView
			SDND2_PathPanel(prefix,i)
			PauseForUser $(winname(0,64)),$graphName
			//
			// Save the x and y cursor positions if the UPDATE button was pressed
			// (the other settings are saved to the path wave in real time)
			// (And don't do this if the user has just turned the VIEW back on)
			if(gUpdate == 1 && (gView==0 && viewBefore!=1) ) 
				for(j=0; j<numObjects; j+=1)
					path[7*j][i] = pcsr($(num2char(65+j)))
					path[7*j+1][i] = qcsr($(num2char(65+j))) 
				endfor
			endif
			//
			// Erase the ovals and rectangles before the next time through the loop
			SetDrawLayer/K ProgFront		
		while(gUpdate)
		//
		// If the user pressed the "ignore" button, then set the path values for this stack to NaN
		if(gIgnore==1)
			path[][i] = NaN				
			gIgnore=0
		else
			prevPath = path[p][i]
		endif
		// If the user pressed the "back" button, then set gJump to the previous stack
		if(gBack==1)
			gJump = i - 1
		endif
		//
		// If the user selected a different stack to go to (instead of the next stack), adjust
		// the loop index (i) accordingly
		// (without changing anything, gJump is just i+1 and there is no change here)
		i += gJump - i - 1
		//
	endfor
	//
	// Get rid of waves made inside this function that we don't keep using
	DoWindow/K $graphName
	KillWaves/Z xFrame, yFrame, zFrame
	KillWaves/Z prevPath
	// 
End // of SDND2_GetPath
//
//
Function/S SDND2_num2strNdigit(inputNumber,N)
	variable inputNumber
	variable N
	//
	// If no extra zeros are needed at the beginning, just use the usual num2str function
	if( inputNumber > (10^(N-1)) )
		return num2str(inputNumber)
	endif
	//
	string outputString = ""
	variable i
	//	
	// Add as many zeros as needed to the beginning of the output string
	for(i=1; i<N; i+=1)	
		if(inputNumber < 10^i)
			outputString += "0"
		endif			
	endfor
	//
	outputString += num2str(inputNumber)	
	return outputString
End // of SDND2_num2strNdigit function
//
// Modified version of FindZProjFraction from SpinDisk4.ipf
Function/C SDND2_FindZProjFrac(prefix)
	string prefix
	//
	//
	// Declare parameter wave, meta waves, and extract relevant numbers from them
	wave paramWave = $(prefix+"params")
	variable numRows = paramWave[3]
	variable numCols = paramWave[4]
	variable numPlanes = paramWave[5]
	variable xyStep = paramWave[13]
	variable zStep = paramWave[14]	
	//
	variable horzZstep =(0.95)*(xyStep*numRows)/(xyStep*numRows + zStep*numPlanes)
	variable vertZstep = (0.95)*(xyStep*numCols)/(xyStep*numCols + zStep*numPlanes)	
	variable/C output = cmplx(horzZstep,vertZstep)
	//
	return output
End // of SDND2_FindZProjFrac function
//
//
// Panel for when user is generating the path:
// Function that creates a panel with buttons that adjusts the settings for each object when
// the user is creating the "path" wave from projection images
Function SDND2_PathPanel(prefix,stackNum) : Panel
	string prefix
	variable stackNum
	string/G gPrefix = prefix
	//
	// Declare parameter wave and the path wave
	wave/Z paramWave = $(prefix+"params")
	wave/Z path = $(prefix+"path")
	if( !waveexists(paramWave) || !waveexists(path) )
		print "either parameter wave or path wave does not exist"
		return 0
	endif
	variable numObjects = paramWave[1]	
	variable numStacks = paramWave[6]
	nvar gJump, gIndex
	gIndex = stackNum
	gJump = stackNum+1
	//
	// Build the panel and write the title at the top (with the current stack number)
	variable panelHt = 255 + (100*numObjects)
	NewPanel/W=(1110,100,1410,(100+panelHt))
	string panelTitle = "PATH FOR\r" + prefix + "\r(stack " + num2str(stackNum) + " of " + num2str(numStacks-1) + ")"
	TitleBox title0, fixedsize=1,fstyle=1,fsize=20,pos={15,8},size={270,82},title=panelTitle
	//
	// Place controls for oval radii, rotation angle, and z-slab properties
	variable i
	for(i=0; i<numObjects; i+=1)
		SetVariable $("radx"+num2str(i)), fsize=20,pos={20,100+100*i},size={125,30}, value = path[7*i+2][stackNum], title=(num2char(65+i) + ": rad x")
		SetVariable $("rady"+num2str(i)), fsize=20,pos={20,130+100*i},size={125,30}, value = path[7*i+3][stackNum], title=(num2char(65+i) + ": rad y")
		SetVariable $("rot"+num2str(i)), fsize=20,pos={20,160+100*i},size={125,30}, value = path[7*i+4][stackNum], title=(num2char(65+i) + ": rotate")
		SetVariable $("posz"+num2str(i)), fsize=20,pos={155,100+100*i},size={125,30}, value = path[7*i+5][stackNum], title=(num2char(65+i) + ": z pos")
		SetVariable $("radz"+num2str(i)), fsize=20,pos={155, 130+100*i},size={125,30}, value = path[7*i+6][stackNum], title=(num2char(65+i) + ": rad z")
	endfor
	//
	// Place buttons for moving the z-positions of all objects together
	SetDrawEnv fsize=20
	DrawText 20,(120+100*numObjects),"ALL OBJECTS: z pos"
	Button button0, fsize=20,pos={220,90+100*numObjects},size={20,20},proc=SDND2_AllZPlusButton,title="+"
	Button button1, fsize=20,pos={220,110+100*numObjects},size={20,20},proc=SDND2_AllZMinusButton,title="-"
	//
	// Place the other buttons, and a new feature that allows jumping to any stack
	Button button2, fstyle=1,pos={20,175+100*numObjects},size={55,30},proc=SDND2_IgnoreButton,title="IGNORE"
	Button button3, fstyle=1,pos={20,210+100*numObjects},size={55,30},proc=SDND2_BackButton,title="BACK"
	Button button4, fstyle=1,pos={85,175+100*numObjects},size={55,30},proc=SDND2_LinesButton,title="LINES"
	Button button5, fstyle=1,pos={85,210+100*numObjects},size={55,30},proc=SDND2_ViewButton,title="VIEW"
	Button button6, fstyle=1,pos={150,210+100*numObjects},size={130,30},proc=SDND2_NextButton,title="NEXT"
	SetVariable sv0,fsize=20,pos={150,175+100*numObjects},size={130,30},value=gJump,title="jump to"
	Button button7, fstyle=1,pos={75,135+100*numObjects},size={150,30},proc=SDND2_UpdateButton,title="UPDATE"
End

// These functions define the behavior of the 7 buttons in the Path Panel
// (all of them have the feature that they do nothing when in view mode (gView=1) )
Function SDND2_IgnoreButton(ctrlName) : ButtonControl
	string ctrlName
	//
	nvar gView, gIgnore, gUpdate
	if(gView==1)
		return 0
	endif
	//
	gIgnore = 1
	gUpdate = 0
	DoWindow/K $(winname(0,64))
End

Function SDND2_BackButton(ctrlName) : ButtonControl
	string ctrlName
	//
	nvar gView,gBack,gUpdate, gIndex
	if(gView==1 || gIndex==0)
		return 0
	endif
	//
	gBack = 1
	gUpdate = 0
	DoWindow/K $(winname(0,64))
End

Function SDND2_LinesButton(ctrlName) : ButtonControl
	string ctrlName
	//
	nvar gView,gLines,gUpdate
	if(gView==1)
		return 0
	endif
	//
	if(gLines==0)
		gLines=1
	else
		gLines=0
	endif
	gUpdate = 1
	DoWindow/K $(winname(0,64))
End

Function SDND2_ViewButton(ctrlName) : ButtonControl 
	string ctrlName
	//
	nvar gView,gUpdate,gIndex
	svar gPrefix
	//
	wave path = $(gPrefix+"path")
	wave/Z pathhold = $(gPrefix+"pathhold")
	if(!waveexists(pathhold))
		Make/O/N=(dimsize(path,0)) $(gPrefix+"pathhold") = 0
		wave pathhold = $(gPrefix+"pathhold")
	endif
	//
	if(gView==0)
		pathhold = path[p][gIndex]
		path[][gIndex] = 0
		gView=1
	else
		path[][gIndex] = pathhold[p]
		gView=0
	endif
	//
	gUpdate = 1
	DoWindow/K $(winname(0,64))
End

Function SDND2_NextButton(ctrlName) : ButtonControl
	string ctrlName
	//
	nvar gView,gUpdate
	if(gView==1)
		return 0
	endif
	//
	gUpdate = 0
	DoWindow/K $(winname(0,64))
End

Function SDND2_AllZPlusButton(ctrlName) : ButtonControl
	string ctrlName
	//
	svar gPrefix
	nvar gIndex, gView
	if(gView==1)
		return 0
	endif
	//
	wave path = $(gPrefix+"path")
	variable numObjects = dimsize(path,0)/7
	//
	variable i
	for(i=0; i<numObjects; i+=1)
		path[7*i+5][gIndex] += 1
	endfor
End

Function SDND2_AllZMinusButton(ctrlName) : ButtonControl
	string ctrlName
	//
	svar gPrefix
	nvar gIndex, gView
	if(gView==1)
		return 0
	endif
	//
	wave path = $(gPrefix+"path")
	variable numObjects = dimsize(path,0)/7
	//
	variable i
	for(i=0; i<numObjects; i+=1)
		path[7*i+5][gIndex] -= 1
	endfor
End

Function SDND2_UpdateButton(ctrlName) : ButtonControl
	string ctrlName
	//
	nvar gView, gUpdate
	if(gView==1)
		return 0
	endif
	//
	gUpdate = 1
	DoWindow/K $(winname(0,64))
End
//
Function SDND2_GetSignals(prefix)
	string prefix
	//
	variable runTime = datetime
	//
	wave paramWave = $(prefix+"params")
	wave path = $(prefix+"path")
	svar gFileName, gFileNameXY, gFileNameYZ, gFileNameXZ
	//
	string dataPathName = "path" + num2str(paramWave[0])
	string dataPathNameXY = dataPathName + "xy"
	string dataPathNameYZ =  dataPathName + "yz"
	string dataPathNameXZ = dataPathName + "xz"
	//
	variable numObjects = paramWave[1]
	variable numToAvg = paramWave[2]
	variable numRows = paramWave[3]
	variable numCols = paramWave[4]
	variable numPlanes = paramWave[5]	
	variable numStacks = paramWave[6]
	variable startStack = paramWave[10]
	variable endstack = paramWave[11]
	variable numStackDigits = strlen(num2str(numStacks))
	//
	// Make waves that will hold intensity signals, and that will hold the 3D mask
	Make/O/N=(numObjects,numStacks) $(prefix+"signal"+num2str(numToAvg)) = NaN
	wave signal = $(prefix+"signal"+num2str(numToAvg))
	Make/B/U/O/N=(numRows,numCols,numPlanes) $(prefix+"mask3D")
	wave mask3D = $(prefix+"mask3D")
	// 
	// Initialize projections and build the graph	
	ImageLoad/N=$"xFrame"/O/P=$dataPathNameYZ/Q/T=tiff (gFileNameYZ+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	ImageLoad/N=$"yFrame"/O/P=$dataPathNameXZ/Q/T=tiff (gFileNameXZ+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	ImageLoad/N=$"zFrame"/O/P=$dataPathNameXY/Q/T=tiff (gFileNameXY+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	variable graphSizeScale = 4
	//variable graphSizeScale = floor(400/numCols)
	NewImage zFrame; ShowInfo
	Appendimage/T=T2/L xFrame
	Appendimage/T/L=L2 yFrame
	ModifyGraph width=(graphSizeScale*numRows),height=(graphSizeScale*numCols)	
	ModifyGraph nticks=10,minor=1,freePos=0, fSize=8,btLen=3, mirror=0, lblPos=13, tlOffset=-2.00
	ModifyGraph tkLblRot(L2)=90
	SetAxis/A/R L2
	variable/C zProjFrac = SDND2_FindZProjFrac(prefix)
	ModifyGraph axisEnab(top)={0.00,real(zProjFrac)},axisEnab(left)={(1-imag(zProjFrac)),1.00}
	ModifyGraph axisEnab(L2)={0.00,(1-imag(zProjFrac)-0.05)},axisEnab(T2)={(real(zProjFrac)+0.05),1.00}
	string graphName = winName(0,1)
	//
	// Loop through the stacks and compute what we want
	variable i, j, k
	for(i=startStack; i<endStack; i+=1)
		// Load the 3D image (stack)
		SDND2_LoadOneStack(prefix,"stackWave",i)
		wave imageWave = $"stackWave"
		//
		// Update the projections on the graph
		ImageLoad/N=$"xFrame"/O/P=$dataPathNameYZ/Q/T=tiff (gFileNameYZ+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		ImageLoad/N=$"yFrame"/O/P=$dataPathNameXZ/Q/T=tiff (gFileNameXZ+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		ImageLoad/N=$"zFrame"/O/P=$dataPathNameXY/Q/T=tiff (gFileNameXY+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		DoUpdate
		//
		// If the path values aren't zero or NaN, do the computation
		if(path[0][i] == 0 || path[0][i]==NaN)
			signal[][i] = NaN
		else
			string command = "ImageGenerateROIMask/E=1/I=0 " + nameofwave(zFrame)
			// Loop through the objects, making a mask for each and analyzing the ROI pixels in it
			for(j=0; j<numObjects; j+=1)
				// Draw the oval for the object of interest
				SetDrawLayer ProgFront
				SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=top,ycoord=left, save	
				SetDrawEnv translate=path[7*j][i],path[7*j+1][i],rotate=path[7*j+4][i],rsabout
				DrawOval path[7*j][i]-path[7*j+2][i], path[7*j+1][i]-path[7*j+3][i], path[7*j][i]+path[7*j+2][i], path[7*j+1][i]+path[7*j+3][i]
				DoUpdate
				//
				// Generate the mask for the oval, then the full 3D mask by looping through the planes 
				execute command
				wave mask = $("M_ROIMask")
				mask3D = 1
				for(k=path[7*j+5][i]-path[7*j+6][i]; k<=path[7*j+5][i]+path[7*j+6][i]; k+=1)
					ImageTransform/P=(k)/D=mask setplane mask3D
				endfor	
				//
				// Use the 3D mask as the ROI to grab pixel values in the oval shape
				ImageTransform/R=mask3D roiTo1D imageWave
				wave W_roi_to_1d
				sort/R W_roi_to_1d W_roi_to_1d
				//
				// Average the brightest ___ pixels in the ROI to get the signal level for the object in this stack
				signal[j][i] = mean(w_roi_to_1d, 0, numToAvg)			
				SetDrawLayer/K ProgFront		
				KillWaves/Z mask			
			endfor			
		endif	
	endfor
	//
	KillWaves/Z xFrame, yFrame, zFrame
	KillWaves/Z mask3D
	KillWaves/Z imageWave
	//
	runTime = datetime - runTime
	print "GetSignals for " + prefix + "took " + num2str(runTime) + " seconds"
	//	
End
//
// Updated version of LoadOneStack, with altered references to the parameter wave
// and using the format of single-image TIFF files saved in Elements
Function SDND2_LoadOneStack(prefix, imageName, stackNum)
	string prefix
	string imageName
	variable stackNum
	//
	// Declare parameter wave and grab the relevant information from it
	wave/Z paramWave = $(prefix+"params")
	if(!waveexists($(prefix+"params")))
		print "tried to use LoadOneStack without a parameter wave"
		return 0
	endif
	variable numRows = paramWave[3]
	variable numCols = paramWave[4]
	variable numSlices = paramWave[5]
	variable numStacks = paramWave[6]
	variable numStackDigits = strlen(num2str(numStacks))
	variable numSliceDigits = strlen(num2str(numSlices))
	string dataPathStr = "path"+num2str(paramWave[0])
	svar gFileName
	variable i
	//
	//
	// Single-image TIFF series case...

		// Make the 3D wave that will hold all the slices we are about to load
		Make/W/U/O/N=(numRows,numCols,numSlices) $"temp" = 0
		wave imageWave = $"temp"		
		// Calculate how many digits are needed to represent the image number, and exactly which
		// files we will be loading
		variable firstIndex = stackNum * numSlices
		variable lastIndex = firstIndex + numSlices
		string nameStr = gFileName + "_t" + SDND2_num2strNdigit(stackNum+1,numStackDigits)
		string fileStr
		//
		// Loop through the individual TIFFs, loading each into Igor and then adding it to the 3D wave
		for(i=0; i<numSlices; i+=1)
			fileStr = nameStr + "z" + SDND2_num2strNdigit(i+1,numSliceDigits) + ".tif"
			//
			ImageLoad/N=$"slice"/O/P=$dataPathStr/Q/T=tiff fileStr
			wave sliceWave = $"slice"
			imageWave[][][i] = sliceWave[p][q]
			//
			KillWaves/Z sliceWave
		endfor	 
	//
	// Rename the loaded wave to be whatever the input to this function says it should be 
	Duplicate/O imageWave $(imageName)
	KillWaves imageWave
End
//
//
// Updated/modified version of GetPath4 from the SpinDisk4.ipf procedure
Function SDND2_GetPathAndSignals(prefix)
	string prefix
	//
	// Global variables
	variable/G gUpdate=1, gBack=0, gLines=0, gIgnore=0, gJump=0, gView=0, gIndex=0
	svar gFileName, gFileNameXY//, gFileNameYZ, gFileNameXZ
	variable defaultRad = 10
	//
	// Parameter wave
	wave paramWave = $(prefix+"params")
	string dataPathName = "path" + num2str(paramWave[0])
	string dataPathNameXY = dataPathName + "xy"	
	variable numObjects = paramWave[1]
	variable numToAvg = paramWave[2]
	variable numRows = paramWave[3]
	variable numCols = paramWave[4]
	variable numPlanes = paramWave[5]	
	variable numStacks = paramWave[6]
	variable startStack = paramWave[10]
	variable endStack = paramWave[11]
	variable numStackDigits = strlen(num2str(numStacks))
	variable numPlaneDigits = strlen(num2str(numPlanes))
	//
	// Make a new path with default values, if it doesn't already exist
	// Path has seven variables per object...
	// [0] = horizontal position
	// [1] = vertical position
	// [2] = horizontal radius
	// [3] = vertical radius
	// [4] = rotation angle
	// [5] = z-position
	// [6] = z "radius" 
	// Then the second object will take up indexes 7 through 13, and so on
	variable i
	if( !waveexists($(prefix+"path")) )
		Make/O/N=(7*numObjects,numStacks) $(prefix+"path") = 0
		wave path = $(prefix+"path")
		//
		for(i=0; i<numObjects; i+=1)
			path[7*i][startStack] = round(numRows/2)
			path[7*i+1][startStack] = round(numCols/2)
			path[7*i+2,7*i+3][startStack] = defaultRad
			path[7*i+4][startStack] = 0
			path[7*i+5][startStack] = round(numPlanes/2)
			path[7*i+6][startStack] = round(numPlanes/10)		
		endfor		
	else
		wave path = $(prefix+"path")
	endif
	//
	// Redimension the path wave if needed
	if(dimsize(path,0)/7 != numObjects)
		Redimension/N=(7*numObjects,-1) path
	endif
	//
	// Wave used to hold onto path values when stacks are skipped or jumped to
	Make/O/N=(7*numObjects,numStacks) $(prefix+"prevPath")
	wave prevPath = $(prefix+"prevPath")
	prevPath = path[p][startStack]
	//
	// Initialize projection images
	SDND2_LoadOneStack(prefix,"stackWave",0)
	wave imageWave = $"stackWave"
	// x-projection:
	ImageTransform/METH=1 xProjection imageWave
	wave M_xProjection
	MatrixOp/O M_xProjection = M_xProjection^t
	// y-projection:
	ImageTransform/METH=1 yProjection imageWave
	wave M_yProjection
	// z-projection:
	ImageLoad/N=$"zFrame"/O/P=$dataPathNameXY/Q/T=tiff (gFileNameXY+"_t"+SDND2_num2strNdigit(1,numStackDigits)+".tif")
	//
	// Build the graph to display / update as path is built
	// *** //
	variable graphSizeScale = 4
	//variable graphSizeScale = floor(400/numCols)
	NewImage zFrame; ShowInfo
	Appendimage/T=T2/L M_xProjection
	Appendimage/T/L=L2 M_yProjection
	ModifyGraph width=(graphSizeScale*numRows),height=(graphSizeScale*numCols)	
	ModifyGraph nticks=10,minor=1,freePos=0, fSize=8,btLen=3, mirror=0, lblPos=13, tlOffset=-2.00
	ModifyGraph tkLblRot(L2)=90
	SetAxis/A/R L2
	variable/C zProjFrac = SDND2_FindZProjFrac(prefix)
	ModifyGraph axisEnab(top)={0.00,real(zProjFrac)},axisEnab(left)={(1-imag(zProjFrac)),1.00}
	ModifyGraph axisEnab(L2)={0.00,(1-imag(zProjFrac)-0.05)},axisEnab(T2)={(real(zProjFrac)+0.05),1.00}
	string graphName = winName(0,1)
	//
	// Make waves that will hold intensity signals, and that will hold the 3D mask
	wave signal = $(prefix+"signal"+num2str(numToAvg))
	Make/B/U/O/N=(numRows,numCols,numPlanes) $(prefix+"mask3D")
	wave mask3D = $(prefix+"mask3D")
	//
	// Loop through the stacks (projection views) and let the user indicate the paths for the objects of interest, then
	// also compute the average intensity of the brightest ___ pixels
	variable j, k
	for(i=startStack; i<=endStack; i+=1)
		//
		// Update the projections on the graph
		SDND2_LoadOneStack(prefix,"stackWave",i)
		wave imageWave = $"stackWave"
		ImageTransform/METH=1 xProjection imageWave
		MatrixOp/O M_xProjection = M_xProjection^t
		ImageTransform/METH=1 yProjection imageWave
		ImageLoad/N=$"zFrame"/O/P=$dataPathNameXY/Q/T=tiff (gFileNameXY+"_t"+SDND2_num2strNdigit(i+1,numStackDigits)+".tif")
		DoUpdate
		//
		// For each object, set the path values for this stack equal to the previous stack 
		// if the current stack's path values are blank (as would be the case for a previously
		// partially-completed path, or empty path, or part of a path with a new object added, etc.)
		// (Note that "previous stack" will mean the previous stack that doesn't have 0's or NaN's)
		for(j=0; j<numObjects; j+=1)
			if(path[7*j][i]==0 || path[7*j][i]==NaN)
				for(k=0; k<7; k+=1)
					path[7*j+k][i] = prevPath[7*j+k]
				endfor
			endif		
		endfor
		//
		// Do-while loop through the panel that draws outlines of the ROIs onto the projection graph
		// (this continues looping as long as gUpdate = 1)
		gUpdate = 1
		gBack = 0
		gLines = 0
		gView = 0
		variable viewBefore
		//
		do
			// Place cursors on the z-projection
			for(j=0; j<numObjects; j+=1)
				cursor/C=(65535,65535,0)/H=(gLines)/I/P/S=1 $(num2char(65+j)) $"zFrame" path[7*j][i], path[7*j+1][i]
			endfor
			//
			// Draw ovals on the z-projection (light blue hollow shapes)
			SetDrawLayer ProgFront
			SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=top,ycoord=left,save
			for(j=0; j<numObjects; j+=1)
				SetDrawEnv translate=path[7*j][i],path[7*j+1][i],rotate=path[7*j+4][i],rsabout
				DrawOval path[7*j][i]-path[7*j+2][i], path[7*j+1][i]-path[7*j+3][i], path[7*j][i]+path[7*j+2][i], path[7*j+1][i]+path[7*j+3][i]			
			endfor
			//
			// Draw rectangles onto y-projection
			SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=top,ycoord=L2,save
			for(j=0; j<numObjects; j+=1)
				DrawRect path[7*j][i]-path[7*j+2][i], path[7*j+5][i]-path[7*j+6][i], path[7*j][i]+path[7*j+2][i], path[7*j+5][i]+path[7*j+6][i]
			endfor
			//
			// Draw rectangles onto x-projection
			SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=T2,ycoord=left,save
			for(j=0; j<numObjects; j+=1)
				DrawRect path[7*j+5][i]-path[7*j+6][i], path[7*j+1][i]-path[7*j+3][i], path[7*j+5][i]+path[7*j+6][i], path[7*j+1][i]+path[7*j+3][i]
			endfor
			//
			// Display the panel, which allows the user to move cursors and change ROI size in all three dimensions
			// (Also, note the value of the "view" global variable before the panel is displayed)
			viewBefore = gView
			SDND2_PathPanel(prefix,i)
			PauseForUser $(winname(0,64)),$graphName
			//
			// Save the x and y cursor positions if the UPDATE button was pressed
			// (the other settings are saved to the path wave in real time)
			// (And don't do this if the user has just turned the VIEW back on)
			if(gUpdate == 1 && (gView==0 && viewBefore!=1) ) 
				for(j=0; j<numObjects; j+=1)
					path[7*j][i] = pcsr($(num2char(65+j)))
					path[7*j+1][i] = qcsr($(num2char(65+j))) 
				endfor
			endif
			//
			// Erase the ovals and rectangles before the next time through the loop
			SetDrawLayer/K ProgFront	
			//
		while(gUpdate)
		//
		// If the user pressed the "ignore" button, then set the path values for this stack to NaN
		if(gIgnore==1)
			path[][i] = NaN				
			gIgnore=0
		else
			prevPath = path[p][i]
		endif
		//
		// GET SIGNALS PART 
		// If the path values aren't zero or NaN, do the computation
		if(path[0][i] == 0 || path[0][i]==NaN)
			signal[][i] = NaN
		else
			string command = "ImageGenerateROIMask/E=1/I=0 " + nameofwave(zFrame)
			// Loop through the objects, making a mask for each and analyzing the ROI pixels in it
			for(j=0; j<numObjects; j+=1)
				// Draw the oval for the object of interest
				SetDrawLayer ProgFront
				SetDrawEnv linefgc=(0,65535,65535),fillpat=0,xcoord=top,ycoord=left, save	
				SetDrawEnv translate=path[7*j][i],path[7*j+1][i],rotate=path[7*j+4][i],rsabout
				DrawOval path[7*j][i]-path[7*j+2][i], path[7*j+1][i]-path[7*j+3][i], path[7*j][i]+path[7*j+2][i], path[7*j+1][i]+path[7*j+3][i]
				DoUpdate
				//
				// Generate the mask for the oval, then the full 3D mask by looping through the planes 
				execute command
				wave mask = $("M_ROIMask")
				mask3D = 1
				for(k=path[7*j+5][i]-path[7*j+6][i]; k<=path[7*j+5][i]+path[7*j+6][i]; k+=1)
					ImageTransform/P=(k)/D=mask setplane mask3D
				endfor	
				//
				// Use the 3D mask as the ROI to grab pixel values in the oval shape
				ImageTransform/R=mask3D roiTo1D imageWave
				wave W_roi_to_1d
				sort/R W_roi_to_1d W_roi_to_1d
				//
				// Average the brightest ___ pixels in the ROI to get the signal level for the object in this stack
				signal[j][i] = mean(w_roi_to_1d, 0, numToAvg)			
				SetDrawLayer/K ProgFront		
				KillWaves/Z mask			 
			endfor			
		endif	
		//
		// If the user pressed the "back" button, then set gJump to the previous stack
		if(gBack==1)
			gJump = i - 1
		endif
		//
		// If the user selected a different stack to go to (instead of the next stack), adjust the loop index (i) accordingly
		// (without changing anything, gJump is just i+1 and there is no change here)
		i += gJump - i - 1
		//
	endfor
	//
	// Get rid of waves made inside this function that we don't keep using
	DoWindow/K $graphName
	KillWaves/Z zFrame
	KillWaves/Z prevPath
	KillWaves/Z mask3D
	KillWaves/Z imageWave
	KillWaves/Z W_roi_to_1d
	// 
End // of SDND2_GetPathAndSignals
//
// Simple "OK" panel...
Function SDND2_OKpanel() : Panel
	//
	NewPanel/W=(1200,300,1500,400)
	//
	Button button0,fsize=20,pos={50,30},size={200,40},proc=SDND2_OKbuttonProc,title="OK"
	//
End
//
// Copy of the function that assigns colors (used for traces on graphs)
Function SDND2_TenColors(colorStr, colorNum)
	string colorStr
	variable colorNum
	//
	// response to invalid color string
	if(cmpstr(colorStr,"red")!=0 && cmpstr(colorStr,"green")!=0 && cmpstr(colorStr,"blue")!=0)
		return 0
	endif
	//
	// choose the r, g, or b component of the color number
	switch(colorNum)
		case 0: // (BLUE)
			if(cmpstr(colorStr,"red")==0)
				return 0
			elseif(cmpstr(colorStr,"green")==0)
				return 15872
			else
				return 65280
			endif
			break
		case 1: // (GREEN)
			if(cmpstr(colorStr,"red")==0)
				return 0
			elseif(cmpstr(colorStr,"green")==0)
				return 52224
			else
				return 0
			endif
			break
		case 2: // (BLACK)
			if(cmpstr(colorStr,"red")==0)
				return 0
			elseif(cmpstr(colorStr,"green")==0)
				return 0
			else
				return 0
			endif
			break
		case 3: // (RED)
			if(cmpstr(colorStr,"red")==0)
				return 65280
			elseif(cmpstr(colorStr,"green")==0)
				return 0
			else
				return 0
			endif
			break
		case 4: // (PURPLE)
			if(cmpstr(colorStr,"red")==0)
				return 36864
			elseif(cmpstr(colorStr,"green")==0)
				return 14592
			else
				return 58880
			endif
			break
		case 5: // (ORANGE)
			if(cmpstr(colorStr,"red")==0)
				return 65280
			elseif(cmpstr(colorStr,"green")==0)
				return 32512
			else
				return 16384
			endif
			break
		case 6: // (YELLOW)
			if(cmpstr(colorStr,"red")==0)
				return 52224
			elseif(cmpstr(colorStr,"green")==0)
				return 52224
			else
				return 0
			endif
			break
		case 7: // (LIGHT BLUE)
			if(cmpstr(colorStr,"red")==0)
				return 16384
			elseif(cmpstr(colorStr,"green")==0)
				return 48896
			else
				return 65280
			endif
			break
		case 8: // (MAGENTA)
			if(cmpstr(colorStr,"red")==0)
				return 65280
			elseif(cmpstr(colorStr,"green")==0)
				return 16384
			else
				return 55552
			endif
			break
		case 9: // (GRAY)
			if(cmpstr(colorStr,"red")==0)
				return 47872
			elseif(cmpstr(colorStr,"green")==0)
				return 47872
			else
				return 47872
			endif
			break
		default: // (just use black if the color number is not 0-9)
			return 0
	endswitch		
End // of TenColors function
//
//
// 
// (12/17/2012)
// Writing function(s) to compute DF/F in an automated fashion, it was getting pretty annoying
// doing it manually
//
Function SDND2_deltaFoverF(prefix, Fwave, BGlevel, delaytime)

	string prefix
	wave Fwave
	variable BGlevel
	variable delaytime // in seconds, and assumes the temperature file starts recording after the movie does)
	
	// waves that will be our sources
	// (for the paramter wave, element [2] holds the number of brightest pixels averaged)
	wave paramWave = $(prefix+"params")
	wave temperatureWave = $(prefix+"temp0")
	wave signalWave = $(prefix+"signal"+num2str(paramWave[2]) )
	wave timeWave = $(prefix+"time0")
	
	// The F-wave needs the right format or quit the function
	if( dimsize(Fwave, 1) != dimsize(signalWave,0)+1 )
		print "F-wave is the wrong size"
		return 0
	endif
	
	// new temperature waves (so we don't mess with the originals)
	make/O/N=(dimsize(temperatureWave,0)) $(prefix+"temptime") = temperatureWave[p][0]
	make/O/N=(dimsize(temperatureWave,0)) $(prefix+"tempSV") = temperatureWave[p][1]
	make/O/N=(dimsize(temperatureWave,0)) $(prefix+"tempPV") = temperatureWave[p][2]
	wave temptime = $(prefix+"temptime")
	wave tempSV = $(prefix+"tempSV")
	wave tempPV = $(prefix+"tempPV")
	
	// waves that will hold the raw values
	variable i
	for(i=0; i<dimsize(signalWave,0); i+=1)
		make/O/N=(dimsize(signalWave,1)) $(nameofwave(signalWave)+"raw"+num2str(i)) = signalWave[i][p]
	endfor
	make/O/N=(dimsize(signalWave,1)) $(nameofwave(signalWave)+"rawtime") = timeWave[p]
	wave timeToGraph = $(nameofwave(signalWave)+"rawtime")
	
	// waves that will hold the F-levels:
	for(i=0; i<dimsize(signalWave,0); i+=1)
		make/O Fx = Fwave[p][0]
		make/O Fy = Fwave[p][i+1]
		make/O/N=(dimsize(signalWave,1)) $(nameofwave(signalWave)+"F"+num2str(i)) = interp(timeWave, Fx,Fy)
	endfor
	
	// waves that will hold the DF/F values:
	for(i=0; i<dimsize(signalWave,0); i+=1)
		wave rawWave =  $(nameofwave(signalWave)+"raw"+num2str(i))
		wave Fwave = $(nameofwave(signalWave)+"F"+num2str(i))
		make/O/N=(dimsize(signalWave,1)) $(nameofwave(signalWave)+"DFF"+num2str(i)) = (rawWave-Fwave)/(Fwave-BGlevel)
	endfor
	
	// re-zero the time and account for the delay-shift-thing
	variable x1, x2
	for(i=1; i<dimsize(tempSV,0); i+=1)
		if(abs(tempSV[i] - tempSV[i-1]) > 0)
			x1 = i
			break
		endif	
	endfor
	//x2 = x1 + 20 // this is currently too specific for T=3 s periods
	//FindPeak/Q/P/R=(x1,x2) tempSV
	//variable shift = temptime[V_PeakLoc] - 0.75 // also specific to T = 3 s
	variable shift = temptime[x1]
	temptime -= shift
	timeToGraph -= shift
	timeToGraph -= delaytime
	
	// plot it?
	Display tempSV, tempPV vs temptime
	ModifyGraph axisEnab(left) = {0,0.4}, nticks(bottom)=20
	Label bottom "Time"
	Label left "Temp. (C)"
	for(i=0; i<dimsize(signalWave,0); i+=1)
		wave waveToGraph = $(nameofwave(signalWave)+"DFF"+num2str(i))
		AppendToGraph/C=(SDND2_TenColors("red",i),SDND2_TenColors("green",i),SDND2_TenColors("blue",i))/L=left2 waveToGraph vs timeToGraph
	endfor
	Label left2 "delta F / F"
	ModifyGraph lsize=2,lblPos(left2)=49,axisEnab(left2)={0.5,1},freePos(left2)=0
	Legend


End // of SDND2_deltaFoverF function
//
//
// (09/14/2014)
// Function for displaying the temperature in a separate graph, and printing the normally-chosen
// "zero time" that is used later for DF/F and averaging of data graphs
Function SDND2_OtherTempGraph(tempWave)

	// input temperature wave, which is generally four columns (time, SV, PV, PV2), and
	// is originally loaded in a text file that is generated in LabVIEW during experiments
	wave tempWave
	
	// determine the index where the SV first starts to change
	variable numTemps = dimsize(tempWave,0)
	variable prevSV, currentSV
	variable zeroIndex = -1
	variable i
	for(i=1; i<numTemps; i+=1)
		prevSV = tempWave[i-1][1]
		currentSV = tempWave[i][1]
		if(currentSV != prevSV)
			zeroIndex = i-1
			break
		endif	
	endfor
	
	// string for reporting the time where SV changes
	string SVchange = "SV change @ " + num2str(tempWave[zeroIndex][0]) + " s"
	
	// display the temperatures (SV=yellow, PV=red)
	display tempWave[][1] vs tempWave[][0]
	modifygraph rgb=(52224,52224,0)
	showinfo
	//Cursor/P A, tempWave, zeroIndex
	appendtograph tempWave[][2] vs tempWave[][0]	
	Legend/K/N=text0
	TextBox/C/N=text0 SVchange
	
	

End // of SDND2_OtherTempGraph




// Functions written 06/29/2012 for just loading raw images (Elements naming format) and storing
// projection waves
Function GetProjectionsSimple(filePrefix)
	string filePrefix
	//
	// Keep track of how long it takes the function to run
	variable runTime = datetime

	variable numStacks = 92
	variable numPlanes = 46
	variable numRows = 512
	variable numCols = 512
	
	variable xyStep = 0.222
	variable zStep = 0.6
	variable horzZstep =(0.95)*(xyStep*numRows)/(xyStep*numRows + zStep*numPlanes)
	variable vertZstep = (0.95)*(xyStep*numCols)/(xyStep*numCols + zStep*numPlanes)	
	variable/C zProjFrac = cmplx(horzZstep,vertZstep)
	//
	// Make the waves to hold the projections
	Make/W/U/O/N=(numPlanes,numRows,numStacks) $(filePrefix+"xProj") = 0
	Make/W/U/O/N=(numRows,numPlanes,numStacks) $(filePrefix+"yProj") = 0
	Make/W/U/O/N=(numRows,numCols,numStacks) $(filePrefix+"zProj") = 0
	wave xProj = $(filePrefix+"xProj")
	wave yProj = $(filePrefix+"yProj")
	wave zProj = $(filePrefix+"zProj")
	//
	// Data path:
	string pathStr = filePrefix+"path"
	NewPath /O/Q $pathStr	
	//
	// Load the first stack, just for setting up the graph and declaring the projection waves
	LoadOneStackNoInfo(filePrefix, "stackWave", 0, numPlanes, numStacks, pathStr)
	wave imageWave = $"stackWave"
	ImageTransform/METH=1 xProjection imageWave
	wave M_xProjection
	MatrixOp/O M_xProjection = M_xProjection^t
	ImageTransform/METH=1 yProjection imageWave
	wave M_yProjection
	ImageTransform/METH=1 zProjection imageWave
	wave M_zProjection
	//
	// Display the the projections and format the graph in a pleasing fashion
	// (this is taken straight from the AskForROI function)
	NewImage M_zProjection; ShowInfo
	string graphName = winname(0,1)
	AppendImage/T=T2/L M_xProjection
	AppendImage/T/L=L2 M_yProjection
	ModifyGraph nticks=10,minor=1,freePos=0,fSize=8,btLen=3,mirror=0,lblPos=13,tlOffset=-2.00
	ModifyGraph tkLblRot(L2)=90
	SetAxis/A/R L2
	//variable/C zProjFrac = FindZProjFraction(prefix)
	ModifyGraph axisEnab(top)={0.00,real(zProjFrac)},axisEnab(left)={(1-imag(zProjFrac)),1.00}
	ModifyGraph axisEnab(L2)={0.00,(1-imag(zProjFrac)-0.05)},axisEnab(T2)={(real(zProjFrac)+0.05),1.00}
	// 
	// Go through the images, loading each stack and projecting it
	variable i
	for(i=0; i<numStacks; i+=1)
		// Load the image
		LoadOneStackNoInfo(filePrefix, "stackWave", i, numPlanes, numStacks, pathStr)
		wave imageWave = $"stackWave"
		//
		// Projections
		ImageTransform/METH=1 zProjection imageWave
		ImageTransform/METH=1 yProjection imageWave
		ImageTransform/METH=1 xProjection imageWave
		MatrixOp/O M_xProjection = M_xProjection^t
		//
		// Update the graph and delete the 3D wave
		DoUpdate
		KillWaves imageWave
		//
		// Save the projection images into the three projection waves
		xProj[][][i] = M_xProjection[p][q]
		yProj[][][i] = M_yProjection[p][q]
		zProj[][][i] = M_zProjection[p][q]
		//
	endfor
	//
	// Delete the graph
	DoWindow/K $graphName
	KillWaves/Z M_xProjection, M_yProjection, M_zProjection
	//
	// Report the elapsed time
	runTime = datetime - runTime
	print "GetProjectionsSimple for " + filePrefix + "took " +num2str(runTime) + " seconds"
End
//
Function LoadOneStackNoInfo(filePrefix, imageName, stackNum, numSlices, numStacks, dataPathStr)
	string filePrefix
	string imageName
	variable stackNum
	
	variable numSlices 
	variable numStacks
	
	string dataPathStr
	
	variable numRows = 512
	variable numCols = 512
	variable numStackDigits = strlen(num2str(numStacks))
	variable numSliceDigits = strlen(num2str(numSlices))
	
	variable i

		// Make the 3D wave that will hold all the slices we are about to load
		Make/W/U/O/N=(numRows,numCols,numSlices) $"temp" = 0
		wave imageWave = $"temp"		
		// Calculate how many digits are needed to represent the image number, and exactly which
		// files we will be loading
		variable firstIndex = stackNum * numSlices
		variable lastIndex = firstIndex + numSlices
		string nameStr = filePrefix + "_t" + SDND2_num2strNdigit(stackNum+1,numStackDigits)
		string fileStr
		//
		// Loop through the individual TIFFs, loading each into Igor and then adding it to the 3D wave
		for(i=0; i<numSlices; i+=1)
			fileStr = nameStr + "z" + SDND2_num2strNdigit(i+1,numSliceDigits) + ".tif"
			//
			ImageLoad/N=$"slice"/O/P=$dataPathStr/Q/T=tiff fileStr
			wave sliceWave = $"slice"
			imageWave[][][i] = sliceWave[p][q]
			//
			KillWaves/Z sliceWave
		endfor	 
	//
	// Rename the loaded wave to be whatever the input to this function says it should be 
	Duplicate/O imageWave $(imageName)
	KillWaves imageWave
End
// 
Function PlayProjSimplePanel(prefix) : Panel
	string prefix
	//
	// Global Variables...
	nvar/Z saveMovieYN
	if(!nvar_exists(saveMovieYN))
		variable/G saveMovieYN = 1
	endif
	//
	nvar/Z currentGraphYN
	if(!nvar_exists(currentGraphYN))
		variable/G currentGraphYN = 0
	endif
	//
	svar/Z colorName
	if(!svar_exists(colorName))
		string/G colorName = "Grays"
	endif
	variable colorMode = 1 + WhichListItem(colorName,CTabList())
	//
	nvar/Z LUTlow
	if(!nvar_exists(LUTlow))
		variable/G LUTlow = 0
	endif
	//
	nvar/Z LUThigh
	if(!nvar_exists(LUThigh))
		variable/G LUThigh = 2^16-1
	endif
	//	
	// Build the panel itself
	NewPanel/W=(1200,100,1550,420)
	// Title:
	string panelTitle = "PLAY MOVIE " + prefix
	TitleBox title0, fixedsize=1,fstyle=1,fsize=20,pos={15,8},size={320,32},title=panelTitle
	// Check boxes:
	CheckBox check0,fsize=20,pos={20,50},size={100,30},variable=saveMovieYN,title="save movie as an AVI?"
	CheckBox check1,fsize=20,pos={20,140},size={100,30},variable=currentGraphYN,title="play on existing graph?"
	// LUT values:
	SetVariable sv0, fsize=20,pos={20,175},size={160,30},value=LUTlow,limits={0,2^(16)-1,1},title="LUT low"
	SetVariable sv1, fsize=20,pos={20,210},size={160,30},value=LUThigh,limits={0,2^(16)-1,1},title="LUT high"
	// Color scheme:	
	PopupMenu pop0, fsize=20,pos={20,245},bodywidth=200,size={200,30},mode=colorMode,value="*COLORTABLEPOP*",proc=SDND2_ColorScheme,title="color"
	// OK button:
	Button button0,fsize=20,pos={115,275},size={120,30},proc=SDND2_okbuttonPROC,title="OK"
End
//
Function SDND2_ColorScheme(ctrlName,popNum,popStr) : PopupMenuControl
	string ctrlName
	variable popNum
	string popStr
	//
	svar colorName
	colorName = popStr	
End
//
Function PlayProjectionsSimple(prefix)
	string prefix
	//
	variable numStacks = 92
	variable numPlanes = 46
	variable numRows = 512
	variable numCols = 512
	//
	variable xyStep = 0.222
	variable zStep = 0.6
	variable horzZstep =(0.95)*(xyStep*numRows)/(xyStep*numRows + zStep*numPlanes)
	variable vertZstep = (0.95)*(xyStep*numCols)/(xyStep*numCols + zStep*numPlanes)	
	variable/C zProjFrac = cmplx(horzZstep,vertZstep)
	//
	// Panel:
	PlayProjSimplePanel(prefix)
	PauseForUser $(winname(0,64))
	nvar saveMovieYN,currentGraphYN, LUTlow, LUThigh
	svar colorName
	//
	// Declare the input projection waves, and end the function if they aren't there
	wave xProj = $(prefix+"xProj")
	wave yProj = $(prefix+"yProj")
	wave zProj = $(prefix+"zProj")
	if(!waveexists(xProj))
		return 0
	endif
	//
	// Declare or make the waves that will be displayed in the movie
	wave/Z xFrame,yFrame,zFrame
	if(!waveexists(xFrame) || currentGraphYN==0)
		Make/O/N=(dimsize(xProj,0),dimsize(xProj,1)) xFrame
		Make/O/N=(dimsize(yProj,0),dimsize(yProj,1)) yFrame
		Make/O/N=(dimsize(zProj,0),dimsize(zProj,1)) zFrame
		currentGraphYN=0	
	endif
	//
	// Initialize the first frame
	xFrame = xProj[p][q][0]
	yFrame = yProj[p][q][0]
	zFrame = zProj[p][q][0]
	//
	// Set up the graph, if not using an existing one
	if(currentGraphYN==0)
		NewImage zFrame; ShowInfo
		Appendimage/T=T2/L xFrame
		Appendimage/T/L=L2 yFrame
		//
		ModifyGraph width=(1.0*dimsize(zFrame,0)),height=(1.0*dimsize(zFrame,1))	
		ModifyGraph nticks=10,minor=1,freePos=0, fSize=8,btLen=3, mirror=0, lblPos=13, tlOffset=-2.00
		ModifyGraph tkLblRot(L2)=90
		SetAxis/A/R L2			
		ModifyGraph axisEnab(top)={0.00,real(zProjFrac)},axisEnab(left)={(1-imag(zProjFrac)),1.00}
		ModifyGraph axisEnab(L2)={0.00,(1-imag(zProjFrac)-0.05)},axisEnab(T2)={(real(zProjFrac)+0.05),1.00}
	endif		
	string graphName = winName(0,1)
	//
	// Scale the image brightness (LUT) and use the color scheme chosen by the user
	ModifyImage xFrame ctab={LUTlow,LUThigh,$colorName,0}
	ModifyImage yFrame ctab={LUTlow,LUThigh,$colorName,0}
	ModifyImage zFrame ctab={LUTlow,LUThigh,$colorName,0}
	DoUpdate
	//
	// If saving a movie, set it up now
	if(saveMovieYN==1)
		NewMovie/A/F=5/I
	endif
	//
	// Loop through the stacks and update the image, effectively playing the movie
	variable i
	for(i=0; i<numStacks; i+=1)
		// Update the images
		xFrame = xProj[p][q][i]
		yFrame = yProj[p][q][i]
		zFrame = zProj[p][q][i]
		// Update the graph
		DoUpdate
		// Save the current image to the movie, if applicable
		if(saveMovieYN==1)
			AddMovieFrame
		endif	
	endfor
	//
	// Close the movie, if applicable (i.e., finish making the AVI file)
	if(saveMovieYN==1)
		CloseMovie
	endif
End