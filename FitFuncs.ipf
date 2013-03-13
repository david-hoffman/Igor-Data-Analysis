#pragma rtGlobals=3		// Use modern global access method.

//Various fitting functions for analyzing FSRS and ISRS data, predominantly various
//combinations of exponentials convoluted with an instrument response function which 
//is modeled as being a gaussian

//******************************************************************
//DISP Gaussian Function
//******************************************************************

Function fDispGauss(pw, yw, xw)
//An all at once style function that simulates dispersive gaussians
//This is useful for fitting RINE ((1) McCamant, D. W.; Kukura, P.; Mathies, R. A.
//The Journal of Physical Chemistry B 2005, 109, 10449Ð10457.) 
//data when the Raman pump is short enough to truncate the FID of the vibrational coherence.
	WAVE pw, yw, xw
	//this function has no offset!
	//Error checking here
	variable length = numpnts(pw) //how many parameters?
	variable i=0
	yw=0//set the target wave to zero, we'll be adding to it later
	//Make a temp wave to hold the intermediate results
	Make/D/O/N=(numpnts(yw)) myTempWave = 0
	for(i=0;i<length;i+=4)
		// cycle through the parameters building up the spectrum
		// from the various peaks.
		Duplicate/O/R=[i,i+3] pw myCoefs
		dispGauss(myCoefs,myTempWave,xw)
		yw+=myTempWave
	EndFor
End

STATIC Function dispGauss(pw, yw, xw)
//This is the function that actually creates a single dispersed gaussian
	//this function has no offset!
	//my passed waves
	WAVE pw, yw, xw
	
	//Make my signal wave and my lorentzian wave
	Make/D/O/N=20000 lor, signal
	//Set my scales, to make things easier
	SetScale/P x 0,0.1,"", lor, signal
	
	//Make my signal which is essentially a gaussian plus its derivative
	Signal = (pw[0]+pw[1]*(x-pw[2]))*exp(-((x-pw[2])/1.3)^2)/1.3
	
	//Make my lorentzian, this determines the final "width" of the line shape
	lor = 1/(((x-1000)/pw[3])^2+1)
	
	//Convolve them
	Duplicate/O lor,lor_conv
	Convolve/A signal, lor_conv
	
	//Make the result wave by sampling at the correct x values.
	yw = lor_conv(xw[p])
End

Function dispLor(w,x)
//A dispersive lorentzian
//w[i] = amplitude
//w[i+1] = center frequency
//w[i+2] =width

	Wave w
	Variable x
	
	Variable i, Val=0
	
	For(i=0;i<numpnts(w);i+=3)
		val += w[i]*(x-w[i+1])*w[i+2]/((x-w[i+1])^2+w[i+2]^2)
	EndFor
	Return val
End

//******************************************************************

//******************************************************************
//Various convoluted exponentials
//******************************************************************

Function convWeird(w,t) 
	Wave w; Variable t
	
		//|  w contains the parameters needed for the convolution and exponential functions.
		//| Equation = A*(1 + Erf((t-t0)/w) +Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2)) (1 + Erf((t-t0)/w - w/(2*t1))))
		//|  w[0]=w the gaussian width(measured by cross-correlation)
		//|  w[1]=t0 t-zero
		//|  w[2]=y0 offset
		//|  w[3]=A0 scaling coefficient
		//|  w[4]=t1, the decay (1/e) time of the first exponential
		//|  w[5]=A1 scaling coefficient
		//|  f(t) = w[5] + (w[3] - w[5])exp(-t/w[6]) - w[3] exp(-t/w[4])
		Variable val= -w[3]*Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[4])))
		val+=(w[3]-w[5])*Exp((-4*(t-w[1])*w[6] + w[0]^2)/(4*w[6]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[6])))
		val+=w[5]*(1 + Erf((t-w[1])/w[0]))
		val/=2
		val+=w[2]
		return val
End

Function IntermediateDecay(w,t) : FitFunc
	WAVE w
	Variable t
	
	//|This function was written by David Hoffman 10/2009
	//| This function models the decay of an intermediate convoluted with an exponential.
	//| For example [A2](t) in [A1]->[A2]->[A3] 
	//| Equation: A*((Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2))*(1 + Erf((t-t0)/w - w/(2*t1)]))/(t1 - t2) + (Exp((-4*(t-t0)*t2 + w^2)/(4*t2^2)) (1 + Erf((t-t0)/w - w/(2*t2)]))/(-t1 + t2))
	//| Coefficients 5
	//| w[0] = w gaussian width of the cross correlation
	//| w[1] = t0 t-zero
	//| w[2] = A Amplitude factor
	//| w[3] = t1 inverse rate for A1->A2
	//| w[4] = t2 inverse rate for A2->A3
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = a
	//CurveFitDialog/ w[3] = t1
	//CurveFitDialog/ w[4] = t2

	variable val=((Exp((-4*(t-w[1])*w[3] + w[0]^2)/(4*w[3]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[3]))))/(w[3] - w[4]) + (Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0]- w[0]/(2*w[4]))))/(-w[3] + w[4]))
	val*=w[2]*w[4]
	val/=2
	return val
End

Function ID1k(w,t) : FitFunc
	WAVE w
	Variable t
	
	//|This function was written by David Hoffman 10/2012
	//| This function models the decay of an intermediate convoluted with an exponential plus one more exponential
	//| For example [A2](t) in [A1]->[A2]->[A3] 
	//| In this case there is only 1 time constant
	//| Equation:  ( (2 t1 w+E^((2 t t1-w^2)^2/(4 t1^2 w^2)) (2 t t1-w^2) (1+Erf[t/w-w/(2 t1)])))/(E^(t^2/w^2))
	//| Coefficients 4
	//| w[0] = w gaussian width of the cross correlation
	//| w[1] = t0 t-zero
	//| w[2] = A Amplitude factor
	//| w[3] = t1 inverse rate for A1->A2 and A2->A3
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = a
	//CurveFitDialog/ w[3] = t1
	
	t=t-w[1]
	
	Variable t1 = w[3]
	Variable width = w[0]
	Variable A = w[2]
	
	variable val=-Exp((-4*t*t1+width^2)/(4*t1^2))
	val*=Sqrt(Pi)
	val*=(-2*t*t1+width^2)
	val*=(1+Erf(t/width-width/(2*t1))*Sign(1-(2*t*t1)/width^2)^2)
	val+=2*t1*width*Exp(-(t/width)^2)
	val*=A*width/(2*t1)^2
	return val
End

Function IntermediateDecayPlus1(w,t) : FitFunc
	WAVE w
	Variable t
	
	//|This function was written by David Hoffman 5/2011
	//| This function models the decay of an intermediate plus an exponential decay convoluted with an exponential
	//| For example [A2](t) in [A1]->[A2]->[A3] 
	//| Equation: A*((Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2))*(1 + Erf((t-t0)/w - w/(2*t1)]))/(t1 - t2) + (Exp((-4*(t-t0)*t2 + w^2)/(4*t2^2)) (1 + Erf((t-t0)/w - w/(2*t2)]))/(-t1 + t2))
	//| Coefficients 7
	//| w[0] = w gaussian width of the cross correlation
	//| w[1] = t0 t-zero
	//| w[2]=y0 offset
	//| w[3] = A Amplitude factor
	//| w[4] = t1 inverse rate for A1->A2
	//| w[5] = t2 inverse rate for A2->A3
	//| w[6] = freq of single exp
	//| w[7] = Amplitude of single exp
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = A2
	//CurveFitDialog/ w[4] = t2
	//CurveFitDialog/ w[5] = t3
	//CurveFitDialog/ w[6] = t1
	//CurveFitDialog/ w[7] = A1

	variable val=w[3]*((Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[4]))))*sign(w[4] - w[5]) + (Exp((-4*(t-w[1])*w[5] + w[0]^2)/(4*w[5]^2))*(1 + Erf((t-w[1])/w[0]- w[0]/(2*w[5]))))*sign(-w[4] + w[5]))
	val+=w[7]*Exp((-4*(t-w[1])*w[6] + w[0]^2)/(4*w[6]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[6])))
	val/=2
	val+=w[2]
	return val
End

Function conv_ExpRise(w,t) : FitFunc
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation=A(1 + Erf[(t-t0)/w] + E^((-4 (t-t0) t1 + w^2)/(4 t1^2)) (-1 - Erf[(t-t0)/w - w/(2 t1)]))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a
	//CurveFitDialog/ w[4] = t
	
	Variable val=-Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[4])))
	val+=1 + Erf((t-w[1])/w[0])
	val*=w[3]
	val/=2
	val+=w[2]
	return val
End

Function Kerr(pw,yw,tw)
//Models the kerr effect, which includes an instantaneous electronic response
//and a slower, exponential, vibrational response
//This function is written as an all at once function
	Wave pw,yw,tw
	//Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation=A*Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2))*(1 + Erf((t-t0)/w - w/(2*t1)))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	//CurveFitDialog/ w[5] = A
	
	yw=Exp((-4*(tw-pw[1])*pw[4]+pw[0]^2)/(4*pw[4]^2))*(1 + Erf((tw-pw[1])/pw[0]-pw[0]/(2*pw[4]), 1e-16))
	yw*=pw[3]
	yw/=2
	yw+= Exp(-((tw-pw[1])/pw[0])^2)*pw[5]
	yw+=pw[2]
	return 0
End

Function KerrFunc(w,t) : FitFunc
//Same as above but written in the traditional style so that if can be used
//in the curve fitting dialog
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation=A*Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2))*(1 + Erf((t-t0)/w - w/(2*t1)))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	//CurveFitDialog/ w[5] = A
	
	variable yw=Exp((-4*(t-w[1])*w[4]+w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0]-w[0]/(2*w[4]), 1e-16))
	yw*=w[3]
	yw/=2
	yw+= Exp(-((t-w[1])/w[0])^2)*w[5]
	yw+=w[2]
	return yw
End

Function ErrorFunction(w,t) : FitFunc
//An error function that is trunctated at t0
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0+A*erf((x-x0)/w)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = A
	//CurveFitDialog/ w[2] = t0
	//CurveFitDialog/ w[3] = w
	

	return w[0]+(t>=w[2])*w[1]*erf((t-w[2])/w[3])
End

Function popdepletionTEST(pw, yw, xw)
	//my passed waves
	WAVE pw, yw, xw
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation=A*Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2))*(1 + Erf((t-t0)/w - w/(2*t1)))
	//|  pw[0]=w the gaussian width(measured by cross-correlation)
	//|  pw[1]=t0 t-zero
	//|  pw[2]=y0 offset
	//|  pw[3]=A0 scaling coefficient
	//|  pw[4]=t1, the decay (1/e) time of the first exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ pw[0] = W
	//CurveFitDialog/ pw[1] = t0
	//CurveFitDialog/ pw[2] = y0
	//CurveFitDialog/ pw[3] = a1
	//CurveFitDialog/ pw[4] = t1
	
	//We will use a tenth of the pw[0] as our resolution
	Variable dT = max(5,abs(round(pw[0]/10)))//fs
	
	//This ensures that there is enough time for the gaussian to decay completely
	Variable maxT = wavemax(xw)*1.1, minT = wavemin(xw)*1.1
	
	Variable numberOfPoints = round((maxT-minT)/dT)
	
	//Make my signal wave and my gaussian wave
	Make/D/O/N=(numberOfPoints) G, signal
	//Set my scales, to make things easier
	SetScale/P x minT, dT,"", G, signal//Needs to be symmetric about zero for later
	
	//Make my signal which is essentially a gaussian plus its derivative
	G = exp(-((x-(minT+MaxT)/2)/pw[0])^2)
	//Normalize gaussian
	G/=sum(G)
	
	//Signal
	Signal = (x>=pw[1])*pw[3]*erf((x-pw[1])/pw[4])
	
	//Convolve them
	Duplicate/O G,G_conv
	Convolve/A signal, G_conv
	
	//Make the result wave by sampling at the correct x values.
	yw = G_conv(xw[p])+pw[2]
End

Function conv_Stretchedexp(pw, yw, xw) : FitFunc
//A stretched exponential function convoluted with a gaussian instrument response
	//my passed waves
	WAVE pw, yw, xw
	
	//| A function of an arbitrary number of exponentials convoluted with a 
	//| gaussian IRF
	//| 
	//| Coefficients 5
	//| pw[0] = W = FWHM/2*sqrt(ln(2)) of the IRF
	//| pw[1] = t0
	//| pw[2] = y0
	//| pw[3] = H offset for any photoproduct formation
	//| pw[4] = a1
	//| pw[5] = t1
	//| ... and so on
	
	//We will use a tenth of the pw[0] as our resolution
	Variable dT = max(1,abs(pw[0]/40))//fs
	
	//This ensures that there is enough time for the gaussian to decay completely
	Variable maxT = wavemax(xw)*1.1
	Variable minT = wavemin(xw)
	
	Variable numberOfPoints = round(2*maxT/dT)+1
	
	//Make my signal wave and my gaussian wave
	Make/D/O/N=(numberOfPoints) G, signal=0
	//Set my scales, to make things easier
	SetScale/P x -maxT, dT,"", G, signal//Needs to be symmetric about zero for later
	
	//Make my instrument response function, in this case modeled by a gaussian
	// with FWHM=2*sqrt(ln(2))*pw[0]
	//We'll make this gaussian in the middle of the wave to help with the convolution
	G = exp(-(x/pw[0])^2)
	//Normalize it
	Variable GaussArea=sum(G)
	
	//now we can make the signal
	
	Signal = (x>=0)*pw[3]
	variable npts= numpnts(pw),i=4
	do
		if( i>=npts )
			break
		endif
		Signal += StretchedExp(pw[i],pw[i+1],pw[i+2],x)
		i+=3
	while(1)
	
	//Convolve them
	Duplicate/O G,G_conv
	Convolve/A signal, G_conv
	
	//Correct t0 due to the effects of the varying resolution (i.e. dT)
	variable t0 = pw[1]+(0.00397793573758223+0.0124675379745083*dT*40)
	//Make the result wave by sampling at the correct x values.
	yw = G_conv(xw[p]-t0)/GaussArea+pw[2]
End

STATIC Function StretchedExp(A,tau,B,t)
//A function that conveniently calculates the stretched exponential function for use above
	Variable A, tau, B, t
	if(t>=0)
		Return A*exp(-(t/tau)^B)
	Else
		Return 0
	EndIf
End

Function conv_Bi(pw, yw, xw) : FitFunc
//This convolutes a gaussian instrument response with a bimolecular decay
//i.e. A+B->C
	//my passed waves
	WAVE pw, yw, xw
	
	//| A function of an arbitrary number of exponentials convoluted with a 
	//| gaussian IRF
	//| 
	//| Coefficients 5
	//| pw[0] = W = FWHM/2*sqrt(ln(2)) of the IRF
	//| pw[1] = t0
	//| pw[2] = y0
	//| pw[3] = H offset for any photoproduct formation
	//| pw[4] = a1
	//| pw[5] = t1
	//| ... and so on
	
	//We will use a tenth of the pw[0] as our resolution
	Variable dT = max(1,abs(pw[0]/40))//fs
	
	//This ensures that there is enough time for the gaussian to decay completely
	Variable maxT = wavemax(xw)*1.1
	Variable minT = wavemin(xw)
	
	Variable numberOfPoints = round(2*maxT/dT)+1
	
	//Make my signal wave and my gaussian wave
	Make/D/O/N=(numberOfPoints) G, signal=0
	//Set my scales, to make things easier
	SetScale/P x -maxT, dT,"", G, signal//Needs to be symmetric about zero for later
	
	//Make my instrument response function, in this case modeled by a gaussian
	// with FWHM=2*sqrt(ln(2))*pw[0]
	//We'll make this gaussian in the middle of the wave to help with the convolution
	G = exp(-(x/pw[0])^2)
	//Normalize it
	Variable GaussArea=sum(G)
	
	//now we can make the signal
	
	Signal = (x>=0)*pw[3]
	variable npts= numpnts(pw),i=4
	do
		if( i>=npts )
			break
		endif
		Signal += myBi(pw[i],pw[i+1],x)
		i+=2
	while(1)
	
	//signal/=2/sqrt(ln(2))
	
	//Convolve them
	Duplicate/O G,G_conv
	Convolve/A signal, G_conv
	
	//Make the result wave by sampling at the correct x values.
	variable t0 = pw[1]+(0.00397793573758223+0.0124675379745083*pw[0])
	yw = G_conv(xw[p]-t0)/GaussArea+pw[2]
End

STATIC Function myBi(A,tau,t)
	Variable A, tau, t
	if(t>=0)
		Return 1/(1/A+t/tau)
	Else
		Return 0
	EndIf
End

Function conv_exp(pw, yw, xw)
//An all-at-once version of a convoluted exponential.
	//my passed waves
	WAVE pw, yw, xw
	
	//| A function of an arbitrary number of exponentials convoluted with a 
	//| gaussian IRF
	//| 
	//| Coefficients 5
	//| pw[0] = W = FWHM/2*sqrt(ln(2)) of the IRF
	//| pw[1] = t0
	//| pw[2] = y0
	//| pw[3] = H offset for any photoproduct formation
	//| pw[4] = a1
	//| pw[5] = t1
	//| ... and so on
	
	//We will use a tenth of the pw[0] as our resolution
	Variable dT = max(1,abs(pw[0]/40))//fs
	
	//This ensures that there is enough time for the gaussian to decay completely
	Variable maxT = wavemax(xw)*1.1
	Variable minT = wavemin(xw)
	
	Variable numberOfPoints = round(2*maxT/dT)+1
	
	//Make my signal wave and my gaussian wave
	Make/D/O/N=(numberOfPoints) G, signal=0
	//Set my scales, to make things easier
	SetScale/P x -maxT, dT,"", G, signal//Needs to be symmetric about zero for later
	
	//Make my instrument response function, in this case modeled by a gaussian
	// with FWHM=2*sqrt(ln(2))*pw[0]
	//We'll make this gaussian in the middle of the wave to help with the convolution
	G = exp(-(x/pw[0])^2)
	//Normalize it
	Variable GaussArea=sum(G)
	
	//now we can make the signal
	
	Signal = (x>=0)*pw[3]
	variable npts= numpnts(pw),i=4
	do
		if( i>=npts )
			break
		endif
		Signal += myExp(pw[i],pw[i+1],x)
		i+=2
	while(1)
	
	//Convolve them
	Duplicate/O G,G_conv
	Convolve/A signal, G_conv
	
	//Make the result wave by sampling at the correct x values.
	variable t0 = pw[1]+(0.00397793573758223+0.0124675379745083*pw[0])
	yw = G_conv(xw[p]-t0)/GaussArea+pw[2]
End

STATIC Function myExp(A,tau,t)
	Variable A, tau, t
	if(t>=0)
		Return A*exp(-t/tau)
	Else
		Return 0
	EndIf
End

Function conv_exp1(w,t) : FitFunc
//A single exponential convoluted with a gaussian instrument response function
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation=A*Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2))*(1 + Erf((t-t0)/w - w/(2*t1)))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	t-=w[1]
	
	Variable val=Exp((-4*t*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf(t/w[0] - w[0]/(2*w[4]), 1e-16))
	val*=w[3]
	val/=2
	val+=w[2]
	return val
End

Function conv_exp2(w,t) : FitFunc
//Two exponentials convoluted with a gaussian instrument response function
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation = A0*Exp((-4*t*t1 + 4*t0*t1 + w^2)/(4*t1^2))*(1+Erf((t -2t0)/w - w/(2*t1)))+
	//|					A1*Exp((-4*t*t2 + 4*t0*t2 + w^2)/(4*t2^2))*(1+Erf((t -2t0)/w - w/(2*t2)))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//|  w[5]=A1scaling coefficient
	//|  w[6]=t2, the decay (1/e) time of the first exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	//CurveFitDialog/ w[5] = a2
	//CurveFitDialog/ w[6] = t2
	
	Variable val= w[3]*Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[4])))
	val+=w[5]*Exp((-4*(t-w[1])*w[6] + w[0]^2)/(4*w[6]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[6])))
	val/=2
	val+=w[2]
	
	return val
End

Function conv_exp3(w,t) : FitFunc
//Three exponentials convoluted with a gaussian instrument response function
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation = A0*Exp((-4*t*t1 + 4*t0*t1 + w^2)/(4*t1^2))*(1+Erf((t -2t0)/w - w/(2*t1)))+
	//|					A1*Exp((-4*t*t2 + 4*t0*t2 + w^2)/(4*t2^2))*(1+Erf((t -2t0)/w - w/(2*t2)))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//|  w[5]=A1 scaling coefficient
	//|  w[6]=t2, the decay (1/e) time of the second exponential
	//|  w[7]=A2 scaling coefficient
	//|  w[8]=t3, the decay (1/e) time of the third exponential
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 9
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	//CurveFitDialog/ w[5] = a2
	//CurveFitDialog/ w[6] = t2
	//CurveFitDialog/ w[7] = a3
	//CurveFitDialog/ w[8] = t3
	
	Variable val= w[3]*Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[4])))
	val+=w[5]*Exp((-4*(t-w[1])*w[6] + w[0]^2)/(4*w[6]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[6])))
	val+=w[7]*Exp((-4*(t-w[1])*w[8] + w[0]^2)/(4*w[8]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[8])))
	val/=2
	val+=w[2]
	return val
	
End

Function conv_expH(w,t) : FitFunc
//three exponentials and a step function convoluted with a gaussian instrument response function
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation = A*(1 + Erf((t-t0)/w) +Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2)) (1 + Erf((t-t0)/w - w/(2*t1))))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//|  w[5]=A1 scaling coefficient
	//|  w[6]=t2, the decay (1/e) time of the second exponential
	//|  w[7]=A2 scaling coefficient
	//|  w[8]=t3, the decay (1/e) time of the third exponential
	//|  w[9]= Heaviside Amplitude
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 10
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	//CurveFitDialog/ w[5] = a2
	//CurveFitDialog/ w[6] = t2
	//CurveFitDialog/ w[7] = a3
	//CurveFitDialog/ w[8] = t3
	//CurveFitDialog/ w[9] = H
	
	Variable val= w[3]*Exp((-4*(t-w[1])*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[4])))
	val+=w[5]*Exp((-4*(t-w[1])*w[6] + w[0]^2)/(4*w[6]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[6])))
	val+=w[7]*Exp((-4*(t-w[1])*w[8] + w[0]^2)/(4*w[8]^2))*(1 + Erf((t-w[1])/w[0] - w[0]/(2*w[8])))
	val+=w[9]*(1 + Erf((t-w[1])/w[0]))
	val/=2
	val+=w[2]
	return val
End

Function conv_expHA(w,t) : FitFunc
//three exponentials, a step function and a delta function
//convoluted with a gaussian instrument response function
	Wave w
	Variable t
	
	//|  w contains the parameters needed for the convolution and exponential functions.
	//| Equation = A*(1 + Erf((t-t0)/w) +Exp((-4*(t-t0)*t1 + w^2)/(4*t1^2)) (1 + Erf((t-t0)/w - w/(2*t1))))
	//|  w[0]=w the gaussian width(measured by cross-correlation)
	//|  w[1]=t0 t-zero
	//|  w[2]=y0 offset
	//|  w[3]=A0 scaling coefficient
	//|  w[4]=t1, the decay (1/e) time of the first exponential
	//|  w[5]=A1 scaling coefficient
	//|  w[6]=t2, the decay (1/e) time of the second exponential
	//|  w[7]=A2 scaling coefficient
	//|  w[8]=t3, the decay (1/e) time of the third exponential
	//|  w[9]= Heaviside Amplitude
	//|  w[9]= amplitude of the instantaneous response
	//CurveFitDialog/ 
	//CurveFitDialog/ Coefficients 11
	//CurveFitDialog/ w[0] = W
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = y0
	//CurveFitDialog/ w[3] = a1
	//CurveFitDialog/ w[4] = t1
	//CurveFitDialog/ w[5] = a2
	//CurveFitDialog/ w[6] = t2
	//CurveFitDialog/ w[7] = a3
	//CurveFitDialog/ w[8] = t3
	//CurveFitDialog/ w[9] = H
	//CurveFitDialog/ w[10] = A
	
	t-=w[1]
	
	Variable val= w[3]*Exp((-4*t*w[4] + w[0]^2)/(4*w[4]^2))*(1 + Erf(t/w[0] - w[0]/(2*w[4])))
	val+=w[5]*Exp((-4*t*w[6] + w[0]^2)/(4*w[6]^2))*(1 + Erf(t/w[0] - w[0]/(2*w[6])))
	val+=w[7]*Exp((-4*t*w[8] + w[0]^2)/(4*w[8]^2))*(1 + Erf(t/w[0] - w[0]/(2*w[8])))
	val+=w[9]*(1 + Erf(t/w[0]))
	val/=2
	val+=w[10]*Exp(-(t/w[0])^2)
	val+=w[2]
	return val
End

Function conv_expMulti(w,t) : FitFunc
//A function that can model an arbitrary number of exponentials and a single step
//function convoluted with a gaussian instrument response function.
	Wave w
	Variable t
	
	//Analytical convoluted exponential function with as many exponentials as you want
	// and an offset, not usuable through the fitting menu
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
	//CurveFitDialog/ w[4] = a1
	//CurveFitDialog/ w[5] = t1
	t-=w[1]
	Variable val=w[3]*(1 + Erf(t/w[0])),npts=numpnts(w),i=4
	do
		if( i>=npts )
			break
		endif
		val+= w[i]*Exp((-4*t*w[i+1] + w[0]^2)/(4*w[i+1]^2))*(1 + Erf(t/w[0] - w[0]/(2*w[i+1])))
		i+=2
	while(1)
	val/=2
	val+=w[2]
	return val
End

Function triExp_XOffSet(w,x) : FitFunc
//Simple triple exponential funtion, trying to mimic the built in dblexp_XOffset
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0 + A1*Exp(-(x-x0)/t1) + A2*Exp(-(x-x0)/t2) + A3*Exp(-(x-x0)/t3)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ w[0] = y0
	//CurveFitDialog/ w[1] = x0
	//CurveFitDialog/ w[2] = A1
	//CurveFitDialog/ w[3] = t1
	//CurveFitDialog/ w[4] = A2
	//CurveFitDialog/ w[5] = t2
	//CurveFitDialog/ w[6] = A3
	//CurveFitDialog/ w[7] = t3

	return w[0] + w[2]*Exp(-(x-w[1])/w[3]) + w[4]*Exp(-(x-w[1])/w[5]) + w[6]*Exp(-(x-w[1])/w[7])
End
