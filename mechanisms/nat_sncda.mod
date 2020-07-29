
COMMENT

nat.mod

Sodium channel, Hodgkin-Huxley style kinetics.  

Kinetics are from Qian...Canavier 2014 J Neurophysiol paper on NaT in DA neurons
Adapted by Yu 2014 and Yu 2015 - claims these fit data more accurately...

This model is used in a few Canavier papers on DAs from 2011-2015

Original fit is to Seutin Engel 2010

Then with additional slow component of inactivation
based on Ding 2011

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu

Edited by Tim Rumbell, IBM Research, 2019, thrumbel@us.ibm.com

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nat
	USEION na READ ena WRITE ina
	RANGE m, h, hs, gna, gbar, ina
    RANGE vhm, vhh, vhhs, vh_shift
	GLOBAL km, kh, khs
    GLOBAL tm0, tm1, tmac, tma0, tma1, tmb0, tmb1
    GLOBAL thf0, thf1, thfa0, thfa1, thfb0, thfb1
    GLOBAL ths0, ths1 :, thsa0, thsa1, thsb0, thsb1
	RANGE minf, hinf, hsinf, mtau, htau, hstau
	GLOBAL q10, tadj, vshift
    GLOBAL Vhalf, taumod, vhm_ax, gdend, gax, vhh_shift, vhhs_shift
}

PARAMETER {
	gbar = 1   	(S/cm2)
	vshift = 0	(mV)		: voltage shift (affects ll)
        vh_shift = 0    (mV) 		: shifts vh values per section
    
    : parmeters for easy tuning
    Vhalf = -30.0907
    taumod = 1
    vhm_ax = 0
    gdend  = 1
    gax    = 1
    vhh_shift = 0
    vhhs_shift = 0
    
    : Parameters from Qian 2014 Table 1

    : activation parameters (m)
    vhm  = -30.0907 (mV)    : v 1/2 for act
    km   = 9.7264   (mV)    : slope for act

    : parameters for tau m
    tm0     = 0.01
    tm1     = 1.0
    : tma0    = 15.6504
    : tma1    = 0.4043
    tmac    = 0.79992
    tma0    = -19.565
    tma1    = -0.50542
    tmb0    = 3.0212
    tmb1    = -7.463e-03

    : fast inactivation parameters (hf)
	vhh  = -54.0289	(mV)	: v 1/2 for inact 	
    kh   = -10.7665 (mV)        : slope for inact 
    
    : parameters for tau hf
    thf0    = 0.4
    thf1    = 1.0
    thfa0   = 5.0754e-04
    thfa1   = -6.3213e-02
    thfb0   = 9.7529
    thfb1   = 0.13442

    : slow inactivation parameters (hs)
    vhhs = -54.8	(mV)	: v 1/2 for inact 	
    khs  = -1.57    (mV)        : slope for inact 

    : parameters for tau hs
    : From Yu 2014:
    ths0    = 20
    ths1    = 580

	temp = 24	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)

    tadj = 1

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(S/um2)
	ena		(mV)	
	minf 		hinf        hsinf
	mtau (ms)	htau (ms)   hstau (ms)
}
 

STATE { m h hs }

INITIAL { 
    : Assume that v has been constant for long enough to reach steady state
    setVhalf(Vhalf)
    setTauMod(taumod)
	m = minf_na(v+vshift)
	h = hinf_na(v+vshift)
    hs = hsinf_na(v+vshift)
    tadj = tadj_na()
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gbar*m*m*m*h*hs
	    ina = gna * (v - ena)
} 

DERIVATIVE states {   :Computes state variables m, h, and hs 
        m' = -(m-minf_na(v+vshift))/(tadj*taum_na(v+vshift))
        h' = -(h-hinf_na(v+vshift))/(tadj*tauh_na(v+vshift))
        hs' = -(hs-hsinf_na(v+vshift))/(tadj*tauhs_na(v+vshift))
}

FUNCTION boltz(x,y,z) {
		boltz = 1/(1+exp(-(x-y)/z))
}

FUNCTION minf_na(v (mV)) (1) {
        minf_na = boltz(v,vhm,km) 
}

FUNCTION taum_na(v (mV)) (1/ms) {
        LOCAL a, b
        a = tmac*(tma0 + (tma1*v)) / ( exp( tma0 + (tma1 * v) ) - 1 )
        b = tmb0 * exp(tmb1*v)
        taum_na = tm0 + ( tm1 / (a+b) )
}

FUNCTION hinf_na(v (mV)) (1) {
        hinf_na = boltz(v,vhh,kh) 
}

FUNCTION tauh_na(v (mV)) (1/ms) {
        LOCAL a, b
        a = thfa0 * exp(thfa1 * v)
        b = thfb0 * exp(thfb1 * v)
        tauh_na = thf0 + ( thf1 / (a+b) )
}

FUNCTION hsinf_na(v (mV)) (1) {
        hsinf_na = boltz(v,vhhs,khs) 
}

FUNCTION tauhs_na(v (mV)) (1/ms) {
        tauhs_na = ths0 + ths1 / (1 + exp(v) )
}

FUNCTION tadj_na() {
        tadj_na =  1/(q10^((celsius - temp)/10))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf + vh_shift
    vhh = Vhalf - 23.9382 + vhh_shift + vh_shift
    vhhs = Vhalf - 24.7093 + vhhs_shift + vh_shift
}

PROCEDURE setTauMod(taumod) {
    tm0 = 0.01 * taumod
    tm1 = 1.0 * taumod

    thf0 = 0.4 * taumod
    thf1 = 1 * taumod

    ths0 = 20 * taumod
    ths1 = 580 * taumod
}

