: kerg_sncda.mod
:
: KERG current model
:
: Current is defined by
:	i = g * o (v-e)
:
: kinetic scheme is:
:       C <-> O <-> I

: Written by Tim Rumbell, IBM Research, 2017, thrumbel@us.ibm.com 
: Based on model of Yu and Canavier 2015
: See also Ficker 1998, Wang 1997, Canavier 2007, Ji 2012

NEURON {
	SUFFIX kerg
	USEION k READ ek WRITE ik
    RANGE o, i, ik
    RANGE gbar, alphaa, betaa, alphai, betai
    EXTERNAL apc_metap, fpc_metap
    GLOBAL Vhalf, taumod
    GLOBAL vshift
    GLOBAL aa0, aac, ba0, bac
    GLOBAL ai0, aic, bi0, bic
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	ek		        (mV)

    gbar   = 0.05 	(S/cm2)
    
    Vhalf   = -40   (mV)
    taumod  = 1

    vshift = 0 (mV)

    ai0 = 91.11
    aic = 0.1189
    bi0 = 12.6
    bic = 0.0733
    aa0 = 0.00236
    aac = 0.0759
    ba0 = 1.2523e-05
    bac = -0.0671
}

STATE {
	o i
}

ASSIGNED {
	v		        (mV)
	ik		(mA/cm2)
	alphaa  (1/ms)  
	alphai  (1/ms)  
	betaa  (1/ms)  
	betai  (1/ms)  
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    ik = gbar * o * (v-ek)
}

DERIVATIVE states{
	rates(v)
    o' = (alphaa * (1 - o - i)) + (betai * i) - ((alphai + betaa) * o)
	i' = (alphai * o) - (betai * i)
}

INITIAL {
    setVhalf(Vhalf)
    setTauMod(taumod)
	rates(v)
	o = 0
    i = 0
}

PROCEDURE rates(v(mV)) {
    alphaa = aa0 * exp(aac * (v + vshift))
    betaa = ba0 * exp(bac * (v + vshift))
    alphai = ai0 * exp(aic * (v + vshift))
    betai = bi0 * exp(bic * (v + vshift))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vshift = Vhalf
}

PROCEDURE setTauMod(taumod) {
    aa0 = 0.00236 / taumod
    ba0 = 1.2523e-05 / taumod
    ai0 = 91.11 / taumod
    bi0 = 12.6 / taumod
}

