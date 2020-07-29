: bk_sncda.mod
:
: Based on the kv4 model of:
: Josh Held	j-held@northwestern.edu
: 3/2003
:
: BK current model
:
: Current is defined by
:	i = g * m^2 * (v-e)
:
: Edited by Tim Rumbell thrumbel@us.ibm.com
: Parameters edited to match those found in
: Kimm et al, J Neuroscience, 2015
: This is the voltage-dependent component of BK, but
: without accounting for Ca changes during the AP, so
: can be considered to have a gating 'fudge factor'
: that accounts for combined V and Ca dependence
: Recordings from SNc dopaminergic neurons


NEURON {
	SUFFIX bk
	USEION k READ ek WRITE ik
    RANGE minf, tm, ik
    RANGE gbar
    GLOBAL Vhalf, taumod
    GLOBAL vhm, vcm
    GLOBAL vhtm, atm, Ctm, tm0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	ek		        (mV)

    gbar   = 0.05 	(mho/cm2)
    
    Vhalf   = -16   (mV)
    taumod  = 1

    vhm     = -16   (mV)
    vcm     = -8.5    (mV)

    vhtm    = -16 (mv)
    atm     = 10.0  (mV)
    Ctm     = 0.13  (mV)
    tm0     = 0.87 (ms)
}

STATE {
	m
}

ASSIGNED {
	v		        (mV)
	ik		(mA/cm2)
	minf
	tm		(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    : ik = gbar * m^2 * (v-ek)
    : Not sure about power, so using 1 so that Vhalf and slope match Kimm 2015
    ik = gbar * m * (v-ek)
}

DERIVATIVE states{
	rates(v)
	m' = (minf - m)/tm
}

INITIAL {
    setVhalf(Vhalf)
    setTauMod(taumod)
	rates(v)
	m = minf
}

PROCEDURE rates(v(mV)) {LOCAL q10
    q10 = 3^((celsius-34)/10)
    minf = 1/(1 + exp((v-vhm)/vcm))
    tm = (1/q10)*(tm0 + Ctm/(1 + exp((v-vhtm)/atm)))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf
    vhtm = Vhalf
}

PROCEDURE setTauMod(taumod) {
    tm0 = 0.87 * taumod    
    Ctm = 0.13 * taumod
}

