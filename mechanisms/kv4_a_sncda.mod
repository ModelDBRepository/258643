: kv4_a_sncda.mod
:
: Josh Held	j-held@northwestern.edu
: 3/2003
:
: Kv4 current (A current) model
:
: Current is defined by
:	i = g * m^3 * h * (v-e)
:
: Data from 9_13_2_1
:
: Edited by Tim Rumbell, IBM Research, 2017, thrumbel@us.ibm.com
: Parameters edited to match those found in
: Amendola et al 2012 (J Neurosci)
: Recordings from SNc dopaminergic neurons


NEURON {
	SUFFIX kv4_a
	USEION k READ ek WRITE ik
    RANGE minf, tm, hinf, th, ik, th
    RANGE gbar
    EXTERNAL apc_metap, fpc_metap
    GLOBAL Vhalf, taumod
    GLOBAL vhm, vcm
    GLOBAL vhh, vch
    GLOBAL vhtm, atm, Ctm, tm0
    GLOBAL vhth, ath, Cth, th0
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

    vhm     = -40   (mV)
    vcm     = -7    (mV)

    vhh     = -73   (mV)
    vch     = 4.9   (mv)

    vhtm    = -56.7 (mv)
    atm     = 6.22  (mV)
    Ctm     = 4.83  (mV)
    tm0     = 1.029 (ms)

    vhth    = -68.5 (mv)
    ath     = 5.95  (mV)
    Cth     = 78.4  (mV)
    th0     = 39.04 (ms)
}

STATE {
	m h
}

ASSIGNED {
	v		        (mV)
	ik		(mA/cm2)
	minf
	tm		(ms)
	hinf
	th		(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    ik = gbar * m^3 * h * (v-ek)
    :ik = gbar * m * h * (v-ek)
}

DERIVATIVE states{
	rates(v)
	m' = (minf - m)/tm
    h' = (hinf - h)/th
}

INITIAL {
    setVhalf(Vhalf)
    setTauMod(taumod)
	rates(v)
	m = minf
    h = hinf
}

PROCEDURE rates(v(mV)) {LOCAL q10
    q10 = 3^((celsius-32)/10)
    minf = (1/(1 + exp((v-vhm)/vcm)))^(1/3) : cube-rooted because we want to match Amendola parameters for Vhalf and slope when we cube this gating variable due to m^3 in BREAKPOINT
    hinf = 1/(1 + exp((v-vhh)/vch))
    tm = (1/q10)*(tm0 + Ctm/(1 + exp((v-vhtm)/atm)))
    th = (1/q10)*(th0 + Cth/(1 + exp((v-vhth)/ath)))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf
    vhh = Vhalf - 33
}

PROCEDURE setTauMod(taumod) {
    tm0 = 1.029 * taumod    
    Ctm = 4.83 * taumod
    
    th0 = 39.04 * taumod
    Cth = 78.4 * taumod
}

