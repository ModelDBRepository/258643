: kv2_dr_sncda.mod
:
: Based on the kv4 model of:
: Josh Held	j-held@northwestern.edu
: 3/2003
:
: Kv2 current (DR current) model
:
: Current is defined by
:	i = g * m^22 * (v-e)
:
: Edited by Tim Rumbell, IBM Research, 2017, thrumbel@us.ibm.com
: Parameters edited to match those found in
: Kimm et al, J Neuroscience, 2015
: Recordings from SNc dopaminergic neurons


NEURON {
	SUFFIX kv2_dr
	USEION k READ ek WRITE ik
    RANGE minf, tm, ik
    RANGE gbar
    : EXTERNAL pc1_metap, pc3_metap, pc24_metap
    GLOBAL Vhalf, taumod
    GLOBAL vhm, vcm
    GLOBAL vhtm, atm, Ctm, tm0
    GLOBAL gax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	ek		        (mV)

    gbar   = 0.05 	(S/um2)
    gax    = 1
    
    Vhalf   = -30.5   (mV)
    taumod  = 1

    vhm     = -30.5   (mV)
    vcm     = -13    (mV)

    vhtm    = -9 (mv)
    atm     = 10.0  (mV)
    Ctm     = 3.6  (mV)
    tm0     = 3.6 (ms)
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
    ik = gbar * m^4 * (v-ek)
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
    vhtm = Vhalf + 21.5
}

PROCEDURE setTauMod(taumod) {
    tm0 = 3.6 * taumod    
    Ctm = 3.6 * taumod
}

