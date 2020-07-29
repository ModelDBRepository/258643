
COMMENT

cat_sncda.mod

T-type Calcium channel, Hodgkin-Huxley style kinetics.  

Kinetics are from Poetschke et al 2015 Scientific Reports paper
on CaT in rodent SNc DA cells

Adapted from Canavier 2014 NaT channel model due to both having
fast and slow components of inactivation - According to the Poetschke 2015 
paper, which provides incomplete data, the equation structure should 
be similar
Added ghk for driving force

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu
Qian...Canavier, 2014
Tim Rumbell, IBM Research, 2016, thrumbel@us.ibm.com

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cat
	USEION ca READ eca,cai,cao WRITE ica
	RANGE m, h, hs, gcat, gbar, ica
    EXTERNAL apc_metap, fpc_metap
	GLOBAL vhm, vhh, vhhs, km, kh, khs
    GLOBAL tm0, th0, ths0, vhtm, vhth, atm, ath, Ctm, Cth
	RANGE minf, hinf, hsinf, mtau, htau, hstau
	GLOBAL q10, tadj
    GLOBAL Vhalf, taumod 
    : Vhalf and taumod are modifiers that can be used as external parameters, varying 
    : several of the interal parameters simultaneously to modulate the channel 
    :properties up or down within the experimentally observed range
}

PARAMETER {
	gbar = 1   	(cm/s) : (S/cm2)
	Vhalf = -40	(mV)		: voltage shift (affects ll)
    taumod = 1
    cai         (mM)
    cao         (mM)
    
    : Parameters from combination of:
    : Poetschke 2015 - activation kinetics, and target for tuning of inactivation
    : McRory 2001 - starting point for inactivation tuning was inactivation of alpha_1H
    : alpha_1H in McRory looks like closest match to kinetic parameters in Poetschke
    : Poetschke shows more expression of Cav3.2 in WT == alpha_1H (I think...)
    : Q10 of 3 used to approximate 0.35 ratio for temp dependence found in Iftinca 2006
    : where Tau recovery from inactivation was 184 ms @ 37C and 527 ms @ 21C

    vhm  = -54.5 (mV)    : v 1/2 for act
    km   = 5   (mV)    : slope for act

	vhh  = -64.5	(mV)	: v 1/2 for inact 	
    kh   = -1.6 (mV)        : slope for inact 
    
    vhhs = -64.5	(mV)	: v 1/2 for inact 	
    khs  = -1.6    (mV)        : slope for inact 

    tm0  = 3.2
    th0  = 76
    ths0 = 600

    vhtm  = -40 (mV)
    vhth  = -46 (mV)

    atm   = 4.6
    ath   =8.85

    Ctm   = 19
    Cth   = 43

	temp = 33	(degC)		: original temp 
	q10  = 3.0			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)

    tadj = 1

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    F = 9.6485e4    (coul)
    R = 8.3145      (joule/degC)
	(S) = (siemens)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ica 		(mA/cm2)
	gca		(cm/s) : (S/um2)
	eca		(mV)	
	minf 		hinf        hsinf
	mtau (ms)	htau (ms)   hstau (ms)
    T   (degC)
    E   (volts)
    z
}
 

STATE { m h hs }

INITIAL { 
    : Assume that v has been constant for long enough to reach steady state
    setVhalf(Vhalf)
    setTauMod(taumod)
	m = minf_cat(v)
	h = hinf_cat(v)
    hs = hsinf_cat(v)
    tadj = tadj_ca_t()
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = gbar*m*m*m*((0.6*h)+(0.4*hs))
	    ica = gca * ghk(v,cai,cao) : (v - eca)
} 

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25.26 (mV) /293.15 (degC) )*(celsius + 273.15 (degC) ))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

DERIVATIVE states {   :Computes state variables m, h, and hs 
        m' = -(m-minf_cat(v))/(tadj*taum_cat(v))
        h' = -(h-hinf_cat(v))/(tadj*tauh_cat(v))
        hs' = -(hs-hsinf_cat(v))/(tadj*tauhs_cat(v))
}

FUNCTION boltz(x,y,z) {
		boltz = 1/(1+exp(-(x-y)/z))
}

FUNCTION minf_cat(v (mV)) (1) {
        minf_cat = boltz(v,vhm,km) 
}

FUNCTION taum_cat(v (mV)) (1/ms) {
        taum_cat = tm0 + Ctm/(1+exp((v-vhtm)/atm))
}

FUNCTION hinf_cat(v (mV)) (1) {
        hinf_cat = (boltz(v,vhh,kh))
}

FUNCTION tauh_cat(v (mV)) (1/ms) {
        tauh_cat = th0 + Cth/(1+exp((v-vhth)/ath))
}

FUNCTION hsinf_cat(v (mV)) (1) {
        hsinf_cat = (boltz(v,vhhs,khs))
}

FUNCTION tauhs_cat(v (mV)) (1/ms) {
        tauhs_cat = ths0
}

FUNCTION tadj_ca_t() {
        tadj_ca_t =  1/(q10^((celsius - temp)/10))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf
    vhh = Vhalf - 10
    vhhs = Vhalf - 10
    vhtm = Vhalf + 14.5
    vhth = Vhalf + 8.5
}

PROCEDURE setTauMod(taumod) {
    tm0 = 3.2 * taumod
    th0 = 76 * taumod
    ths0 = 600 * taumod
    Ctm = 19 * taumod
    Cth = 43 * taumod
}

