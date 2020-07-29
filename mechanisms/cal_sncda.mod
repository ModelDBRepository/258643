NEURON {
	SUFFIX cal
	USEION ca READ eca,cai,cao WRITE ica
    RANGE  gbar, ica, g
    RANGE m, h, f
    EXTERNAL apc_metap, fpc_metap
    GLOBAL kf, tauf
	GLOBAL vhm, km, Ctm, tm0, vhtm, atm
    GLOBAL vhh, kh, th0
    GLOBAL q10, tadj
    RANGE minf, taum, hinf, tauh, finf
    GLOBAL Vhalf, taumod : modifiers that can be used externally to shift the channel around
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    F = 9.6485e4    (coul)
    R = 8.3145      (joule/degC)
	(S) = (siemens)
	(mM) = (milli/liter)
}

PARAMETER {
	v		(mV)
	celsius	(degC)
	gbar = .003 (cm/s) :  	(S/cm2)
	cai 		(mM)
	cao 		(mM)

	: ki = .001 	(mM)
    kf = .001 	(mM)
    tauf = 30   (ms)

	vhm = -35	(mV)
    km = 5.5    (mV)
	Ctm = 0.2		(ms)
	atm = 12    	(mV)

	vhtm = -35	(mV)

    tm0 = 0.2   (ms)

	vhh = -55	(mV)
    kh = -8     (mV)
	th0 = 300		(ms)

    Vhalf = -35 (mV)
    taumod = 1  

    temp = 30
    q10 = 2
    
    tadj = 1
}

STATE { m h f }

ASSIGNED {
	vf	(mV)
	ica	(mA/cm2)
    g	(cm/s) : (S/cm2)
    eca (mV)
    minf
    taum	(ms)
    hinf
    tauh    (ms)
	a	(1/ms)
    T   (degC)
    E   (volts)
    z
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    g = gbar*m*m*h*f
	ica  = g*ghk(v,cai,cao) : (v-eca)
}

INITIAL {
    setVhalf(Vhalf)
    setTauMod(taumod)
	m = minf_cal(v)
    h = hinf_cal(v)
    f = finf_cal(cai)
    tadj = tadj_ca_l()
}

FUNCTION fca(cai(mM)) {
    : fca is the CDI function
	fca = kf/(kf+cai)
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

DERIVATIVE states { :Computes state variables m and h
        m' = (minf_cal(v) - m)/(tadj*taum_cal(v))
        h' = (hinf_cal(v) - h)/(tadj*tauh_cal(v))
        f' = (finf_cal(cai) - f)/(tadj*tauf)
}

FUNCTION boltz(x,y,z) {
    boltz = 1/(1+exp(-(x-y)/z))
}

FUNCTION finf_cal(cai (mM)) (1) {
    finf_cal = 1/(1+(cai/kf))
}

FUNCTION minf_cal(v (mV)) (1) {
    minf_cal = boltz(v,vhm,km)
}

FUNCTION taum_cal(v (mV)) (1/ms) {
    taum_cal = tm0 + Ctm/(1+exp((v-vhtm)/atm))
}

FUNCTION hinf_cal(v (mV)) (1) {
    hinf_cal = 0.2 + (0.8*boltz(v,vhh,kh))
}

FUNCTION tauh_cal(v (mV)) (1/ms) {
    : tauh_cal = th0 + Cth/(1+exp((v-vhth)/ath))
    tauh_cal = th0
}

FUNCTION tadj_ca_l() {
    tadj_ca_l = 1/(q10^((celsius - temp)/10))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf
    vhh = Vhalf - 21 : determined by difference between V50 m and h in Koschak 2001, as per method of Tuckwell 2012
    vhtm = vhm
}

PROCEDURE setTauMod(taumod) {
    tm0 = 0.2 * taumod
    Ctm = 0.2 * taumod
    
    th0 = 300 * taumod
    
    tauf = 30 * taumod
}

