: Model and most parameters from Wang, Chen, Nolan and Siegelbaum, Neuron, 2002
: Some parameters tuned to match findings of Gambardella, Pignatelli and Belluzi 2012
: The h-current in the Substantia Nigra pars Compacta Neurons : A Re-examination
: PLoS ONE December 2012 7:12 e52329
: Adapted by Tim Rumbell, 2017, thrumbel@us.ibm.com

NEURON {
	SUFFIX hcn
	NONSPECIFIC_CURRENT i
	RANGE i, ehcn, g, gbar
    EXTERNAL apc_metap, fpc_metap
	GLOBAL a0, b0, ah, bh, ac, bc, aa0, ba0
	GLOBAL aa0, ba0, aah, bah, aac, bac
	GLOBAL kon, koff, b, bf, gca, shift
    GLOBAL Vhalf, vh1, vh2, vh_shift, avh1, avh2, avh_shift
	RANGE ai
}

UNITS {
	(mV)	= (millivolt)
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA)	= (milliamp)
	(S)	= (siemens)
}

PARAMETER {
	gbar    = 0.0		(S/cm2)
    ehcn    = -44       (mV) : parameter from Gambardella 2012
    : : Params established using least squares non linear fitting from orig params +- 25%
    : a0      = .0023		(/ms)		: parameters for alpha and beta
    : b0      = .0022		(/ms)
    : ah      = -93.9898	(mV)
    : bh      = -62.7310		(mV)
    : ac      = -0.1137		(/mV)
    : bc      = 0.1203		(/mV)
    : aa0     = 0.1203	(/ms)		: parameters for alphaa and betaa
    : ba0     = 0.0042		(/ms)
    : aah     = -109.7611		(mV)
    : bah     = -27.6119		(mV)
    : aac     = -0.0678		(/mV)
    : bac     = 0.07		(/mV)
    : kon     = 1.327		(/mM-ms)	: cyclic AMP binding parameters
    : koff    = 1.217e-05	(/ms)
    : b       = 101.1513
    : bf      = 17.5464
    : Params established from using least squares non linear fitting
    : Used original params +- 100%
    : Used 26 different values for Vhalf
    : Boltzmann slope 7.5
    : cAMP effect on half act: +3 mV @ 1e-5 mM; +10 mV @ 0.01 mM
    : tau tuned to Amendola 2012 fig 7b, Vhalf = half act +4.9 mV (from Amendola)
    : cAMP binding rates left as in Siegelbaum

    :::These are set to tuned values    
    
    a0      = .00032743		(/ms)		: parameters for alpha and beta
    b0      = .00029334		(/ms)
    ac      = -0.1103		(/mV)
    bc      = 0.1025	    	(/mV)

    aa0     = 0.0011    	(/ms)		: parameters for alphaa and betaa
    ba0     = 0.0164 		(/ms)
    aac     = -0.0774		(/mV)
    bac     = 0.1486   		(/mV)

    ::: These are set according to the Vhalf kinetic parameter:
    Vhalf = -90.0           (mV)
    vh1   = 1.057
    vh2   = 79.3
    avh1   = 1.886        
    avh2 = 164.1       

    ::: And the rest are left at default
    kon     = 3.086		    (/mM-ms)	: cyclic AMP binding parameters
    koff    = 4.4857e-05	(/ms)
    b       = 80
    bf      = 8.94

	ai	= 1e-05		(mM)		: concentration cyclic AMP
	gca     = 1				: relative conductance of the bound state
	shift   = 0		(mV)		: shift in voltage dependence
	q10v    = 4				: q10 value from Magee 1998
	q10a    = 1.5				: estimated q10 for the cAMP binding reaction
	celsius			(degC)
}

ASSIGNED {
	v	(mV)
	g	(S/cm2)
	i	(mA/cm2)
	alpha	(/ms)
	beta    (/ms)
	alphaa	(/ms)
	betaa	(/ms)

    vh_shift
    avh_shift

    ah  (mV)
    bh  (mV)
    aah (mV)
    bah (mV)
}

STATE {
	c
	cac
	o
	cao
}

INITIAL {
    setVhalf(Vhalf)
    SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*(o + cao*gca)
	i = g*(v-ehcn)
}

KINETIC kin {
	LOCAL qa
	qa = q10a^((celsius-22 (degC))/10 (degC)) : original
	: qa = q10a^((celsius-35 (degC))/10 (degC))
	rates(v)
	~ c <-> o       (alpha, beta)
	~ c <-> cac     (kon*qa*ai/bf,koff*qa*b/bf)
	~ o <-> cao     (kon*qa*ai, koff*qa)
	~ cac <-> cao   (alphaa, betaa)
	CONSERVE c + cac + o + cao = 1
}

PROCEDURE rates(v(mV)) {
	LOCAL qv
	qv = q10v^((celsius-22 (degC))/10 (degC)) : original
	: qv = q10v^((celsius-37 (degC))/10 (degC))
	if (v > -200) {
		alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))
	} else {
		alpha = a0*qv / (1 + exp(-((-200)-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-((-200)-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-((-200)-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-((-200)-bah-shift)*bac))
	}
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vh_shift = Vhalf*vh1+vh2
    :avh_shift = vh_shift+avh1
    avh_shift = Vhalf*avh1+avh2

    ah      = -87.7 + vh_shift
    bh      = -51.7 + vh_shift
    aah     = -94.2 + avh_shift
    bah     = -35.5 + avh_shift
}

