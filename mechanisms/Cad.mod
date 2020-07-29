: Calcium handler using parameters from Table 1 in
: Wilson and Callaway, 2000, J Neurophysiol

NEURON {
	SUFFIX cad
	USEION ca READ cai, ica WRITE cai
	GLOBAL Pmax, beta
	RANGE cai_prime, CaCurr, CaDep
    EXTERNAL apc_metap, fpc_metap
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (milli/liter)
	(um) = (micrometer)
}

PARAMETER {
	Pmax = 2	(um/ms)
	beta = 0.001
	z = 2
	F = 9.64846e4	(coulomb)
}

ASSIGNED {
	diam		(um)
	: cai			(mM)
	ica			(mA/cm2)
	icalcium	(mA/m2)
	dt			(ms)
	cai_prime	(mM/ms)
    CaCurr      (mM/ms)
    CaDep       (mM/ms)
}

STATE {
    cai (mM) : intracellular calcium
}

BREAKPOINT {
    SOLVE CaUpdate METHOD cnexp
}

DERIVATIVE CaUpdate {
    rates(ica, cai)
    cai' = CaCurr - CaDep
}

INITIAL {
    cai = 5e-5
}

PROCEDURE rates(ica (mA/cm2), cai (mM)) {
	icalcium = ica * -1e4 (cm2/m2)
    CaCurr = (icalcium*4*beta)/(z*F*diam)
    CaDep = (Pmax*4*beta*cai)/diam
	: cai_prime = ((icalcium*4*beta)/(z*F*diam)) - ((Pmax*4*beta*cai)/diam)
}

