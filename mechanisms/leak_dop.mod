NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
    EXTERNAL apc_metap, fpc_metap
	RANGE i, e, g, gbar
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gbar = 1	(S/cm2)
	e = 0		(mV)
}

ASSIGNED {
	g	(S/cm2)
	i 	(mA/cm2)
	v	(mV)
}

INITIAL {
    g = gbar
	i = g*(v-e)

}

BREAKPOINT {
	g = gbar
	i = g*(v-e)
}

