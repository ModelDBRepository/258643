TITLE SK

COMMENT
    SK small conductance potassium channel
    Implementation taken from:
        Forrest MD (2014) Two Compartment model of the cerebellar purkinje neuron [modelDB 180789]
    Purely calcium dependent
    km is half-activation calcium concentration
    n is Hill number
    Edited by Tim Rumbell, IBM Research, 2019, thrumbel@us.ibm.com
ENDCOMMENT

NEURON {
    SUFFIX sk
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gbar, oinf, n, ik
    GLOBAL km:, taumod
    EXTERNAL apc_metap, fpc_metap
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)    
    (molar) = (1/liter)
    (mM) = (millimolar)
    (S) = (siemens)
}

PARAMETER {
    gbar = 1e-05    (S/cm2)
    n = 4 : from many previous DA models :4.2 : from Xia 98 (including calmodulin)
    km = 0.0002 (mM) : many previous DA models between 180 and 250 nM : 0.00035 (mM) : from Xia 98 (incl. calmodulin)
    cai (mM)
    dt  (ms)
    ek  (mV)

    : to = 15 (ms)    
    : taumod = 1

    v   (mV)
}
ASSIGNED {
    ik  (mA/cm2)
    oinf
}

BREAKPOINT {
    oinf = 1/(1 + pow(km/cai,n))
    ik = oinf*gbar*(v-ek)
}

