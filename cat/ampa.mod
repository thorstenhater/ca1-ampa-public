COMMENT
Approximate AMPA/NMDA/GABA synapse model after NEST

References
==========

* Probabilistic Decision Making by Slow Reverberation in Cortical Circuits (Wang 2002)
* A flexible and precise approximation to the Wang (2002) NMDA Model (Skaar and Plesser 2024)
  https://nest-simulator.readthedocs.io/en/master/model_details/Brunel_Wang_2001_Model_Approximation.html

ENDCOMMENT

NEURON {
    POINT_PROCESS ampa
    NONSPECIFIC_CURRENT i
    RANGE k0, k1, e_ex, e_in, gbar_ampa, gbar_gaba, gbar_nmda, mg, alpha, tau_gaba, tau_ampa, tau_nmda, tau_nmda_rise
}

STATE { s_ampa s_gaba s_nmda }

PARAMETER {
  : builtins to read
  v

  : decay times
  tau_gaba      =   5.0      (ms)
  tau_ampa      =   2.0      (ms)
  tau_nmda      = 100.0      (ms)

  : NMDA paramters
  alpha         =   0.5      (kHz)
  tau_nmda_rise =   2.0      (ms)

  : maximum conductances from Wang 2002 Experimental Procedures for Pyramidal cells
  gbar_gaba     =   0.0013   (uS)
  gbar_ampa     =   0.0021   (uS)
  gbar_nmda     =   0.000165 (uS)

  : external Mg concentration
  mg            =   1.0      (mM)

  : reversal potentials
  e_ex          =   0.0      (mV)
  e_in          = -70.0      (mV)
}

ASSIGNED { k0 k1 }

INITIAL {
  LOCAL at, rt

  s_ampa = 0
  s_gaba = 0
  s_nmda = 0

  at = alpha*tau_nmda_rise
  rt = tau_nmda_rise/tau_nmda
  k0 = at^rt * linc_gamma(1 - rt, at)
  k1 = exp(-at) - 1
}

BREAKPOINT {
  SOLVE dS METHOD cnexp

  LOCAL i_ampa, g_ampa, i_gaba, g_gaba, i_nmda, g_nmda

  : AMPA current
  g_ampa = gbar_ampa * s_ampa
  i_ampa = g_ampa * (v - e_ex)

  : NMDA current
  g_nmda = gbar_nmda * s_nmda / (1 + mg/3.57*exp(-0.062*v))
  i_nmda = g_nmda * (v - e_ex)

  : GABA current
  g_gaba = gbar_gaba * s_gaba
  i_gaba = g_gaba * (v - e_in)

  : i (nA)
  i = i_ampa + i_gaba + i_nmda
}

DERIVATIVE dS {
  s_ampa' = -s_ampa/tau_ampa
  s_gaba' = -s_gaba/tau_gaba
  s_nmda' = -s_nmda/tau_nmda
}

NET_RECEIVE(weight) {
  if (weight > 0) {
     s_gaba = s_gaba + weight
     s_ampa = s_ampa + weight
  }
  else if (weight < 0) {
    s_nmda = s_nmda + weight*(k0 + k1*s_nmda)
  }
}

: Approximate series upto four terms. Might be bad!
: Based on:
:                 --- oo           s + k
:                 \          k    z
:    gamma(s, z) = >     (-1)   ------------
:                 /              k! (s + k)
:                 --- k = 0
: https://en.wikipedia.org/wiki/Incomplete_gamma_function#Evaluation_formulae
FUNCTION linc_gamma(s, z) { linc_gamma = (z^s)*(1/s - z/(s + 1) + z^2/(2*(s + 2)) - z^3/(6*(s + 3))) }
