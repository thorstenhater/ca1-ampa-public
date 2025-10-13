NEURON {
  SUFFIX Ca_HVA
  USEION ca WRITE ica READ eca
  RANGE conductance
}

PARAMETER {
  conductance = 0.00001 (uS)
}

STATE { gates_h_q gates_m_q }

INITIAL {
  LOCAL gates_m_reverseRate_r, gates_m_forwardRate_x, gates_m_forwardRate_r, gates_m_inf, gates_h_forwardRate_r, gates_h_reverseRate_r, gates_h_inf

  gates_m_reverseRate_r = 0.9399999976158142 * exp(-0.058823529411764705 * (75 + v))
  gates_m_forwardRate_x = 0.26315789803903855 * (27 + v)
  gates_m_forwardRate_r = 0.20900000631809235 * exprelr(-gates_m_forwardRate_x)

  gates_h_forwardRate_r = 0.0004569999873638153 * exp(-0.02 * (13 + v))
  gates_h_reverseRate_r = 0.006500000134110451 / (1 + exp(-0.03571428571428571 * (15 + v)))

  gates_h_q = gates_h_forwardRate_r / (gates_h_forwardRate_r + gates_h_reverseRate_r)
  gates_m_q = gates_m_forwardRate_r / (gates_m_forwardRate_r + gates_m_reverseRate_r)
}

DERIVATIVE dstate {
  LOCAL gates_m_reverseRate_r, gates_m_forwardRate_x, gates_m_forwardRate_r, gates_m_inf, gates_m_tau, gates_h_forwardRate_r, gates_h_reverseRate_r, gates_h_inf, gates_h_tau

  gates_m_reverseRate_r = 0.9399999976158142 * exp(-0.058823529411764705 * (75 + v))
  gates_m_forwardRate_x = 0.26315789803903855 * (27 + v)
  gates_m_forwardRate_r = 0.20900000631809235 * exprelr(-gates_m_forwardRate_x)

  gates_h_forwardRate_r = 0.0004569999873638153 * exp(-0.02 * (13 + v))
  gates_h_reverseRate_r = 0.006500000134110451 / (1 + exp(-0.03571428571428571 * (15 + v)))

  gates_h_q' = gates_h_forwardRate_r - gates_h_q * (gates_h_forwardRate_r + gates_h_reverseRate_r)
  gates_m_q' = gates_m_forwardRate_r - gates_m_q * (gates_m_forwardRate_r + gates_m_reverseRate_r)
}

BREAKPOINT {
  SOLVE dstate METHOD cnexp
  LOCAL g
  g = conductance * gates_h_q * gates_m_q * gates_m_q
  ica = g * (v -eca)
}

