NEURON {
  SUFFIX NaTa_t
  USEION na WRITE ina READ ena
  RANGE conductance
}

PARAMETER {
  conductance = 0.00001 (uS)
}

STATE { gates_h_q gates_m_q }

INITIAL {
  LOCAL gates_m_reverseRate_x, gates_m_reverseRate_r, gates_m_forwardRate_r, gates_h_reverseRate_x, gates_h_reverseRate_r, gates_h_forwardRate_r

  gates_m_reverseRate_x = 0.16666666666666666 * (38 + v)
  gates_m_reverseRate_r = 0.7440000176429749 * exprelr( gates_m_reverseRate_x)
  gates_m_forwardRate_r = 1.0920000076293945 * exprelr(-gates_m_reverseRate_x)

  gates_h_reverseRate_x = 0.16666666666666666 * (66 + v)
  gates_h_reverseRate_r = 0.09000000357627869 * exprelr(-gates_h_reverseRate_x)
  gates_h_forwardRate_r = 0.09000000357627869 * exprelr( gates_h_reverseRate_x)

  gates_h_q = gates_h_forwardRate_r / (gates_h_forwardRate_r + gates_h_reverseRate_r)
  gates_m_q = gates_m_forwardRate_r / (gates_m_forwardRate_r + gates_m_reverseRate_r)
}

DERIVATIVE dstate {
  LOCAL gates_m_reverseRate_x, gates_m_reverseRate_r, gates_m_forwardRate_r, gates_h_reverseRate_x, gates_h_reverseRate_r, gates_h_forwardRate_r

  gates_m_reverseRate_x = 0.16666666666666666 * (38 + v)
  gates_m_reverseRate_r = 0.7440000176429749 * exprelr( gates_m_reverseRate_x)
  gates_m_forwardRate_r = 1.0920000076293945 * exprelr(-gates_m_reverseRate_x)

  gates_h_reverseRate_x = 0.16666666666666666 * (66 + v)
  gates_h_reverseRate_r = 0.09000000357627869 * exprelr(-gates_h_reverseRate_x)
  gates_h_forwardRate_r = 0.09000000357627869 * exprelr( gates_h_reverseRate_x)

  gates_h_q' = (gates_h_forwardRate_r - (gates_h_forwardRate_r + gates_h_reverseRate_r) * gates_h_q) * 2.9528825283050537
  gates_m_q' = (gates_m_forwardRate_r - (gates_m_forwardRate_r + gates_m_reverseRate_r) * gates_m_q) * 2.9528825283050537
}

BREAKPOINT {
  SOLVE dstate METHOD cnexp
  LOCAL g

  g = conductance * gates_h_q * gates_m_q * gates_m_q * gates_m_q
  ina = g * (v - ena)
}

