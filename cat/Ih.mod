NEURON {
  SUFFIX Ih
  NONSPECIFIC_CURRENT ihcn
  RANGE conductance
  GLOBAL ehcn
}

PARAMETER {
  conductance = 0.00001 (uS)
  ehcn = -45 (mV)
}

STATE { gates_m_q }

INITIAL {
  LOCAL gates_m_reverseRate_r, gates_m_forwardRate_x, gates_m_forwardRate_r

  gates_m_reverseRate_r = 0.19300000369548798 * exp(0.030211481755258694 * v)
  gates_m_forwardRate_x = -0.08403361613918324 * (154.89999389648438 + v)
  gates_m_forwardRate_r = 0.0765170007944107 * exprelr(-gates_m_forwardRate_x)

  gates_m_q =  gates_m_forwardRate_r / (gates_m_forwardRate_r + gates_m_reverseRate_r)
}

DERIVATIVE dstate {
  LOCAL gates_m_reverseRate_r, gates_m_forwardRate_x, gates_m_forwardRate_r

  gates_m_reverseRate_r = 0.19300000369548798 * exp(0.030211481755258694 * v)
  gates_m_forwardRate_x = -0.08403361613918324 * (154.89999389648438 + v)
  gates_m_forwardRate_r = 0.0765170007944107 * exprelr(-gates_m_forwardRate_x)

  gates_m_q' = gates_m_forwardRate_r - gates_m_q * (gates_m_forwardRate_r + gates_m_reverseRate_r)
}

BREAKPOINT {
  SOLVE dstate METHOD cnexp
  LOCAL g

  g = conductance * gates_m_q
  ihcn = g * (v -ehcn)
}

