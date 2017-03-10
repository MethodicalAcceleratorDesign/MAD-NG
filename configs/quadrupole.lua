return DEFAULTS {
  el_type = "quadrupole",
  el_args = T[[k1:=${k1}, k1s:=${k1s}, tilt=${tilt}, l=1, at=0]],
  k1      = 0.01,
  k1s     = 0.0,
  tilt    = 0.0,
  -- output:
  prefix = 'quadrupole/',
  studies = DEFAULTS.studies {
    k1   = {stop = 0.01, k1s = 0},
    k1s  = {stop = 0.01, k1 = 0},
    tilt = {stop = "pi/4"},
    -- TODO: l
  },
}
