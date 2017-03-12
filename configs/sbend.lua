local DEFAULTS = R 'defaults'
return DEFAULTS {
  el_type = "sbend",
  el_args = T[[angle:=${angle}, k0:=${k0}, tilt:=${tilt}, l=${el_len}, at=0]],
  angle   = 0.1,
  tilt    = 0.0,
  k0      = T"${angle}/${el_len}",
  el_len  = 1,
  -- output:
  prefix = 'sbend/',
  studies = DEFAULTS.studies {
    angle = {stop = "pi/4", start=0.01},
    tilt  = {stop = "pi/4"},
  },
}
