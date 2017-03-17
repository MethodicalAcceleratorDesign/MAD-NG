local DEFAULTS = R 'defaults'
return DEFAULTS {
  prefix = 'sbend/',
  el_type = "sbend",
  el_args = T[[angle:=${angle}, k0:=${k0}, tilt:=${tilt}, l=${el_len}, at=0.5]],
  sad_el_type = 'BEND',
  sad_el_args = T'ANGLE=${angle} L=${el_len} ROTATE=${tilt}',
  angle   = 0.1,
  tilt    = 0.0,
  k0      = T"${angle}/${el_len}",
  el_len  = 1,
  -- engines = DEFAULTS.engines { thick = false },
  -- output:
  studies = DEFAULTS.studies {
    angle = {stop = "pi/4", start=0.01, sad_attr='ANGLE' },
    tilt  = {stop = "pi/4",             sad_attr='ROTATE'},
  },
}
