local DEFAULTS = R 'defaults'
return DEFAULTS {
  prefix = 'sextupole/',
  el_type = "sextupole",
  el_args = T[[k2:=${k2}, k2s:=${k2s}, l=1, at=0.5]],
  sad_el_type = 'SEXT',
  sad_el_args = T'L=1 K2=${k2} ROTATE=${tilt}',
  k2      = 0.01,
  k2s     = 0.00,
  tilt    = 0,
  nst = 3,
  engines = DEFAULTS.engines { thick = false },
  -- output:
  studies = DEFAULTS.studies {
    k2  = {stop=0.1,         sad_attr = 'K2', tilt=0     },
    k2s = {stop=0.1, k2=0.0, sad_attr = 'K2', tilt='pi/4'},
  },
}
