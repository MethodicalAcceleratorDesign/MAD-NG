return DEFAULTS {
  el_type = "sextupole",
  el_args = T[[k2:=${k2}, k2s:=${k2s}, l=1, at=0]],
  k2      = 0.01,
  k2s     = 0.00,
  makethin = 'makethin, sequence=seq;',
  -- output:
  prefix = 'sextupole/',
  studies = DEFAULTS.studies {
    k2  = {stop=0.1},
    k2s = {stop=0.1, k2=0.0},
  },
}
