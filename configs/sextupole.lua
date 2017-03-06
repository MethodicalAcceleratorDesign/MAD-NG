return DEFAULTS {
  el_type = "sextupole",
  el_args = [[k2:=${k2}, l=1, at=0]],
  k2      = 0.01,
  -- output:
  prefix = 'sextupole/',
  studies = DEFAULTS.studies {
    k2 = {stop=0.1, x=0.01, y=0.01},
  },
}
