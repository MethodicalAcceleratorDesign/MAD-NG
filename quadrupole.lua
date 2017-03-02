return DEFAULTS {
  el_type = "quadrupole",
  el_args = [[k1:=${k1}, l=1, at=0]],
  k1      = 0.1,
  -- output:
  prefix = 'quadrupole/',
  studies = merge(DEFAULTS.studies, {
    k1 = {stop=0.01, x=0.01, y=0.01},
  }),
}
