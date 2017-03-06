return DEFAULTS {
  el_type = "multipole",
  el_args = [[knl:={k0l_, k1l_, k2l_}, at=0.5]],
  k0l     = 0,
  k1l     = 0,
  k2l     = 0,
  -- output:
  prefix = 'multipole/',
  studies = DEFAULTS.studies {
    k0l_ = {stop=0.01},
    k1l_ = {stop=0.01},
    k2l_ = {stop=0.01},
  },
}
