return DEFAULTS {
  el_type = "multipole",
  el_args = [[knl:={k0l, k1l, k2l}, at=0.5]],
  k0l     = 0,
  k1l     = 0,
  k2l     = 0,
  -- output:
  prefix = 'multipole/',
  studies = DEFAULTS.studies {
    k0l = {stop=0.01, px=0.001},
    k1l = {stop=0.01, px=0.001},
    k2l = {stop=0.01, px=0.001},
  },
}
