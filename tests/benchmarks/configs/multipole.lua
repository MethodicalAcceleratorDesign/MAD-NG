local DEFAULTS = R 'thinelement'
return DEFAULTS {
  prefix = 'multipole/',
  el_type = "multipole",
  el_args = T'knl:={${k0l}, ${k1l}, ${k2l}}, at=${at}',
  sad_el_type = 'MULT',
  sad_el_args = T'K0=${k0l} K1=${k1l} K2=${k2l} L=0',
  k0l     = 0,
  k1l     = 0,
  k2l     = 0,
  -- output:
  studies = DEFAULTS.studies {
    k0l = {stop=0.01, px=0.001, sad_attr='K0'},
    k1l = {stop=0.01, px=0.001, sad_attr='K1'},
    k2l = {stop=0.01, px=0.001, sad_attr='K2'},
  },
}
