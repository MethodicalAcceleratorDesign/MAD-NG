local DEFAULTS = R 'defaults'
return DEFAULTS {
  prefix = 'multipole/',
  el_type = "multipole",
  el_args = T[[knl:={${k0l}, ${k1l}, ${k2l}}, at=0.5]],
  sad_el_type = 'MULT',
  sad_el_args = T'L=1 K0=${k0l} K1=${k1l} K2=${k2l}',
  mad_sequence = T[[
local sequence = MAD.sequence "${seq_name}" {
  refer="entry", l=${seq_len},
  MAD.element.drift "dr" {l=0.5, at=0},
  MAD.element.${el_type} "${el_name}" { ${el_args} },
}]],
  k0l     = 0,
  k1l     = 0,
  k2l     = 0,
  observe = '#e',
  -- output:
  studies = DEFAULTS.studies {
    k0l = {stop=0.01, px=0.001, sad_attr='K0'},
    k1l = {stop=0.01, px=0.001, sad_attr='K1'},
    k2l = {stop=0.01, px=0.001, sad_attr='K2'},
  },
}
