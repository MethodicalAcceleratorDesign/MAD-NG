local DEFAULTS = R 'defaults'
return DEFAULTS {
  mad_sequence = T[[
local sequence = MAD.sequence "${seq_name}" {
  refer="entry", l=${seq_len},
  MAD.element.drift "dr" {l=0.5, at=0.25},
  MAD.element.${el_type} "${el_name}" { ${el_args} },
}]],
  el_type = "multipole",
  el_args = [[knl:={k0l, k1l, k2l}, at=0.5]],
  k0l     = 0,
  k1l     = 0,
  k2l     = 0,
  observe = '#e',
  -- output:
  prefix = 'multipole/',
  studies = DEFAULTS.studies {
    k0l = {stop=0.01, px=0.001},
    k1l = {stop=0.01, px=0.001},
    k2l = {stop=0.01, px=0.001},
  },
}
