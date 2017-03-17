-- TODO: consistency check with SBEND: E1 E2
local DEFAULTS = R 'defaults'
return DEFAULTS {
  prefix = 'dipedge_e1/',
  el_type = "dipedge",

  at = 0,
  sb_len = 1,

  madx_sequence = T[[
${seq_name}: sequence, refer=center, l=${seq_len};
  ${el_name}: DIPEDGE, ${el_args};
  BND: SBEND, ${sb_args};
  MOBS: MARKER, AT=${at_exit};
endsequence;]],
  observe = 'MOBS',
  el_args = T'HGAP=${HGAP}, TILT=${tilt}, H=${h}, E1=${e1}, FINT=${fint}, AT=${at}',
  sb_args = T'HGAP=${hgap}, TILT=${tilt}, ANGLE:=${angle}, L=${sb_len}, AT=${at}+${sb_len}/2',
  h       = T"${angle}/${sb_len}",
  el_len  = 1,
  angle   = 0.1,
  tilt    = 0.0,
  e1      = 0,
  fint    = 0,
  engines = DEFAULTS.engines { sad = false, mad = false },
  -- output:
  studies = DEFAULTS.studies {
    angle = {stop = "pi/4", start=0.01, },
    tilt  = {stop = "pi/4",             },
    e1    = {stop = "pi/6",             },
    fint  = {                           },
  },
}
