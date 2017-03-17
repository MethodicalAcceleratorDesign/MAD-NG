local DEFAULTS = R 'defaults'
return DEFAULTS {
  at      = 0.5,
  at_exit = 1.0,
  madx_sequence = T[[
${seq_name}: sequence, refer=center, l=${seq_len};
  ${el_name}: ${el_type}, ${el_args};
  MOBS: MARKER, AT=${at_exit};
endsequence;]],
  mad_sequence = T[[
local sequence = MAD.sequence "${seq_name}" {
  refer="entry", l=${seq_len},
  MAD.element.drift "dr0" {l=${at}},
  MAD.element.${el_type} "${el_name}" { ${el_args} },
  MAD.element.marker "MOBS" {at=${at_exit}},
}]],
  sad_sequence = T[[
MARK M1 = ();
DRIFT DR0 = (L=${at});
DRIFT DR1 = (L=${at_exit}-${at});
${sad_el_type} ${el_name} = (${sad_el_args});
MARK MOBS = ();
LINE ${seq_name} = (M1 DR0 EL DR1 MOBS);]],
  observe = 'MOBS',
}
