local DEFAULTS = R 'defaults'
return DEFAULTS {
  prefix = 'drift/',
  el_type = "drift",
  el_args = [[l=1, at=0.5]],
  sad_el_type = 'DRIFT',
  sad_el_args = 'L=1',
  studies = DEFAULTS.studies {
    x_px = {px=0.01, x=T"${varname}"},
  },
}
