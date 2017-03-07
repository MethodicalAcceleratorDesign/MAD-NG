return DEFAULTS {
  el_type = "drift",
  el_args = [[l=1, at=0]],
  prefix = 'drift/',
  studies = DEFAULTS.studies {
    x_px = {px=0.01, x="${varname}"},
  },
}

