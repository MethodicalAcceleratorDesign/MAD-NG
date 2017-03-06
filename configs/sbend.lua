return DEFAULTS {
  el_type = "sbend",
  el_args = [[angle:=${angle}, l=1, at=0]],
  angle   = 0.1,
  -- output:
  prefix = 'sbend/',
  studies = DEFAULTS.studies {
    angle = {stop=1.5707963267948966},
  },
}
