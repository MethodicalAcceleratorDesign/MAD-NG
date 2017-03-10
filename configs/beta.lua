return DEFAULTS {
  el_type = "drift",
  el_args = [[l=1, at=0]],
  prefix = 'beta/',
  studies = MAD.Object {
    energy = {start=1, stop=2000, count=50, varfunc="logrange"},
  },
}
