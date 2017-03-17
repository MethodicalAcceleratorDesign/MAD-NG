local DEFAULTS = R 'thinelement'
return DEFAULTS {
  prefix = 'hkicker/',
  el_type = "hkicker",
  el_args = T'kick:=${kick}, at=${at}, tilt:=${tilt}',
  sad_el_type = 'BEND',
  sad_el_args = T'ANGLE=0 L=0 K0=${kick} ROTATE=${tilt}',
  kick    = 0.1,
  tilt    = 0.0,
  -- output:
  studies = DEFAULTS.studies {
    kick  = {stop = "pi/4", sad_attr='K0' },
    tilt  = {stop = "pi/4", sad_attr='ROTATE'},
  },
}
