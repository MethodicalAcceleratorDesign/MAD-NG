local DEFAULTS = R 'defaults'
return DEFAULTS {
  prefix = 'sbend/',
  el_type = "sbend",
  el_args = T"\
angle:=${angle}, tilt:=${tilt}, \
k0:=${k0}, k1:=${k1}, k2:=${k2}, \
e1:=${e1}, e2:=${e2}, \
h1:=${h1}, h2:=${h2}, \
fint=${fint}, fintx=${fintx}, \
hgap=${hgap}, \
l=${el_len}, at=0.5",
  sad_el_type = 'BEND',
  sad_el_args = T'ANGLE=${angle} L=${el_len} ROTATE=${tilt} E1=${e1} E2=${e2}',
  el_len  = 1,
  -- cannot use ANGLE=0 for MAD-NG:
  angle   = 0.1,
  tilt    = 0.0,
  k0      = T"${angle}/${el_len}",
  k1 = 0, k2 = 0,
  e1 = 0, e2 = 0,
  h1 = 0, h2 = 0,
  fint    = 0,
  fintx   = 0,
  hgap    = 0,
  -- engines = DEFAULTS.engines { thick = false },
  -- output:
  studies = DEFAULTS.studies {
    -- cannot use ANGLE=0 for MAD-NG:
    angle = {stop = "pi/4", start=0.01, sad_attr='ANGLE'  },
    tilt  = {stop = "pi/4",             sad_attr='ROTATE' },
    e1    = {stop = "pi/6",             sad_attr='E1'     },
    e2    = {stop = "pi/6",             sad_attr='E2'     },
    k1    = {                           sad_attr='K1'     },
    -- MISSING (for now):
    -- k0 : ignored in MAD-X
    -- k2  : not in SAD
    -- hgap: not in SAD?
    -- h1, h2: ?
    -- F1 ~ fint=fintx?
  },
}
