local DEFAULTS = R 'sbend'
return DEFAULTS {
  prefix = 'sbend_nosad/',
  engines = DEFAULTS.engines { sad = false },
  -- output:
  studies = {
    -- k0 : ignored in MAD-X
    k1    = {                           },
    k2    = {                           },
    h1    = {                           },
    h2    = {                           },
    fint  = {                           },
    fintx = {                           },
    hgap  = {                           },
  },
}
