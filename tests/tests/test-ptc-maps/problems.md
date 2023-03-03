  MAD-X:
    - KNL[1] and ksl[1] are never communicated to PTC (for bend reasons?)
    - If solenoid strength is zero, PTC ignores it and therefore knl and ksl are not communicated to PTC
    - In a quadrupole, MAD-X removes any input of knl[2] and ksl[2]
    - MAD-X does not communicate tilt in ELSEPARATOR to PTC
    - MAD-X does not allow to turn off fringe in RBEND
    - MAD-X defines a k0 for the crab cavity based on volt
    - MAD-X does not support totalpath in RFMULTIPOLE
    - MAD-X does not support totalpath in CRABCAVITY
  
  Other:
    - MAD-NG does not support the same syntax for the rf cavity as MAD-X for totalpath
    - MAD-NG does not support the same syntax for the elseparator as MAD-X for exl and eyl