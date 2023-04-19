  MAD-X:

- KNL[1] and ksl[1] are never communicated to PTC (for bend reasons?)
- If solenoid strength is zero, PTC ignores it and therefore knl and ksl are not communicated to PTC
- In a quadrupole, MAD-X removes any input of knl[2] and ksl[2]
- MAD-X does not communicate tilt in ELSEPARATOR to PTC
- MAD-X does not allow to turn off fringe in RBEND
- MAD-X defines a k0 for the crab cavity based on volt
- MAD-X does not support totalpath in RFMULTIPOLE
- MAD-X does not support totalpath in CRABCAVITY
- PTC seems to ignore tilt in RBEND (MAD-X does not seem to be the culprit.)
- RBEND e1 ~= 0 -> 

```
  The true parallel face bend 
  only accepts the total angle and e1 as an input 
  if e1=0, then the pipe angle to the entrance face is 
  angle/2. It is a normal rbend.
  If e1/=0, then the pipe angle to the entrance face is 
  angle/2+e1 and the exit pipe makes an angle "angle/2-e1" 
  with the exit face.
  The offending non-zero t2 = (e2 - angle/2) is set to zero! 
  Make sure that this is what you want!!! 
```
- RBEND fint and hgap ~= 0 is ignored by MAD-X-PTC
- RBEND k1 is ignored by MAD-X-PTC

  Other:

- MAD-NG does not support the same syntax for the rf cavity as MAD-X for totalpath
- MAD-NG does not support the same syntax for the elseparator as MAD-X for exl and eyl
- Crab cavity is not identical as a rfmultipole, instead k0 and lag are redefined.
- If volt = 0 on an rf multipole, the effects on the multipole (thick only) is not the same as if it went through rfcav_thick map

What needs to be tested for MAD-NG vs MAD-NG:
- General
  - tilt + skew
  - backtracking
  - multipole part vs non multipole part
  - lots of dkd -> tkt
  - negation

- Multipole  
  - Multiple at same point with ksi, ksl, knl seperately and together
