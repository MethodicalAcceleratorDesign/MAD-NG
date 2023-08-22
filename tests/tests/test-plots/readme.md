Run the following to run everything successfully:
```bash
../mad test-all-plots.mad -s --ptc --ng --cpp --conv -x Deca -x Elsep
../mad test-all-plots.mad -s --ng --cpp --conv -p Deca
../mad test-all-plots.mad -s --ng --cpp --conv -p Elsep
```

Reasons for the three commands:
- The first command runs everything but:
    - Decapole
    - Dodecapole
    - Electrostatic separators
- The second command runs the decapole and dodecapole, without PTC
- The third command runs the electrostatic separators, without PTC 

For the decapole and dodecapole, the elements are not implemented in PTC.
For electrostatic separators, PTC does include any fringes or slice the element, so (dbgmap fails).

