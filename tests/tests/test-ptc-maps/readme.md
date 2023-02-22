Building a test for an element
------------------------------

To create a test for a single element, there are two main components that need to be created; the sequence and the configuration specific to the element.

The sequence is a string that contains the MAD-X script to create the element that you desire to test, with all the necessary parameters, the parameters can be assigned through formatting with the configuration table. For example, a quadrupole with the parameters `k1`, `k1s` `fringe` and `tilt` can be created with the following string: 
```
local elm_str = "QUADRUPOLE, at=0.75, l=1.5, k1=${k1}, k1s=${k1s}, tilt=${tilt}*pi/8, fringe=${fringe}"
```

The configuration table is a table that contains the names of the parameters and therefore keys in the array part of the table, and the values attached to the parameter name in the hash part of the table, see the example below;

```
quad_cfg = {
    "tilt", "fringe", "k1", "k1s",
    tilt   = 0  ..4,
    fringe = 0  ..3 ..3,
    k1     = -0.2..0.2..0.2,
    k1s    = -0.2..0.2..0.2,
}
```

Then this can be combined with the reference table using ``tbl_cat(ref_cfg, quad_cfg)`` to create a table that contains all the parameters for the element. The reference table contains the parameters that are common to all elements, such as the energy, the number of slices, the initial conditions, etc. To change the reference table, you must copy the table using ``tblcopy`` and then change the values in the copy, as the reference table is used by all the tests. 

Finally, the test is run by calling ``run_test(element, elm_str, cfg, tol_)``, where `element` is the element object that would like to be tested, `elm_str` is the string that contains the MAD-X script to create the element, `cfg` is the total configuration table, and `tol_` is the (optional) tolerance for the test. This test will first create the config mtable then loop through the config, creating a file called ``element.seq`` that is called by MAD-X and MAD-NG. The final map is then gathered from the mflow in MAD-NG and from the debug output of MAD-X-PTC to be compared and stored into a max dif mtable, which is then saved to a file called ``element_name_max_err.tfs`` at the end of the test.

Running the tests
-----------------

To run the test, just call the element function with a optionally a tolerance as the only argument. The tolerance should be an integer, as the number that the difference between the maps is compared to is this tolerance multiplied by your machine epsilon (`eps`). 

For example, ``testQUAD(1000)`` will run the test `testQUAD` (the quadrupole) with a tolerance of `1000*eps`.

There are four flags that can be used to control how your tests run and the reported results:

* `doprnt` - If set to `true`, the test will print the differences between ptc and MAD-NG to the console for each order at each config.
* `dorun` - If set to `true`, the test will run the test, rather than just reading the results from the file and outputting them.
* `dosave` - If set to `true`, the test will save the results and the config to the files, even if a file with the same name already exists.
* `dodbg` - If set to `true`, the test will check the result against a tolerance, given by the `tol` argument, if a failure occurs, the test will stop, create files called `ref.mad`, `element.seq`, `element_name_n.txt`, `element_name_p.txt`, `element_name_d.txt`. Then to run the same scenario again with the failing config, you can run a MAD-NG script (below) where `element_name` is the name of the element you are testing. Also if `dosave==true`, the test will save the results and the config to the files, with the results stopping at the failed config.

    >```
    >os.execute("../mad  ref.mad  > element_name_n.txt")
    >os.execute("../madx ref.madx > element_name_p.txt")
    >require "madl_dbgmap".cmpmdump("element_name")
    >```

With the `doprnt=true`, you will also receive a final output. This output is the sum of all the failures for each order, for each column of your config. The idea of this sum is to give you an idea of what is causing the failure. For example, see the output below for a quadrupole (truncated at order 0), run on linux with ``testQUAD(1000)`` as the input command.

```
cfg = {
    "model", "energy", "method", "nslice", "x0i", "order", "icase", "tilt", "fringe", "k1", "k1s",
    model = 1..2,
    energy = {1, 6500}, 
    method = 2..2..2,
    nslice = 1..3..1, 
    x0i = 1..4,
    order = {2},
    icase = {56},
    tilt   = 0  ..4,
    fringe = 0  ..3 ..3,
    k1     = -0.2..0.2..0.2,
    k1s    = -0.2..0.2..0.2,
}

order 0:
nslice  = {360, 360, 360}
x0i     = {[2]=360, [3]=360, [4]=360}
energy  = {540, [6500]=540}
icase   = {[56]=1080}
k1      = {[0]=360, [-0.2]=360, [0.2]=360}
k1s     = {[-0.2]=540, [0.2]=540}
fringe  = {[3]=1080}
tilt    = {216, 216, 216, 216, [0]=216}
method  = {[2]=1080}
model   = {540, 540}
order   = {[2]=1080}
```

In this list, you can see situations where the failure never occured and situations where changing the variable had no effect. So the variables that had no effect here are ``nslice, energy, icase, k1, tilt, method, model, order``. While the specifics that are missing are `x0i = 1`, `k1s = 0`, `fringe = 0`, indicating the issue only occurs when the initial conditions are `x0i > 1` (not the zero orbit) and when there is a fringe field, and a `k1s` value.

This can then be further investigated by setting ``dodbg`` to `true`, and running the test again, then using the MAD-NG script above to run the failing config. 