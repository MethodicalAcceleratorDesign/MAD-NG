Building a test for an element
------------------------------

To create a test for a single element, there is a single component that is required; an object containing the configuration of the test. This object must hold the following information: 

* A string, named `elm`, containing the MAD-X element definition, with the parameters to be tested as part of the string, ready to be formatted with the configuration table, see below.
* An array part, named `alist` for attribute list, within the object, containing the names of the parameters to be tested.
* The hash part of the object, containing an *iterable* for each parameter, see below.

Since the array part and the hash part also need to contain the reference parameters, the configuration object specific to your test must inherit from ref_cfg, which is a table that contains the reference parameters. To take from reference `alist`, you must use the `tbl_cat` function, see below.  

The sequence is a string that contains the MAD-X script to create the element that you desire to test, with all the necessary parameters, the parameters can be assigned through formatting with the configuration table. For example, a quadrupole with the parameters `k1`, `k1s` `fringe` and `tilt` can be created with the following string: 
```
"QUADRUPOLE, at=0.75, l=1.5, k1=${k1}, k1s=${k1s}, tilt=${tilt}*pi/8, fringe=${fringe}"
```

Below is an example of a configuration object for a quadrupole, with the parameters `k1`, `k1s` `fringe` and `tilt`:
```
local cfg = ref_cfg "quadrupole" {
    elm = " QUADRUPOLE, at=0.75, l=1.5, k1=${k1}, k1s=${k1s}, tilt=${tilt}*pi/8, fringe=${fringe}",
    tol = 1000,

    alist = tblcat(ref_cfg.alist, {"tilt", "fringe", "k1", "k1s"}),
    tilt   = 0  ..4,
    fringe = 0  ..3 ..3,
    k1     = -0.2..0.2..0.2,
    k1s    = -0.2..0.2..0.2,
}
run_test(cfg)
```

Above shows the required parts of the object, the element string, tolerance, an extended attribute list and the parameters with their ranges. The tolerance is the number that the difference between the maps is compared to, it is multiplied by your machine epsilon (`eps`) to get the number to compare to. The tolerance is optional as the default is 1000. If the tolerance is set to a string then the test will try to read the tolerance from a file with the same name as the tolerance string

Finally, the test is run by calling ``run_test(cfg)``. 

How the tests work
------------------

The entry function to running the tests is inventively called `run_test`. This function takes a configuration object as the only argument. Within this function the following steps are taken:

1. Checks if cfg.dorun is not nil or false, if it is, then the function prints the previous results from the previous run and returns.
2. An mtable with 2 columns is created, the first column will store the configuration snapshots, while the second column will store the results of each configuration snapshot.
3. The configuration snapshots are created using a recursive function that takes the configuration object and creates the snapshots and runs the function `run_cfg` on each snapshot.  
    The function `run_cfg` takes a configuration snapshot and does the following:
    1. Calls `do_trck`, which runs track and ptc to get the difference between MAD and PTC, then stores the results in a table and returns it.
    2. Stores the configuration snapshot and the results into the mtable.
    3. If cfg.doprnt then the results are printed to the console in a single line per config.
    4. If cfg.dodbg then the results are compared to the tolerance, if the results are outside the tolerance then the test stops and the results are printed to the console.
4. Now the mtable has been filled with the results of each configuration snapshot, we create the generator columns for the mtable, which link the tables of each configuration snapshot and the results to a column in the mtable.
5. If cfg.dosave or the file does not exist then the mtable is saved to a file with the name of the element appended with `_cfg` or `_res` for the configuration snapshots and the results respectively.
6. If cfg.doprnt then a summary of the results is printed to the console. (See below for an example of the summary)
7. If the test has not been stopped by cfg.dodbg then several excess files are removed.


Running the tests
-----------------

To run the test, just call the element function For example, ``testQUAD()`` will run the test `testQUAD` (the quadrupole).

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

With the `doprnt=true`, you will also receive a final output. This output is the sum of all the failures for each order, for each column of your config. The idea of this sum is to give you an idea of what is causing the failure. For example, see the output below for a quadrupole (truncated at order 0), run on linux with ``testQUAD()`` as the input command.

```
cfg = {
    elm = " QUADRUPOLE, at=0.75, l=1.5, k1=${k1}, k1s=${k1s}, tilt=${tilt}*pi/8, fringe=${fringe}",
    tol = 1000,
    
    alist = {"model", "energy", "method", "nslice", "x0i", "order", "icase", "tilt", "fringe", "k1", "k1s",},
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