# Test Suite for PTC vs NG maps

## Purpose

The purpose of this test suite is to test MAD-NG maps against PTC. 

## Prerequisites

A binary of MAD-NG and MAD-X must be in the parent directory of the test directory. The binaries must be named `mad` and `madx` respectively.

## Commands to run the tests

To run all the tests:

```shell
../mad test-all-maps.mad -t
```
or 
```shell
../mad test-all-maps.mad --test
```

To run a category of tests, like in unit tests:

```shell
../mad test-all-maps.mad --test TestModule
```

To run an individual test: 

```shell
../mad test-category-maps.mad --test TestModule.TestName
```

If you would like to run a unit test, you must ensure that you have generated the files for the unit tests (these take a bit of time, could be more than 3 hours for them all), to generate these run the following, this can be done while running the test as well:
```shell 
../mad test-all-maps.mad -g
```

Once the files have all been generated, now you can run unit tests. To run a unit test, it is identical to how the unit tests in the folder utest are called:
```shell
../mad test-all-maps.mad TestModule.TestName
```
where TestName is optional, if it is not provided, then all the tests in the module will be run. Also TestModule is optional, if it is not provided, then all the tests in the folder will be run.

## Details on Running the Tests

To run the tests you must put the first argument after the file name as `-test` or the file will assume that you are running the file as a unit test. 


There are five flags that can be used to control how your tests run and the reported results:

* `doprnt` - If set to `true`, the test will print the differences between PTC and MAD-NG to the console for each order at each config (default: `false`).
* `dodbg` - If set to `true`, the test will stop if the difference between PTC and MAD-NG is outside the tolerance `tol` (default: `false`).
* `dosave` - If set to `true`, the test will save the results to a file (default: `false`).
* `dorun` - If set to `true`, the test will run, if set to `false`, the test will immediately do the save, plot and print steps (default: `true`).
* `doplot` - If set to `true`, the test will plot the results (default: `false`).
* `gen_utest` - If set to `true`, the files required for the unit tests are generated (default: `false`).


These flags can be set when you call the file. After the file, you can set values to the opposite of the default (shown above and in the file).

Run the following to get the help message, describing the flags and how to use them:
```shell
../mad test-all-maps.mad -h
```
or 
```shell
../mad test-all-maps.mad --help
```

The following example runs the crab cavity test and saves the cfg and results files. The next line then does the plot, without running the test again by retrieving the results from the cfg and res files.
```shell
../mad test-electric-maps.mad -t TestElectric.testCCAV -s
../mad test-electric-maps.mad -t TestElectric.testCCAV -f -n
```
But of course you could just save and plot in the same run:
```shell
../mad test-electric-maps.mad -t -s -f TestElectric.testCCAV
```

## Defining the Tests

Within each test, we define an object, inherited from the reference config (which is essentially identical across all files, with changes to order and icase depending on the files need). 

This object must contain the following information:
- The element string, which is the MAD-X element definition, with the parameters to be tested as part of the string, ready to be formatted with the configuration table, see below.
- The following attributes must be define as an iterable:
    - model (1 = DKD, 2 = TKT)
    - method (2, 4, 6, 8, "'teapot2'", "'teapot4'")
    - nslice (Number of slices)
    - energy (Energy of the beam)


Then optionally contains:
- The tolerance, which is the number that the difference between the maps is compared to, it is multiplied by your machine epsilon (`eps`) to get the number to compare to. The tolerance is optional as the default is 10. 
- Additional parameters that are in the element string to be tested, these must be defined as an iterable. If you add parameters, such as k0, this must be added to the attribute list `alist`.
- The information for the plot. This is a list that will be sent directly to the plot function, with a couple adjustments:
    - x_axis_attribute is the value that will be plotted on the x axis, this must be a string and will be created by each row of the configuration. If you had `x_axis_attribute = "${k1}"`, then the x axis would be the k1 value for each row. If this is left blanl, then it is `"${cfgid}"` by default, which is the configuration id.
    - Series is a list of strings that evaluate to a logical, once string interpolation is done using the configuration row. This will be used to create the series' on the plot. If this is left blank, then it is `{"${cfgid} > 0"}` (therefore all the configurations will be plotted in the same series). 
    - If filename is ommited, the plot will not be saved and instead will be displayed on the screen.

Once the configuration object is defined, the test is run by calling `run_test(cfg)`.

All the configurations are placed inside functions, such as `local cfg = ref_cfg "rbend"{ ... }` is the configuration for the rbend and within the function `testRBEND()`, which is called later in the file.

## How the Tests Work

Using the configuration object, the function `run_test` will do the following:
1. Checks if cfg.dorun is not nil or false, if it is, then the function skips to step 5.
2. An mtable with 2 columns is created, the first column will store the configuration snapshots, while the second column will store the results of each configuration snapshot.
3. The configuration snapshots are created using a recursive function that takes the configuration object and runs PTC and MAD-NG on each snapshot in a function called `run_cfg`.
    The function `run_cfg` takes a configuration snapshot and does the following:
    1. Calls `do_trck`, which runs track and ptc to get the difference between MAD and PTC, then stores the results in a table and returns it.
    2. Stores the configuration snapshot and the results into the mtable.
    3. If cfg.doprnt then the results are printed to the console in a single line per config.
    4. If cfg.dodbg then the results are compared to the tolerance, if the results are outside the tolerance then the test stops and the results are printed to the console.
4. Now the mtable has been filled with the results of each configuration snapshot, we create the generator columns for the mtable, which link the tables of each configuration snapshot and the results to a column in the mtable.
5. If cfg.dosave or the file does not exist then the mtable is saved to a file with the name of the element appended with `_cfg` or `_res` for the configuration snapshots and the results respectively.
6. If cfg.doprnt then a summary of the results is printed to the console. (See below for an example of the summary)
7. If cfg.doplot then the results are plotted and saved to a file if cfg.filename is not nil.
8. If the test has not been stopped by cfg.dodbg then several excess files are removed.

## Debugging the Tests

When you run the tests, you can check the results in two ways:
- Check the plot for any particularly high differences.
- Check the results in the results file.

The next steps are up to you, however the best way I found to debug the tests is to do the following:
1. Set `dorun` to `false` (`-n`) and `doprnt` to `true` (`-v`).
2. Run the specific test you want to debug. This will give you an output like so: 
    ```shell
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
    In this list, you can see situations where the failure never occured and situations where changing the variable had no effect. So the variables that had no effect here are ``nslice, energy, icase, k1, tilt, method, model, order``. While the specifics that are missing are `x0i = 1`, `k1s = 0`, `fringe = 0`, indicating the issue only occurs when the initial conditions are `x0i > 1` (not the zero orbit), there is a fringe field, and a `k1s` value.
    
    This step gives you an idea of where the problem may have occured, so you can next move to debugging the specific maps.

3. Set `dorun` to `true` and `dodbg` to `true` (`-d`). This will run the test until your tolerance is reached. Once this happens, the test will stop and run MAD-NG and PTC in debug mode to files called `cfg.name .. "_n"` and `cfg.name .. "_p"` respectively and then a diff of the two files to `cfg.name .. "_d"` (all in your output folder), created by madl_dbgmap.mad. This will give you the differences between every map that is run, so you can see exactly where the problem begins. 
    
    The sequence, mad and madx files that are used to run everything can be found in the input folder.

4. Make changes and repeat step 3 until you have found/resolved the issue.

5. Once you have resolved the issue, you can turn off `dodbg` and `doprnt` and run the test again to get the final results. (Not a necessary step, but it makes the output cleaner)


## The file trackvsptc.mad
This file does the test or unittest, once given a configuration object.

On loading of this file, this file reads in the following files:
- `input/ref.mad` - This file contains the run that loads the sequence (through `MADX`), runs track and returns the results. The string is placed into `mad_ref`.
- `input/ref.madx` - This file contains the MAD-X file that loads the sequence and runs ptc_normal in debug mode. The string is placed into `madx_ref`.

### Functions

#### do_trck (cfg)
This function calls `create_madx_seq` from `track-tool.mad`, then runs the function created by `loadstring(mad_ref%cfg)`. If it is a unittest and unittest generation is off, then the function will read the results from the unittest reference file, otherwise, it will run `run_madx(cfg.name)` from `track-tool.mad` and grab the final map from PTC. If the unittest generation is on, then the results are saved to the unittest reference file. Then the difference between the final maps from MAD-NG and PTC are returned, using `get_diff` from `track-tool.mad`.

#### run_cfg (cfg, results)
This function first gets the diffs by running `do_trck(cfg)`, then runs `store_results` and `prnt_results` from `track-tool.mad`. If `cfg.dodbg` is on, then the function will run `debug_chk` from `track-tool.mad`, which will stop the test if the difference is outside the tolerance and create the file `cfg.name .. "_n.txt"` and run dbgmap to create `cfg.name .. "_d.txt"`. If it is a unittest, then the function will assert `chk_tol` from `track-tool returns true.

#### run_test (cfg)
This function runs the test, first running `args_to_cfg(cfg)` from `test-tool.mad`, to get the test config. PTC model is then turned on before creating the madx reference file and loading the mad reference string specific to this test. 

The function then creates the mtable, writes the printing header, runs `gen_cfg(cfg, 1, \-> run_cfg(cfg, results))`, adds generator columns to the mtable, then saves, prints and plots the results. If the program was not stopped by `debug_chk`, then the function will cleanup the files created by the test.

## The file testvsptc.mad
This file overwrites the default help of the test suite, to include the help for the additional options. While also including the application of the additional options.