# Test Suite for NG vs NG maps

## Running the Tests
For help on how to run the tests, run the following:
```shell
../mad test-all-maps.mad -h
```

### Types of tests
Most of the details of this test suite is identical to the test suite for ng vs ptc. The main difference is the type of test that is run.

By default, the test will run:
- chg vs -chg
- edir vs -edir
- sdir vs -sdir
- lua vs cmap


### Creating your element
Creating your element in these tests are almost identical to that in test-ptc-maps, however there is a requirement that you must add `${bdir}, ${tdir}` to the element, following the convention:

Usage of directions:
  - bending angles are multiplied by `tdir`
  - strengths      are multiplied by `bdir`

Note: `sdir` has no effect on the element, as this is never changed in the tests. This is because, instead of comparing the element to itself, with `sdir=-1`, we backtrack through the result of going forwards through the element, therefore never changing the sequence.

### Custom Tests
You can also add more tests, such as rotation of an element, converting the skew component of an element to a normal component. To do this, when you call `run_test`, you can add a table of the following form, as the second argument to `run_test`:
```lua
local equiv = object "quad" {
    rotate = object { -- Requires deferred expression
      "tilt", "k1", "k1s", -- These are the parameters that will be varied
      tilt = {"math.pi/4", "-math.pi/4"}, -- Must be the same length as the other parameters, equal to n
      k1  := {-cfg.cur_cfg.k1s,  cfg.cur_cfg.k1s}, -- cfg is the configuration table, cur_cfg is the current configuration
      k1s := { cfg.cur_cfg.k1 , -cfg.cur_cfg.k1 },
      
      n = 2, -- This is the number of tests to run
    },

    alist = {"rotate"}, -- List of parameters to be tested (must be in the overall table)
  }
```

In this object, we define a custom test called rotate, this tells the program that three parameters will be varied, `tilt, k1, k1s`, and that the number of tests to run is 2. 
The program will look into alist, take the table, and loop `n` times, each time setting the parameters to the values in the table, so for the first loop, it will set `tilt = math.pi/4, k1 = -cfg.cur_cfg.k1s, k1s = cfg.cur_cfg.k1`, and for the second loop, it will set `tilt = -math.pi/4, k1 = cfg.cur_cfg.k1s, k1s = -cfg.cur_cfg.k1`. 

The program will compare each of these tests to the reference test, which is the test with the parameters set to the default values (that you set) in the configuration table. Therefore the purpose of custom tests is to create a test that is equivalent to the reference test, but with a different configuration.

## The file trackvsng.mad
This file does the test or unittest, once given a configuration object.

### Variables

#### Configuration strings

```lua
local track_str = [[
-- Constants
local models = {'DKD', 'TKT'}
local X0s = 
]] .. tostring(X0s, ",\n") .. [[
local X0 = X0s[${x0i}]
if ${order} > 0 then 
  X0 = MAD.damap{nv = 6, mo = ${order}}
  for i, c in ipairs({"x", "px", "y", "py", "t", "pt"}) do 
    X0[c]:set0(X0s[${x0i}][c])
  end
end

local drift, quadrupole, sextupole, octupole, sbend, rbend, solenoid, rfcavity,
      kicker, elseparator, crabcavity, multipole, xrotation, yrotation, 
      srotation, translate, changeref, decapole, dodecapole, rfmultipole         in MAD.element

local particles = {}
for i = 1, ${npar} do particles[i] = X0 end

local seq = MAD.sequence "seq" {dir = ${edir}, l=${seql}, ${elm}}
return MAD.track { -- see ref_cfg for list of values
  dir      = ${sdir},
  beam     = MAD.beam {energy = ${energy}, charge=${chg}},
  sequence = seq,
  X0       = particles,
! mapdef   = ${order},
  model    = models[${model}],
  method   = ${method},
  nslice   = ${nslice},
  debug    = ${debug},
  cmap     = ${cmap},
  aperture = {kind='circle', 10},
}
]]

local backtrack_str = [[
local _, mflw = MAD.track { -- see ref_cfg for list of values
  dir      = ${sdir},
  beam     = MAD.beam {energy = ${energy}, charge=${chg}},
  sequence = seq,
  X0       = mflw0,
  model    = models[${model}],
  method   = ${method},
  nslice   = ${nslice},
  cmap     = ${cmap},
  debug    = ${debug},
}
return mflw, particles
]]
```

The first string is used in all the tests, it is the string that is used to track through the element and return the results. The second string is only used for the backtracking tests, where we track through the element, then track backwards through the element, returning the results and the initial particle coordinates.

### Functions

#### create_run (cfg, cur_cfg, track_str_)
This function creates the string that is loaded by mad to run the test. It takes the configuration object, the current configuration, interpolating `track_str` (see above) with `cfg` and `cur_cfg`, then returning the string. 

If `track_str_` is given, then this string is used instead of `track_str`. This is used for the backtracking tests, where the following is sent to the function:
```lua
track_str_ = cfg.ref_script:gsub("return", "local _, mflw0 = ") .. backtrack_str 
```

`ref_script` is the string that is created at the very beginning of the test, which is the reference test, with all attributes set to the default values. See above for `backtrack_str`.

#### debug_chk (cfg, dif, script2) (Could be improved)
This function is used to check if the difference between the two maps is within the tolerance (if debug mode is on). If the difference is greater than the tolerance, then the function will cause the test to stop, creating a two files in the output directory, `setup1.mad` and `setup2.mad`, which are the two scripts that were used to create the two maps that disagree.

#### do_unit_test(cfg, dif)
If the test is a unit test, then this function will assert that `chk_tol(dif, cfg.tol, cfg.order)` (from track_tool.mad) results in true. 

#### dif_save_prnt_dbg (cfg, exp, res, result_mtbl, script) 
This function runs through all of the particles that were just tracked, setting the test  type based on which particle and if cmap is on, running the following:

```lua
local dif_tbl = get_diff(exp[i], res[i], cfg.order) -- get_diff is from track_tool.mad
store_results(cfg, dif_tbl, result_mtbl) -- store_results is from track_tool.mad
prnt_results( cfg, dif_tbl) -- prnt_results is from track_tool.mad
debug_chk    (cfg, dif_tbl, script) -- debug_chk is defined above
do_unit_test (cfg, dif_tbl) -- do_unit_test is defined above
```

`test_type` is used in the configuration mtable, so that you can look at the results and see which test it is. 

#### backtrack (cfg, results)
If the backtracking is on, then this function will run forwards and backwards through the same configuration, checking if the output is the same as the input coordinates of the particles. 

#### reverse_attr (cfg, results, attr)
This function compares the results from the default configuration to the results from the configuration after performing `cfg[attr] = -cfg[attr]`. 

#### reverse_chg (cfg, results)
If reversing charge is on, then this function calls `reverse_attr` with `attr = "chg"`.

#### reverse_edir (cfg, results)
If reversing edir is on, then this function calls `reverse_attr` with `attr = "edir"`.

#### cmaps (cfg, results)
This function sets `cfg.cmap` to `true`, creating the results from the default configuration with `cfg.cmap = true`. Then if reversing cmap is on, it will save the results, comparing to the results from the default configuration with `cfg.cmap = false`. 

Then, while `cfg.cmap` is true, it will run the following before resetting `cfg.cmap` to false:
```lua
backtrack   (cfg, results) -- Backtrack
reverse_chg (cfg, results) -- Change sign of chg
reverse_edir(cfg, results) -- Change sign of edir
```

#### run_cfg (cfg, equiv, results)
This function runs the test for the given configuration, first generating the reference script and map (for comparison later), then runs `cmaps`, `backtrack`, `reverse_chg` and `reverse_edir`. Then if custom tests are on, it will run the custom tests, using the table `equiv` (see above).

#### gen_plot_cfg (cfg)
This function generates the defualt plot configuration, which is used for all the plots, if `plot_info.series` or `plot_info.legend` is not defined. This function is required, since series is defined through an iterable, by the legend is defined through a mappable that needs to match the order of the series. Therefore if you have an empty series, it may take the wrong legend.

#### run_test (cfg, equiv)
This is the entry function for all the tests, where the user inputs the configuration object and the optional table of custom tests. This function first creates the default values for the configuration object (`edir, sdir, chg, tdir, bdir, cmap`), then inherits the global configuration created from the command line (by calling `args_to_cfg`).

The function then creates the mtable, writes the printing header, runs `gen_cfg(cfg, 1, \-> run_cfg(cfg, equiv or {}, results))`, adds generator columns to the mtable, then saves, prints and plots the results. If the program was not stopped by `debug_chk`, then the function will cleanup the files created by the test.

## The file testvsng.mad
This file overwrites the default help of the test suite, to include the help for the additional options. While also including the application of the additional options.