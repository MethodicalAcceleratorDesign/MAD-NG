# What each tool does 

## Track Tool

On loading this file, a new folder is created in the current directory called `output`.  

### Functions

#### run_madx (fnam)
This function runs a madx file called `input/fnam_ref.madx`, where `fnam` is the input argument converted to a string, sending the output to `output/fnam_p.txt`, where `fnam` is the input argument converted to a string.

#### create_madx_seq (cfg)
This function creates a madx sequence file called `"input/"..cfg.name.."_seq.seq"`, where `cfg.name` is the name of the configuration object. 

The contents of the file is determined with string interpolation of the configuration object in the following string:

```lua
local seq_ctx = [[
x0  = ${x};
px0 = ${px};
y0  = ${y};
py0 = ${py};
t0  = ${t};
pt0 = ${pt};
model    = ${model};
method   = ${method};
nslice   = ${nslice};
energy   = ${energy};
snm      = ${snm};
order    = ${order};
icase    = ${icase};
debug    = ${debug};
ptcmodel = ${ptcmodel};
cmap     = ${cmap};
seq: sequence, l=${seql} ;
  elm: ${elm};
endsequence ;
]]
```

The values of `x`, `px`, `y`, `py`, `t` and `pt` are determined by the `x0i` attribute of the configuration object, which is an index to the table `X0s` (all possible configurations). This table is also defined in this file as:

```lua
local X0s = {{x=0   , px=0    , y=0    , py=0   , t=0   , pt=0   }, -- zero
             {x=3e-3, px=-2e-4, y=-2e-3, py=3e-4, t=0   , pt=0   }, -- 4D
             {x=3e-3, px=-2e-4, y=-2e-3, py=3e-4, t=0   , pt=2e-5}, -- 5D weak
             {x=3e-3, px=-2e-4, y=-2e-3, py=3e-4, t=0   , pt=3e-2}, -- 5D strong
             {x=3e-3, px=-2e-4, y=-2e-3, py=3e-4, t=1e-5, pt=2e-5}} -- 6D
```
All other values are determined by the attributes of the configuration object.

#### get_last_ptc_map (filename)
This function reads the output of PTC from the file `filename`, and returns the final DA map as a `damap`

#### get_diff (ref, res, order)
This function uses the `dif` method of the `damap` to find the difference between the reference map `ref` and the result map `res` up to order `order`. 

Two variables are returned, the first a table containing the difference by coordinate and order, and the second the actual difference returned by the `dif` method.

The table is of the form (for `order=4`):

```lua
dif_tbl = {
    x_eps = 3, px_eps = 2, y_eps = 1, py_eps = 2, t_eps = 3, pt_eps = 2,
    order0_eps = 0, order1_eps = 1, order2_eps = 2, order3_eps = 3, order4_eps = 3,
}
```

#### chk_tol (res, tol, order)
This function loops through the table `res` (of the same form as `dif_tbl` above) and checks if the difference is less than the tolerance `tol` (in units of `eps`) up to order `order`. If the difference is less than the tolerance, the function returns `true`, otherwise it returns `false`.

#### store_results (cfg, res, results)
This function adds a new row to the mtable `results`, placing a copy of the current configuration `cfg` of the test into a column called `__cfg`, and the results table `res` (of the same form as `dif_tbl` above) in a column called `__res`.

#### prnt_results (cfg, res)
If `cfg.doprnt`, this function prints the results of the test to the terminal. \
During the test, this is intended to print per single configuration. The first number is the config id, and then the difference in each order is printed (in units of `eps`). The configuration is then printed in the form `key=value,` for each key in the table `alist` in the configuration object (see below).

```
cfgid   order 0 order 1 order 2 order 3 order 4
1       0       1       2       2       3       model=1,  energy=1,  method=2,  nslice=1,  x0i=1,  k0=0.5,  no_fringe=false
```

#### gen_cfg (cfg, idx, gen_fun)
This function generates many single configurations using the table `alist` in the configuration object, and the corresponding iterables placed into the configuration object.

This is a recursive function, continuously calling itself until all combinations of the iterables are exhausted. To initialise the function, the `idx` argument should be set to `1`, and the `gen_fun` argument should be set a function in which you would like to be called for each configuration. The function `gen_fun` should not take any arguments.

#### add_gen_cols (results, cfg)
This function adds generator columns to the mtable `results`. 

One set of columns columns is created from the table `alist` in the configuration object, creating a column for each string in the table that refers to the corresponding key in the column `__cfg` (see above).

The other set of columns is created from the table `res_cols` in the mtable `results`, creating a column for each string in the table that refers to the corresponding key in the column `__res` (see above).

#### add_trk_gen_cols (results, cfg)
This is a wrapper around the above function by first adding to the table `res_cols` the corresponding keys to the table `dif_tbl` (see above).

#### get_lower_bnds (res, tol) (Consistency)
This function returns a function `(order, row_index)` that returns the tolerance based on the order and row index, using the mtable `res` and `tol`. 

The input `tol` can be a file name with columns `orderi_eps` (where i is the order from 0 - the order of the test), or a number, in which case the same tolerance is used for all orders, or an table where the index is the order+1 (0 order is the first index) and the value is the tolerance.

#### show_res (res, cfg, attr_cols, tol) (Consistency)
This function prints the final results of the test to the terminal, using the mtable `res` and `cfg`. The variable `attr_cols` should be the table `cfg.alist`, while `tol` is a function that takes order and row_index to return a tolerance \

The output is of the form:

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
Where each line is an attribute from the table `attr_cols`. 
The keys of the table are the values of the attribute, while the values are the number of configurations that fail the tolerance for that attribute. This can be used to determine which attributes are causing the test to fail.

#### get_prev_res (test_name, out_dir) (out_dir not needed)
This function reads the configuration and results files from the test `out_dir/test_name` and returns the mtables `cfg` and `res`.

#### save_res (cfg, res, cfg_cols_, res_cols_, hdr_list_)
This function saves the configuration and results tables to the files `"output/"..res.name.."_cfg.txt"` and `"output/"..res.name.."_res.txt"` respectively, where `res.name` is the name of the results .

The variables `cfg_cols_`, `res_cols_` are used to determine which columns in the config and results files respectively to save. If these are not defined, then `cfg.alist` and `res.res_cols` are used respectively. 

The variable `hdr_list_` is used to determine what is output to the header of the mtables. If this is not defined, then `{"name", "date", "time", "origin", "max_order", "run_tol"}` is stored. Note, `run_tol` is the tolerance used for the test and is used in show_res, so is important to store.

#### tbl2da (tbl)
This function converts a table of the form:

```lua
{x=3e-3, px=-2e-4, y=-2e-3, py=3e-4, t=0, pt=0}
```

to a `damap`, with an `mo` of 1 and `nv` of 6, setting the 0 order values to each coordinate in the table.

#### in_dir (s)
This function returns a string with the `s` appended to the `input/` directory.

#### out_dir (s)
This function returns a string with the `s` appended to the `output/` directory.

### Variables

#### ptc_strs 
This table contains the strings that are used to run PTC, and is used in the function `run_cfg`. The strings are:

```lua
local coord_str = [[
x0  = ${x};
px0 = ${px};
y0  = ${y};
py0 = ${py};
t0  = ${t};
pt0 = ${pt};
]]

local seq_ctx = [[
seq: sequence, l=${seql} ;
  elm: ${elm};
endsequence ;
]]
local madx_script = [[
${coords}
model    = ${model};
method   = ${method};
nslice   = ${nslice};
energy   = ${energy};
snm      = ${snm};
order    = ${order};
icase    = ${icase};
debug    = ${debug};
ptcmodel = ${ptcmodel};
cmap     = ${cmap};
${seq_ctx}
]]

local ptc_strs = {
  coord_str   = coord_str,
  seq_ctx     = seq_ctx,
  madx_script = madx_script,
}
```

This means the sequence from `create_madx_seq` can be adjusted by changing the one of the strings in this table.

### Return Values

The following functions and variables are returned by this file:

run_madx        ,
get_diff        ,
chk_tol         ,
store_results   ,
prnt_results    ,
gen_cfg         ,
add_gen_cols    ,
add_trk_gen_cols,
show_res        ,
get_prev_res    ,
get_lower_bnds  ,
X0s             ,
create_madx_seq ,
get_last_ptc_map,
save_res        ,
tbl2da          ,
in_dir          ,
out_dir         ,
ptc_strs.

## Plot Tool

On loading this file, a new folder is created in the current directory called `plots` in the output directory imported from the track tool. 


### Functions

#### newSID ()
This function returns a new SID for the plot, which allows for multiple plots to be created at the same time.

#### plot_res (get_y_val, res_cfg_tbl, cfg, plt_dir, cfg_tbl_) (plot_dir not needed?)
This function plots the results of the test. 

The variable `get_y_val` is a function that takes a row of the results mtable and `cfg`, returning the y value for the plot.

The variable `res_cfg_tbl` is a either an mtable that contains just the results or an mtable that contains both the results and the configuration. If it just contains the results, then `cfg_tbl_` must be defined as the mtable that contains the configuration.

The variable `cfg` is the configuration object of the test.

The variable `plt_dir` is the directory to save the plot to. 

Within the configuration object, a table called `plot_jnfo` can be defined. \
`series` is a list of strings that evaluate to a logical, once string interpolation is done using the configuration row. This will be used to create the series' on the plot. If this is left blank, then it is `{"${cfgid} > 0"}` (therefore all the configurations will be plotted in the same series).
`x_axis_attribute` is the value that will be plotted on the x axis, this must be a string and will be created by each row of the configuration. If you had `x_axis_attribute = "${k1}"`, then the x axis would be the k1 value for each row. If this is left blank, then it is `"${cfgid}"` by default, which is the configuration id. \
All other values are passed directly to the plot function.

#### get_trk_y_val (res, cfg)
This function returns the maximum `eps` error in the row `res` of the results mtable multiplied by the machine epsilon `eps`.

#### plot_trk_res (res_cfg_tbl, cfg, plt_dir, cfg_tbl_)
This function is a wrapper around `plot_res` with the function `get_trk_y_val` as the first argument.

#### do_norun (cfg)
This function completes the actions of the test without running the test by reading the configuration and results files of the test.

#### plt_dir (s)
This function returns a string with the `s` appended to the `plots/` directory, inside the output directory.

### Return Values
The following functions and variables are returned by this file:

plot_res     ,
plot_trk_res ,
get_trk_y_val,
plt_dir      ,
newSID       ,
do_norun.


## Test Tool

### Variables

#### arg_fun
This table contains the default functions to be called based on an argument input on the command line.

#### args_dict
This table contains the default conversion of the arguments from the command line to a string that indicates which function call to make, if it exists in the table `arg_fun`, otherwise flip the boolean value of the argument in the configuration object.

For example:

```lua
local args_dict = { 
  ["--save"   ] = "dosave" ,
  ["-s"       ] = "dosave" ,

  ["-t"       ] = "test"   ,
  ["--test"   ] = "test"   ,
}
local arg_fun = {
  test    = fnil,
}
```

This means that if the argument `--save` or `-s` is input, then the flag `dosave` in the configuration object will be set to the opposite of the default value. If the argument `--test` or `-t` is input, then the function `fnil` will be called.

#### global_cfg
This table contains the configuration set by the command line arguments. By default, it is set to:

```lua
local global_cfg = {dorun = true}
```

Thereforel, from the command line, the only argument in the configuration object that can be set to `false` is `dorun`.

### Functions

#### remove_arg (is_unittest, i, arg)
To prevent the lua unittest from throwing errors, this function removes the arguments that are not needed for the unittests. A table called `args_in_utest` is defined in this file, which contains the arguments that are needed for the unittests and is the determining factor for which arguments are removed (along with is_unittest).

#### parse_cmd_args (is_unittest)
This function parses the command line arguments and sets the corresponding flags into a global configuration object called `global_cfg` or runs the function corresponding to the argument. It also deals with inputs that do not begin with a `-` or `--` by assuming that the input is a test name or class name, storing this information in `tests_to_run`.

This returns a boolean if the test is not a unittest (so returns `true` if its just a test).

### Functions (only for non unittests)

#### chk_test (mod_name, test_name)
This function checks if the test should be run. If the test name is not input, then it checks if a keep or exclusion pattern has been input on the command line. 

The rules for patterns and exclusions are as follows:
- If the test name does not match every pattern, then the test is not run.
- If the test name matches any exclusion, then the test is not run.
- Note: should this affect module names?

If the test name is input, then it checks the module name and test name if they were input on the command line.

#### chk_unknown_tests (ran_tests)
This function checks if any tests were input on the command line that do not exist. If there are, then the function prints to the terminal.

#### is_test_name (s)
If the input is a string and the first four characters is `test` independent of case, then this function returns `true`, otherwise it returns `false`.

#### run_tests ()
This function loops through the global environment finds all the classes that are tests, and runs all functions that pass `is_test_name` and `chk_test`.

### Return Values
The following functions and variables are returned by this file:

`parse_cmd_args,
args_to_cfg,
run_tests,
args_dict,
arg_fun,
global_cfg.
`