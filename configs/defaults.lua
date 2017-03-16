return MAD.Object {
  -- methods
  _prepare  = \ctx\ make_dirs(dirname(ctx.madfile)),
  _generate = \ctx\ save_file(ctx.madfile, ctx.whole),
  _run      = \ctx\ run(ctx.command, ctx.madfile, ctx.logfile),
  -- engines
  engines = MAD.Object {
    thick = { whole = LT'basic_track.madx', command = 'madx', makethin=0  },
    madx  = { whole = LT'basic_track.madx', command = 'madx', makethin=1  },
    ptc   = { whole = LT'basic_track.ptc',  command = 'madx'              },
    mad   = { whole = LT'basic_track.mad',  command = 'mad',
              varname = \ctx -> ctx.parent.varname:gsub('->', '.'),       },
    sad   = { whole = LT'basic_track.sad',  command = 'sad'               },
  },
  -- sequence
  madx_sequence  = T[[
${seq_name}: sequence, refer=center, l=${seq_len};
  ${el_name}: ${el_type}, ${el_args};
endsequence;]],
  ptc_sequence = T"${madx_sequence}",
  mad_sequence = T[[
local sequence = MAD.sequence "${seq_name}" {
  refer="entry", l=${seq_len},
  MAD.element.${el_type} "${el_name}" { ${el_args} },
}]],
  sad_sequence = T[[
MARK M1 = ();
${sad_el_type} ${el_name} = (${sad_el_args});
LINE ${seq_name} = (M1 EL);]],
  seq_name  = "SEQ",          -- must be uppercase for SAD
  el_name   = "EL",
  observe   = T"${el_name}",
  seq_len   = 1,
  -- loop:
  varname   = "loopvalue",    -- no underscores in SAD!
  varfunc   = "linrange",
  sad_var   = \ctx ctx.sad_attr and T'Element["${sad_attr}", "${el_name}"]'(ctx) or ctx.varname,
  sad_attr  = false,
  start     = 0,
  stop      = 0.5,
  count     = 20,
  -- particle:
  x = 0, px = 0,
  y = 0, py = 0,
  t = 0, pt = 0,
  energy = 450,
  beam_args = T[[particle="proton", energy=${energy}]],
  mass      = 0.9382720813,   -- particle mass (proton), needed for SAD
  -- control
  studies = MAD.Object {
    x = {}, px = {},
    y = {}, py = {},
    t = {}, pt = {},
    energy = {
      start=1.0, stop=1000,
      varfunc="logrange",
    }
  },
  nst = 1,
  -- output:
  madfile = T"${prefix}.${engine}",
  tfsfile = T"${prefix}.${engine}.tfs",
  logfile = T"${prefix}.${engine}.log",
}
