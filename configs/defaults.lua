return MAD.Object {
    -- methods
    _generate = \ctx\ save_file(ctx.madfile, ctx.whole),
    _run      = \ctx\ run(ctx.command, ctx.madfile, ctx.logfile),
    -- engines
    engines = MAD.Object {
      madx = { whole = LT'basic_track.madx', command = 'madx' },
      ptc =  { whole = LT'basic_track.ptc',  command = 'madx' },
      mad =  { whole = LT'basic_track.mad',  command = 'mad'  },
    },
    -- sequence
    madx_sequence  = T[[
${seq_name}: sequence, refer=center, l=${seq_len};
    ${el_name}: ${el_type}, ${el_args};
endsequence;
    ]],
    ptc_sequence = T"${madx_sequence}",
    mad_sequence = T[[
local sequence = MAD.sequence "${seq_name}" {
  refer="entry", l=${seq_len},
  MAD.element.${el_type} "${el_name}" { ${el_args} },
}
    ]],
    seq_name  = "seq",
    el_name   = "el",
    beam_args = T[[particle="proton", energy=${energy}]],
    observe   = T"${el_name}",
    seq_len   = 1,
    -- loop:
    varname   = "loop_value",
    varfunc   = "linrange",
    start     = 0,
    stop      = 0.5,
    count     = 20,
    -- particle:
    x = 0, px = 0,
    y = 0, py = 0,
    t = 0, pt = 0,
    energy = 450,
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
