import struct, os, subprocess, sys, select
from typing import Union, Callable, Any
import numpy as np

__all__ = ["mad_process"]


def is_not_private(varname):
  if varname[:2] == "__" and varname[:8] != "__last__":
    return False
  return True


class mad_process:
  def __init__(self, mad_path: str, py_name: str = "py", debug: bool = False) -> None:
    self.py_name = py_name

    # Create the pipes for communication
    self.from_mad, mad_write = os.pipe()
    mad_read, self.to_mad = os.pipe()

    # Open the pipes for communication to MAD (the stdin of MAD)
    self.fto_mad = os.fdopen(self.to_mad, "wb", buffering=0)

    # Create a chunk of code to start the process
    startupChunk = (
      f"MAD.pymad '{py_name}' {{_dbg = {str(debug).lower()}}} :__ini({mad_write})"
    )

    # Start the process
    self.process = subprocess.Popen(
      [mad_path, "-q", "-e", startupChunk],
      bufsize=0,
      stdin=mad_read,  # Set the stdin of MAD to the read end of the pipe
      stdout=sys.stdout.fileno(),  # Forward stdout
      preexec_fn=os.setpgrp,  # Don't forward signals
      pass_fds=[
        mad_write,
        sys.stdout.fileno(),
        sys.stderr.fileno(),
      ],  # Don't close these (python closes all fds by default)
    )

    # Close the ends of the pipes that are not used by the process
    os.close(mad_write)
    os.close(mad_read)

    # Create a global variable dictionary for the exec function (could be extended to include more variables)
    self.globalVars = {"np": np}

    # Open the pipe from MAD (this is where MAD will no longer hang)
    self.ffrom_mad = os.fdopen(self.from_mad, "rb")

    # stdout should be line buffered by default, but for jupyter notebook,
    # stdout is redirected and not line buffered by default
    self.send(
      f"""io.stdout:setvbuf('line')
    {self.py_name}:send(1)"""
    )

    # Check if MAD started successfully using select
    checker = select.select([self.ffrom_mad], [], [], 1)  # May not work on windows
    if not checker[0] or self.recv() != 1:  # Need to check number?
      raise OSError(f"Unsuccessful starting of {mad_path} process")

  def send_rng(self, start: float, stop: float, size: int):
    """Send a numpy array as a rng to MAD"""
    self.fto_mad.write(b"rng_")
    send_grng(self, start, stop, size)

  def send_lrng(self, start: float, stop: float, size: int):
    """Send a numpy array as a logrange to MAD"""
    self.fto_mad.write(b"lrng")
    send_grng(self, start, stop, size)

  def send_tpsa(self, monos: np.ndarray, coefficients: np.ndarray):
    """Send the monomials and coeeficients of a TPSA to MAD, creating a table representing the TPSA object"""
    self.fto_mad.write(b"tpsa")
    send_gtpsa(self, monos, coefficients, send_num)

  def send_ctpsa(self, monos: np.ndarray, coefficients: np.ndarray):
    """Send the monomials and coeeficients of a complex TPSA to MAD, creating a table representing the complex TPSA object"""
    self.fto_mad.write(b"ctpa")
    send_gtpsa(self, monos, coefficients, send_cpx)

  def send(self, data: Union[str, int, float, np.ndarray, bool, list]):
    """Send data to MAD, returns self for chaining"""
    try:
      typ = type_str[get_typestr(data)]
      self.fto_mad.write(typ.encode("utf-8"))
      type_fun[typ]["send"](self, data)
      return self
    except KeyError:  # raise not in exception to reduce error output
      raise TypeError(
        f"Unsupported data type, expected a type in: \n{list(type_str.keys())}, got {type(data)}"
      ) from None

  def psend(self, string: str):
    """Perform a protected send to MAD, by first enabling error handling, so that if an error occurs, an error is returned"""
    return self.send(f"{self.py_name}:__err(true); {string}; {self.py_name}:__err(false);")

  def precv(self, name: str):
    """Perform a protected send receive to MAD, by first enabling error handling, so that if an error occurs, an error is received"""
    return self.send(f"{self.py_name}:__err(true):send({name}):__err(false)").recv(name)

  def errhdlr(self, on_off: bool):
    """Enable or disable error handling"""
    self.send(f"{self.py_name}:__err({str(on_off).lower()})")

  def recv(self, varname: str = None):
    """Receive data from MAD, if a function is returned, it will be executed with the argument mad_communication"""
    typ = self.ffrom_mad.read(4).decode("utf-8")
    self.varname = varname  # For mad reference
    return type_fun[typ]["recv"](self)

  def recv_and_exec(self, env: dict = {}):
    """Read data from MAD and execute it"""
     # Check if user has already defined mad (madp_object will have mad defined), otherwise define it
    try:             env["mad"] 
    except KeyError: env["mad"] = self
    exec(compile(self.recv(), "ffrom_mad", "exec"), self.globalVars, env)
    return env

  # ----------------- Dealing with communication of variables ---------------- #
  def send_vars(self, **vars):
    for name, var in vars.items():
      if isinstance(var, mad_ref):
        self.send(f"{name} = {var.__name__}")
      else:
        self.send(f"{name} = {self.py_name}:recv()").send(var)

  def recv_vars(self, *names):
    if len(names) == 1:
      if is_not_private(names[0]):
        return self.precv(names[0])
    else:
      return tuple(self.precv(name) for name in names if is_not_private(name))

  # -------------------------------------------------------------------------- #

  def __del__(self):
    self.send(f"{self.py_name}:__fin()")
    self.ffrom_mad.close()
    self.process.terminate()  # In case user left mad waiting
    self.fto_mad.close()
    self.process.wait()


class mad_ref(object):
  def __init__(self, name: str, mad_proc: mad_process):
    assert name is not None, "Reference must have a variable to reference to. Did you forget to put a name in the receive functions?"
    self.__name__ = name
    self.__mad__ = mad_proc

  def __getattr__(self, item):
    if item[0] != "_":
      try:
        return self[item]
      except (IndexError, KeyError):
        pass
    raise AttributeError(item)  # For python

  def __getitem__(self, item: Union[str, int]):
    if isinstance(item, int):
      result = self.__mad__.precv(f"{self.__name__}[{item+1}]")
      if result is None:
        raise IndexError(item)  # For python
    elif isinstance(item, str):
      result = self.__mad__.precv(f"{self.__name__}['{item}']")
      if result is None:
        raise KeyError(item)  # For python
    else:
      raise TypeError("Cannot index type of ", type(item))

    return result

  def eval(self):
    return self.__mad__.recv_vars(self.__name__)


# data transfer -------------------------------------------------------------- #

# Data ----------------------------------------------------------------------- #


def send_dat(self: mad_process, dat_fmt: str, *dat: Any):
  self.fto_mad.write(struct.pack(dat_fmt, *dat))


def recv_dat(self: mad_process, dat_sz: int, dat_typ: np.dtype):
  return np.frombuffer(self.ffrom_mad.read(dat_sz), dtype=dat_typ)


# None ----------------------------------------------------------------------- #

send_nil = lambda self, input: None
recv_nil = lambda self: None

# Boolean -------------------------------------------------------------------- #

send_bool = lambda self, input: self.fto_mad.write(struct.pack("?", input))
recv_bool = lambda self: recv_dat(self, 1, np.bool_)[0]

# int32 ---------------------------------------------------------------------- #

send_int = lambda self, input: send_dat(self, "i", input)
recv_int = lambda self: recv_dat(self, 4, np.int32)[
  0
]  # Should it be a python int or a numpy int32?

# String --------------------------------------------------------------------- #


def send_str(self: mad_process, input: str):
  send_int(self, len(input))
  self.fto_mad.write(input.encode("utf-8"))


def recv_str(self: mad_process) -> str:
  return self.ffrom_mad.read(recv_int(self)).decode("utf-8")


# number (float64) ----------------------------------------------------------- #

send_num = lambda self, input: send_dat(self, "d", input)
recv_num = lambda self: recv_dat(self, 8, np.float64)[0]

# Complex (complex128) ------------------------------------------------------- #

send_cpx = lambda self, input: send_dat(self, "dd", input.real, input.imag)
recv_cpx = lambda self: recv_dat(self, 16, np.complex128)[0]

# Range ---------------------------------------------------------------------- #

send_grng = lambda self, start, stop, size: send_dat(self, "ddi", start, stop, size)


def recv_rng(self: mad_process) -> np.ndarray:
  return np.linspace(*struct.unpack("ddi", self.ffrom_mad.read(20)))


def recv_lrng(self: mad_process) -> np.ndarray:
  return np.geomspace(*struct.unpack("ddi", self.ffrom_mad.read(20)))


# irange --------------------------------------------------------------------- #

send_irng = lambda self, rng: send_dat(self, "iii", rng.start, rng.stop, rng.step)

def recv_irng(self: mad_process) -> range:
  start, stop, step = recv_dat(self, 12, np.int32)
  return range(start, stop + 1, step)  # MAD is inclusive at both ends


# matrix --------------------------------------------------------------------- #


def send_gmat(self: mad_process, mat: np.ndarray):
  assert len(mat.shape) == 2, "Matrix must be of two dimensions"
  send_dat(self, "ii", *mat.shape)
  self.fto_mad.write(mat.tobytes())


def recv_gmat(self: mad_process, dtype: np.dtype) -> str:
  shape = recv_dat(self, 8, np.int32)
  return recv_dat(self, shape[0] * shape[1] * dtype.itemsize, dtype).reshape(shape)


recv_mat = lambda self: recv_gmat(self, np.dtype("float64"))
recv_cmat = lambda self: recv_gmat(self, np.dtype("complex128"))
recv_imat = lambda self: recv_gmat(self, np.dtype("int32"))

# monomial ------------------------------------------------------------------- #


def send_mono(self: mad_process, mono: np.ndarray):
  send_int(self, mono.size)
  self.fto_mad.write(mono.tobytes())

recv_mono = lambda self: recv_dat(self, recv_int(self), np.ubyte)

# TPSA ----------------------------------------------------------------------- #


def send_gtpsa(
  self: mad_process,
  monos: np.ndarray,
  coefficients: np.ndarray,
  send_num: Callable[[mad_process, Union[float, complex]], None],
):
  assert len(monos.shape) == 2, "The list of monomials must have two dimensions"
  assert len(monos) == len(coefficients), "The number of monomials must be equal to the number of coefficients"
  assert monos.dtype == np.uint8, "The monomials must be of type 8-bit unsigned integer "
  send_dat(self, "ii", len(monos), len(monos[0]))
  for mono in monos:
    self.fto_mad.write(mono.tobytes())
  for coefficient in coefficients:
    send_num(self, coefficient)


def recv_gtpsa(self: mad_process, dtype: np.dtype) -> np.ndarray:
  num_mono, mono_len = recv_dat(self, 8, np.int32)
  mono_list = np.reshape(
    recv_dat(self, mono_len * num_mono, np.ubyte),
    (num_mono, mono_len),
  )
  coefficients = recv_dat(self, num_mono * dtype.itemsize, dtype)
  return mono_list, coefficients


recv_ctpa = lambda self: recv_gtpsa(self, np.dtype("complex128"))
recv_tpsa = lambda self: recv_gtpsa(self, np.dtype("float64"))

# lists ---------------------------------------------------------------------- #


def send_list(self: mad_process, lst: list):
  send_int(self, len(lst))
  for item in lst:
    self.send(item)


def recv_list(self: mad_process) -> list:
  varname = self.varname  # cache
  haskeys = recv_bool(self)
  lstLen = recv_int(self)
  vals = [self.recv(varname and varname + f"[{i+1}]") for i in range(lstLen)]
  self.varname = varname  # reset
  if haskeys and lstLen == 0:
    return type_fun["ref_"]["recv"](self)
  elif haskeys:
    return vals, type_fun["ref_"]["recv"](self)
  else:
    return vals


# object (table with metatable are treated as pure reference) ---------------- #

recv_ref = lambda self: mad_ref(self.varname, self)
send_ref = lambda self, obj: send_str(self, f"return {obj.__name__}")

# error ---------------------------------------------------------------------- #


def recv_err(self: mad_process):
  self.errhdlr(False)
  raise RuntimeError("MAD Errored (see the MAD error output)")


# ---------------------------- dispatch tables ------------------------------- #
type_fun = {
  "nil_": {"recv": recv_nil , "send": send_nil },
  "bool": {"recv": recv_bool, "send": send_bool},
  "str_": {"recv": recv_str , "send": send_str },
  "tbl_": {"recv": recv_list, "send": send_list},
  "ref_": {"recv": recv_ref , "send": send_ref },
  "fun_": {"recv": recv_ref , "send": send_ref },
  "obj_": {"recv": recv_ref , "send": send_ref },
  "int_": {"recv": recv_int , "send": send_int },
  "num_": {"recv": recv_num , "send": send_num },
  "cpx_": {"recv": recv_cpx , "send": send_cpx },
  "mat_": {"recv": recv_mat , "send": send_gmat},
  "cmat": {"recv": recv_cmat, "send": send_gmat},
  "imat": {"recv": recv_imat, "send": send_gmat},
  "rng_": {"recv": recv_rng ,                  },
  "lrng": {"recv": recv_lrng,                  },
  "irng": {"recv": recv_irng, "send": send_irng},
  "mono": {"recv": recv_mono, "send": send_mono},
  "tpsa": {"recv": recv_tpsa,                  },
  "ctpa": {"recv": recv_ctpa,                  },
  "err_": {"recv": recv_err ,                  },
}


def get_typestr(a: Union[str, int, float, np.ndarray, bool, list]):
  if isinstance(a, np.ndarray):
    return a.dtype
  elif type(a) is int:  # Check for signed 32 bit int
    if a.bit_length() < 31:
      return int
    else:
      return float
  else:
    return type(a)


type_str = {
  type(None)              : "nil_",
  bool                    : "bool",
  str                     : "str_",
  list                    : "tbl_",
  tuple                   : "tbl_",
  mad_ref                 : "ref_",
  int                     : "int_",
  np.int32                : "int_",
  float                   : "num_",
  np.float64              : "num_",
  complex                 : "cpx_",
  np.complex128           : "cpx_",
  range                   : "irng",
  np.dtype("float64")     : "mat_",
  np.dtype("complex128")  : "cmat",
  np.dtype("int32")       : "imat",
  np.dtype("ubyte")       : "mono",
}
# ---------------------------------------------------------------------------- #
