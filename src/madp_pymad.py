import struct, os, subprocess, sys, select
from typing import Union, Callable, Any
import numpy as np

__all__ = ["mad_process"]

def get_typestring(a: Union[str, int, float, np.ndarray, bool, list]):
  if isinstance(a, np.ndarray):
    return a.dtype
  elif type(a) is int: # Check for signed 32 bit int
    if a.bit_length() < 31: return int
    else:                   return float
  else:
    return type(a)
  
data_types = {
  type(None)              : "nil_",
  str                     : "str_",
  int                     : "int_",
  np.int32                : "int_",
  float                   : "num_",
  np.float64              : "num_",
  complex                 : "cpx_",
  np.complex128           : "cpx_",
  bool                    : "bool",
  list                    : "tbl_",
  range                   : "irng",
  np.dtype("float64")     : "mat_",
  np.dtype("complex128")  : "cmat",
  np.dtype("int32")       : "imat",
  np.dtype("ubyte")       : "mono",
}
class mad_process:
  def __init__(self, mad_path: str, py_name: str = "py", debug: bool = False) -> None:
    self.py_name = py_name

    self.from_mad, mad_write = os.pipe()
    mad_read, self.to_mad = os.pipe()
    self.fto_mad = os.fdopen(self.to_mad, "wb", buffering=0) # Sensible to not buffer stdin?
    
    startupChunk = f"MAD.pymad '{py_name}' {{_dbg = {str(debug).lower()}}} :__ini({mad_write})"

    self.process = subprocess.Popen(
      [mad_path, "-q", "-e", startupChunk],
      bufsize=0,
      stdin=mad_read,
      stdout=sys.stdout.fileno(),
      preexec_fn=os.setpgrp,  # Don't forward signals
      pass_fds=[mad_write, sys.stdout.fileno(), sys.stderr.fileno()],
    )
    os.close(mad_write)
    os.close(mad_read)

    self.globalVars = {"np" : np} 
    self.ffrom_mad = os.fdopen(self.from_mad, "rb")

    # stdout should be line buffered by default, but for jupyter notebook, 
    # stdout is redirected and not line buffered by default
    self.send(
      f"""io.stdout:setvbuf('line')
        {self.py_name}:send(1)"""
    )
    checker = select.select([self.ffrom_mad], [], [], 1) # May not work on windows
    if not checker[0] or self.recv() != 1: # Need to check number?
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
      typ = data_types[get_typestring(data)]
      self.fto_mad.write(typ.encode("utf-8"))
      str_to_fun[typ]["send"](self, data)
      return self
    except KeyError:  # raise not in exception to reduce error output
      pass
    raise TypeError(
      f"Unsupported data type, expected a type in: \n{list(data_types.keys())}, got {type(data)}"
    )

  def psend(self, string: str):
    """Perform a protected send to MAD, by first enabling error handling, so that if an error occurs, an error is returned"""
    return self.send(f"{self.py_name}:__err(true); {string}; {self.py_name}:__err(false);")
  
  def precv(self, name: str):
    """Perform a protected send receive to MAD, by first enabling error handling, so that if an error occurs, an error is received"""
    return self.send(f"{self.py_name}:__err(true):send({name}):__err(false)").recv(name)

  def errhdlr(self, on_off: bool):
    """Enable or disable error handling"""
    self.send(f"{self.py_name}:__err({str(on_off).lower()})")

  def recv(
    self, varname: str = None
  ) -> Union[str, int, float, np.ndarray, bool, list]:
    """Receive data from MAD, if a function is returned, it will be executed with the argument mad_communication"""
    typ = self.ffrom_mad.read(4).decode("utf-8")
    self.varname = varname  # For mad reference
    return str_to_fun[typ]["recv"](self)

  def recv_and_exec(self, env: dict = {}) -> dict:
    """Read data from MAD and execute it"""
    try:             env["mad"]        # Check if user has already defined mad (madp_object will have mad defined)
    except KeyError: env["mad"] = self # If not, define it
    exec(compile(self.recv(), "ffrom_mad", "exec"), self.globalVars, env)
    return env

  # ----------------- Dealing with communication of variables ---------------#
  def send_vars(self, names, vars):
    if isinstance(names, str): 
      names = [names]
      vars = [vars]
    else:
      assert isinstance(vars, list), "A list of names must be matched with a list of variables"
      assert len(vars) == len(names), "The number of names must match the number of variables"
    for i, var in enumerate(vars):
      if isinstance(vars[i], mad_ref):
        self.send(f"{names[i]} = {var.__name__}")
      else:
        self.send(f"{names[i]} = {self.py_name}:recv()").send(var)

  def recv_vars(self, names) -> Any:
    if isinstance(names, str): 
      names = [names]
      cnvrt = lambda rtrn: rtrn[0]
    else: 
      cnvrt = lambda rtrn: tuple(rtrn)

    rtrn_vars = []
    for name in names:
      if name[:2] != "__" or name[:8] == "__last__":  # Check for private variables
        rtrn_vars.append(self.precv(name))        
    return cnvrt(rtrn_vars)

  # -------------------------------------------------------------------------#

  def __del__(self):
    self.send(f"{self.py_name}:__fin()")
    self.ffrom_mad.close()
    self.process.terminate() #In case user left mad waiting
    self.fto_mad.close()
    self.process.wait()
    
class mad_ref(object):    
  def __init__(self, name: str, mad_proc: mad_process) -> None:
    assert name is not None, "Reference must have a variable to reference to. Did you forget to put a name in the receive functions?"
    self.__name__ = name
    self.__mad__ = mad_proc

  def __getattr__ (self, item):
    if item[0] != "_":
      try:
        return self[item]
      except (IndexError, KeyError):
        pass
    raise AttributeError (item)  # For python

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

data_types[mad_ref] = "ref_" # Add mad_ref to the datatypes

# -------------------------------- Sending data -------------------------------#
send_nil = lambda self, input: None

def send_int(self: mad_process, input: int) -> None:
  self.fto_mad.write(struct.pack("i", input))

def send_str(self: mad_process, input: str) -> None:
  send_int(self, len(input))
  self.fto_mad.write(input.encode("utf-8"))

def send_ref(self: mad_process, obj: mad_ref) -> None:
  send_str(self, f"return {obj.__name__}")

def send_num(self: mad_process, input: float) -> None:
  self.fto_mad.write(struct.pack("d", input))

def send_cpx(self: mad_process, input: complex) -> None:
  self.fto_mad.write(struct.pack("dd", input.real, input.imag))

def send_bool(self: mad_process, input: bool) -> None:
  self.fto_mad.write(struct.pack("?", input))

def send_gmat(self: mad_process, mat: np.ndarray) -> None:
  assert len(mat.shape) == 2, "Matrix must be of two dimensions"
  send_int(self, mat.shape[0])
  send_int(self, mat.shape[1])
  self.fto_mad.write(mat.tobytes())

def send_list(self: mad_process, lst: list) -> None:
  n = len(lst)
  send_int(self, n)
  for item in lst:
    self.send(item)  # deep copy
  return self

def send_grng(self: mad_process, start: float, stop: float, size: int) -> None:
  self.fto_mad.write(struct.pack("ddi", start, stop, size))

def send_irng(self: mad_process, rng: range) -> None:
  self.fto_mad.write(struct.pack("iii", rng.start, rng.stop, rng.step))

def send_mono(self: mad_process, mono: np.ndarray) -> None:
  send_int(self, mono.size)
  self.fto_mad.write(mono.tobytes())

def send_gtpsa(
  self: mad_process,
  monos: np.ndarray,
  coefficients: np.ndarray,
  fsendNum: Callable[[mad_process, Union[float, complex]], None],
) -> None:
  assert len(monos.shape) == 2, "The list of monomials must have two dimensions"
  assert len(monos) == len(coefficients), "The number of monomials must be equal to the number of coefficients"
  assert monos.dtype == np.uint8, "The monomials must be of type 8-bit unsigned integer "
  send_int(self, len(monos))  # Num monomials
  send_int(self, len(monos[0]))  # Monomial length
  for mono in monos:
    self.fto_mad.write(mono.tobytes())
  for coefficient in coefficients:
    fsendNum(self, coefficient)
# -----------------------------------------------------------------------------#

# --------------------------- Receiving data ----------------------------------#

recv_nil = lambda self: None

def recv_ref(self: mad_process) -> mad_ref:
  return mad_ref(self.varname, self)

def recv_str(self: mad_process) -> str:
  return self.ffrom_mad.read(recv_int(self)).decode("utf-8")

def recv_int(self: mad_process) -> int:  # Must be int32
  return int.from_bytes(self.ffrom_mad.read(4), sys.byteorder, signed=True)

def recv_num(self: mad_process) -> float:
  return np.frombuffer(self.ffrom_mad.read(8), dtype=np.float64)[0]

def recv_cpx(self: mad_process) -> complex:
  return np.frombuffer(self.ffrom_mad.read(16), dtype=np.complex128)[0]

def recv_bool(self: mad_process) -> str:
  return np.frombuffer(self.ffrom_mad.read(1), dtype=np.bool_)[0]

def recv_gmat(self: mad_process, dtype: np.dtype) -> str:
  shape = np.frombuffer(self.ffrom_mad.read(8), dtype=np.int32)
  arraySize = shape[0] * shape[1] * dtype.itemsize
  return np.frombuffer(self.ffrom_mad.read(arraySize), dtype=dtype).reshape(shape)

def recv_mat(self: mad_process) -> str:
  return recv_gmat(self, np.dtype("float64"))

def recv_cmat(self: mad_process) -> str:
  return recv_gmat(self, np.dtype("complex128"))

def recv_imat(self: mad_process) -> str:
  return recv_gmat(self, np.dtype("int32"))

def recv_list(self: mad_process) -> list:
  varname = self.varname  # cache
  haskeys = recv_bool(self)
  lstLen = recv_int(self)
  vals = [self.recv(varname and varname + f"[{i+1}]") for i in range(lstLen)]
  self.varname = varname  # reset
  if haskeys and lstLen == 0:
    return recv_ref(self)
  elif haskeys:
    return vals, recv_ref(self)
  else:
    return vals

def recv_irng(self: mad_process) -> range:
  start, stop, step = np.frombuffer(self.ffrom_mad.read(12), dtype=np.int32)
  return range(start, stop + 1, step)  # MAD is inclusive at both ends

def recv_rng(self: mad_process) -> np.ndarray:
  return np.linspace(*struct.unpack("ddi", self.ffrom_mad.read(20)))

def recv_lrng(self: mad_process) -> np.ndarray:
  return np.geomspace(*struct.unpack("ddi", self.ffrom_mad.read(20)))

def recv_mono(self: mad_process) -> np.ndarray:
  mono_len = recv_int(self)
  return np.frombuffer(self.ffrom_mad.read(mono_len), dtype=np.ubyte)

def recv_gtpsa(self: mad_process, dtype: np.dtype) -> np.ndarray:
  num_mono, mono_len = np.frombuffer(self.ffrom_mad.read(8), dtype=np.int32)
  mono_list = np.reshape(
    np.frombuffer(self.ffrom_mad.read(mono_len * num_mono), dtype=np.ubyte),
    (num_mono, mono_len),
  )
  coefficients = np.frombuffer(
    self.ffrom_mad.read(num_mono * dtype.itemsize), dtype=dtype
  )
  return mono_list, coefficients

def recv_ctpa(self: mad_process):
  return recv_gtpsa(self, np.dtype("complex128"))

def recv_tpsa(self: mad_process):
  return recv_gtpsa(self, np.dtype("float64"))

def recv_err(self: mad_process):
  self.errhdlr(False)
  raise RuntimeError("MAD Errored (see the MAD error output)")

# -----------------------------------------------------------------------------#

# ----------------------------- Dict for data ---------------------------------#
str_to_fun = {
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
# ---------------------------------------------------------------------------- #