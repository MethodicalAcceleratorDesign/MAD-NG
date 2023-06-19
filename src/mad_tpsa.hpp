#ifndef MAD_TPSA_HPP
#define MAD_TPSA_HPP

/*
 o-----------------------------------------------------------------------------o
 |
 | Simple C++ wrapper to GTPSA
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - Provide memory management and operator overloading for GTPSA.

 o-----------------------------------------------------------------------------o
 */

#include <memory>

// --- types ------------------------------------------------------------------o

namespace mad {

struct tpsa_del {
  void operator()(tpsa_t* t) { mad_tpsa_del(t); }
};

typedef std::unique_ptr<tpsa_t, tpsa_del> tpsa;

}

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_HPP
