/*
 o-----------------------------------------------------------------------------o
 |
 | String module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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
*/

#include <ctype.h>
#include <assert.h>

#include "mad_str.h"

// --- implementation ---------------------------------------------------------o

const char*
mad_str_trim (const char *str, ssz_t *len)
{
  assert(str && len);

  while (*len > 0 && isspace(str[   0  ])) --*len, ++str;
  while (*len > 0 && isspace(str[*len-1])) --*len;

  return str;
}
