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

str_t
mad_str_trim (str_t str, ssz_t len[2])
{
  assert(str && len);

  while (len[1] > 0 && isspace(str[len[0]         ])) --len[1], ++len[0];
  while (len[1] > 0 && isspace(str[len[0]+len[1]-1])) --len[1];

  return str;
}

str_t
mad_str_split (str_t str, ssz_t len[5], str_t sep)
{
  assert(str && len && sep);

  ssz_t i = len[0], k = len[0]+len[1];
  while (i < k && str[i] != *sep) ++i;

  if (i == k) // no sep found
    len[2] = len[3] = len[4] = -1;
  else {
    len[1] = i-len[0];
    len[2] = i+1;
    len[3] = k-(i+1);
    len[4] = i;
    mad_str_trim(str, len+2);
  }

  return mad_str_trim(str, len);
}

str_t
mad_str_split_bracket (str_t str, ssz_t len[6])
{
  assert(str && len);

  ssz_t i = len[0], k = len[0]+len[1];
  while (i < k && str[i] != '[' && str[i] != '{'
               && str[i] != ']' && str[i] != '}') ++i;

  ssz_t j = k-1;
       if (str[i] == '[') while (j > i && str[j] != ']') --j;
  else if (str[i] == '{') while (j > i && str[j] != '}') --j;
  else j = i;

       if (i == k) // no bracket found
    len[2] = len[3] = len[4] = len[5] = -1;
  else if (i == j) // bad pairs of brackets
    len[2] = len[3] = len[5] = -1, len[4] = i;
  else {
    len[1] = i-len[0];
    len[2] = i+1;
    len[3] = j-(i+1);
    len[4] = i;
    len[5] = j;
    mad_str_trim(str, len+2);
  }

  return mad_str_trim(str, len);
}
