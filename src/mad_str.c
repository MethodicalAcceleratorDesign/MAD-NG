/*
 o-----------------------------------------------------------------------------o
 |
 | String module implementation
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
*/

#include <ctype.h>
#include <assert.h>

#include "mad_str.h"

// --- implementation ---------------------------------------------------------o

static inline void
mad_str_trim_front (str_t str, ssz_t arg[2])
{
  while (arg[1] > 0 && isspace(str[arg[0]])) --arg[1], ++arg[0];
}

static inline void
mad_str_trim_back (str_t str, ssz_t arg[2])
{
  while (arg[1] > 0 && isspace(str[arg[0]+arg[1]-1])) --arg[1];
}

str_t
mad_str_trim (str_t str, ssz_t arg[2])
{
  assert(str && arg);
  mad_str_trim_front(str, arg);
  mad_str_trim_back (str, arg);
  return str;
}

str_t
mad_str_quote (str_t str, ssz_t arg[5], log_t sq)
{
  assert(str && arg);
  mad_str_trim_front(str, arg);

  if (str[arg[0]] != '"' && (str[arg[0]] != '\'' || !sq)) { // no quote found
    arg[2] = -1, arg[3] = arg[4] = 0;
    return str;
  }

  idx_t i = arg[0], j = i+1, k = arg[0]+arg[1], q = 0;

  if (str[i] == '"')
    while (j < k && str[j] != '"' )
      j += (str[j] == '\\' && str[j+1] == '"' ) ? ++q, 2 : 1;
  else if (sq)
    while (j < k && str[j] != '\'')
      j += (str[j] == '\\' && str[j+1] == '\'') ? ++q, 2 : 1;

  if (j == k) return NULL; // error: no closing quote found

  arg[0] = i+1;
  arg[1] = j-(i+1);
  arg[2] = j;
  arg[3] = (str[i] == '\'' && sq) + 1;
  arg[4] = q;

  return str;
}

str_t
mad_str_bracket (str_t str, ssz_t arg[6])
{
  assert(str && arg);
  mad_str_trim_front(str, arg);

  idx_t i = arg[0], k = arg[0]+arg[1];

  while (i < k && str[i] != '[' && str[i] != '{'
               && str[i] != ']' && str[i] != '}') ++i;

  if (i == k) { // no bracket found
    arg[2] = arg[3] = arg[4] = -1, arg[5] = 0;
    mad_str_trim_back(str, arg);
    return str;
  }

  if (str[i] == ']' || str[i] == '}') // error: no opening bracket
    return NULL;

  idx_t j = i+1;
  while (j < k && str[j] != '[' && str[j] != '{'
               && str[j] != ']' && str[j] != '}') ++j;

  if (str[j] == '[' || str[j] == '{') // error: no closing bracket
    return NULL;

  arg[1] = i-arg[0];
  arg[2] = j;
  arg[3] = i+1;
  arg[4] = j-(i+1);
  arg[5] = (str[i] == '{') + 1;
  mad_str_trim_back(str, arg  );
  mad_str_trim     (str, arg+3);

  return str;
}

str_t
mad_str_split (str_t str, ssz_t arg[4], str_t sep)
{
  assert(str && arg && sep);

  idx_t i = arg[0], j = -1, k = arg[0]+arg[1], l = arg[2];

  switch(l) {
    case 1: while (i < k && str[i] != sep[j=0]) ++i; break;
    case 2: while (i < k && str[i] != sep[j=0] && str[i] != sep[j=1]) ++i; break;
    case 3: while (i < k && str[i] != sep[j=0] && str[i] != sep[j=1] && str[i] != sep[j=2]) ++i; break;
    case 4: while (i < k && str[i] != sep[j=0] && str[i] != sep[j=1] && str[i] != sep[j=2] && str[i] != sep[j=3]) ++i; break;
    default: for (;i < k; ++i) for (j=0; j<l; ++j) if (str[i] == sep[j]) goto found;
  }

found:
  if (i == k) { // no sep found
    arg[2] = arg[3] = -1;
  } else {
    arg[1] = i-arg[0];
    arg[2] = i;
    arg[3] = j;
  }
  return mad_str_trim(str, arg);
}

str_t
mad_str_num (str_t str, ssz_t arg[5])
{
  assert(str && arg);
  mad_str_trim_front(str, arg);

  idx_t i = arg[0], d = -1, e = -1, n = 0;

  // sign
  if (str[i] == '-' || str[i] == '+') ++i;

  // inf or nan ?
  if (isalpha(str[i])) {
    if (tolower(str[i  ]) == 'i' &&  tolower(str[i+1]) == 'n' &&
        tolower(str[i+2]) == 'f' && !isalpha(str[i+3])) {
      i += 3; goto fini;
    }
    if (tolower(str[i  ]) == 'n' &&  tolower(str[i+1]) == 'a' &&
        tolower(str[i+2]) == 'n' && !isalpha(str[i+3])) {
      i += 3; goto fini;
    }
    arg[1] = 0, arg[2] = arg[3] = arg[4] = -1; // no number found
    return str;
  }

  // integer part
  while(isdigit(str[i])) ++i, ++n;

  // decimal part
  if (str[i] == '.') {
    d = i++;

    // concat ..
    if (str[i] == '.') {
      if (n) { i -= 2, d = -1; goto fini; }
      arg[1] = 0, arg[2] = arg[3] = arg[4] = -1; // no number found
      return str;
    }

    while(isdigit(str[i])) ++i, ++n;
  }

  // ensure at least ±# or ±#. or ±.#
  if(!n && d >= 0) {
    arg[1] = 0, arg[2] = arg[3] = arg[4] = -1; // no number found
    return str;
  }

  // exponent part
  if (str[i] == 'e' || str[i] == 'E') {
    e = i++;

    // sign
    if (str[i] == '-' || str[i] == '+') ++i;

    // digits
    while(isdigit(str[i])) ++i;

    // ensure e# or e±# otherwise backtrack
    if (!isdigit(str[i-1])) { i = e-1, e = -1; goto fini; }
  }

fini:
  arg[1] = i-arg[0]; // len
  arg[2] = i; // index right after
  arg[3] = d;
  arg[4] = e;

  return str;
}

str_t
mad_str_ident (str_t str, ssz_t arg[4])
{
  assert(str && arg);
  mad_str_trim_front(str, arg);

  idx_t i = arg[0];

  if (!isalpha(str[i]) && str[i] != '_') {
    arg[1] = 0, arg[2] = arg[3] = -1; // no identifier found
    return str;
  }

  ++i; while (isalnum(str[i]) || str[i] == '_') ++i;

  arg[1] = i-arg[0]; // len
  arg[2] = i; // index right after identifier

  while (isspace(str[i])) ++i;

  arg[3] = i; // index right after leading spaces

  return str;
}
