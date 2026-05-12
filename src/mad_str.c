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

#include <assert.h>
#include "mad_str.h"

// --- implementation ---------------------------------------------------------o

static inline int mad_isspace(chr_t c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

static inline int mad_isdigit(chr_t c) {
  return (u32_t)(c - '0') <= 9u;
}

static inline int mad_islower(chr_t c) {
  return (u32_t)(c - 'a') <= 25u;
}

static inline int mad_isupper(chr_t c) {
  return (u32_t)(c - 'A') <= 25u;
}

static inline int mad_isalpha(chr_t c) {
  return mad_islower(c) || mad_isupper(c);
}

static inline int mad_isalnum(chr_t c) {
  return mad_isalpha(c) || mad_isdigit(c);
}

static inline chr_t mad_tolower(chr_t c) {
  return mad_isupper(c) ? (chr_t)(c + ('a' - 'A')) : c;
}

static inline chr_t mad_toupper(chr_t c) {
  return mad_islower(c) ? (chr_t)(c - ('a' - 'A')) : c;
}

static inline void
mad_str_trim_front (str_t str, ssz_t arg[2])
{
  while (arg[1] > 0 && mad_isspace(str[arg[0]])) --arg[1], ++arg[0];
}

static inline void
mad_str_trim_back (str_t str, ssz_t arg[2])
{
  while (arg[1] > 0 && mad_isspace(str[arg[0]+arg[1]-1])) --arg[1];
}

str_t
mad_str_trim (str_t str, ssz_t arg[2])
{
  assert(str && arg);
  assert(arg[0] >= 0);
  assert(arg[1] >= 0);

  mad_str_trim_front(str, arg);
  mad_str_trim_back (str, arg);
  return str;
}

str_t
mad_str_quote (str_t str, ssz_t arg[5], log_t sq)
{
  assert(str && arg);
  assert(arg[0] >= 0);
  assert(arg[1] >= 0);

  mad_str_trim_front(str, arg);

  if (!arg[1] || (str[arg[0]] != '"' && (str[arg[0]] != '\'' || !sq))) {
    arg[2] = -1, arg[3] = arg[4] = 0; // no quote found
    return str;
  }

  idx_t i = arg[0], j = i+1, k = arg[0]+arg[1], q = 0;

  if (str[i] == '"')
    while (j < k && str[j] != '"' )
      j += (str[j] == '\\' && (j+1) < k && str[j+1] == '"' ) ? ++q, 2 : 1;
  else if (sq)
    while (j < k && str[j] != '\'')
      j += (str[j] == '\\' && (j+1) < k && str[j+1] == '\'') ? ++q, 2 : 1;

  if (j == k) {
    arg[2] = -1, arg[3] = arg[4] = 0;
    return NULL; // error: no closing quote found
  }

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
  assert(arg[0] >= 0);
  assert(arg[1] >= 0);

  mad_str_trim_front(str, arg);

  if (!arg[1]) {
    arg[2] = arg[3] = arg[4] = -1, arg[5] = 0; // no bracket found
    return str;
  }

  idx_t i = arg[0], k = arg[0]+arg[1];

  while (i < k && str[i] != '[' && str[i] != '{'
               && str[i] != ']' && str[i] != '}') ++i;

  if (i == k) { // no bracket found
    arg[2] = arg[3] = arg[4] = -1, arg[5] = 0;
    mad_str_trim_back(str, arg);
    return str;
  }

  if (str[i] == ']' || str[i] == '}') { // error: no opening bracket
    arg[2] = arg[3] = arg[4] = -1, arg[5] = 0;
    return NULL;
  }

  idx_t j = i+1;
  while (j < k && str[j] != '[' && str[j] != '{'
               && str[j] != ']' && str[j] != '}') ++j;

  if (j == k || (str[i] == '[' && str[j] != ']')    // error: no closing bracket
             || (str[i] == '{' && str[j] != '}')) { // or invalid nested opening
    arg[2] = arg[3] = arg[4] = -1, arg[5] = 0;
    return NULL;
  }

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
  assert(arg[0] >= 0);
  assert(arg[1] >= 0);
  assert(arg[2] >  0);

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
  assert(arg[0] >= 0);
  assert(arg[1] >= 0);

  mad_str_trim_front(str, arg);

  if (!arg[1]) {
    arg[2] = arg[3] = arg[4] = -1; // no number found
    return str;
  }

  idx_t i = arg[0], k = arg[0]+arg[1], d = -1, e = -1, n = 0;

  // sign
  if (i < k && (str[i] == '-' || str[i] == '+')) ++i;

  // inf or nan ?
  if (i < k && mad_isalpha(str[i])) {
    if (i+2 < k && mad_tolower(str[i  ]) == 'i' &&  mad_tolower(str[i+1]) == 'n' &&
                   mad_tolower(str[i+2]) == 'f' && (i+3 == k || !mad_isalpha(str[i+3]))) {
      i += 3; goto fini;
    }
    if (i+2 < k && mad_tolower(str[i  ]) == 'n' &&  mad_tolower(str[i+1]) == 'a' &&
                   mad_tolower(str[i+2]) == 'n' && (i+3 == k || !mad_isalpha(str[i+3]))) {
      i += 3; goto fini;
    }
    arg[1] = 0, arg[2] = arg[3] = arg[4] = -1; // no number found
    return str;
  }

  // integer part
  while(i < k && mad_isdigit(str[i])) ++i, ++n;

  // decimal part
  if (i < k && str[i] == '.') {
    d = i++;

    // concat ..
    if (i < k && str[i] == '.') {
      if (n) { i -= 2, d = -1; goto fini; }
      arg[1] = 0, arg[2] = arg[3] = arg[4] = -1; // no number found
      return str;
    }

    while(i < k && mad_isdigit(str[i])) ++i, ++n;
  }

  // ensure at least ±# or ±#. or ±.#
  if(!n && d >= 0) {
    arg[1] = 0, arg[2] = arg[3] = arg[4] = -1; // no number found
    return str;
  }

  // exponent part
  if (i < k && (str[i] == 'e' || str[i] == 'E')) {
    e = i++;

    // sign
    if (i < k && (str[i] == '-' || str[i] == '+')) ++i;

    // digits
    while(i < k && mad_isdigit(str[i])) ++i;

    // ensure e# or e±# otherwise backtrack
    if (!mad_isdigit(str[i-1])) { i = e, e = -1; goto fini; }
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
  assert(arg[0] >= 0);
  assert(arg[1] >= 0);

  mad_str_trim_front(str, arg);

  if (!arg[1]) {
    arg[2] = arg[3] = -1; // no identifier found
    return str;
  }

  idx_t i = arg[0], k = arg[0]+arg[1];

  if (i >= k || (!mad_isalpha(str[i]) && str[i] != '_')) {
    arg[1] = 0, arg[2] = arg[3] = -1; // no identifier found
    return str;
  }

  ++i; while (i < k && (mad_isalnum(str[i]) || str[i] == '_')) ++i;

  arg[1] = i-arg[0]; // len
  arg[2] = i; // index right after identifier

  while (i < k && mad_isspace(str[i])) ++i;

  arg[3] = i; // index right after trailing spaces

  return str;
}
