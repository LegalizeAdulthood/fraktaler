// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

enum number_type
{
  nt_none = 0,
  nt_float = 1,
  nt_double = 2,
  nt_longdouble = 3,
  nt_floatexp = 4
};

extern const char *nt_string[];

extern number_type nt_current;

void delete_bla();
count_t getM(number_type nt);
bool convert_reference(const number_type to, const number_type from);
bool convert_bla(const number_type to, const number_type from);
void reference_thread(stats &sta, const formula *form, const param &par, progress_t *progress, bool *running, bool *ended);
void subframe_thread(map &out, stats &sta, const formula *form, const param &par, const count_t subframe, progress_t *progress, bool *running, bool *ended);