// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

struct phybrid;

extern const char *nt_string[
#ifdef HAVE_FLOAT128
  7
#else
  6
#endif
];

extern number_type nt_current;

extern std::string pref_path;

void populate_number_type_wisdom(void);
void delete_ref();
void delete_bla();
count_t getM(number_type nt, count_t phase);
void reference_thread(stats &sta, param &par, progress_t *progress, bool *running, bool *ended);
void subframe_thread(map &out, stats &sta, const param &par, const count_t subframe, progress_t *progress, bool *running, bool *ended);
