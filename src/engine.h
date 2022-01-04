// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

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
void delete_bla();
count_t getM(number_type nt);
bool convert_reference(const number_type to, const number_type from);
bool convert_bla(const number_type to, const number_type from);
void reference_thread(stats &sta, const formula *form, param &par, progress_t *progress, bool *running, bool *ended);
void subframe_thread(map &out, stats &sta, const formula *form, const param &par, const count_t subframe, progress_t *progress, bool *running, bool *ended);
