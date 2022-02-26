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

#define newton_action_period 0
#define newton_action_center 1
#define newton_action_zoom 2
#define newton_action_transform 3

#define newton_zoom_minibrot_relative 0
#define newton_zoom_minibrot_absolute 1
#define newton_zoom_domain_relative 2
#define newton_zoom_domain_absolute 3

void populate_number_type_wisdom(void);
void delete_ref();
void delete_bla();
count_t getM(number_type nt, count_t phase);
void reference_thread(stats &sta, param &par, bool just_did_newton, progress_t *progress, volatile bool *running, volatile bool *ended);
void subframe_thread(map &out, stats &sta, const param &par, const count_t subframe, progress_t *progress, volatile bool *running, volatile bool *ended);
void newton_thread(param &out, bool &ok, const param &par, const complex<floatexp> &c, const floatexp &r, volatile progress_t *progress, volatile bool *running, volatile bool *ended);
