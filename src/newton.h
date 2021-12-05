// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "complex.h"
#include "types.h"

template <typename real>
count_t find_period(const complex<real> *Zp, const count_t M, const complex<real> c, const count_t N, const real r, progress_t *progress, bool *running);
