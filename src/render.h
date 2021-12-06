// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

template <typename real>
void render(map &out, stats &sta, const param &par, const real Zoom, const count_t M, const complex<real> *Zp, const formula *form, progress_t *progress, bool *running);
