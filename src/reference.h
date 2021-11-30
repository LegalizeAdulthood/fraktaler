// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

#include <mpfr.h>

template <typename real>
count_t reference(complex<real> *Z, count_t MaximumReferenceIterations, const mpfr_t Cx, const mpfr_t Cy, progress_t *progress, bool *running);
