// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "types.h"

#include <mpfr.h>

template <typename real>
count_t run_reference(complex<real> *Z, const count_t MaximumReferenceIterations, reference *ref, progress_t *progress, bool *running);
