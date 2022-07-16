// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "param.h"

extern param par;

int batch_cli(int verbosity);
int batch_cl(int verbosity);
int gui(const char *progname, int verbosity, const char *persistence);

extern std::string pref_path;
