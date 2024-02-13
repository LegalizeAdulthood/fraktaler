// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include "gles2.h"

struct colour;
struct tile;

colour *colour_new();
void colour_delete(colour *c);
void colour_set_image_size(colour *c, int width, int height);
void colour_set_zoom_log_2(colour *c, float zoom_log_2);
void colour_set_time(colour *c, float time);
void colour_tile(colour *c, int x, int y, int subframe, tile *data);
