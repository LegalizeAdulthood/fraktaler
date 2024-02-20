// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2024 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <string>

#include "types.h"

// true on success
bool save_jpeg_rgb8(const std::string &filename, const unsigned char *data, coord_t width, coord_t height, const std::string &comment, int quality);
bool save_jpeg_yuv8(const std::string &filename, const unsigned char *data, coord_t width, coord_t height, const std::string &comment, int quality);

// "" on failure
std::string read_jpeg_comment(const std::string &filename);
