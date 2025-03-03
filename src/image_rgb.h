// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#pragma once

#include <string>

#include "image.h"

struct image_rgb : public image
{
  float *RGBA;
  std::mutex mutex;
  image_rgb(coord_t width, coord_t height);
  image_rgb(image_rgb &source, bool vflip = false);
  virtual ~image_rgb();
  virtual void clear() override;
  virtual void blit(coord_t tx, coord_t ty, const struct tile *t) override;
  virtual bool save_exr(const std::string &filename, int threads, const std::string &metadata, const std::string &kf2plus_metadata = "");
};
