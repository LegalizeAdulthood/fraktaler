// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <cstdio>

#include "source.h"

const char *fraktaler_3_source_7z = "fraktaler-3-" FRAKTALER_3_VERSION_STRING ".7z";

#if DEFAULT_LINKAGE_HAS_UNDERSCORE
extern const char  binary_fraktaler_3_source_7z_start[],  binary_fraktaler_3_source_7z_end[];
#else
extern const char _binary_fraktaler_3_source_7z_start[], _binary_fraktaler_3_source_7z_end[];
#endif

int write_file(const char *name, const char *start, const char *end)
{
  std::fprintf(stderr, "writing '%s'... ", name);
  FILE *file = std::fopen(name, "wxb");
  if (! file)
  {
    std::fprintf(stderr, "FAILED\n");
    return 0;
  }
  int ok = (std::fwrite(start, end - start, 1, file) == 1);
  std::fclose(file);
  if (ok)
  {
    std::fprintf(stderr, "ok\n");
  }
  else
  {
    std::fprintf(stderr, "FAILED\n");
  }
  return ok;
}

bool write_source(const std::filesystem::path &file)
{
#if DEFAULT_LINKAGE_HAS_UNDERSCORE
  return write_file(file.string().c_str(),  binary_fraktaler_3_source_7z_start,  binary_fraktaler_3_source_7z_end);
#else
  return write_file(file.string().c_str(), _binary_fraktaler_3_source_7z_start, _binary_fraktaler_3_source_7z_end);
#endif
}
