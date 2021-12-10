// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <mpreal.h>
#ifdef HAVE_EXR
#include ...
#endif
#ifdef HAVE_GUI
#include <SDL2/SDL.h>
#include "imgui.h"
#endif

#include "version.h"

std::string version()
{
  std::ostringstream out;
  out << "fraktaler-3 version " << FRAKTALER_3_VERSION_STRING << "\n";
#ifdef __GNUC__
  out << "g++ version " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
#else
#ifdef __clang__
  out << "clang version " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << "\n";
#else
  out << "compiler version unknown\n";
#endif
#endif
  out << "gmp version " << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << "\n";
  out << "mpfr version " << MPFR_VERSION_MAJOR << "." << MPFR_VERSION_MINOR << "." << MPFR_VERSION_PATCHLEVEL << "\n";
  out << "mpfrc++ version " << MPREAL_VERSION_MAJOR << "." << MPREAL_VERSION_MINOR << "." << MPREAL_VERSION_PATCHLEVEL << "\n";
#ifdef HAVE_EXR
  out << "zlib version " << zlib_version << "\n";
  out << "ilmbase version " << ILMBASE_VERSION_STRING << "\n";
  out << "openexr version " << OPENEXR_VERSION_STRING << "\n";
#endif
#ifdef HAVE_GUI
  out << "sdl2 version " << SDL_MAJOR_VERSION << "." << SDL_MINOR_VERSION << "." << SDL_PATCHLEVEL << "\n";
  out << "glew version " << GLEW_VERSION_STRING << "\n";
  out << "imgui version " << IMGUI_VERSION << " (" << IMGUI_GIT_VERSION_STRING << ")\n";
//  out << "imgui-filebrowser version git (" << IMGUI_FILEBROWSER_GIT_VERSION_STRING << ")\n";
#endif
//  out << "tomlplusplus version " << TOML_LANG_MAJOR << "." << TOML_LANG_MINOR << "." << TOML_LANG_PATCH << " (" << TOMLPLUSPLUS_GIT_VERSION_STRING << ")\n";
//  out << "miniaudio version " << MA_VERSION_STRING << " (" << MINIAUDIO_GIT_VERSION_STRING << ")\n";
//  out << "fftw version " << fftw_version << "\n";
  return out.str();
}
