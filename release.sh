#!/bin/sh
# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

make SYSTEM=i686-w64-mingw32
make clean SYSTEM=i686-w64-mingw32
make SYSTEM=x86_64-w64-mingw32
make clean SYSTEM=x86_64-w64-mingw32
make SYSTEM=aarch64-w64-mingw32
make clean SYSTEM=aarch64-w64-mingw32
make SYSTEM=emscripten
make clean SYSTEM=emscripten
make SYSTEM=native-clang
make clean SYSTEM=native-clang
make SYSTEM=docs
make clean SYSTEM=docs
