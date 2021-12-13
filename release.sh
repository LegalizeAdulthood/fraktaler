#!/bin/sh
# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

make -j $(nproc) SYSTEM=i686-w64-mingw32
make clean SYSTEM=i686-w64-mingw32
make -j $(nproc) SYSTEM=x86_64-w64-mingw32
make clean SYSTEM=x86_64-w64-mingw32
make -j $(nproc) SYSTEM=aarch64-w64-mingw32
make clean SYSTEM=aarch64-w64-mingw32
make -j $(nproc) SYSTEM=emscripten
make clean SYSTEM=emscripten
make -j $(nproc) SYSTEM=native-gcc
make clean SYSTEM=native-gcc
make -j $(nproc) SYSTEM=native-clang
make clean SYSTEM=native-clang
make -j $(nproc) SYSTEM=docs
make clean SYSTEM=docs
