# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021-2023 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

EMSCRIPTEN ?= $(HOME)/opt/emscripten/pthreads-exceptions
COMPILER = em++
CFLAGS += -std=$(STDCXX) -Wall -Wextra -pedantic -fexceptions -O3 -MMD -I$(EMSCRIPTEN)/include -s USE_SDL=2 -s USE_PTHREADS -DHAVE_GUI -DIMGUI_IMPL_OPENGL_ES2
LDFLAGS += -fexceptions -L$(EMSCRIPTEN)/lib -lgmp -lmpfr -lidbfs.js -s USE_SDL=2 -s ALLOW_MEMORY_GROWTH=1 -s USE_PTHREADS -s PTHREAD_POOL_SIZE="(navigator.hardwareConcurrency+2)" -s MAX_WEBGL_VERSION=2 -s ASYNCIFY=1 -s EXPORTED_RUNTIME_METHODS=ccall
OEXT = .emscripten.o
EXEEXT =
TARGETS = live/$(VERSION)/index.html
