# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

EMSCRIPTEN ?= $(HOME)/opt/emscripten
COMPILER = em++
CFLAGS += -std=c++20 -Wall -Wextra -pedantic -O3 -MMD -I$(EMSCRIPTEN)/include -s USE_SDL=2 -s USE_PTHREADS -DHAVE_GLEW
LDFLAGS += -L$(EMSCRIPTEN)/lib -lgmp -lmpfr -s USE_SDL=2 -s ALLOW_MEMORY_GROWTH=1 -s USE_PTHREADS -s PTHREAD_POOL_SIZE=4 -s MAX_WEBGL_VERSION=2
OEXT = .emscripten.o
EXEEXT =
TARGETS = live/$(VERSION)/index.html
