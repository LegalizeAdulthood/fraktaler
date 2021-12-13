# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

COMPILER = clang-11
CFLAGS += -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -fopenmp -DHAVE_OMP -MMD
LDFLAGS += -lstdc++ -lstdc++fs -lm
OEXT = .native-clang.o
EXEEXT =
TARGETS = fraktaler-3-$(VERSION)-cli$(EXEEXT) fraktaler-3-$(VERSION)-gui$(EXEEXT)
