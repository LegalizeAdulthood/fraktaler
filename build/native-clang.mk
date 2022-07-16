# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021,2022 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

COMPILER = clang-11
CFLAGS += -std=c++20 -Wall -Wextra -pedantic -O3 -mtune=native -MMD -DHAVE_GLEW -DHAVE_GLDEBUG -DHAVE_EXR
LDFLAGS += -lstdc++ -lstdc++fs -lm
LIBS_IMGUI += -ldl
LIBS_GUI += glew
LIBS_CL += OpenCL
OEXT = .native-clang.o
EXEEXT = .clang
TARGETS = fraktaler-3-$(VERSION)$(EXEEXT)
