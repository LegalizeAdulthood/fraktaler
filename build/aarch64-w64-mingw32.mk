# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021-2023 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

WIN ?= $(HOME)/opt/windows/posix/aarch64
CLEWPREFIX ?= ../clew
CLFLAGS = -I$(CLEWPREFIX)/include -Dclew_STATIC $(CLEWPREFIX)/src/clew.c -Wno-cast-function-type
PKG_CONFIG_PATH = $(WIN)/lib/pkgconfig
PKG_CONFIG_FLAGS = --static
PKG_CONFIG_SED = s/-pthread/-Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic/g
COMPILER = aarch64-w64-mingw32-g++
STRIP = aarch64-w64-mingw32-strip
CFLAGS += -std=$(STDCXX) -Wall -Wextra -pedantic -O3 -MMD
CPPFLAGS += -D__USE_MINGW_ANSI_STDIO=1 -DWINVER=0x501 -D_WIN32_WINNT=0x501 -I$(WIN)/include -I$(WIN)/include/OpenEXR -DHAVE_EXR=$(EXR) -I$(CLEWPREFIX)/include -Dclew_STATIC -DHAVE_CLEW $(CL)
LDFLAGS += -static -static-libgcc -static-libstdc++ -static -L$(WIN)/lib -lpthread
LIBS_IMGUI +=
SOURCES_GUI_C +=
OEXT = .win32-aarch64.o
EXEEXT = .aarch64.exe
TARGETS = fraktaler-3-$(VERSION)$(EXEEXT)
