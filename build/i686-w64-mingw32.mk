# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

WIN ?= $(HOME)/win/posix/i686
PKG_CONFIG_PATH = $(WIN)/lib/pkgconfig
PKG_CONFIG_FLAGS = --static
PKG_CONFIG_SED = s/-pthread/-Wl,-Bstatic -lstdc++ -lstdc++fs -lpthread -Wl,-Bdynamic/g
COMPILER = i686-w64-mingw32-g++
CFLAGS += -std=c++20 -Wall -Wextra -pedantic -O3 -MMD
CPPFLAGS += -D__USE_MINGW_ANSI_STDIO=1 -DWINVER=0x501 -D_WIN32_WINNT=0x501 -I$(WIN)/include -I$(WIN)/include/OpenEXR -I$(WIN)/src/glew-2.1.0/include -DGLEW_STATIC -DHAVE_GLEW -DHAVE_GLDEBUG -DHAVE_EXR
LDFLAGS += -static -static-libgcc -static-libstdc++ -L$(WIN)/lib
LIBS_IMGUI += -lopengl32
SOURCES_GUI_C += $(WIN)/src/glew-2.1.0/src/glew.c
OEXT = .win32-i686.o
EXEEXT = .i686.exe
TARGETS = fraktaler-3-$(VERSION)-cli$(EXEEXT) fraktaler-3-$(VERSION)-gui$(EXEEXT)
