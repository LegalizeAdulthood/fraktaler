# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

SOURCES_CC = \
src/bla.cc \
src/display.cc \
src/main.cc \
src/map.cc \
src/param.cc \
src/reference.cc \
src/render.cc \

SOURCES_H = \
src/bla.h \
src/complex.h \
src/display.h \
src/floatexp.h \
src/main.h \
src/map.h \
src/param.h \
src/reference.h \
src/render.h \
src/types.h \

SOURCES = $(SOURCES_CC) $(SOURCES_H)

all: fraktaler-3-glfw fraktaler-3-sdl2 fraktaler-3.pdf index.html

fraktaler-3-glfw: $(SOURCES) src/main_glfw.cc
	g++ -std=c++20 -Wall -Wextra -pedantic -Og -march=native -fopenmp -o fraktaler-3-glfw $(SOURCES_CC) src/main_glfw.cc `pkg-config --cflags --libs glew glfw3 mpfr OpenEXR` -ggdb

fraktaler-3-sdl2: $(SOURCES) src/main_sdl2.cc
	g++ -std=c++20 -Wall -Wextra -pedantic -Og -march=native -fopenmp -o fraktaler-3-sdl2 $(SOURCES_CC) src/main_sdl2.cc `pkg-config --cflags --libs glew mpfr OpenEXR sdl2` -ggdb

fraktaler-3.pdf: README.md
	pandoc README.md -o fraktaler-3.pdf

index.html: README.md
	pandoc README.md -o index.html
