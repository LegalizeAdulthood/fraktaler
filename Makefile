# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

SOURCES_CC = \
src/bla.cc \
src/main.cc \
src/map.cc \
src/reference.cc \
src/render.cc \

SOURCES_H = \
src/bla.h \
src/complex.h \
src/floatexp.h \
src/map.h \
src/param.h \
src/reference.h \
src/render.h \
src/types.h \

SOURCES = $(SOURCES_CC) $(SOURCES_H)

fraktaler-3: $(SOURCES)
	g++ -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -fopenmp -o fraktaler-3 $(SOURCES_CC) `pkg-config --cflags --libs OpenEXR` -lmpfr

fraktaler-3.pdf: README.md
	pandoc README.md -o fraktaler-3.pdf

index.html: README.md
	pandoc README.md -o index.html
