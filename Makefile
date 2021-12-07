# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

CFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -fopenmp -MMD
LIBS = glew glm mpfr OpenEXR sdl2

CFLAGS_IMGUI = -I../imgui -I../imgui/backends -I../imgui/misc/cpp
LIBS_IMGUI = -ldl

COMPILE = g++ $(CFLAGS) `pkg-config --cflags $(LIBS)` $(CFLAGS_IMGUI)
LINK = g++ $(CFLAGS)
LINK_FLAGS = `pkg-config --libs $(LIBS)` $(LIBS_IMGUI)

SOURCES_CC = \
src/bla.cc \
src/colour.cc \
src/display.cc \
src/formula.cc \
src/main.cc \
src/main_sdl2.cc \
src/map.cc \
src/param.cc \
src/stats.cc \

SOURCES_H = \
src/bla.h \
src/colour.h \
src/colour_rainbow.h \
src/complex.h \
src/display.h \
src/dual.h \
src/floatexp.h \
src/formula.h \
src/formula_burningship.h \
src/formula_mandelbrot.h \
src/main.h \
src/map.h \
src/matrix.h \
src/param.h \
src/stats.h \
src/types.h \

SOURCES_IMGUI_CC = \
../imgui/imgui.cpp \
../imgui/imgui_demo.cpp \
../imgui/imgui_draw.cpp \
../imgui/imgui_tables.cpp \
../imgui/imgui_widgets.cpp \
../imgui/backends/imgui_impl_sdl.cpp \
../imgui/backends/imgui_impl_opengl3.cpp \
../imgui/misc/cpp/imgui_stdlib.cpp \

OBJECTS = \
$(patsubst %.cc,%.o,$(SOURCES_CC)) \
$(patsubst %.cpp,%.o,$(SOURCES_IMGUI_CC)) \

DEPENDS = \
$(patsubst %.o,%.d,$(OBJECTS)) \

SOURCES = $(SOURCES_CC) $(SOURCES_H)

all: fraktaler-3 fraktaler-3.pdf index.html

#fraktaler-3-glfw: $(SOURCES) src/main_glfw.cc
#	g++ -std=c++20 -Wall -Wextra -pedantic -Og -march=native -fopenmp -o fraktaler-3-glfw $(SOURCES_CC) src/main_glfw.cc `pkg-config --cflags --libs glew glfw3 mpfr OpenEXR` -ggdb

clean:
	-rm $(OBJECTS)
	-rm $(DEPENDS)

fraktaler-3: $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LINK_FLAGS)

%.o: %.cc
	$(COMPILE) -o $@ -c $<

%.o: %.cpp
	$(COMPILE) -o $@ -c $<

fraktaler-3.pdf: README.md
	pandoc README.md -o fraktaler-3.pdf

index.html: README.md
	pandoc README.md -o index.html

-include \
$(DEPENDS) \
