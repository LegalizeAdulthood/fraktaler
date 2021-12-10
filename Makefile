# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

CFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -fopenmp -MMD
LIBS = glm mpfr OpenEXR
LIBS_GUI = glew sdl2

CFLAGS_IMGUI = -I../imgui -I../imgui/backends -I../imgui/misc/cpp
LIBS_IMGUI = -ldl

COMPILE_CLI = g++ $(CFLAGS) `pkg-config --cflags $(LIBS)`
COMPILE_GUI = g++ $(CFLAGS) `pkg-config --cflags $(LIBS) $(LIBS_GUI)` $(CFLAGS_IMGUI)

LINK = g++ $(CFLAGS)
LINK_FLAGS_CLI = `pkg-config --libs $(LIBS)`
LINK_FLAGS_GUI = `pkg-config --libs $(LIBS) $(LIBS_GUI)` $(LIBS_IMGUI)

EMSCRIPTEN=$(HOME)/opt/emscripten
COMPILE_WEB = em++ -std=c++20 -Wall -Wextra -pedantic -O3 -MMD -I$(EMSCRIPTEN)/include -s USE_SDL=2 $(CFLAGS_IMGUI) -s USE_PTHREADS
LINK_WEB = em++ -L$(EMSCRIPTEN)/lib
LINK_FLAGS_WEB = -lgmp -lmpfr -s USE_SDL=2 -s ALLOW_MEMORY_GROWTH=1 -s USE_PTHREADS -s PTHREAD_POOL_SIZE=4

SOURCES_CC = \
src/bla.cc \
src/colour.cc \
src/engine.cc \
src/formula.cc \
src/map.cc \
src/param.cc \
src/stats.cc \

SOURCES_H = \
src/bla.h \
src/colour.h \
src/colour_monochrome.h \
src/colour_rainbow.h \
src/complex.h \
src/display.h \
src/dual.h \
src/engine.h \
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

SOURCES_CLI_CC = \
src/cli.cc \
src/display_cpu.cc \

SOURCES_GUI_CC = \
src/display_gl.cc \
src/glutil.cc \
src/gui.cc \

SOURCES_WEB_CC = \
src/display_cpu.cc \
src/display_web.cc \
src/glutil.cc \
src/gui.cc \

SOURCES_IMGUI_CC = \
../imgui/imgui.cpp \
../imgui/imgui_demo.cpp \
../imgui/imgui_draw.cpp \
../imgui/imgui_tables.cpp \
../imgui/imgui_widgets.cpp \
../imgui/backends/imgui_impl_sdl.cpp \
../imgui/backends/imgui_impl_opengl3.cpp \
../imgui/misc/cpp/imgui_stdlib.cpp \

OBJECTS_CLI = \
$(patsubst %.cc,%.cli.o,$(SOURCES_CC)) \
$(patsubst %.cc,%.cli.o,$(SOURCES_CLI_CC)) \

OBJECTS_GUI = \
$(patsubst %.cc,%.gui.o,$(SOURCES_CC)) \
$(patsubst %.cc,%.gui.o,$(SOURCES_GUI_CC)) \
$(patsubst %.cpp,%.gui.o,$(SOURCES_IMGUI_CC)) \

OBJECTS_WEB = \
$(patsubst %.cc,%.web.o,$(SOURCES_CC)) \
$(patsubst %.cc,%.web.o,$(SOURCES_WEB_CC)) \
$(patsubst %.cpp,%.web.o,$(SOURCES_IMGUI_CC)) \

DEPENDS = \
$(patsubst %.o,%.d,$(OBJECTS_CLI)) \
$(patsubst %.o,%.d,$(OBJECTS_GUI)) \

all: cli gui web doc
cli: fraktaler-3-cli
gui: fraktaler-3-gui
web: live/fraktaler-3.html
doc: fraktaler-3.pdf index.html

clean:
	-rm $(OBJECTS_CLI)
	-rm $(OBJECTS_GUI)
	-rm $(OBJECTS_WEB)
	-rm $(DEPENDS)

fraktaler-3-cli: $(OBJECTS_CLI)
	$(LINK) -o $@ $(OBJECTS_CLI) $(LINK_FLAGS_CLI)

fraktaler-3-gui: $(OBJECTS_GUI)
	$(LINK) -o $@ $(OBJECTS_GUI) $(LINK_FLAGS_GUI)

live/fraktaler-3.html: $(OBJECTS_WEB)
	$(LINK_WEB) -o $@ $(OBJECTS_WEB) $(LINK_FLAGS_WEB)

%.cli.o: %.cc
	$(COMPILE_CLI) -o $@ -c $<

%.gui.o: %.cc
	$(COMPILE_GUI) -o $@ -c $<

%.cli.o: %.cpp
	$(COMPILE_CLI) -o $@ -c $<

%.gui.o: %.cpp
	$(COMPILE_GUI) -o $@ -c $<

%.web.o: %.cc
	$(COMPILE_WEB) -o $@ -c $<

%.web.o: %.cpp
	$(COMPILE_WEB) -o $@ -c $<

fraktaler-3.pdf: README.md
	pandoc README.md -o fraktaler-3.pdf

index.html: README.md
	pandoc README.md -o index.html

-include \
$(DEPENDS) \
