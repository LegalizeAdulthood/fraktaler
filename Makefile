# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

VERSION ?= $(shell test -d .git && git describe --always --dirty=+ || (cat VERSION.txt | head -n 1))
DATE ?= $(shell test -d .git && date --iso || (cat VERSION.txt | tail -n+1 | head -n 1))

SOURCE := $(shell cat INDEX.txt)
RELEASE := \
fraktaler-3-$(VERSION)-cli \
fraktaler-3-$(VERSION)-gui \
live/$(VERSION)/fraktaler-3.html \
fraktaler-3-$(VERSION).pdf \
fraktaler-3-$(VERSION).html.gz \
fraktaler-3-$(VERSION).css.gz \
fraktaler-3-$(VERSION).png \
fraktaler-3-$(VERSION).7z \

EMBEDSOURCE = -Wl,--format=binary -Wl,fraktaler-3-source.7z -Wl,--format=default

VERSIONS = \
-DFRAKTALER_3_VERSION_STRING="\"$(VERSION)\"" \
-DIMGUI_GIT_VERSION_STRING="\"$(shell cd ../imgui && git describe --always)\"" \
-DGLEW_VERSION_STRING="\"2.2.0\"" \

CFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -fopenmp -MMD
LIBS = glm mpfr OpenEXR
LIBS_GUI = glew sdl2

CFLAGS_IMGUI = -I../imgui -I../imgui/backends -I../imgui/misc/cpp
LIBS_IMGUI = -ldl

COMPILE_CLI = g++ $(CFLAGS) `pkg-config --cflags $(LIBS)` $(VERSIONS)
COMPILE_GUI = g++ $(CFLAGS) `pkg-config --cflags $(LIBS) $(LIBS_GUI)` $(CFLAGS_IMGUI) $(VERSIONS)

LINK = g++ $(CFLAGS)
LINK_FLAGS_CLI = `pkg-config --libs $(LIBS)`
LINK_FLAGS_GUI = `pkg-config --libs $(LIBS) $(LIBS_GUI)` $(LIBS_IMGUI)

EMSCRIPTEN=$(HOME)/opt/emscripten
COMPILE_WEB = em++ -std=c++20 -Wall -Wextra -pedantic -O3 -MMD $(CFLAGS_IMGUI) $(VERSIONS) -I$(EMSCRIPTEN)/include -s USE_SDL=2 -s USE_PTHREADS
LINK_WEB = em++ -L$(EMSCRIPTEN)/lib
LINK_FLAGS_WEB = -lgmp -lmpfr -s USE_SDL=2 -s ALLOW_MEMORY_GROWTH=1 -s USE_PTHREADS -s PTHREAD_POOL_SIZE=4

SOURCES_CC = \
src/bla.cc \
src/colour.cc \
src/engine.cc \
src/formula.cc \
src/map.cc \
src/param.cc \
src/source.cc \
src/stats.cc \
src/version.cc \

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

default: gui

cli: fraktaler-3-cli
gui: fraktaler-3-gui
web: live/latest/fraktaler-3.html

release: $(RELEASE)

clean:
	-rm $(OBJECTS_CLI)
	-rm $(OBJECTS_GUI)
	-rm $(OBJECTS_WEB)
	-rm $(DEPENDS)

VERSION.txt:
	echo "$(VERSION)" > VERSION.txt
	date --iso >> VERSION.txt
	touch -c -d '@0' VERSION.txt

# distribution

fraktaler-3-$(VERSION)-cli: fraktaler-3-cli
	cp -avf $< $@
	strip --strip-unneeded $@

fraktaler-3-$(VERSION)-gui: fraktaler-3-gui
	cp -avf $< $@
	strip --strip-unneeded $@

fraktaler-3-$(VERSION).7z: fraktaler-3-source.7z
	cp -avf $< $@

fraktaler-3-source.7z: $(SOURCE)
	-rm -rf fraktaler-3-source.7z "fraktaler-3-$(VERSION)"
	mkdir fraktaler-3-$(VERSION)
	cat INDEX.txt |	cpio -pdv "fraktaler-3-$(VERSION)"
	7zr a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on "fraktaler-3-source.7z" "fraktaler-3-$(VERSION)/"

fraktaler-3-$(VERSION).html: README.md fraktaler-3-$(VERSION).css
	pandoc README.md --metadata="title=fraktaler-3-$(VERSION)" --metadata="date=$(DATE)" --standalone -c "fraktaler-3-$(VERSION).css" --toc -o "fraktaler-3-$(VERSION).html"
	sed -i "s/src=.fraktaler-3.png./src='fraktaler-3-$(VERSION).png'/g" "fraktaler-3-$(VERSION).html"
	sed -i "s/href=.fraktaler-3.css./href='fraktaler-3-$(VERSION).css'/g" "fraktaler-3-$(VERSION).html"

fraktaler-3-$(VERSION).css: fraktaler-3.css
	cp -avf $< $@

fraktaler-3-$(VERSION).png: fraktaler-3.png
	cp -avf $< $@

fraktaler-3-$(VERSION).pdf: README.md fraktaler-3.png
	pandoc README.md --metadata="title=fraktaler-3-$(VERSION)" --metadata="date=$(DATE)" -o fraktaler-3-$(VERSION).pdf

# link

fraktaler-3-cli: $(OBJECTS_CLI)
	$(LINK) -o $@ $(OBJECTS_CLI) $(LINK_FLAGS_CLI) $(EMBEDSOURCE)

fraktaler-3-gui: $(OBJECTS_GUI)
	$(LINK) -o $@ $(OBJECTS_GUI) $(LINK_FLAGS_GUI) $(EMBEDSOURCE)

live/$(VERSION)/fraktaler-3.html: $(OBJECTS_WEB)
	mkdir -p live/$(VERSION)
	cp -avi live/latest/index.html live/$(VERSION)
	$(LINK_WEB) -o $@ $(OBJECTS_WEB) $(LINK_FLAGS_WEB)
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.js
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.wasm
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.worker.js
	gzip -9 -k -f live/$(VERSION)/index.html

%.gz: %
	gzip -k -9 -f $<

# compile

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

# dependencies

.PHONY: default release clean VERSION.txt

-include \
$(DEPENDS) \
