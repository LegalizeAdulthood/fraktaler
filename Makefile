# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

VERSION ?= $(shell test -d .git && git describe --always --dirty=+ || (cat VERSION.txt | head -n 1))
DATE ?= $(shell test -d .git && date --iso || (cat VERSION.txt | tail -n+1 | head -n 1))

SYSTEM ?= native-gcc
COMPILER = false
OEXT = .o
include build/$(SYSTEM).mk

STRIP ?= strip

SOURCE := $(shell cat INDEX.txt)

VERSIONS += \
-DFRAKTALER_3_VERSION_STRING="\"$(VERSION)\"" \
-DIMGUI_GIT_VERSION_STRING="\"$(shell test -d ../imgui && cd ../imgui && git describe --always || echo none)\"" \
-DGLEW_VERSION_STRING="\"2.1.0\"" \

LIBS += glm mpfr OpenEXR zlib
LIBS_GUI += sdl2

CFLAGS_IMGUI += -I../imgui -I../imgui/backends -I../imgui/misc/cpp -DHAVE_GUI
LIBS_IMGUI +=

CFLAGS += -ggdb
LDFLAGS += -ggdb

COMPILE_CLI = $(COMPILER) $(CPPFLAGS) $(CFLAGS) `PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config $(PKG_CONFIG_FLAGS) --cflags $(LIBS) | sed "$(PKG_CONFIG_SED)"` $(VERSIONS)
COMPILE_GUI = $(COMPILER) $(CPPFLAGS) $(CFLAGS) `PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config $(PKG_CONFIG_FLAGS) --cflags $(LIBS) $(LIBS_GUI) | sed "$(PKG_CONFIG_SED)"` $(CFLAGS_IMGUI) $(VERSIONS)
COMPILE_WEB = $(COMPILER) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_IMGUI) $(VERSIONS)

LINK = $(COMPILER) $(CFLAGS)
LINK_FLAGS_CLI = $(LDFLAGS) `PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config $(PKG_CONFIG_FLAGS) --libs $(LIBS) | sed "$(PKG_CONFIG_SED)"`
LINK_FLAGS_GUI = $(LDFLAGS) `PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config $(PKG_CONFIG_FLAGS) --libs $(LIBS) $(LIBS_GUI) | sed "$(PKG_CONFIG_SED)"` $(LIBS_IMGUI)

SOURCES_CC += \
src/bla.cc \
src/colour.cc \
src/engine.cc \
src/formula.cc \
src/map.cc \
src/param.cc \
src/source.cc \
src/stats.cc \
src/version.cc \

SOURCES_CLI_CC += \
src/cli.cc \
src/display_cpu.cc \

SOURCES_GUI_CC += \
src/display_gl.cc \
src/glutil.cc \
src/gui.cc \

SOURCES_GUI_C += \

SOURCES_WEB_CC += \
src/display_cpu.cc \
src/display_web.cc \
src/glutil.cc \
src/gui.cc \

SOURCES_IMGUI_CC += \
../imgui/imgui.cpp \
../imgui/imgui_draw.cpp \
../imgui/imgui_tables.cpp \
../imgui/imgui_widgets.cpp \
../imgui/backends/imgui_impl_sdl.cpp \
../imgui/backends/imgui_impl_opengl3.cpp \
../imgui/misc/cpp/imgui_stdlib.cpp \

OBJECTS_CLI = \
$(patsubst %.cc,%.cli$(OEXT),$(SOURCES_CC)) \
$(patsubst %.cc,%.cli$(OEXT),$(SOURCES_CLI_CC)) \

OBJECTS_GUI = \
$(patsubst %.cc,%.gui$(OEXT),$(SOURCES_CC)) \
$(patsubst %.cc,%.gui$(OEXT),$(SOURCES_GUI_CC)) \
$(patsubst %.cpp,%.gui$(OEXT),$(SOURCES_IMGUI_CC)) \
$(patsubst %.c,%.gui$(OEXT),$(SOURCES_GUI_C)) \

OBJECTS_WEB = \
$(patsubst %.cc,%.web$(OEXT),$(SOURCES_CC)) \
$(patsubst %.cc,%.web$(OEXT),$(SOURCES_WEB_CC)) \
$(patsubst %.cpp,%.web$(OEXT),$(SOURCES_IMGUI_CC)) \

DEPENDS = \
$(patsubst %.o,%.d,$(OBJECTS_CLI)) \
$(patsubst %.o,%.d,$(OBJECTS_GUI)) \

default: $(TARGETS)

cli: fraktaler-3-$(VERSION)-cli$(EXEEXT)
gui: fraktaler-3-$(VERSION)-gui$(EXEEXT)
web: live/$(VERSION)/index.html

clean:
	-rm -f src/fraktaler-3-source.7z.h
	-rm -f $(OBJECTS_CLI)
	-rm -f $(OBJECTS_GUI)
	-rm -f $(OBJECTS_WEB)
	-rm -f $(DEPENDS)

VERSION.txt:
	echo "$(VERSION)" > VERSION.txt
	date --iso >> VERSION.txt
	touch -c -d '@0' VERSION.txt

# distribution

fraktaler-3-$(VERSION)-cli$(EXEEXT): fraktaler-3-cli$(EXEEXT)
	cp -avf $< $@
	$(STRIP) --strip-unneeded $@

fraktaler-3-$(VERSION)-gui$(EXEEXT): fraktaler-3-gui$(EXEEXT)
	cp -avf $< $@
	$(STRIP) --strip-unneeded $@

fraktaler-3-$(VERSION).7z: fraktaler-3-source.7z
	cp -avf $< $@

fraktaler-3-source.7z: $(SOURCE)
	-rm -rf fraktaler-3-source.7z "fraktaler-3-$(VERSION)"
	mkdir fraktaler-3-$(VERSION)
	cat INDEX.txt |	cpio -pdv "fraktaler-3-$(VERSION)"
	7zr a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on "fraktaler-3-source.7z" "fraktaler-3-$(VERSION)/"

fraktaler-3-$(VERSION).html: README.md fraktaler-3-$(VERSION).css
	pandoc README.md --metadata="title=fraktaler-3-$(VERSION)" --metadata="date=$(DATE)" --standalone -c "fraktaler-3-$(VERSION).css" --toc -o "fraktaler-3-$(VERSION).html"
	sed -i "s/<head>/<head profile='http://www.w3.org/2005/10/profile'>/g" "fraktaler-3-$(VERSION).html"
	sed -i "s/src=.fraktaler-3.png./src='fraktaler-3-$(VERSION).png'/g" "fraktaler-3-$(VERSION).html"
	sed -i "s/href=.fraktaler-3.css./href='fraktaler-3-$(VERSION).css'/g" "fraktaler-3-$(VERSION).html"

fraktaler-3-$(VERSION).css: fraktaler-3.css
	cp -avf $< $@

fraktaler-3-$(VERSION).png: fraktaler-3.png
	cp -avf $< $@

fraktaler-3-$(VERSION).pdf: README.md fraktaler-3.png
	pandoc README.md --metadata="title=fraktaler-3-$(VERSION)" --metadata="date=$(DATE)" -o fraktaler-3-$(VERSION).pdf

LICENSE.pdf: LICENSE.md
	pandoc LICENSE.md --toc --toc-depth=4 -V geometry="margin=1in" --metadata="author=Free Software Foundation, Inc." --metadata="title=GNU Affero General Public License" --metadata="date=19 November 2007" -o LICENSE.pdf

src/fraktaler-3-source.7z.h: fraktaler-3-source.7z
	xxd -i $< | sed "s/unsigned/const unsigned/g" > $@

# link

fraktaler-3-cli$(EXEEXT): $(OBJECTS_CLI)
	$(LINK) -o $@ $(OBJECTS_CLI) $(LINK_FLAGS_CLI)

fraktaler-3-gui$(EXEEXT): $(OBJECTS_GUI)
	$(LINK) -o $@ $(OBJECTS_GUI) $(LINK_FLAGS_GUI)

live/$(VERSION)/index.html: $(OBJECTS_WEB) fraktaler-3-$(VERSION).7z fraktaler-3.ico
	mkdir -p live/$(VERSION)
	cp -avf src/index.html fraktaler-3-$(VERSION).7z fraktaler-3.ico live/$(VERSION)
	sed -i "s/VERSION/$(VERSION)/g" "live/$(VERSION)/index.html"
	$(LINK) -o live/$(VERSION)/fraktaler-3.html $(OBJECTS_WEB) $(LDFLAGS)
	rm live/$(VERSION)/fraktaler-3.html
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.ico
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.js
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.wasm
	gzip -9 -k -f live/$(VERSION)/fraktaler-3.worker.js
	gzip -9 -k -f live/$(VERSION)/index.html

%.gz: %
	gzip -k -9 -f $<

# compile

%.cli$(OEXT): %.cc
	$(COMPILE_CLI) -o $@ -c $<

%.cli$(OEXT): %.cpp
	$(COMPILE_CLI) -o $@ -c $<

%.cli$(OEXT): %.c
	$(COMPILE_CLI) -o $@ -c $<

%.gui$(OEXT): %.cc
	$(COMPILE_GUI) -o $@ -c $<

%.gui$(OEXT): %.cpp
	$(COMPILE_GUI) -o $@ -c $<

%.gui$(OEXT): %.c
	$(COMPILE_GUI) -o $@ -c $<

%.web$(OEXT): %.cc
	$(COMPILE_WEB) -o $@ -c $<

%.web$(OEXT): %.cpp
	$(COMPILE_WEB) -o $@ -c $<

release:
	mkdir -p fraktaler-3-$(VERSION)-windows
	cp -avft fraktaler-3-$(VERSION)-windows fraktaler-3-$(VERSION)-*.exe fraktaler-3-$(VERSION).pdf LICENSE.pdf
	7zr a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on fraktaler-3-$(VERSION)-windows.7z fraktaler-3-$(VERSION)-windows/*

# dependencies

.PHONY: default clean VERSION.txt
.SUFFIXES:

source.cli$(OEXT): src/fraktaler-3-source.7z.h
source.gui$(OEXT): src/fraktaler-3-source.7z.h
source.web$(OEXT): src/fraktaler-3-source.7z.h

-include \
$(DEPENDS) \
