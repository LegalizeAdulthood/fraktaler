#!/bin/bash
# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

make src/fraktaler-3-source.7z.cc
cd ..
wget -c https://www.libsdl.org/release/SDL2-2.0.18.tar.gz
tar xaf SDL2-2.0.18.tar.gz
cd SDL2-2.0.18/build-scripts
./androidbuild.sh uk.co.mathr.fraktaler-3 ../../fraktaler-3/src/main.cc
cd ../build/uk.co.mathr.fraktaler-3/app/jni
rm -r src
ln -s ../../../../../fraktaler-3/src/
cd ../..
./gradlew installDebug
