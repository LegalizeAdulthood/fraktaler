#!/bin/bash
# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021,2022 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only

TOP="$(pwd)"
make src/fraktaler-3-source.7z.h
VERSION="$(cat VERSION.txt | head -n 1)"
if [[ "$1" =~ "prepare" ]]
then
mkdir -p "${TOP}/android/src"
cd "${TOP}/src"
ln -fs ../android/arm64-v8a/
ln -fs ../android/armeabi-v7a/
ln -fs ../android/x86/
ln -fs ../android/x86_64/
ln -fs ../../imgui/
ln -fs ../../imgui-filebrowser/
ln -fs ../../toml11/
ln -fs ../../fraktaler-3/
cd "${TOP}/android/src"
git clone https://code.mathr.co.uk/android-build-scripts.git
wget -c https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
wget -c https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.xz
wget -c https://github.com/advanpix/mpreal/archive/refs/tags/mpfrc++-3.6.8.tar.gz
wget -c https://github.com/g-truc/glm/releases/download/0.9.9.8/glm-0.9.9.8.7z
wget -c https://www.libsdl.org/release/SDL2-2.0.18.tar.gz
tar xaf gmp-6.2.1.tar.lz
cd gmp-6.2.1
NDK="${ANDROID_NDK_HOME}" DESTBASE="${TOP}/android" "${TOP}/android/src/android-build-scripts/build/compile-gmp-x86.sh"
NDK="${ANDROID_NDK_HOME}" DESTBASE="${TOP}/android" "${TOP}/android/src/android-build-scripts/build/compile-gmp-arm.sh"
#NDK="${ANDROID_NDK_HOME}" DESTBASE="${TOP}/android" "${TOP}/android/src/android-build-scripts/build/compile-gmp-mips.sh"
cd ..
tar xaf mpfr-4.1.0.tar.xz
cd mpfr-4.1.0
NDK="${ANDROID_NDK_HOME}" DESTBASE="${TOP}/android" "${TOP}/android/src/android-build-scripts/build/compile-mpfr-x86.sh"
NDK="${ANDROID_NDK_HOME}" DESTBASE="${TOP}/android" "${TOP}/android/src/android-build-scripts/build/compile-mpfr-arm.sh"
#NDK="${ANDROID_NDK_HOME}" DESTBASE="${TOP}/android" "${TOP}/android/src/android-build-scripts/build/compile-mpfr-mips.sh"
cd ..
tar xaf mpfrc++-3.6.8.tar.gz
cp -avf mpreal-mpfrc-3.6.8/mpreal.h "${TOP}/android/armeabi-v7a/include"
cp -avf mpreal-mpfrc-3.6.8/mpreal.h "${TOP}/android/arm64-v8a/include"
cp -avf mpreal-mpfrc-3.6.8/mpreal.h "${TOP}/android/x86/include"
cp -avf mpreal-mpfrc-3.6.8/mpreal.h "${TOP}/android/x86_64/include"
7zr x glm-0.9.9.8.7z
ln -fs "${TOP}/android/src/glm/glm" "${TOP}/android/armeabi-v7a/include/glm"
ln -fs "${TOP}/android/src/glm/glm" "${TOP}/android/arm64-v8a/include/glm"
ln -fs "${TOP}/android/src/glm/glm" "${TOP}/android/x86/include/glm"
ln -fs "${TOP}/android/src/glm/glm" "${TOP}/android/x86_64/include/glm"
tar xaf SDL2-2.0.18.tar.gz
cd SDL2-2.0.18/build-scripts
./androidbuild.sh uk.co.mathr.fraktaler.v3 ../../../../src/main.cc
cd ../build/uk.co.mathr.fraktaler.v3/
rm -rf app
ln -fs "${TOP}/app" app
mkdir -p app/jni/SDL
ln -fs "${TOP}/src" app/jni/src
ln -fs "${TOP}/android/src/SDL2-2.0.18/include" app/jni/SDL/include
ln -fs "${TOP}/android/src/SDL2-2.0.18/src" app/jni/SDL/src
ln -fs "${TOP}/android/src/SDL2-2.0.18/android-project/app/src/main/java/org" app/src/main/java/org
else
cd "${TOP}/android/src/SDL2-2.0.18/build/uk.co.mathr.fraktaler.v3"
sed "s|VERSION|${VERSION}|g" < app/build.gradle.in > app/build.gradle
if [[ "$1" =~ "release" ]]
then
  ./gradlew assembleRelease
  cd app/build/outputs/apk/release/
  zipalign -v -p 4 app-release-unsigned.apk app-release-unsigned-aligned.apk
  apksigner sign --ks ~/.fraktaler-3.ks --out "uk.co.mathr.fraktaler.v3-${VERSION}.apk" app-release-unsigned-aligned.apk
  adb install -r "uk.co.mathr.fraktaler.v3-${VERSION}.apk"
  cp -avi "uk.co.mathr.fraktaler.v3-${VERSION}.apk" "${TOP}"
else
  ./gradlew installDebug
fi
fi
