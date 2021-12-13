#!/bin/bash
# Fraktaler 3 -- fast deep escape time fractals
# Copyright (C) 2021 Claude Heiland-Allen
# SPDX-License-Identifier: AGPL-3.0-only
set -e
NCPUS="$(( $(nproc) * 2 ))"
export CPPFLAGS="-D__USE_MINGW_ANSI_STDIO=1 -DWINVER=0x501 -D_WIN32_WINNT=0x501"
export LDFLAGS="-static-libgcc -static-libstdc++ -static -Wl,-Bstatic -lstdc++ -Wl,-Bdynamic"
ALL_ARCH="x86_64 i686 aarch64 armv7 emscripten"
ALL_LIBS="gmp mpfr mpreal zlib glm openexr clew sdl2"
THREADMODEL="posix"
if [[ "x$1" = "x" ]]
then
  ACTION="dl ${ALL_ARCH}"
else
  ACTION="$1"
fi
if [[ "x$2" = "x" ]]
then
  PREPARE="${ALL_LIBS}"
else
  PREPARE="$2"
fi
if [[ "x$3" = "x" ]]
then
  COMPILER="gcc"
else
  COMPILER="$3"
fi
if [[ "${ACTION}" =~ "-h" ]]
then
  echo "usage:"
  echo "  $0             # download and build everything"
  echo "  $0 dl          # download sources"
  echo "  $0 \$arch       # build all libraries for one architecture"
  echo "    # supported architectures:"
  echo "    ${ALL_ARCH}"
  echo "  $0 \$arch \$lib  # build one library for one architecture"
  echo "    # supported libraries:"
  echo "    ${ALL_LIBS}"
  echo "  $0 \$arch \$lib \$compiler # build with a specific compiler"
  echo "    # supported compilers:"
  echo "    gcc llvm"
  exit 0
fi
if [[ "${ACTION}" =~ "dl" ]]
then
  mkdir -p ~/win/src
  #cp -avft ~/win/src *.patch
  # download
  cd ~/win/src
  wget -c https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
  wget -c https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.xz
  #wget -c https://www.mpfr.org/mpfr-current/allpatches
  wget -c https://github.com/advanpix/mpreal/archive/refs/tags/mpfrc++-3.6.8.tar.gz
  wget -c https://zlib.net/zlib-1.2.11.tar.xz
  wget -c https://jpegclub.org/support/files/jpegsrc.v6b2.tar.gz
  #wget -c https://download.sourceforge.net/libpng/libpng-1.6.37.tar.xz
  wget -c https://download.osgeo.org/libtiff/tiff-4.3.0.tar.gz
  wget -c https://github.com/g-truc/glm/releases/download/0.9.9.8/glm-0.9.9.8.7z
  wget -c https://github.com/AcademySoftwareFoundation/openexr/archive/refs/tags/v2.5.7.tar.gz -O openexr-2.5.7.tar.gz
  wget -c https://libsdl.org/release/SDL2-2.0.18.tar.gz
  git clone https://github.com/meganz/mingw-std-threads.git || ( cd mingw-std-threads && git pull )
  git clone https://github.com/martijnberger/clew.git || ( cd clew && git pull )
fi
if [[ "${ACTION}" =~ "x86_64" ]]
then
  if [[ "${PREPARE}" =~ "gmp" ]]
  then
    # gmp 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/gmp-6.2.1.tar.lz
    cd gmp-6.2.1/
    CC_FOR_BUILD="gcc" CPP_FOR_BUILD="gcc -E" ./configure --build=x86_64-pc-linux-gnu --host=x86_64-w64-mingw32 --enable-fat --prefix=$HOME/win/${THREADMODEL}/x86_64
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpfr" ]]
  then
    # mpfr 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/mpfr-4.1.0.tar.xz
    cd mpfr-4.1.0/
    #patch -N -Z -p1 < ~/win/src/allpatches
    ./configure --host=x86_64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/x86_64 --with-gmp-build=../gmp-6.2.1 --enable-static --disable-shared
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpreal" ]]
  then
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/mpfrc++-3.6.8.tar.gz
    cp -avf mpreal-mpfrc-3.6.8/mpreal.h ~/win/${THREADMODEL}/x86_64/include
  fi
  if [[ "${PREPARE}" =~ "zlib" ]]
  then
    # zlib 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/zlib-1.2.11.tar.xz
    cd zlib-1.2.11/
    CC=x86_64-w64-mingw32-gcc ./configure --static --prefix=$HOME/win/${THREADMODEL}/x86_64
    CC=x86_64-w64-mingw32-gcc make -j $NCPUS
    CC=x86_64-w64-mingw32-gcc make install
  fi
  if [[ "${PREPARE}" =~ "png" ]]
  then
    # png 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/libpng-1.6.37.tar.xz
    cd libpng-1.6.37/
    ./configure --disable-shared --host=x86_64-w64-mingw32 CPPFLAGS="$CPPFLAGS -I$HOME/win/${THREADMODEL}/x86_64/include" LDFLAGS="$LDFLAGS -L$HOME/win/${THREADMODEL}/x86_64/lib" --prefix=$HOME/win/${THREADMODEL}/x86_64
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "jpeg" ]]
  then
    # jpeg 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/jpegsrc.v6b2.tar.gz
    cd jpeg-6b2/
    ./configure --disable-shared --host=x86_64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/x86_64
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "tiff" ]]
  then
    # tiff 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/tiff-4.3.0.tar.gz
    cd tiff-4.3.0/
    ./configure --disable-shared --host=x86_64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/x86_64
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "gsl" ]]
  then
    # gsl 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/gsl-2.7.tar.gz
    cd gsl-2.7/
    ./configure --disable-shared --host=x86_64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/x86_64
    make -j $NCPUS
    make install
    rm -f gsl-histogram
    ln -s gsl-histogram.exe gsl-histogram # hack for test suite
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "pixman" ]]
  then
    # pixman 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    rm -rf pixman
    cp -avf ~/win/src/pixman pixman
    cd pixman/
    NOCONFIGURE=true ./autogen.sh
    CC=x86_64-w64-mingw32-gcc LDFLAGS=-L$HOME/win/${THREADMODEL}/x86_64/lib ./configure --disable-shared --disable-openmp --prefix=$HOME/win/${THREADMODEL}/x86_64
    make SUBDIRS="pixman" -j $NCPUS
    make SUBDIRS="pixman" install
    #make SUBDIRS="pixman test" check -k || echo "expected 1 FAIL (thread-test)"
  fi
  if [[ "${PREPARE}" =~ "boost" ]]
  then
    # boost 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    7zr x ~/win/src/boost_1_76_0.7z
    cd ~/win/${THREADMODEL}/x86_64/include
    rm -f boost
    ln -s ../src/boost_1_76_0/boost/
  fi
  if [[ "${PREPARE}" =~ "glm" ]]
  then
    # glm 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    7zr x ~/win/src/glm-0.9.9.8.7z
    cd ~/win/${THREADMODEL}/x86_64/include
    rm -f glm
    ln -s ../src/glm/glm/
  fi
  if [[ "${PREPARE}" =~ "mingw-std-threads" ]]
  then
    # mingw-std-threads 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/include
    cd ~/win/${THREADMODEL}/x86_64/include
    rm -f mingw-std-threads
    ln -s ~/win/src/mingw-std-threads
  fi
  if [[ "${PREPARE}" =~ "openexr" ]]
  then
    # openexr 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/openexr-2.5.7.tar.gz
    cd openexr-2.5.7/
    #if [[ "x${COMPILER}" = "xgcc" ]]
    #then
    #  patch -p1 < ~/win/src/openexr-2.4.0.patch
    #fi
    sed -i "s/#ifdef _WIN32/#if 0/g" OpenEXR/IlmImf/ImfStdIO.cpp
    mkdir -p build
    cd build
    cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain-mingw.cmake -DCMAKE_CXX_FLAGS=-I$HOME/win/${THREADMODEL}/x86_64/include -DZLIB_INCLUDE_DIR=$HOME/win/${THREADMODEL}/x86_64/include -DZLIB_LIBRARY=$HOME/win/${THREADMODEL}/x86_64/lib/libz.a -DCMAKE_INSTALL_PREFIX=$HOME/win/${THREADMODEL}/x86_64 ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "glfw" ]]
  then
    # glfw 64
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    unzip ~/win/src/glfw-3.3.4.zip
    cd glfw-3.3.4/
    mkdir -p build
    cd build
    cmake -DCMAKE_TOOLCHAIN_FILE=../CMake/x86_64-w64-mingw32.cmake -DCMAKE_INSTALL_PREFIX=${HOME}/win/${THREADMODEL}/x86_64/ -DBUILD_SHARED_LIBS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF -DGLFW_BUILD_DOCS=OFF -DGLFW_USE_HYBRID_HPG=ON ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "sdl2" ]]
  then
    # sdl2 32
    mkdir -p ~/win/${THREADMODEL}/x86_64/src
    cd ~/win/${THREADMODEL}/x86_64/src
    tar xaf ~/win/src/SDL2-2.0.18.tar.gz
    cd SDL2-2.0.18/
    ./configure --prefix=${HOME}/win/${THREADMODEL}/x86_64 --host=x86_64-w64-mingw32
    make -j $NCPUS
    make install
  fi
  # clew 64
  # nop
fi
if [[ "${ACTION}" =~ "i686" ]]
then
  if [[ "${PREPARE}" =~ "gmp" ]]
  then
    # gmp 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/gmp-6.2.1.tar.lz
    cd gmp-6.2.1/
    CC_FOR_BUILD="gcc" CPP_FOR_BUILD="gcc -E" ./configure --build=x86_64-pc-linux-gnu --host=i686-w64-mingw32 --enable-fat --prefix=$HOME/win/${THREADMODEL}/i686
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpfr" ]]
  then
    # mpfr 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/mpfr-4.1.0.tar.xz
    cd mpfr-4.1.0/
    #patch -N -Z -p1 < ~/win/src/allpatches
    ./configure --host=i686-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/i686 --with-gmp-build=../gmp-6.2.1 --enable-static --disable-shared
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpreal" ]]
  then
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/mpfrc++-3.6.8.tar.gz
    cp -avf mpreal-mpfrc-3.6.8/mpreal.h ~/win/${THREADMODEL}/i686/include
  fi
  if [[ "${PREPARE}" =~ "zlib" ]]
  then
    # zlib 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/zlib-1.2.11.tar.xz
    cd zlib-1.2.11/
    CC=i686-w64-mingw32-gcc ./configure --static --prefix=$HOME/win/${THREADMODEL}/i686
    CC=i686-w64-mingw32-gcc make -j $NCPUS
    CC=i686-w64-mingw32-gcc make install
  fi
  if [[ "${PREPARE}" =~ "png" ]]
  then
    # png 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/libpng-1.6.37.tar.xz
    cd libpng-1.6.37/
    ./configure --disable-shared --host=i686-w64-mingw32 CPPFLAGS="$CPPFLAGS -I$HOME/win/${THREADMODEL}/i686/include" LDFLAGS="$LDFLAGS -L$HOME/win/${THREADMODEL}/i686/lib" --prefix=$HOME/win/${THREADMODEL}/i686
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "jpeg" ]]
  then
    # jpeg 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/jpegsrc.v6b2.tar.gz
    cd jpeg-6b2/
    ./configure --disable-shared --host=i686-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/i686
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "tiff" ]]
  then
    # tiff 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/tiff-4.3.0.tar.gz
    cd tiff-4.3.0/
    ./configure --disable-shared --host=i686-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/i686
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "gsl" ]]
  then
    # gsl 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/gsl-2.7.tar.gz
    cd gsl-2.7/
    ./configure --disable-shared --host=i686-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/i686
    make -j $NCPUS
    make install
    rm -f gsl-histogram
    ln -s gsl-histogram.exe gsl-histogram # hack for test suite
    make -k check || echo "expected 2 FAILs (multifit_nlinear, multilarge_nlinear)"
  fi
  if [[ "${PREPARE}" =~ "pixman" ]]
  then
    # pixman 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    rm -rf pixman
    cp -avf ~/win/src/pixman pixman
    cd pixman/
    #patch -p1 < ~/win/src/pixman-0.38.4-i686.patch
    NOCONFIGURE=true ./autogen.sh
    CC=i686-w64-mingw32-gcc LDFLAGS=-L$HOME/win/${THREADMODEL}/i686/lib ./configure --disable-shared --disable-openmp --prefix=$HOME/win/${THREADMODEL}/i686
    make SUBDIRS="pixman" -j $NCPUS
    make SUBDIRS="pixman" install
    #make SUBDIRS="pixman test" check -k || echo "expected 1 FAIL (thread-test)"
  fi
  if [[ "${PREPARE}" =~ "boost" ]]
  then
    # boost 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    7zr x ~/win/src/boost_1_76_0.7z
    cd ~/win/${THREADMODEL}/i686/include
    rm -f boost
    ln -s ../src/boost_1_76_0/boost/
  fi
  if [[ "${PREPARE}" =~ "glm" ]]
  then
    # glm 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    7zr x ~/win/src/glm-0.9.9.8.7z
    cd ~/win/${THREADMODEL}/i686/include
    rm -f glm
    ln -s ../src/glm/glm/
  fi
  if [[ "${PREPARE}" =~ "mingw-std-threads" ]]
  then
    # mingw-std-threads 32
    mkdir -p ~/win/${THREADMODEL}/i686/include
    cd ~/win/${THREADMODEL}/i686/include
    rm -f mingw-std-threads
    ln -s ~/win/src/mingw-std-threads
  fi
  if [[ "${PREPARE}" =~ "openexr" ]]
  then
    # openexr 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xf ~/win/src/openexr-2.5.7.tar.gz
    cd openexr-2.5.7/
    #if [[ "x${COMPILER}" = "xgcc" ]]
    #then
    #  patch -p1 < ~/win/src/openexr-2.4.0.patch
    #fi
    sed -i "s/#ifdef _WIN32/#if 0/g" OpenEXR/IlmImf/ImfStdIO.cpp
    sed -i "s/x86_64/i686/g" cmake/Toolchain-mingw.cmake
    mkdir -p build
    cd build
    cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain-mingw.cmake -DCMAKE_CXX_FLAGS=-I$HOME/win/${THREADMODEL}/i686/include -DZLIB_INCLUDE_DIR=$HOME/win/${THREADMODEL}/i686/include -DZLIB_LIBRARY=$HOME/win/${THREADMODEL}/i686/lib/libz.a -DCMAKE_INSTALL_PREFIX=$HOME/win/${THREADMODEL}/i686 ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "glfw" ]]
  then
    # glfw 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    unzip ~/win/src/glfw-3.3.4.zip
    cd glfw-3.3.4/
    mkdir -p build
    cd build
    cmake -DCMAKE_TOOLCHAIN_FILE=../CMake/i686-w64-mingw32.cmake -DCMAKE_INSTALL_PREFIX=${HOME}/win/${THREADMODEL}/i686/ -DBUILD_SHARED_LIBS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF -DGLFW_BUILD_DOCS=OFF -DGLFW_USE_HYBRID_HPG=ON ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "sdl2" ]]
  then
    # sdl2 32
    mkdir -p ~/win/${THREADMODEL}/i686/src
    cd ~/win/${THREADMODEL}/i686/src
    tar xaf ~/win/src/SDL2-2.0.18.tar.gz
    cd SDL2-2.0.18/
    ./configure --prefix=${HOME}/win/${THREADMODEL}/i686 --host=i686-w64-mingw32
    make -j $NCPUS
    make install
  fi
  # clew 32
  #nop
fi
if [[ "${ACTION}" =~ "aarch64" ]]
then
  if [[ "${PREPARE}" =~ "gmp" ]]
  then
    # gmp 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/gmp-6.2.1.tar.lz
    cd gmp-6.2.1/
    CC_FOR_BUILD="gcc" CPP_FOR_BUILD="gcc -E" ./configure --build=aarch64-pc-linux-gnu --host=aarch64-w64-mingw32 --enable-fat --prefix=$HOME/win/${THREADMODEL}/aarch64
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpfr" ]]
  then
    # mpfr 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/mpfr-4.1.0.tar.xz
    cd mpfr-4.1.0/
    #patch -N -Z -p1 < ~/win/src/allpatches
    ./configure --host=aarch64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/aarch64 --with-gmp-build=../gmp-6.2.1 --enable-static --disable-shared
    make -j $NCPUS
    make install
    # make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpreal" ]]
  then
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/mpfrc++-3.6.8.tar.gz
    cp -avf mpreal-mpfrc-3.6.8/mpreal.h ~/win/${THREADMODEL}/aarch64/include
  fi
  if [[ "${PREPARE}" =~ "zlib" ]]
  then
    # zlib 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/zlib-1.2.11.tar.xz
    cd zlib-1.2.11/
    CC=aarch64-w64-mingw32-gcc AR=aarch64-w64-mingw32-ar RANLIB=aarch64-w64-mingw32-ranlib ./configure --static --prefix=$HOME/win/${THREADMODEL}/aarch64
    CC=aarch64-w64-mingw32-gcc AR=aarch64-w64-mingw32-ar RANLIB=aarch64-w64-mingw32-ranlib make -j $NCPUS -k || echo
    CC=aarch64-w64-mingw32-gcc AR=aarch64-w64-mingw32-ar RANLIB=aarch64-w64-mingw32-ranlib make install
  fi
  if [[ "${PREPARE}" =~ "png" ]]
  then
    # png 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/libpng-1.6.37.tar.xz
    cd libpng-1.6.37/
    ./configure --disable-shared --host=aarch64-w64-mingw32 CPPFLAGS="$CPPFLAGS -I$HOME/win/${THREADMODEL}/aarch64/include" LDFLAGS="$LDFLAGS -L$HOME/win/${THREADMODEL}/aarch64/lib" --prefix=$HOME/win/${THREADMODEL}/aarch64
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "jpeg" ]]
  then
    # jpeg 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/libpng-1.6.37.tar.xz
    tar xaf ~/win/src/jpegsrc.v6b2.tar.gz
    cp -avf libpng-1.6.37/config.sub jpeg-6b2/config.sub
    cd jpeg-6b2/
    ./configure --disable-shared --host=aarch64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/aarch64
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "tiff" ]]
  then
    # tiff 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/tiff-4.3.0.tar.gz
    cd tiff-4.3.0/
    ./configure --disable-shared --host=aarch64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/aarch64
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "gsl" ]]
  then
    # gsl 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/gsl-2.7.tar.gz
    cd gsl-2.7/
    ./configure --disable-shared --host=aarch64-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/aarch64
    make -j $NCPUS
    make install
    rm -f gsl-histogram
    ln -s gsl-histogram.exe gsl-histogram # hack for test suite
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "pixman" ]]
  then
    # pixman 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    rm -rf pixman
    cp -avf ~/win/src/pixman pixman
    cd pixman/
    NOCONFIGURE=true ./autogen.sh
    LDFLAGS=-L$HOME/win/${THREADMODEL}/aarch64/lib ./configure --host aarch64-w64-mingw32 --disable-shared --disable-openmp --prefix=$HOME/win/${THREADMODEL}/aarch64
    make SUBDIRS="pixman" -j $NCPUS
    make SUBDIRS="pixman" install
    #make SUBDIRS="pixman test" check -k || echo "expected 1 FAIL (thread-test)"
  fi
  if [[ "${PREPARE}" =~ "boost" ]]
  then
    # boost 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    mkdir -p ~/win/${THREADMODEL}/aarch64/include
    cd ~/win/${THREADMODEL}/aarch64/src
    7zr x ~/win/src/boost_1_76_0.7z
    cd ~/win/${THREADMODEL}/aarch64/include
    rm -f boost
    ln -s ../src/boost_1_76_0/boost/
  fi
  if [[ "${PREPARE}" =~ "glm" ]]
  then
    # glm 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    mkdir -p ~/win/${THREADMODEL}/aarch64/include
    cd ~/win/${THREADMODEL}/aarch64/src
    7zr x ~/win/src/glm-0.9.9.8.7z
    cd ~/win/${THREADMODEL}/aarch64/include
    rm -f glm
    ln -s ../src/glm/glm/
  fi
  if [[ "${PREPARE}" =~ "mingw-std-threads" ]]
  then
    # mingw-std-threads 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/include
    cd ~/win/${THREADMODEL}/aarch64/include
    rm -f mingw-std-threads
    ln -s ~/win/src/mingw-std-threads
  fi
  if [[ "${PREPARE}" =~ "openexr" ]]
  then
    # openexr 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/openexr-2.5.7.tar.gz
    cd openexr-2.5.7/
    #if [[ "x${COMPILER}" = "xgcc" ]]
    #then
    #  patch -p1 < ~/win/src/openexr-2.4.0.patch
    #fi
    sed -i "s/#ifdef _WIN32/#if 0/g" OpenEXR/IlmImf/ImfStdIO.cpp
    sed -i "s/x86_64/aarch64/g" cmake/Toolchain-mingw.cmake
    mkdir -p build
    cd build
    cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain-mingw.cmake -DCMAKE_CXX_FLAGS=-I$HOME/win/${THREADMODEL}/aarch64/include -DZLIB_INCLUDE_DIR=$HOME/win/${THREADMODEL}/aarch64/include -DZLIB_LIBRARY=$HOME/win/${THREADMODEL}/aarch64/lib/libz.a -DCMAKE_INSTALL_PREFIX=$HOME/win/${THREADMODEL}/aarch64 ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "glfw" ]]
  then
    # glfw 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    unzip ~/win/src/glfw-3.3.4.zip
    cd glfw-3.3.4/
    sed -i "s/x86_64/aarch64/g" CMake/x86_64-w64-mingw32.cmake
    mkdir -p build
    cd build
    cmake -DCMAKE_TOOLCHAIN_FILE=../CMake/x86_64-w64-mingw32.cmake -DCMAKE_INSTALL_PREFIX=${HOME}/win/${THREADMODEL}/aarch64/ -DBUILD_SHARED_LIBS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF -DGLFW_BUILD_DOCS=OFF -DGLFW_USE_HYBRID_HPG=ON ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "sdl2" ]]
  then
    # sdl2 64
    mkdir -p ~/win/${THREADMODEL}/aarch64/src
    cd ~/win/${THREADMODEL}/aarch64/src
    tar xaf ~/win/src/SDL2-2.0.18.tar.gz
    cd SDL2-2.0.18/
    ./configure --prefix=${HOME}/win/${THREADMODEL}/aarch64 --host=aarch64-w64-mingw32
    make -j $NCPUS
    make install
  fi
  # clew 64
  # nop
fi
if [[ "${ACTION}" =~ "armv7" ]]
then
  if [[ "${PREPARE}" =~ "gmp" ]]
  then
    # gmp 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/gmp-6.2.1.tar.lz
    cd gmp-6.2.1/
    CC_FOR_BUILD="gcc" CPP_FOR_BUILD="gcc -E" ./configure --build=x86_64-pc-linux-gnu --host=armv7-w64-mingw32 --disable-assembly --prefix=$HOME/win/${THREADMODEL}/armv7
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpfr" ]]
  then
    # mpfr 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/mpfr-4.1.0.tar.xz
    cd mpfr-4.1.0/
    #patch -N -Z -p1 < ~/win/src/allpatches
    ./configure --host=armv7-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/armv7 --with-gmp-build=../gmp-6.2.1 --enable-static --disable-shared
    make -j $NCPUS
    make install
    make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpreal" ]]
  then
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/mpfrc++-3.6.8.tar.gz
    cp -avf mpreal-mpfrc-3.6.8/mpreal.h ~/win/${THREADMODEL}/armv7/include
  fi
  if [[ "${PREPARE}" =~ "zlib" ]]
  then
    # zlib 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/zlib-1.2.11.tar.xz
    cd zlib-1.2.11/
    CC=armv7-w64-mingw32-gcc AR=armv7-w64-mingw32-ar RANLIB=armv7-w64-mingw32-ranlib ./configure --static --prefix=$HOME/win/${THREADMODEL}/armv7
    CC=armv7-w64-mingw32-gcc AR=armv7-w64-mingw32-ar RANLIB=armv7-w64-mingw32-ranlib make -j $NCPUS
    CC=armv7-w64-mingw32-gcc AR=armv7-w64-mingw32-ar RANLIB=armv7-w64-mingw32-ranlib make install
  fi
  if [[ "${PREPARE}" =~ "png" ]]
  then
    # png 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/libpng-1.6.37.tar.xz
    cd libpng-1.6.37/
    ./configure --disable-shared --host=armv7-w64-mingw32 CPPFLAGS="$CPPFLAGS -I$HOME/win/${THREADMODEL}/armv7/include" LDFLAGS="$LDFLAGS -L$HOME/win/${THREADMODEL}/armv7/lib" --prefix=$HOME/win/${THREADMODEL}/armv7
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "jpeg" ]]
  then
    # jpeg 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/libpng-1.6.37.tar.xz
    tar xaf ~/win/src/jpegsrc.v6b2.tar.gz
    cp -avf libpng-1.6.37/config.sub jpeg-6b2/config.sub
    cd jpeg-6b2/
    ./configure --disable-shared --host=armv7-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/armv7
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "tiff" ]]
  then
    # tiff 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/tiff-4.3.0.tar.gz
    cd tiff-4.3.0/
    ./configure --disable-shared --host=armv7-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/armv7
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "gsl" ]]
  then
    # gsl 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/gsl-2.7.tar.gz
    cd gsl-2.7/
    ./configure --disable-shared --host=armv7-w64-mingw32 --prefix=$HOME/win/${THREADMODEL}/armv7
    make -j $NCPUS
    make install
    rm -f gsl-histogram
    ln -s gsl-histogram.exe gsl-histogram # hack for test suite
    make -k check || echo "expected 2 FAILs (multifit_nlinear, multilarge_nlinear)"
  fi
  if [[ "${PREPARE}" =~ "pixman" ]]
  then
    # pixman 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    rm -rf pixman
    cp -avf ~/win/src/pixman pixman
    cd pixman/
    NOCONFIGURE=true ./autogen.sh
    LDFLAGS=-L${HOME}/win/${THREADMODEL}/armv7/lib ./configure --host armv7-w64-mingw32 --disable-shared --disable-openmp --prefix=$HOME/win/${THREADMODEL}/armv7
    make SUBDIRS="pixman" -j $NCPUS
    make SUBDIRS="pixman" install
    #make SUBDIRS="pixman test" check -k || echo "expected 1 FAIL (thread-test)"
  fi
  if [[ "${PREPARE}" =~ "boost" ]]
  then
    # boost 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    7zr x ~/win/src/boost_1_76_0.7z
    cd ~/win/${THREADMODEL}/armv7/include
    rm -f boost
    ln -s ../src/boost_1_76_0/boost/
  fi
  if [[ "${PREPARE}" =~ "glm" ]]
  then
    # glm 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    7zr x ~/win/src/glm-0.9.9.8.7z
    cd ~/win/${THREADMODEL}/armv7/include
    rm -f glm
    ln -s ../src/glm/glm/
  fi
  if [[ "${PREPARE}" =~ "mingw-std-threads" ]]
  then
    # mingw-std-threads 32
    mkdir -p ~/win/${THREADMODEL}/armv7/include
    cd ~/win/${THREADMODEL}/armv7/include
    rm -f mingw-std-threads
    ln -s ~/win/src/mingw-std-threads
  fi
  if [[ "${PREPARE}" =~ "openexr" ]]
  then
    # openexr 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xf ~/win/src/openexr-2.5.7.tar.gz
    cd openexr-2.5.7/
    #if [[ "x${COMPILER}" = "xgcc" ]]
    #then
    #  patch -p1 < ~/win/src/openexr-2.4.0.patch
    #fi
    sed -i "s/#ifdef _WIN32/#if 0/g" OpenEXR/IlmImf/ImfStdIO.cpp
    sed -i "s/x86_64/armv7/g" cmake/Toolchain-mingw.cmake
    mkdir -p build
    cd build
    cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain-mingw.cmake -DCMAKE_CXX_FLAGS=-I$HOME/win/${THREADMODEL}/armv7/include -DZLIB_INCLUDE_DIR=$HOME/win/${THREADMODEL}/armv7/include -DZLIB_LIBRARY=$HOME/win/${THREADMODEL}/armv7/lib/libz.a -DCMAKE_INSTALL_PREFIX=$HOME/win/${THREADMODEL}/armv7 ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "glfw" ]]
  then
    # glfw 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    unzip ~/win/src/glfw-3.3.4.zip
    cd glfw-3.3.4/
    sed -i "s/i686/armv7/g" CMake/i686-w64-mingw32.cmake
    mkdir -p build
    cd build
    cmake -DCMAKE_TOOLCHAIN_FILE=../CMake/i686-w64-mingw32.cmake -DCMAKE_INSTALL_PREFIX=${HOME}/win/${THREADMODEL}/armv7/ -DBUILD_SHARED_LIBS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF -DGLFW_BUILD_DOCS=OFF -DGLFW_USE_HYBRID_HPG=ON ..
    make -j $NCPUS
    make install
  fi
  if [[ "${PREPARE}" =~ "sdl2" ]]
  then
    # sdl2 32
    mkdir -p ~/win/${THREADMODEL}/armv7/src
    cd ~/win/${THREADMODEL}/armv7/src
    tar xaf ~/win/src/SDL2-2.0.18.tar.gz
    cd SDL2-2.0.18/
    ./configure --prefix=${HOME}/win/${THREADMODEL}/armv7 --host=armv7-w64-mingw32
    make -j $NCPUS
    make install
  fi
  # clew 32
  #nop
fi
if [[ "${ACTION}" =~ "emscripten" ]]
then
  if [[ "${PREPARE}" =~ "emsdk" ]]
  then
    mkdir -p ~/opt/emscripten
    cd ~/opt/emscripten
    git clone https://github.com/emscripten-core/emsdk.git
    cd emsdk
    ./emsdk install latest
    ./emsdk activate latest
  fi
  source ~/opt/emscripten/emsdk/emsdk_env.sh
  export EMMAKEN_CFLAGS="-s USE_PTHREADS=1"
  if [[ "${PREPARE}" =~ "gmp" ]]
  then
    mkdir -p ~/opt/emscripten/src
    cd ~/opt/emscripten/src
    tar xaf ~/win/src/gmp-6.2.1.tar.lz
    cd gmp-6.2.1/
    emconfigure ./configure --host none --prefix ${HOME}/opt/emscripten --disable-shared --enable-static --disable-assembly --enable-cxx
    emmake make -j $NCPUS
    emmake make install
    # emmake make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpfr" ]]
  then
    mkdir -p ~/opt/emscripten/src
    cd ~/opt/emscripten/src
    tar xaf ~/win/src/mpfr-4.1.0.tar.xz
    cd mpfr-4.1.0/
    #patch -N -Z -p1 < ~/win/src/allpatches
    emconfigure ./configure --host none --prefix=${HOME}/opt/emscripten --with-gmp=${HOME}/opt/emscripten
    emmake make -j $NCPUS
    emmake make install
    # emmake make check -k || echo
  fi
  if [[ "${PREPARE}" =~ "mpreal" ]]
  then
    cd ~/opt/emscripten/src
    tar xaf ~/win/src/mpfrc++-3.6.8.tar.gz
    cp -avf mpreal-mpfrc-3.6.8/mpreal.h ${HOME}/opt/emscripten/include
  fi
  if [[ "${PREPARE}" =~ "glm" ]]
  then
    cd ~/opt/emscripten/src
    7zr x glm-0.9.9.8.7z
    cp -avf glm/glm ~/opt/emscripten/include/glm
  fi
fi
