---
author: Claude Heiland-Allen
toc: true
geometry:
- margin=1in
...

# Fraktaler 3

Fast deep escape time fractals.

<https://mathr.co.uk/f3>

## Try It Online

<https://mathr.co.uk/f3/live/latest>

Requires support for `SharedArrayBuffer`, among other web APIs.

## Run

### Run GUI

```
./fraktaler-3-gui
```

### Run CLI

```
./fraktaler-3-cli
```

### Run Web

Configure web server with headers:

```
Cross-Origin-Embedder-Policy: require-corp
Cross-Origin-Resource-Policy: same-origin
Cross-Origin-Opener-Policy: same-origin
```

Make sure `*.wasm` is served with MIME type `application/wasm`

Serve the `live/` sub-folder.  Needs httpS for non-localhost domains.

You must serve the corresponding source code, for legal reasons.

## Build

### Build GUI

```
sudo apt install \
  build-essential \
  git \
  libglew-dev \
  libglm-dev \
  libmpfr-dev \
  libmpfrc++-dev \
  libopenexr-dev \
  libsdl2-dev \
  p7zip \
  pkg-config
git clone https://github.com/ocornut/imgui.git
git clone https://code.mathr.co.uk/fraktaler-3.git
cd fraktaler-3
make gui -j $(nproc)
```

### Build CLI

```
sudo apt install \
  build-essential \
  git \
  libglm-dev \
  libmpfr-dev \
  libmpfrc++-dev \
  libopenexr-dev \
  p7zip \
  pkg-config
git clone https://code.mathr.co.uk/fraktaler-3.git
cd fraktaler-3
make cli -j $(nproc)
```

### Build Web

```
sudo apt install \
  build-essential \
  git \
  libglew-dev \
  libglm-dev \
  libmpfr-dev \
  libmpfrc++-dev \
  libopenexr-dev \
  libsdl2-dev \
  p7zip \
  pkg-config
# download
wget "https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz"
wget "https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.xz"
wget "https://github.com/advanpix/mpreal/archive/refs/tags/mpfrc++-3.6.8.tar.gz"
wget "https://github.com/g-truc/glm/releases/download/0.9.9.8/glm-0.9.9.8.7z"
git clone https://github.com/emscripten-core/emsdk.git
git clone https://github.com/ocornut/imgui.git
git clone https://code.mathr.co.uk/fraktaler-3.git
# emscripten
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
cd ..
mkdir -p ${HOME}/opt/emscripten
export EMMAKEN_CFLAGS="-s USE_PTHREADS=1"
# gmp
tar xaf gmp-6.2.1.tar.lz
cd gmp-6.2.1
emconfigure ./configure \
  --host none \
  --prefix ${HOME}/opt/emscripten \
  --disable-shared \
  --enable-static \
  --disable-assembly \
  --enable-cxx
emmake make -j $(nproc)
emmake make install
cd ..
# mpfr
tar xaf mpfr-4.1.0.tar.xz
cd mpfr-4.1.0
emconfigure ./configure \
  --host none \
  --prefix=${HOME}/opt/emscripten \
  --with-gmp=${HOME}/opt/emscripten
emmake make -j $(nproc)
emmake make install
cd ..
# mpfrc++
tar xaf mpfrc++-3.6.8.tar.gz
cp mpreal-mpfrc-3.6.8/mpreal.h ${HOME}/opt/emscripten/include
# glm
7zr x glm-0.9.9.8.7z
cp -a glm/glm ~/opt/emscripten/include/glm
cd fraktaler-3
make web -j $(nproc)
```

### Build Documentation

Needs `pandoc`.  Built as part of release.

### Build Release

```
make release -j $(nproc)
```

## Legal

Fraktaler 3 -- Fast deep escape time fractals

Copyright (C) 2021  Claude Heiland-Allen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

---
<https://mathr.co.uk>
