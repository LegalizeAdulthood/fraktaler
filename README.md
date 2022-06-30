---
author: Claude Heiland-Allen
toc: true
geometry:
- margin=1in
keywords:
- escape time
- fractal
- deep zoom
- Mandelbrot set
- Burning Ship
header-includes:
- <link rel="shortcut icon" href="fraktaler-3.ico" />
...

# Fraktaler 3

![Screenshot of Fraktaler 3](fraktaler-3.png)

Fast deep escape time fractals.

<https://fraktaler.mathr.co.uk>

Fraktaler 3 is a cross-platform program (Linux, Windows, Android, Web)
for fast deep zooming of hybrid escape-time 2D fractals. It has a
graphical explorer using SDL2, OpenGL and Dear ImGUI, and a command line
version for high resolution images and zoom sequences, with export of
raw data in EXR format compatible with
[Kalles Fraktaler 2 +](https://mathr.co.uk/kf/kf.html) and
[zoomasm](https://mathr.co.uk/zoomasm).

Other fractal deep zoom software that also uses bilinear approximation
(BLA) for acceleration includes:

[fractalshades](https://gbillotey.github.io/Fractalshades/)

:   Fractalshades is a Python package for creating static and
    interactive visualisations of 2d fractals. It targets Windows and
    Unix operating systems and implements efficient algorithms for
    very-deep exploration of the Mandelbrot and the Burning_Ship sets
    (1.e-2000 scale and beyond).

...

:   Get in touch if you know of other software (closed or open source,
    payware or gratis) that is comparable and I'll add it to the list!

## Live

Try Fraktaler 3 live online in your web browser.

<https://fraktaler.mathr.co.uk/live/latest>

Requires support for `SharedArrayBuffer`, among other web APIs.
(This rules out Firefox/Fennec on Android at the moment.)

Performance is significantly slower than native versions, which are
available for download below.

## Download

<https://fraktaler.mathr.co.uk/download>

## Run

### Run GUI

You need support for recent OpenGL.  If you don't have it, the program
window may appear briefly before closing without any error messages
visible.

On Microsoft Windows, if your GPU drivers do not support it you can
install Mesa 3D and the Vulkan Runtime from:

- <https://github.com/pal1000/mesa-dist-win/releases>
- <https://vulkan.lunarg.com/sdk/home#windows>

Use the `mesa-dist-win` per-app deployment script.

```
./fraktaler-3-gui
```

### Run CLI

```
./fraktaler-3-cli parameter.f3.toml
```

or if you have OpenCL, edit the parameter to add

```
opencl.platform = 0
opencl.device = 0
```

(replaced with the indices of your desired platform and device) and run

```
./fraktaler-3-cl parameter.f3.toml
```

OpenCL is only good until about 1e300 zoom (double precision), and
progress reporting is less frequent.  Throughput may be higher depending
on hardware.  OpenCL with CPU is typically faster than the regular CPU
code.  OpenCL with GPU may or may not be faster than OpenCL with CPU.

### Run Web

Configure web server with headers:

```
Cross-Origin-Embedder-Policy: require-corp
Cross-Origin-Resource-Policy: same-origin
Cross-Origin-Opener-Policy: same-origin
```

Make sure `*.wasm` is served with MIME type `application/wasm`

Serve the `live/` sub-folder.  Needs httpS for non-localhost domains.

You must serve the corresponding source code to comply with the license.

## User Manual

### Navigation

The F10 key toggles the graphical user interface windows, so you can
focus on exploring the fractal.  If you do not have a keyboard, you can
manually close all the windows by deselecting their checkboxes in the
main Fraktaler 3 window, which can also be collapsed by clicking on the
small triangle in the top left.

The fractal can be explored with a mouse.  Left mouse down on the
desired new image center and drag to set the new image size; a rectangle
is highlighted during the gesture.  Release the left mouse button to
confirm the new view, or press the right mouse button (while the left
button is still held) to cancel the action.  Alternatively use the
scroll wheel to zoom in and out around the mouse cursor position.  The
middle mouse button centers the view on the click location.

The fractal can be explored with a keyboard.  Numeric keypad keys 1-9
zoom to different quadrants of the view (1 is bottom left, 9 is top
right, 4 is middle left, and so on, as per usual layout).  The 0 key
zooms out.  Numeric keypad keys + and - adjust the maximum iteration
count (doubling and halving respectively), which can also be set in the
Bailout window.

The fractal can be explored with multi-touch.  One touch translates the
view.  Two touches zoom and rotate.  Three touches enables stretching or
skewing the image.  If you have no multi-touch device, but do have a
mouse, you can use multi-touch emulation.  Hold Ctrl+Shift and the left
mouse button to add or move a touch point.  Hold Ctrl+Shift and press
the right mouse button to delete a touch point.  Delete all touch points
to finish the gesture and confirm the action.

### Fraktaler 3 Window

This window has toggles to open/close all the other subwindows.

### Input/Output Window

This has a Home button to zoom out to the original view.  You must set
the checkbox to the left to enable this to avoid accidents.  There are
also buttons to Load and Save, which can be as parameter file text
(suggested extension .f3.toml) or images (EXR format, extension .exr).
Clicking the Load or Save buttons opens a file browser dialog.  Note:
saved images do not yet include parameter metadata, so be sure to save
a parameter file too if you want to return to the location later.

### Formula Window

The default formula is the Mandelbrot set, with one line with |X|, |Y|,
-X, -Y all unchecked and P=2.  This corresponds to the familiar formula
(X+iY)^2 + C.  If you check both |X| and |Y| then you get the Burning
Ship (|X|+i|Y|)^2 + C, if instead you check -Y you get the Mandelbar
(aka Tricorn) (X-iY)^2 + C.  The + button on the right lets you add more
than one formula, which can be edited independently.  These are iterated
in an interleaved fashion, one line after the other in a loop, creating
hybrid escape time fractals.  Note: reference orbit processing and
memory requirements increase with each line (N lines need N times the
amounts total as 1 line);.

### Colour Window

This window lets you adjust the colouring algorithm.  Currently you can
choose between monochrome (black on white) or rainbow colourings.

### Status Window

Shows various progress bars to show how rendering is proceeding.  There
is also a timer.

### Location Window

Shows the coordinates and magnification of the view.

### Reference Window

Shows the coordinates and period (if any) of the reference (which is
usually the image center).

### Bailout Window

Adjust maximum iteration count.  The first two items should usually be
the same, and should be increased if there are solid regions that look
out of place.  The third item can be increased for complex images if
increasing the first two does not fix the issue.  Use the information
window to diagnose the necessary iteration counts.

The escape radius is adjusted at the bottom, decrease it for high power
formulas if unsightly rings appear around the fractal.

### Transform Window

Adjust image transformation, including reflection (useful if your
Burning Ship is upside down), rotation, and stretch.  The exponential
map feature is not so useful in the graphical program, but can be used
in the command line version for rendering a zoom out sequence for later
assembly into a video using zoomasm (<https://mathr.co.uk/zoomasm>).

### Algorithm Window

Contains advanced algorithm tuning options.  Be careful if you adjust
these as sometimes bad images can result.

### Information Window

Displays statistics on various aspects of the calculations, including a
speedup factor for the bilinear approximation algorithm vs doing regular
perturbation iterations for each pixel.

### Quality Window

Control image quality.  Increasing top slider decreases quality (but
increases speed) by subsampling the image.  Increasing the bottom slider
increases quality by computing many versions of the image and averaging
them.  Setting the bottom slider to 0 will compute more subframes
indefinitely, allowing you to stop when the quality gets high enough for
you.

### Newton Zooming Window

Zooms automatically to mini-sets or embedded Julia sets deep in the
fractal.  Set the options (each action includes the ones above), then
select the activate checkbox and left-click in the image where you want
to zoom.  Remember to deselect the activate checkbox if you want to use
the left mouse zooming feature.

### About Window

Displays version information and software licenses.

## Source

You can browse the source code repository at:

<https://code.mathr.co.uk/fraktaler-3>

## Build

### Source Dependencies

```
git clone https://github.com/ocornut/imgui.git
git clone https://github.com/AirGuanZ/imgui-filebrowser.git
git clone https://github.com/ToruNiina/toml11.git
git clone https://code.mathr.co.uk/fraktaler-3.git
```

### Debian Dependencies

Bullseye or newer is required.

```
sudo apt install \
  build-essential \
  clang-11 \
  git \
  libglew-dev \
  libglm-dev \
  libmpfr-dev \
  libmpfrc++-dev \
  libomp-11-dev \
  libopenexr-dev \
  libsdl2-dev \
  p7zip \
  pkg-config \
  xxd
```

### Windows Dependencies

For cross-compilation from Debian.

```
sudo dpkg --add-architecture i386
sudo apt update
sudo apt install \
  build-essential \
  git \
  mingw-w64 \
  p7zip \
  wine32 \
  wine64 \
  wine-binfmt \
  xxd
sudo update-alternatives --set x86_64-w64-mingw32-g++ /usr/bin/x86_64-w64-mingw32-g++-posix
sudo update-alternatives --set x86_64-w64-mingw32-gcc /usr/bin/x86_64-w64-mingw32-gcc-posix
sudo update-alternatives --set x86_64-w64-mingw32-gfortran /usr/bin/x86_64-w64-mingw32-gfortran-posix
sudo update-alternatives --set x86_64-w64-mingw32-gnat /usr/bin/x86_64-w64-mingw32-gnat-posix
sudo update-alternatives --set i686-w64-mingw32-g++ /usr/bin/i686-w64-mingw32-g++-posix
sudo update-alternatives --set i686-w64-mingw32-gcc /usr/bin/i686-w64-mingw32-gcc-posix
sudo update-alternatives --set i686-w64-mingw32-gfortran /usr/bin/i686-w64-mingw32-gfortran-posix
sudo update-alternatives --set i686-w64-mingw32-gnat /usr/bin/i686-w64-mingw32-gnat-posix
```

Use the `prepare.sh` script to download and build dependencies for your
architecture.  For help:

```
./build/prepare.sh -h
```

#### Windows i686

```
make SYSTEM=i686-w64-mingw32
```

#### Windows x86_64

```
make SYSTEM=x86_64-w64-mingw32
```

#### Windows armv7

You need `llvm-mingw` because `gcc-mingw` does not support Windows on
ARM: <https://github.com/mstorsjo/llvm-mingw>

Note: `-lopengl32` is not supported upstream yet, so the GUI won't
compile.

Note: Wine is untested.  Microsoft Windows is untested.

```
make SYSTEM=armv7-w64-mingw32
```

#### Windows aarch64

You need `llvm-mingw` because `gcc-mingw` does not support Windows on
ARM: <https://github.com/mstorsjo/llvm-mingw>

Note: `-lopengl32` is not supported upstream yet, so the GUI won't
compile.

Note: Wine does not yet support `__C_specific_handler`, so it won't run
in Wine.  Microsoft Windows is untested.

```
make SYSTEM=aarch64-w64-mingw32
```

### Emscripten Dependencies

Use the `prepare.sh` script to download and build dependencies for the
`emscripten` architecture.  For help:

```
./build/prepare.sh -h
```

### Build For Android

Use the `android.sh` script to download and build dependencies for
Android.  Needs Android command line tools, SDK, NDK.  Set environment
variables to configure, for example:

```
ANDROID_HOME=${HOME}/opt/android
ANDROID_NDK_HOME=${ANDROID_HOME}/ndk/23.1.7779620
PATH="${ANDROID_HOME}/tools:$PATH"
PATH="${ANDROID_HOME}/platform-tools:$PATH"
PATH="${ANDROID_NDK_HOME}:$PATH"
./build/android.sh prepare
./build/android.sh
```

Default is a debug build (runs slow).  Release build requires signing.

### Build Documentation

Needs `pandoc`.  Built as part of release.

### Build Release

Builds all architectures and documentation ready for release.  Does not
yet include Android.

```
./build/release.sh
```

## Theory

References:

perturbation technique

:   <http://www.science.eclipse.co.uk/sft_maths.pdf>

rebasing and BLA

:   <https://fractalforums.org/f/28/t/4360>

distance estimates

:   <https://mathr.co.uk/helm>

interior detection

:   <https://fractalforums.org/f/28/t/4802>

### The Mandelbrot Set

High precision reference orbit:

$$Z_{m+1} = Z_m^2 + C$$

$m$ starts at $0$ with $Z_0 = 0$.

### Perturbation

Low precision deltas relative to high precision orbit.
Pixel orbit $Z_m + z_n$, $C + c$.

$$z_{n+1} = 2 Z_m z_n + z_n^2 + c$$

$m$ and $n$ start at $0$ with $z_0 = 0$.

### Rebasing

Rebasing to avoid glitches: when
$$|Z_m + z_n| < |z_n|$$
replace $z_n$ with $Z_m + z_n$ and
reset the reference iteration count $m$ to $0$.

### Bivariate Linear Approximation

When $Z$ is large and $z$ is small, the iterations can be approximated
by bivariate linear function;

$$z_{n+l} = A_{n,l} z_n + B_{n,l} c$$

This is valid when the non-linear part of the full perturbation
iterations is so small that omitting it would cause fewer problems than
the rounding error of the low precision data type.

### Single Step BLA

Approximation of a single step by bilinear form is valid when
$$\begin{aligned}
|z_n^2| &<< |2 Z_n z_n + c| \\
&\Uparrow \quad \text{ definition of $A_{n,1}, B_{n,1}$ for single step } \\
|z_n^2| &<< |A_{n,1} z_n + B_{n,1} c| \\
&\Uparrow \quad \text{ definition of $\epsilon$ (for example, $\epsilon = 2^{-24}$) } \\
|z_n^2| &< \epsilon |A_{n,1} z_n + B_{n,1} c| \\
&\Uparrow \quad \text{ triangle inequality } \\
|z_n^2| &< \epsilon |A_{n,1} z_n| - \epsilon |B_{n,1} c| \\
&\Uparrow \quad \text{ algebra } \\
|z_n|^2 - \epsilon |A_{n,1}| |z_n| + \epsilon |B_{n,1}| |c_n| &< 0 \\
&\Uparrow \quad \text{ quadratic formula } \\
|z_n| &< \frac{\epsilon |A_{n,1}| + \sqrt{ (\epsilon |A_{n,1}|)^2 - 4 \epsilon |B_{n,1}| |c| }}{2} \\
&\Uparrow \quad \text{ linear Taylor polynomial (**approximation**) } \\
|z_n| &< \epsilon |A_{n,1}| - \frac{|B_{n,1}|}{|A_{n,1}|} |c| =: R_{n,1}
\end{aligned}$$

For single step of Mandelbrot set:
$$\begin{aligned}
A_{m,1} &= \frac{\partial Z_{m+1}}{\partial Z_m} = 2 Z_m \\
B_{m,1} &= \frac{\partial Z_{m+1}}{\partial C} = 1 \\
R_{m,1} &= \max\left\{ 0, \epsilon 2 |Z_m| - \frac{|c|}{2 |Z_m|} \right\}
\end{aligned}$$

Note: this is different to the formulas suggested by Zhuoran on Fractal
Forums, but I couldn't get them to work, and this version does seem to
work fine.

### Merging BLA Steps

If $T_x$ skips $l_x$ iterations from iteration $m_x$ when $|z| < R_x$
and $T_y$ skips $l_y$ iterations from iteration $m_x + l_x$ when $|z| < R_y$
then $T_z = T_y \circ T_x$ skips $l_x + l_y$ iterations from iteration $m_x$ when $|z| < R_z$:
$$\begin{aligned}
z_{m_x + l_x + l_y} &= A_y (A_x z_{m_x} + B_x c) + B_y c = A_z z_{m_x} + B_z c \\
A_{m_x, l_x + l_y} = A_z &= A_y A_x \\
B_{m_x, l_x + l_y} = B_z &= A_y B_x + B_y \\
R_{m_x, l_x + l_y} = R_z &= \max\left\{ 0, \min\left\{ R_x, \frac{R_y - |B_x| |c|}{|A_x|} \right\} \right\}
\end{aligned}$$

### BLA Table Construction

Suppose the reference has $M$ iterations.  Create $M$ BLAs each skipping
$1$ iteration (this can be done in parallel).  Then merge neighbours
without overlap to create $\left\lceil \frac{M}{2} \right\rceil$ each
skipping $2$ iterations (except for perhaps the last which skips less).
Repeat until there is only $1$ BLA skipping $M-1$ iterations: it's best
to start the merge from iteration $1$ because reference iteration $0$
always corresponds to a non-linear perturbation step as $Z = 0$.

The resulting table has $O(M)$ elements.

### BLA Table Lookup

Find the BLA starting from iteration $m$ that has the largest skip $l$
satisfying $|z| < R$.  If there is none, do a perturbation iteration.
Check for rebasing opportunities after each BLA application or
perturbation step.

### Non-Conformal BLA

The Mandelbrot set is conformal (angles are preserved).  This means
complex numbers can be used for derivatives.  Some other formulas are
not conformal, for example the Tricorn aka Mandelbar, defined by:
$$ X + i Y \to (X - i Y)^2 + C $$

For non-conformal formulas, replace complex numbers by $2 \times 2$ real
matrices for $A, B$.  Dual numbers with two dual parts can be used to
calculate the derivatives.

Be careful finding norms.  Define $\sup|M|$ and $\inf|M|$ as the largest
and smallest singular values of $M$.  Then single step BLA radius
becomes
$$R = \epsilon \inf|A| - \frac{\sup|B|}{\inf|A|} |c|$$
and merging BLA steps radius becomes
$$R_z = \max\left\{ 0, \min\left\{ R_x, \frac{R_y - \sup|B_x| |c|}{\sup|A_x|} \right\} \right\}$$

### ABS Variation BLA

The only problem with the Mandelbrot set is the non-linearity, but some
other formulas have other problems, for example the Burning Ship,
defined by:
$$ X + i Y \to (|X| + i |Y|)^2 + C $$
The absolute value folds the plane when $X$ or $Y$ are near $0$, so the
single step BLA radius becomes the minimum of the non-linearity radius
and the folding radii:
$$ R = \max\left\{ 0, \min\left\{ \epsilon \inf|A| - \frac{\sup|B|}{\inf|A|} |c|, |X|, |Y| \right\} \right\} $$
Currently Fraktaler 3 uses a fudge factor for paranoia, dividing $|X|$
and $|Y|$ by $2$.  The merged BLA step radius is unchanged.

### Hybrid BLA

For a hybrid loop with multiple phases, you need multiple references,
one starting at each phase in the loop.  Rebasing switches to the
reference for the current phase.  You need one BLA table per reference.

### Multiple Critical Points

Some formulas (but none among those implemented in Fraktaler 3) have
multiple critical points.  In this case some modifications need to be
made: you need a reference per critical point, and rebasing needs to
switch to the nearest orbit among all critical points.  There needs to
be a separate BLA table for each reference.  This also applies to
hybrids, you need one reference and BLA table per critical point per
phase.

### Distance Estimation

Keep track of derivatives of $Z+z$ wrt. pixel coordinates $k$.  As $Z$
is constant for the whole image, you just need $\frac{dz}{dk}$.  An easy
way to do this is with dual numbers for automatic numeric
differentiation.  Set up the pixel coordinates as dual numbers with dual
part $1+0i$, then transform them to the complex C+c plane of the fractal
iterations.  At the end you plug the complex derivative into the
(directional) distance estimate formula, it is already prescaled by the
pixel spacing (this also helps to avoid overflow during iteration).

For non-complex-analytic formulas (like Mandelbar and Burning Ship), you
can use dual numbers with two dual parts, for each of the real and
imaginary components.  At the end they can be combined into a Jacobian
matrix and used in the (directional) distance estimate formula for
general iterations.

### Interior Detection

Keep track of derivatives of $Z+z$ wrt. $Z_1+z_1$ (where $Z_0+z_0$ is at
the critical point, usually $0$).  When the absolute value of the
derivative drops below a threshold such as $0.001$, classify it as
interior and stop iterating.  For non-complex-analytic formulas, dual
numbers with four dual parts can be used (two for distance estimation
and two for interior detection), along with matrix operator norm.

Using $\frac{dz}{dz1}$ works because:

$$\begin{aligned}
 &\frac{d(Z+x)}{d(Z_1+z_1)} \\
=& \frac{dZ}{d(Z_1+z_1)} + \frac{dz}{d(Z_1+z_1)} \\
=& \frac{1}{\frac{dZ_1}{dZ} + \frac{dz_1}{dZ}} + \frac{1}{\frac{dZ_1}{dz} + \frac{dz_1}{dz}} \\
=& \frac{1}{\frac{dZ_1}{dZ} + 0} + \frac{1}{0 + \frac{dz_1}{dz}} \\
=& \frac{dZ}{dZ_1} + \frac{dz}{dz_1} \\
=& 0 + \frac{dz}{dz_1} \\
=& \frac{dz}{dz_1}
\end{aligned}$$

where the last two lines hold when $C$ is periodic with $Z = 0$ in the
orbit which happens precisely when the formula has a critical point at
$0$ and $C$ is the nucleus of a hyperbolic component.

## TODO

- fix (or remove from GUI) image IO
  - GUI should save EXR if `*.exr`, currently saves TOML
  - should save metadata to EXR
  - should load metadata from images
  - CLI should have an option to save TOML from argument (which could
    be an image)
  - CLI should have an option to merge TOML parameters (for example to
    keep palette or OpenCL settings separate from image location)
  - implement EXR channel output filters (to save disk space and time)
- implement low + high bailout
  - ensure BLA doesn't escape past low bailout
  - don't use BLA between low bailout and high bailout
  - store cooked values at low bailout
  - store cooked values at high bailout
  - option to rename channels to avoid clashes
  - channel filters to save memory and calculation time (no-DE mode?)
- optimize MPFR memory allocation
  - reference orbit
  - period detection
  - root finding
  - size calculation
- optimize conformal formulas
  - use complex numbers instead of matrices
  - Mandelbrot set / multibrot only
- extend OpenCL to other number types
  - float
  - floatexp
  - softfloat
  - ensure zoom out sequence works as expected across type boundaries,
    even with reuse reference enabled
  - ensure long double, float128 are never used for references or bla
- support OpenCL in regular CLI renderer
- support OpenCL in GUI
- high resolution rendering dialog
  - dimensions in inches and dots per inch
  - automatically translated to/from pixels
  - exports to toml for command line renderer
  - option to select OpenCL platform and device
  - option to enable reuse reference and zoom out sequence
- extend colouring algorithms
  - parameterize
  - allow custom OpenCL source for colouring snippet (no parameters)
  - allow custom GLSL source with dynamically generated UI for uniforms
  - use OpenCL/OpenGL interop to do colouring with custom GLSL with UI
- compat with other software
  - KFR location import, including metadata from image files
  - KFP palette import (with default GLSL implementation copied from KF)
  - KF custom GLSL import mode (see zoomasm)
  - custom GLSL export for zoomasm
- Windows
  - support OpenCL via CLEW
- Android
  - mechanisms to access SD card or clipboard or Share With or something
  - put copy/paste buttons in IO window so touchscreen can be used?
- Web
  - fix copy/paste from host OS into ImGUI dialog boxes
  - export/import parameters to/from host clipboard or via up/download
  - export/import parameters to/from URL hash (base64)
  - export image as download (browser canvas right click is captured)

## Legal

Fraktaler 3 -- Fast deep escape time fractals

Copyright (C) 2021,2022  Claude Heiland-Allen

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
