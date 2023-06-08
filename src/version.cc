// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021-2023 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <mpreal.h>
#include <toml.hpp>
#ifdef HAVE_EXR
#if HAVE_EXR == 0
#undef HAVE_EXR
#else
#include <zlib.h>
#include <OpenEXR/OpenEXRConfig.h>
#endif
#endif
#include <SDL.h>
#include <glm/glm.hpp>
#ifdef HAVE_GUI
#include "gles2.h"
#include <imgui.h>
#ifdef HAVE_FS
#include <imfilebrowser.h>
#endif
#include <implot.h>
#endif

#include "source.h"
#include "version.h"

std::string version(const char *gl_version)
{
  (void) gl_version;
  std::ostringstream out;
  out << "fraktaler-3 version " << fraktaler_3_version_string << "\n";
  out << "\n";
#ifdef __clang__
  out << "clang version " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << "\n";
#else
#ifdef __GNUC__
  out << "g++ version " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
#else
  out << "compiler version unknown\n";
#endif
#endif
  out << "gmp version " << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << "\n";
  out << "mpfr version " << MPFR_VERSION_MAJOR << "." << MPFR_VERSION_MINOR << "." << MPFR_VERSION_PATCHLEVEL << "\n";
  out << "mpfrc++ version " << MPREAL_VERSION_MAJOR << "." << MPREAL_VERSION_MINOR << "." << MPREAL_VERSION_PATCHLEVEL << "\n";
  out << "toml11 version " << TOML11_VERSION_MAJOR << "." << TOML11_VERSION_MINOR << "." << TOML11_VERSION_PATCH << " (" << TOML11_GIT_VERSION_STRING << ")\n";
#ifdef HAVE_EXR
  out << "zlib version " << zlib_version << "\n";
  out << "openexr version " << OPENEXR_VERSION_STRING << "\n";
#endif
  out << "sdl2 version " << SDL_MAJOR_VERSION << "." << SDL_MINOR_VERSION << "." << SDL_PATCHLEVEL << "\n";
#ifdef HAVE_CLEW
  out << "clew version " << CLEW_GIT_VERSION_STRING << "\n";
#endif
  out << "glm version " << GLM_VERSION_MAJOR << "." << GLM_VERSION_MINOR << "." << GLM_VERSION_PATCH << "." << GLM_VERSION_REVISION << "\n";
#ifdef HAVE_GUI
  out << "glad version " << GLAD_GENERATOR_VERSION << "\n";
  out << "imgui version " << IMGUI_VERSION << " (" << IMGUI_GIT_VERSION_STRING << ")\n";
#ifdef HAVE_FS
  out << "imgui-filebrowser version " << IMGUI_FILEBROWSER_GIT_VERSION_STRING << "\n";
#endif
  out << "implot version " << IMPLOT_VERSION << " (" << IMPLOT_GIT_VERSION_STRING << ")\n";
  out << "\n";
  out << "opengl version " << (gl_version ? gl_version : "unknown") << "\n";
#endif
  return out.str();
}

std::string license()
{
  return
"fraktaler-3 license\n"
"\n"
"\n"
"Fraktaler 3 -- Fast deep escape time fractals\n"
"Copyright (C) 2021-2023 Claude Heiland-Allen\n"
"\n"
"This program is free software: you can redistribute it and/or modify\n"
"it under the terms of the GNU Affero General Public License as\n"
"published by the Free Software Foundation, version 3.\n"
"\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"GNU Affero General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU Affero General Public License\n"
"along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"\n"
"\n"
"\n"
"\n"

"gmp license\n"
"\n"
"\n"
"The GNU MP Library\n"
"Copyright 1991, 1996, 1999-2016, 2018-2020 Free Software Foundation, Inc.\n"
"\n"
"The GNU MP Library is free software; you can redistribute it and/or modify\n"
"it under the terms of either:\n"
"\n"
"  * the GNU Lesser General Public License as published by the Free\n"
"    Software Foundation; either version 3 of the License, or (at your\n"
"    option) any later version.\n"
"\n"
"or\n"
"\n"
"  * the GNU General Public License as published by the Free Software\n"
"    Foundation; either version 2 of the License, or (at your option) any\n"
"    later version.\n"
"\n"
"or both in parallel, as here.\n"
"\n"
"The GNU MP Library is distributed in the hope that it will be useful, but\n"
"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY\n"
"or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License\n"
"for more details.\n"
"\n"
"You should have received copies of the GNU General Public License and the\n"
"GNU Lesser General Public License along with the GNU MP Library.  If not,\n"
"see https://www.gnu.org/licenses/.\n"
"\n"
"\n"
"\n"
"\n"

"mpfr license\n"
"\n"
"\n"
"The GNU MPFR Library\n"
"Copyright 2000-2023 Free Software Foundation, Inc.\n"
"Contributed by the AriC and Caramba projects, INRIA.\n"
"\n"
"The GNU MPFR Library is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU Lesser General Public License as published by\n"
"the Free Software Foundation; either version 3 of the License, or (at your\n"
"option) any later version.\n"
"\n"
"The GNU MPFR Library is distributed in the hope that it will be useful, but\n"
"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY\n"
"or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public\n"
"License for more details.\n"
"\n"
"You should have received a copy of the GNU Lesser General Public License\n"
"along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see\n"
"https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,\n"
"51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.\n"
"\n"
"\n"
"\n"
"\n"

"mpfrc++ license\n"
"\n"
"\n"
"MPFR C++: Multi-precision floating point number class for C++.\n"
"Based on MPFR library:    http://mpfr.org\n"
"\n"
"Project homepage:    http://www.holoborodko.com/pavel/mpfr\n"
"Contact e-mail:      pavel@holoborodko.com\n"
"\n"
"Copyright (c) 2008-2020 Pavel Holoborodko\n"
"\n"
"Contributors:\n"
"Dmitriy Gubanov, Konstantin Holoborodko, Brian Gladman,\n"
"Helmut Jarausch, Fokko Beekhof, Ulrich Mutze, Heinz van Saanen,\n"
"Pere Constans, Peter van Hoof, Gael Guennebaud, Tsai Chia Cheng,\n"
"Alexei Zubanov, Jauhien Piatlicki, Victor Berger, John Westwood,\n"
"Petr Aleksandrov, Orion Poplawski, Charles Karney, Arash Partow,\n"
"Rodney James, Jorge Leitao, Jerome Benoit, Michal Maly.\n"
"\n"
"Licensing:\n"
"(A) MPFR C++ is under GNU General Public License (\"GPL\").\n"
"\n"
"(B) Non-free licenses may also be purchased from the author, for users who\n"
"    do not want their programs protected by the GPL.\n"
"\n"
"    The non-free licenses are for users that wish to use MPFR C++ in\n"
"    their products but are unwilling to release their software\n"
"    under the GPL (which would require them to release source code\n"
"    and allow free redistribution).\n"
"\n"
"    Such users can purchase an unlimited-use license from the author.\n"
"    Contact us for more details.\n"
"\n"
"GNU General Public License (\"GPL\") copyright permissions statement:\n"
"**************************************************************************\n"
"This program is free software: you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License as published by\n"
"the Free Software Foundation, either version 3 of the License, or\n"
"(at your option) any later version.\n"
"\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"GNU General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU General Public License\n"
"along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
"\n"
"\n"
"\n"
"\n"

"toml11 license\n"
"\n"
"\n"
"The MIT License (MIT)\n"
"\n"
"Copyright (c) 2017-2022 Toru Niina\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in\n"
"all copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
"THE SOFTWARE.\n"
"\n"
"\n"
"\n"
"\n"

#ifdef HAVE_EXR

"zlib license\n"
"\n"
"\n"
"Copyright (C) 1995-2022 Jean-loup Gailly and Mark Adler\n"
"\n"
"This software is provided 'as-is', without any express or implied\n"
"warranty.  In no event will the authors be held liable for any damages\n"
"arising from the use of this software.\n"
"\n"
"Permission is granted to anyone to use this software for any purpose,\n"
"including commercial applications, and to alter it and redistribute it\n"
"freely, subject to the following restrictions:\n"
"\n"
"1. The origin of this software must not be misrepresented; you must not\n"
"   claim that you wrote the original software. If you use this software\n"
"   in a product, an acknowledgment in the product documentation would be\n"
"   appreciated but is not required.\n"
"2. Altered source versions must be plainly marked as such, and must not be\n"
"   misrepresented as being the original software.\n"
"3. This notice may not be removed or altered from any source distribution.\n"
"\n"
"Jean-loup Gailly        Mark Adler\n"
"jloup@gzip.org          madler@alumni.caltech.edu\n"
"\n"
"\n"
"\n"
"\n"

"openexr license\n"
"\n"
"\n"
"Copyright Contributors to the OpenEXR Project. All rights reserved.\n"
"\n"
"Redistribution and use in source and binary forms, with or without\n"
"modification, are permitted provided that the following conditions are\n"
"met:\n"
"\n"
"1. Redistributions of source code must retain the above copyright notice,\n"
"   this list of conditions and the following disclaimer.\n"
"\n"
"2. Redistributions in binary form must reproduce the above copyright notice,\n"
"   this list of conditions and the following disclaimer in the documentation\n"
"   and/or other materials provided with the distribution.\n"
"\n"
"3. Neither the name of the copyright holder nor the names of its\n"
"   contributors may be used to endorse or promote products derived from\n"
"   this software without specific prior written permission.\n"
"\n"
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n"
"AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n"
"IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE\n"
"ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE\n"
"LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR\n"
"CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF\n"
"SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS\n"
"INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN\n"
"CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)\n"
"ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n"
"POSSIBILITY OF SUCH DAMAGE.\n"
"\n"
"\n"
"\n"
"\n"

#endif

"sdl2 license\n"
"\n"
"\n"
"Copyright (C) 1997-2023 Sam Lantinga <slouken@libsdl.org>\n"
"\n"
"This software is provided 'as-is', without any express or implied\n"
"warranty.  In no event will the authors be held liable for any damages\n"
"arising from the use of this software.\n"
"\n"
"Permission is granted to anyone to use this software for any purpose,\n"
"including commercial applications, and to alter it and redistribute it\n"
"freely, subject to the following restrictions:\n"
"\n"
"1. The origin of this software must not be misrepresented; you must not\n"
"   claim that you wrote the original software. If you use this software\n"
"   in a product, an acknowledgment in the product documentation would be\n"
"   appreciated but is not required. \n"
"2. Altered source versions must be plainly marked as such, and must not be\n"
"   misrepresented as being the original software.\n"
"3. This notice may not be removed or altered from any source distribution.\n"
"\n"
"\n"
"\n"
"\n"

#ifdef HAVE_CLEW

"clew license\n"
"\n"
"\n"
"(c) 2009-2011   Organic Vectory B.V.\n"
"Written by George van Venrooij  <george (at) organicvectory (dot) com>\n"
"\n"
"Boost Software License - Version 1.0 - August 17th, 2003\n"
"\n"
"Permission is hereby granted, free of charge, to any person or organization\n"
"obtaining a copy of the software and accompanying documentation covered by\n"
"this license (the \"Software\") to use, reproduce, display, distribute,\n"
"execute, and transmit the Software, and to prepare derivative works of the\n"
"Software, and to permit third-parties to whom the Software is furnished to\n"
"do so, all subject to the following:\n"
"\n"
"The copyright notices in the Software and this entire statement, including\n"
"the above license grant, this restriction and the following disclaimer,\n"
"must be included in all copies of the Software, in whole or in part, and\n"
"all derivative works of the Software, unless such copies or derivative\n"
"works are solely in the form of machine-executable object code generated by\n"
"a source language processor.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT\n"
"SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE\n"
"FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,\n"
"ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER\n"
"DEALINGS IN THE SOFTWARE.\n"
"\n"
"\n"
"\n"
"\n"

#endif

"glm licence\n"
"\n"
"\n"
"OpenGL Mathematics (GLM)\n"
"GLM is licensed under The Happy Bunny License and MIT License\n"
"\n"
"\n"
"The Happy Bunny License (Modified MIT License)\n"
"\n"
"Copyright (c) 2005 - 2014 G-Truc Creation\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in\n"
"all copies or substantial portions of the Software.\n"
"\n"
"Restrictions:\n"
" By making use of the Software for military purposes, you choose to make a\n"
" Bunny unhappy.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
"THE SOFTWARE.\n"
"\n"
"\n"
"The MIT License\n"
"\n"
"Copyright (c) 2005 - 2014 G-Truc Creation\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in\n"
"all copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
"THE SOFTWARE.\n"
"\n"
"\n"
"\n"
"\n"

#ifdef HAVE_GUI

"glad license\n"
"\n"
"\n"
"The glad source code:\n"
"\n"
"    The MIT License (MIT)\n"
"\n"
"    Copyright (c) 2013-2020 David Herberth\n"
"\n"
"    Permission is hereby granted, free of charge, to any person obtaining a copy of\n"
"    this software and associated documentation files (the \"Software\"), to deal in\n"
"    the Software without restriction, including without limitation the rights to\n"
"    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of\n"
"    the Software, and to permit persons to whom the Software is furnished to do so,\n"
"    subject to the following conditions:\n"
"\n"
"    The above copyright notice and this permission notice shall be included in all\n"
"    copies or substantial portions of the Software.\n"
"\n"
"    THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS\n"
"    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR\n"
"    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER\n"
"    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN\n"
"    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n"
"\n"
"\n"
"The Khronos Specifications:\n"
"\n"
"    Copyright (c) 2013-2020 The Khronos Group Inc.\n"
"\n"
"    Licensed under the Apache License, Version 2.0 (the \"License\");\n"
"    you may not use this file except in compliance with the License.\n"
"    You may obtain a copy of the License at\n"
"\n"
"        http://www.apache.org/licenses/LICENSE-2.0\n"
"\n"
"    Unless required by applicable law or agreed to in writing, software\n"
"    distributed under the License is distributed on an \"AS IS\" BASIS,\n"
"    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
"    See the License for the specific language governing permissions and\n"
"    limitations under the License.\n"
"\n"
"\n"
"The EGL Specification and various headers:\n"
"\n"
"    Copyright (c) 2007-2016 The Khronos Group Inc.\n"
"\n"
"    Permission is hereby granted, free of charge, to any person obtaining a\n"
"    copy of this software and/or associated documentation files (the\n"
"    \"Materials\"), to deal in the Materials without restriction, including\n"
"    without limitation the rights to use, copy, modify, merge, publish,\n"
"    distribute, sublicense, and/or sell copies of the Materials, and to\n"
"    permit persons to whom the Materials are furnished to do so, subject to\n"
"    the following conditions:\n"
"\n"
"    The above copyright notice and this permission notice shall be included\n"
"    in all copies or substantial portions of the Materials.\n"
"\n"
"    THE MATERIALS ARE PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,\n"
"    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF\n"
"    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.\n"
"    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY\n"
"    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,\n"
"    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE\n"
"    MATERIALS OR THE USE OR OTHER DEALINGS IN THE MATERIALS.\n"
"\n"
"\n"
"\n"
"\n"

"imgui license\n"
"\n"
"\n"
"The MIT License (MIT)\n"
"\n"
"Copyright (c) 2014-2022 Omar Cornut\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in all\n"
"copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
"SOFTWARE.\n"
"\n"
"\n"
"\n"
"\n"

#ifdef HAVE_FS

"imgui-filebrowser license\n"
"\n"
"\n"
"MIT License\n"
"\n"
"Copyright (c) 2019-2022 Zhuang Guan\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in all\n"
"copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
"SOFTWARE.\n"
"\n"
"\n"
"\n"
"\n"

#endif

"implot license\n"
"\n"
"\n"
"MIT License\n"
"\n"
"Copyright (c) 2020 Evan Pezent\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in all\n"
"copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
"SOFTWARE.\n"
"\n"
"\n"
"\n"
"\n"

#endif
  ;
}
