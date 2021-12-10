// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include "map.h"

#ifndef __EMSCRIPTEN__
#include <ImfNamespace.h>
#include <ImfOutputFile.h>
#include <ImfHeader.h>
#include <ImfChannelList.h>
#include <ImfIntAttribute.h>
#include <ImfStringAttribute.h>
#include <ImfArray.h>
#include <ImfFrameBuffer.h>

namespace IMF = OPENEXR_IMF_NAMESPACE;
using namespace IMF;
using namespace IMATH_NAMESPACE;

const char kf2plus[] = "KallesFraktaler2+";
const char fraktaler3[] = "Fraktaler3";
#endif

void map::saveEXR(const std::filesystem::path &filename, const channel_mask_t channels, const int threads, const std::string &metadata, const std::string &kf2plus_metadata) const
{
#ifndef __EMSCRIPTEN__
  setGlobalThreadCount(threads);
  // prepare preview image
  Header header(width, height);
  // insert metadata
  if (metadata != "") header.insert(fraktaler3, StringAttribute(metadata));
  if (kf2plus_metadata != "") header.insert(kf2plus, StringAttribute(kf2plus_metadata));
  if (maxiters + Nbias < INT_MAX - 1)
  {
    header.insert("Iterations", IntAttribute(maxiters));
  }
  else
  {
    std::ostringstream s;
    s << maxiters;
    header.insert("Iterations", StringAttribute(s.str()));
  }
  if (Nbias < INT_MAX - 1)
  {
    header.insert("IterationsBias", IntAttribute(Nbias));
  }
  else
  {
    std::ostringstream s;
    s << Nbias;
    header.insert("IterationsBias", StringAttribute(s.str()));
  }
  // write image
  if (RGB && (channels & (1 << Channel_R))) header.channels().insert("R", Channel(IMF::HALF));
  if (RGB && (channels & (1 << Channel_G))) header.channels().insert("G", Channel(IMF::HALF));
  if (RGB && (channels & (1 << Channel_B))) header.channels().insert("B", Channel(IMF::HALF));
  if (N0 && (channels & (1 << Channel_N0)) &&
      N1 && (channels & (1 << Channel_N1)))
  {
    header.channels().insert("N0",  Channel(IMF::UINT));
    header.channels().insert("N1",  Channel(IMF::UINT));
  }
  else if (N0 && (channels & (1 << Channel_N0)))
  {
    header.channels().insert("N",  Channel(IMF::UINT));
  }
  if (NF  && (channels & (1 << Channel_NF ))) header.channels().insert("NF",  Channel(IMF::FLOAT));
  if (T   && (channels & (1 << Channel_T  ))) header.channels().insert("T",   Channel(IMF::FLOAT));
  if (DEX && (channels & (1 << Channel_DEX))) header.channels().insert("DEX", Channel(IMF::FLOAT));
  if (DEY && (channels & (1 << Channel_DEY))) header.channels().insert("DEY", Channel(IMF::FLOAT));
  OutputFile of(filename.c_str(), header);
  FrameBuffer fb;
  if (RGB && (channels & (1 << Channel_R))) fb.insert("R", Slice(IMF::HALF, (char *)(RGB + 0), sizeof(*RGB) * 3, sizeof(*RGB) * 3 * width));
  if (RGB && (channels & (1 << Channel_G))) fb.insert("G", Slice(IMF::HALF, (char *)(RGB + 1), sizeof(*RGB) * 3, sizeof(*RGB) * 3 * width));
  if (RGB && (channels & (1 << Channel_B))) fb.insert("B", Slice(IMF::HALF, (char *)(RGB + 2), sizeof(*RGB) * 3, sizeof(*RGB) * 3 * width));
  if (N0 && (channels & (1 << Channel_N0)) && N1 && (channels & (1 << Channel_N1)))
  {
    fb.insert("N0", Slice(IMF::UINT, (char *)(N0), sizeof(*N0), sizeof(*N0) * width));
    fb.insert("N1", Slice(IMF::UINT, (char *)(N1), sizeof(*N1), sizeof(*N1) * width));
  }
  else if (N0 && (channels & (1 << Channel_N0)))
  {
    fb.insert("N", Slice(IMF::UINT, (char *)(N0), sizeof(*N0), sizeof(*N0) * width));
  }

  if (NF  && (channels & (1 << Channel_NF ))) fb.insert("NF",  Slice(IMF::FLOAT, (char *)(NF ), sizeof(*NF ), sizeof(*NF ) * width));
  if (T   && (channels & (1 << Channel_T  ))) fb.insert("T",   Slice(IMF::FLOAT, (char *)(T  ), sizeof(*T  ), sizeof(*T  ) * width));
  if (DEX && (channels & (1 << Channel_DEY))) fb.insert("DEX", Slice(IMF::FLOAT, (char *)(DEX), sizeof(*DEX), sizeof(*DEX) * width));
  if (DEY && (channels & (1 << Channel_DEY))) fb.insert("DEY", Slice(IMF::FLOAT, (char *)(DEY), sizeof(*DEY), sizeof(*DEY) * width));
  of.setFrameBuffer(fb);
  of.writePixels(height);
#endif
}
