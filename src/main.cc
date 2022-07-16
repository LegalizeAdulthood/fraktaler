// Fraktaler 3 -- fast deep escape time fractals
// Copyright (C) 2021,2022 Claude Heiland-Allen
// SPDX-License-Identifier: AGPL-3.0-only

#include <getopt.h>

#include <SDL.h>

#include "colour.h"
#include "engine.h"
#include "main.h"
#include "source.h"
#include "version.h"

param par;

std::string pref_path, default_source_path, default_persistence_path, default_wisdom_path;

void initialize_paths()
{
  pref_path = "";
  if (! SDL_Init(0))
  {
    char *p = SDL_GetPrefPath("uk.co.mathr", "fraktaler-3");
    if (p)
    {
      pref_path = std::string(p);
      SDL_free(p);
    }
    SDL_Quit();
  }
  default_source_path = source_filename;
  default_persistence_path = pref_path + "persistence.f3.toml";
  default_wisdom_path = pref_path + "fraktaler-3.wisdom";
}

int load_wisdom(const char *wisdom) // FIXME
{
  return 0;
}

int generate_wisdom(const char *wisdom) // FIXME
{
  return 0;
}

int print_mode_error(const char *progname)
{
  std::fprintf(stderr, "%s: error: exactly one mode of operation must be specified\n", progname);
  return 1;
}

int print_help(const char *progname)
{
  std::fprintf
    ( stdout
    , "usage:\n"
      "  %s [mode] [flags ...] [inputfile [inputfile ...]]\n"
      "modes of operation:\n"
      "  -h, --help                print this message and exit\n"
      "  -V, --version             print version information and exit\n"
      "  -i, --interactive         interactive graphical user interface\n"
      "  -b, --batch               command line batch processing\n"
      "  -W, --generate-wisdom     benchmark hardware for optimal efficiency\n"
      "  -S, --export-source       export this program's source code\n"
      "                            default location: %s\n"
      "flags:\n"
      "  -v, --verbose             increase verbosity\n"
      "  -q, --quiet               decrease verbosity\n"
      "  -p, --persistence file    path to persist state\n"
      "                            default location: %s\n"
      "  -P, --no-persistence      don't persist state\n"
      "  -w, --wisdom file         path to wisdom\n"
      "                            default location: %s\n"
      "input files are merged in command line order\n"
    , progname
    , default_source_path.c_str()
    , default_persistence_path.c_str()
    , default_wisdom_path.c_str()
    );
  return 0;
}

int print_version()
{
  std::fprintf(stdout, "%s", version().c_str());
  return 0;
}

int export_source()
{
  if (write_source(std::filesystem::path(source_filename)))
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

int print_generate_wisdom_error(const char *progname, const char *wisdom)
{
  std::fprintf(stderr, "%s: error: could not generate wisdom file '%s'\n", progname, wisdom);
  return 1;
}

int print_load_wisdom_error(const char *progname, const char *wisdom)
{
  std::fprintf(stderr, "%s: error: could not load wisdom file '%s'\n", progname, wisdom);
  return 1;
}

int print_load_parameter_error(const char *progname, const char *parameter)
{
  std::fprintf(stderr, "%s: error: could not load parameter file '%s'\n", progname, parameter);
  return 1;
}

int main(int argc, char **argv)
{
  // initialize environment
  initialize_paths();
  // modes of operation  
  bool
    do_help = false,
    do_version = false,
    do_interactive = false,
    do_batch = false,
    do_wisdom = false,
    do_source = false;
  int verbosity = 1;
  const char *wisdom = ""; // empty string -> use default app path
  const char *persistence = ""; // empty string -> use default app path
  // parse arguments
  while (1)
  {
    static struct option long_options[] =
      { { "version", no_argument, 0, 'V' }
      , { "help", no_argument, 0, 'h' }
      , { "verbose", no_argument, 0, 'v' }
      , { "quiet", no_argument, 0, 'q' }
      , { "interactive", no_argument, 0, 'i' }
      , { "batch", no_argument, 0, 'b' }
      , { "persistence", required_argument, 0, 'p' }
      , { "no-persistence", no_argument, 0, 'P' }
      , { "wisdom", required_argument, 0, 'w' }
      , { "generate-wisdom", no_argument, 0, 'W' }
      , { "export-source", no_argument, 0, 'S' }
      , { 0, 0, 0, 0 }
      };
    int option_index = 0;
    int c = getopt_long(argc, argv, "Vhvqibp:Pw:W:S", long_options, &option_index);
    if (c == -1)
    {
      break;
    }
    switch (c)
    {
      case 'V':
        do_version = true;
        break;
      case 'h':
        do_help = true;
        break;
      case 'v':
        verbosity++;
        break;
      case 'q':
        verbosity--;
        break;
      case 'i':
        do_interactive = true;
        break;
      case 'b':
        do_batch = true;
        break;
      case 'p':
        persistence = optarg;
        break;
      case 'P':
        persistence = nullptr;
        break;
      case 'w':
        wisdom = optarg;
        break;
      case 'W':
        wisdom = optarg;
        do_wisdom = true;
        break;
      case 'S':
        do_source = true;
        break;
    }
  }
  int count = int(do_help) + int(do_version) + int(do_source) + int(do_wisdom) + int(do_interactive) + int(do_batch);
  if (count == 0)
  {
    do_interactive = true;
    count++;
  }
  if (count != 1)
  {
    return print_mode_error(argv[0]);
  }
  if (std::string(wisdom) == std::string(""))
  {
    wisdom = default_wisdom_path.c_str();
  }
  if (persistence && std::string(persistence) == std::string(""))
  {
    persistence = default_persistence_path.c_str();
  }
  if (do_help)
  {
    return print_help(argv[0]);
  }
  if (do_version)
  {
    return print_version();
  }
  if (do_source)
  {
    return export_source();
  }
  // initialize rendering engine
  populate_number_type_wisdom();
  colours_init();
  if (do_wisdom)
  {
    return generate_wisdom(wisdom);
  }
  if (do_batch || do_interactive)
  {
    // load wisdom
    if (load_wisdom(wisdom))
    {
      // load failed, try to generate and retry
      if (generate_wisdom(wisdom))
      {
        return print_generate_wisdom_error(argv[0], wisdom);
      }
      if (load_wisdom(wisdom))
      {
        return print_load_wisdom_error(argv[0], wisdom);
      }
    }
    // load persistence
    if (persistence)
    {
      try
      {
        par.load_toml(persistence);
      }
      catch (...)
      {
        // ignore persistence load failure
      }
    }
    // load parameters
    for (int arg = optind; arg < argc; ++arg)
    {
      try
      {
        par.load_toml(argv[arg]);
      }
      catch (...)
      {
        return print_load_parameter_error(argv[0], argv[arg]);
      }
    }
    if (do_batch)
    {
      // FIXME unify batch modes
#ifdef HAVE_OPENCL
      if (par.p.opencl.platform != -1)
      {
#ifdef HAVE_CLEW
        if (clewInit())
        {
          std::fprintf(stderr, "%s: error: clewInit() failed\n", argv[0]);
          return 1;
        }
#endif
        return batch_cl(verbosity);
      }
      else
#endif
      {
        return batch_cli(verbosity);
      }
    }
    else
    {
      return gui(argv[0], verbosity, persistence);
    }
  }
  // unreachable
  return 1;
}
