#ifndef CUT_TETGEN_H
#define CUT_TETGEN_H

#ifdef QHULL
#define TETLIBRARY
#undef PI
#ifdef TETGENINCLUDED
#include "../tetgen/tetgen.h"
#else
#include <tetgen.h>
#endif
#undef DOT
#endif

#endif
