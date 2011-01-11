#ifndef CUT_TETGEN_H
#define CUT_TETGEN_H

#ifdef QHULL
#define TETLIBRARY
#undef PI
#include <tetgen.h>
#undef DOT
#undef REAL
#endif

#endif
