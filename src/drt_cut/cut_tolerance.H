#ifndef CUT_TOLERANCE_H
#define CUT_TOLERANCE_H

// Our input files enjoy lax accuracies. Thus, we need to have a large catch
// point radius during setup. Later on, however, the catch radius determines
// the numerical exactness we can reach. We need to use a rather small catch
// radius during the actual cutting.
//#define SETUPNODECATCHTOL 1e-6

// Tolerances that build on each other. We assume the numbers we compare with
// these are about 1. These tolerances build on each other. The largest one
// determines the exactness we can reach in out simulations.
#define MINIMALTOL  1e-11
#define TOLERANCE   1e-12
#define LINSOLVETOL 1e-13

// Special tolerance to distinguish planar facets from curved ones.
#define PLANARTOL   1e-6

#define BOXOVERLAP  1e-7

// Scale all points that go inside tetgen, so the tetgen internal tolerance
// does not merge nodes.
#define TETGENPOINTSCALE 1e5

#endif
