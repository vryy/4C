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
#define MINIMALTOL  1e-10
#define TOLERANCE   1e-12
#define LINSOLVETOL 1e-13
#define VOLUMETOL   1e-14  // the tolerance when created integration cells are discarded, the boundary cell however will be created

// special tolerances for some special cut situations

// tolerance in cut_kernel to decide if there is a cut between side and line (dependent on det of system)
// or if these lie on each other or are parallel
#define PARALLEL_DET_TOL 1e-16 // has to be smaller than 1e-14 -> cut_test

// tolerance to identify points of volumecells between processors
#define PARALLEL_COORD_TOL 1e-10

// Special tolerance to distinguish planar facets from curved ones.
#define PLANARTOL   1e-6

#define BOXOVERLAP  1e-7

// if cross product constructed from any 3 points are less than this,
// then these points are on the same line
#define TOL_POINTS_ON_LINE 1e-11

// In equation of plane ax+by+cz=d for a facet, this tolerance is used to check
// whether any of the four coefficients are zero
#define TOL_EQN_PLANE 1e-12

// For moment fitting, during line integration of base functions, if normal componenet of
// lines of less than this, they are considered zero
#define TOL_LINE_NORMAL 1e-12

// Obtaining  internal Gauss points in directdivergence involves division by
// coefficent 'a' in equation of plane. A facet whose 'a' is more than this tolerance is
// chosen as reference plane to avoid a large number (and hence large error)
#define REF_PLANE_DIRDIV 1e-3

// We may get volumecells with all facets on a plane
// in such case directdivergence may produce negative volume, because it works only for closed geometries
// if volume of cell is less than this, we do not generate gauss rules
#define REF_VOL_DIRDIV 1e-6

#endif
