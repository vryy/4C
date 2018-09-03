/*---------------------------------------------------------------------*/
/*!
\file cut_tolerance.H

\brief Central file for !all! Tolerances used in the cut

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_TOLERANCE_H
#define CUT_TOLERANCE_H

#define BASICTOL 1e-14  // basic tolerance for the cut, limited by double precicion!

// derived tolerances
//--> do not change the relation of these tolerances from the basictolerance, if there isn't a
//reasonable argument to do that!!!

// to check l2-norm of vector with 3-entries of magnitude 1.0!
#define LINSOLVETOL (BASICTOL / sqrt(3) * 10.0)
// sqrt(3) because of l2 norm with 3 values in vector
// 10 because scaled with norminf(x) || if norminf(x)=1.0 --> for double precicion 1.(15 digits)
// REMARK: using this tolerance for increments in reference coordinate system norminf(x)=1!!!

// to define general tolerance for position represented by a vector of doubles!
#define POSITIONTOL (BASICTOL * 10.0)


#define REFERENCETOL \
  (BASICTOL * 10.0)  // tolerance to check coordinates in the relative coordinate system (assumtion
                     // is that abs(max(value)) = 1.0)!

// Not in use anymore.
#define VOLUMETOL \
  1e-14  // JUST FOR TESSELATION: the tolerance when created integration cells are discarded, the
         // boundary cell however will be created

// special tolerances for some special cut situations --> we should define them relative to our
// BASICTOL!!!

// tolerance to identify points of volumecells between processors
#define PARALLEL_COORD_TOL 1e-10

// Special tolerance to distinguish planar facets from curved ones.
#define PLANARTOL 1e-6

#define BOXOVERLAP 1e-7

// if cross product constructed from any 3 points are less than this,
// then these points are on the same line
#define TOL_POINTS_ON_LINE 1e-14

#define TOL_POINTS_ON_LINE_FOR_DELETING 1e-11

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
// in such case directdivergence may produce negative volume, because it works only for closed
// geometries if volume of cell is less than this, we do not generate gauss rules
#define REF_VOL_DIRDIV 1e-26

#define SELF_CUT_POS_TOL 1e-10

// If the area of boundarycell is less than this value, then it is eliminated from boundary
// integration Such small area can lead to zero or a very small negative value in the boundary cell
// mapping's Jacobian This is only for the case of Direct Divergence!
#define REF_AREA_BCELL 1e-13  // Make it viable for GLOBAL!!!

// If the distance of the line between the main Gauss points and its projected point is less than
// this value, then this main Gauss point is not considered
#define JAC_LINE_TOL 1e-14

// Defines the order of the Gaussrule used for the direct divergence.
//  To make it fit also the DivergenceCells have to be fully triangulated for the mesh as is done in
//  the LS. See the DirectDiv algorithm in facetintegration.
#define DIRECTDIV_GAUSSRULE 7  // 7 for test cases

// Construct the integration rules directly on the local coordinate system of the element
#define LOCAL

// utilized quadrature compression (for Details see paper by Sudhakar et al.)
//#define QUADCOMP


#endif
