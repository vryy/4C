/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_surface_evaluation_data.H"


/**
 *
 */
template <typename scalar_type>
void GEOMETRYPAIR::LineToSurfaceEvaluationData<scalar_type>::Reset()
{
  // Call reset on the base method.
  LineTo3DEvaluationData::Reset();
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::LineToSurfaceEvaluationData<double>;
