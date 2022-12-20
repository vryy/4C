/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameter for geometry pairs.

\level 3
*/
// End doxygen header.


#ifndef INPAR_GEOMETRYPAIR_H
#define INPAR_GEOMETRYPAIR_H


#include "utils_integration.H"

// Forward declaration.
namespace Teuchos
{
  class ParameterList;

}  // namespace Teuchos


namespace INPAR
{
  namespace GEOMETRYPAIR
  {
    /**
     * Method to be used for line to 3D pairs.
     */
    enum class LineTo3DStrategy
    {
      none,                                                  ///< Default value.
      gauss_point_projection_without_boundary_segmentation,  ///< Project the Gauss points.
      gauss_point_projection_boundary_segmentation,  ///< Project the Gauss points, segment lines
                                                     ///< that poke out of the volume.
      gauss_point_projection_cross_section,  ///< Project the Gauss points on the beam surface.
      segmentation                           ///< Segment the line with the volume.
    };

    /**
     * Method to be used for normals in surface pairs.
     */
    enum class SurfaceNormals
    {
      standard,        ///< Use normals on the actual surface.
      extended_volume  ///< Use normals from extended volume projection.
    };

    /**
     * \brief Map number of gauss points to 1D gauss rule
     */
    DRT::UTILS::GaussRule1D IntToGaussRule1D(const int n_gauss_points);

    /**
     * \brief Set valid input parameters for line to 3D geometry pairs.
     *
     * This function adds all available input parameters to a given section. Call this function with
     * the parameter list of the input section you want the line to 3D parameters to be added to.
     *
     * @param (out) Parameter list to add the line to 3D parameters to.
     */
    void SetValidParametersLineTo3D(Teuchos::ParameterList& list);

    /**
     * \brief Set valid input parameters for line to surface geometry pairs.
     *
     * This function adds all available input parameters for line to surface pairs, that are not yet
     * added in the SetValidParametersLineTo3D function.
     *
     * @param (out) Parameter list to add the line to surface parameters to.
     */
    void SetValidParametersLineToSurface(Teuchos::ParameterList& list);

  }  // namespace GEOMETRYPAIR

}  // namespace INPAR


#endif
