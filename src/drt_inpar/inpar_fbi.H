/*-----------------------------------------------------------*/
/*! \file
\file inpar_fbi.H

\brief input parameter for Fluid-Beam Interaction

\maintainer Nora Hagmeyer

\level 3

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_FBI_H
#define INPAR_FBI_H

#include "inpar_parameterlist_utils.H"
#include <boost/algorithm/string/predicate.hpp>

#include "../drt_lib/drt_dserror.H"

#include "../drt_fem_general/drt_utils_integration.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace FBI
  {
    /// Constraint enforcement for beam to fluid meshtying.
    enum class BeamToFluidConstraintEnforcement
    {
      //! Default value.
      none,
      //! Penalty method.
      penalty
    };

    /// Discretization approach for beam to fluid meshtying.
    enum class BeamToFluidDiscretization
    {
      none,                    //< Default value
      gauss_point_to_segment,  //< Gauss point to segment approach
      //! Gauss point on cross section method.
      //      gauss_point_cross_section,
      //! Mortar type coupling.
      //      mortar
    };

    //! Map number of gauss points to 1D gauss rule
    static inline DRT::UTILS::GaussRule1D IntToGaussRule1D(const int n_gauss_points)
    {
      switch (n_gauss_points)
      {
        case 1:
          return DRT::UTILS::GaussRule1D::intrule_line_1point;
        case 2:
          return DRT::UTILS::GaussRule1D::intrule_line_2point;
        case 3:
          return DRT::UTILS::GaussRule1D::intrule_line_3point;
        case 4:
          return DRT::UTILS::GaussRule1D::intrule_line_4point;
        case 5:
          return DRT::UTILS::GaussRule1D::intrule_line_5point;
        case 6:
          return DRT::UTILS::GaussRule1D::intrule_line_6point;
        case 7:
          return DRT::UTILS::GaussRule1D::intrule_line_7point;
        case 8:
          return DRT::UTILS::GaussRule1D::intrule_line_8point;
        case 9:
          return DRT::UTILS::GaussRule1D::intrule_line_9point;
        case 10:
          return DRT::UTILS::GaussRule1D::intrule_line_10point;
        case 16:
          return DRT::UTILS::GaussRule1D::intrule_line_16point;
        case 20:
          return DRT::UTILS::GaussRule1D::intrule_line_20point;
        case 32:
          return DRT::UTILS::GaussRule1D::intrule_line_32point;
        case 50:
          return DRT::UTILS::GaussRule1D::intrule_line_50point;
        default:
        {
          dserror("No Gauss rule defined for %d points", n_gauss_points);
          return DRT::UTILS::GaussRule1D::intrule1D_undefined;
        }
      }
    };

    /// set the beam interaction parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set beam interaction specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace FBI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif
