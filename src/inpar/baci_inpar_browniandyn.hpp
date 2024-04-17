/*-----------------------------------------------------------*/
/*! \file
\brief input parameter for Brownian dynamics simulation


\level 2

*/
/*-----------------------------------------------------------*/
#ifndef FOUR_C_INPAR_BROWNIANDYN_HPP
#define FOUR_C_INPAR_BROWNIANDYN_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace BROWNIANDYN
  {
    /// the way how damping coefficient values for beams are specified
    enum BeamDampingCoefficientSpecificationType
    {
      cylinder_geometry_approx,
      input_file,
      vague
    };

    /// set the brownian dynamic parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace BROWNIANDYN

}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
