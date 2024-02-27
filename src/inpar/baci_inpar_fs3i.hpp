/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for fs3i

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef BACI_INPAR_FS3I_HPP
#define BACI_INPAR_FS3I_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{
  namespace FS3I
  {
    /// Type of coupling strategy for SSI problems
    enum SolutionSchemeOverFields
    {
      fs3i_SequStagg,
      //    fs3i_SequStagg_SolidToScatra,
      fs3i_IterStagg,
      //    fs3i_IterStaggFixedRel_ScatraToSolid,
      //    fs3i_IterStaggFixedRel_SolidToScatra,
      //    fs3i_IterStaggAitken_ScatraToSolid,
      //    fs3i_IterStaggAitken_SolidToScatra
      //    fs3i_IterStaggAitkenIrons,
      //    fs3i_Monolithic
    };

    /// Type of coupling strategy between the two fields of the SSI problems
    enum VolumeCoupling
    {
      coupling_match,
      coupling_nonmatch
    };

    /// set the fs3i parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace FS3I

}  // namespace INPAR

BACI_NAMESPACE_CLOSE

#endif  // INPAR_FS3I_H
