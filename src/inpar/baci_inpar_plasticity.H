/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for plasticity
\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_PLASTICITY_H
#define BACI_INPAR_PLASTICITY_H

#include "baci_inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace PLASTICITY
  {
    /// set the plasticity parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PLASTICITY

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // INPAR_PLASTICITY_H
