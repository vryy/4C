/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for plasticity
\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_PLASTICITY_H
#define INPAR_PLASTICITY_H

#include "inpar_parameterlist_utils.H"

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
#endif  // #ifndef INPAR_PLASTICITY_H
