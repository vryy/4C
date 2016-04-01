/*----------------------------------------------------------------------*/
/*!
\file inpar_plasticity.H
\maintainer Alexander Seitz
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_PLASTICITY_H
#define INPAR_PLASTICITY_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for potential-based beam interaction
  namespace PLASTICITY
  {
    /// set the plasticity parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PLASTICITY

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_PLASTICITY_H
