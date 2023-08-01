/*---------------------------------------------------------------------*/
/*! \file

\brief tutorial input parameters


\level 2

*/
/*---------------------------------------------------------------------*/

// All header files should have #define guards to prevent multiple inclusion.
// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef BACI_TUTORIAL_INPAR_H
#define BACI_TUTORIAL_INPAR_H

#include "baci_inpar_parameterlist_utils.H"

namespace INPAR
{
  namespace TUTORIAL
  {
    /*----------------------------------------------------------------------*
     | Coupling Methods                                                     |
     *----------------------------------------------------------------------*/
    enum _TUTORIAL_TYPE_
    {
      nonlinear_truss,
      partitioned_fixed_point_scheme,
      struct_ale_coupling
    };


    /// set the tutorial parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace TUTORIAL
}  // namespace INPAR
/*----------------------------------------------------------------------*/
#endif  // TUTORIAL_INPAR_H
