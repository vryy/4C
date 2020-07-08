/*---------------------------------------------------------------------*/
/*! \file

\brief tutorial input parameters


\level 2

*/
/*---------------------------------------------------------------------*/

// All header files should have #define guards to prevent multiple inclusion.
// The format of the symbol name should be <PROJECT>_<PATH>_<FILE>_H_.
#ifndef INPAR_TUTORIAL_H_
#define INPAR_TUTORIAL_H_

#include "../drt_inpar/inpar_parameterlist_utils.H"

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
#endif /* INPAR_TUTORIAL_H_ */
