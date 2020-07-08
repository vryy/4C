/*--------------------------------------------------------------------------*/
/*! \file

\brief utility class for  elastohydrodynamic lubrication (lubrication structure interaction)

\level 3


*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_EHL_EHL_UTILS_H_
#define SRC_DRT_EHL_EHL_UTILS_H_


#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

namespace EHL
{
  namespace Utils
  {
    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    int CheckTimeStepping(double dt1, double dt2);

    // Modification of time parameter list for problem with different time step size
    void ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ehlparams,
        Teuchos::ParameterList& lubricationdyn, Teuchos::ParameterList& sdyn);

  };  // namespace Utils

  //! prints the BACI EHL-logo on the screen
  void printlogo();

}  // namespace EHL


#endif /* SRC_DRT_EHL_EHL_UTILS_H_ */
