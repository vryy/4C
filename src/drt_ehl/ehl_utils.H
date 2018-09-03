/*--------------------------------------------------------------------------*/
/*!
\file ehl_utils.H

\brief utility class for  elastohydrodynamic lubrication (lubrication structure interaction)

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
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
}  // namespace EHL


#endif /* SRC_DRT_EHL_EHL_UTILS_H_ */
