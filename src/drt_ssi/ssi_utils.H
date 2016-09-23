/*!------------------------------------------------------------------------------------------------*
 \file ssi_utils.H

 \brief Utility methods for SSI

\level 1

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_UTILS_H_
#define SSI_UTILS_H_

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

namespace SSI{
namespace Utils{
/// Function for checking that the different time steps are a
/// multiplicative of each other
int CheckTimeStepping(double dt1, double dt2);

//! Modification of time parameter list for problem with different time step size
void ChangeTimeParameter(const Epetra_Comm& comm,
    Teuchos::ParameterList& ssiparams,
    Teuchos::ParameterList& scatradyn,
    Teuchos::ParameterList& sdyn);

};// namespace Utils
}// namespace SSI

#endif /* SSI_UTILS_H_ */
