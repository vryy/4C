/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_UTILS_H_
#define SSI_UTILS_H_

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

namespace ADAPTER
{
  class Coupling;
}

namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace SSI
{
  namespace Utils
  {
    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    int CheckTimeStepping(double dt1, double dt2);

    //! Modification of time parameter list for problem with different time step size
    void ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
        Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn);

    void CheckConsistencyWithS2IMeshtyingCondition(
        std::vector<DRT::Condition*> conditionsToBeTested,
        Teuchos::RCP<DRT::Discretization>& structdis);

    Teuchos::RCP<ADAPTER::Coupling> SetupInterfaceCouplingAdapterStructure(
        Teuchos::RCP<DRT::Discretization> structdis);

  };  // namespace Utils
}  // namespace SSI

#endif /* SSI_UTILS_H_ */
