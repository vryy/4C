/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_UTILS_H_
#define SSI_UTILS_H_

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
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

namespace LINALG
{
  class SparseOperator;
}

namespace SSI
{
  class SSI_Mono;

  namespace Utils
  {
    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    int CheckTimeStepping(double dt1, double dt2);

    //! Modification of time parameter list for problem with different time step size
    void ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
        Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn);

    //! check, if nodes of input conditions equal s2i conditions
    void CheckConsistencyWithS2IMeshtyingCondition(
        std::vector<DRT::Condition*>& conditionsToBeTested,
        Teuchos::RCP<DRT::Discretization>& structdis);

    //! setup coupling adapter for structure interfaces from s2i condition
    Teuchos::RCP<ADAPTER::Coupling> SetupInterfaceCouplingAdapterStructure(
        Teuchos::RCP<DRT::Discretization> structdis);

    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds all subblocks of system matrices and system matrix for SSTI simulations
    class SSISubMatrices
    {
     public:
      SSISubMatrices(
          const SSI::SSI_Mono& ssi_mono_algorithm, Teuchos::RCP<Epetra_Map> interface_map_scatra);

      Teuchos::RCP<LINALG::SparseOperator> SystemMatrix() { return systemmatrix_; };

      //! return subblocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureDomain()
      {
        return scatrastructuredomain_;
      };
      Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureInterface()
      {
        return scatrastructureinterface_;
      };
      Teuchos::RCP<LINALG::SparseOperator> StructureScaTraDomain()
      {
        return structurescatradomain_;
      };
      //@}

     private:
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;
      //! subblocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;
      Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;
      Teuchos::RCP<LINALG::SparseOperator> structurescatradomain_;
      //@}
    };

  };  // namespace Utils
}  // namespace SSI

#endif /* SSI_UTILS_H_ */
