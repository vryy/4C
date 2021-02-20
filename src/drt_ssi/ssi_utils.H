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
  class CouplingSlaveConverter;
}  // namespace ADAPTER

namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace LINALG
{
  class BlockSparseMatrixBase;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace SSI
{
  class SSIMono;

  namespace UTILS
  {
    //! Modification of time parameter list for problem with different time step size
    void ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
        Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn);

    //! check, if nodes of input conditions equal s2i conditions
    void CheckConsistencyWithS2IMeshtyingCondition(
        const std::vector<DRT::Condition*>& conditionsToBeTested,
        Teuchos::RCP<DRT::Discretization>& structdis);

    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    int CheckTimeStepping(double dt1, double dt2);

    //! setup coupling adapter for structure interfaces from ssi condition for line coupling
    Teuchos::RCP<ADAPTER::Coupling> SetupInterfaceCouplingAdapterStructure3DomainIntersection(
        Teuchos::RCP<DRT::Discretization> structdis,
        const std::string& conditionname_3_domain_intersection);

    //! setup coupling adapter for structure interfaces from ssi condition for surface coupling
    Teuchos::RCP<ADAPTER::Coupling> SetupInterfaceCouplingAdapterStructure(
        Teuchos::RCP<DRT::Discretization> structdis, bool meshtying_3_domain_intersection,
        const std::string& conditionname_coupling,
        const std::string& conditionname_3_domain_intersection);

    //! Build a vector of condition with slave-master pairs ordered by coupling ID
    std::vector<std::pair<DRT::Condition* const, DRT::Condition* const>> BuildSlaveMasterPairing(
        const std::vector<DRT::Condition*>& conditions);

    //! clone scatra specific parameters for solver of manifold. Add manifold specific parameters
    Teuchos::ParameterList CloneScaTraManifoldParams(const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& sublist_manifold_params, const Epetra_Comm& comm);

    //! create MultiMapExtractor with dofs on SSISurfaceManifold and others
    Teuchos::RCP<LINALG::MultiMapExtractor> CreateManifoldMultiMapExtractor(
        Teuchos::RCP<DRT::Discretization> dis);


    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds all subblocks of system matrices and system matrix for SSTI simulations
    class SSIMatrices
    {
     public:
      SSIMatrices(
          const SSI::SSIMono& ssi_mono_algorithm, Teuchos::RCP<Epetra_Map> interface_map_scatra);

      Teuchos::RCP<LINALG::SparseOperator> SystemMatrix() { return systemmatrix_; }

      //! return subblocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureDomain()
      {
        return scatrastructuredomain_;
      }
      Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureInterface()
      {
        return scatrastructureinterface_;
      }
      Teuchos::RCP<LINALG::SparseOperator> StructureScaTraDomain()
      {
        return structurescatradomain_;
      }
      Teuchos::RCP<LINALG::SparseOperator> ScaTraManifoldStructureDomain()
      {
        return scatramanifoldstructuredomain_;
      }
      //@}

     private:
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> SetupBlockMatrix(
          Teuchos::RCP<const LINALG::MultiMapExtractor> row_map,
          Teuchos::RCP<const LINALG::MultiMapExtractor> col_map);

      Teuchos::RCP<LINALG::SparseMatrix> SetupSparseMatrix(
          const Teuchos::RCP<const Epetra_Map> row_map);

      Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;
      //! subblocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;
      Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;
      Teuchos::RCP<LINALG::SparseOperator> structurescatradomain_;
      Teuchos::RCP<LINALG::SparseOperator> scatramanifoldstructuredomain_;
      //@}
    };

    class SSISlaveSideConverter
    {
     public:
      SSISlaveSideConverter(Teuchos::RCP<ADAPTER::Coupling> icoup_structure,
          Teuchos::RCP<ADAPTER::Coupling> icoup_structure_3_domain_intersection,
          bool meshtying_3_domain_intersection);

      //! Slave structure converter for structure coupling
      ADAPTER::CouplingSlaveConverter& InterfaceCouplingAdapterStructureSlaveConverter() const
      {
        return *icoup_structure_slave_converter_;
      }

      //! Slave structure converter for structure coupling with 3 domain intersection
      ADAPTER::CouplingSlaveConverter&
      InterfaceCouplingAdapterStructureSlaveConverter3DomainIntersection() const
      {
        return *icoup_structure_slave_converter_3_domain_intersection_;
      }

     private:
      //! Slave structure converter for structure coupling
      Teuchos::RCP<ADAPTER::CouplingSlaveConverter> icoup_structure_slave_converter_;

      //! Slave structure converter for structure coupling with 3 domain intersection
      Teuchos::RCP<ADAPTER::CouplingSlaveConverter>
          icoup_structure_slave_converter_3_domain_intersection_;
    };

  }  // namespace UTILS
}  // namespace SSI

#endif /* SSI_UTILS_H_ */
