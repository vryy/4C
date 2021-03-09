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

    //! check for a consistent input file definition of the SSIInterfaceContact condition
    void CheckConsistencyOfSSIInterfaceContactCondition(
        const std::vector<DRT::Condition*>& conditionsToBeTested,
        Teuchos::RCP<DRT::Discretization>& structdis);

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
    //! sets up and holds all sub blocks of system matrices and system matrix for SSI simulations
    class SSIMatrices
    {
     public:
      SSIMatrices(const SSI::SSIMono& ssi_mono_algorithm,
          Teuchos::RCP<const Epetra_Map> interface_map_scatra);

      //! method that clears all ssi matrices
      void ClearMatrices();

      //! return the system matrix
      Teuchos::RCP<LINALG::SparseOperator> SystemMatrix() { return systemmatrix_; }

      //! return sub blocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> ScaTraManifoldStructureDomain()
      {
        return scatramanifoldstructuredomain_;
      }
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
      Teuchos::RCP<LINALG::SparseMatrix> StructureMatrix() { return structurematrix_; }
      //@}

      /*!
       * @brief set up a pointer to a block matrix
       *
       * @param[in] row_map  row map the block matrix is based on
       * @param[in] col_map  column map the block matrix is based on
       * @return pointer to block matrix
       */
      static Teuchos::RCP<LINALG::BlockSparseMatrixBase> SetupBlockMatrix(
          Teuchos::RCP<const LINALG::MultiMapExtractor> row_map,
          Teuchos::RCP<const LINALG::MultiMapExtractor> col_map);

      /*!
       * @brief set up a pointer to a sparse matrix
       *
       * @param[in] row_map  row map the sparse matrix is based on
       * @return pointer to sparse matrix
       */
      static Teuchos::RCP<LINALG::SparseMatrix> SetupSparseMatrix(
          const Teuchos::RCP<const Epetra_Map> row_map);

     private:
      /*!
       * @brief set up and return the scatra interface block map
       *
       * @param[in] ssi_mono_algorithm    scatra-structure interaction monolithic time integrator
       * @param[in] interface_map_scatra  slave- and master-side of scatra interface map
       * @return scatra interface block map
       */
      Teuchos::RCP<const LINALG::MultiMapExtractor> GetScaTraInterfaceBlockMap(
          const SSI::SSIMono& ssi_mono_algorithm,
          Teuchos::RCP<const Epetra_Map> interface_map_scatra) const;

      /*!
       * @brief initialize the scatra-structure interaction off-diagonal matrices
       *
       * @param[in] ssi_mono_algorithm    scatra-structure interaction monolithic time integrator
       * @param[in] interface_map_scatra  slave- and master-side of scatra interface map
       */
      void InitializeOffDiagMatrices(const SSI::SSIMono& ssi_mono_algorithm,
          Teuchos::RCP<const Epetra_Map> interface_map_scatra);

      /*!
       * @brief initialize the system matrix
       *
       * @param[in] ssi_mono_algorithm  scatra-structure interaction monolithic time integrator
       */
      void InitializeSystemMatrix(const SSI::SSIMono& ssi_mono_algorithm);

      //! flag indicating if we have a scatra manifold
      const bool is_scatra_manifold_;
      //! flag indicating if we have ssi interface mesh tying
      const bool is_ssi_interface_meshtying_;
      //! system matrix
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;
      //! sub blocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> scatramanifoldstructuredomain_;
      Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;
      Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;
      Teuchos::RCP<LINALG::SparseOperator> structurescatradomain_;
      Teuchos::RCP<LINALG::SparseMatrix> structurematrix_;
      //@}
    };

    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds the system residuals and increment for SSI simulations
    class SSIVectors
    {
     public:
      //! constructor
      explicit SSIVectors(const SSI::SSIMono& ssi_mono_algorithm);

      //! clear the increment vector
      void ClearIncrement();

      //! clear all residual vectors
      void ClearResiduals();

      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> Increment() { return increment_; }

      //! residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }

      //! residual vector on right-hand side of structure system
      Teuchos::RCP<Epetra_Vector> StructureResidual() { return structure_residual_; }

     private:
      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> increment_;

      //! residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> residual_;

      //! residual vector on right-hand side of structure system
      Teuchos::RCP<Epetra_Vector> structure_residual_;
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
