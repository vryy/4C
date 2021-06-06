/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_UTILS_H_
#define SSI_UTILS_H_

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include "unordered_map"
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
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
  enum class MatrixType;
}  // namespace LINALG

namespace SSI
{
  class SSIMono;
  enum class Subproblem;

  namespace UTILS
  {
    // forward declaration
    class SSIMaps;

    //! Modification of time parameter list for problem with different time step size
    void ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
        Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn);

    //! check for a consistent input file definition of the SSIInterfaceContact condition
    void CheckConsistencyOfSSIInterfaceContactCondition(
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
        const Teuchos::ParameterList& sublist_manifold_params);

    //! create MultiMapExtractor with dofs on SSISurfaceManifold and others
    Teuchos::RCP<LINALG::MultiMapExtractor> CreateManifoldMultiMapExtractor(
        Teuchos::RCP<DRT::Discretization> dis);

    //! modify scatra parameters for ssi specific values
    Teuchos::ParameterList ModifyScaTraParams(const Teuchos::ParameterList& scatraparams);


    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds all sub blocks of system matrices and system matrix for SSI simulations
    class SSIMatrices
    {
     public:
      /*!
       * @brief constructor
       *
       * @param[in] ssi_maps            pointer to the ssi maps object containing all relevant maps
       * @param[in] ssi_matrixtype      the ssi matrix type
       * @param[in] scatra_matrixtype   the scalar transport matrix type
       * @param[in] is_scatra_manifold  flag indicating if a scatra manifold is used
       */
      SSIMatrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps,
          LINALG::MatrixType ssi_matrixtype, LINALG::MatrixType scatra_matrixtype,
          bool is_scatra_manifold);

      //! call complete on the scalar transport manifold - structure off-diagonal matrix
      void CompleteScaTraManifoldStructureMatrix();

      //! call complete on the scalar transport - structure off-diagonal matrix
      void CompleteScaTraStructureMatrix();

      //! call complete on the structure - scalar transport off-diagonal matrix
      void CompleteStructureScaTraMatrix();

      //! method that clears all ssi matrices
      void ClearMatrices();

      //! return the system matrix
      Teuchos::RCP<LINALG::SparseOperator> SystemMatrix() { return system_matrix_; }

      //! return sub blocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> ScaTraMatrix() { return scatra_matrix_; }
      Teuchos::RCP<LINALG::SparseOperator> ScaTraManifoldStructureMatrix()
      {
        return scatramanifold_structure_matrix_;
      }
      Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureMatrix()
      {
        return scatra_structure_matrix_;
      }
      Teuchos::RCP<LINALG::SparseOperator> StructureScaTraMatrix()
      {
        return structure_scatra_matrix_;
      }
      Teuchos::RCP<LINALG::SparseMatrix> StructureMatrix() { return structure_matrix_; }
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
       * @brief initialize the scatra-structure interaction main-diagonal matrices
       *
       * @param[in] ssi_maps            pointer to the ssi maps object containing all relevant maps
       */
      void InitializeMainDiagMatrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

      /*!
       * @brief initialize the scatra-structure interaction off-diagonal matrices
       *
       * @param[in] ssi_maps            pointer to the ssi maps object containing all relevant maps
       */
      void InitializeOffDiagMatrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps);

      /*!
       * @brief initialize the system matrix
       *
       * @param[in] ssi_maps         pointer to the ssi maps object containing all relevant maps
       * @param[in] ssi_matrixtype   the ssi matrix type
       */
      void InitializeSystemMatrix(
          Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, LINALG::MatrixType ssi_matrixtype);

      //! flag indicating if we have a scatra manifold
      const bool is_scatra_manifold_;

      //! matrix type of scatra matrix
      const LINALG::MatrixType scatra_matrixtype_;

      //! the scalar transport dof row map
      Teuchos::RCP<const Epetra_Map> scatra_dofrowmap_;

      //! the scalar transport manifold dof row map
      Teuchos::RCP<const Epetra_Map> scatramanifold_dofrowmap_;

      //! the structure dof row map
      Teuchos::RCP<const Epetra_Map> structure_dofrowmap_;

      //! system matrix
      Teuchos::RCP<LINALG::SparseOperator> system_matrix_;
      //! sub blocks of system matrix
      //@{
      Teuchos::RCP<LINALG::SparseOperator> scatra_matrix_;
      Teuchos::RCP<LINALG::SparseOperator> scatramanifold_structure_matrix_;
      Teuchos::RCP<LINALG::SparseOperator> scatra_structure_matrix_;
      Teuchos::RCP<LINALG::SparseOperator> structure_scatra_matrix_;
      Teuchos::RCP<LINALG::SparseMatrix> structure_matrix_;
      //@}
    };

    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    //! sets up and holds the system residuals and increment for SSI simulations
    class SSIVectors
    {
     public:
      //! constructor
      explicit SSIVectors(const SSI::SSIMono& ssi_mono);

      //! clear the increment vector
      void ClearIncrement();

      //! clear all residual vectors
      void ClearResiduals();

      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> Increment() { return increment_; }

      //! residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }

      //! residual vector on right-hand side of scalar transport system
      Teuchos::RCP<Epetra_Vector> ScatraResidual() { return scatra_residual_; }

      //! residual vector on right-hand side of structure system
      Teuchos::RCP<Epetra_Vector> StructureResidual() { return structure_residual_; }

     private:
      //! global increment vector for Newton-Raphson iteration
      Teuchos::RCP<Epetra_Vector> increment_;

      //! residual vector on right-hand side of global system of equations
      Teuchos::RCP<Epetra_Vector> residual_;

      //! residual vector on right-hand side of scalar transport system
      Teuchos::RCP<Epetra_Vector> scatra_residual_;

      //! residual vector on right-hand side of structure system
      Teuchos::RCP<Epetra_Vector> structure_residual_;
    };

    /*---------------------------------------------------------------------------------*
     *---------------------------------------------------------------------------------*/
    class SSIMaps
    {
     public:
      //! constructor
      explicit SSIMaps(const SSI::SSIMono& ssi_mono_algorithm);

      //! get position within global dof map for specific sub problem
      static int GetProblemPosition(Subproblem subproblem);

      //! block maps of all sub problems organized in std map
      std::unordered_map<Subproblem, Teuchos::RCP<const LINALG::MultiMapExtractor>>
      BlockMapsSubProblems() const
      {
        return block_maps_sub_problems_;
      }

      //! the multi map extractor of the scalar transport field
      Teuchos::RCP<const LINALG::MultiMapExtractor> MultiMapExtractorScaTra() const;

      //! the multi map extractor of the scalar transport on manifold field
      Teuchos::RCP<const LINALG::MultiMapExtractor> MultiMapExtractorScaTraManifold() const;

      //! the multi map extractor of the structure field
      Teuchos::RCP<const LINALG::MultiMapExtractor> MultiMapExtractorStructure() const;

      //! map extractor associated with blocks of global system matrix
      Teuchos::RCP<const LINALG::MultiMapExtractor> BlockMapsSystemMatrix() const
      {
        return block_maps_systemmatrix_;
      }

      //! all dofs of the SSI algorithm
      Teuchos::RCP<const Epetra_Map> MapSystemMatrix() const { return map_system_matrix_; }

      //! global map extractor (0: scalar transport, 1: structure, [2: scatra manifold])
      Teuchos::RCP<const LINALG::MultiMapExtractor> MapsSubProblems() const
      {
        return maps_sub_problems_;
      }

      //! the scalar transport dof row map
      Teuchos::RCP<const Epetra_Map> ScaTraDofRowMap() const;

      //! the scalar transport on manifolds dof row map
      Teuchos::RCP<const Epetra_Map> ScaTraManifoldDofRowMap() const;

      //! the structure dof row map
      Teuchos::RCP<const Epetra_Map> StructureDofRowMap() const;

      //! condensed structure dofs (inner and master side): other map
      //! slave side structure dofs: cond map
      Teuchos::RCP<const LINALG::MapExtractor> SlaveCondensedStructureMaps() const
      {
        return slave_condensed_structure_maps_;
      }

      //! slave side map in case of 3 domain mesh tying
      Teuchos::RCP<const Epetra_Map> SlaveSide3DomainMeshytyingMap() const
      {
        return slave_side_3_domain_meshtying_map_;
      }

     private:
      //! create and check the block maps of all sub problems
      void CreateAndCheckBlockMapsSubProblems(const SSIMono& ssi_mono_algorithm);

      //! block maps of all sub problems organized in std map
      std::unordered_map<Subproblem, Teuchos::RCP<const LINALG::MultiMapExtractor>>
          block_maps_sub_problems_;

      //! map extractor associated with blocks of global system matrix
      Teuchos::RCP<const LINALG::MultiMapExtractor> block_maps_systemmatrix_;

      //! all dofs of the SSI algorithm
      Teuchos::RCP<const Epetra_Map> map_system_matrix_;

      //! global map extractor (0: scalar transport, 1: structure, [2: scatra manifold])
      Teuchos::RCP<const LINALG::MultiMapExtractor> maps_sub_problems_;

      //! condensed structure dofs (inner and master side): other map
      //! slave side structure dofs: cond map
      Teuchos::RCP<const LINALG::MapExtractor> slave_condensed_structure_maps_;

      //! slave side map in case of 3 domain mesh tying
      Teuchos::RCP<const Epetra_Map> slave_side_3_domain_meshtying_map_;
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
