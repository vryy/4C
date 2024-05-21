/*----------------------------------------------------------------------*/
/*! \file
 \brief Provides utilities to incorporate transport on manifolds into SSI, mainly
 - flux between ScaTra and ScaTra on manifolds incl. coupling matrices
 - mesh tying between manifold fields
 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_MANIFOLD_UTILS_HPP
#define FOUR_C_SSI_MANIFOLD_UTILS_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_ssi_utils.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>
#include <set>

FOUR_C_NAMESPACE_OPEN

namespace CORE::ADAPTER
{
  class Coupling;
  class CouplingSlaveConverter;
  class CouplingMasterConverter;
}  // namespace CORE::ADAPTER

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
}

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace CORE::LINALG
{
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
  enum class MatrixType;
}  // namespace CORE::LINALG

namespace SSI
{

  class SSIMono;

  //! type of coupling/system matrix block
  enum class BlockMatrixType
  {
    ManifoldScaTra,     //! derivative of flux into/out of manifold w.r.t. scatra
    ManifoldStructure,  //! derivative of flux into/out of manifold w.r.t. structure
    ScaTraManifold,     //! derivative of flux into/out of scatra w.r.t. manifold
    ScaTraStructure,    //! derivative of flux into/out of scatra w.r.t. structure
    SysMatManifold,     //! derivative of flux into/out of manifold w.r.t. manifold
    SysMatScaTra        //! derivative of flux into/out of scatra w.r.t. scatra
  };

  //! holds everything to evaluate coupling between scatra and scatra manifold
  class ManifoldScaTraCoupling
  {
   public:
    ManifoldScaTraCoupling(Teuchos::RCP<DRT::Discretization> manifolddis,
        Teuchos::RCP<DRT::Discretization> scatradis,
        CORE::Conditions::Condition* condition_manifold,
        CORE::Conditions::Condition* condition_kinetics, int ndof_per_node);

    //! Check if graph of matrix has changed compared to last evaluation. Afterwards, store new size
    //! of graph
    //! \param block  graph of which matrix should be checked
    //! \param size   new size of graph
    //! \return       does old size of graph match new size?
    bool CheckAndSetSizeOfMatrixGraph(BlockMatrixType block, int size);

    //! Kinetics condition on scatra dis
    CORE::Conditions::Condition* ConditionKinetics() const { return condition_kinetics_; }

    //! manifold condition
    CORE::Conditions::Condition* ConditionManifold() const { return condition_manifold_; }

    //! coupling adapter between manifold (slave) and scatra (master)
    Teuchos::RCP<CORE::ADAPTER::Coupling> CouplingAdapter() const { return coupling_adapter_; }

    //! inverse of thickness of manifold
    double InvThickness() const { return inv_thickness_; }

    //! condition ID of manifold condition
    int ManifoldConditionID() const { return manifold_condition_id_; }

    //! from master to slave side
    Teuchos::RCP<CORE::ADAPTER::CouplingMasterConverter> MasterConverter() const
    {
      return master_converter_;
    }

    //! condition ID of kinetics condition
    int KineticsConditionID() const { return kinetics_condition_id_; }

    //! Map exctractor for dofs in this manifold condition
    Teuchos::RCP<CORE::LINALG::MapExtractor> ManifoldMapExtractor() const
    {
      return manifold_map_extractor_;
    }

    //! Map exctractor for dofs in this kinetics condition
    Teuchos::RCP<CORE::LINALG::MapExtractor> ScaTraMapExtractor() const
    {
      return scatra_map_extractor_;
    }

   private:
    //! Kinetics condition on scatra dis
    CORE::Conditions::Condition* condition_kinetics_;

    //! manifold condition
    CORE::Conditions::Condition* condition_manifold_;

    //! coupling adapter between manifold (slave) and scatra (master)
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupling_adapter_;

    //! inverse of thickness of manifold
    const double inv_thickness_;

    //! condition ID of manifold condition
    const int manifold_condition_id_;

    //! condition ID of kinetics condition
    const int kinetics_condition_id_;

    //! Map exctractor for dofs in this manifold condition
    Teuchos::RCP<CORE::LINALG::MapExtractor> manifold_map_extractor_;

    //! Master converter for scatra - manifold coupling
    Teuchos::RCP<CORE::ADAPTER::CouplingMasterConverter> master_converter_;

    //! Map exctractor for dofs in this kinetics condition
    Teuchos::RCP<CORE::LINALG::MapExtractor> scatra_map_extractor_;

    //! map with the size of the graph of each matrix
    std::map<BlockMatrixType, int> size_matrix_graph_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //!
  class ScaTraManifoldScaTraFluxEvaluator
  {
   public:
    explicit ScaTraManifoldScaTraFluxEvaluator(const SSI::SSIMono& ssi_mono);

    //! call complete on system matrices and coupling matrices
    //@{
    void CompleteMatrixManifoldScaTra();
    void CompleteMatrixManifoldStructure();
    void CompleteMatrixScaTraManifold();
    void CompleteMatrixScaTraStructure();
    void CompleteSystemMatrixManifold();
    void CompleteSystemMatrixScaTra();
    //@}

    //! write inflow fluxes to csv file
    bool DoOutput() const { return do_output_; }

    //! Evaluate everything including coupling
    void Evaluate();

    //! Evaluate inflow into manifold field from coupling with scatra field
    void EvaluateScaTraManifoldInflow();

    //! get all RHS
    //@{
    Teuchos::RCP<Epetra_Vector> RHSManifold() { return rhs_manifold_; }
    Teuchos::RCP<Epetra_Vector> RHSScaTra() { return rhs_scatra_; }
    //@}

    //! get all matrices
    //@{
    Teuchos::RCP<CORE::LINALG::SparseOperator> SystemMatrixManifold()
    {
      return systemmatrix_manifold_;
    }
    Teuchos::RCP<CORE::LINALG::SparseOperator> SystemMatrixScaTra() { return systemmatrix_scatra_; }

    Teuchos::RCP<CORE::LINALG::SparseOperator> MatrixManifoldScatra()
    {
      return matrix_manifold_scatra_;
    }
    Teuchos::RCP<CORE::LINALG::SparseOperator> MatrixScaTraManifold()
    {
      return matrix_scatra_manifold_;
    }
    Teuchos::RCP<CORE::LINALG::SparseOperator> MatrixManifoldStructure()
    {
      return matrix_manifold_structure_;
    }
    Teuchos::RCP<CORE::LINALG::SparseOperator> MatrixScaTraStructure()
    {
      return matrix_scatra_structure_;
    }
    //@}

    //! write coupling fluxes (inflow into manifold) to csv file
    void Output();

    //! return all scatra-scatra manifold coupling operators
    std::vector<Teuchos::RCP<SSI::ManifoldScaTraCoupling>>& ScaTraManifoldCouplings()
    {
      return scatra_manifold_couplings_;
    }

   private:
    //! Add to global matrices and rhs
    void AddConditionContribution();

    //! Copy and scale (-1.0) to manifold side
    void CopyScaTraScaTraManifoldSide(
        Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling);

    //! Evaluate flux and linearizations on bulk side
    void EvaluateBulkSide(Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling);

    //! Evaluate integral on scatra manifold over 1.0
    void EvaluateScaTraManifoldDomainIntegral(
        Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling);

    //! Evaluate integral on scatra manifold over positive fluxes
    void EvaluateScaTraManifoldInflowIntegral(
        Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling);

    //! prepare evaluation of coupling condition: set elemental data
    void PreEvaluate(Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling);

    //! uncomplete all global matrices if any matrices holding the condition contributions have
    //! updated graphs (i.e. zeros become non-zeros or vice versa). In this case the graph of the
    //! global matrices needs to be updated as well to be able to add the local matrices to the
    //! global matrices
    void UnCompleteMatricesIfNecessary(
        Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling);

    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_scatra_;

    //! map extractor associated with all degrees of freedom inside scatra manifold field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_scatra_manifold_;

    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure_;

    //! write inflow fluxes to csv file
    const bool do_output_;

    //! integral of manifold domain
    std::map<int, double> domainintegral_;

    //! map of all scatra manifold dofs
    Teuchos::RCP<const Epetra_Map> full_map_manifold_;

    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const Epetra_Map> full_map_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! integrated flux for each scalar into manifold
    std::map<int, std::vector<double>> inflow_;

    //! coupling matrices
    //@{
    Teuchos::RCP<CORE::LINALG::SparseOperator> matrix_manifold_scatra_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix_manifold_scatra_cond_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> matrix_manifold_structure_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix_manifold_structure_cond_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> matrix_scatra_manifold_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix_scatra_manifold_cond_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> matrix_scatra_structure_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix_scatra_structure_cond_;
    //@}

    //! rhs for manifold and scatra
    //@{
    Teuchos::RCP<Epetra_Vector> rhs_manifold_;
    Teuchos::RCP<Epetra_Vector> rhs_manifold_cond_;
    Teuchos::RCP<Epetra_Vector> rhs_scatra_;
    Teuchos::RCP<Epetra_Vector> rhs_scatra_cond_;
    //@}

    // writes evaluated data to output
    std::optional<IO::RuntimeCsvWriter> runtime_csvwriter_;

    //! scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! scatra manifold problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_manifold_;

    //! all scatra-scatra manifold coupling operators
    std::vector<Teuchos::RCP<SSI::ManifoldScaTraCoupling>> scatra_manifold_couplings_;

    //! SSI structure mesh tying object containing coupling adapters, converters and maps
    Teuchos::RCP<SSI::UTILS::SSIMeshTying> ssi_structure_meshtying_;

    //! system matrices of scatra and manifold
    //@{
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix_manifold_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix_manifold_cond_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix_scatra_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix_scatra_cond_;
    //@}
  };

  //! data types for mesh tying handlers for sparse case: coupling adapter and MultiMapExtractor
  //! splitting dofs into interior, master, and slave
  using meshtying_handler_type = std::pair<Teuchos::RCP<CORE::ADAPTER::Coupling>,
      Teuchos::RCP<CORE::LINALG::MultiMapExtractor>>;
  //! data types for mesh tying handlers for block case: standard mesh tying handler and
  //! MultiMapExtractor splitting slave dofs into blocks
  using meshtying_block_handler_type = std::pair<std::vector<Teuchos::RCP<CORE::ADAPTER::Coupling>>,
      std::vector<Teuchos::RCP<const Epetra_Map>>>;

  //! Base class to handle mesh tying between manifold fields
  class ManifoldMeshTyingStrategyBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~ManifoldMeshTyingStrategyBase() = default;

    explicit ManifoldMeshTyingStrategyBase(Teuchos::RCP<DRT::Discretization> scatra_manifold_dis,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps, bool is_manifold_meshtying);

    //! apply mesh tying to right hand side
    void ApplyMeshTyingToManifoldRHS(Teuchos::RCP<Epetra_Vector> rhs_manifold);

    //! apply mesh tying to manifold system matrix
    virtual void ApplyMeshtyingToManifoldMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_matrix) = 0;

    //! apply mesh tying to scatra manifold - scatra matrix
    virtual void ApplyMeshtyingToManifoldScatraMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_scatra_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_scatra_matrix) = 0;

    //! apply mesh tying to scatra manifold - structure matrix
    virtual void ApplyMeshtyingToManifoldStructureMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_structure_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) = 0;

    //! apply mesh tying to scatra manifold - scatra matrix
    virtual void ApplyMeshtyingToScatraManifoldMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_scatra_manifold_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> scatra_manifold_matrix,
        bool do_uncomplete) = 0;

    //! coupling adpaters, maps, and extractors for mesh tying
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> SSIMeshTying() const { return ssi_meshtying_; }

   protected:
    //! should meshtying between manifold fields be applied?
    const bool is_manifold_meshtying_;

    //! all interior and master dofs
    Teuchos::RCP<const Epetra_Map> condensed_dof_map_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps_;

   private:
    //! coupling adpaters, maps, and extractors for mesh tying
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_meshtying_;
  };

  class ManifoldMeshTyingStrategySparse : public ManifoldMeshTyingStrategyBase
  {
   public:
    explicit ManifoldMeshTyingStrategySparse(Teuchos::RCP<DRT::Discretization> scatra_manifold_dis,
        Teuchos::RCP<UTILS::SSIMaps> ssi_maps, bool is_manifold_meshtying);

    void ApplyMeshtyingToManifoldMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_matrix) override;

    void ApplyMeshtyingToManifoldScatraMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_scatra_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_scatra_matrix) override;

    void ApplyMeshtyingToManifoldStructureMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_structure_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraManifoldMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_scatra_manifold_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> scatra_manifold_matrix,
        bool do_uncomplete) override;
  };

  class ManifoldMeshTyingStrategyBlock : public ManifoldMeshTyingStrategyBase
  {
   public:
    explicit ManifoldMeshTyingStrategyBlock(Teuchos::RCP<DRT::Discretization> scatra_manifold_dis,
        Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps, bool is_manifold_meshtying);

    void ApplyMeshtyingToManifoldMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_matrix) override;

    void ApplyMeshtyingToManifoldScatraMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_scatra_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_scatra_matrix) override;

    void ApplyMeshtyingToManifoldStructureMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_structure_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_structure_matrix,
        bool do_uncomplete) override;

    void ApplyMeshtyingToScatraManifoldMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_scatra_manifold_matrix,
        Teuchos::RCP<const CORE::LINALG::SparseOperator> scatra_manifold_matrix,
        bool do_uncomplete) override;

    //! handles all coupling adapters and multimap extractors with additional MultiMapExtractor
    //! splitting the slave dof map into blocks
    std::vector<meshtying_block_handler_type> MeshTyingBlockHandler() const
    {
      return meshtying_block_handler_;
    }

   private:
    //! intersects @p block_map with @p intersecting_map. Extracts also dofs from @p perm_map if the
    //! corresponding dof in @p intersecting map intersects wirh @p block_map. Returns both in a
    //! tuple.
    std::tuple<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map>>
    IntersectCouplingMapsBlockMap(Teuchos::RCP<const Epetra_Map> block_map,
        Teuchos::RCP<const Epetra_Map> intersecting_map,
        Teuchos::RCP<const Epetra_Map> permuted_map, const Epetra_Comm& comm);

    //! all interior and master dofs split into blocks
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> condensed_block_dof_map_;

    //! handles all coupling adapters and multimap extractors with additional multimap extractors
    //! that split the slave dof map into blocks
    std::vector<meshtying_block_handler_type> meshtying_block_handler_;
  };

  //! build specific mesh tying strategy
  Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> BuildManifoldMeshTyingStrategy(
      Teuchos::RCP<DRT::Discretization> scatra_manifold_dis, Teuchos::RCP<UTILS::SSIMaps> ssi_maps,
      bool is_manifold_meshtying, CORE::LINALG::MatrixType matrixtype_manifold);

}  // namespace SSI

FOUR_C_NAMESPACE_CLOSE

#endif
