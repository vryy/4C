/*----------------------------------------------------------------------*/
/*! \file
 \brief Evaluates flux between ScaTra and ScaTra on manifolds incl. coupling matrices

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_MANIFOLD_FLUX_EVALUATOR_H_
#define SSI_MANIFOLD_FLUX_EVALUATOR_H_

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

namespace ADAPTER
{
  class Coupling;
  class CouplingSlaveConverter;
  class ScaTraBaseAlgorithm;
}  // namespace ADAPTER

namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace LINALG
{
  class MapExtractor;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace SSI
{
  class SSIMono;

  //! holds everything to evaluate coupling between scatra and scatra manifold
  class ManifoldScaTraCoupling
  {
   public:
    ManifoldScaTraCoupling(Teuchos::RCP<DRT::Discretization> manifolddis,
        Teuchos::RCP<DRT::Discretization> scatradis, DRT::Condition* condition_manifold,
        DRT::Condition* condition_kinetics);

    //! Kinetics condition on scatra dis
    DRT::Condition* ConditionKinetics() const { return condition_kinetics_; }

    //! manifold condition
    DRT::Condition* ConditionManifold() const { return condition_manifold_; }

    //! coupling adapter between manifold (master) and scatra (slave)
    Teuchos::RCP<ADAPTER::Coupling> CouplingAdapter() const { return coupling_adapter_; }

    //! condition ID of manifold condition
    int ManifoldConditionID() const { return manifold_conditionID_; }

    //! Map exctractor for dofs in this manifold condition
    Teuchos::RCP<LINALG::MapExtractor> ManifoldMapExtractor() const
    {
      return manifold_map_extractor_;
    }

    //! Map exctractor for dofs in this kinetics condition
    Teuchos::RCP<LINALG::MapExtractor> ScaTraMapExtractor() const { return scatra_map_extractor_; }

    //! Slave converter for scatra - manifold coupling
    Teuchos::RCP<ADAPTER::CouplingSlaveConverter> SlaveConverter() const
    {
      return slave_converter_;
    }

   private:
    //! Kinetics condition on scatra dis
    DRT::Condition* condition_kinetics_;

    //! manifold condition
    DRT::Condition* condition_manifold_;

    //! coupling adapter between manifold (master) and scatra (slave)
    Teuchos::RCP<ADAPTER::Coupling> coupling_adapter_;

    //! condition ID of manifold condition
    const int manifold_conditionID_;

    //! Map exctractor for dofs in this manifold condition
    Teuchos::RCP<LINALG::MapExtractor> manifold_map_extractor_;

    //! Map exctractor for dofs in this kinetics condition
    Teuchos::RCP<LINALG::MapExtractor> scatra_map_extractor_;

    //! Slave converter for scatra - manifold coupling
    Teuchos::RCP<ADAPTER::CouplingSlaveConverter> slave_converter_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //!
  class ScaTraManifoldScaTraFluxEvaluator
  {
   public:
    ScaTraManifoldScaTraFluxEvaluator(const SSI::SSIMono& ssi_mono,
        Teuchos::RCP<const Epetra_Map> interface_map_scatra,
        Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure);

    //! Evaluate everything including couling
    void Evaluate();

    //! get all RHS
    //@{
    Teuchos::RCP<Epetra_Vector> RHSManifold() { return rhs_manifold_; }
    Teuchos::RCP<Epetra_Vector> RHSScaTra() { return rhs_scatra_; }
    //@}

    //! get all matrices
    //@{
    Teuchos::RCP<LINALG::SparseOperator> SystemMatrixManifold() { return systemmatrix_manifold_; }
    Teuchos::RCP<LINALG::SparseOperator> SystemMatrixScaTra() { return systemmatrix_scatra_; }

    Teuchos::RCP<LINALG::SparseOperator> MatrixManifoldScatra() { return matrix_manifold_scatra_; }
    Teuchos::RCP<LINALG::SparseOperator> MatrixScaTraManifold() { return matrix_scatra_manifold_; }
    Teuchos::RCP<LINALG::SparseOperator> MatrixManifoldStructure()
    {
      return matrix_manifold_structure_;
    }
    Teuchos::RCP<LINALG::SparseOperator> MatrixScaTraStructure()
    {
      return matrix_scatra_structure_;
    }
    //@}

   private:
    //! Evaluate flux and linearizations on manifold side
    void EvaluateManifoldSide(Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling,
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_manifold_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_manifold_scatra_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_manifold_structure_cond,
        Teuchos::RCP<Epetra_Vector> rhs_manifold_cond);

    //! Copy and scale (-1.0) to scatra side
    void CopyScaTraManifoldScaTraMasterSide(
        Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling,
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_manifold_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_manifold_scatra_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_manifold_structure_cond,
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_scatra_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_scatra_manifold_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_scatra_structure_cond,
        Teuchos::RCP<Epetra_Vector> rhs_manifold_cond, Teuchos::RCP<Epetra_Vector> rhs_scatra_cond);

    //! Add to global matrices and rhs
    void AddConditionContribution(Teuchos::RCP<LINALG::SparseMatrix> symstemmatrix_manifold_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_manifold_scatra_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_manifold_structure_cond,
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_scatra_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_scatra_manifold_cond,
        Teuchos::RCP<LINALG::SparseMatrix> matrix_scatra_structure_cond,
        Teuchos::RCP<Epetra_Vector> rhs_manifold_cond, Teuchos::RCP<Epetra_Vector> rhs_scatra_cond);

    //! map extractor associated with all degrees of freedom on scatra interfaces
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_scatra_interface_;

    //! map extractor associated with all degrees of freedom inside structure field
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure_;

    //! map of all scatra manifold dofs
    Teuchos::RCP<const Epetra_Map> full_map_manifold_;

    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const Epetra_Map> full_map_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Epetra_Map> full_map_structure_;

    //! scatra-scatra interface coupling adapter for structure field
    Teuchos::RCP<const ADAPTER::Coupling> icoup_structure_;

    //! map associated with all degrees of freedom on scatra interface
    Teuchos::RCP<const Epetra_Map> interface_map_scatra_;

    //! coupling matrices
    //@{
    Teuchos::RCP<LINALG::SparseOperator> matrix_manifold_scatra_;
    Teuchos::RCP<LINALG::SparseOperator> matrix_manifold_structure_;
    Teuchos::RCP<LINALG::SparseOperator> matrix_scatra_manifold_;
    Teuchos::RCP<LINALG::SparseOperator> matrix_scatra_structure_;
    //@}

    //! rhs for manifold and scatra
    //@{
    Teuchos::RCP<Epetra_Vector> rhs_manifold_;
    Teuchos::RCP<Epetra_Vector> rhs_scatra_;
    //@}

    //! scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! scatra manifold problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_manifold_;

    //! all scatra-scatra manifold coupling operators
    std::vector<Teuchos::RCP<SSI::ManifoldScaTraCoupling>> scatra_manifold_couplings_;

    //! system matrices of scatra and manifold
    //@{
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix_manifold_;
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix_scatra_;
    //@}
  };
}  // namespace SSI

#endif /* SSI_MANIFOLD_FLUX_EVALUATOR_H_ */
