/*----------------------------------------------------------------------*/
/*! \file
 \brief routines for coupling between 1D arterial network and 2D/3D
        scatra-algorithm

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_ARTERY_HPP
#define FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_ARTERY_HPP

#include "baci_config.hpp"

#include "baci_linalg_utils_sparse_algebra_assemble.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_print.hpp"
#include "baci_scatra_timint_meshtying_strategy_base.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace DRT
{
  class Discretization;
}
namespace ADAPTER
{
  class Coupling;
  class ArtNet;
}  // namespace ADAPTER
namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScaTraArtCouplBase;
}

namespace SCATRA
{
  class MeshtyingStrategyArtery : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyArtery(
        SCATRA::ScaTraTimIntImpl* scatratimint  //!< scalar transport time integrator
    );

    //! return global map of degrees of freedom
    const Epetra_Map& DofRowMap() const override;

    //! return global map of degrees of freedom
    Teuchos::RCP<const Epetra_Map> ArtScatraDofRowMap() const;

    //! evaluate mesh-tying
    //! \note  nothing is done here
    //!        actual coupling (meshtying) is evaluated in Solve
    //!        reason for that is that we need the system matrix of the continuous scatra
    //!        problem with DBCs applied which is performed directly before calling solve
    void EvaluateMeshtying() override{};

    //! init
    void InitMeshtying() override;

    bool SystemMatrixInitializationNeeded() const override { return false; }

    Teuchos::RCP<CORE::LINALG::SparseOperator> InitSystemMatrix() const override
    {
      dserror(
          "This meshtying strategy does not need to initialize the system matrix, but relies "
          "instead on the initialization of the field. If this changes, you also need to change "
          "'SystemMatrixInitializationNeeded()' to return true");
      // dummy return
      return Teuchos::null;
    }

    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> InterfaceMaps() const override
    {
      dserror("InterfaceMaps() is not implemented in MeshtyingStrategyArtery.");
      return Teuchos::null;
    }

    //! setup
    void SetupMeshtying() override;

    //! solver
    const CORE::LINALG::Solver& Solver() const override;

    //! init the convergence check
    void InitConvCheckStrategy() override;

    //! solve resulting linear system of equations
    void Solve(const Teuchos::RCP<CORE::LINALG::Solver>& solver,         //!< solver
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
        const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
        const Teuchos::RCP<Epetra_Vector>& phinp,  //!< state vector at time n+1
        const int iteration,                       //!< number of current Newton-Raphson iteration
        CORE::LINALG::SolverParams& solver_params) const override;

    void SetupSystem(
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& residual                      //!< residual vector
    ) const;

    //! init the convergence check
    void SetArteryScatraTimeIntegrator(Teuchos::RCP<SCATRA::ScaTraTimIntImpl> artscatratimint);

    //! set the artery time integrator
    void SetArteryTimeIntegrator(Teuchos::RCP<ADAPTER::ArtNet> arttimint);

    //! set the element pairs that are close as found by search algorithm
    void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs);

    //! prepare a time step
    void PrepareTimeStep() const;

    //! set the artery pressure
    void SetArteryPressure() const;

    //! apply mesh movement
    void ApplyMeshMovement();

    //! block systemmatrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> CombinedSystemMatrix()
    {
      return comb_systemmatrix_;
    }

    //! get the combined rhs
    Teuchos::RCP<Epetra_Vector> CombinedRHS() const { return rhs_; }

    //! get the combined increment
    Teuchos::RCP<Epetra_Vector> CombinedIncrement() const { return comb_increment_; }

    //! access to time integrator
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ArtScatraField() { return artscatratimint_; }

    //! check if initial fields match
    void CheckInitialFields() const;

    //! update increment of 1D discretization
    void UpdateArtScatraIter(Teuchos::RCP<const Epetra_Vector> combined_inc);

    /*!
     * extract single field vectors
     * @param[i] globalvec combined 1D-3D vector
     * @param[o] vec_cont  3D vector
     * @param[o] vec_art   1D vector
     */
    void ExtractSingleFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
        Teuchos::RCP<const Epetra_Vector>& vec_cont,
        Teuchos::RCP<const Epetra_Vector>& vec_art) const;

   private:
    //! initialize the linear solver
    void InitializeLinearSolver(const Teuchos::ParameterList& scatraparams);
    //! time integrators
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> artscatratimint_;
    Teuchos::RCP<ADAPTER::ArtNet> arttimint_;

    //! mesh tying object
    Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase> arttoscatracoupling_;

    //! the two discretizations
    Teuchos::RCP<DRT::Discretization> artscatradis_;
    Teuchos::RCP<DRT::Discretization> scatradis_;

    //! block systemmatrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> comb_systemmatrix_;

    //! combined rhs
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! combined rhs
    Teuchos::RCP<Epetra_Vector> comb_increment_;

  };  // class MeshtyingStrategyArtery

}  // namespace SCATRA

BACI_NAMESPACE_CLOSE

#endif
