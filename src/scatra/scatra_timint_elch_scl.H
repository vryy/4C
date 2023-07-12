/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for elch

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_ELCH_SCL_H
#define SCATRA_TIMINT_ELCH_SCL_H

#include "inpar_elch.H"

#include "scatra_timint_elch.H"
#include "io_runtime_csv_writer.H"
#include "coupling_adapter.H"

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace FLD
{
  class Meshtying;
}

namespace ADAPTER
{
  class Coupling;
  class ScaTraBaseAlgorithm;
}  // namespace ADAPTER

namespace SCATRA
{

  class ScaTraTimIntElchSCL : public virtual ScaTraTimIntElch
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScaTraTimIntElchSCL(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! create result test for micro field
    Teuchos::RCP<DRT::ResultTest> CreateMicroFieldTest();

    //! get time integration of micro problem
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> MicroScaTraField();

    void NonlinearSolve() override;

    void Output(int num) override;

    void PrepareTimeLoop() override;

    void PrepareTimeStep() override;

    void ReadRestartProblemSpecific(int step, IO::DiscretizationReader& reader) override;

    void Setup() override;

    void TestResults() override;

    void Update(int num) override;

   protected:
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params) override;

    void CalcInitialPotentialField() override;

    void CreateMeshtyingStrategy() override;

   private:
    //! assemble micro and macro and apply mesh tying between micro and macro model
    void AssembleAndApplyMeshTying();

    //! stop Netwon loop on convergence and print L2-Norm of increments and residuals
    bool BreakNewtonLoopAndPrintConvergence();

    //! copy solution from coupling nodes from macro discretization to micro discretization
    void CopySolutionToMicroField();

    //! redistribute micro discretization to minimize processor interfaces
    void RedistributeMicroDiscretization();

    //! scale micro problem with associated area of macro field
    void ScaleMicroProblem();

    //! setup coupling between micro and macro field
    void SetupCoupling();

    //! update increments in micro and macro field
    void UpdateIterMicroMacro();

    //! write coupled nodes and node coordinates to csv file
    //! \param glob_micro_macro_coupled_node_gids      coupled micro macro nodes of all procs
    //! \param glob_macro_slave_node_master_node_gids  macro slave/master nodes of all procs
    void WriteCouplingToCSV(const std::map<int, int>& glob_micro_macro_coupled_node_gids,
        const std::map<int, int>& glob_macro_slave_node_master_node_gids);

    //! the micro problem is split into sub discretiations. This map relates all nodes in the sub
    //! problem (key) to the coupled node of each sub problem (value)
    std::map<int, int> coupled_micro_nodes_;

    //! DBC maps for coupled problem
    Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps_elch_scl_;

    //! block map of coupled ELCH-SCL problem
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> full_block_map_elch_scl_;

    //! map of coupled ELCH-SCL problem
    Teuchos::RCP<const Epetra_Map> full_map_elch_scl_;

    //! increment of coupled ELCH-SCL problem
    Teuchos::RCP<Epetra_Vector> increment_elch_scl_;

    //! map extractor to get the coupled dofs from macro discretization (CondMap) out of all macro
    //! dofs
    Teuchos::RCP<CORE::LINALG::MapExtractor> macro_coupling_dofs_;

    //! coupling adapter between micro (slave) and macro discretization (master).
    Teuchos::RCP<CORE::ADAPTER::Coupling> macro_micro_coupling_adapter_;

    //! map extractor to get micro and macro dofs from global vector
    //! cond. map: micro, other map: macro
    Teuchos::RCP<CORE::LINALG::MapExtractor> macro_micro_dofs_;

    //! type of system matrix of coupled ELCH-SCL problem
    const CORE::LINALG::MatrixType matrixtype_elch_scl_;

    //! map extractor to get the coupled dofs from micro discretization (CondMap) out of all micro
    //! dofs
    Teuchos::RCP<CORE::LINALG::MapExtractor> micro_coupling_dofs_;

    //! time integrator for micro problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> micro_timint_;

    //! residual of coupled ELCH-SCL problem
    Teuchos::RCP<Epetra_Vector> residual_elch_scl_;

    //! solver for coupled ELCH-SCL problem
    Teuchos::RCP<CORE::LINALG::Solver> solver_elch_scl_;

    //! system matrix of coupled ELCH-SCL problem
    Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix_elch_scl_;
  };
}  // namespace SCATRA
#endif