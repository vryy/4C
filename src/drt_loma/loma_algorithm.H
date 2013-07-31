/*----------------------------------------------------------------------*/
/*!
\file loma_algorithm.H

\brief Basis of all LOMA algorithms

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef LOMA_ALGORITHM_H
#define LOMA_ALGORITHM_H

#include <Epetra_Comm.h>
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "../linalg/linalg_mapextractor.H"


/*!
\brief LOMA: namespace of the module for low-Mach-number flow problems in BACI

*/
namespace LOMA
{

  /// LOMA algorithm base
  /*!

    Base class of LOMA algorithms. Derives from FluidBaseAlgorithm
    and ScatraBaseAlgorithm.
    There can (and will) be different subclasses that implement
    different coupling schemes.

    \author vg
    \date 08/08
   */
  class Algorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm
  {
  public:

    /// constructor
    explicit Algorithm(
        const Epetra_Comm& comm,             ///< communicator
        const Teuchos::ParameterList& prbdyn,///< problem-specific parameter list
        const Teuchos::ParameterList& solverparams
        );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// LOMA time loop
    void TimeLoop();

    /// read restart for preceding turbulent inflow simulation
    void ReadInflowRestart(int restart);

  protected:

    /// do initial calculations
    void InitialCalculations();

    /// prepare time step
    void PrepareTimeStep();

    /// do (partitioned) outer iteration loop
    void OuterLoop();

    /// do monolithic iteration loop
    void MonoLoop();

    /// set fluid values required in scatra
    void SetFluidValuesInScaTra();

    /// set scatra values required in fluid
    void SetScaTraValuesInFluid();

    /// set up right-hand-side for monolithic low-Mach-number system
    void SetupMonoLomaMatrix();

    /// evaluate off-diagonal block with fluid weighting functions
    void EvaluateLomaODBlockMatFluid(Teuchos::RCP<LINALG::SparseMatrix> mat_fs);

    /// evaluate off-diagonal block with scatra weighting functions
    //void EvaluateLomaODBlockMatScaTra(Teuchos::RCP<LINALG::SparseMatrix> mat_sf);

    /// set up right-hand-side for monolithic low-Mach-number system
    void SetupMonoLomaRHS();

    /// solve monolithic low-Mach-number system
    void MonoLomaSystemSolve();

    /// update for next iteration step for monolithic low-Mach-number system
    void IterUpdate();

    /// convergence Check for present iteration step
    bool ConvergenceCheck(int itnum);

    /// update for next time step
    void TimeUpdate();

    /// write output
    void Output();

    /// flag for monolithic solver
    bool monolithic_;

    /// dof row map splitted in (field) blocks for monolithic solver
    LINALG::MultiMapExtractor lomablockdofrowmap_;

    /// combined Dirichlet boundary condition map for monolithic solver
    /// (unique map of all dofs with Dirichlet boundary conditions)
    Teuchos::RCP<Epetra_Map> lomadbcmap_;

    /// incremental vector for monolithic solver
    Teuchos::RCP<Epetra_Vector> lomaincrement_;

    /// rhs vector for monolithic solver
    Teuchos::RCP<Epetra_Vector> lomarhs_;

    /// vector of zeros for Dirichlet boundary conditions for monolithic solver
    Teuchos::RCP<Epetra_Vector> zeros_;

    /// block matrix for monolithic solver
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> lomasystemmatrix_;

    /// monolithic solver
    Teuchos::RCP<LINALG::Solver> lomasolver_;

    /// time-step length, maximum time and maximum number of steps
    double dt_;
    double maxtime_;
    int stepmax_;

    /// (preliminary) maximum number of iterations and tolerance for outer iteration
    int     itmax_;
    int     itmaxpre_;
    int     itmaxbs_;
    double  ittol_;

    /// flag for constant thermodynamic pressure
    std::string consthermpress_;

    /// flag for special flow
    std::string special_flow_;

    /// start of sampling period
    int samstart_;

    /// flag for turbulent inflow
    bool turbinflow_;
    /// number of inflow steps
    int numinflowsteps_;

  private:

  };

} // namespace LOMA

#endif // LOMA_ALGORITHM_H
