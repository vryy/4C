/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all LOMA algorithms

\level 2

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef LOMA_ALGORITHM_H
#define LOMA_ALGORITHM_H

#include <Epetra_Comm.h>
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "../linalg/linalg_mapextractor.H"

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
    Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn,
        const Teuchos::ParameterList& solverparams);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Init(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            scatradyn,  ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
        const std::string& disname = "scatra",       ///< name of scalar transport discretization
        const bool isale = false                     ///< ALE flag
    );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

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
    // void EvaluateLomaODBlockMatScaTra(Teuchos::RCP<LINALG::SparseMatrix> mat_sf);

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
    int itmax_;
    int itmaxpre_;
    int itmaxbs_;
    double ittol_;

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
    //! problem dynamic parameters
    Teuchos::ParameterList probdyn_;
  };

}  // namespace LOMA

#endif  // LOMA_ALGORITHM_H
