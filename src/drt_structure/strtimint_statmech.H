/*!----------------------------------------------------------------------
\file strtimint_statmech.H

\brief time integration for structural problems with statistical mechanics

\level 3

\maintainer Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276

*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_STATMECH_H
#define STRTIMINT_STATMECH_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_ost.H"
#include <iostream>
#include <iomanip>
#include <sys/stat.h>

/*----------------------------------------------------------------------*/
/*forward declarations*/
namespace LINALG
{
  class Solver;
  class SparseOperator;
}

namespace STATMECH
{
  class StatMechManager;
  class StatMechManagerBilayer;
}

namespace CONTACT
{
  class Beam3cmanager;
}
/*----------------------------------------------------------------------*/


/* belongs to structural dynamics namespace */
namespace STR
{
  /*!
  \brief implicit time integration scheme for structural problems involving statistical mechanics

  \author mueller (mueller@lnm.mw.tum.de)
  \date 03/12
  */
  class TimIntStatMech : public TimIntOneStepTheta
  {
  public:
    //! \brief Standard Constructor
    explicit TimIntStatMech(const Teuchos::ParameterList&          params,
                            const Teuchos::ParameterList&          sdynparams,
                            const Teuchos::ParameterList&          xparams,
                            Teuchos::RCP<DRT::Discretization>      actdis,
                            Teuchos::RCP<LINALG::Solver>           solver,
                            Teuchos::RCP<LINALG::Solver>           contactsolver,
                            Teuchos::RCP<IO::DiscretizationWriter> output);

    //! \brief Destructor
    virtual ~TimIntStatMech() {};

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
    \return bool
    \date 08/16
    \author rauch  */
    virtual void Init(
        const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::RCP<LINALG::Solver> solver );

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

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const { return INPAR::STR::dyna_statmech;}

    //! \brief print BC type to screen
    void StatMechPrintBCType();

    //! \brief retrieve random numbers per element
    void RandomNumbersPerElement();

    //! \brief retrieve random numbers per element
    void RandomNumbersPerNode();

    //! \brief suppresses output to screen/lof file in some special cases
    void SuppressOutput();

    //! \brief initialize beam contact
    void InitializeBeamContact();

//    //! \brief Do time integration
//    int Integrate();

    //! \brief do predictor step
    void Predict();

    //! \brief predict displacements and velocities
    void PredictConstDisConsistVel();

    //! nonlinear solve
    virtual INPAR::STR::ConvergenceStatus Solve();

    //! error action in case of non-convergence of nonlinear solver
    bool PerformErrorAction();

    //! \brief evaluate residual
    void EvaluateForceStiffResidual(Teuchos::ParameterList& params);

    //! \brief Apply Dirichlet Boundary values
    void ApplyDirichletBC(const double                time,
                          Teuchos::RCP<Epetra_Vector> dis,
                          Teuchos::RCP<Epetra_Vector> vel);

    //! \brief Apply external forces
    void ApplyForceExternal(const double                          time,
                            const Teuchos::RCP<Epetra_Vector>     dis,
                            const Teuchos::RCP<Epetra_Vector>     disn,
                            const Teuchos::RCP<Epetra_Vector>     vel,
                            Teuchos::RCP<Epetra_Vector>&          fext,
                            Teuchos::RCP<LINALG::SparseOperator>& fextlin);

    //! \brief internal forces and stiffness matrix
    void ApplyForceStiffInternal(const double                         time,
                                 const double                         dt,
                                 const Teuchos::RCP<Epetra_Vector>    dis,  //!<  displacement state
                                 const Teuchos::RCP<Epetra_Vector>    disi,  //!<  residual displacements
                                 const Teuchos::RCP<Epetra_Vector>    vel,  //!<  velocity state
                                 Teuchos::RCP<Epetra_Vector>          fint,  //!<  internal force
                                 Teuchos::RCP<LINALG::SparseOperator> stiff,  //!<  stiffness matrix
                                 double                               t_eval=0.0);  //!<  evaluation time

    //! \brief evaluate Mid State for one-step-theta scheme
    void EvaluateMidState();

    //! \brief build residual
    void BuildResidual();

    //!\brief do full Newton iteration
    void NewtonFull();

    //!\brief solve non-linear problem by means of PTC method
    void PTC();

    //! \brief calculates ref norms of displacement and force
    void CalcRefNorms();

    //!\brief initialize Newton for 2nd, 3rd, ... Uzawa iteration
    void InitializeNewtonUzawa();

    //!\brief does what it says
    void UpdateAndOutput();

    //!\brief does what it says
    bool IsConverged(){ return isconverged_; }

    //! \brief preparation for the current time step relevant to Statistical Mechanics
    void StatMechPrepareStep();

    //! \brief evaluates beam contact according to options chosen in the input file
    void BeamContactNonlinearSolve();

    //! \brief Reset relevant values, vectors, discretization state, etc. before repeating the time step (if unconverged)
    void StatMechRestoreConvState();

    //! \brief Reset relevant values, vectors, discretization state, etc. before repeating the time step (if unconverged)
    void PrintTotalSurfAreaperIter();

    void PrintTotalVolperIter();

    void PrintTotalEdgeLenperIter();

    void GmshOutputEveryIter();

    //! \brief call Update() method of statmechmanager_ according to options chosen
    void StatMechUpdate(bool newrandomnumbers);

    //! \brief call Output() method of statmechmanager_ according to options chosen
    void StatMechOutput();

    //! Do Line search iteration
    //!
    //! This routines expects a prepared negative residual force #fres_
    //! and associated effective stiffness matrix #stiff_
    int NewtonLS();

    //! Solver call (line search)
    int LsSolveNewtonStep();

    //! Update structural RHS and stiffness matrix (line search)
    void LsUpdateStructuralRHSandStiff(bool& isexcept,double& merit_fct);

    //! Evaluate the specified merit function (line search)
    //! (for pure structural problems this routine is rather short.
    //! However, to incorporate contact problems more easily the
    //! evaluation of the merit function is performed here.)
    //! return 0 if successful
    int LsEvalMeritFct(double& merit_fct);

    //! Print information concerning the last line search step (line search)
    void LsPrintLineSearchIter(double* mf_value, int iter_ls, double step_red);

    //! Check the inner linesearch loop for convergence (line search)
    bool LsConverged(double* mf_value, double step_red);

    //!\brief writing output including restart information
    //virtual void Output();

//    //!\brief read restart information
//    virtual void ReadRestart(int step);

  private:

    // don't want = operator and cctor
    TimIntStatMech operator = (const TimIntStatMech& old);
    TimIntStatMech(const TimIntStatMech& old);

    //! \brief print updated statmech values to screen
    void StatMechPrintUpdate(const double& t_admin);

    //! \brief check convergence status of Newton scheme for Augmented Lagrange method
    bool BeamContactExitUzawaAt(int& maxuzawaiter);

    //! \brief beam contact penalty approach
    void BeamContactPenalty();

    //! \brief beam contact augmented lagrange approach
    void BeamContactAugLag();

    //! \brief Compute Statistical Forces and Damping Matrix for Brownian Dynamics Applications
    void PTCBrownianForcesAndDamping(double& dt, double& crotptc, double& ctransptc, double& sumptc);

    //! \brief Statistical Mechanics PTC Update
    void PTCStatMechUpdate(double& ctransptc, double& crotptc, double& nc, double& resinit, double& alphaptc);

    void UpdateIter(const int iter);

    void UpdateIterIncrementally();

    void UpdateIterIteratively();

    //! \brief set PTC status to unconverged
    void ConvergenceStatusUpdate(bool converged);

    //! \brief spatial dimension
    int ndim_;

    //! \brief maximal number of random numbers to be generated in each time step per element
    int maxrandomnumbersperglobalelement_;

    //! \brief maximal number of random numbers to be generated in each time step per node
    int maxrandomnumbersperglobalnode_;

    //! \brief MultiVector holding random numbers
    Teuchos::RCP<Epetra_MultiVector> randomnumbers_;

    //! \brief false if in current time step the iterative solution of the nonlinear system of equations has not converged
    bool   isconverged_;

    //! total number of Newton iterations in simulation
    int iterges_;

    //! \brief flag signaling a new octree search for contact detection (especially important if a time step has to be repeated)
    bool buildoctree_;

    //! \brief object for managing beam contact
    Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager_;

    //! \brief enum beam contact solution strategy
    enum beamcstrategy
    {
        contact_none,
        contact_penalty,
        contact_auglag,
        contact_lagmult
    };
  };
}

#endif  // #ifndef STRTIMINT_STATMECH_H
