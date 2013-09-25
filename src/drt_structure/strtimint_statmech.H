/*!----------------------------------------------------------------------
\file statmech_time.H

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
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

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const { return INPAR::STR::dyna_statmech;}

    //! \brief print BC type to screen
    void StatMechPrintBCType();

    //! \brief retrieve random numbers per element
    void RandomNumbersPerElement();

    //! \brief suppresses output to screen/lof file in some special cases
    void SuppressOutput();

    //! \brief initialize beam contact
    void InitializeBeamContact();

    //! \brief Do time integration
    int Integrate();

    //! \brief do predictor step
    void Predict();

    //! \brief predict displacements and velocities
    void PredictConstDisConsistVel();

    //! \brief evaluate residual
    void EvaluateForceStiffResidual(bool predict=false);

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

    //!\brief writing output including restart information
    //virtual void Output();

//    //!\brief read restart information
//    virtual void ReadRestart(int step);

  private:

    // don't want = operator and cctor
    TimIntStatMech operator = (const TimIntStatMech& old);
    TimIntStatMech(const TimIntStatMech& old);

    //! \brief preparation for the current time step relevant to Statistical Mechanics
    void StatMechPrepareStep();

    //! \brief call Update() method of statmechmanager_ according to options chosen
    void StatMechUpdate();

    //! \brief print updated statmech values to screen
    void StatMechPrintUpdate(const double& t_admin);

    //! \brief call Output() method of statmechmanager_ according to options chosen
    void StatMechOutput();

    //! \brief Reset relevant values, vectors, discretization state, etc. before repeating the time step (if unconverged)
    void StatMechRestoreConvState();

    //! \brief preparation for the current time step relevant to Beam Contact
    void BeamContactPrepareStep();

    //! \brief check convergence status of Newton scheme for Augmented Lagrange method
    bool BeamContactExitUzawaAt(int& maxuzawaiter);

    //! \brief evaluates beam contact according to options chosen in the input file
    void BeamContactNonlinearSolve();

    //! \brief beam contact penalty approach
    void BeamContactPenalty();

    //! \brief beam contact augmented lagrange approach
    void BeamContactAugLag();

    //! \brief Compute Statistical Forces and Damping Matrix for Brownian Dynamics Applications
    void PTCBrownianForcesAndDamping(double& dt, double& crotptc, double& ctransptc, double& sumptc);

    //! \brief Statistical Mechanics PTC Update
    void PTCStatMechUpdate(double& ctransptc, double& crotptc, double& nc, double& resinit, double& alphaptc);

    //! \brief evaluate whether to repeat the time step or just the uzawa step with altered penalty parameter
    void PTCConvergenceStatus(int& numiter, int& maxiter, bool fresmnormdivergent);

    void UpdateIter(const int iter);

    void UpdateIterIncrementally();

    void UpdateIterIteratively();

    //! \brief set PTC status to unconverged
    void ConvergenceStatusUpdate(bool converged = false, bool increasestepcount = true);

    //! \brief spatial dimension
    int ndim_;

    //! \brief maximal number of random numbers to be generated in each time step per element
    int maxrandomnumbersperglobalelement_;

    //! \brief MultiVector holding random numbers
    Teuchos::RCP<Epetra_MultiVector> randomnumbers_;

    //! \brief old style inverted dirichtoggle vector
    Teuchos::RCP<Epetra_Vector> invtoggle_;

    //! \brief false if in current time step the iterative solution of the nonlinear system of equations has not converged
    bool   isconverged_;

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
