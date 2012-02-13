/*!----------------------------------------------------------------------
\file statmech_time.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STATMECH_TIME_H
#define STATMECH_TIME_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_structure/strugenalpha.H"
/*!
\brief time integration for structural problems with statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/

// forward declarations
class StatMechManager;

namespace CONTACT
{
  class Beam3cmanager;
}


class StatMechTime : public virtual StruGenAlpha
{
public:

  //! \brief Standard Constructor
  explicit StatMechTime(Teuchos::ParameterList& params,
                   DRT::Discretization& dis,
                   LINALG::Solver& solver,
                   IO::DiscretizationWriter& output);

  //! \brief Destructor
  virtual ~StatMechTime() {};

  //! \brief Do time integration
  void Integrate();

  //! \brief do consistent predictor step
  virtual void ConsistentPredictor(Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //!\brief do full Newton iteration
  virtual void FullNewton(Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //!\brief solve non-linear problem by means of PTC method
  virtual void PTC(Teuchos::RCP<Epetra_MultiVector> randomnumbers, int& istep, bool* redoptc=NULL);

  //!\brief initialize Newton for 2nd, 3rd, ... Uzawa iteration
  void InitializeNewtonUzawa(Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //!\brief writing output including restart information
  virtual void Output();

  //!\brief read restart information
  virtual void ReadRestart(int step);

protected:

  // don't want = operator and cctor
  StatMechTime operator = (const StatMechTime& old);
  StatMechTime(const StatMechTime& old);

  //! \brief Initialize StatMechManager and BeamContactManager
  void InitializeManagers(int& currstep, int& startstep, int& ndim, double& dt, bool& buildoctree);

  //! \brief check convergence status of Newton scheme for Augmented Lagrange method
  bool ExitUzawaUpon(int& maxuzawaiter);

  //! \brief switch time step size at given point in time
  void StatMechManagerTimeUpdate(double& time, double& dt);

  //! \brife call Update() method of statmechmanager_ according to options chosen
  void StatMechManagerUpdate(int& istep, int& ndim, double& dt, bool buildoctree);

  //! \brief evaluates beam contact according to options chosen in the input file
  void BeamContactNonlinearSolve(int& istep, int& ndim, Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //! \brief beam contact penalty approach
  void SolveContactPenalty(int& istep, int& ndim, Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //! \brief beam contact augmented lagrange approach
  void SolveContactAugLag(int& istep, int& ndim, Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //! \brief In case an Uzawa step has to be repeated, this method restores the preivous state
  void RestorePrevious(Teuchos::RCP<LINALG::SparseOperator> stiffprev, Epetra_Vector& fresmprev, Epetra_Vector& disiprev,
                       Epetra_Vector& dismprev, Epetra_Vector& disnprev, Epetra_Vector& velmprev, Teuchos::RCP<Epetra_MultiVector> randomnumbers);

  //! \brief call Output() method of statmechmanager_ according to options chosen
  void StatMechManagerOutput(int& istep, int& ndim, double& time, double& dt);

  //! \brief Reset relevant values, vectors, discretization state, etc. before repeating the time step (if unconverged)
  void RestoreOldConvergedState(bool& buildoctree);

  //! \brief Compute Statistical Forces and Damping Matrix for Brownian Dynamics Applications
  void PTCBrownianForcesAndDamping(double& dt, double& crotptc, double& ctransptc);

  //! \brief evaluate whether to repeat the time step or just the uzawa step with altered penalty parameter
  bool PTCRepeat(int& numiter, int& maxiter, bool fresmnormdivergent);

  //! \brief set PTC status to unconverged
  void ConvergenceStatusUpdate(bool converged = false, bool increasestepcount = true);

  //! \brief manager for administration of statistical mechanics stuff
  Teuchos::RCP<StatMechManager> statmechmanager_;

  //! \brief maximal number of random numbers to be generated in each time step per element
  int maxrandomnumbersperglobalelement_;

  //! \brief zero if in current time step the iterative solution of the nonlinear system of equations has not converged
  int   isconverged_;

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

  StatMechTime::beamcstrategy contactstrategy_;

}; // class StatMechTime

#endif  // #ifndef STATMECH_H
#endif  // #ifdef CCADISCRET
