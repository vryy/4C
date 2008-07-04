/*----------------------------------------------------------------------*/
/*!
\file strutimint.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_H
#define STRUTIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/uzawasolver.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_surfstress/drt_potential_manager.H"
#include "../drt_contact/drt_contact_manager.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include <NOX.H>
#include <NOX_Epetra.H>
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"

#include "strutimint_state.H"

/*======================================================================*/
/*!
\brief Front-end for structural dynamics by structure time integration. 
\author bborn
\date 06/08
*/
class StruTimInt
{

public:

  //! Type of stress
  enum StressEnum
  {
    stress_none,  //!< undetermined
    stress_cauchy,  //!< Cauchy stress
    stress_pk1,  //!< 1st Piola-Kirchhoff stress
    stress_pk2  //! 2nd Piola-Kirchhoff stress
  };

  //! Map stress input string to enum
  static enum StressEnum MapStressStringToEnum
  (
    const std::string name  //!< identifier
  );

  //! Type of strain
  enum StrainEnum
  {
    strain_none,  //!< undetermined
    strain_ea,  //!< Euler-Almansi strain
    strain_gl  //!< Green-Lagrange strain
  };

  //! Map strain input string to enum
  static enum StrainEnum MapStrainStringToEnum
  (
    const std::string name  //!< identifier
  );

  //! List of all available time integration schemes
  enum NameEnum
  {
    name_genalpha,  //!< generalised-alpha
    name_onesteptheta,  //!< one-step-theta
    name_ab2  //!< Adams-Bashforth 2nd order
  };

  //! Constructor
  StruTimInt
  (
    const Teuchos::ParameterList& ioparams,  //!< ioflags
    const Teuchos::ParameterList& sdynparams,  //!< input parameters
    const Teuchos::ParameterList& xparams,  //!< extra flags
    DRT::Discretization& actdis,  //!< current discretisation
    LINALG::Solver& solver,  //!< the solver
    IO::DiscretizationWriter& output  //!< the output
  );

  //! Destructor
  virtual ~StruTimInt() { ; }

  //! @name Actions
  //@{

  //! Do time integration of multiple steps
  void Integrate();

  //! Do time integration of single step
  virtual void IntegrateStep() = 0;

  //! Update configuration after time step
  //!
  //! Thus the 'last' converged is lost and a reset of the time step becomes
  //! impossible. We are ready and keen awating the next time step.
  virtual void UpdateStep() = 0;

  //! print summary after step
  virtual void PrintStep() = 0;

  //! Output to file
  virtual void Output() = 0;

  //@}

  //! @name Attributes
  //@{

  //! Provide Name
  virtual enum NameEnum Name() = 0;

  //! Return true, if time integrator is implicit
  virtual bool Implicit() = 0;

  //! Return true, if time integrator is explicit
  bool Explicit() { return (not Implicit()); }

  //! Provide number of steps, e.g. a single-step method returns 1,
  //! a m-multistep method returns m
  virtual int MethodSteps() = 0;

  //! Access Dicretisation
  virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

  //@}

  //! @name Phase state
  //@{

  //! Displacements
  virtual Teuchos::RCP<Epetra_Vector> Disp() = 0;

  //! Velocities
  virtual Teuchos::RCP<Epetra_Vector> Vel() = 0;

  //!  Accelerations
  virtual Teuchos::RCP<Epetra_Vector> Acc() = 0;

  //@}

protected:

  //! @name General purpose algorithm parameters
  //@{
  DRT::Discretization& discret_;  //!< attached discretisation
  int myrank_;  //!< ID of actual processor in parallel
  const Epetra_Map* dofrowmap_;  //!< vector layout from the discretization 
                                 //!< to construct matching
                                 //!< vectors and matrices
  LINALG::Solver& solver_;  //!< linear algebraic solver
  bool solveradapttol_;  //!< adapt solver tolerance
  double solveradaptolbetter_;  //!< tolerance to which is adapated ????
  //@}

  //! @name Printing and output
  //@{
  IO::DiscretizationWriter& output_;  //!< binary output
  bool printscreen_;  //!< print infos to standard out
  FILE* errfile_;  //!< error file handle
  bool printerrfile_;  //!< print infos to error file #errfile_
  bool printiter_;  //!< print intermediate iterations during solution
  int writerestartevery_;  //!< write restart every given step;
                           //!< if 0, restart is not written
  bool writestate_;  //!< write state on/off
  int writestateevery_;  //!< write state every given step
  int writestrevery_;  //!< write stress/strain every given step
  enum StressEnum writestress_;  //!< stress output type
  enum StrainEnum writestrain_;  //!< strain output type
  //@}

  //! @name Rayleigh damping
  //!       \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
  //@{
  bool damping_;  //!< damping active/inactive
  double dampk_;  //!< damping factor for stiffness \f$c_\text{K}\f$
  double dampm_;  //!< damping factor for mass \f$c_\text{M}\f$
  //@}

  //! @name General control parameters
  //@{
  double time_;  //!< time \f$t_{n}\f$ of last converged step
  double timen_;  //!< target time \f$t_{n+1}\f$
  double dt_;  //!< time step size \f$\Delta t\f$
  double timemax_;  //!< final time \f$t_\text{fin}\f$
  int stepmax_;  //!< final step \f$N\f$
  int step_;  //!< time step index \f$n\f$
  int stepn_;  //!< time step index \f$n+1\f$
  //@}

  //! @name Global vectors
  //@{
  Teuchos::RCP<Epetra_Vector> dirichtoggle_;  //!< At a free DOF this vector
                                     //!< holds a 0; at a supported DOF
                                     //!< it is a 1.
  Teuchos::RCP<Epetra_Vector> invtoggle_;  //!< Opposite of 
                                           //!< #dirichtoggle_ vector:
                                           //!< i.e.\ for each component 
                                           //!< a 0 indicates a supported 
                                           //!< corresponding DOF; a 1 indicates
                                           //!< a free/unsupported DOF.
  Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
  //@}


  //! @name Phase states
  //@{
  Teuchos::RCP<StruTimIntState> state_;  //!< converged states 
                                         //!< at t_{n}, t_{n-1}, ...
  Teuchos::RCP<StruTimIntState> staten_;  //!< current state at t_{n+1}
  //@}

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_H
#endif  // #ifdef CCADISCRET
