/*----------------------------------------------------------------------*/
/*!
\file strutimint_ab2.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_AB2_H
#define STRUTIMINT_AB2_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_expl.H"

/*======================================================================*/
/*!
\brief Front-end for structural dynamics by structure time integration. 
\author bborn
\date 06/08
*/
class StruTimIntAB2 : public StruTimIntExpl
{

public:

  //! Constructor
  StruTimIntAB2
  (
    const Teuchos::ParameterList& ioparams,  //!< ioflags
    const Teuchos::ParameterList& sdynparams,  //!< input parameters
    const Teuchos::ParameterList& xparams,  //!< extra flags
    //const Teuchos::ParameterList& ab2params,  //!< AB2 flags
    DRT::Discretization& actis,  //!< current discretisation
    LINALG::Solver& solver,  //!< the solver
    IO::DiscretizationWriter& output  //!< the output
  );

  //! Destructor
  // ....

  //! @name Actions
  //@{

  //! Do time integration of single step
  void IntegrateStep();

  //! Update configuration after time step
  //!
  //! Thus the 'last' converged is lost and a reset of the time step 
  //! becomes impossible. We are ready and keen awating the next time step.
  void UpdateStep();

  //@}

  //! @name Attribute access functions
  //@{

  //! Return time integrator name
  enum NameEnum MethodName() { return name_ab2; }

  //! Provide number of steps, e.g. a single-step method returns 1,
  //! a m-multistep method returns m
  int MethodSteps() { return 2; }

  //! Give local order of accuracy of displacement part
  int MethodOrderOfAccuracyDis() { return 2; }

  //! Give local order of accuracy of velocity part
  int MethodOrderOfAccuracyVel() { return 2; }

  //@}

protected:

  //! @name Global convenience state vectors at \f$t_{n-1}\f$
  //@{
  Teuchos::RCP<Epetra_Vector> diso_;  //!< mid-displacements
                                      //!< \f$D_{n-1}\f$
  Teuchos::RCP<Epetra_Vector> velo_;  //!< mid-velocities
                                      //!< \f$V_{n-1}\f$
  Teuchos::RCP<Epetra_Vector> acco_;  //!< mid-accelerations
                                      //!< \f$A_{n-1}\f$
  //@}

  //! @name Global forces at \f$t_{n+1}\f$
  //@{
  Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force
                                       //!< \f$F_{int;n+1}\f$
  Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force
                                       //!< \f$F_{int;n+1}\f$
  Teuchos::RCP<Epetra_Vector> fviscn_;  //!< Rayleigh viscous forces
                                        //!< \f$C \cdot V_{n+1}\f$
  Teuchos::RCP<Epetra_Vector> frimpn_;  //!< time derivative of
                                        //!< linear momentum
                                        //!< (temporal rate of impulse)
                                        //!< \f$\dot{P}_{n+1} = M \cdot \dot{V}_{n+1}\f$
  //@}

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_AB2_H
#endif  // #ifdef CCADISCRET
