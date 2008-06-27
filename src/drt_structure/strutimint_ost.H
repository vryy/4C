/*----------------------------------------------------------------------*/
/*!
\file strutimint_ost.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_OST_H
#define STRUTIMINT_OST_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_impl.H"


/*======================================================================*/
/*!
\brief One-step-theta time integration
\author bborn
\date 06/08
*/
class StruTimIntOneStepTheta : public StruTimIntImpl
{

public:

  //! Constructor
  StruTimIntOneStepTheta
  (
    const Teuchos::ParameterList& ioparams,  //!< ioflags
    const Teuchos::ParameterList& sdynparams,  //!< input parameters
    const Teuchos::ParameterList& xparams,  //!< extra flags
    const Teuchos::ParameterList& onestepthetaparams,  //!< onesteptheta flags
    DRT::Discretization& actis,  //!< current discretisation
    LINALG::Solver& solver,  //!< the solver
    IO::DiscretizationWriter& output  //!< the output
  );

  //! Destructor
  // ....

  //! @name Pure virtual methods which have to be implemented
  //@{

  //! Consistent predictor with constant displacements
  //! and consistent velocities and displacements
  void PredictConstDisConsistVelAcc();

  //! Create force residual #fres_ and ist stiffness #stiff_
  void EvaluateForceStiffResidual();

  //! Determine characteristic norm for displacements
  //! \author lw (originally)
  double CalcRefNormDisplacement();

  //! Determine characteristic norm for force
  //! \author lw (originally)
  double CalcRefNormForce();

  //! Update iteration
  void UpdateIteration();

  //! Update step
  void UpdateStep();

  //@}

  //! @name Generalised-alpha specific methods
  //@{

  //! Evaluate mid-state vectors by averaging end-point vectors
  void EvaluateMidState();

  //@}

protected:

  //! equal operator is NOT wanted
  StruTimIntOneStepTheta operator=(const StruTimIntOneStepTheta& old);

  //! copy constructor is NOT wanted
  StruTimIntOneStepTheta(const StruTimIntOneStepTheta& old);

  //! @name Key coefficients
  //@{
  double theta_;  //!< factor (0,1]
  //@}

  //! @name Global mid-state vectors
  //@{
  Teuchos::RCP<Epetra_Vector> dism_;  //!< mid-displacements
                                      //!< \f$D_m = D_{n+1-\alpha_f}\f$
  Teuchos::RCP<Epetra_Vector> velm_;  //!< mid-velocities
                                      //!< \f$V_m = V_{n+1-\alpha_f}\f$
  Teuchos::RCP<Epetra_Vector> accm_;  //!< mid-accelerations
                                      //!< \f$A_m = A_{n+1-\alpha_m}\f$
  //@}

  //! @name Global force vectors
  //! Residual \c fres_ exists already in base class
  //@{
  Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$
  Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force
  Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

  Teuchos::RCP<Epetra_Vector> fext_;  //!< external force at \f$t_n\f$
  Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force
  Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

  Teuchos::RCP<Epetra_Vector> finertm_;  //!< interial force
  Teuchos::RCP<Epetra_Vector> fviscm_;  //!< viscous force

  Teuchos::RCP<Epetra_Vector> frobin_;  //!< Robin force due to FSI
  //@}

};  // end class StruTimIntOneStepTheta

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_OST_H
#endif  // #ifdef CCADISCRET
