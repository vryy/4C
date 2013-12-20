/*!----------------------------------------------------------------------
\file fluid_timint_ost.H
\brief

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUID_TIMINT_OST_H
#define FLUID_TIMINT_OST_H

#include "fluidimplicitintegration.H"



#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntOneStepTheta : public virtual FluidImplicitTimeInt
{

public:

  /// Standard Constructor
  TimIntOneStepTheta(      const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntOneStepTheta();

  /*!
  \brief Print information about current time step to screen (reimplementation)

  */
  void PrintTimeStepInfo();

  /*!
  \brief Set the part of the righthandside belonging to the last
         timestep for incompressible or low-Mach-number flow

     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

     Stationary/af-generalized-alpha:

                   mom: hist_ = 0.0
                  (con: hist_ = 0.0)

     One-step-Theta:

                   mom: hist_ = veln_  + dt*(1-Theta)*accn_
                  (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

     BDF2: for constant time step:

                   mom: hist_ = 4/3 veln_  - 1/3 velnm_
                  (con: hist_ = 4/3 densn_ - 1/3 densnm_)


  */
  void SetOldPartOfRighthandside();

  /*!
  \brief Set states in the time integration schemes: differs between GenAlpha and the others

  */
  void SetStateTimInt();

  /*!
  \brief Calculate time derivatives for
         stationary/one-step-theta/BDF2/af-generalized-alpha time integration
         for incompressible and low-Mach-number flow
  */
  void CalculateAcceleration(
      const Teuchos::RCP<const Epetra_Vector>    velnp,    ///< velocity at     n+1
      const Teuchos::RCP<const Epetra_Vector>    veln,     ///< velocity at     n
      const Teuchos::RCP<const Epetra_Vector>    velnm,    ///< velocity at     n-1
      const Teuchos::RCP<const Epetra_Vector>    accn ,    ///< acceleration at n-1
      const Teuchos::RCP<Epetra_Vector>          accnp     ///< acceleration at n+1
  );

  /*!
  \brief Set gamma to a value

  */
  void SetGamma(Teuchos::ParameterList& eleparams);

  /*!
  \brief Scale separation

  */
  void Sep_Multiply();

  /*!
  \brief Output of filtered velocity

  */
  void OutputofFilteredVel(Teuchos::RCP<Epetra_Vector> outvec,
                           Teuchos::RCP<Epetra_Vector> fsoutvec);

  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementTimeParameter();

  /*!
  \brief return scheme-specific time integration parameter

  */
  const double TimIntParam() const;

  /*!
  \brief filtered quantities for classical LES models

  */
  void ApplyScaleSeparationForLES();

  /*!
  \brief return scaling factor for the residual

  */
  double ResidualScaling() const
  {
    return 1.0/(theta_*dta_);
  }

protected:

  /// don't want = operator and cctor
  TimIntOneStepTheta operator = (const TimIntOneStepTheta& old);

  /// copy constructor
  TimIntOneStepTheta (const TimIntOneStepTheta& old);


private:


}; // class TimIntOneStepTheta

} // namespace FLD

#endif  // #ifndef FLD_TIMINT_OST_H
