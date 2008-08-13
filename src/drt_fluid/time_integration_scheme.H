/*!----------------------------------------------------------------------
\file time_integration_scheme.H
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef TIME_INTEGRATION_SCHEME_H
#define TIME_INTEGRATION_SCHEME_H

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_geometry/vector_definitions.H"


namespace FLD
{

//! namespace for all time discretization routines of the fluid field
namespace TIMEINT_THETA_BDF2
{
  /*!
  \brief do explicit predictor step to start nonlinear iteration from
         a better initial value
                       +-                                      -+
                       | /     dta \          dta  veln_-velnm_ |
   velnp_ =veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
                       | \     dtp /          dtp     dtp       |
                       +-                                      -+
  */
  void ExplicitPredictor(
      const Teuchos::RCP<Epetra_Vector>&   veln,
      const Teuchos::RCP<Epetra_Vector>&   velnm,
      const Teuchos::RCP<Epetra_Vector>&   accn,
      const double                         dta,
      const double                         dtp,
      Teuchos::RCP<Epetra_Vector>&         velnp
  );

  /*!
  \brief Set the part of the righthandside belonging to the last
         timestep.

   set part of the rhs vector belonging to the old timestep

           One-step-Theta:

                   hist_ = veln_ + dta*(1-Theta)*accn_


           BDF2: for constant time step:

                   hist_ = 4/3 veln_ - 1/3 velnm_
  */
  void SetOldPartOfRighthandside(
      const Teuchos::RCP<Epetra_Vector>&   veln,
      const Teuchos::RCP<Epetra_Vector>&   velnm,
      const Teuchos::RCP<Epetra_Vector>&   accn,
      const FLUID_TIMEINTTYPE              timealgo,
      const double                         dta,
      const double                         theta,
      Teuchos::RCP<Epetra_Vector>&         hist
  );
  
  //! Calculate time derivatives for Stationary/OST/BDF2 time integration
  void CalculateAcceleration(
      const Teuchos::RCP<Epetra_Vector>&   velnp,    ///< velocity at     n+1
      const Teuchos::RCP<Epetra_Vector>&   veln,     ///< velocity at     n
      const Teuchos::RCP<Epetra_Vector>&   velnm,    ///< velocity at     n-1
      const Teuchos::RCP<Epetra_Vector>&   accn,     ///< acceleration at n
      const FLUID_TIMEINTTYPE              timealgo, ///< time integration id
      const int                            step,     ///< time step
      const double                         theta,    ///< factor for one-step-theta method
      const double                         dta,      ///< time step size (delta t)
      const double                         dtp,      
      Teuchos::RCP<Epetra_Vector>&         accnp     ///< acceleration at n+1
  );
  
}
}

#endif  // #ifndef TIME_INTEGRATION_SCHEME_H
#endif  // #ifdef CCADISCRET
