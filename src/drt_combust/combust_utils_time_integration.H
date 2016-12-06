/*!----------------------------------------------------------------------
\file combust_utils_time_integration.H
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o af-generalized-alpha time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#ifndef COMBUST_UTILS_TIME_INTEGRATION_H
#define COMBUST_UTILS_TIME_INTEGRATION_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H" // should go away

//forward declarations

namespace LINALG
{
  class MapExtractor;
}

namespace COMBUST
{

//! namespace for all time discretization routines of the fluid field
namespace UTILS
{
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
  void SetOldPartOfRighthandside(
      const Teuchos::RCP<Epetra_Vector>&   veln,
      const Teuchos::RCP<Epetra_Vector>&   velnm,
      const Teuchos::RCP<Epetra_Vector>&   accn,
      const INPAR::FLUID::TimeIntegrationScheme timealgo,
      const double                         dta,
      const double                         theta,
      Teuchos::RCP<Epetra_Vector>&         histmom
  );

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
      const INPAR::FLUID::TimeIntegrationScheme  timealgo, ///< time integration id
      const int                                  step,     ///< time step
      const double                               theta,    ///< factor for OST method
      const double                               dta,      ///< time step size (delta t) n
      const double                               dtp,      ///< previous time step size (delta t) n-1
      const Teuchos::RCP<Epetra_Vector>          accnp     ///< acceleration at n+1
  );



}
}

#endif  // #ifndef
