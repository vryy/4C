/*!----------------------------------------------------------------------
\file time_integration_scheme.H
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and ale displacement.

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

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_globalproblem.H"

#include "Epetra_Vector.h"

using namespace std;
using namespace Teuchos;

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
      const RCP<Epetra_Vector>&   veln,
      const RCP<Epetra_Vector>&   velnm,
      const RCP<Epetra_Vector>&   accn,
      const double                dta,
      const double                dtp,
      RCP<Epetra_Vector>&         velnp
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
      const RCP<Epetra_Vector>&   veln,
      const RCP<Epetra_Vector>&   velnm,
      const RCP<Epetra_Vector>&   accn,
      const FLUID_TIMEINTTYPE     timealgo,
      const double                dta,
      const double                theta,
      RCP<Epetra_Vector>&         hist
  );

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration.
  */
  void CalculateAcceleration(
      const RCP<Epetra_Vector>&   velnp,
      const RCP<Epetra_Vector>&   veln,
      const RCP<Epetra_Vector>&   velnm,
      const FLUID_TIMEINTTYPE     timealgo,
      const int                   step,
      const double                theta,
      const double                dta,
      const double                dtp,
      RCP<Epetra_Vector>&         accn,
      RCP<Epetra_Vector>&         accnm
  );
}
}

#endif  // #ifndef TIME_INTEGRATION_SCHEME_H
#endif  // #ifdef CCADISCRET
