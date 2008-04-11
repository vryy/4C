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
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/drt_globalproblem.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>


using namespace std;
using namespace Teuchos;

//! namespace for all time discretization routines of the fluid field
namespace FLUIDTIMEINTEGRATION
{
  /*!
  \brief do explicit predictor step to start nonlinear iteration from
         a better initial value

  */
  void ExplicitPredictor(
          RCP<Epetra_Vector>&         veln,
          RCP<Epetra_Vector>&         velnm,
          RCP<Epetra_Vector>&         accn,
          const double                dta,
          const double                dtp,
          RCP<Epetra_Vector>&         velnp
          );

  /*!
  \brief Set the part of the righthandside belonging to the last
         timestep.
  */
  void SetOldPartOfRighthandside(
          RCP<Epetra_Vector>&         veln,
          RCP<Epetra_Vector>&         velnm,
          RCP<Epetra_Vector>&         accn,
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
          RCP<Epetra_Vector>&         velnp,
          RCP<Epetra_Vector>&         veln,
          RCP<Epetra_Vector>&         velnm,
          const FLUID_TIMEINTTYPE     timealgo,
          const int                   step,
          const double                theta,
          const double                dta,
          const double                dtp,
          RCP<Epetra_Vector>&         accn,
          RCP<Epetra_Vector>&         accnm
          );


  /*!
  \brief determine grid velocity

  */
  void UpdateGridv(
          RCP<Epetra_Vector>&         dispnp,
          RCP<Epetra_Vector>&         dispn,
          RCP<Epetra_Vector>&         dispnm,
          const int                   order,
          const int                   step,
          const double                theta,
          const double                dta,
          const double                dtp,
          RCP<Epetra_Vector>&         gridv
          );
 }


#endif  // #ifndef TIME_INTEGRATION_SCHEME_H
#endif  // #ifdef CCADISCRET
