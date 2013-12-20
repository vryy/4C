/*!----------------------------------------------------------------------
\file fluid_utils_time_integration.H
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o af-generalized-alpha time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUID_UTILS_TIME_INTEGRATION_H
#define FLUID_UTILS_TIME_INTEGRATION_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include "../drt_inpar/inpar_fluid.H"
#include "fluid_utils_mapextractor.H" // should go away

//forward declarations

namespace LINALG
{
  class MapExtractor;
}

namespace FLD
{

//! namespace for all time discretization routines of the fluid field
namespace UTILS
{

  /*!
  \brief do explicit predictor step to start nonlinear iteration from
         a better initial value
                        +-                                      -+
                        | /     dta \          dta  veln_-velnm_ |
   velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
                        | \     dtp /          dtp     dtp       |
                        +-                                      -+
  */
  void ExplicitPredictor(
      const std::string                               predictor,
      const Teuchos::RCP<const Epetra_Vector>    veln,
      const Teuchos::RCP<const Epetra_Vector>    velnm,
      const Teuchos::RCP<const Epetra_Vector>    accn,
      const LINALG::MapExtractor&                velpressplitter,
      const double                               timealgo_constant,
      const double                               dta,
      const double                               dtp,
      const Teuchos::RCP<Epetra_Vector>          velnp,
      const Epetra_Comm&                         comm
  );

}
}

#endif  // #ifndef TIME_INTEGRATION_SCHEME_H
