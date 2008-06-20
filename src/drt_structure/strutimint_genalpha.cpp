/*----------------------------------------------------------------------*/
/*!
\file strudyn_direct.cpp
\brief Structural time integration

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_genalpha.H"

/*----------------------------------------------------------------------*/
/* constructor */
StruTimIntGenAlpha::StruTimIntGenAlpha
(
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& genalphaparams,
  DRT::Discretization& actis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
: StruTimIntImpl
  (
    sdynparams,
    actis,
    solver,
    output
  ),
  beta_(genalphaparams.get<double>("BETA")),
  gamma_(genalphaparams.get<double>("GAMMA")),
  alphaf_(genalphaparams.get<double>("ALPHA_F")),
  alpham_(genalphaparams.get<double>("ALPHA_M"))
{
  return;
}



/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
