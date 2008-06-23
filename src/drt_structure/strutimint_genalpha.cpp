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
#include <Teuchos_StringToIntMap.hpp>

/*----------------------------------------------------------------------*/
/* converts string to enum */
enum StruTimIntGenAlpha::MidAverageEnum StruTimIntGenAlpha::MidAvgStringToEnum
(
  const std::string instring
)
{
  if (instring == "Vague")
  {
    return midavg_vague;
  }
  else if (instring == "ImrLike")
  {    
    return midavg_imrlike;
  }
  else if (instring == "TrLike")
  {
    return midavg_trlike;
  }
  else
  {
    return midavg_vague;
  }
}

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
  midavg_(MidAvgStringToEnum(genalphaparams.get<string>("GENAVG"))),
  beta_(genalphaparams.get<double>("BETA")),
  gamma_(genalphaparams.get<double>("GAMMA")),
  alphaf_(genalphaparams.get<double>("ALPHA_F")),
  alpham_(genalphaparams.get<double>("ALPHA_M"))
{
  return;
}



/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
