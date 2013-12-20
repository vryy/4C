/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_red_stat.cpp
\brief TimIntRedAirwaysStat

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_red_stat.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntRedAirwaysStat::TimIntRedAirwaysStat(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      TimIntStationary(actdis,solver,params,output,alefluid),
      TimIntRedAirways(actdis,solver,params,output,alefluid)
{

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntRedAirwaysStat::~TimIntRedAirwaysStat()
{
  return;
}

