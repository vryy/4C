/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_two_phase_ost.cpp
\brief TimIntTwoPhaseOst

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_two_phase_ost.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntTwoPhaseOst::TimIntTwoPhaseOst(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      TimIntOneStepTheta(actdis,solver,params,output,alefluid),
      TimIntTwoPhase(actdis,solver,params,output,alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhaseOst::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  TimIntTwoPhase::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntTwoPhaseOst::~TimIntTwoPhaseOst()
{
  return;
}

