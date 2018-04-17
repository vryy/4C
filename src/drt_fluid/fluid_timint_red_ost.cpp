/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_red_ost.cpp
\brief TimIntRedModelsOst

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_red_ost.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntRedModelsOst::TimIntRedModelsOst(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      TimIntOneStepTheta(actdis,solver,params,output,alefluid),
      TimIntRedModels(actdis,solver,params,output,alefluid)
{

  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModelsOst::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  TimIntRedModels::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  read restart data                                   rasthofer 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModelsOst::ReadRestart(int step)
{
  // call of base classes
  TimIntOneStepTheta::ReadRestart(step);
  TimIntRedModels::ReadRestart(step);

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntRedModelsOst::~TimIntRedModelsOst()
{
  return;
}

