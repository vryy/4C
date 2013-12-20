/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_loma_ost.cpp
\brief TimIntLomaOst

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_loma_ost.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLomaOst::TimIntLomaOst(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      TimIntOneStepTheta(actdis,solver,params,output,alefluid),
      TimIntLoma(actdis,solver,params,output,alefluid)
{

  std::cout << "\nWARNING: Loma has never been tested with BDF2 time integration!!\n" << std::endl;
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntLomaOst::~TimIntLomaOst()
{
  return;
}

