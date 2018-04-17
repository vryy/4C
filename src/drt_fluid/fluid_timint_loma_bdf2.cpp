/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_loma_bdf2.cpp
\brief TimIntLomaBDF2

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_loma_bdf2.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLomaBDF2::TimIntLomaBDF2(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      TimIntBDF2(actdis,solver,params,output,alefluid),
      TimIntLoma(actdis,solver,params,output,alefluid)
{

  std::cout << "\nWARNING: Loma has never been tested with BDF2 time integration!!\n" << std::endl;
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLomaBDF2::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntBDF2::Init();
  TimIntLoma::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntLomaBDF2::~TimIntLomaBDF2()
{
  return;
}

