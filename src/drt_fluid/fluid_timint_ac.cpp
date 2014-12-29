/*!----------------------------------------------------------------------
\file fluid_timint_ac.cpp
\brief

<pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_ac.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     Thon 12/14 |
 *----------------------------------------------------------------------*/
FLD::TimIntAC::TimIntAC(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid)
{
  return;
}

/*----------------------------------------------------------------------*
| Destructor (public)                                        Thon 12/14 |
*-----------------------------------------------------------------------*/
FLD::TimIntAC::~TimIntAC()
{
  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);

//  reader.ReadVector(gridv_,"gridv");
  reader.ReadVector(trueresidual_,"residual");

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::Output()
{
  FluidImplicitTimeInt::Output();

  if ( uprestart_ > 0 and step_%uprestart_ == 0 ) //iff we write an restartable output
  {
    output_->WriteVector("residual", trueresidual_); //we need this to be able to calculate wss
  }
  return;
}
