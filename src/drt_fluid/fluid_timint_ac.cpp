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
#include "../drt_lib/drt_globalproblem.H"
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
  const Teuchos::ParameterList& fs3idynac = DRT::Problem::Instance()->FS3IDynamicParams().sublist("AC");
  const bool restartfromfsi = DRT::INPUT::IntegralValue<int>(fs3idynac,"RESTART_FROM_PART_FSI"); //fs3idynac.get<int>("RESTART_FROM_PART_FSI");

  if (not restartfromfsi) //standard restart
  {
    IO::DiscretizationReader reader(discret_,step);

    reader.ReadVector(trueresidual_,"residual");
  }
  //else //nothing to do in this case since we start from a partitioned Fsi problem where the trueresidual has not been written!

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
    output_->WriteVector("residual", trueresidual_); //we need this to be able to calculate WSS when restarting
  }
  return;
}
