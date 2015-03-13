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
  //AC-FSI specific output
  if ( DRT::INPUT::IntegralValue<INPAR::FLUID::WSSType>(DRT::Problem::Instance()->FluidDynamicParams() ,"WSS_TYPE") == INPAR::FLUID::wss_mean )
  {
    IO::DiscretizationReader reader(discret_,step);

    Teuchos::RCP<Epetra_Vector> SumWss =Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap(0),true));
    reader.ReadVector(SumWss,"wss_fluid");

    double SumDtWss = 0.0;
    SumDtWss = reader.ReadDouble("wss_time");

    SumWss->Scale(SumDtWss);
    StressManager()->RestartWss(SumWss,SumDtWss);

  }

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::Output()
{
  FluidImplicitTimeInt::Output();

  //AC-FSI specific output
  if ( DRT::INPUT::IntegralValue<INPAR::FLUID::WSSType>(DRT::Problem::Instance()->FluidDynamicParams() ,"WSS_TYPE") == INPAR::FLUID::wss_mean )
  {
    if ( uprestart_ > 0 and step_%uprestart_ == 0 ) //iff we write an restartable output
    {
      output_->WriteVector("wss_fluid",stressmanager_->GetPreCalcWallShearStresses(trueresidual_)); //we need WSS when restarting an AC FSI simulation
      output_->WriteDouble("wss_time",stressmanager_->GetSumDtWss()); //we need WSS when restarting an AC FSI simulation
    }
  }
  return;
}
