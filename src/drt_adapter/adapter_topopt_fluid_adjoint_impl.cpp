/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_adjointTimeInt_adjointTimeInt_impl.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifdef CCADISCRET

#include "adapter_topopt_fluid_adjoint_impl.H"
#include "../linalg/linalg_utils.H"
#include "../drt_opti/topopt_fluidAdjointImplTimeIntegration.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidAdjointImpl::FluidAdjointImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output)
  : dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  adjointTimeInt_ = rcp(new TOPOPT::ADJOINT::ImplicitTimeInt(dis, *solver, *params, *output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdjointImpl::Velnp()
{
  return adjointTimeInt_->Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdjointImpl::Velaf()
{
  return adjointTimeInt_->Velaf();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdjointImpl::Veln()
{
  return adjointTimeInt_->Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdjointImpl::Velnm()
{
  return adjointTimeInt_->Velnm();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidAdjointImpl::Discretization()
{
  return adjointTimeInt_->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::FluidAdjointImpl::GetDBCMapExtractor()
{
  return adjointTimeInt_->DirichMaps();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::TimeLoop()
{
   adjointTimeInt_->Integrate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
{
  if (stepinc!=Teuchos::null)
  {
    adjointTimeInt_->Evaluate(stepinc);
  }
  else
  {
    adjointTimeInt_->Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::Output()
{
  adjointTimeInt_->Output();
}


IO::DiscretizationWriter& ADAPTER::FluidAdjointImpl::DiscWriter()
{
  return adjointTimeInt_->DiscWriter();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::ReadRestart(int step)
{
  adjointTimeInt_->ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::SetRestart(const int step,
                const double time,
                Teuchos::RCP<const Epetra_Vector> readvelnp,
                Teuchos::RCP<const Epetra_Vector> readveln,
                Teuchos::RCP<const Epetra_Vector> readvelnm,
                Teuchos::RCP<const Epetra_Vector> readaccnp,
                Teuchos::RCP<const Epetra_Vector> readaccn)
{
  adjointTimeInt_->SetRestart(step,time,readvelnp,readveln,readvelnm,readaccnp,readaccn);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const INPAR::FLUID::TimeIntegrationScheme ADAPTER::FluidAdjointImpl::TimIntScheme() const
{
  return adjointTimeInt_->TimIntScheme();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidAdjointImpl::CreateFieldTest()
{
  dserror("not implemented until now");
  return Teuchos::null;
//  return Teuchos::rcp(new FLD::FluidResultTest(adjointTimeInt_));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::SetTopOptData(RCP<Epetra_Vector> porosity,RCP<TOPOPT::Optimizer> optimizer)
{
  adjointTimeInt_->SetTopOptData(porosity,optimizer);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno)
{
   adjointTimeInt_->SetInitialFlowField(initfield,startfuncno);
   return;
}

#endif  // #ifdef CCADISCRET
