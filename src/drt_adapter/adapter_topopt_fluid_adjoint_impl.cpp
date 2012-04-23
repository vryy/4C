/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_fluid_adjoint_impl.cpp

\brief adapter for element routines of fluid adjoint equations in topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "adapter_topopt_fluid_adjoint_impl.H"
#include "../drt_opti/topopt_fluidAdjointImplTimeIntegration.H"
#include "../drt_opti/topopt_fluidAdjointResulttest.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidAdjointImpl::FluidAdjointImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output)
  : dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  adjointTimeInt_ = Teuchos::rcp(new TOPOPT::ADJOINT::ImplicitTimeInt(dis, *solver, *params, *output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdjointImpl::Velnp()
{
  return adjointTimeInt_->Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidAdjointImpl::Veln()
{
  return adjointTimeInt_->Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidAdjointImpl::Discretization()
{
  return adjointTimeInt_->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::TimeLoop()
{
   adjointTimeInt_->Integrate();
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
  return Teuchos::rcp(new TOPOPT::ADJOINT::FluidAdjointResultTest(adjointTimeInt_));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::SetTopOptData(
    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > fluidvelocities,
    Teuchos::RCP<Epetra_Vector> porosity,
    Teuchos::RCP<TOPOPT::Optimizer> optimizer)
{
  adjointTimeInt_->SetTopOptData(fluidvelocities,porosity,optimizer);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAdjointImpl::SetInitialAdjointField(const INPAR::TOPOPT::InitialAdjointField initfield,const int startfuncno)
{
   adjointTimeInt_->SetInitialAdjointField(initfield,startfuncno);
   return;
}

