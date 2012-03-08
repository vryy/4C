/*!------------------------------------------------------------------------------------------------*
\file adapter_topopt_fluid_adjoint_impl.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifdef CCADISCRET


#include "../drt_opti/topopt_fluidAdjointImplTimeIntegration.H"
#include "../drt_lib/drt_dserror.H"

#include "adapter_topopt_fluid_adjoint_impl.H"


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
void ADAPTER::FluidAdjointImpl::SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno)
{
   adjointTimeInt_->SetInitialFlowField(initfield,startfuncno);
   return;
}

#endif  // #ifdef CCADISCRET
