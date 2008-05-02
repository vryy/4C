#ifdef CCADISCRET

#include "adapter_condif.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ConDifTimeIntegration::~ConDifTimeIntegration()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ConDifImplicit::ConDifImplicit(
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output
    ):
    condif_(dis, *solver, *params, *output),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::PrepareTimeStep()
{
  condif_.PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::Update()
{
  condif_.Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::Output()
{
  condif_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::SetVelocityField(int veltype, int velfuncno)
{
  condif_.SetVelocityField(veltype, velfuncno);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel)
{
  condif_.SetVelocityField(veltype, extvel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::ReadRestart(int step)
{
  condif_.ReadRestart(step);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::Integrate()
{
  condif_.Integrate();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifImplicit::Solve()
{
  condif_.Solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ConDifImplicit::CreateFieldTest()
{
  return Teuchos::rcp(new ConDifResultTest(condif_));
}

#endif  // #ifdef CCADISCRET
