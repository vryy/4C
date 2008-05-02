#ifdef CCADISCRET

#include "adapter_condif_genalpha.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ConDifGenAlpha::ConDifGenAlpha(
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
void ADAPTER::ConDifGenAlpha::PrepareTimeStep()
{
  dserror("PrepareTimeStep() for GenAlpha not implemented.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::Update()
{
  condif_.GenAlphaTimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::Output()
{
  condif_.GenAlphaOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::SetVelocityField(int veltype, int velfuncno)
{
  condif_.SetVelocityField(veltype, velfuncno);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel)
{
  condif_.SetVelocityField(veltype, extvel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::ReadRestart(int step)
{
  condif_.ReadRestart(step);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::Integrate()
{
  condif_.TimeLoop();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ConDifGenAlpha::Solve()
{
  dserror("Gen. Alpha of convection-diffusion not yet supported");
}

#endif  // #ifdef CCADISCRET
