#ifdef CCADISCRET

#include "adapter_condif_genalpha.H"

// further includes for FluidBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;


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
