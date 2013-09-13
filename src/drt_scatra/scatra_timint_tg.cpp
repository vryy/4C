/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_tg.cpp

\brief  explicit Taylor Galerkin time integration scheme (TG) and implicit Characteristic Galerkin scheme (ICG)
        just for the pure 1st order transport equation, not yet implemented for convections-diffusion equation

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_tg.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"



/*----------------------------------------------------------------------*
 |  Constructor (public)                                   schott 05/11 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntTaylorGalerkin::TimIntTaylorGalerkin(
    Teuchos::RCP<DRT::Discretization>      actdis,
    Teuchos::RCP<LINALG::Solver>           solver,
    Teuchos::RCP<Teuchos::ParameterList>   params,
    Teuchos::RCP<Teuchos::ParameterList>   extraparams,
    Teuchos::RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output)
{
  if(scatratype_ != INPAR::SCATRA::scatratype_levelset)
    dserror("Taylor Galerkin timeintegration scheme should be used only for the 1D level-set transport equation");


  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  //solution at time n-1, for level set problems
  // only used for 2-step 3rd order taylor galerkin method
  phinm_  = LINALG::CreateVector(*dofrowmap,true);

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                 schott 05/11 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntTaylorGalerkin::~TimIntTaylorGalerkin()
{
  return;
}

/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ExplicitPredictor()
{
  return;
}

/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::PredictThermPressure()
{
  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::SetTimeForNeumannEvaluation(
  Teuchos::ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::AddNeumannToResidual()
{
  if (neumanninflow_) ComputeNeumannInflowTG(sysmat_,residual_);

  return;
}


/*----------------------------------------------------------------------*
 | compute Neumann inflow terms                            schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ComputeNeumannInflowTG(
    RCP<LINALG::SparseOperator> matrix,
    RCP<Epetra_Vector>          rhs)
{
  // time measurement: evaluate condition 'Neumann inflow'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'TaylorGalerkinNeumannInflow'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("scatratype",scatratype_);
  condparams.set<int>("timealgo",timealgo_);
  condparams.set("incremental solver",incremental_);

  condparams.set("time-step length", dta_);

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(condparams,"convective velocity field",convel_);
  AddMultiVectorToParameterList(condparams,"velocity field",vel_);

  //provide displacement field in case of ALE
  condparams.set("isale",isale_);
  if (isale_) AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();

  // add element parameters according to time-integration scheme
  AddSpecificTimeIntegrationParameters(condparams);

  std::string condstring("TaylorGalerkinNeumannInflow");
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  return;
}



/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::AVM3Separation()
{
  dserror("not available for TaylorGalerkin time integration");
  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme     schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::AddSpecificTimeIntegrationParameters(
  Teuchos::ParameterList& params)
{
  params.set("using stationary formulation",false);
  params.set("using generalized-alpha time integration",false);
  params.set("total time",time_);
  params.set("alpha_F",1.0);


  discret_->SetState("phinp",phinp_);
  discret_->SetState("phin",phin_);
  discret_->SetState("phinm",phinm_);

  if(reinitswitch_) discret_->SetState("phistart", phistart_);

  params.set<int>("timealgo",timealgo_);

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ComputeThermPressure()
{
  dserror("not available for TaylorGalerkin time integration");
}


/*----------------------------------------------------------------------*
 | compute time derivative                                     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ComputeTimeDerivative()
{
  dserror("not available for TaylorGalerkin time integration");
  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative of thermodynamic pressure           vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ComputeThermPressureTimeDerivative()
{
  dserror("not available for TaylorGalerkin time integration");
  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                         schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::Update(const int num)
{
  // phinm is needed for 2-step Taylor Galerkin methods as well as for restart of level set problems
  phinm_ ->Update(1.0,*phin_,0.0);

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  return;
}

/*----------------------------------------------------------------------*
 | update level set after reinitialization                  schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::UpdateReinit()
{

  //phinm is needed for restart of level set problems
  phinm_ ->Update(1.0,*phin_,0.0);

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  return;
}

/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::UpdateThermPressure()
{
  dserror("not available for TaylorGalerkin time integration");
  return;
}


/*----------------------------------------------------------------------*
 | update density at n for ELCH natural convection         schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::UpdateDensityElch()
{
  dserror("not available for TaylorGalerkin time integration");
  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart               schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::OutputRestart()
{
  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phin", phin_);

  // phinm is needed for 2-step Taylor Galerkin methods and to reconstruct the interface
  output_->WriteVector("phinm", phinm_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                         schott 05/11 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed
  reader.ReadVector(phinp_, "phinp");
  reader.ReadVector(phin_,  "phin");

  // phinm is needed for restart of level set problems
  reader.ReadVector(phinm_,  "phinm");

  return;
}

/*--------------------------------------------------------------------------------------------*
 | Redistribute the scatra discretization and vectors according to nodegraph   wichmann 10/11 |
 *--------------------------------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  // let the base class do the basic redistribution and transfer of the base class members
  ScaTraTimIntImpl::Redistribute(nodegraph);

  // now do all the tg specfic steps
  const Epetra_Map* newdofrowmap = discret_->DofRowMap();
  Teuchos::RCP<Epetra_Vector> old;

  if (fsphinp_ != Teuchos::null)
  {
    old = fsphinp_;
    fsphinp_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *fsphinp_);
  }

  if (phinm_ != Teuchos::null)
  {
    old = phinm_;
    phinm_ = LINALG::CreateVector(*newdofrowmap,true);
    LINALG::Export(*old, *phinm_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step     schott 05/11 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::PrepareFirstTimeStep()
{
  // evaluate Dirichlet boundary conditions at time t=0
  // the values should match your initial field at the boundary!
  //ApplyDirichletBC(time_,phin_,phidtn_);
  ApplyDirichletBC(time_,phin_,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | update of time-dependent variables for electrode kinetics  gjb 11/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ElectrodeKineticsTimeUpdate(const bool init)
{
	dserror("do not call this function");
  return;
}


/*----------------------------------------------------------------------*
 | set old part of RHS for galvanostatic equation             gjb 04/10 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntTaylorGalerkin::ElectrodeKineticsSetOldPartOfRHS()
{
	dserror("do not call this function");
  return;
}



