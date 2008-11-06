/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_ost.cpp
\brief One-Step-Theta time integration scheme

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "scatra_timint_ost.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,output)
{
    // -------------------------------------------------------------------
    // get a vector layout from the discretization to construct matching
    // vectors and matrices
    //                 local <-> global dof numbering
    // -------------------------------------------------------------------
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // Vectors passed to the element
    // -----------------------------

    // temporal solution derivative at time n
    phidtn_       = LINALG::CreateVector(*dofrowmap,true);

    return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::~TimIntOneStepTheta()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  // hist_ = phin_ + dt*(1-Theta)*phidtn_
  hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ExplicitPredictor()
{
  phinp_->Update(dta_, *phidtn_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Update()
{
  // update time derivative of phi
  double fact1 = 1.0/(theta_*dta_);
  double fact2 = (-1.0/theta_) +1.0;

  // phidt(n) = (phi(n)-phi(n-1)) / (Theta*dt(n)) - (1/Theta -1)*phidt(n-1)
  phidtn_->Update( fact1,*phinp_,-fact1,*phin_ ,fact2);

  // we know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  ApplyDirichletBC(time_,Teuchos::null,phidtn_);

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::OutputRestart()
{
  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin", phin_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed for One-Step-Theta restart
  reader.ReadVector(phinp_,"phinp");
  reader.ReadVector(phin_, "phin");
  reader.ReadVector(phidtn_, "phidtn");

  return;
}

/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step        gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrepareFirstTimeStep()
{
  ApplyDirichletBC(time_, phin_,phidtn_);
  CalcInitialPhidt();
  return;
}

/*----------------------------------------------------------------------*
 | calculate initial time derivative of phi at t=t_0           gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::CalcInitialPhidt()
{
  // time measurement:
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + calc inital phidt");
  if (myrank_ == 0)
  cout<<"SCATRA: calculating initial time derivative of phi\n"<<endl;

  // are we really at step 0?
  dsassert(step_==0,"Step counter is not 0");

  // call elements to calculate matrix and right-hand-side
  {
    // zero out matrix entries
    sysmat_->Zero();

    // create the parameters for the discretization
    ParameterList eleparams;

    // action for elements
    eleparams.set("action","initialize_one_step_theta");
    // other parameters that are needed by the elements
    eleparams.set("total time",time_);
    eleparams.set("thsl",theta_*dta_);
    eleparams.set("problem type",prbtype_);
    eleparams.set("fs subgrid diffusivity",fssgd_);
    if (prbtype_=="elch")
    {
      // get ELCH-specific paramter F/RT (default value for the temperature is 298K)
      const double frt = 96485.3399/(8.314472 * params_->get<double>("TEMPERATURE",298.0));
      eleparams.set("frt",frt); // factor F/RT
    }

    //provide velocity field (export to column map necessary for parallel evaluation)
    //SetState cannot be used since this Multivector is nodebased and not dofbased
    const Epetra_Map* nodecolmap = discret_->NodeColMap();
    RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
    LINALG::Export(*convel_,*tmp);
    eleparams.set("velocity field",tmp);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phi0",phin_);
    discret_->SetState("dens0",densnp_);
    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,residual_);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();
  }

  // apply Dirichlet boundary conditions to system matrix
  LINALG::ApplyDirichlettoSystem(sysmat_,phidtn_,residual_,phidtn_,*(dbcmaps_->CondMap()));

  // solve for phidtn
  solver_->Solve(sysmat_->EpetraOperator(),phidtn_,residual_,true,true);

  // reset the matrix (and its graph!) since we solved
  // a very special problem here that has a different sparsity pattern_
  if (params_->get<int>("BLOCKPRECOND") )
    ; //how to reset a block matrix ??
  else
    SystemMatrix()->Reset();

  return;
}

#endif /* CCADISCRET */
