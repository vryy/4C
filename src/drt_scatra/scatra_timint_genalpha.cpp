/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_genalpha.cpp
\brief Generalized-alpha time-integration scheme

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "scatra_timint_genalpha.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       vg 11/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntGenAlpha::TimIntGenAlpha(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,output),
  alphaM_(params_->get<double>("alpha_M")),
  alphaF_(params_->get<double>("alpha_F")),
  gamma_ (params_->get<double>("gamma"))
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // Vectors passed to the element
  // -----------------------------

  // scalar at time n+alpha_F
  phiaf_ = LINALG::CreateVector(*dofrowmap,true);

  // temporal solution derivative at times n+1 and n
  phidtnp_ = LINALG::CreateVector(*dofrowmap,true);
  phidtn_  = LINALG::CreateVector(*dofrowmap,true);

  // only required for low-Mach-number flow
  if (prbtype_ == "loma")
  {
   // density at times n, n+alpha_M and n+alpha_F
    densn_  = LINALG::CreateVector(*dofrowmap,true);
    densam_ = LINALG::CreateVector(*dofrowmap,true);
    densaf_ = LINALG::CreateVector(*dofrowmap,true);

    // time derivative of density at times n+1 and n
    densdtnp_ = LINALG::CreateVector(*dofrowmap,true);
    densdtn_  = LINALG::CreateVector(*dofrowmap,true);

    // time derivative of thermodynamic pressure at n+alpha_F and n
    // (computed if not constant, otherwise remaining zero)
    thermpressdtaf_ = 0.0;
    thermpressdtn_  = 0.0;
  }

  // compute specific time factor for generalized-alpha time integration:
  // genalphatimefac = gamma*alpha_F/alpha_M
  if (alphaM_ < EPS12) dserror("factor alpha_M lower than or equal zero");
  genalphafac_ = gamma_/alphaM_;

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     vg 11/08 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntGenAlpha::~TimIntGenAlpha()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                             vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::SetOldPartOfRighthandside()
{
  // calculation of history vector only for non-incremental formulation:
  // (History vector is used in bot cases, but in incremental case, it
  // contains time derivatives of (density times) scalar, see below.)
  if (not incremental_)
  {
    // For conservative formulation of low-Mach-number flow:
    // hist_ = densn_*phin_ + dt*(1-(gamma/alpha_M))*densn_*phidtn_
    //                      + dt*(1-(gamma/alpha_M))*phin_*densdtn_
    //       = densn_*phin_ + dt*(1-genalphafac)*densn_*phidtn_
    //                      + dt*(1-genalphafac)*phin_*densdtn_
    if (prbtype_ == "loma" and convform_ =="conservative")
    {
      hist_->Multiply(1.0, *phin_, *densn_, 0.0);
      hist_->Update(dta_*(1.0-genalphafac_),*phidtn_,1.0);
    }
    // hist_ = phin_ + dt*(1-(gamma/alpha_M))*phidtn_
    //       = phin_ + dt*(1-genalphafac)*phidtn_
    else hist_->Update(1.0, *phin_, dta_*(1.0-genalphafac_), *phidtn_, 0.0);
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute initial time derivative of density field            vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeInitialDensityDerivative()
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // define auxiliary vectors
  Teuchos::RCP<Epetra_Vector> invphi = LINALG::CreateVector(*dofrowmap,true);
  Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*dofrowmap,true);

  // densdtn_ = densn_*thermpressdtn_/thermpressn_ - densn_*phidtn_/phin_
  invphi->Reciprocal(*phin_);
  tmp->Multiply(1.0, *invphi, *densn_, 0.0);
  densdtn_->Multiply(-1.0, *tmp, *phidtn_, 0.0);
  densdtn_->Update((thermpressdtn_/thermpressn_), *densn_, 1.0);

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                          vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ExplicitPredictor()
{
  // constant predictor
  phinp_->Update(1.0,*phin_,0.0);
  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)
  // prediction of time derivative:
  double fact = (gamma_-1.0)/gamma_;
  thermpressdtnp_ = fact*thermpressdtn_;

  // same-thermodynamic-pressure-derivative predictor (currrently not used)
  //thermpressnp_ += dta_*thermpressdtn_;
  // prediction of time derivative not required (would also not be required
  // to be performed, since we just updated the time derivatives of density,
  // and thus, thermpressdtnp_ = thermpressdtn_)

  // time derivative of thermodynamic pressure at n+alpha_F
  // -> required as right-hand-side contribution to temperature equation
  thermpressdtaf_ = alphaF_*thermpressdtnp_ + (1.0-alphaF_)*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | predict density for next time step for low-Mach-number flow vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PredictDensity()
{
  // same-density predictor (not required to be performed, since we just
  // updated the density field, and thus, densnp_ = densn_)
  // consistent prediction of time derivative
  // (not required here, for the time being, since it will be calculated
  //  before being required by FLUID solver)
  //double fact = (gamma_-1.0)/gamma_;
  //densdtnp_->Update(fact,*densdtn_,0.0);

  // same-density-derivative predictor
  //densnp_->Update(dta_,*densdtn_,1.0);
  // prediction of time derivative not required (would also not be required
  // to be performed, since we just updated the time derivatives of density,
  // and thus, densdtnp_ = densdtn_)

  return;
}


/*----------------------------------------------------------------------*
 | compute values at intermediate time steps                   vg 02/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeIntermediateValues()
{
  // computations for incremental case: phi at time n+alpha_F and time deriv.
  if (incremental_)
  {
    // computation of time derivative of phi:
    // phidt(n+1) = (phi(n+1)-phi(n)) / (gamma*dt) - (1/gamma -1)*phidt(n)
    // in conservative low-Mach-number case: density times phi instead of phi
    const double fact1 = 1.0/(gamma_*dta_);
    const double fact2 = (-1.0/gamma_) +1.0;
    phidtnp_->Update(1.0,*phidtn_,fact2);
    if (prbtype_ == "loma" and convform_ =="conservative")
    {
      phidtnp_->Multiply( fact1,*phinp_,*densnp_,1.0);
      phidtnp_->Multiply(-fact1,*phin_,*densn_,1.0);
    }
    else phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

    // calculation of time derivative of phi at n+alpha_M, stored on
    // history vector for comfortable later transport to element routine
    hist_->Update(alphaM_,*phidtnp_,(1.0-alphaM_),*phidtn_,0.0);

    // compute phi at n+alpha_F
    phiaf_->Update(alphaF_,*phinp_,(1.0-alphaF_),*phin_,0.0);
  }

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_-(1-alphaF_)*dta_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::AddNeumannToResidual()
{
  residual_->Update(genalphafac_*dta_,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::AddSpecificTimeIntegrationParameters(
  ParameterList& params)
{
  params.set("using stationary formulation",false);
  params.set("using generalized-alpha time integration",true);
  params.set("total time",time_-(1-alphaF_)*dta_);
  params.set("time factor",genalphafac_*dta_);
  params.set("alpha_F",alphaF_);

  if (prbtype_ == "loma")
    params.set("time derivative of thermodynamic pressure",thermpressdtaf_);

  if (incremental_) discret_->SetState("phinp",phiaf_);
  else              discret_->SetState("phinp",phin_);

  if (prbtype_ == "loma" and convform_ != "conservative")
       discret_->SetState("densnp",densam_);
  else discret_->SetState("densnp",densnp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
double SCATRA::TimIntGenAlpha::ComputeThermPressure()
{
  // compute temperature at n+alpha_F
  phiaf_->Update(alphaF_,*phinp_,(1.0-alphaF_),*phin_,0.0);

  // set scalar and density vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phiaf_);
  discret_->SetState("densnp",densaf_);

  // define element parameter list
  ParameterList eleparams;

  // provide velocity field (export to column map necessary for parallel evaluation)
  // SetState cannot be used since this Multivector is nodebased and not dofbased
  const Epetra_Map* nodecolmap = discret_->NodeColMap();
  RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
  LINALG::Export(*convel_,*tmp);
  eleparams.set("velocity field",tmp);

  // set action for elements
  eleparams.set("action","calc_domain_and_bodyforce");
  eleparams.set("total time",time_-(1-alphaF_)*dta_);

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(2));

  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint  = (*scalars)[0];
  double parbofint  = (*scalars)[1];

  // evaluate domain integral
  // set action for elements
  eleparams.set("action","calc_therm_press");

  // variables for integrals of velocity-divergence and diffusive flux
  double divuint = 0.0;
  double diffint = 0.0;
  eleparams.set("velocity-divergence integral",divuint);
  eleparams.set("diffusive-flux integral",     diffint);

  // evaluate velocity-divergence and rhs on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for 
  // thermodynamic pressure, since it is usually at the same boundary.
  vector<std::string> condnames;
  condnames.push_back("FluxCalculation");
  for (unsigned int i=0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condnames[i]);
  }

  // get integral values on this proc
  divuint = eleparams.get<double>("velocity-divergence integral");
  diffint = eleparams.get<double>("diffusive-flux integral");

  // get integral values in parallel case
  double pardivuint = 0.0;
  double pardiffint = 0.0;
  discret_->Comm().SumAll(&divuint,&pardivuint,1);
  discret_->Comm().SumAll(&diffint,&pardiffint,1);

  // clean up
  discret_->ClearState();

  // compute thermodynamic pressure (with specific heat ratio fixed to be 1.4)
  const double shr  = 1.4;
  const double divt = shr*pardivuint/pardomint;
  const double lhs  = alphaF_*genalphafac_*dta_*divt;
  const double rhs  = genalphafac_*dta_*(shr-1.0)*(pardiffint+parbofint)/pardomint;
  const double hist = thermpressn_
                     - (1.0 - alphaF_)*genalphafac_*dta_*divt*thermpressn_
                     + (1.0 - genalphafac_)*dta_*thermpressdtn_;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

  // compute time derivative of thermodynamic pressure
  // tpdt(n+1) = (tp(n+1)-tp(n))/(gamma*dt)+((gamma-1)/gamma)*tpdt(n)
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = (gamma_-1.0)/gamma_;
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  // generalized-alpha version for time derivative of thermodynamic pressure
  /*const double lhs  = alphaF_*gamma_*dta_*divt;
  const double rhs  = (shr-1.0)*(pardiffint+parbofint)/pardomint;
  const double hist = (alphaF_-1.0)*divt*thermpressn_
                     + (alphaM_-1.0+alphaF_*(gamma_-1.0)*dta_*divt)*thermpressdtn_;
  thermpressdtnp_ = (rhs + hist)/(alphaM_ + lhs);
  const double fact1 = gamma_*dta_;
  constdouble fact2 = (1.0-gamma_)*dta_;
  thermpressnp_ = thermpressn_ + fact1*thermpressdtnp_ + fact2*thermpressdtn_;*/

  // backward-Euler version for thermodynamic pressure
  /*const double lhs  = dta_*divt;
  const double rhs  = dta_*(shr-1.0)*(pardiffint+parbofint)/pardomint;
  const double hist = thermpressn_;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);
  thermpressdtnp_ = (thermpressnp_-thermpressn_)/dta_;*/

  // print out thermodynamic pressure
  if (myrank_ == 0) 
  {
    cout << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
    cout << "Data output for instationary thermodynamic pressure:" << endl;
    cout << "Velocity in-/outflow at indicated boundary: " << pardivuint << endl;
    cout << "Diffusive flux at indicated boundary: "       << pardiffint << endl;
    cout << "Thermodynamic pressure: "                     << thermpressnp_ << endl;
    cout << "+--------------------------------------------------------------------------------------------+" << endl;
  }

  // time derivative of thermodynamic pressure at n+alpha_F
  // -> required as right-hand-side contribution to temperature equation
  thermpressdtaf_ = alphaF_*thermpressdtnp_ + (1.0-alphaF_)*thermpressdtn_;

  return thermpressnp_;
}


/*----------------------------------------------------------------------*
 | compute time derivative of density for low-Mach-number flow vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeDensityDerivative()
{
  // densdt(n+1) = (dens(n+1)-dens(n))/(gamma*dt)+((gamma-1)/gamma)*densdt(n)
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = (gamma_-1.0)/gamma_;
  densdtnp_->Update(fact1,*densnp_,-fact1,*densn_ ,0.0);
  densdtnp_->Update(fact2,*densdtn_,1.0);

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 | and computation of acceleration for next time step          vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::Update()
{
  // compute time derivative of phi for non-incremental case:
  // phidt(n) = (phi(n)-phi(n-1)) / (gamma*dt(n)) - (1/gamma -1)*phidt(n-1)
  if (not incremental_)
  {
    const double fact1 = 1.0/(gamma_*dta_);
    const double fact2 = (-1.0/gamma_) +1.0;
    phidtn_->Update( fact1,*phinp_,-fact1,*phin_ ,fact2);
  }
  // time deriv. of this step becomes most recent time derivative of
  // last step for incremental solver
  else phidtn_->Update(1.0,*phidtnp_,0.0);

  // solution of this step becomes most recent solution of last step
  phin_->Update(1.0,*phinp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::UpdateThermPressure()
{
  thermpressn_   = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 | update density at n for low-Mach-number flow                vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::UpdateDensity()
{
  densn_->Update(1.0,*densnp_,0.0);
  densdtn_->Update(1.0,*densdtnp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                  vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::OutputRestart()
{
  // additional state vectors that are needed for generalized-alpha restart
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin", phin_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                             vg 11/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed for generalized-alpha restart
  reader.ReadVector(phinp_,  "phinp");
  reader.ReadVector(phin_,   "phin");
  reader.ReadVector(phidtn_, "phidtn");

  return;
}

/*----------------------------------------------------------------------*
 | Initialization procedure before the first time step         vg 11/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PrepareFirstTimeStep()
{
  ApplyDirichletBC(time_, phin_,phidtn_);
  CalcInitialPhidt();
  return;
}

/*----------------------------------------------------------------------*
 | calculate initial time derivative of phi at t=t_0            vg 11/08|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::CalcInitialPhidt()
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
    eleparams.set("action","calc_initial_time_deriv");
    // other parameters that are needed by the elements
    eleparams.set("using generalized-alpha time integration",true);
    eleparams.set("total time",time_);
    eleparams.set("time-step length",dta_);
    eleparams.set("time factor",genalphafac_*dta_);
    eleparams.set("alpha_F",alphaF_);
    eleparams.set("problem type",prbtype_);
    eleparams.set("incremental solver",incremental_);
    eleparams.set("form of convective term",convform_);
    eleparams.set("fs subgrid diffusivity",fssgd_);
    if (prbtype_=="elch")
    {
      // get ELCH-specific paramter F/RT (default value for the temperature is 298K)
      const double frt = 96485.3399/(8.314472 * params_->get<double>("TEMPERATURE",298.0));
      eleparams.set("frt",frt); // factor F/RT
    }
    else if (prbtype_ == "loma")
      eleparams.set("time derivative of thermodynamic pressure",thermpressdtn_);

    //provide velocity field (export to column map necessary for parallel evaluation)
    //SetState cannot be used since this Multivector is nodebased and not dofbased
    const Epetra_Map* nodecolmap = discret_->NodeColMap();
    RefCountPtr<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(*nodecolmap,3));
    LINALG::Export(*convel_,*tmp);
    eleparams.set("velocity field",tmp);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

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

  // copy values to phidtnp
  phidtnp_->Update(1.0,*phidtn_,0.0);

  // reset the matrix (and its graph!) since we solved
  // a very special problem here that has a different sparsity pattern_
  if (params_->get<int>("BLOCKPRECOND") )
    ; //how to reset a block matrix ??
  else
    SystemMatrix()->Reset();

  return;
}

/*----------------------------------------------------------------------*
 | set velocity field for low-Mach-number flow                 vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::SetLomaVelocity(RCP<const Epetra_Vector> extvel,
    RCP<DRT::Discretization> fluiddis)
{
  // store temperature and velocity of previous iteration for convergence check
  tempincnp_->Update(1.0,*phinp_,0.0);
  //velincnp_->Update(1.0,*convel_,0.0);

  // calculation of density fields at intermediate time steps
  densam_->Update(alphaM_,*densnp_,(1.0-alphaM_),*densn_,0.0);
  densaf_->Update(alphaF_,*densnp_,(1.0-alphaF_),*densn_,0.0);

  // check vector compatibility and determine space dimension
  int numdim =-1;
  if (extvel->MyLength()<= (4* convel_->MyLength()) and
      extvel->MyLength() > (3* convel_->MyLength()))
    numdim = 3;
  else if (extvel->MyLength()<= (3* convel_->MyLength()))
    numdim = 2;
  else
    dserror("fluid velocity vector too large");

  // get noderowmap of scatra discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // get dofrowmap of fluid discretization
  const Epetra_Map* dofrowmap = fluiddis->DofRowMap();

  // loop over local nodes of scatra discretization
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();lnodeid++)
  {
    // first of all, assume the present node is not a slavenode
    bool slavenode=false;

    // get the processor-local scatra node
    DRT::Node*  scatralnode = discret_->lRowNode(lnodeid);

    // get the processor-local fluid node
    DRT::Node*  fluidlnode = fluiddis->lRowNode(lnodeid);

    // the set of degrees of freedom associated with the fluid node
    vector<int> nodedofset = fluiddis->Dof(fluidlnode);

    // check whether we have a pbc condition on this scatra node
    vector<DRT::Condition*> mypbc;
    scatralnode->GetCondition("SurfacePeriodic",mypbc);

    // yes, we have a periodic boundary condition on this scatra node
    if (mypbc.size()>0)
    {
      // get master and list of all his slavenodes
      map<int, vector<int> >::iterator master = pbcmapmastertoslave_->find(scatralnode->Id());

      // check whether this is a slavenode
      if (master == pbcmapmastertoslave_->end())
      {
        // indeed a slavenode
        slavenode = true;
      }
      else
      {
        // we have a masternode: set values for all slavenodes
        vector<int>::iterator i;
        for(i=(master->second).begin();i!=(master->second).end();++i)
        {
          // global and processor-local scatra node ID for slavenode
          int globalslaveid = *i;
          int localslaveid  = noderowmap->LID(globalslaveid);

          // get the processor-local fluid slavenode
          DRT::Node*  fluidlslavenode = fluiddis->lRowNode(localslaveid);

          // the set of degrees of freedom associated with the node
          vector<int> slavenodedofset = fluiddis->Dof(fluidlslavenode);

          for(int index=0;index<numdim;++index)
          {
            // global and processor-local fluid dof ID
            int gid = slavenodedofset[index];
            int lid = dofrowmap->LID(gid);

            // get density for this processor-local scatra node
            double dens  = (*densaf_)[localslaveid];
            // get velocity for this processor-local fluid dof
            double velocity =(*extvel)[lid];
            // insert velocity*density-value in vector
            convel_->ReplaceMyValue(localslaveid, index, velocity*dens);
          }
        }
      }
    }

    // do this for all nodes other than slavenodes
    if (slavenode == false)
    {
      for(int index=0;index<numdim;++index)
      {
        // global and processor-local fluid dof ID
        int gid = nodedofset[index];
        int lid = dofrowmap->LID(gid);

        // get density for this processor-local scatra node
        double dens  = (*densaf_)[lnodeid];
        // get velocity for this processor-local fluid dof
        double velocity = (*extvel)[lid];
        // insert velocity*density-value in vector
        convel_->ReplaceMyValue(lnodeid, index, velocity*dens);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | return the right time-scaling-factor for the true residual gjb 02/09 |
 *----------------------------------------------------------------------*/
double SCATRA::TimIntGenAlpha::ResidualScaling() const
{
  if (incremental_) 
  return (alphaF_/(genalphafac_*dta_));
else 
  return (alphaF_ / ((1.0-alphaF_)*genalphafac_*dta_));
}

#endif /* CCADISCRET */
