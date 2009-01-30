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
: ScaTraTimIntImpl(actdis,solver,params,output),
  theta_(params_->get<double>("theta"))
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // temporal solution derivative at time n
  phidtn_ = LINALG::CreateVector(*dofrowmap,true);

  // only required for low-Mach-number flow
  if (prbtype_ == "loma")
  {
    // density at time n
    densn_  = LINALG::CreateVector(*dofrowmap,true);
    // density at time n-1 (only defined since required by fluid OST solver)
    densnm_ = LINALG::CreateVector(*dofrowmap,true);

    // time derivative of density at time n
    // (required for same-density-derivative predictor and conservative
    // formulation)
    densdtn_  = LINALG::CreateVector(*dofrowmap,true);

    // time derivative of thermodynamic pressure at n+1 and n
    // (computed if not constant, otherwise remaining zero)
    thermpressdtnp_ = 0.0;
    thermpressdtn_ = 0.0;
  }

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
  // low-Mach-number flow: hist_ = densn_*phin_+dt*(1-Theta)*densn_*phidtn_
  if (prbtype_ == "loma")
  {
    hist_->Multiply(1.0, *phin_, *densn_, 0.0);
    hist_->Multiply(dta_*(1.0-theta_), *phidtn_, *densn_, 1.0);

    // for conservative formulation:
    // hist_ = hist_+dt*(1-Theta)*phin_*densdtn_
    if (convform_ =="conservative")
      hist_->Multiply(dta_*(1.0-theta_), *densdtn_, *phin_, 1.0);
  }
  // else: hist_=phin_+dt*(1-Theta)*phidtn_
  else hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | compute initial time derivative of density field            vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeInitialDensityDerivative()
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
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ExplicitPredictor()
{
  phinp_->Update(dta_, *phidtn_,1.0);

  // for the electric potential we just use the 'old' value of last time step
  Teuchos::RCP<Epetra_Vector> onlypot = conpotsplitter_.ExtractCondVector(phin_);
  conpotsplitter_.InsertCondVector(onlypot, phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)

  // same-thermodynamic-pressure-derivative predictor (currently not used)
  //thermpressnp_ += dta_*thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | predict density for next time step for low-Mach-number flow vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PredictDensity()
{
  // same-density predictor (not required to be performed, since we just
  // updated the density field, and thus, densnp_ = densn_)

  // same-density-derivative predictor (currently not used)
  //densnp_->Update(dta_, *densdtn_,1.0);

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | reset the residual vector and add actual Neumann loads               |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddNeumannToResidual()
{
  residual_->Update(theta_*dta_,*neumann_loads_,0.0);
  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddSpecificTimeIntegrationParameters(
  ParameterList& params)
{
  params.set("using stationary formulation",false);
  params.set("using generalized-alpha time integration",false);
  params.set("total time",time_);
  params.set("time factor",theta_*dta_);
  params.set("alpha_F",1.0);

  if (prbtype_ == "loma")
    params.set("time derivative of thermodynamic pressure",thermpressdtnp_);

  discret_->SetState("densnp",densnp_);
  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
double SCATRA::TimIntOneStepTheta::ComputeThermPressure()
{
  // compute "history" part
  const double hist = thermpressn_ + (1.0-theta_)*dta_*thermpressdtn_;

  // set scalar and density vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  discret_->SetState("densnp",densnp_);

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
  eleparams.set("total time",time_);

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
  const double shr = 1.4;
  const double lhs = theta_*dta_*shr*pardivuint/pardomint;
  const double rhs = theta_*dta_*(shr-1.0)*(pardiffint+parbofint)/pardomint;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

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

  // compute time derivative of thermodynamic pressure at n+1
  // tpdt(n+1) = (tp(n+1)-tp(n))/(theta*dt)+((theta-1)/theta)*tpdt(n)
  double fact1 = 1.0/(theta_*dta_);
  double fact2 = (theta_-1.0)/theta_;
  thermpressdtnp_ = fact1*(thermpressnp_-thermpressn_) + fact2*thermpressdtn_;

  return thermpressnp_;
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
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateThermPressure()
{
  thermpressn_   = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 | update density at n for low-Mach-number flow                vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::UpdateDensity()
{
  // compute density derivative at time n if required for
  // same-density-derivative predictor or conservative formulation
  // densdt(n) = (dens(n)-dens(n-1))/(theta*dt)+((theta-1)/theta)*densdt(n-1)
  if (convform_ =="conservative")
  {
    double fact1 = 1.0/(theta_*dta_);
    double fact2 = (theta_-1.0)/theta_;
    densdtn_->Update(fact1,*densnp_,-fact1,*densn_ ,fact2);
  }

  densn_->Update(1.0,*densnp_,0.0);

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
    eleparams.set("action","calc_initial_time_deriv");
    // other parameters that are needed by the elements
    eleparams.set("using generalized-alpha time integration",false);
    eleparams.set("total time",time_);
    eleparams.set("time-step length",dta_);
    eleparams.set("time factor",theta_*dta_);
    eleparams.set("alpha_F",1.0);
    eleparams.set("problem type",prbtype_);
    eleparams.set("form of convective term",convform_);
    eleparams.set("fs subgrid diffusivity",fssgd_);
    if (prbtype_=="elch")
    {
      // get ELCH-specific paramter F/RT (already set in base class constructor)
      eleparams.set("frt",frt_);
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
void SCATRA::TimIntOneStepTheta::SetLomaVelocity(RCP<const Epetra_Vector> extvel,
    RCP<DRT::Discretization> fluiddis)
{
  // store temperature and velocity of previous iteration for convergence check
  tempincnp_->Update(1.0,*phinp_,0.0);
  //velincnp_->Update(1.0,*convel_,0.0);

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
            double dens  = (*densnp_)[localslaveid];
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
        double dens  = (*densnp_)[lnodeid];
        // get velocity for this processor-local fluid dof
        double velocity = (*extvel)[lid];
        // insert velocity*density-value in vector
        convel_->ReplaceMyValue(lnodeid, index, velocity*dens);
      }
    }
  }

  return;
}

#endif /* CCADISCRET */
