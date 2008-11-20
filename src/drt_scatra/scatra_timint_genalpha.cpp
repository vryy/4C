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

  // temporal solution derivative at time n
  phidtn_  = LINALG::CreateVector(*dofrowmap,true);

  if (prbtype_ == "loma")
  {
    // density at times n, n+alpha_M and n+alpha_F
    densn_  = LINALG::CreateVector(*dofrowmap,true);
    densam_ = LINALG::CreateVector(*dofrowmap,true);
    densaf_ = LINALG::CreateVector(*dofrowmap,true);

    // time derivative of density at times n+1 and n
    densdtnp_ = LINALG::CreateVector(*dofrowmap,true);
    densdtn_  = LINALG::CreateVector(*dofrowmap,true);
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
  // hist_ = phin_ + dt*(1-(gamma/alpha_M))*phidtn_
  //       = phin_ + dt*(1-genalphafac)*phidtn_
  hist_->Update(1.0, *phin_, dta_*(1.0-genalphafac_), *phidtn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                          vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ExplicitPredictor()
{
  // constant predictor
  phinp_ ->Update(1.0,*phin_,0.0);
  return;
}


/*----------------------------------------------------------------------*
 | predict density for next time step for low-Mach-number flow vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::PredictDensity()
{
  // same-density predictor (not required to be performed, since we just
  // updated the density field, and thus, densnp_ = densn_)
  // prediction of time derivative not required, for the time being
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
 | reset the residual vector and add actual Neumann loads               |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::AddNeumannToResidual()
{
  residual_->Update(alphaF_*genalphafac_*dta_,*neumann_loads_,0.0);
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
  params.set("time factor",genalphafac_*dta_);
  params.set("alpha_F",alphaF_);

  if (prbtype_ == "loma") discret_->SetState("densnp",densam_);
  else                    discret_->SetState("densnp",densnp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative of density for low-Mach-number flow vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntGenAlpha::ComputeDensityDerivative()
{
  // densdt(n+1) = (dens(n+1)-dens(n))/(gamma*dt)+((gamma-1)/gamma)*densdt(n)
  double fact1 = 1.0/(gamma_*dta_);
  double fact2 = (gamma_-1.0)/gamma_;
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
  // update time derivative of phi
  double fact1 = 1.0/(gamma_*dta_);
  double fact2 = (-1.0/gamma_) +1.0;

  // phidt(n) = (phi(n)-phi(n-1)) / (gamma*dt(n)) - (1/gamma -1)*phidt(n-1)
  phidtn_->Update( fact1,*phinp_,-fact1,*phin_ ,fact2);

  // we know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  ApplyDirichletBC(time_,Teuchos::null,phidtn_);

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

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
    eleparams.set("action","initialize_one_step_theta");
    // other parameters that are needed by the elements
    eleparams.set("using generalized-alpha time integration",true);
    eleparams.set("total time",time_);
    eleparams.set("time-step length",dta_);
    eleparams.set("time factor",genalphafac_*dta_);
    eleparams.set("alpha_F",alphaF_);
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
void SCATRA::TimIntGenAlpha::SetLomaVelocity(RCP<const Epetra_Vector> extvel,
    RCP<DRT::Discretization> fluiddis)
{
  // for generalized-alpha time integration, at first, density fields at
  // intermediate time steps have to be calculated
  densam_->Update((alphaM_),*densnp_,(1.0-alphaM_),*densn_,0.0);
  densaf_->Update((alphaF_),*densnp_,(1.0-alphaF_),*densn_,0.0);


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

#endif /* CCADISCRET */
