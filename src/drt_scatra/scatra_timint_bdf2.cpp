/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_bdf2.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "scatra_timint_bdf2.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntBDF2::TimIntBDF2(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,output),
  theta_(0.0)
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // state vector for solution at time t_{n-1}
  phinm_      = LINALG::CreateVector(*dofrowmap,true);

  // only required for low-Mach-number flow
  if (prbtype_ == "loma")
  {
    // density at times n and n-1
    densn_  = LINALG::CreateVector(*dofrowmap,true);
    densnm_ = LINALG::CreateVector(*dofrowmap,true);

    // time derivative of thermodynamic pressure at n+1
    // (computed if not constant, otherwise remaining zero)
    thermpressdtnp_ = 0.0;
  }

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntBDF2::~TimIntBDF2()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetOldPartOfRighthandside()
{
  /*
  BDF2: for variable time step:

                 hist_ = (1+omega)^2/(1+ 2*omega) * phin_
                           - omega^2/(1+ 2*omega) * phinm_

  BDF2: for constant time step:

                 hist_ = 4/3 phin_ - 1/3 phinm_

  For low-Mach-number flow, densn_*phin_ and densnm_*phinm_ are
  used instead of phin_ and phinm_, respectively.
  */
  if (step_>1)
  {
    double omega = dta_/dtp_;
    double fact1 = (1.0 + omega)*(1.0 + omega)/(1.0 + (2.0*omega));
    double fact2 = -(omega*omega)/(1+ (2.0*omega));

    // low-Mach-number flow: multiply by density
    if (prbtype_ == "loma")
    {
      hist_->Multiply(fact1, *phin_, *densn_, 0.0);
      hist_->Multiply(fact2, *phinm_, *densnm_, 1.0);
    }
    else hist_->Update(fact1, *phin_, fact2, *phinm_, 0.0);

    // for BDF2 theta is set by the timestepsizes, 2/3 for const. dt
    theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
  }
  else   // for start-up of BDF2 we do one step with backward Euler
  {
    // low-Mach-number flow: multiply by density
    if (prbtype_ == "loma") hist_->Multiply(1.0, *phin_, *densn_, 0.0);
    else                    hist_->Update(1.0, *phin_, 0.0);

    // backward Euler => use theta=1.0
    theta_=1.0;
  }

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ExplicitPredictor()
{
  if (step_>1) phinp_->Update(-1.0, *phinm_,2.0);
  // for step == 1 phinp_ is already correctly initialized with the
  // initial field phin_

  return;
}


/*----------------------------------------------------------------------*
 | predict density for next time step for low-Mach-number flow vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::PredictDensity()
{
  // same-density predictor (not required to be performed, since we just
  // updated the density field, and thus, densnp_ = densn_)

  // same-density-increment predictor (currently not used)
  //if (step_>1) densnp_->Update(-1.0,*densnm_,2.0);

  return;
}


/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)

  // same-thermodynamic-pressure-increment predictor (currently not used)
  //if (step_>1) thermpressnp_ = 2.0*thermpressn_ - thermpressnm_;

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | reset the residual vector and add actual Neumann loads               |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::AddNeumannToResidual()
{
  residual_->Update(theta_*dta_,*neumann_loads_,0.0);
  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::AddSpecificTimeIntegrationParameters(
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
double SCATRA::TimIntBDF2::ComputeThermPressure()
{
  // compute "history" part (start-up of BDF2: one step backward Euler)
  double hist = 0.0;
  if (step_>1)
  {
    double omega = dta_/dtp_;
    double fact1 = (1.0 + omega)*(1.0 + omega)/(1.0 + (2.0*omega));
    double fact2 = -(omega*omega)/(1+ (2.0*omega));

    hist = fact1*thermpressn_ + fact2*thermpressnm_;
  }
  else hist = thermpressn_;

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
  eleparams.set("action","calc_therm_press");

  // variables for integrals of velocity-divergence, rhs and domain
  double divuint = 0.0;
  double rhsint  = 0.0;
  double domint  = 0.0;
  eleparams.set("velocity-divergence integral",divuint);
  eleparams.set("rhs integral",                rhsint);
  eleparams.set("domain integral",             domint);

  // evaluate integrals of velocity-divergence, rhs and domain
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // get integral values on this proc
  divuint = eleparams.get<double>("velocity-divergence integral");
  rhsint = eleparams.get<double>("rhs integral");
  domint  = eleparams.get<double>("domain integral");

  // get integral values in parallel case
  double pardivuint = 0.0;
  double parrhsint = 0.0;
  double pardomint  = 0.0;
  discret_->Comm().SumAll(&divuint,&pardivuint,1);
  discret_->Comm().SumAll(&rhsint,&parrhsint,1);
  discret_->Comm().SumAll(&domint,&pardomint,1);

  // clean up
  discret_->ClearState();

  // compute thermodynamic pressure (with specific heat ratio fixed to be 1.4)
  const double shr = 1.4;
  const double lhs = theta_*dta_*shr*divuint/domint;
  const double rhs = theta_*dta_*(shr-1.0)*rhsint/domint;
  thermpressnp_ = (rhs + hist)/(1.0 + lhs);

  // print out thermodynamic pressure
  if (myrank_ == 0) cout << "Thermodynamic pressure: " << thermpressnp_ << endl;

  // compute time derivative of thermodynamic pressure at n+1
  // (start-up of BDF2: one step backward Euler)
  if (step_>1)
       thermpressdtnp_ = (3.0*thermpressnp_-4.0*thermpressn_+thermpressnm_)/(2.0*dta_);
  else thermpressdtnp_ = (thermpressnp_-thermpressn_)/dta_;

  return thermpressnp_;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::Update()
{
  // solution of this step becomes most recent solution of the last step
  phinm_->Update(1.0,*phin_ ,0.0);
  phin_ ->Update(1.0,*phinp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::UpdateThermPressure()
{
  thermpressnm_ = thermpressn_;
  thermpressn_  = thermpressnp_;

  return;
}


/*----------------------------------------------------------------------*
 | update density at n-1 and n for low-Mach-number flow        vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::UpdateDensity()
{
  densnm_->Update(1.0,*densn_ ,0.0);
  densn_->Update(1.0,*densnp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::OutputRestart()
{
  // additional state vectors that are needed for BDF2 restart
  output_->WriteVector("phin", phin_);
  output_->WriteVector("phinm", phinm_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/0 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed for BDF2 restart
  reader.ReadVector(phinp_,"phinp");
  reader.ReadVector(phin_, "phin");
  reader.ReadVector(phinm_,"phinm");

  return;
}


/*----------------------------------------------------------------------*
 | initialization procedure before the first time step         gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::PrepareFirstTimeStep()
{
  ApplyDirichletBC(time_, phin_,Teuchos::null);
  return;
}


/*----------------------------------------------------------------------*
 | set velocity field for low-Mach-number flow                 vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetLomaVelocity(RCP<const Epetra_Vector> extvel,
    RCP<DRT::Discretization> fluiddis)
{
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
