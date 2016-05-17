/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_hdg.cpp
\brief HDG time-integration scheme

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_hdg.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_scatra_ele/scatra_ele_hdg.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

#include "../drt_lib/drt_assemblestrategy.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntHDG::TimIntHDG(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<Teuchos::ParameterList>&   extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>        output
)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output),
  TimIntGenAlpha(actdis,solver,params,extraparams,output),
  nds_intvar_(2),
  intphinp_(Teuchos::null),
  intphin_(Teuchos::null),
  interpolatedPhinp_(Teuchos::null),
  timealgoset_(INPAR::SCATRA::timeint_gen_alpha),
  startalgo_(true),
  theta_(-1)

{
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::Init()
{
  DRT::DiscretizationHDG* hdgdis = dynamic_cast<DRT::DiscretizationHDG*>(discret_.get());
  if (hdgdis == NULL)
    dserror("Did not receive an HDG discretization");

  int elementndof = hdgdis->NumMyRowElements() > 0 ?
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG *>(hdgdis->lRowElement(0))->NumDofPerElementAuxiliary() : 0;


  // add proxy for interior degrees of freedom to scatra discretization
  if (discret_->BuildDofSetAuxProxy(0,elementndof,0,false) != 2)
    dserror("Scatra discretization has illegal number of dofsets!");
  discret_->FillComplete();


  // HDG vectors passed to the element
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(nds_intvar_);
  intphinp_  = LINALG::CreateVector(*intdofrowmap,true);
  intphin_  = LINALG::CreateVector(*intdofrowmap,true);

  // write number of degrees of freedom for hdg and interior variables to screen output
  if (discret_->Comm().MyPID()==0)
  {
    std::cout << "Number of degrees of freedom in HDG system: "
              << discret_->DofRowMap(0)->NumGlobalElements() << std::endl;
    std::cout << "Number of degrees of freedom of interior variables: "
            << discret_->DofRowMap(nds_intvar_)->NumGlobalElements() << std::endl;
  }

  // implement ost and bdf2 through gen-alpha facilities
  // TO DO: implement other time integration schemes, at the moment only one-step-theta implemented
  if (timealgo_ == INPAR::SCATRA::timeint_bdf2)
  { dserror("At the moment only one step theta implemented");
    alphaM_ = 1.5;
    alphaF_ = 1.0;
    gamma_ = 1.0;
  }
  else if (timealgo_ == INPAR::SCATRA::timeint_one_step_theta)
  {
    alphaM_ = 1.0;
    alphaF_ = 1.0;
    gamma_ = params_->get<double>("THETA");
  }
  else if (timealgo_ == INPAR::SCATRA::timeint_stationary)
    dserror("Stationary case not implemented for HDG");
  else
    dserror("At the moment only one step theta implemented");

  timealgoset_ = timealgo_;
  timealgo_ = INPAR::SCATRA::timeint_gen_alpha;

  // call Init()-functions of base classes
  // note: this order is important
  SCATRA::TimIntGenAlpha::Init();

  // create dofsets for concentration at nodes for output
  interpolatedPhinp_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);

}


/*------------------------------------------------------------------------*
| calculate pseudo-theta for startalgo_, modified for HDG  hoermann 09/15 |
*-------------------------------------------------------------------------*/
void SCATRA::TimIntHDG::SetTheta()
{
  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  // starting algorithm
  if (startalgo_ || (step_ <= 2 && timealgoset_ == INPAR::SCATRA::timeint_bdf2))
  {
    // use backward-Euler-type parameter combination
    if (step_ <= 1 && timealgoset_ == INPAR::SCATRA::timeint_bdf2)
    {
      if (myrank_==0)
      {
        std::cout<<"Starting algorithm for Af_GenAlpha active. "
            //<<"Performing step "<<step_ <<" of "<<numstasteps_
            <<" Backward Euler starting steps"<<std::endl;
      }
      alphaM_ = 1.0;
      alphaF_ = 1.0;
      gamma_  = 1.0;
    }
    else
    {
      // recall original user wish
      if (timealgoset_ == INPAR::SCATRA::timeint_one_step_theta)
      {
        alphaM_ = alphaF_ = 1.0;
        gamma_  = params_->get<double>("THETA");
      }
      else if (timealgoset_ == INPAR::SCATRA::timeint_bdf2)
      {
        alphaF_ = gamma_ = 1.0;
        alphaM_ = 3./2.;
      }
      else
      {
        alphaM_ = params_->get<double>("alpha_M");
        alphaF_ = params_->get<double>("alpha_F");
        gamma_  = params_->get<double>("gamma");
      }

      // do not enter starting algorithm section in the future
      startalgo_ = false;
    }
  }

  // compute "pseudo-theta" for af-generalized-alpha scheme
  theta_ = alphaF_*gamma_/alphaM_;
}


/*----------------------------------------------------------------------*
| set HDG state vectors                                  hoermann 09/15 |
*-----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
{
  // set hdg vector and interior variables vector
  discret_->SetState(0, "phin", phin_);
  discret_->SetState(0, "phiaf",phinp_);
  discret_->SetState(nds_intvar_, "intphinp",intphinp_);
  discret_->SetState(nds_intvar_, "intphin",intphin_);
}//AddTimeIntegrationSpecificVectors

/*----------------------------------------------------------------------*
 | compute values at intermediate time steps for gen.-alpha  hoer 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::GenAlphaIntermediateValues()
{
  // set intermediate values for concentration derivatives
  //
  //       n+alphaM                n+1                      n
  //  dtphi       = alpha_M * dtphi    + (1-alpha_M) * dtphi
  //       (i)                     (i)
  phidtam_->Update((alphaM_),*phidtnp_,(1.0-alphaM_),*phidtn_,0.0);

  // set intermediate values for concentration, concentration gradient
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  phiaf_->Update((alphaF_),*phinp_,(1.0-alphaF_),*phin_,0.0);

  phiam_->Update(alphaM_,*phinp_,(1.0-alphaM_),*phin_,0.0);

}//GenAlphaIntermediateValues

/*----------------------------------------------------------------------*
| set old part of right hand side                        hoermann 09/15 |
*-----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::SetOldPartOfRighthandside()
{
  SetTheta();
  hist_->PutScalar(0.0);

  // This code is entered at the beginning of the nonlinear iteration, so
  // store that the assembly to be done next is going to be the first one
  // (without combined vector update) for HDG.
//  SCATRA::TimIntGenAlpha::SetOldPartOfRighthandside();
}

/*----------------------------------------------------------------------*
 * Update
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::Update(const int num)
{
  SCATRA::TimIntGenAlpha::Update();

  // concentrations of this step become most recent
  // concentrations of the last step
  intphin_->Update(1.0,*intphinp_,0.0);


}//Update

namespace
{
  // internal helper function for output
  void getNodeVectorsHDG (DRT::Discretization &                 dis,
                          const Teuchos::RCP<Epetra_Vector> &   interiorValues,
                          const Teuchos::RCP<Epetra_Vector> &   traceValues,
                          const int                             ndim,
                          Teuchos::RCP<Epetra_Vector> &         phi,
                          Teuchos::RCP<Epetra_MultiVector> &    gradphi,
                          Teuchos::RCP<Epetra_Vector> &         tracephi,
                          int                                   nds_intvar_,
                          int                                   ndofs
                          )
  {
    dis.ClearState(true);

    // create dofsets for concentration at nodes
    tracephi.reset(new Epetra_Vector(phi->Map()));
    gradphi.reset(new Epetra_MultiVector(*dis.NodeRowMap(),ndim));

    // call element routine to interpolate HDG to elements
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::interpolate_hdg_to_node);
    dis.SetState(0,"phiaf",traceValues);
    dis.SetState(nds_intvar_,"intphinp",interiorValues);
    DRT::Element::LocationArray la(ndofs);
    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(dis.NumMyRowNodes());

    phi->PutScalar(0.);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      interpolVec.Size(ele->NumNode()*(2+ndim));

      ele->Evaluate(eleparams,dis,la,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i=0; i<ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = dis.NodeRowMap()->LID(node->Id());
        if (localIndex < 0)
          continue;
        touchCount[localIndex]++;
        (*phi)[localIndex] += interpolVec(i);
        (*tracephi)[localIndex] += interpolVec(i+ele->NumNode());
        for (int d=0; d<ndim; ++d)
          (*gradphi)[d][localIndex] += interpolVec(i+(d+2)*ele->NumNode());
      }
    }

    //build average of values
    for (int i=0; i<phi->MyLength(); ++i)
    {
        (*phi)[i] /= touchCount[i];
        (*tracephi)[i] /= touchCount[i];
        for (int d=0; d<ndim; ++d)
          (*gradphi)[d][i] /= touchCount[i];
    }
    dis.ClearState(true);
  }
}


/*----------------------------------------------------------------------*
 | output of solution vector to binio                     hoermann 09/15|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::OutputState()
{
  // output of solution

  Teuchos::RCP<Epetra_MultiVector> interpolatedGradPhi;
  Teuchos::RCP<Epetra_Vector> interpolatedtracePhi;
  // get (averaged) values at element nodes
  getNodeVectorsHDG(*discret_, intphinp_, phinp_, DRT::Problem::Instance()->NDim(),
      interpolatedPhinp_, interpolatedGradPhi, interpolatedtracePhi, nds_intvar_, discret_->NumDofSets());

  // write vector to output file
  output_->WriteVector("phi_hdg",interpolatedPhinp_, IO::DiscretizationWriter::nodevector);
  output_->WriteVector("gradphi_hdg",interpolatedGradPhi, IO::DiscretizationWriter::nodevector);
  output_->WriteVector("tracephi_hdg",interpolatedtracePhi, IO::DiscretizationWriter::nodevector);

  WriteProblemSpecificOutput(interpolatedPhinp_);

  SCATRA::TimIntGenAlpha::OutputState();

}//OutputState

/*----------------------------------------------------------------------*
 | output of solution vector to binio for restart         hoermann 09/15|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::OutputRestart()
{
  SCATRA::TimIntGenAlpha::OutputRestart();
  output_->WriteVector("intphinp",intphinp_);

}

/*----------------------------------------------------------------------*
 | read restart                                          hoermann 09/15 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::ReadRestart(const int step)
{
  dserror("Restart functionality not implemented yet");

  SCATRA::TimIntGenAlpha::ReadRestart(step);

  IO::DiscretizationReader reader(discret_,step);
  // read state vectors that are needed for hdg
  reader.ReadVector(intphinp_, "intphinp");

  intphin_->Update(1.0,*intphinp_,0.0);
  phin_->Update(1.0,*phinp_,0.0);

  return;
}

/*----------------------------------------------------------------------*
 |  set initial field for phi                            hoermann 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::SetInitialField(
    const INPAR::SCATRA::InitialField init,
    const int startfuncno)
{
  switch(init)
  {
  case INPAR::SCATRA::initfield_zero_field:
  {
    // set initial field to zero
    phin_->PutScalar(0.0);
    phinp_->PutScalar(0.0);
    intphin_->PutScalar(0.0);
    intphinp_->PutScalar(0.0);
    break;
  }
  case INPAR::SCATRA::initfield_field_by_function:
  {
    // set initial field defined by function
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::set_initial_field);
    eleparams.set<int>("funct",startfuncno);

    discret_->SetState("phiaf",phinp_);
    discret_->SetState("phin",phin_);
    discret_->SetState(nds_intvar_, "intphin",intphin_);
    discret_->SetState(nds_intvar_, "intphinp",intphinp_);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector updateVec1, updateVec2, dummyVec;
    DRT::Element::LocationArray la(discret_->NumDofSets());

    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    const Epetra_Map* intdofrowmap = discret_->DofRowMap(nds_intvar_);
    double error = 0;

    for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
    {
      updateVec2 = Teuchos::null;
      DRT::Element *ele = discret_->lColElement(iele);
      ele->LocationVector(*discret_,la,false);
      if (static_cast<std::size_t>(updateVec1.M()) != la[0].lm_.size())
        updateVec1.Shape(la[0].lm_.size(), 1);
      if (updateVec2.M() != discret_->NumDof(nds_intvar_,ele))
        updateVec2.Shape(discret_->NumDof(nds_intvar_,ele), 1);
      ele->Evaluate(eleparams,*discret_,la,dummyMat,dummyMat,updateVec1,updateVec2,dummyVec);

      if (ele->Owner() == discret_->Comm().MyPID())
      {
        std::vector<int> localDofs = discret_->Dof(nds_intvar_, ele);
        dsassert(localDofs.size() == static_cast<std::size_t>(updateVec2.M()), "Internal error");
        for (unsigned int i=0; i<localDofs.size(); ++i)
          localDofs[i] = intdofrowmap->LID(localDofs[i]);
        intphinp_->ReplaceMyValues(localDofs.size(), updateVec2.A(), &localDofs[0]);
      }

      // now fill the element vector into the discretization
      for (unsigned int i=0; i<la[0].lm_.size(); ++i)
      {
        const int lid = dofrowmap->LID(la[0].lm_[i]);
        if (lid >= 0)
        {
          // safety check if initial value for trace dof is set for all elements the same (interior face)
          if ((*phinp_)[lid] != 0)
            error += std::abs((*phinp_)[lid]-updateVec1(i));
          (*phinp_)[lid] = updateVec1(i);
          (*phin_)[lid] = updateVec1(i);
        }
      }
    }

    double globerror = 0;
    discret_->Comm().SumAll(&error, &globerror, 1);
    if (discret_->Comm().MyPID() == 0)
      std::cout << "Error project when setting face twice: " << globerror << std::endl;

    // initialize also the solution vector. These values are a pretty good guess for the
    // solution after the first time step (much better than starting with a zero vector)
    intphin_->Update(1.0,*intphinp_ ,0.0);

    break;
  }

  default:
    dserror("Option for initial field not implemented: %d", init); break;
  } // switch(init)

}//SetInitialField


/*----------------------------------------------------------------------*
 | calculate intermediate solution                        hoermann 09/15|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::ComputeIntermediateValues()
{
  // time derivatives are not independent but rather have to be computed
  // from phinp_, phin_ and phidtn_
  GenAlphaComputeTimeDerivative();
  // compute values at intermediate time steps
  GenAlphaIntermediateValues();

  return;
}

/*----------------------------------------------------------------------*
 | compute values at the interior of the elements         hoermann 09/15|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::ComputeInteriorValues()
{
  //Update the interior variables
  UpdateInteriorVariables(intphinp_);
  return;
}

/*----------------------------------------------------------------------*
 | update time derivative for gen-alpha time integration hoermann 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::GenAlphaComputeTimeDerivative()
{

  //                                  n+1     n
  //                               phi   - phi
  //       n+1      n  gamma-1.0      (i)
  //phidt    = phidt * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // compute factors
  const double fact1 = 1.0/(gamma_*dta_);
  const double fact2 = 1.0 - (1.0/gamma_);

  phidtnp_->Update(fact2,*phidtn_,0.0);
  phidtnp_->Update(fact1,*phinp_,-fact1,*phin_,1.0);

}//GenAlphaComputeTimeDerivative


/*----------------------------------------------------------------------*
 | update interior variables                             hoermann 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::UpdateInteriorVariables(
    Teuchos::RCP<Epetra_Vector>          updatevector)
{

  discret_->ClearState(true);
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::update_interior_variables);
  discret_->SetState("phiaf",phinp_);
  discret_->SetState("phin",phin_);
  discret_->SetState(nds_intvar_, "intphin",intphin_);
  discret_->SetState(nds_intvar_, "intphinp",intphinp_);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;
  Epetra_SerialDenseVector updateVec;
  DRT::Element::LocationArray la(discret_->NumDofSets());
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(nds_intvar_);

  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::Element *ele = discret_->lColElement(iele);
    ele->LocationVector(*discret_,la,false);
      updateVec.Shape(discret_->NumDof(nds_intvar_,ele),1);

    ele->Evaluate(eleparams,*discret_,la,dummyMat,dummyMat,updateVec,dummyVec,dummyVec);

    if (ele->Owner() == discret_->Comm().MyPID())
    {
      std::vector<int> localDofs = discret_->Dof(nds_intvar_, ele);
      dsassert(localDofs.size() == static_cast<std::size_t>(updateVec.M()), "Internal error");
      for (unsigned int i=0; i<localDofs.size(); ++i)
      {
        localDofs[i] = intdofrowmap->LID(localDofs[i]);
      }
      updatevector->ReplaceMyValues(localDofs.size(), updateVec.A(), &localDofs[0]);
    }
  }

  discret_->ClearState(true);

}


/*-------------------------------------------------------------------------------------*
 | finite difference check for system matrix (for debugging only)       hoermann 09/15 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::TimIntHDG::FDCheck()
{
  // create mat for finite difference approximation
//  Epetra_SerialDenseMatrix fdmat(phinp_->MyLength(),phinp_->MyLength());

  for (int i = 0; i<16; ++i)
  {
      double eps = 0.1;
      for(int j=0; j<i; ++j)
         eps *= 0.1;


    // make a copy of state variables to undo perturbations later
    Teuchos::RCP<Epetra_Vector> phinp_original = Teuchos::rcp(new Epetra_Vector(*phinp_));

    discret_->ClearState(true);

    const Epetra_Map* dofrowmap = discret_->DofRowMap(0);
    const Epetra_Map* intdofrowmap = discret_->DofRowMap(nds_intvar_);

    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, systemmatrix2;
    Teuchos::RCP<Epetra_Vector>          systemvector1, systemvector2, systemvector3;

    // create matrix and vector for calculation of sysmat and assemble
    systemmatrix1 = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()),27));
    systemvector1 = LINALG::CreateVector(*dofrowmap,true);
    DRT::AssembleStrategy strategy( 0, 0, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3 );

    // fill state vector with original state variables
    phinp_->Update(1.,*phinp_original,0.);

    // make temporary vector for interior variables for the update of the last step and use this vector
    // for the calculation of the original residual
    // the temporary vector is necessary because afterwards we need to calculate also the interior vectors
    // with the state vector with perturbation without influence for this update
    Teuchos::RCP<Epetra_Vector> intphitemp;
    intphitemp = LINALG::CreateVector(*intdofrowmap,true);

    strategy.Zero();

    // calculate of residual vector
    UpdateInteriorVariables(intphitemp);

    discret_->ClearState(true);
    discret_->SetState("phiaf",phinp_);
    discret_->SetState(nds_intvar_, "intphin",intphin_);
    discret_->SetState(0, "phin", phin_);
    discret_->SetState(nds_intvar_, "intphinp",intphitemp);
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_mat_and_rhs);
    DRT::Element::LocationArray la(discret_->NumDofSets());

    // loop over elements
    for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
    {
      DRT::Element *ele = discret_->lColElement(iele);
      ele->LocationVector(*discret_,la,false);

      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      // evaluate
      ele->Evaluate(eleparams,*discret_,la,
                       strategy.Elematrix1(),
                       strategy.Elematrix2(),
                       strategy.Elevector1(),
                       strategy.Elevector2(),
                       strategy.Elevector3());
      int eid = ele->Id();
      strategy.AssembleMatrix1( eid, la[0].lm_, la[0].lm_, la[0].lmowner_, la[0].stride_ );
      strategy.AssembleVector1( la[0].lm_, la[0].lmowner_ );
    }
    strategy.Complete();

    // make a copy of system matrix as Epetra_CrsMatrix
    Teuchos::RCP<Epetra_CrsMatrix> sysmatcopy = Teuchos::null;
    sysmatcopy =  (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::SparseMatrix>(systemmatrix1))))->EpetraMatrix();
    sysmatcopy->FillComplete();

    // make a copy of system right-hand side vector
    Teuchos::RCP<Epetra_Vector> residualVec = Teuchos::rcp(new Epetra_Vector(*systemvector1));
    Teuchos::RCP<Epetra_Vector> fdvec = LINALG::CreateVector(*dofrowmap,true);

    // initialize tracking variable for maximum absolute and relative errors
    double maxabserr(0.);
    double maxrelerr(0.);

    // calculate fd matrix
    for (int i= 0; i<phinp_->MyLength(); ++i)
    {
      strategy.Zero();

      // fill state vector with original state variables
      phinp_->Update(1.,*phinp_original,0.);

      // impose perturbation and update interior variables
      phinp_->SumIntoGlobalValue(i,0,eps);
      UpdateInteriorVariables(intphitemp);

      discret_->ClearState(true);

      discret_->SetState("phiaf",phinp_);
      discret_->SetState(nds_intvar_, "intphin",intphin_);
      discret_->SetState(0, "phin", phin_);
      discret_->SetState(nds_intvar_, "intphinp",intphitemp);

      Teuchos::ParameterList eleparams;
      eleparams.set<int>("action",SCATRA::calc_mat_and_rhs);

      DRT::Element::LocationArray la(discret_->NumDofSets());

      for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
      {
        DRT::Element *ele = discret_->lColElement(iele);
        ele->LocationVector(*discret_,la,false);

        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );
        ele->Evaluate(eleparams,*discret_,la,
                         strategy.Elematrix1(),
                         strategy.Elematrix2(),
                         strategy.Elevector1(),
                         strategy.Elevector2(),
                         strategy.Elevector3());
        int eid = ele->Id();
        strategy.AssembleMatrix1( eid, la[0].lm_, la[0].lm_, la[0].lmowner_, la[0].stride_ );
        strategy.AssembleVector1( la[0].lm_, la[0].lmowner_ );
      }
      strategy.Complete();
      fdvec->PutScalar(0.0);

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      for (int j=0; j<phinp_->MyLength(); ++j)
        (*fdvec)[j] = (*systemvector1)[j]/eps-(*residualVec)[j]/eps;

      for(int row=0; row<discret_->DofRowMap()->NumMyElements(); ++row)
      {
        // get current entry in original system matrix
        double entry(0.);
        int length = sysmatcopy->NumMyEntries(row);
        int numentries;
        std::vector<double> values(length);
        std::vector<double> valnew(length);
        std::vector<int> indices(length);
        sysmatcopy->ExtractMyRowCopy(row,length,numentries,&values[0],&indices[0]);
        for(int ientry=0; ientry<length; ++ientry)
        {
          if(sysmatcopy->ColMap().GID(indices[ientry]) == i)
          {
            entry = values[ientry];

            // absolute and relative errors in first comparison
            const double abserr1 = entry + (*fdvec)[row];
            double relerr1 = 0;
//            abserrmat(col,row) = abs(abserr1);
            if(abs(entry) > 1.e-17)
              relerr1 = abserr1 / abs(entry);
            else if(abs((*fdvec)[row]) > 1.e-17)
              relerr1 = abserr1 / abs((*fdvec)[row]);
            // store max abs and rel error
            if(abs(abserr1) > maxabserr)
              maxabserr = abs(abserr1);
            if(abs(relerr1) > maxrelerr)
              maxrelerr = abs(relerr1);
            break;
          }
        }
      }
    }
    // end calculate fd matrix

    // screen output
    if(myrank_ == 0)
    {
      std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SCATRA HDG SYSTEM MATRIX" << std::endl;
      std::cout << "EPS:        " << eps << std::endl;
      // FD output: sysmat, FD matrix and difference
  //    std::cout << std::endl << "SYSTEM MATRIX" << std::endl;
  //    (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix1)()->EpetraMatrix())->Print(std::cout);

  //    std::cout << std::endl << "FD MATRIX" << std::endl;
  //    std::cout << fdmat << std::endl;

  //    std::cout << std::endl << "ABSOLUT ERROR BETWEEN SYSTEM MATRIX AND FD MATRIX" << std::endl;
  //    std::cout << abserrmat << std::endl;

  //    std::cout << std::endl << "RELATIVE ERROR BETWEEN SYSTEM MATRIX AND FD MATRIX" << std::endl;
  //    std::cout << relerrmat << std::endl;

  //    std::cout << std::endl << "MAXIMUM ABSOLUTE AND RELATIVE ERRORS" << std::endl;
      std::cout << "ABSOLUT: " << maxabserr << std::endl;
      std::cout << "RELATIVE: " << maxrelerr << std::endl;
    }


  }
  dserror("FD check END");

}


/*----------------------------------------------------------------------*
 | prepare time loop                                     hoermann 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::PrepareTimeLoop()
{
  // call base class routine
  ScaTraTimIntImpl::PrepareTimeLoop();

  // check validity of material and element formulation
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_mat_initial);

  discret_->SetState("phiaf",phinp_);
  discret_->SetState("phin",phin_);
  discret_->SetState(nds_intvar_, "intphin",intphin_);
  discret_->SetState(nds_intvar_, "intphinp",intphinp_);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;
  DRT::Element::LocationArray la(discret_->NumDofSets());

  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::Element *ele = discret_->lColElement(iele);
    ele->LocationVector(*discret_,la,false);
    ele->Evaluate(eleparams,*discret_,la,dummyMat,dummyMat,dummyVec,dummyVec,dummyVec);
  }

  return;
} // SCATRA::TimIntHDG::PrepareTimeLoop
