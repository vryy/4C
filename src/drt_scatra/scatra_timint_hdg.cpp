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

#include "../drt_lib/drt_dofset.H"
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

#include "../drt_binstrategy/binning_strategy.H"

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
  theta_(-1),
  hdgdis_(NULL),
  padaptivity_(DRT::INPUT::IntegralValue<bool>(*params,"PADAPTIVITY")),
  padapterrortol_(params->get<double>("PADAPTERRORTOL")),
  padapterrorbase_(params->get<double>("PADAPTERRORBASE")),
  padaptdegreemax_(params->get<int>("PADAPTDEGREEMAX")),
  elementdegree_(Teuchos::null)

{
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::Setup()
{
  hdgdis_ = dynamic_cast<DRT::DiscretizationHDG*>(discret_.get());
  if (hdgdis_ == NULL)
    dserror("Did not receive an HDG discretization");

  // vector to store the dofs per element
  const Teuchos::RCP<Epetra_IntVector> eledofs = Teuchos::rcp(new Epetra_IntVector(*discret_->ElementColMap()));

  // loop over elements
  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::ELEMENTS::ScaTraHDG *hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG *>(discret_->lColElement(iele));
    (*eledofs)[iele] = hdgele->NumDofPerElementAuxiliary();
  }

  // add proxy for interior degrees of freedom to scatra discretization
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0,eledofs,0,false));
  if (discret_->AddDofSet(dofsetaux) != 2)
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
  SCATRA::TimIntGenAlpha::Setup();

  // create vector for concentration at nodes for output
  interpolatedPhinp_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);

  // vector to store the elementdegree at each time step
  elementdegree_ = LINALG::CreateVector(*(discret_->ElementRowMap()),true);

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

  if(padaptivity_)
    AdaptDegree();

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
  output_->WriteVector("phi_hdg",interpolatedPhinp_, IO::nodevector);
  output_->WriteVector("gradphi_hdg",interpolatedGradPhi, IO::nodevector);
  output_->WriteVector("tracephi_hdg",interpolatedtracePhi, IO::nodevector);

  WriteProblemSpecificOutput(interpolatedPhinp_);

  output_->WriteVector("elementdegree",elementdegree_, IO::elementvector);

}//OutputState

/*----------------------------------------------------------------------*
 | output of solution vector to binio for restart         hoermann 09/15|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::OutputRestart() const
{
  SCATRA::TimIntGenAlpha::OutputRestart();
  output_->WriteVector("intphinp",intphinp_);
  output_->WriteVector("phinp_trace",phinp_);
  output_->WriteVector("intphin",intphin_);

  output_->WriteMesh(step_,time_); // add info to control file for reading all variables in restart

}

/*----------------------------------------------------------------------*
 | read restart                                          hoermann 09/15 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input)
{

  IO::DiscretizationReader reader(discret_,step);

  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadHistoryData(step); // Read all saved data in nodes and elements and call nodal and element Unpacking each global variable has to be read

  if(padaptivity_)
  {
    // redistribute discr. with help of binning strategy
    if(discret_->Comm().NumProc()>1)
    {
      // create vector of discr.
      std::vector<Teuchos::RCP<DRT::Discretization> > dis;
      dis.push_back(discret_);

      //binning strategy for parallel redistribution
      Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy;

      std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
      std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

      // binning strategy is created and parallel redistribution is performed
      binningstrategy = Teuchos::rcp( new BINSTRATEGY::BinningStrategy() );
      binningstrategy->Init(dis);
      binningstrategy->WeightedPartitioning(dis,stdelecolmap,stdnodecolmap);
    }
  }

  // vector to store the dofs per element
  const Teuchos::RCP<Epetra_IntVector> eledofs = Teuchos::rcp(new Epetra_IntVector(*discret_->ElementColMap()));

  // build new maps for face dofs with adapted element order
  hdgdis_->BuildFaces();
  hdgdis_->BuildFaceRowMap();
  hdgdis_->BuildFaceColMap();

  // assign the degrees of freedom to the adapted dofsets
  hdgdis_->AssignDegreesOfFreedom(0);

  // replace all ghosted element with the original thus the correct polynomial degree is used
  discret_->ExportColumnElements(*discret_->ElementColMap(),false, false);

  hdgdis_->FillComplete();

  // store the number of dofs per element on vector
  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::ELEMENTS::ScaTraHDG *hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG *>(discret_->lColElement(iele));
    // store the number of dofs for the element
    (*eledofs)[iele] = hdgele->NumDofPerElementAuxiliary();
  }

  // create new local dofset for the new interior element dofs with adapted element order
  Teuchos::RCP<DRT::DofSetPredefinedDoFNumber> eledofs_new = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0,eledofs,0,false));
  // replace old interior element dofs with the new created dofset
  discret_->ReplaceDofSet(nds_intvar_,eledofs_new,false);

  hdgdis_->AssignDegreesOfFreedom(0);

  // clear map cache since after every FillComplete() / AssignDegreesOfFreedom() old maps are stored in the mapstack
  output_->ClearMapCache();

  // reset the residual, increment and sysmat to the size
  residual_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  increment_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  neumann_loads_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  sysmat_ = Teuchos::null;
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()),27));

  // reset the state vectors
  intphinp_.reset(new Epetra_Vector(*(discret_->DofRowMap(nds_intvar_)),true));
  intphin_.reset(new Epetra_Vector(*(discret_->DofRowMap(nds_intvar_)),true));
  phinp_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  phin_.reset(new Epetra_Vector(*(discret_->DofRowMap())));

  // read state vectors that are needed for hdg
  reader.ReadVector(phinp_, "phinp_trace");
  reader.ReadVector(intphinp_, "intphinp");

  intphin_->Update(1.0,*intphinp_,0.0);
  phin_->Update(1.0,*phinp_,0.0);

  // reset vector
  interpolatedPhinp_.reset(new Epetra_Vector(*(discret_->NodeRowMap())));
  elementdegree_.reset(new Epetra_Vector(*(discret_->ElementRowMap())));

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
      DRT::Element *ele = discret_->lColElement(iele);
      ele->LocationVector(*discret_,la,false);
      if (static_cast<std::size_t>(updateVec1.M()) != la[0].lm_.size())
        updateVec1.Shape(la[0].lm_.size(), 1);
      else
        memset(updateVec1.Values(),0.,la[0].lm_.size()*sizeof(double));
      if (updateVec2.M() != discret_->NumDof(nds_intvar_,ele))
        updateVec2.Shape(discret_->NumDof(nds_intvar_,ele), 1);
      else
        memset(updateVec2.Values(),0.,discret_->NumDof(nds_intvar_,ele)*sizeof(double));
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
    if (ele->Owner() != discret_->Comm().MyPID()) continue;

    ele->LocationVector(*discret_,la,false);
      updateVec.Shape(discret_->NumDof(nds_intvar_,ele),1);

    ele->Evaluate(eleparams,*discret_,la,dummyMat,dummyMat,updateVec,dummyVec,dummyVec);

    std::vector<int> localDofs = discret_->Dof(nds_intvar_, ele);
    dsassert(localDofs.size() == static_cast<std::size_t>(updateVec.M()), "Internal error");
    for (unsigned int i=0; i<localDofs.size(); ++i)
    {
      localDofs[i] = intdofrowmap->LID(localDofs[i]);
    }
    updatevector->ReplaceMyValues(localDofs.size(), updateVec.A(), &localDofs[0]);
  }

  discret_->ClearState(true);

}


/*-------------------------------------------------------------------------------------*
 | finite difference check for system matrix (for debugging only)       hoermann 09/15 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::TimIntHDG::FDCheck()
{
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

  for (int k = 0; k<16; ++k)
  {
    double eps = 1000;
    for(int j=0; j<k; ++j)
       eps *= 0.1;

    // initialize tracking variable for maximum absolute and relative errors
    double maxabserr(0.);
    double maxrelerr(0.);

    // calculate fd matrix
    for (int colgid=0; colgid<=sysmatcopy->ColMap().MaxAllGID(); ++colgid)
    {
      // check whether current column index is a valid global column index and continue loop if not
      int collid(sysmatcopy->ColMap().LID(colgid));
      int maxcollid(-1);
      discret_->Comm().MaxAll(&collid,&maxcollid,1);
      if(maxcollid < 0)
        continue;

      strategy.Zero();

      // fill state vector with original state variables
      phinp_->Update(1.,*phinp_original,0.);

      // impose perturbation and update interior variables
      if(phinp_->Map().MyGID(colgid))
        if(phinp_->SumIntoGlobalValue(colgid,0,eps))
          dserror("Perturbation could not be imposed on state vector for finite difference check!");
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
      fdvec->Scale(0.0);

      // finite difference suggestion (first divide by epsilon and then subtract for better conditioning)
      for (int j=0; j<phinp_->MyLength(); ++j)
        (*fdvec)[j] = - (*systemvector1)[j]/eps + (*residualVec)[j]/eps;

      for(int rowlid=0; rowlid<discret_->DofRowMap()->NumMyElements(); ++rowlid)
      {
        // get global index of current matrix row
        const int rowgid = sysmatcopy->RowMap().GID(rowlid);
        if(rowgid < 0)
          dserror("Invalid global ID of matrix row!");

        // get current entry in original system matrix
        double entry(0.);
        int length = sysmatcopy->NumMyEntries(rowlid);
        int numentries;
        std::vector<double> values(length);
        std::vector<int> indices(length);
        sysmatcopy->ExtractMyRowCopy(rowlid,length,numentries,&values[0],&indices[0]);

        for(int ientry=0; ientry<length; ++ientry)
        {
          if(sysmatcopy->ColMap().GID(indices[ientry]) == colgid)
          {
            entry = values[ientry];
            break;
          }
        }

        // absolute and relative errors in first comparison
        const double abserr1 = entry - (*fdvec)[rowlid];
        double relerr1 = 0;
        if(abs(entry) > 1.e-17)
          relerr1 = abserr1 / abs(entry);
        else if(abs((*fdvec)[rowlid]) > 1.e-17)
          relerr1 = abserr1 / abs((*fdvec)[rowlid]);
        // store max abs and rel error
        if(abs(abserr1) > maxabserr)
          maxabserr = abs(abserr1);
        if(abs(relerr1) > maxrelerr)
          maxrelerr = abs(relerr1);
      }
    }
    // end calculate fd matrix

    // screen output
    if(myrank_ == 0)
    {
      std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SCATRA HDG SYSTEM MATRIX" << std::endl;
      std::cout << "EPS:        " << eps << std::endl;
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

  // calculate matrices on element
  CalcMatInitial();

} // SCATRA::TimIntHDG::PrepareTimeLoop


/*----------------------------------------------------------------------*
 | calculate matrices on element                        hoermann 07/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::CalcMatInitial()
{

  TEUCHOS_FUNC_TIME_MONITOR("SCATRA::TimIntHDG::CalcMat");

  discret_->ClearState(true);

  // check validity of material and element formulation
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_mat_initial);

  discret_->SetState("phiaf",phinp_);
  discret_->SetState("phin",phin_);
  discret_->SetState(nds_intvar_, "intphin",intphin_);
  discret_->SetState(nds_intvar_, "intphinp",intphinp_);

  Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, systemmatrix2;
  Teuchos::RCP<Epetra_Vector>          systemvector1, systemvector2, systemvector3;

  // create matrix and vector for calculation of sysmat and assemble
  systemmatrix1 = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()),27));
  DRT::AssembleStrategy strategy( 0, 0, sysmat_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null );

  strategy.Zero();
  DRT::Element::LocationArray la(discret_->NumDofSets());

//    // get cpu time
//    const double tcmatinit = Teuchos::Time::wallTime();

  // loop over elements
  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::Element *ele = discret_->lColElement(iele);

    //if the element has only ghosted nodes it will not assemble -> skip evaluation
    if(ele->HasOnlyGhostNodes(discret_->Comm().MyPID())) continue;
    ele->LocationVector(*discret_,la,false);

    strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

    // evaluate
    int err = ele->Evaluate(eleparams,*discret_,la,
                     strategy.Elematrix1(),
                     strategy.Elematrix2(),
                     strategy.Elevector1(),
                     strategy.Elevector2(),
                     strategy.Elevector3());
    if (err) dserror("Proc %d: Element %d returned err=%d",discret_->Comm().MyPID(),ele->Id(),err);

    int eid = ele->Id();
    strategy.AssembleMatrix1( eid, la[0].lm_, la[0].lm_, la[0].lmowner_, la[0].stride_ );
  }
  sysmat_->Complete();

//    // end time measurement for element
//    double dtmatinit=Teuchos::Time::wallTime()-tcmatinit;
//    std::cout << "Time measurement evaluate: " << dtmatinit << std::endl;

  // Output of non-zeros in system matrix
  if(step_==0 and discret_->Comm().MyPID()==0)
  {
    int numglobalnonzeros = SystemMatrix()->EpetraMatrix()->NumGlobalNonzeros();
    std::cout << "Number of non-zeros in system matrix: " <<  numglobalnonzeros << std::endl;
  }

  return;
} // SCATRA::TimIntHDG::CalcMatIntitial


/*----------------------------------------------------------------------*
 | adapt degree of test function on element               hoermann 07/16|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::AdaptDegree()
{

  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + adapt degree");

  // get cpu time
  const double tcadapt = Teuchos::Time::wallTime();

  // cast and check if hdg discretization is provided
  DRT::DiscretizationHDG* hdgdis = dynamic_cast<DRT::DiscretizationHDG*>(discret_.get());
  if (hdgdis == NULL)
    dserror("Did not receive an HDG discretization");

  // vector to store the dofs per single element
  const Teuchos::RCP<Epetra_IntVector> eledofs = Teuchos::rcp(new Epetra_IntVector(*discret_->ElementColMap()));

  // vector to store the location array of the dofsets before the adaption with the new order
  std::vector<DRT::Element::LocationArray> la_old;

  // copy the old face dof map and the old interior element dof map
  Teuchos::RCP<Epetra_Map> facedofs_old = Teuchos::rcp( new Epetra_Map(*discret_->DofColMap(0)) );
  Teuchos::RCP<Epetra_Map> eledofs_old = Teuchos::rcp( new Epetra_Map(*discret_->DofColMap(nds_intvar_)) );

  // set action
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_padaptivity);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;

  discret_->SetState("phiaf",phinp_);
  discret_->SetState(nds_intvar_, "intphinp",intphinp_);

  // get cpu time
//  const double tccalcerr = Teuchos::Time::wallTime();


  // store if degree changes
  int degchange(0);

  // loop over elements
  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    // add new location array in vector for each element
    la_old.push_back(DRT::Element::LocationArray(discret_->NumDofSets()));

    DRT::Element *ele = discret_->lColElement(iele);

    // fill location array and store it for later use
    ele->LocationVector(*discret_,la_old[iele],false);

    if(ele->Owner() == discret_->Comm().MyPID())
    {
      DRT::ELEMENTS::ScaTraHDG *hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG *>(discret_->lColElement(iele));

      // call routine on elements to calculate error on element
      ele->Evaluate(eleparams,*discret_,la_old[iele],dummyMat,dummyMat,dummyVec,dummyVec,dummyVec);

      double error = eleparams.get<double>("error");
      double errorlog=0;

      if(error < 0)
        dserror("Error is negative!");

      if(error > 0)
        errorlog = log(error/padapterrortol_);
      else
        errorlog = 0.;

      int deg = hdgele->Degree() + ceil(errorlog/padapterrorbase_);

      if(deg < 0)
        deg = 0;
      else if(deg > padaptdegreemax_)
        deg = padaptdegreemax_;

      if (hdgele->Degree() != deg)
        degchange = 1;

        // set degree on element
      hdgele->SetDegree(deg);

      // store element degree (only for output)
      const int eleIndex = discret_->ElementRowMap()->LID(ele->Id());
      if (eleIndex >= 0)
        (*elementdegree_)[eleIndex] = deg;
    }
  }

  int degchangeall;
  discret_->Comm().SumAll(&degchange,&degchangeall,1);

  if (!degchangeall)
    return;

  PackMaterial();

//  // end time measurement for element
//  double dtcalcerr=Teuchos::Time::wallTime()-tccalcerr;
//  std::cout << "Time measurement for error calculation: " << dtcalcerr << std::endl;

//  // get cpu time
//  const double tcfillcomplete = Teuchos::Time::wallTime();

  // number of dofset in location array
  int nds_intvar_old(discret_->NumDofSets()+2);
  int nds_var_old(discret_->NumDofSets());

  // build new maps for face dofs with adapted element order
  hdgdis_->BuildFaces();
  hdgdis_->BuildFaceRowMap();
  hdgdis_->BuildFaceColMap();

  // assign the degrees of freedom to the adapted dofsets
  hdgdis_->AssignDegreesOfFreedom(0);

  // replace all ghosted element with the original thus the correct polynomial degree is used
  discret_->ExportColumnElements(*discret_->ElementColMap(),false, false);

  hdgdis_->FillComplete();

  // store the number of dofs per element on vector
  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::ELEMENTS::ScaTraHDG *hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG *>(discret_->lColElement(iele));
    // store the number of dofs for the element
    (*eledofs)[iele] = hdgele->NumDofPerElementAuxiliary();
  }

  // create new local dofset for the new interior element dofs with adapted element order
  Teuchos::RCP<DRT::DofSetPredefinedDoFNumber> eledofs_new = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0,eledofs,0,false));
  // replace old interior element dofs with the new created dofset
  discret_->ReplaceDofSet(nds_intvar_,eledofs_new,false);

  hdgdis_->AssignDegreesOfFreedom(0);

  // clear map cache since after every FillComplete() / AssignDegreesOfFreedom() old maps are stored in the mapstack
  output_->ClearMapCache();

  // copy old values of the state vectors phi and intphi into vectors, which are then used for the projection
  Teuchos::RCP<Epetra_Vector> phinp_old = LINALG::CreateVector(*facedofs_old,true);
  LINALG::Export(*phinp_,*phinp_old);

  Teuchos::RCP<Epetra_Vector> intphinp_old = LINALG::CreateVector(*eledofs_old,true);
  LINALG::Export(*intphinp_,*intphinp_old);

  // reset the residual, increment and sysmat to the size of the adapted new dofset
  residual_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  increment_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  neumann_loads_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  sysmat_ = Teuchos::null;
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()),27));

  // reset the state vectors
  intphinp_.reset(new Epetra_Vector(*(discret_->DofRowMap(nds_intvar_)),true));
  intphin_.reset(new Epetra_Vector(*(discret_->DofRowMap(nds_intvar_)),true));
  phinp_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  phin_.reset(new Epetra_Vector(*(discret_->DofRowMap())));

//  // end time measurement for element
//  double dtfillcomplete=Teuchos::Time::wallTime()-tcfillcomplete;
//  std::cout << "Time measurement fill complete: " << dtfillcomplete << std::endl;

//  // get cpu time
//  const double tcproject = Teuchos::Time::wallTime();

  // unpack material data
  UnpackMaterial();

  AdaptVariableVector(
      phinp_,
      phinp_old,
      intphinp_,
      intphinp_old,
      nds_var_old,
      nds_intvar_old,
      la_old);

  intphin_->Update(1.0,*intphinp_,0.0);
  phin_->Update(1.0,*phinp_,0.0);

  ProjectMaterial();

//  // end time measurement for element
//  double dtproject=Teuchos::Time::wallTime()-tcproject;
//  std::cout << "Time measurement for projection: " << dtproject << std::endl;
//
//  // get cpu time
//  const double tcmatinit = Teuchos::Time::wallTime();

  CalcMatInitial();

//  // end time measurement for element
//  double dtmatinit=Teuchos::Time::wallTime()-tcmatinit;
//  std::cout << "Time measurement calc mat initial: " << dtmatinit << std::endl;

  // end time measurement for element
  double dtadapt=Teuchos::Time::wallTime()-tcadapt;

  if(myrank_ == 0)
    std::cout << "Time measurement for adaption of element degree: " << dtadapt << std::endl;


  return;
}

/*----------------------------------------------------------------------*
 | adapt trace vector and interior variables when adapting element      |
 | degrees                                                hoermann 07/16|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::AdaptVariableVector(
    Teuchos::RCP<Epetra_Vector>                phi_new,
    Teuchos::RCP<Epetra_Vector>                phi_old,
    Teuchos::RCP<Epetra_Vector>                intphi_new,
    Teuchos::RCP<Epetra_Vector>                intphi_old,
    int                                        nds_var_old,
    int                                        nds_intvar_old,
    std::vector<DRT::Element::LocationArray>   la_old
    )
{

  // set action
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::project_field);

  // set number of dofset for the old dofsets on the parameter list to extract the correct location array
  eleparams.set<int>("nds_var_old",nds_var_old);
  eleparams.set<int>("nds_intvar_old",nds_intvar_old);

  // dof row map for adapted dofset
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(nds_intvar_);
  const Epetra_Map* dofrowmap = discret_->DofRowMap(0);


  // set old state vector on parameter list
  eleparams.set<Teuchos::RCP<Epetra_Vector> >("phi", phi_old);
  eleparams.set<Teuchos::RCP<Epetra_Vector> >("intphi", intphi_old);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector intphi_ele, phi_ele, dummyVec;

  // create location array for new and old dofsets (old ones are already filled and only copied to the location array)
  DRT::Element::LocationArray la(2*discret_->NumDofSets());
  // create location array for new dofsets
  DRT::Element::LocationArray la_temp(discret_->NumDofSets());

  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::Element *ele = discret_->lColElement(iele);

    //if the element has only ghosted nodes it will not assemble -> skip evaluation
    if(ele->HasOnlyGhostNodes(discret_->Comm().MyPID())) continue;

    // fill location array for adapted dofsets
    ele->LocationVector(*discret_,la_temp,false);

    for(int i=0; i<discret_->NumDofSets(); i++)
    {
      // copy old and new location arrays to global location array la
      la[i] = la_temp[i];
      la[discret_->NumDofSets()+i] = la_old[iele][i];
    }

    const unsigned size = la_temp[0].lm_.size();

    if (static_cast<std::size_t>(phi_ele.M()) != size)
      phi_ele.Shape(la[0].lm_.size(), 1);
    else
      memset(phi_ele.Values(),0.,size*sizeof(double));
    if (intphi_ele.M() != discret_->NumDof(nds_intvar_,ele))
      intphi_ele.Shape(discret_->NumDof(nds_intvar_,ele), 1);
    else
      memset(intphi_ele.Values(),0.,discret_->NumDof(nds_intvar_,ele)*sizeof(double));

    // call routine on elements to project values from old to new element vector
    ele->Evaluate(eleparams,*discret_,la,dummyMat,dummyMat,phi_ele,intphi_ele,dummyVec);

    // store projected values of the element on the new state vector for the interior variables
    if (ele->Owner() == discret_->Comm().MyPID())
    {
      std::vector<int> localDofs = discret_->Dof(nds_intvar_, ele);
      dsassert(localDofs.size() == static_cast<std::size_t>(intphi_ele.M()), "Internal error");
      for (unsigned int i=0; i<localDofs.size(); ++i)
        localDofs[i] = intdofrowmap->LID(localDofs[i]);
      (intphi_new)->ReplaceMyValues(localDofs.size(), intphi_ele.A(), &localDofs[0]);
    }

    // now fill the element vector into the new state vector for the trace values
    for (unsigned int i=0; i<la[0].lm_.size(); ++i)
    {
      const int lid = dofrowmap->LID(la[0].lm_[i]);

      if (lid >= 0)
        (*phi_new)[lid] = phi_ele(i);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | chooses the assembly process for matrix and rhs       hoermann 06/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::AssembleMatAndRHS()
{
  if(!DRT::INPUT::IntegralValue<int>(*params_,"SEMIIMPLICIT"))
    SCATRA::ScaTraTimIntImpl::AssembleMatAndRHS();
  else // in semi-implicit evaluation matrix does not change, thus only rhs is assembled in every step
    AssembleRHS();

  return;
}// TimIntHDG::AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | contains the assembly process only for rhs            hoermann 06/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntHDG::AssembleRHS()
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  Teuchos::RCP<LINALG::SparseOperator> dummyMat;

  // reset the residual vector
  residual_->PutScalar(0.0);

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",SCATRA::calc_mat_and_rhs);

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  eleparams.set<int>("ndsvel",nds_vel_);

  // set vector values needed by elements
  discret_->ClearState();

  // add state vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // add problem specific time-integration parameters
  AddProblemSpecificParametersAndVectors(eleparams);

  DRT::AssembleStrategy strategy( 0, 0, Teuchos::null, Teuchos::null, residual_, Teuchos::null, Teuchos::null );

  strategy.Zero();

  DRT::Element::LocationArray la(discret_->NumDofSets());

  // loop over elements
  for (int iele=0; iele<discret_->NumMyColElements(); ++iele)
  {
    DRT::Element *ele = discret_->lColElement(iele);

    //if the element has only ghosted nodes it will not assemble -> skip evaluation
    if(ele->HasOnlyGhostNodes(discret_->Comm().MyPID())) continue;

    ele->LocationVector(*discret_,la,false);

    strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

    // evaluate
    ele->Evaluate(eleparams,*discret_,la,
                     strategy.Elematrix1(),
                     strategy.Elematrix2(),
                     strategy.Elevector1(),
                     strategy.Elevector2(),
                     strategy.Elevector3());
    strategy.AssembleVector1( la[0].lm_, la[0].lmowner_ );
  }

  discret_->ClearState();

  // potential residual scaling and potential addition of Neumann terms
  ScalingAndNeumann();

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpuele;

  return;
} // TimIntHDG::AssembleRHS
