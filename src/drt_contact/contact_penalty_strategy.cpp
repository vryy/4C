/*---------------------------------------------------------------------*/
/*!
\file contact_penalty_strategy.cpp

\brief Penalty contact solving strategy: The contact constrains are enforced
       by a penalty formulation.

\level 2

\maintainer Philipp Farah, Alexander Seitz

*/
/*---------------------------------------------------------------------*/

#include "contact_penalty_strategy.H"
#include "contact_interface.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_utils.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>
#include "../drt_lib/drt_utils.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                              popp 05/09|
 *----------------------------------------------------------------------*/
CONTACT::CoPenaltyStrategy::CoPenaltyStrategy(
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
    int dim,
    Teuchos::RCP<Epetra_Comm> comm,
    double alphaf,
    int maxdof)
    : CoAbstractStrategy(Teuchos::rcp(new CONTACT::AbstractStratDataContainer()),
        DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof),
      constrnorm_(0.0),
      constrnormtan_(0.0),
      initialpenalty_(Params().get<double>("PENALTYPARAM")),
      initialpenaltytan_(Params().get<double>("PENALTYPARAMTAN"))
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoPenaltyStrategy::CoPenaltyStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
    int dim,
    Teuchos::RCP<Epetra_Comm> comm,
    double alphaf,
    int maxdof)
    : CoAbstractStrategy(data_ptr,DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof),
      constrnorm_(0.0),
      constrnormtan_(0.0),
      initialpenalty_(Params().get<double>("PENALTYPARAM")),
      initialpenaltytan_(Params().get<double>("PENALTYPARAMTAN"))
{
  // empty constructor
}


/*----------------------------------------------------------------------*
 |  save the gap-scaling kappa from reference config          popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis)
{
  // initialize the displacement field
  SetState(MORTAR::state_new_displacement, *dis);

  // kappa will be the shape function integral on the slave sides
  // (1) build the nodal information
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // interface needs to be complete
    if (!interface_[i]->Filled() && Comm().MyPID()==0)
      dserror("ERROR: FillComplete() not called on interface %", i);

    // get out of here if not participating in interface
    if (!interface_[i]->lComm()) continue;

    // do the computation of nodal shape function integral
    // (for convenience, the results will be stored in nodal gap)

    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int j=0; j<interface_[i]->SlaveColElements()->NumMyElements(); ++j)
    {
      int gid1 = interface_[i]->SlaveColElements()->GID(j);
      DRT::Element* ele1 = interface_[i]->Discret().gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CoElement* selement = dynamic_cast<CoElement*>(ele1);

      interface_[i]->IntegrateKappaPenalty(*selement);
    }

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // get nodal weighted gap
      // (this is where we stored the shape function integrals)
      double gap = cnode->CoData().Getg();

      // store kappa as the inverse of gap
      // (this removes the scaling introduced by weighting the gap!!!)
      cnode->CoData().Kappa() = 1.0/gap;

      //std::cout << "S-NODE #" << gid << " kappa=" << cnode->CoData().Kappa() << std::endl;
    }
  }
}

/*----------------------------------------------------------------------*
 | evaluate relative movement in predictor step               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::EvaluateRelMovPredict()
{
  // only for frictional contact
  if (friction_ == false) return;

  // call evaluation method of base class
  EvaluateRelMov();

  return;
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step   popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::Initialize()
{
  // (re)setup global matrices containing fc derivatives
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  linmmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // (re)setup global vector containing lagrange multipliers
  z_ = LINALG::CreateVector(*gsdofrowmap_, true);

  // (re)setup global matrix containing lagrange multiplier derivatives
  linzmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));

  return;
}

/*----------------------------------------------------------------------*
 | evaluate contact and create linear system                  popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                 Teuchos::RCP<Epetra_Vector>& feff)
{
  // in the beginning of this function, the regularized contact forces
  // in normal and tangential direction are evaluated from geometric
  // measures (gap and relative tangential velocity). Here, also active and
  // slip nodes are detected. Then, the insertion of the according stiffness
  // blocks takes place.

  bool isincontact = false;
  bool activesetchange = false;

  for (int i=0; i<(int)interface_.size(); ++i)
  {
    bool localisincontact = false;
    bool localactivesetchange = false;

    // evaluate lagrange multipliers (regularized forces) in normal direction
    // and nodal derivz matrix values, store them in nodes
    interface_[i]->AssembleRegNormalForces(localisincontact, localactivesetchange);

    // evaluate lagrange multipliers (regularized forces) in tangential direction
    INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");

    if(friction_ and soltype==INPAR::CONTACT::solution_penalty)
      interface_[i]->AssembleRegTangentForcesPenalty();

    if(friction_ and soltype==INPAR::CONTACT::solution_uzawa)
      interface_[i]->AssembleRegTangentForcesUzawa();

    isincontact = isincontact || localisincontact;
    activesetchange = activesetchange || localactivesetchange;
  }

  // broadcast contact status & active set change
  int globalcontact, globalchange = 0;
  int localcontact = isincontact;
  int localchange = activesetchange;

  Comm().SumAll(&localcontact, &globalcontact, 1);
  Comm().SumAll(&localchange, &globalchange, 1);

  if (globalcontact>=1)
  {
    isincontact_=true;
    wasincontact_=true;
  }
  else
    isincontact_=false;

  if( (Comm().MyPID()==0) && (globalchange>=1) )
    std::cout << "ACTIVE CONTACT SET HAS CHANGED..." << std::endl;

  // (re)setup active global Epetra_Maps
  // the map of global active nodes is needed for the penalty case, too.
  // this is due to the fact that we want to monitor the constraint norm
  // of the active nodes
  gactivenodes_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gactivedofs_=Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
  }

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // since we will modify the graph of kteff by adding additional
  // meshtyong stiffness entries, we have to uncomplete it
  kteff->UnComplete();

  // assemble contact quantities on all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // assemble global lagrangian multiplier vector
    interface_[i]->AssembleLM(*z_);
    // assemble global derivatives of lagrangian multipliers
    interface_[i]->AssembleLinZ(*linzmatrix_);
    // assemble global derivatives of mortar D and M matrices
    interface_[i]->AssembleLinDM(*lindmatrix_, *linmmatrix_);
  }

  // FillComplete() global matrices LinD, LinM, LinZ
  lindmatrix_->Complete(*gsmdofrowmap_, *gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofrowmap_, *gmdofrowmap_);
  linzmatrix_->Complete(*gsmdofrowmap_, *gsdofrowmap_);

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify lindmatrix_ and dmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*invtrafo_,true,*lindmatrix_,false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    lindmatrix_ = temp1;
     dmatrix_    = temp2;
  }

#ifdef CONTACTFDPENALTYTRAC
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(Params(),"FRICTION");

  // check derivatives of penalty traction
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if( IsInContact() )
    {
      if(ftype==INPAR::CONTACT::friction_coulomb )
      {
        std::cout << "LINZMATRIX" << *linzmatrix_ << std::endl;
        interface_[i]->FDCheckPenaltyTracFric();
      }
      else if (ftype==INPAR::CONTACT::friction_none)
      {
        std::cout << "-- CONTACTFDDERIVZ --------------------" << std::endl;
        interface_[i]->FDCheckPenaltyTracNor();
        std::cout << "-- CONTACTFDDERIVZ --------------------" << std::endl;
      }
      else
        dserror("Error: FD Check for this friction type not implemented!");
    }
  }
#endif

  // **********************************************************************
  // Build Contact Stiffness #1
  // **********************************************************************
  // involving contributions of derivatives of D and M:
  //  Kc,1 = delta[ 0 -M(transpose) D] * LM

  // transform if necessary
  if (ParRedist())
  {
    lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
    linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
  }

  // add to kteff
  kteff->Add(*lindmatrix_, false, 1.0-alphaf_, 1.0);
  kteff->Add(*linmmatrix_, false, 1.0-alphaf_, 1.0);

  // **********************************************************************
  // Build Contact Stiffness #2
  // **********************************************************************
  // involving contributions of derivatives of lagrange multipliers:
  //  Kc,2= [ 0 -M(transpose) D] * deltaLM

  // multiply Mortar matrices D and M with LinZ
  Teuchos::RCP<LINALG::SparseMatrix> dtilde = LINALG::MLMultiply(*dmatrix_,true,*linzmatrix_,false,false,false,true);
  Teuchos::RCP<LINALG::SparseMatrix> mtilde = LINALG::MLMultiply(*mmatrix_,true,*linzmatrix_,false,false,false,true);

  // transform if necessary
  if (ParRedist())
  {
    dtilde = MORTAR::MatrixRowTransform(dtilde,pgsdofrowmap_);
    mtilde = MORTAR::MatrixRowTransform(mtilde,pgmdofrowmap_);
  }

  // add to kteff
  kteff->Add(*dtilde, false, 1.0-alphaf_, 1.0);
  kteff->Add(*mtilde, false, -(1.0-alphaf_), 1.0);

  // **********************************************************************
  // Build RHS
  // **********************************************************************
  // feff += -alphaf * fc,n - (1-alphaf) * fc,n+1,k

  {
    // we initialize fcmdold with dold-rowmap instead of gsdofrowmap
    // (this way, possible self contact is automatically included)

    Teuchos::RCP<Epetra_Vector> fcmdold = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
    dold_->Multiply(true, *zold_, *fcmdold);
    Teuchos::RCP<Epetra_Vector> fcmdoldtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fcmdold, *fcmdoldtemp);
    feff->Update(-alphaf_, *fcmdoldtemp, 1.0);
  }

  {
    // we initialize fcmmold with mold-domainmap instead of gmdofrowmap
    // (this way, possible self contact is automatically included)

    Teuchos::RCP<Epetra_Vector> fcmmold = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
    mold_->Multiply(true, *zold_, *fcmmold);
    Teuchos::RCP<Epetra_Vector> fcmmoldtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fcmmold, *fcmmoldtemp);
    feff->Update(alphaf_, *fcmmoldtemp, 1.0);
  }

  {
    Teuchos::RCP<Epetra_Vector> fcmd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true, *z_, *fcmd);
    Teuchos::RCP<Epetra_Vector> fcmdtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fcmd, *fcmdtemp);
    feff->Update(-(1-alphaf_), *fcmdtemp, 1.0);
  }

  {
    Teuchos::RCP<Epetra_Vector> fcmm = LINALG::CreateVector(*gmdofrowmap_, true);
    mmatrix_->Multiply(true, *z_, *fcmm);
    Teuchos::RCP<Epetra_Vector> fcmmtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fcmm, *fcmmtemp);
    feff->Update(1-alphaf_, *fcmmtemp, 1.0);
  }

#ifdef CONTACTFDGAP
   // FD check of weighted gap g derivatives (non-penetr. condition)

  std::cout << "-- CONTACTFDGAP -----------------------------" << std::endl;
  interface_[0]->FDCheckGapDeriv();
  std::cout << "-- CONTACTFDGAP -----------------------------" << std::endl;

#endif // #ifdef CONTACTFDGAP

  return;
}

/*----------------------------------------------------------------------*
 | evaluate frictional contact and create linear system gitterle   10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                  Teuchos::RCP<Epetra_Vector>& feff)
{
  // this is almost the same as in the frictionless contact
  // whereas we chose the EvaluateContact routine with
  // one difference

  // check if friction should be applied
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(Params(),"FRICTION");

  // coulomb friction case
  if (ftype == INPAR::CONTACT::friction_coulomb ||
      ftype == INPAR::CONTACT::friction_stick)
  {
    EvaluateContact(kteff,feff);
  }
  else if (ftype == INPAR::CONTACT::friction_tresca)
  {
    dserror("Error in AbstractStrategy::Evaluate: Penalty Strategy for"
           " Tresca friction not yet implemented");
  }
  else
    dserror("Error in AbstractStrategy::Evaluate: Unknown friction type");

  return;
}

/*----------------------------------------------------------------------*
 | reset penalty parameter to intial value                    popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::ResetPenalty()
{
  // reset penalty parameter in strategy
  Params().set<double>("PENALTYPARAM",InitialPenalty());
  Params().set<double>("PENALTYPARAMTAN",InitialPenaltyTan());

  // reset penalty parameter in all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->IParams().set<double>("PENALTYPARAM",InitialPenalty());
    interface_[i]->IParams().set<double>("PENALTYPARAMTAN",InitialPenaltyTan());
  }

  return;
}

/*----------------------------------------------------------------------*
 | modify penalty parameter to intial value                    mhv 03/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::ModifyPenalty()
{

  // generate random number between 0.95 and 1.05
  double randnum = ((double) rand()/(double) RAND_MAX) * 0.1 + 0.95;
  double pennew = randnum * InitialPenalty();

  // modify penalty parameter in strategy
  Params().set<double>("PENALTYPARAM",pennew);
  Params().set<double>("PENALTYPARAMTAN",pennew);

  // modify penalty parameter in all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->IParams().set<double>("PENALTYPARAM",pennew);
    interface_[i]->IParams().set<double>("PENALTYPARAMTAN",pennew);
  }

  return;
}

/*----------------------------------------------------------------------*
 | intialize second, third,... Uzawa step                     popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                 Teuchos::RCP<Epetra_Vector>& feff)
{
  // remove old stiffness terms
  // (FIXME: redundant code to EvaluateContact(), expect for minus sign)

  // since we will modify the graph of kteff by adding additional
  // meshtying stiffness entries, we have to uncomplete it
  kteff->UnComplete();

  // remove contact stiffness #1 from kteff
  kteff->Add(*lindmatrix_, false, -(1.0-alphaf_), 1.0);
  kteff->Add(*linmmatrix_, false, -(1.0-alphaf_), 1.0);

  // multiply Mortar matrices D and M with LinZ
  Teuchos::RCP<LINALG::SparseMatrix> dtilde = LINALG::MLMultiply(*dmatrix_, true, *linzmatrix_, false,false,false,true);
  Teuchos::RCP<LINALG::SparseMatrix> mtilde = LINALG::MLMultiply(*mmatrix_, true, *linzmatrix_, false,false,false,true);

  // transform if necessary
  if (ParRedist())
  {
    dtilde = MORTAR::MatrixRowTransform(dtilde,pgsdofrowmap_);
    mtilde = MORTAR::MatrixRowTransform(mtilde,pgmdofrowmap_);
  }

  // remove contact stiffness #2 from kteff
  kteff->Add(*dtilde, false, -(1.0-alphaf_), 1.0);
  kteff->Add(*mtilde, false, (1.0-alphaf_), 1.0);

  // remove old force terms
  // (FIXME: redundant code to EvaluateContact(), expect for minus sign)

  Teuchos::RCP<Epetra_Vector> fcmdold = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
  dold_->Multiply(true, *zold_, *fcmdold);
  Teuchos::RCP<Epetra_Vector> fcmdoldtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcmdold, *fcmdoldtemp);
  feff->Update(alphaf_, *fcmdoldtemp, 1.0);

  Teuchos::RCP<Epetra_Vector> fcmmold = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
  mold_->Multiply(true, *zold_, *fcmmold);
  Teuchos::RCP<Epetra_Vector> fcmmoldtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcmmold, *fcmmoldtemp);
  feff->Update(-alphaf_, *fcmmoldtemp, 1.0);

  Teuchos::RCP<Epetra_Vector> fcmd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(true, *z_, *fcmd);
  Teuchos::RCP<Epetra_Vector> fcmdtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcmd, *fcmdtemp);
  feff->Update(1-alphaf_, *fcmdtemp, 1.0);

  Teuchos::RCP<Epetra_Vector> fcmm = LINALG::CreateVector(*gmdofrowmap_, true);
  mmatrix_->Multiply(true, *z_, *fcmm);
  Teuchos::RCP<Epetra_Vector> fcmmtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcmm, *fcmmtemp);
  feff->Update(-(1-alphaf_), *fcmmtemp, 1.0);

  // reset some matrices
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  linmmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // reset nodal derivZ values
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    for (int j=0;j<interface_[i]->SlaveColNodesBound()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveColNodesBound()->GID(i);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      for (int k=0; k<(int)((cnode->CoData().GetDerivZ()).size()); ++k)
        (cnode->CoData().GetDerivZ())[k].clear();
      (cnode->CoData().GetDerivZ()).resize(0);
    }
  }

  // now redo Initialize()
  Initialize();

  // and finally redo Evaluate()
  Teuchos::RCP<Epetra_Vector> nullvec = Teuchos::null;
  Evaluate(kteff,feff,nullvec);

  // complete stiffness matrix
  kteff->Complete();

  return;
}

/*----------------------------------------------------------------------*
 | evaluate L2-norm of active constraints                     popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::UpdateConstraintNorm(int uzawaiter)
{
  // initialize parameters
  double cnorm = 0.0;
  double cnormtan = 0.0;
  bool updatepenalty = false;
  bool updatepenaltytan = false;
  double ppcurr = Params().get<double>("PENALTYPARAM");
  double ppcurrtan = Params().get<double>("PENALTYPARAMTAN");

  // gactivenodes_ is undefined
  if (gactivenodes_==Teuchos::null)
  {
    constrnorm_=0;
    constrnormtan_=0;
  }

  // gactivenodes_ has no elements
  else if (gactivenodes_->NumGlobalElements()==0)
  {
    constrnorm_=0;
    constrnormtan_=0;
  }

  // gactivenodes_ has at least one element
  else
  {
    // export weighted gap vector to gactiveN-map
    Teuchos::RCP<Epetra_Vector> gact;
    if (constr_direction_==INPAR::CONTACT::constr_xyz)
    {
      gact = LINALG::CreateVector(*gactivedofs_,true);
      LINALG::Export(*g_,*gact);
    }
    else
    {
      gact = LINALG::CreateVector(*gactivenodes_,true);
      if (gact->GlobalLength()) LINALG::Export(*g_,*gact);
    }

    // compute constraint norm
    gact->Norm2(&cnorm);

    // Evaluate norm in tangential direction for frictional contact
    if (friction_)
    {
      for (int i=0; i<(int)interface_.size(); ++i)
        interface_[i]->EvaluateTangentNorm(cnormtan);

      cnormtan = sqrt(cnormtan);
    }

    //********************************************************************
    // adaptive update of penalty parameter
    // (only for Uzawa Augmented Lagrange strategy)
    //********************************************************************
    INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");

    if (soltype==INPAR::CONTACT::solution_uzawa)
    {
      // check convergence of cnorm and update penalty parameter
      // only do this for second, third, ... Uzawa iteration
      // cf. Wriggers, Computational Contact Mechanics, 2nd edition (2006), p. 340
      if ((uzawaiter >= 2) && (cnorm > 0.25 * ConstraintNorm()))
      {
        updatepenalty = true;

        // update penalty parameter in strategy
        Params().set<double>("PENALTYPARAM",10*ppcurr);

        // update penalty parameter in all interfaces
        for (int i=0; i<(int)interface_.size(); ++i)
        {
          double ippcurr = interface_[i]->IParams().get<double>("PENALTYPARAM");
          if (ippcurr != ppcurr) dserror("Something wrong with penalty parameter");
          interface_[i]->IParams().set<double>("PENALTYPARAM",10*ippcurr);
        }
        // in the case of frictional contact, the tangential penalty
        // parameter is also dated up when this is done for the normal one
        if (friction_)
        {
          updatepenaltytan = true;

          // update penalty parameter in strategy
          Params().set<double>("PENALTYPARAMTAN",10*ppcurrtan);

          // update penalty parameter in all interfaces
          for (int i=0; i<(int)interface_.size(); ++i)
          {
            double ippcurrtan = interface_[i]->IParams().get<double>("PENALTYPARAMTAN");
            if (ippcurrtan != ppcurrtan) dserror("Something wrong with penalty parameter");
            interface_[i]->IParams().set<double>("PENALTYPARAMTAN",10*ippcurrtan);
          }
        }
      }
    }
    //********************************************************************

    // update constraint norm
    constrnorm_ = cnorm;
    constrnormtan_ = cnormtan;
  }

  // output to screen
  if (Comm().MyPID()==0)
  {
    std::cout << "********************************************\n";
    std::cout << "Normal Constraint Norm: " << cnorm << "\n";
    if (friction_)
      std::cout << "Tangential Constraint Norm: " << cnormtan << "\n";
    if (updatepenalty)
      std::cout << "Updated normal penalty parameter: " << ppcurr << " -> " << Params().get<double>("PENALTYPARAM") << "\n";
    if (updatepenaltytan == true && friction_)
     std::cout << "Updated tangential penalty parameter: " << ppcurrtan << " -> " << Params().get<double>("PENALTYPARAMTAN") << "\n";
    std::cout << "********************************************\n";
  }

  return;
}

/*----------------------------------------------------------------------*
 | store Lagrange multipliers for next Uzawa step             popp 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoPenaltyStrategy::UpdateUzawaAugmentedLagrange()
{
  // store current LM into Uzawa LM
  // (note that this is also done after the last Uzawa step of one
  // time step and thus also gives the guess for the initial
  // Lagrange multiplier lambda_0 of the next time step)
  zuzawa_ = Teuchos::rcp(new Epetra_Vector(*z_));
  StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);

  return;
}
