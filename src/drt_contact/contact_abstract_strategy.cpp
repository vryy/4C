/*!----------------------------------------------------------------------
\file contact_abstract_strategy.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "Epetra_SerialComm.h"
#include "contact_abstract_strategy.H"
#include "contact_defines.H"
#include "contact_interface.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 05/09 |
 *----------------------------------------------------------------------*/
CONTACT::CoAbstractStrategy::CoAbstractStrategy(RCP<Epetra_Map> problemrowmap, Teuchos::ParameterList params,
                                                vector<RCP<CONTACT::CoInterface> > interface,
                                                int dim, RCP<Epetra_Comm> comm, double alphaf, int maxdof) :
MORTAR::StrategyBase(problemrowmap,params,dim,comm,alphaf,maxdof),
interface_(interface),
isincontact_(false),
wasincontact_(false),
wasincontactlts_(false),
isselfcontact_(false),
friction_(false),
dualquadslave3d_(false),
tsi_(false),
wear_(false)
{
  // set potential global self contact status
  // (this is TRUE if at least one contact interface is a self contact interface)
  bool selfcontact = 0;
  for (int i=0;i<(int)interface_.size();++i)
    if (interface_[i]->SelfContact()) ++selfcontact;

  if (selfcontact) isselfcontact_=true;

  // set frictional contact status
  INPAR::CONTACT::FrictionType ftype = Teuchos::getIntegralValue<INPAR::CONTACT::FrictionType>(Params(),"FRICTION");
  if (ftype != INPAR::CONTACT::friction_none)
    friction_ = true;
  
  // set wear contact status
  INPAR::CONTACT::WearType wtype = Teuchos::getIntegralValue<INPAR::CONTACT::WearType>(Params(),"WEAR");
  if (wtype != INPAR::CONTACT::wear_none)
    wear_ = true;

  // set thermo-structure-interaction with contact  
  if (DRT::Problem::Instance()->ProblemType()=="tsi")
    tsi_ = true;

  // check for infeasible self contact combinations
  if (isselfcontact_ && ftype != INPAR::CONTACT::friction_none)
    dserror("ERROR: Self contact only implemented for frictionless contact!");

  // call setup method with flag redistributed=FALSE, init=TRUE
  Setup(false,true);

  // store interface maps with parallel distribution of underlying
  // problem discretization (i.e. interface maps before parallel
  // redistribution of slave and master sides)
  if (ParRedist())
  {
    pglmdofrowmap_ = rcp(new Epetra_Map(*glmdofrowmap_));
    pgsdofrowmap_  = rcp(new Epetra_Map(*gsdofrowmap_));
    pgmdofrowmap_  = rcp(new Epetra_Map(*gmdofrowmap_));
    pgsmdofrowmap_ = rcp(new Epetra_Map(*gsmdofrowmap_));
  }
  
  // intialize storage fields for parallel redistribution
  tunbalance_.clear();
  eunbalance_.clear();

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CoAbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RedistributeContact(RCP<Epetra_Vector> dis)
{
  // get out of here if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (!ParRedist() || Comm().NumProc()==1) return;

  // decide whether redistribution should be applied or not
  double taverage = 0.0;
  double eaverage = 0;
  bool doredist = false;
  double max_balance = Params().get<double>("MAX_BALANCE");

  //**********************************************************************
  // (1) static redistribution: ONLY at time t=0 or after restart
  // (both cases can be identified via empty unbalance vectors)
  //**********************************************************************
  if (WhichParRedist()==INPAR::MORTAR::parredist_static)
  {
    // this is the first time step (t=0) or restart
    if ((int)tunbalance_.size()==0 && (int)eunbalance_.size()==0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k=0;k<(int)tunbalance_.size();++k) taverage += tunbalance_[k];
      taverage/=(int)tunbalance_.size();
      for (int k=0;k<(int)eunbalance_.size();++k) eaverage += eunbalance_[k];
      eaverage/=(int)eunbalance_.size();

      // delete balance factors of last time step
      tunbalance_.resize(0);
      eunbalance_.resize(0);

      // no redistribution
      doredist = false;
    }
  }

  //**********************************************************************
  // (2) dynamic redistribution: whenever system is out of balance
  //**********************************************************************
  else if (WhichParRedist()==INPAR::MORTAR::parredist_dynamic)
  {
    // this is the first time step (t=0) or restart
    if ((int)tunbalance_.size()==0 && (int)eunbalance_.size()==0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k=0;k<(int)tunbalance_.size();++k) taverage += tunbalance_[k];
      taverage/=(int)tunbalance_.size();
      for (int k=0;k<(int)eunbalance_.size();++k) eaverage += eunbalance_[k];
      eaverage/=(int)eunbalance_.size();
      
      // delete balance factors of last time step
      tunbalance_.resize(0);
      eunbalance_.resize(0);

      // decide on redistribution
      // -> (we allow a maximum value of the balance measure in the
      // system as defined in the input parameter MAX_BALANCE, i.e.
      // the maximum local processor workload and the minimum local
      // processor workload for mortar evaluation of all interfaces
      // may not differ by more than (MAX_BALANCE - 1.0)*100%)
      // -> (moreover, we redistribute if in the majority of iteration
      // steps of the last time step there has been an unbalance in
      // element distribution, i.e. if eaverage >= 0.5)
      if (taverage >= max_balance || eaverage >= 0.5)
        doredist = true;
    }
  }

  // print balance information to screen
  if (Comm().MyPID()==0)
  {
    cout << "**********************************************************" << endl;
    if (taverage>0)
    {
      printf("Parallel balance (time): %e (limit %e) \n",taverage,max_balance);
      printf("Parallel balance (eles): %e (limit %e) \n",eaverage,0.5);
    }
    else
      printf("Parallel balance: t=0/restart \n");
    cout << "**********************************************************" << endl;
  }

  // get out of here if simulation is still in balance
  if (!doredist) return;

  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // set old and current displacement state
  // (needed for search within redistribution)
  SetState("displacement",dis);
  SetState("olddisplacement",dis);

  // global flag for redistribution
  bool anyinterfacedone = false;

  // parallel redistribution of all interfaces
  for (int i=0; i<(int)interface_.size();++i)
  {
    // redistribute optimally among procs
    bool done = interface_[i]->Redistribute(i+1);

    // if redistribution has really been performed
    // (the previous method might have found that there
    // are no "close" slave elements and thud redistribution
    // might not be necessary ->indicated by boolean)
    if (done)
    {
      // call fill complete again
      interface_[i]->FillComplete(maxdof_);

      // print new parallel distribution
      interface_[i]->PrintParallelDistribution(i+1);

      // re-create binary search tree
      interface_[i]->CreateSearchTree();

      // set global flag to TRUE
      anyinterfacedone = true;
    }
  }

  // re-setup strategy with redistributed=TRUE, init=FALSE
  if (anyinterfacedone) Setup(true,false);

  // time measurement
  Comm().Barrier();
  double t_end = Teuchos::Time::wallTime()-t_start;
  if (Comm().MyPID()==0) cout << "\nTime for parallel redistribution.........." << t_end << " secs\n" << endl;

  return;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                                popp 08/10 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Setup(bool redistributed, bool init)
{
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  gsdofrowmap_  = Teuchos::null;
  gmdofrowmap_  = Teuchos::null;
  gsmdofrowmap_ = Teuchos::null;
  glmdofrowmap_ = Teuchos::null;
  gdisprowmap_  = Teuchos::null;
  gsnoderowmap_ = Teuchos::null;
  gactivenodes_ = Teuchos::null;
  gactivedofs_  = Teuchos::null;
  gactiven_     = Teuchos::null;
  gactivet_     = Teuchos::null;
  if (!redistributed) gndofrowmap_= Teuchos::null;

  if (friction_)
  {
    gslipnodes_ = Teuchos::null;
    gslipdofs_  = Teuchos::null;
    gslipt_     = Teuchos::null;
  }

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // build Lagrange multiplier dof map
    interface_[i]->UpdateLagMultSets(offset_if);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_, interface_[i]->LagMultDofs());
    offset_if = glmdofrowmap_->NumGlobalElements();
    if (offset_if < 0) offset_if = 0;
        
    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_, interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_, interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_, interface_[i]->MasterRowDofs());
  
    // merge active sets and slip sets of all interfaces
    // (these maps are NOT allowed to be overlapping !!!)
    interface_[i]->BuildActiveSet(init);
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(), false);
    
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the dicretization dof map)
  // (no need to rebuild this map after redistribution)
  if (!redistributed)
  {
    gndofrowmap_ = LINALG::SplitMap(*problemrowmap_, *gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);
  }

  // setup combined global slave and master dof map
  // setup global displacement dof map
  gsmdofrowmap_ = LINALG::MergeMap(*gsdofrowmap_,*gmdofrowmap_,false);
  gdisprowmap_  = LINALG::MergeMap(*gndofrowmap_,*gsmdofrowmap_,false);

  // initialize flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    IsInContact()=true;
    WasInContact()=true;
    WasInContactLastTimeStep()=true;
  }
  
  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // initialize vectors and matrices
  if (!redistributed)
  {
    // setup Lagrange multiplier vectors
    z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));

    // setup global Mortar matrices Dold and Mold
    dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
    dold_->Zero();
    dold_->Complete();
    mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_, *gsdofrowmap_);
  }

  // In the redistribution case, first check if the vectors and
  // matrices have already been defined, If yes, transform them
  // to the new redistributed maps. If not, initialize them.
  // Moreover, store redistributed quantities into nodes!!!
  else
  {
    // setup Lagrange multiplier vectors
    if (z_ == Teuchos::null)
      z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    else
    {
      RCP<Epetra_Vector> newz = rcp(new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*z_,*newz);
      z_ = newz;
    }

    if (zold_ == Teuchos::null)
      zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    else
    {
      RCP<Epetra_Vector> newzold = rcp(new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*zold_,*newzold);
      zold_ = newzold;
    }

    if (zuzawa_ == Teuchos::null)
      zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    else
    {
      RCP<Epetra_Vector> newzuzawa = rcp(new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*zuzawa_,*newzuzawa);
      zuzawa_ = newzuzawa;
    }

    // setup global Mortar matrices Dold and Mold
    if (dold_ == Teuchos::null)
    {
      dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
      dold_->Zero();
      dold_->Complete();
    }
    else
      dold_ = MORTAR::MatrixRowColTransform(dold_,gsdofrowmap_,gsdofrowmap_);

    if (mold_ == Teuchos::null)
    {
      mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_));
      mold_->Zero();
      mold_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    }
    else
      mold_ = MORTAR::MatrixRowColTransform(mold_,gsdofrowmap_,gmdofrowmap_);
  }

  // output contact stress vectors
  stressnormal_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  
  // output wear
  if (wear_) wearoutput_ = rcp(new Epetra_Vector(*gsdofrowmap_));
     
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // These matrices need to be applied to the slave displacements
  // in the cases of dual LM interpolation for tet10/hex20 meshes
  // in 3D. Here, the displacement basis functions have been modified
  // in order to assure positivity of the D matrix entries and at
  // the same time biorthogonality. Thus, to scale back the modified
  // discrete displacements \hat{d} to the nodal discrete displacements
  // {d}, we have to apply the transformation matrix T and vice versa
  // with the transformation matrix T^(-1).
  //----------------------------------------------------------------------
  INPAR::MORTAR::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
  if (shapefcn == INPAR::MORTAR::shape_dual)
    for (int i=0; i<(int)interface_.size(); ++i)
      dualquadslave3d_ += interface_[i]->Quadslave3d();

  //----------------------------------------------------------------------
  // IF SO, COMPUTE TRAFO MATRIX AND ITS INVERSE
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    trafo_    = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    invtrafo_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));

    // set of already processed nodes
    // (in order to avoid double-assembly for N interfaces)
    set<int> donebefore;

    // for all interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
      interface_[i]->AssembleTrafo(*trafo_,*invtrafo_,donebefore);

    // FillComplete() transformation matrices
    trafo_->Complete();
    invtrafo_->Complete();
  }

  return;
}

/*----------------------------------------------------------------------*
 | set current and old deformation state                     popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::SetState(const string& statename, const RCP<Epetra_Vector> vec)
{
  if( (statename=="displacement") || statename=="olddisplacement")
  {
    // set state on interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
      interface_[i]->SetState(statename, vec);
  }

  return;
}

/*----------------------------------------------------------------------*
 | update global master and slave sets (public)               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::UpdateMasterSlaveSetsGlobal()
{
  // reset global slave / master Epetra Maps
  gsnoderowmap_  = rcp(new Epetra_Map(0,0,Comm()));
  gsdofrowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  gmdofrowmap_   = rcp(new Epetra_Map(0,0,Comm()));
  glmdofrowmap_  = rcp(new Epetra_Map(0,0,Comm()));

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i=0;i<(int)interface_.size();++i)
  {
    // build Lagrange multiplier dof map
    interface_[i]->UpdateLagMultSets(offset_if);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_, interface_[i]->LagMultDofs());
    offset_if = glmdofrowmap_->NumGlobalElements();
    if (offset_if < 0) offset_if = 0;

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,interface_[i]->MasterRowDofs());
  }

  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate interface for next Newton step       popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitEvalInterface()
{
  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // initialize / reset interfaces
    interface_[i]->Initialize();

    // evaluate interfaces
    interface_[i]->Evaluate();
  }

  //**********************************************************************
  // PARALLEL REDISTRIBUTION
  //**********************************************************************
  // get out of here if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (!ParRedist() || Comm().NumProc()==1) return;

  // collect information about participation in coupling evaluation
  // and in parallel distribution of the individual interfaces
  vector<int> numloadele((int)interface_.size());
  vector<int> numcrowele((int)interface_.size());
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->CollectDistributionData(numloadele[i],numcrowele[i]);
  
  // time measurement (on each processor)
  double t_end_for_minall = Teuchos::Time::wallTime()-t_start;
  double t_end_for_maxall = t_end_for_minall;

  // restrict time measurement to procs that own at least some part
  // of the "close" slave interface section(s) on the global level,
  // i.e. restrict to procs that actually have to do some work
  int gnumloadele = 0;
  for (int i=0; i<(int)numloadele.size(); ++i)
    gnumloadele += numloadele[i];
  
  // for non-loaded procs, set time measurement to values 0.0 / 1.0e12,
  // which do not affect the maximum and minimum identification
  if (gnumloadele==0)
  {
    t_end_for_minall =  1.0e12;
    t_end_for_maxall =  0.0;
  }
  
  // store time indicator for parallel redistribution
  // (indicator is the maximum local processor time
  // divided by the minimum local processor time)
  double maxall = 0.0;
  double minall = 0.0;
  Comm().MaxAll(&t_end_for_maxall,&maxall,1);
  Comm().MinAll(&t_end_for_minall,&minall,1);
  
  // check for plausibility before storing
  if (maxall==0.0 && minall==1.0e12) tunbalance_.push_back(1.0);
  else                               tunbalance_.push_back(maxall/minall);
  
  // obtain info whether there is an unbalance in element distribution
  bool eleunbalance = false;
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // find out how many close slave elements in total
    int totrowele = 0;
    Comm().SumAll(&numcrowele[i],&totrowele,1);
    
    // find out how many procs have work on this interface
    int lhascrowele = 0;
    int ghascrowele = 0;
    if (numcrowele[i]>0) lhascrowele=1;
    Comm().SumAll(&lhascrowele,&ghascrowele,1);

    // minimum number of elements per proc
    int minele = Params().get<int>("MIN_ELEPROC");
    int numproc = Comm().NumProc();
    
    //--------------------------------------------------------------------
    // check if there is an element unbalance
    //--------------------------------------------------------------------
    // CASE 0: if minimum number of elements per proc is zero, but
    // further procs are still available and more than numproc elements
    if ( (minele == 0) && (totrowele > numproc) && (ghascrowele < numproc) )
      eleunbalance = true;
    
    // CASE 1: in total too few close slave elements but more than one
    // proc is active (otherwise, i.e. if interface small, we have no choice)
    if ( (minele > 0) && (totrowele < ghascrowele * minele) && (ghascrowele > 1) )
      eleunbalance = true;
    
    // CASE 2: in total too many close slave elements, but further procs
    // are still available for redsitribution
    if ( (minele > 0) && (totrowele >= (ghascrowele+1)*minele) && (ghascrowele < numproc) )
      eleunbalance = true;
  }
  
  // obtain global info on element unbalance
  int geleunbalance = 0;
  int leleunbalance = (int)(eleunbalance);
  Comm().SumAll(&leleunbalance,&geleunbalance,1);
  if (geleunbalance>0) eunbalance_.push_back(1);
  else                 eunbalance_.push_back(0);
  
  // debugging output
  //cout << "PROC: " << Comm().MyPID() << "\t LOADELE: " << numloadele[0] << "\t ROWELE: " << numcrowele[0]
  //     << "\t MIN: " << minall << "\t MAX: " << maxall
  //     << "\t tmin: " << t_end_for_minall << "\t tmax: " << t_end_for_maxall
  //     << "\t TUNBALANCE: " << tunbalance_[(int)tunbalance_.size()-1]
  //     << "\t EUNBALANCE: " << eunbalance_[(int)eunbalance_.size()-1] << endl;
  
  //**********************************************************************
  
  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate mortar stuff for next Newton step    popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitEvalMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  if (IsSelfContact()) UpdateMasterSlaveSetsGlobal();

  // initialize Dold and Mold if not done already
  if (dold_==null)
  {
    dold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    dold_->Zero();
    dold_->Complete();
  }
  if (mold_==null)
  {
    mold_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  }

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  mmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  g_       = LINALG::CreateVector(*gsnoderowmap_, true);
  if (friction_) jump_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  // wear
  if (wear_) wearvector_ = LINALG::CreateVector(*gsnoderowmap_, true);
  
  // matrix A for tsi problems
  if (tsi_) amatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  
  // in the case of frictional dual quad 3D, also the modified D matrices are setup
  if (friction_ && Dualquadslave3d())
  {
    // initialize Dold and Mold if not done already
    if (doldmod_==null)
    {
      doldmod_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  }
  
  // (re)setup global matrices containing fc derivatives
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  linmmatrix_ = rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    interface_[i]->AssembleDM(*dmatrix_,*mmatrix_);
    interface_[i]->AssembleG(*g_);
    
    // assemble wear vector
    if (wear_) interface_[i]->AssembleWear(*wearvector_);
    
    // assemble matrix A for tsi with contact
    if (tsi_)  interface_[i]->AssembleA(*amatrix_);
    
#ifdef CONTACTFDNORMAL
    // FD check of normal derivatives
    cout << " -- CONTACTFDNORMAL- -----------------------------------" << endl;
    interface_[i]->FDCheckNormalDeriv();
    cout << " -- CONTACTFDNORMAL- -----------------------------------" << endl;
#endif // #ifdef CONTACTFDNORMAL
  
#ifdef CONTACTFDMORTARD
    // FD check of Mortar matrix D derivatives
    cout << " -- CONTACTFDMORTARD -----------------------------------" << endl;
    dmatrix_->Complete();
    if( dmatrix_->NormOne() )
      interface_[i]->FDCheckMortarDDeriv();
    dmatrix_->UnComplete();
    cout << " -- CONTACTFDMORTARD -----------------------------------" << endl;
#endif // #ifdef CONTACTFDMORTARD
  
#ifdef CONTACTFDMORTARM
    // FD check of Mortar matrix M derivatives
    cout << " -- CONTACTFDMORTARM -----------------------------------" << endl;
    mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    if( mmatrix_->NormOne() )
        interface_[i]->FDCheckMortarMDeriv();
    mmatrix_->UnComplete();
    cout << " -- CONTACTFDMORTARM -----------------------------------" << endl;
#endif // #ifdef CONTACTFDMORTARM
  }

  // modify gap vector towards wear
  if (wear_) g_->Update(1.0,*wearvector_,1.0);

  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  
  if (tsi_) amatrix_->Complete();

  return;
}

/*----------------------------------------------------------------------*
 | evaluate reference state                                gitterle 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateReferenceState(int step,const RCP<Epetra_Vector> vec)
{
#ifndef CONTACTFORCEREFCONFIG 
  
  // only do something for frictional case
  if (!friction_) return;
  
  // only before the first time step
  if(step!=0) return;
  
  // set state and do mortar calculation
  SetState("displacement",vec);
  InitEvalInterface();
  InitEvalMortar();
  
  // store contact state to contact nodes (active or inactive)
  StoreNodalQuantities(MORTAR::StrategyBase::activeold);

  // store D and M to old ones
  StoreDM("old");

  // store nodal entries from D and M to old ones
  StoreToOld(MORTAR::StrategyBase::dm);
  
  // transform dold_ in the case of dual quadratic 3d
  if (Dualquadslave3d())
  {
    RCP<LINALG::SparseMatrix> tempold = LINALG::MLMultiply(*dold_,false,*invtrafo_,false,false,false,true);
    doldmod_ = tempold;
  }  

  // evaluate relative movement
  // needed because it is not called in the predictor of the  
  // lagrange multiplier strategy 
  EvaluateRelMov();

  // reset unbalance factors for redistribution
  // (during EvalRefState the interface has been evaluated once)
  tunbalance_.resize(0);
  eunbalance_.resize(0);

#else  
  
  if(step!=0) return;
  
  // set state and do mortar calculation
  SetState("displacement",vec);
  InitEvalInterface();
  InitEvalMortar();
  
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    dmatrix_ = temp;
  }

  // dump mortar Matrix D
  // this is always the matrix from the reference configuration,
  // also in the restart case
  dmatrix_->Dump("DREF");
  
  // only continue for frictional case
  if (!friction_) return;
  
  // store contact state to contact nodes (active or inactive)
  StoreNodalQuantities(MORTAR::StrategyBase::activeold);

  // store D and M to old ones
  StoreDM("old");

  // store nodal entries from D and M to old ones
  StoreToOld(MORTAR::StrategyBase::dm);

  // evaluate relative movement
  // needed because it is not called in the predictor of the  
  // lagrange multiplier strategy 
  EvaluateRelMov();
  
  // reset unbalance factors for redistribution
  // (during EvalRefState the interface has been evaluated once)
  tunbalance_.resize(0);
  eunbalance_.resize(0);
#endif
  
  return;
}

/*----------------------------------------------------------------------*
 | evaluate relative movement of contact bodies            gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateRelMov()
{
  // only for fricional contact
  if (!friction_) return;
  
  // transformation of slave displacement dofs
  // Dmod       ---->   D * T^(-1)
  if (Dualquadslave3d())
  {
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    dmatrixmod_ = temp;
  }
  
  // vector of slave coordinates xs
  RCP<Epetra_Vector> xsmod = rcp(new Epetra_Vector(*gsdofrowmap_));
 
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->AssembleSlaveCoord(xsmod);
 
  // ATTENTION: for EvaluateRelMov() we need the vector xsmod in
  // fully overlapping layout. Thus, export here. First, allreduce
  // slave dof row map to obtain fully overlapping slave dof map.
  RCP<Epetra_Map> fullsdofs = LINALG::AllreduceEMap(*gsdofrowmap_);
  RCP<Epetra_Vector> xsmodfull = rcp(new Epetra_Vector(*fullsdofs));
  LINALG::Export(*xsmod,*xsmodfull);
  xsmod = xsmodfull;

  // in case of 3D dual quadratic case, slave coordinates xs are modified
  if (Dualquadslave3d())
    invtrafo_->Multiply(false,*xsmod,*xsmod);
  
  // do the evaluation on the interface
  // loop over all slave row nodes on the current interface
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->EvaluateRelMov(xsmod,dmatrixmod_,doldmod_);
    interface_[i]->AssembleRelMov(*jump_);
  }
  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Evaluate(RCP<LINALG::SparseOperator>& kteff,
                                           RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis)
{
  // treat frictional and frictionless cases differently
  if (friction_) EvaluateFriction(kteff,feff);
  else           EvaluateContact(kteff,feff);

  return;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreNodalQuantities: Double active node check needed for n interfaces!");

    // get global quantity to be stored in nodes
    RCP<Epetra_Vector> vectorglobal = null;
    switch (type)
    {
      case MORTAR::StrategyBase::lmcurrent:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmold:
      {
        vectorglobal = LagrMultOld();
        break;
      }
      case MORTAR::StrategyBase::lmupdate:
      {
        vectorglobal = LagrMult();
        break;
      }
      case MORTAR::StrategyBase::lmuzawa:
      {
        vectorglobal = LagrMultUzawa();
        break;
      }
      case MORTAR::StrategyBase::activeold:
      {
        break;
      }
      case MORTAR::StrategyBase::jump:
      {
        vectorglobal = Jump();
        break;
      }
      case MORTAR::StrategyBase::wear:
      {
        break;
      }
      default:
        dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
    } // switch

    // slave dof and node map of the interface
    // columnmap for current or updated LM
    // rowmap for remaining cases 
    RCP<Epetra_Map> sdofmap, snodemap;
    if (type == MORTAR::StrategyBase::lmupdate or type == MORTAR::StrategyBase::lmcurrent)
    {  
      sdofmap = interface_[i]->SlaveColDofs();
      snodemap = interface_[i]->SlaveColNodes();
    }  
    else
    {  
      sdofmap = interface_[i]->SlaveRowDofs();
      snodemap = interface_[i]->SlaveRowNodes();
    }
      
    // export global quantity to current interface slave dof map (column or row)
    RCP<Epetra_Vector> vectorinterface = rcp(new Epetra_Vector(*sdofmap));
    
    if (vectorglobal != null) // necessary for case "activeold" and wear
      LINALG::Export(*vectorglobal, *vectorinterface);

    // loop over all slave nodes (column or row) on the current interface
    for (int j=0; j<snodemap->NumMyElements(); ++j)
    {
      int gid = snodemap->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      vector<int> locindex(dim);

      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(cnode->Dofs()[dof]);
        if (locindex[dof]<0) dserror("ERROR: StoreNodalQuantites: Did not find dof in map");

        switch(type)
        {
        case MORTAR::StrategyBase::lmcurrent:
        {
          cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmold:
        {
          cnode->MoData().lmold()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmuzawa:
        {
          cnode->MoData().lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::lmupdate:
        {
#ifndef CONTACTPSEUDO2D
          // throw a dserror if node is Active and DBC
          if (cnode->IsDbc() && cnode->Active())
            dserror("ERROR: Slave node %i is active AND carries D.B.C.s!", cnode->Id());

          // explicity set global Lag. Mult. to zero for D.B.C nodes
          if (cnode->IsDbc()) (*vectorinterface)[locindex[dof]] = 0.0;
#endif // #ifndef CONTACTPSEUDO2D

          // explicity set global Lag. Mult. to zero for inactive nodes
          // (this is what we wanted to enforce anyway before condensation)
          if (cnode->Active()==false)
            (*vectorinterface)[locindex[dof]] = 0.0;

          // store updated LM into node
          cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::activeold:
        {
          if (!friction_)
            dserror("ERROR: This should not be called for contact without friction");
          CONTACT::FriNode* frinode = static_cast<FriNode*>(cnode);
          frinode->FriData().ActiveOld() = frinode->Active();
          break;
        }
        case MORTAR::StrategyBase::jump:
        {
          if(!friction_)
            dserror("ERROR: This should not be called for contact without friction");
          FriNode* frinode = static_cast<FriNode*>(cnode);
          frinode->FriData().jump()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::wear:
        {
          if(!friction_)
            dserror("ERROR: This should not be called for contact without friction");
            // update wear only once
            if(dof==0)
            {
              FriNode* frinode = static_cast<FriNode*>(cnode);
              double wearcoeff = Params().get<double>("WEARCOEFF", 0.0);
              frinode->FriData().Wear() += wearcoeff*frinode->FriData().DeltaWear();
            }
          break;
        }
        default:
          dserror("ERROR: StoreNodalQuantities: Unknown state string variable!");
        } // switch
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Output vector of normal/tang. contact stresses        gitterle 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::OutputStresses ()
{
  // reset contact stress class variables
  stressnormal_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: OutputStresses: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      double nn[3];
      double nt1[3];
      double nt2[3];
      double lmn = 0.0;
      double lmt1 = 0.0;
      double lmt2 = 0.0;

      for (int j=0;j<3;++j)
      {
        nn[j]=cnode->MoData().n()[j];
        nt1[j]=cnode->CoData().txi()[j];
        nt2[j]=cnode->CoData().teta()[j];
        lmn +=  nn[j]* cnode->MoData().lm()[j];
        lmt1 += nt1[j]* cnode->MoData().lm()[j];
        lmt2 += nt2[j]* cnode->MoData().lm()[j];
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      vector<int> locindex(dim);

      // normal stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (stressnormal_->Map()).LID(cnode->Dofs()[dof]);
        (*stressnormal_)[locindex[dof]] = -lmn*nn[dof];
      }

      // tangential stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (stresstangential_->Map()).LID(cnode->Dofs()[dof]);
        (*stresstangential_)[locindex[dof]] = -lmt1*nt1[dof]-lmt2*nt2[dof];
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------*
|  Output de-weighted wear vector                          gitterle 10/10|
*-----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::OutputWear()
{

  // vectors
  RCP<Epetra_Vector> wear_vector = rcp(new Epetra_Vector(*gsdofrowmap_));
  RCP<Epetra_Vector> real_wear = rcp(new Epetra_Vector(*gsdofrowmap_,true));

  // solver  
  RCP<Teuchos::ParameterList> solver_params = rcp(new Teuchos::ParameterList());
  solver_params->set("solver","umfpack");
  solver_params->set("symmetric",false);
  LINALG::Solver solver(solver_params, Comm(), NULL);

  // multiply the wear with its normal direction and store in wear_vector
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // FIRST: get the wear values and the normal directions for the interface
    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* frinode = static_cast<FriNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = frinode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // nodal normal vector and wear
      double nn[3];
      double wear = 0.0;

      for (int j=0;j<3;++j)
        nn[j]=frinode->MoData().n()[j];
      wear = frinode->FriData().Wear();

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs
      vector<int> locindex(dim);

      for (int dof=0;dof<dim;++dof)
      {
         locindex[dof] = (wear_vector->Map()).LID(frinode->Dofs()[dof]);
        (*wear_vector)[locindex[dof]] = wear * nn[dof];
      }
    }
  }
 
  // approx. undo the weighting of the wear by solving
  // D * w = w~
  // dmatrix_ * real_wear = wear_
  solver.Solve(dmatrix_->EpetraMatrix(), real_wear, wear_vector, true);

  // copy the local part of real_wear into wearoutput_
  for (int i=0; i<(int)gsdofrowmap_->NumMyElements(); ++i)
  {
   int gid = gsdofrowmap_->MyGlobalElements()[i];
   double tmp = (*real_wear)[real_wear->Map().LID(gid)];
   (*wearoutput_)[wearoutput_->Map().LID(gid)] = tmp;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreDirichletStatus: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k=0;k<cnode->NumDof();++k)
      {
        int currdof = cnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid>=0)
        {
          cnode->SetDbc() = true;
          break;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDM(const string& state)
{
  //store Dold and Mold matrix in D and M
  if (state=="current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  }

  // store D and M matrix in Dold and Mold
  else if (state=="old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
    if (friction_) doldmod_ = dmatrixmod_;
  }

  // unknown conversion
  else
    dserror("ERROR: StoreDM: Unknown conversion requested!");

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal quant. to old ones   (last conv. time step)gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreToOld(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreDMToNodes: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      switch (type)
      {
        case MORTAR::StrategyBase::dm:
        {
          // store D and M entries
          cnode->StoreDMOld();
          break;
        }
        case MORTAR::StrategyBase::pentrac:
        {
          // store penalty tractions to old ones
          cnode->StoreTracOld();
          break;
        }
        default:
          dserror("ERROR: StoreDMToNodes: Unknown state string variable!");
      } // switch
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step             popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Update(int istep, RCP<Epetra_Vector> dis)
{
  // store Lagrange multipliers, D and M
  // (we need this for interpolation of the next generalized mid-point)
  // in the case of self contact, the size of z may have changed
  if (IsSelfContact()) zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));

  zold_->Update(1.0,*z_,0.0);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);
  StoreDM("old");

  // old displacements in nodes
  // (this is NOT only needed for friction but also for calculating
  // the auxiliary positions in binarytree contact search)
  SetState("olddisplacement",dis);
  
#ifdef MORTARGMSH1  
  VisualizeGmsh(istep);
#endif // #ifdef MORTARGMSH1

  // double-check if active set is really converged
  // (necessary e.g. for monolithic FSI with Lagrange multiplier contact,
  // because usually active set convergence check has been integrated into
  // structure Newton scheme, but now the monolithic FSI Newton scheme decides)
  // not for thermo-structure-interaction because of partitioned scheme
  // so long
  std::string probtype = DRT::Problem::Instance()->ProblemType();
  if (!tsi_ and (!ActiveSetConverged() || !ActiveSetSemiSmoothConverged()))
    dserror("ERROR: Active set not fully converged!");
  
  // reset active set status for next time step
  ResetActiveSet();

  // update flag for global contact status of last time step
  if (gactivenodes_->NumGlobalElements())
  {
    WasInContact()=true;
    WasInContactLastTimeStep()=true;
  }
  else
  {
    WasInContact()=false;
    WasInContactLastTimeStep()=false;
  }

  //----------------------------------------friction: store history values
  // in the case of frictional contact we have to store several
  // information and quantities at the end of a time step (converged
  // state) which is needed in the next time step as history
  // information / quantities. 
  if (friction_)
  {
    // store contact state to contact nodes (active or inactive)
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);
  
    // store nodal entries of D and M to old ones
    StoreToOld(MORTAR::StrategyBase::dm);

    // store nodal entries form penalty contact tractions to old ones
    StoreToOld(MORTAR::StrategyBase::pentrac);
  }
  
  // update wear
  if(wear_)
    StoreNodalQuantities(MORTAR::StrategyBase::wear);
  
  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoWriteRestart(RCP<Epetra_Vector>& activetoggle,
                                                 RCP<Epetra_Vector>& sliptoggle,
                                                 RCP<Epetra_Vector>& weightedwear)
{
  // initalize
  activetoggle = rcp(new Epetra_Vector(*gsnoderowmap_));
  if (friction_)
    sliptoggle = rcp(new Epetra_Vector(*gsnoderowmap_));
  if (wear_)
    weightedwear = rcp(new Epetra_Vector(*gsnoderowmap_));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: DoWriteRestart: Double active node check needed for n interfaces!");
        
    // loop over all slave nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = static_cast<CoNode*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      // set value active / inactive in toggle vector
      if (cnode->Active()) (*activetoggle)[dof]=1;
      
      // set value slip / stick in the toggle vector 
      if (friction_)
      {
        CONTACT::FriNode* frinode = static_cast<CONTACT::FriNode*>(cnode);
        if (frinode->FriData().Slip()) (*sliptoggle)[dof]=1;
        if (wear_)  (*weightedwear)[dof] = frinode->FriData().Wear();
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoReadRestart(IO::DiscretizationReader& reader,
                                                RCP<Epetra_Vector> dis)
{
  // set restart displacement state
  SetState("displacement", dis);
  SetState("olddisplacement", dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  InitEvalInterface();
  InitEvalMortar();

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_
    RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    dmatrix_ = temp;
  }

  // read restart information on actice set and slip set
  RCP<Epetra_Vector> activetoggle =rcp(new Epetra_Vector(*gsnoderowmap_));
  reader.ReadVector(activetoggle,"activetoggle");
  
  // friction
  RCP<Epetra_Vector> sliptoggle;
  RCP<Epetra_Vector> weightedwear;
  if(friction_)
  {  
    sliptoggle =rcp(new Epetra_Vector(*gsnoderowmap_));
    reader.ReadVector(sliptoggle,"sliptoggle");
  }

  // wear
  if (wear_)
  {
    weightedwear = rcp(new Epetra_Vector(*gsnoderowmap_));
    reader.ReadVector(weightedwear, "weightedwear");
  }
  
  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: DoReadRestart: Double active node check needed for n interfaces!");
        
    // loop over all slave nodes on the current interface
    for (int j=0; j<(interface_[i]->SlaveRowNodes())->NumMyElements(); ++j)
    {
      int gid = (interface_[i]->SlaveRowNodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof]==1)
      {
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = static_cast<CoNode*>(node);

        // set value active / inactive in cnode
        cnode->Active()=true;
        
        if (friction_)
        {
          // set value stick / slip in cnode
          // set wear value
          if ((*sliptoggle)[dof]==1) static_cast<CONTACT::FriNode*>(cnode)->FriData().Slip()=true;
          if (wear_) static_cast<CONTACT::FriNode*>(cnode)->FriData().Wear() = (*weightedwear)[dof];
        }
      }
    }
  }

  // read restart information on Lagrange multipliers
  z_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = rcp(new Epetra_Vector(*gsdofrowmap_));
  reader.ReadVector(LagrMult(),"lagrmultold");
  reader.ReadVector(LagrMultOld(),"lagrmultold");

  // store restart information on Lagrange multipliers into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // only for Augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  INPAR::CONTACT::SolvingStrategy st = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");
  if (st == INPAR::CONTACT::solution_auglag)
  {
    zuzawa_ = rcp(new Epetra_Vector(*gsdofrowmap_));
    reader.ReadVector(LagrMultUzawa(),"lagrmultold");
    StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);
  }

  // store restart Mortar quantities
  StoreDM("old");
  
  if (friction_)
  {
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);
    StoreToOld(MORTAR::StrategyBase::dm);
  }

  // (re)setup active global Epetra_Maps
  gactivenodes_ = null;
  gactivedofs_ = null;
  gactiven_ = null;
  gactivet_ = null;
  gslipnodes_ = null;
  gslipdofs_ = null;
  gslipt_ = null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(), false);
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }  
  }

  // update flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    IsInContact()=true;
    WasInContact()=true;
    WasInContactLastTimeStep()=true;
  }
  
  // evaluate relative movement (jump)
  // needed because it is not called in the predictor of the  
  // lagrange multiplier strategy 
  EvaluateRelMov();
  
  // reset unbalance factors for redistribution
  // (during restart the interface has been evaluated once)
  tunbalance_.resize(0);
  eunbalance_.resize(0);

  return;
}


/*----------------------------------------------------------------------*
 |  Compute interface forces (for debugging only)             popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InterfaceForces(bool output)
{
  // Note that we ALWAYS use a TR-like approach to compute the contact
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}

  // check chosen output option
  INPAR::CONTACT::EmOutputType emtype =
    Teuchos::getIntegralValue<INPAR::CONTACT::EmOutputType>(Params(),"EMOUTPUT");
  
  // get out of here if no output wanted
  if (emtype==INPAR::CONTACT::output_none) return;
  
  // compute two subvectors of fc each via Lagrange multipliers z_n+1, z_n
  RCP<Epetra_Vector> fcslavetemp = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> fcmastertemp = rcp(new Epetra_Vector(mmatrix_->DomainMap()));
  RCP<Epetra_Vector> fcslavetempend = rcp(new Epetra_Vector(dold_->RowMap()));
  RCP<Epetra_Vector> fcmastertempend = rcp(new Epetra_Vector(mold_->DomainMap()));

  // for self contact, slave and master sets may have changed,
  // thus we have to export z and zold to new D and M dimensions
  if (IsSelfContact())
  {
    RCP<Epetra_Vector> zexp = rcp(new Epetra_Vector(dmatrix_->RowMap()));
    RCP<Epetra_Vector> zoldexp = rcp(new Epetra_Vector(dold_->RowMap()));
    if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_,*zexp);
    if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    dmatrix_->Multiply(true,*zexp,*fcslavetemp);
    mmatrix_->Multiply(true,*zexp,*fcmastertemp);
    dold_->Multiply(true,*zoldexp,*fcslavetempend);
    mold_->Multiply(true,*zoldexp,*fcmastertempend);
  }
  // if there is no self contact everything is ok
  else
  {
    dmatrix_->Multiply(true, *z_, *fcslavetemp);
    mmatrix_->Multiply(true, *z_, *fcmastertemp);
    dold_->Multiply(true, *zold_, *fcslavetempend);
    mold_->Multiply(true, *zold_, *fcmastertempend);
  }

  // export the contact forces to full dof layout
  RCP<Epetra_Vector> fcslave = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmaster = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcslaveend = rcp(new Epetra_Vector(*problemrowmap_));
  RCP<Epetra_Vector> fcmasterend = rcp(new Epetra_Vector(*problemrowmap_));
  LINALG::Export(*fcslavetemp, *fcslave);
  LINALG::Export(*fcmastertemp, *fcmaster);
  LINALG::Export(*fcslavetempend, *fcslaveend);
  LINALG::Export(*fcmastertempend, *fcmasterend);

  // contact forces and moments
  vector<double> gfcs(3);
  vector<double> ggfcs(3);
  vector<double> gfcm(3);
  vector<double> ggfcm(3);
  vector<double> gmcs(3);
  vector<double> ggmcs(3);
  vector<double> gmcm(3);
  vector<double> ggmcm(3);

  vector<double> gmcsnew(3);
  vector<double> ggmcsnew(3);
  vector<double> gmcmnew(3);
  vector<double> ggmcmnew(3);

  // weighted gap vector
  RCP<Epetra_Vector> gapslave  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> gapmaster = rcp(new Epetra_Vector(mmatrix_->DomainMap()));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      vector<double> nodeforce(3);
      vector<double> position(3);

      // forces and positions
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      vector<double> nodemoment(3);
      nodemoment[0] = position[1]*nodeforce[2]-position[2]*nodeforce[1];
      nodemoment[1] = position[2]*nodeforce[0]-position[0]*nodeforce[2];
      nodemoment[2] = position[0]*nodeforce[1]-position[1]*nodeforce[0];
      for (int d=0;d<3;++d)
        gmcs[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      for (int d=0;d<Dim();++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapslave,posnode,lm,lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j=0;j<interface_[i]->MasterRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      vector<double> nodeforce(3);
      vector<double> position(3);

      // forces and positions
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcmastertemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: ContactForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      vector<double> nodemoment(3);
      nodemoment[0] = position[1]*nodeforce[2]-position[2]*nodeforce[1];
      nodemoment[1] = position[2]*nodeforce[0]-position[0]*nodeforce[2];
      nodemoment[2] = position[0]*nodeforce[1]-position[1]*nodeforce[0];
      for (int d=0;d<3;++d)
        gmcm[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      vector<int> lm(Dim());
      vector<int> lmowner(Dim());
      for (int d=0;d<Dim();++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapmaster,posnode,lm,lmowner);
    }
  }

  // weighted gap
  RCP<Epetra_Vector> gapslavefinal  = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  RCP<Epetra_Vector> gapmasterfinal = rcp(new Epetra_Vector(mmatrix_->RowMap()));
  dmatrix_->Multiply(false,*gapslave,*gapslavefinal);
  mmatrix_->Multiply(false,*gapmaster,*gapmasterfinal);
  RCP<Epetra_Vector> gapfinal = rcp(new Epetra_Vector(dmatrix_->RowMap()));
  gapfinal->Update(1.0,*gapslavefinal,0.0);
  gapfinal->Update(-1.0,*gapmasterfinal,1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      vector<double> lm(3);
      vector<double> nodegaps(3);
      vector<double> nodegapm(3);

      // LMs and gaps
      for (int d=0;d<Dim();++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid<0) dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = cnode->MoData().lm()[d];
      }

      // moments
      vector<double> nodemoments(3);
      vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1]*lm[2]-nodegaps[2]*lm[1];
      nodemoments[1] = nodegaps[2]*lm[0]-nodegaps[0]*lm[2];
      nodemoments[2] = nodegaps[0]*lm[1]-nodegaps[1]*lm[0];
      nodemomentm[0] = nodegapm[1]*lm[2]-nodegapm[2]*lm[1];
      nodemomentm[1] = nodegapm[2]*lm[0]-nodegapm[0]*lm[2];
      nodemomentm[2] = nodegapm[0]*lm[1]-nodegapm[1]*lm[0];
      for (int d=0;d<3;++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }
    }
  }

  // summing up over all processors
  for (int i=0;i<3;++i)
  {
    Comm().SumAll(&gfcs[i],&ggfcs[i],1);
    Comm().SumAll(&gfcm[i],&ggfcm[i],1);
    Comm().SumAll(&gmcs[i],&ggmcs[i],1);
    Comm().SumAll(&gmcm[i],&ggmcm[i],1);
    Comm().SumAll(&gmcsnew[i],&ggmcsnew[i],1);
    Comm().SumAll(&gmcmnew[i],&ggmcmnew[i],1);
  }

  // print interface results to file
  if (emtype == INPAR::CONTACT::output_file ||
      emtype == INPAR::CONTACT::output_both)
  {
    // do this at end of time step only (output==true)!
    // processor 0 does all the work
    if (output && Comm().MyPID()==0)
    {
      FILE* MyFile = NULL;
      MyFile = fopen("o/scilab_output/OutputInterface.txt", "at+");
      
      if (MyFile)
      {
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggfcs[i]);
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggfcm[i]);
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggmcs[i]);
        for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", ggmcm[i]);
        //for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsfgh[i]);
        //for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsmgh[i]);
        fprintf(MyFile, "\n");
        fclose(MyFile);
      }
      else
        dserror("ERROR: File for writing meshtying forces could not be opened.");
    } 
  }
  
  // print interface results to screen
  if (emtype == INPAR::CONTACT::output_screen ||
      emtype == INPAR::CONTACT::output_both)
  {
    // do this during Newton steps only (output==false)!
    // processor 0 does all the work
    if (!output && Comm().MyPID()==0)
    {
      double snorm = sqrt(ggfcs[0]*ggfcs[0]+ggfcs[1]*ggfcs[1]+ggfcs[2]*ggfcs[2]);
      double mnorm = sqrt(ggfcm[0]*ggfcm[0]+ggfcm[1]*ggfcm[1]+ggfcm[2]*ggfcm[2]);
      printf("Slave Contact Force:   % e  % e  % e \tNorm: % e\n",ggfcs[0],ggfcs[1],ggfcs[2], snorm);
      printf("Master Contact Force:  % e  % e  % e \tNorm: % e\n",ggfcm[0],ggfcm[1],ggfcm[2], mnorm);
      printf("Slave Contact Moment:  % e  % e  % e\n",ggmcs[0],ggmcs[1],ggmcs[2]);
      //printf("Slave Contact Moment:  % e  % e  % e\n",ggmcsnew[0],ggmcsnew[1],ggmcsnew[2]);
      printf("Master Contact Moment: % e  % e  % e\n",ggmcm[0],ggmcm[1],ggmcm[2]);
      //printf("Master Contact Moment: % e  % e  % e\n",ggmcmnew[0],ggmcmnew[1],ggmcmnew[2]);
      fflush(stdout);
    }
  }
  
  return;
}

/*----------------------------------------------------------------------*
 | evaluate contact forces w.r.t. reference configuration     mgit 07/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ForceRefConfig()
{
  // multiply current D matrix with current LM
  RCP<Epetra_Vector> forcecurr = rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*z_,*forcecurr);
  
  // read mortar matrix D of reference configuration 
  RCP<LINALG::SparseMatrix> dref = rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  std::string inputname = "DREF";
  dref->Load(Comm(),inputname);
  dref->Complete();
  
  // LM in reference / current configuration
  RCP<Epetra_Vector> zref  = rcp(new Epetra_Vector(*gsdofrowmap_));
  RCP<Epetra_Vector> zcurr  = rcp(new Epetra_Vector(*z_));

  // solve with default solver
  LINALG::Solver solver(Comm());
  solver.Solve(dref->EpetraOperator(),zref,forcecurr,true);

  // store reference LM into global vector and nodes
  z_ = zref;
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
  
  // print message
  if (Comm().MyPID()==0)
    cout << "\n**** First output is w.r.t. the reference configuration! ****" << endl;

  // print active set
  PrintActiveSet();

  // restore current LM into global vector and nodes
  z_ = zcurr;
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "--------------------------------- CONTACT::CoAbstractStrategy\n"
       << "Contact interfaces: " << (int)interface_.size() << endl
       << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    cout << *(interface_[i]);
  }
  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | print active set information                               popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::PrintActiveSet()
{
  // output message
  Comm().Barrier();
  if (Comm().MyPID()==0)
  {
    printf("\nActive contact set--------------------------------------------------------------\n");
    fflush(stdout);
  }

  //**********************************************************************
  // detailed active set output
  //**********************************************************************

#ifdef CONTACTASOUTPUT
  // create storage for local and global data
  vector<int>    lnid, gnid;
  vector<double> llmn, glmn;
  vector<double> lgap, ggap;

  // introduce integer variable status
  // (0=inactive, 1=active, 2=slip, 3=stick)
  // (this is necessary as all data will be written by proc 0, but
  // the knowledge of the above status ONLY exists on the owner
  // processor of the respective node. Thus this information must
  // also be communicated to proc 0 in addition to the actual data!)
  vector<int>    lsta, gsta;

  // some more storage for local and global friction data
  vector<double> llmt, glmt;
  vector<double> ljtx, gjtx;
  vector<double> ljte, gjte;
  vector<double> lwear, gwear;

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      // gid of current node
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);

      //--------------------------------------------------------------------
      // FRICTIONLESS CASE
      //--------------------------------------------------------------------
      if (!friction_)
      {
        // cast to CoNode
        CoNode* cnode = static_cast<CoNode*>(node);

        // compute weighted gap
        double wgap = (*g_)[g_->Map().LID(gid)];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k=0;k<3;++k)
          nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active()) lsta.push_back(1);
        else                 lsta.push_back(0);
      }

      //--------------------------------------------------------------------
      // FRICTIONAL CASE
      //--------------------------------------------------------------------
      else
      {
        // cast to CoNode and FriNode
        CoNode* cnode = static_cast<CoNode*>(node);
        FriNode* frinode = static_cast<FriNode*>(cnode);

        // compute weighted gap
        double wgap = (*g_)[g_->Map().LID(gid)];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k=0;k<3;++k)
          nz += frinode->MoData().n()[k] * frinode->MoData().lm()[k];

        // compute tangential parts of Lagrange multiplier and jumps and wear
        double txiz = 0.0;
        double tetaz = 0.0;
        double jumptxi = 0.0;
        double jumpteta = 0.0;
        double wear = 0.0;
        
        for (int k=0;k<Dim();++k)
        {
          txiz += frinode->CoData().txi()[k] * frinode->MoData().lm()[k];
          tetaz += frinode->CoData().teta()[k] * frinode->MoData().lm()[k];
          jumptxi += frinode->CoData().txi()[k] * frinode->FriData().jump()[k];
          jumpteta += frinode->CoData().teta()[k] * frinode->FriData().jump()[k];
        }

        // total tangential component
        double tz = sqrt(txiz*txiz+tetaz*tetaz);

        // check for dimensions
        if (Dim()==2 && abs(jumpteta)>0.0001)
          dserror("Error: Jumpteta should be zero for 2D");
        
        // compute weighted wear
        if (wear_) wear = frinode->FriData().Wear();

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);
        llmt.push_back(tz);
        ljtx.push_back(jumptxi);
        ljte.push_back(jumpteta);
        lwear.push_back(wear);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active())
        {
          if (frinode->FriData().Slip()) lsta.push_back(2);
          else                           lsta.push_back(3);
        }
        else
        {
          lsta.push_back(0);
        }
      }
    }
  }

  // we want to gather data from on all procs
  vector<int> allproc(Comm().NumProc());
  for (int i=0; i<Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid,gnid,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(llmn,glmn,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(lgap,ggap,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<int>(lsta,gsta,(int)allproc.size(),&allproc[0],Comm());

  // communicate some more data to proc 0 for friction
  if (friction_)
  {
    LINALG::Gather<double>(llmt,glmt,(int)allproc.size(),&allproc[0],Comm());
    LINALG::Gather<double>(ljtx,gjtx,(int)allproc.size(),&allproc[0],Comm());
    LINALG::Gather<double>(ljte,gjte,(int)allproc.size(),&allproc[0],Comm());
    LINALG::Gather<double>(lwear,gwear,(int)allproc.size(),&allproc[0],Comm());
 }

  // output is solely done by proc 0
  if (Comm().MyPID()==0)
  {
    //--------------------------------------------------------------------
    // FRICTIONLESS CASE
    //--------------------------------------------------------------------
    if (!friction_)
    {
      // loop over all nodes
      for (int k=0;k<(int)gnid.size();++k)
      {
        // print nodes of inactive set *************************************
        if (gsta[k]==0)
        {
          printf("INACTIVE: %d \t wgap: % e \t lm: % e \n",gnid[k],ggap[k],glmn[k]);
          fflush(stdout);
        }

        // print nodes of active set ***************************************
        else if (gsta[k]==1)
        {
          printf("ACTIVE:   %d \t wgap: % e \t lm: % e \n",gnid[k],ggap[k],glmn[k]);
          fflush(stdout);
        }

        // invalid status **************************************************
        else dserror("ERROR: Invalid node status %i for frictionless case",gsta[k]);
      }
    }

    //--------------------------------------------------------------------
    // FRICTIONAL CASE
    //--------------------------------------------------------------------
    else
    {
      // loop over all nodes
      for (int k=0;k<(int)gnid.size();++k)
      {
        // print nodes of slip set **************************************
        if (gsta[k]==2)
        {
          printf("SLIP:  %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",gnid[k],glmn[k],glmt[k],gjtx[k],gjte[k],gwear[k]);
          fflush(stdout);
        }

        // print nodes of stick set *************************************
        else if (gsta[k]==3)
        {
          printf("STICK: %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",gnid[k],glmn[k],glmt[k],gjtx[k],gjte[k],gwear[k]);
          fflush(stdout);
        }

        // print nodes of inactive set *************************************
        else if (gsta[k]==0)
        {
          // do nothing
        }

        // invalid status **************************************************
        else dserror("ERROR: Invalid node status %i for frictional case",gsta[k]);
      }
    }
  }

#else
  //**********************************************************************
  // reduced active set output
  //**********************************************************************

  // counters
  int activenodes = 0;
  int gactivenodes = 0;
  int inactivenodes = 0;
  int ginactivenodes = 0;
  int slipnodes = 0;
  int gslipnodes = 0;

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);

      // increase active counters
      CoNode* cnode = static_cast<CoNode*>(node);
      if (cnode->Active()) activenodes   += 1;
      else                 inactivenodes += 1;

      // increase friction counters
      if (friction_)
      {
        FriNode* frinode = static_cast<FriNode*>(cnode);
        if (cnode->Active() && frinode->FriData().Slip()) slipnodes += 1;
      }
    }
  }

  // sum among all processors
  Comm().SumAll(&activenodes,&gactivenodes,1);
  Comm().SumAll(&inactivenodes,&ginactivenodes,1);
  Comm().SumAll(&slipnodes,&gslipnodes,1);

  // print active set information
  if (Comm().MyPID()==0)
  {
    if (friction_)
    {
      cout << BLUE2_LIGHT  << "Total     SLIP nodes:\t" << gslipnodes << END_COLOR << endl;
      cout << BLUE2_LIGHT  << "Total    STICK nodes:\t" << gactivenodes-gslipnodes << END_COLOR << endl;
      cout << RED_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes << END_COLOR << endl;
    }
    else
    {
      cout << BLUE2_LIGHT <<  "Total   ACTIVE nodes:\t" << gactivenodes << END_COLOR << endl;
      cout << RED_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes << END_COLOR << endl;
    }
  }
#endif // #ifdef CONTACTASOUTPUT

  // output line
  Comm().Barrier();
  if (Comm().MyPID()==0)
  {
    printf("--------------------------------------------------------------------------------\n\n");
    fflush(stdout);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of contact segments with gmsh                popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::VisualizeGmsh(const int step, const int iter)
{
  // visualization with gmsh
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->VisualizeGmsh(step, iter);
}

#endif // CCADISCRET
