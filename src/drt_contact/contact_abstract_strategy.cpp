/*!----------------------------------------------------------------------
\file contact_abstract_strategy.cpp

\maintainer Philipp Farah, Alexander Seitz

*-----------------------------------------------------------------------*/
#include "Epetra_SerialComm.h"
#include "contact_abstract_strategy.H"
#include "contact_defines.H"
#include "contact_interface.H"
#include "../drt_contact_aug/contact_augmented_interface.H"
#include "friction_node.H"

#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_utils.H"

#include "../drt_inpar/inpar_contact.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_sparsematrix.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 05/09 |
 *----------------------------------------------------------------------*/
CONTACT::CoAbstractStrategy::CoAbstractStrategy(
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
    int dim,
    Teuchos::RCP<Epetra_Comm> comm,
    double alphaf,
    int maxdof) :
MORTAR::StrategyBase(DofRowMap,NodeRowMap,params,dim,comm,alphaf,maxdof),
interface_(interface),
step_(0),
iter_(0),
iterls_(-1),
isincontact_(false),
wasincontact_(false),
wasincontactlts_(false),
isselfcontact_(false),
friction_(false),
regularized_(false),
dualquadslave3d_(false),
stype_(DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(params,"STRATEGY")),
constr_direction_(DRT::INPUT::IntegralValue<INPAR::CONTACT::ConstraintDirection>(params,"CONSTRAINT_DIRECTIONS"))
{
  // set potential global self contact status
  // (this is TRUE if at least one contact interface is a self contact interface)
  bool selfcontact = false;
  for (int i = 0; i < (int) interface_.size(); ++i)
    if (interface_[i]->SelfContact())
      selfcontact = true;

  if (selfcontact)
    isselfcontact_ = true;

  INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::FrictionType>(Params(), "FRICTION");

  // set frictional contact status
  if (ftype != INPAR::CONTACT::friction_none)
    friction_ = true;

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::Regularization>(Params(), "CONTACT_REGULARIZATION")
      != INPAR::CONTACT::reg_none)
    regularized_ = true;

  // call setup method with flag redistributed=FALSE, init=TRUE
  Setup(false, true);

  // store interface maps with parallel distribution of underlying
  // problem discretization (i.e. interface maps before parallel
  // redistribution of slave and master sides)
  if (ParRedist())
  {
    pglmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*glmdofrowmap_));
    pgsdofrowmap_  = Teuchos::rcp(new Epetra_Map(*gsdofrowmap_));
    pgmdofrowmap_  = Teuchos::rcp(new Epetra_Map(*gmdofrowmap_));
    pgsmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gsmdofrowmap_));
  }

  // intialize storage fields for parallel redistribution
  tunbalance_.clear();
  eunbalance_.clear();

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator <<(std::ostream& os,
    const CONTACT::CoAbstractStrategy& strategy)
{
  strategy.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::RedistributeContact(
    Teuchos::RCP<Epetra_Vector> dis)
{
  // get out of here if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (!ParRedist() || Comm().NumProc() == 1)
    return;

  // decide whether redistribution should be applied or not
  double taverage = 0.0;
  double eaverage = 0;
  bool doredist = false;
  double max_balance = Params().get<double>("MAX_BALANCE");

  //**********************************************************************
  // (1) static redistribution: ONLY at time t=0 or after restart
  // (both cases can be identified via empty unbalance vectors)
  //**********************************************************************
  if (WhichParRedist() == INPAR::MORTAR::parredist_static)
  {
    // this is the first time step (t=0) or restart
    if ((int) tunbalance_.size() == 0 && (int) eunbalance_.size() == 0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k = 0; k < (int) tunbalance_.size(); ++k)
        taverage += tunbalance_[k];
      taverage /= (int) tunbalance_.size();
      for (int k = 0; k < (int) eunbalance_.size(); ++k)
        eaverage += eunbalance_[k];
      eaverage /= (int) eunbalance_.size();

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
  else if (WhichParRedist() == INPAR::MORTAR::parredist_dynamic)
  {
    // this is the first time step (t=0) or restart
    if ((int) tunbalance_.size() == 0 && (int) eunbalance_.size() == 0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k = 0; k < (int) tunbalance_.size(); ++k)
        taverage += tunbalance_[k];
      taverage /= (int) tunbalance_.size();
      for (int k = 0; k < (int) eunbalance_.size(); ++k)
        eaverage += eunbalance_[k];
      eaverage /= (int) eunbalance_.size();

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
  if (Comm().MyPID() == 0)
  {
    std::cout << "**********************************************************"
        << std::endl;
    if (taverage > 0)
    {
      printf("Parallel balance (time): %e (limit %e) \n", taverage,
          max_balance);
      printf("Parallel balance (eles): %e (limit %e) \n", eaverage, 0.5);
    }
    else
      printf("Parallel balance: t=0/restart \n");
    std::cout << "**********************************************************"
        << std::endl;
  }

  // get out of here if simulation is still in balance
  if (!doredist)
    return;

  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // set old and current displacement state
  // (needed for search within redistribution)
  SetState("displacement", dis);
  SetState("olddisplacement", dis);

  // global flag for redistribution
  bool anyinterfacedone = false;

  // parallel redistribution of all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // redistribute optimally among procs
    bool done = interface_[i]->Redistribute(i + 1);

    // if redistribution has really been performed
    // (the previous method might have found that there
    // are no "close" slave elements and thus redistribution
    // might not be necessary ->indicated by boolean)
    if (done)
    {
      // call fill complete again
      interface_[i]->FillComplete(maxdof_);

      // print new parallel distribution
      interface_[i]->PrintParallelDistribution(i + 1);

      // re-create binary search tree
      interface_[i]->CreateSearchTree();

      // set global flag to TRUE
      anyinterfacedone = true;
    }
  }

  // re-setup strategy with redistributed=TRUE, init=FALSE
  if (anyinterfacedone)
    Setup(true, false);

  // time measurement
  Comm().Barrier();
  double t_end = Teuchos::Time::wallTime() - t_start;
  if (Comm().MyPID() == 0)
    std::cout << "\nTime for parallel redistribution.........." << t_end
        << " secs\n" << std::endl;

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
  gsdofrowmap_ = Teuchos::null;
  gmdofrowmap_ = Teuchos::null;
  gsmdofrowmap_ = Teuchos::null;
  glmdofrowmap_ = Teuchos::null;
  gdisprowmap_ = Teuchos::null;
  gsnoderowmap_ = Teuchos::null;
  gmnoderowmap_ = Teuchos::null;
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  if (!redistributed)
    gndofrowmap_ = Teuchos::null;
  if (init)
    initial_elecolmap_.clear();
  initial_elecolmap_.resize(0);

  if (friction_)
  {
    gslipnodes_ = Teuchos::null;
    gslipdofs_ = Teuchos::null;
    gslipt_ = Teuchos::null;
  }

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // build Lagrange multiplier dof map
    interface_[i]->UpdateLagMultSets(offset_if);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_,
        interface_[i]->LagMultDofs());
    offset_if = glmdofrowmap_->NumGlobalElements();
    if (offset_if < 0)
      offset_if = 0;

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,
        interface_[i]->SlaveRowNodes());
    gmnoderowmap_ = LINALG::MergeMap(gmnoderowmap_,
        interface_[i]->MasterRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,
        interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,
        interface_[i]->MasterRowDofs());

    // merge active sets and slip sets of all interfaces
    // (these maps are NOT allowed to be overlapping !!!)
    interface_[i]->BuildActiveSet(init);
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,
        interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(),
        false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(),
        false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(),
        false);

    // store initial element col map for binning strategy
    initial_elecolmap_.push_back(
        Teuchos::rcp(
            new Epetra_Map(*interface_[i]->Discret().ElementColMap())));

    // ****************************************************
    // friction
    // ****************************************************
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(),
          false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(),
          false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the discretization dof map)
  // (no need to rebuild this map after redistribution)
  if (!redistributed)
  {
    gndofrowmap_ = LINALG::SplitMap(*(ProblemDofs()),
        *gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *gmdofrowmap_);
  }

  // setup combined global slave and master dof map
  // setup global displacement dof map
  gsmdofrowmap_ = LINALG::MergeMap(*gsdofrowmap_, *gmdofrowmap_, false);
  gdisprowmap_ = LINALG::MergeMap(*gndofrowmap_, *gsmdofrowmap_, false);

  // initialize flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
  }

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // initialize vectors and matrices
  if (!redistributed)
  {
    // setup Lagrange multiplier vectors
    z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

    // setup global Mortar matrices Dold and Mold
    dold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,1,true,false));
    dold_->Zero();
    dold_->Complete();
    mold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,1,true,false));
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
    {
      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> newz = Teuchos::rcp(
          new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*z_, *newz);
      z_ = newz;
    }

    if (zincr_ == Teuchos::null)
    {
      zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> newzincr = Teuchos::rcp(
          new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*zincr_, *newzincr);
      zincr_ = newzincr;
    }

    if (zold_ == Teuchos::null)
      zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    else
    {
      Teuchos::RCP<Epetra_Vector> newzold = Teuchos::rcp(
          new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*zold_, *newzold);
      zold_ = newzold;
    }

    if (zuzawa_ == Teuchos::null)
      zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    else
    {
      Teuchos::RCP<Epetra_Vector> newzuzawa = Teuchos::rcp(
          new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*zuzawa_, *newzuzawa);
      zuzawa_ = newzuzawa;
    }

    // setup global Mortar matrices Dold and Mold
    if (dold_ == Teuchos::null)
    {
      dold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,1,true,false));
      dold_->Zero();
      dold_->Complete();
    }
    else
      dold_ = MORTAR::MatrixRowColTransform(dold_, gsdofrowmap_, gsdofrowmap_);

    if (mold_ == Teuchos::null)
    {
      mold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,1,true,false));
      mold_->Zero();
      mold_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    }
    else
      mold_ = MORTAR::MatrixRowColTransform(mold_, gsdofrowmap_, gmdofrowmap_);
  }

  // output contact stress vectors
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

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
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");
  if ((shapefcn == INPAR::MORTAR::shape_dual
      || shapefcn == INPAR::MORTAR::shape_petrovgalerkin) && Dim() == 3)
    for (int i = 0; i < (int) interface_.size(); ++i)
      dualquadslave3d_ += interface_[i]->Quadslave();

  //----------------------------------------------------------------------
  // IF SO, COMPUTE TRAFO MATRIX AND ITS INVERSE
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
    invtrafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));

    // set of already processed nodes
    // (in order to avoid double-assembly for N interfaces)
    std::set<int> donebefore;

    // for all interfaces
    for (int i = 0; i < (int) interface_.size(); ++i)
      interface_[i]->AssembleTrafo(*trafo_, *invtrafo_, donebefore);

    // FillComplete() transformation matrices
    trafo_->Complete();
    invtrafo_->Complete();
  }

  // transform modified old D-matrix in case of friction
  // (ony necessary after parallel redistribution)
  if (redistributed && friction_ && Dualquadslave3d())
  {
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,1,true,false));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    else
      doldmod_ = MORTAR::MatrixRowColTransform(doldmod_, gsdofrowmap_, gsdofrowmap_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | global evaluation method called from time integrator      popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::ApplyForceStiffCmt(
    Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt,
    Teuchos::RCP<Epetra_Vector>& f, const int step, const int iter,
    bool predictor)
{
  // update step and iteration counters
  step_ = step;
  iter_ = iter;

  /******************************************/
  /*     VERSION WITH TIME MEASUREMENT      */
  /******************************************/

#ifdef CONTACTTIME
  // mortar initialization and evaluation
  Comm().Barrier();
  const double t_start1 = Teuchos::Time::wallTime();
  SetState("displacement",dis);
  Comm().Barrier();
  const double t_end1 = Teuchos::Time::wallTime()-t_start1;
  if (Comm().MyPID()==0) std::cout << "    -->SetState:\t" << t_end1 << " seconds\n";

  Comm().Barrier();
  const double t_start2 = Teuchos::Time::wallTime();
  InitMortar();
  Comm().Barrier();
  const double t_end2 = Teuchos::Time::wallTime()-t_start2;
  if (Comm().MyPID()==0) std::cout << "    -->MortarInit  :\t" << t_end2 << " seconds\n";

  Comm().Barrier();
  const double t_start3 = Teuchos::Time::wallTime();
  InitEvalInterface();
  Comm().Barrier();
  const double t_end3 = Teuchos::Time::wallTime()-t_start3;
  if (Comm().MyPID()==0)
  {
    std::cout << "    -->Interface:\t" << t_end3 << " seconds";
    if ((int)tunbalance_.size()==0 && (int)eunbalance_.size()==0) std::cout << "\n";
    else std::cout << " (BALANCE: " << tunbalance_.back() << " " << eunbalance_.back() << ")\n";
  }

  Comm().Barrier();
  const double t_start4 = Teuchos::Time::wallTime();
  AssembleMortar();
  Comm().Barrier();
  const double t_end4 = Teuchos::Time::wallTime()-t_start4;
  if (Comm().MyPID()==0) std::cout << "    -->AssembleMortar  :\t" << t_end4 << " seconds\n";

  // evaluate relative movement for friction
  Comm().Barrier();
  const double t_start5 = Teuchos::Time::wallTime();
  if (predictor) EvaluateRelMovPredict();
  else EvaluateRelMov();
  Comm().Barrier();
  const double t_end5 = Teuchos::Time::wallTime()-t_start5;
  if (Comm().MyPID()==0) std::cout << "    -->RelMov  :\t" << t_end5 << " seconds\n";

  // update active set
  Comm().Barrier();
  const double t_start6 = Teuchos::Time::wallTime();
  if (!predictor) UpdateActiveSetSemiSmooth();
  Comm().Barrier();
  const double t_end6 = Teuchos::Time::wallTime()-t_start6;
  if (Comm().MyPID()==0) std::cout << "    -->ActivSet:\t" << t_end6 << " seconds\n";

  // apply contact forces and stiffness
  Comm().Barrier();
  const double t_start7 = Teuchos::Time::wallTime();
  Initialize();
  Comm().Barrier();
  const double t_end7 = Teuchos::Time::wallTime()-t_start7;
  if (Comm().MyPID()==0) std::cout << "    -->Initial :\t" << t_end7 << " seconds\n";

  Comm().Barrier();
  const double t_start8 = Teuchos::Time::wallTime();
  Evaluate(kt,f,dis);
  Comm().Barrier();
  const double t_end8 = Teuchos::Time::wallTime()-t_start8;
  if (Comm().MyPID()==0) std::cout << "    -->Evaluate:\t" << t_end8 << " seconds\n";

  Comm().Barrier();
  const double t_start9 = Teuchos::Time::wallTime();
  InterfaceForces();
  Comm().Barrier();
  const double t_end9 = Teuchos::Time::wallTime()-t_start9;
  if (Comm().MyPID()==0) std::cout << "    -->IfForces:\t" << t_end9 << " seconds\n";

#else

  /******************************************/
  /*   VERSION WITHOUT TIME MEASUREMENT     */
  /******************************************/

  // mortar initialization and evaluation
  SetState("displacement", dis);

  //---------------------------------------------------------------
  // For selfcontact the master/slave sets are updated within the -
  // contact search, see SelfBinaryTree.                          -
  // Therefore, we have to initialize the mortar matrices after   -
  // interface evaluations.                                       -
  //---------------------------------------------------------------
  if (IsSelfContact())
  {
    InitEvalInterface(); // evaluate mortar terms (integrate...)
    InitMortar(); // initialize mortar matrices and vectors
    AssembleMortar(); // assemble mortar terms into global matrices
  }
  else
  {
    InitMortar(); // initialize mortar matrices and vectors
    InitEvalInterface(); // evaluate mortar terms (integrate...)
    AssembleMortar(); // assemble mortar terms into global matrices
  }

  // evaluate relative movement for friction
  if (predictor)
    EvaluateRelMovPredict();
  else
    EvaluateRelMov();

  // update active set
  if (!predictor and (stype_ != INPAR::CONTACT::solution_augmented))
    UpdateActiveSetSemiSmooth();

  // apply contact forces and stiffness
  Initialize();         // init lin-matrices
  Evaluate(kt, f, dis); // assemble lin. matrices, condensation ...

  // Evaluate structural and constraint rhs. This is also necessary, if the rhs did not
  // change during the predictor step, but a redistribution was executed!
  UpdateStructuralRHS(f);
  UpdateStructuralStiff(kt);
  EvalConstrRHS();

  //only for debugging:
  InterfaceForces();

#endif // #ifdef CONTACTTIME
  return;
}
/*----------------------------------------------------------------------*
 | set current and old deformation state                     popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::SetState(const std::string& statename,
    const Teuchos::RCP<Epetra_Vector> vec)
{
  if ((statename == "displacement") || statename == "olddisplacement")
  {
    // set state on interfaces
    for (int i = 0; i < (int) interface_.size(); ++i)
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
  gsnoderowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  gsdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  gmdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
  glmdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // build Lagrange multiplier dof map
    interface_[i]->UpdateLagMultSets(offset_if);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_,
        interface_[i]->LagMultDofs());
    offset_if = glmdofrowmap_->NumGlobalElements();
    if (offset_if < 0)
      offset_if = 0;

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ = LINALG::MergeMap(gsnoderowmap_,
        interface_[i]->SlaveRowNodes());
    gsdofrowmap_ = LINALG::MergeMap(gsdofrowmap_,
        interface_[i]->SlaveRowDofs());
    gmdofrowmap_ = LINALG::MergeMap(gmdofrowmap_,
        interface_[i]->MasterRowDofs());
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate mean. vel. for bin size                         farah 11/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::CalcMeanVelforBinning(
    Teuchos::RCP<Epetra_Vector> vel)
{
  ivel_.clear();
  ivel_.resize(0);

  // for dynamic problems
  if (alphaf_ != 0.0)
  {
    // create vector of interface velocities
    for (int i = 0; i < (int) interface_.size(); ++i)
    {
      // interface node map
      Teuchos::RCP<Epetra_Vector> velidofs = Teuchos::rcp(
          new Epetra_Vector(*(interface_[i]->Discret().DofRowMap())));
      LINALG::Export(*vel, *velidofs);

      double mean = 0.0;

      int err = velidofs->MeanValue(&mean);
      if (err)
        dserror("error for meanvalue calculation");
      mean = abs(mean);

      ivel_.push_back(mean);
    }
  }
  // static problems
  else
  {
    // TODO: should be fixed for static problems!
    dserror("Binning Strategy is only recommended for dynamic problems! Please use a different Parallel Strategy!");
  }
  return;
}

/*----------------------------------------------------------------------*
 | Prepare binning                                           farah 11/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitBinStrategyforTimestep(
    Teuchos::RCP<Epetra_Vector> vel)
{
  // calc mean value of contact interfaces
  CalcMeanVelforBinning(vel);

  // create bins and perform ghosting for each interface
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // init interface
    interface_[i]->Initialize();

    // call binning strategy
    interface_[i]->BinningStrategy(initial_elecolmap_[i], ivel_[i]);

    //build new search tree or do nothing for bruteforce
    if (DRT::INPUT::IntegralValue<INPAR::MORTAR::SearchAlgorithm>(Params(),
        "SEARCH_ALGORITHM") == INPAR::MORTAR::search_binarytree)
      interface_[i]->CreateSearchTree();
    else if (DRT::INPUT::IntegralValue<INPAR::MORTAR::SearchAlgorithm>(Params(),
        "SEARCH_ALGORITHM") != INPAR::MORTAR::search_bfele)
      dserror("ERROR: Invalid search algorithm");
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

  // get type of parallel strategy
  INPAR::MORTAR::ParallelStrategy strat = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::ParallelStrategy>(Params(), "PARALLEL_STRATEGY");
  INPAR::MORTAR::RedundantStorage redundant = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::RedundantStorage>(Params(), "REDUNDANT_STORAGE");
  // check solving strategy
  if (stype_ == INPAR::CONTACT::solution_augmented && strat!=INPAR::MORTAR::ghosting_redundant)
      dserror("ERROR: Augmented Lagrange strategy supports only redundant ghosting of "
          "master and slave.");

  // Evaluation for all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // initialize / reset interfaces
    interface_[i]->Initialize();

    //store required integration time
    inttime_ += interface_[i]->Inttime();

    /************************************************************
     *  Round Robin loop with evaluations after each iteration  *
     ************************************************************/
    if (strat == INPAR::MORTAR::roundrobinevaluate)
    {
      // check redundant input
      if (redundant != INPAR::MORTAR::redundant_none)
        dserror("Round-Robin-Loop only for none-redundant storage of interface!");

      // this contains the evaluation as well as the rr loop
      interface_[i]->RoundRobinEvaluate();
    }
    /************************************************************
     *  Round Robin loop only for ghosting                      *
     ************************************************************/
    else if (strat == INPAR::MORTAR::roundrobinghost)
    {
      // check redundant input
      if (redundant != INPAR::MORTAR::redundant_none)
        dserror("Round-Robin-Loop only for none-redundant storage of interface!");

      // first perform rrloop to detect the required ghosting
      interface_[i]->RoundRobinDetectGhosting();

      // second step --> evaluate
      interface_[i]->Evaluate(0, step_, iter_);
    }
    /************************************************************
     *  Creating bins and ghost all mele within adjacent bins   *
     ************************************************************/
    else if (strat == INPAR::MORTAR::binningstrategy)
    {
      // check redundant input
      if (redundant != INPAR::MORTAR::redundant_none)
        dserror("Binning strategy only for none-redundant storage of interface!");

      // required master elements are already ghosted (preparestepcontact) !!!
      // call evaluation
      interface_[i]->Evaluate(0, step_, iter_);
    }
    /************************************************************
     *  Fully redundant ghosting of master side                 *
     ************************************************************/
    else //std. evaluation for redundant ghosting
    {
      if (stype_ == INPAR::CONTACT::solution_augmented)
      {
        Teuchos::RCP<CONTACT::AugmentedInterface> augInterface =
            Teuchos::rcp_dynamic_cast<CONTACT::AugmentedInterface>(interface_[i]);
        // evaluate averaged weighted gap
        augInterface->Evaluate(0,step_,iter_);
        // Calculate weighted gap
        augInterface->WGap();
        // Do the augmented active set decision
        BuildGlobalAugActiveSet(i);
        // Calculate averaged weighted gap linearization for all active nodes
        augInterface->AWGapLin();
        // evaluate remaining entities and linearization
        augInterface->RedEvaluate();
      }
      else
        interface_[i]->Evaluate(0,step_,iter_);
    }
  } // end interface loop

  //**********************************************************************
  // PARALLEL REDISTRIBUTION
  //**********************************************************************
  // don't do this if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (ParRedist() && Comm().NumProc() > 1)
  {
    // collect information about participation in coupling evaluation
    // and in parallel distribution of the individual interfaces
    std::vector<int> numloadele((int) interface_.size());
    std::vector<int> numcrowele((int) interface_.size());
    for (int i = 0; i < (int) interface_.size(); ++i)
      interface_[i]->CollectDistributionData(numloadele[i], numcrowele[i]);

    // time measurement (on each processor)
    double t_end_for_minall = Teuchos::Time::wallTime() - t_start;
    double t_end_for_maxall = t_end_for_minall;

    // restrict time measurement to procs that own at least some part
    // of the "close" slave interface section(s) on the global level,
    // i.e. restrict to procs that actually have to do some work
    int gnumloadele = 0;
    for (int i = 0; i < (int) numloadele.size(); ++i)
      gnumloadele += numloadele[i];

    // for non-loaded procs, set time measurement to values 0.0 / 1.0e12,
    // which do not affect the maximum and minimum identification
    if (gnumloadele == 0)
    {
      t_end_for_minall = 1.0e12;
      t_end_for_maxall = 0.0;
    }

    // store time indicator for parallel redistribution
    // (indicator is the maximum local processor time
    // divided by the minimum local processor time)
    double maxall = 0.0;
    double minall = 0.0;
    Comm().MaxAll(&t_end_for_maxall, &maxall, 1);
    Comm().MinAll(&t_end_for_minall, &minall, 1);

    // check for plausibility before storing
    if (maxall == 0.0 && minall == 1.0e12)
      tunbalance_.push_back(1.0);
    else
      tunbalance_.push_back(maxall / minall);

    // obtain info whether there is an unbalance in element distribution
    bool eleunbalance = false;
    for (int i = 0; i < (int) interface_.size(); ++i)
    {
      // find out how many close slave elements in total
      int totrowele = 0;
      Comm().SumAll(&numcrowele[i], &totrowele, 1);

      // find out how many procs have work on this interface
      int lhascrowele = 0;
      int ghascrowele = 0;
      if (numcrowele[i] > 0)
        lhascrowele = 1;
      Comm().SumAll(&lhascrowele, &ghascrowele, 1);

      // minimum number of elements per proc
      int minele = Params().get<int>("MIN_ELEPROC");
      int numproc = Comm().NumProc();

      //--------------------------------------------------------------------
      // check if there is an element unbalance
      //--------------------------------------------------------------------
      // CASE 0: if minimum number of elements per proc is zero, but
      // further procs are still available and more than numproc elements
      if ((minele == 0) && (totrowele > numproc) && (ghascrowele < numproc))
        eleunbalance = true;

      // CASE 1: in total too few close slave elements but more than one
      // proc is active (otherwise, i.e. if interface small, we have no choice)
      if ((minele > 0) && (totrowele < ghascrowele * minele)
          && (ghascrowele > 1))
        eleunbalance = true;

      // CASE 2: in total too many close slave elements, but further procs
      // are still available for redsitribution
      if ((minele > 0) && (totrowele >= (ghascrowele + 1) * minele)
          && (ghascrowele < numproc))
        eleunbalance = true;
    }

    // obtain global info on element unbalance
    int geleunbalance = 0;
    int leleunbalance = (int) (eleunbalance);
    Comm().SumAll(&leleunbalance, &geleunbalance, 1);
    if (geleunbalance > 0)
      eunbalance_.push_back(1);
    else
      eunbalance_.push_back(0);

    // debugging output
    //std::cout << "PROC: " << Comm().MyPID() << "\t LOADELE: " << numloadele[0] << "\t ROWELE: " << numcrowele[0]
    //     << "\t MIN: " << minall << "\t MAX: " << maxall
    //     << "\t tmin: " << t_end_for_minall << "\t tmax: " << t_end_for_maxall
    //     << "\t TUNBALANCE: " << tunbalance_[(int)tunbalance_.size()-1]
    //     << "\t EUNBALANCE: " << eunbalance_[(int)eunbalance_.size()-1] << std::endl;
  }
  //**********************************************************************

  //**********************************************************************
  // OVERVIEW OF PARALLEL MORTAR COUPLING STATUS
  //**********************************************************************
#ifdef CONTACTSTATUS
  // total numbers per processor
  std::vector<int> smpairs(1);
  std::vector<int> smintpairs(1);
  std::vector<int> intcells(1);

  // add numbers of all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    smpairs[0] += interface_[i]->SlaveMasterPairs();
    smintpairs[0] += interface_[i]->SlaveMasterIntPairs();
    intcells[0] += interface_[i]->IntegrationCells();
  }

  // vector containing all proc ids
  const int numproc = Comm().NumProc();
  std::vector<int> allproc(numproc);
  for (int i=0; i<numproc; ++i) allproc[i] = i;

  // global numbers
  std::vector<int> gsmpairs, gsmintpairs, gintcells;
  LINALG::Gather<int>(smpairs,gsmpairs,numproc,&allproc[0],Comm());
  LINALG::Gather<int>(smintpairs,gsmintpairs,numproc,&allproc[0],Comm());
  LINALG::Gather<int>(intcells,gintcells,numproc,&allproc[0],Comm());

  // output to screen
  if (Comm().MyPID()==0)
  {
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << std::setw(10) << "proc ID" << std::setw(16) << "# s/m pairs"
    << std::setw(16) << "# s/m intpairs" << std::setw(16) << "# intcells" << std::endl;
    for (int i=0; i<numproc; ++i)
    {
      std::cout << std::setw(10) << i << std::setw(16) << gsmpairs[i] << std::setw(16)
      << gsmintpairs[i] << std::setw(16) << gintcells[i] << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
  }
#endif // #ifdef CONTACTSTATUS
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 | initialize mortar stuff for next Newton step               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::InitMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  if (IsSelfContact())
    UpdateMasterSlaveSetsGlobal();

  // check solving strategy
  if (stype_ == INPAR::CONTACT::solution_augmented) return;

  // initialize Dold and Mold if not done already
  if (dold_ == Teuchos::null)
  {
    dold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
    dold_->Zero();
    dold_->Complete();
  }
  if (mold_ == Teuchos::null)
  {
    mold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 100));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_, *gsdofrowmap_);
  }

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
  mmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 100));

  if (constr_direction_==INPAR::CONTACT::constr_xyz)
    g_ = LINALG::CreateVector(*gsdofrowmap_, true);
  else if (constr_direction_==INPAR::CONTACT::constr_ntt)
    g_ = LINALG::CreateVector(*gsnoderowmap_, true);
  else
    dserror("unknown contact constraint direction");

  // in the case of frictional dual quad 3D, also the modified D matrices are setup
  if (friction_ && Dualquadslave3d())
  {
    // initialize Dold and Mold if not done already
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble mortar stuff for next Newton step                 popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::AssembleMortar()
{
  // check solving strategy
  if (stype_ == INPAR::CONTACT::solution_augmented) return;

  // for all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    interface_[i]->AssembleDM(*dmatrix_, *mmatrix_);
    interface_[i]->AssembleG(*g_);

#ifdef CONTACTFDNORMAL
    // FD check of normal derivatives
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
//    interface_[i]->FDCheckNormalDeriv();
    interface_[i]->FDCheckNormalCPPDeriv();
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDNORMAL
#ifdef CONTACTFDMORTARD
    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
    dmatrix_->Complete();
    if( dmatrix_->NormOne() )
    interface_[i]->FDCheckMortarDDeriv();
    dmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARD
#ifdef CONTACTFDMORTARM
    // FD check of Mortar matrix M derivatives
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
    mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    if( mmatrix_->NormOne() )
    interface_[i]->FDCheckMortarMDeriv();
    mmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARM
  }

  // FillComplete() global Mortar matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate reference state                               gitterle 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateReferenceState(const Teuchos::RCP<Epetra_Vector> vec)
{
  // flag for initualization of contact with nodal gaps
  bool initcontactbygap = DRT::INPUT::IntegralValue<int>(Params(),"INITCONTACTBYGAP");

  // only do something for frictional case
  // or for initialization of initial contact set with nodal gap
  if (!friction_ and !initcontactbygap) return;

  // set state and do mortar calculation
  SetState("displacement", vec);
  InitMortar();
  InitEvalInterface();
  AssembleMortar();

  // (1) GAP INITIALIZATION CASE
  // initialize init contact with nodal gap
  if (initcontactbygap)
  {
    // merge interface maps to global maps
    for (int i = 0; i < (int) interface_.size(); ++i)
    {
      // merge active sets and slip sets of all interfaces
      // (these maps are NOT allowed to be overlapping !!!)
      interface_[i]->BuildActiveSet(true);
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

    // initialize flags for global contact status
    if (gactivenodes_->NumGlobalElements())
    {
      isincontact_ = true;
      wasincontact_ = true;
      wasincontactlts_ = true;
    }

    // error if no nodes are initialized to active
    if (gactivenodes_->NumGlobalElements() == 0)
      dserror("ERROR: No active nodes: Choose bigger value for INITCONTACTGAPVALUE!");
  }

  // (2) FRICTIONAL CONTACT CASE
  // do some friction stuff
  if (friction_)
  {
    // store contact state to contact nodes (active or inactive)
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);

    // store D and M to old ones
    StoreDM("old");

    // store nodal entries from D and M to old ones
    StoreToOld(MORTAR::StrategyBase::dm);

    // transform dold_ in the case of dual quadratic 3d
    if (Dualquadslave3d())
    {
      Teuchos::RCP<LINALG::SparseMatrix> tempold = LINALG::MLMultiply(*dold_,
        false, *invtrafo_, false, false, false, true);
      doldmod_ = tempold;
    }

    // evaluate relative movement
    // needed because it is not called in the predictor of the
    // lagrange multiplier strategy
    EvaluateRelMov();
  }

  // reset unbalance factors for redistribution
  // (since the interface has been evaluated once above)
  tunbalance_.resize(0);
  eunbalance_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate relative movement of contact bodies           gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::EvaluateRelMov()
{
  // only for fricional contact
  if (!friction_)
    return;

  // transformation of slave displacement dofs
  // Dmod       ---->   D * T^(-1)
  if (Dualquadslave3d())
  {
    Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,
        false, *invtrafo_, false, false, false, true);
    dmatrixmod_ = temp;
  }

  // vector of slave coordinates xs
  Teuchos::RCP<Epetra_Vector> xsmod = Teuchos::rcp(
      new Epetra_Vector(*gsdofrowmap_));

  for (int i = 0; i < (int) interface_.size(); ++i)
    interface_[i]->AssembleSlaveCoord(xsmod);

  // ATTENTION: for EvaluateRelMov() we need the vector xsmod in
  // fully overlapping layout. Thus, export here. First, allreduce
  // slave dof row map to obtain fully overlapping slave dof map.
  Teuchos::RCP<Epetra_Map> fullsdofs = LINALG::AllreduceEMap(*gsdofrowmap_);
  Teuchos::RCP<Epetra_Vector> xsmodfull = Teuchos::rcp(
      new Epetra_Vector(*fullsdofs));
  LINALG::Export(*xsmod, *xsmodfull);
  xsmod = xsmodfull;

  // in case of 3D dual quadratic case, slave coordinates xs are modified
  if (Dualquadslave3d())
    invtrafo_->Multiply(false, *xsmod, *xsmod);

  // evaluation of obj. invariant slip increment
  // do the evaluation on the interface
  // loop over all slave row nodes on the current interface
  if (DRT::INPUT::IntegralValue<int>(Params(), "GP_SLIP_INCR") == false)
    for (int i = 0; i < (int) interface_.size(); ++i)
      interface_[i]->EvaluateRelMov(xsmod, dmatrixmod_, doldmod_);

  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Evaluate(
    Teuchos::RCP<LINALG::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
{
  // treat frictional and frictionless cases differently
  if (friction_)
    EvaluateFriction(kteff, feff);
  else
    EvaluateContact(kteff, feff);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate matrix of normals (for VelocityUpdate)            popp 10/11|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::CoAbstractStrategy::EvaluateNormals(
    Teuchos::RCP<Epetra_Vector> dis)
{
  // set displacement state and evaluate nodal normals
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    interface_[i]->SetState("displacement", dis);
    interface_[i]->EvaluateNodalNormals();
  }

  // create empty global matrix
  // (rectangular: rows=snodes, cols=sdofs)
  Teuchos::RCP<LINALG::SparseMatrix> normals = Teuchos::rcp(
      new LINALG::SparseMatrix(*gsnoderowmap_, 3));

  // assemble nodal normals
  for (int i = 0; i < (int) interface_.size(); ++i)
    interface_[i]->AssembleNormals(*normals);

  // complete global matrix
  // (rectangular: rows=snodes, cols=sdofs)
  normals->Complete(*gsdofrowmap_, *gsnoderowmap_);

  return normals;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreNodalQuantities(
    MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // get global quantity to be stored in nodes
    Teuchos::RCP<Epetra_Vector> vectorglobal = Teuchos::null;

    // start type switch
    switch (type)
    {
    case MORTAR::StrategyBase::lmold:
    {
      vectorglobal = LagrMultOld();
      break;
    }
    case MORTAR::StrategyBase::lmcurrent:
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
    case MORTAR::StrategyBase::slipold:
    {
      break;
    }
    default:
      dserror("ERROR: StoreNodalQuantities: Unknown state std::string variable!");
      break;
    } // switch

    // slave dof and node map of the interface
    // columnmap for current or updated LM
    // rowmap for remaining cases
    Teuchos::RCP<Epetra_Map> sdofmap, snodemap;
    if (type == MORTAR::StrategyBase::lmupdate  or
        type == MORTAR::StrategyBase::lmcurrent)
    {
      sdofmap  = interface_[i]->SlaveColDofs();
      snodemap = interface_[i]->SlaveColNodes();
    }
    else
    {
      sdofmap  = interface_[i]->SlaveRowDofs();
      snodemap = interface_[i]->SlaveRowNodes();
    }

    // export global quantity to current interface slave dof map (column or row)
    Teuchos::RCP<Epetra_Vector> vectorinterface = Teuchos::null;
    vectorinterface = Teuchos::rcp(new Epetra_Vector(*sdofmap));
    if (vectorglobal != Teuchos::null) // necessary for case "activeold" and wear
      LINALG::Export(*vectorglobal, *vectorinterface);

    // loop over all slave nodes (column or row) on the current interface
    for (int j = 0; j < snodemap->NumMyElements(); ++j)
    {
      int gid = snodemap->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      const int dim = Dim();
      const int numdof = cnode->NumDof();
      if (dim != numdof)
        dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      std::vector<int> locindex(dim);

      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(cnode->Dofs()[dof]);
        if (locindex[dof] < 0)
          dserror("ERROR: StoreNodalQuantites: Did not find dof in map");

        switch (type)
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
            dserror("ERROR: Slave node %i is active AND carries D.B.C.s!",cnode->Id());
#endif // #ifndef CONTACTPSEUDO2D

          // store updated LM into node
          cnode->MoData().lm()[dof] = (*vectorinterface)[locindex[dof]];
          break;
        }
        case MORTAR::StrategyBase::activeold:
        {
          cnode->CoData().ActiveOld() = cnode->Active();
          break;
        }
        case MORTAR::StrategyBase::slipold:
        {
          if (!friction_)
            dserror("ERROR: Slip just for friction problems!");

          FriNode* fnode = dynamic_cast<FriNode*>(cnode);
          fnode->FriData().SlipOld() = fnode->FriData().Slip();
          break;
        }
        default:
          dserror("ERROR: StoreNodalQuantities: Unknown state std::string variable!");
          break;
        } // switch
      }
    } // end slave loop
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Output vector of normal/tang. contact stresses        gitterle 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::OutputStresses()
{
  // reset contact stress class variables
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: OutputStresses: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim != numdof)
        dserror("ERROR: Inconsisteny Dim <-> NumDof");

      double nn[3];
      double nt1[3];
      double nt2[3];
      double lmn = 0.0;
      double lmt1 = 0.0;
      double lmt2 = 0.0;

      for (int j = 0; j < 3; ++j)
      {
        nn[j] = cnode->MoData().n()[j];
        nt1[j] = cnode->CoData().txi()[j];
        nt2[j] = cnode->CoData().teta()[j];
        lmn += nn[j] * cnode->MoData().lm()[j];
        lmt1 += nt1[j] * cnode->MoData().lm()[j];
        lmt2 += nt2[j] * cnode->MoData().lm()[j];
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(dim);

      // normal stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (stressnormal_->Map()).LID(cnode->Dofs()[dof]);
        if (DRT::INPUT::IntegralValue<int>(Params(), "LM_NODAL_SCALE") == false
            || cnode->MoData().GetScale() == 0.)
          (*stressnormal_)[locindex[dof]] = -lmn * nn[dof];
        else
          (*stressnormal_)[locindex[dof]] = -lmn * nn[dof]
              / cnode->MoData().GetScale();
      }

      // tangential stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (stresstangential_->Map()).LID(cnode->Dofs()[dof]);
        if (DRT::INPUT::IntegralValue<int>(Params(), "LM_NODAL_SCALE") == false
            || cnode->MoData().GetScale() == 0.)
          (*stresstangential_)[locindex[dof]] = -lmt1 * nt1[dof]
              - lmt2 * nt2[dof];
        else
          (*stresstangential_)[locindex[dof]] = -lmt1 * nt1[dof]
              - lmt2 * nt2[dof] / cnode->MoData().GetScale();
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDirichletStatus(
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: StoreDirichletStatus: Double active node check needed for n interfaces!");

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // check if this node's dofs are in dbcmap
      for (int k = 0; k < cnode->NumDof(); ++k)
      {
        int currdof = cnode->Dofs()[k];
        int lid = (dbcmaps->CondMap())->LID(currdof);

        // store dbc status if found
        if (lid >= 0 && cnode->DbcDofs()[k] == false)
          cnode->SetDbc() = true;

        // check compatibility of contact symmetry condition and displacement dirichlet conditions
        if (lid<0 && cnode->DbcDofs()[k]==true)
        {
          std::cout << "node " << cnode->Id() << " at: " << cnode->X()[0] << " " << cnode->X()[1] << " " << cnode->X()[2] << std::endl;
          std::cout << "dbcdofs: " << cnode->DbcDofs()[0] << cnode->DbcDofs()[1] << cnode->DbcDofs()[2] << std::endl;
          dserror("Inconsistency in structure Dirichlet conditions and Mortar symmetry conditions");
        }
      }
    }
  }
  // create old style dirichtoggle vector (supposed to go away)
  pgsdirichtoggle_ = LINALG::CreateVector(*gsdofrowmap_,true);
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp,*pgsdirichtoggle_);

  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreDM(const std::string& state)
{
  //store Dold and Mold matrix in D and M
  if (state == "current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  }

  // store D and M matrix in Dold and Mold
  else if (state == "old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
    if (friction_ && Dualquadslave3d())
      doldmod_ = dmatrixmod_;
  }

  // unknown conversion
  else
  {
    dserror("ERROR: StoreDM: Unknown conversion requested!");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Store nodal quant. to old ones (last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::StoreToOld(
    MORTAR::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      FriNode* cnode = dynamic_cast<FriNode*>(node);

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
        dserror("ERROR: StoreDMToNodes: Unknown state std::string variable!");
        break;
      } // switch
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step             popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Update(Teuchos::RCP<Epetra_Vector> dis)
{
  // store Lagrange multipliers, D and M
  // (we need this for interpolation of the next generalized mid-point)
  // in the case of self contact, the size of z may have changed
  if (IsSelfContact())
    zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  zold_->Update(1.0, *z_, 0.0);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);
  StoreDM("old");

  // store contact state to contact nodes (active or inactive)
  StoreNodalQuantities(MORTAR::StrategyBase::activeold);

  // old displacements in nodes
  // (this is NOT only needed for friction but also for calculating
  // the auxiliary positions in binarytree contact search)
  SetState("olddisplacement", dis);

  // double-check if active set is really converged
  // (necessary e.g. for monolithic FSI with Lagrange multiplier contact,
  // because usually active set convergence check has been integrated into
  // structure Newton scheme, but now the monolithic FSI Newton scheme decides)
  if (!ActiveSetConverged() || !ActiveSetSemiSmoothConverged())
    dserror("ERROR: Active set not fully converged!");

  // reset active set status for next time step
  ResetActiveSet();

  // update flag for global contact status of last time step
  if (gactivenodes_->NumGlobalElements())
  {
    wasincontact_ = true;
    wasincontactlts_ = true;
  }
  else
  {
    wasincontact_ = false;
    wasincontactlts_ = false;
  }

  //----------------------------------------friction: store history values
  // in the case of frictional contact we have to store several
  // information and quantities at the end of a time step (converged
  // state) which is needed in the next time step as history
  // information / quantities.
  if (friction_)
  {
    // store contact state to friction nodes (slip or stick)
    StoreNodalQuantities(MORTAR::StrategyBase::slipold);

    // store nodal entries of D and M to old ones
    StoreToOld(MORTAR::StrategyBase::dm);

    // store nodal entries form penalty contact tractions to old ones
    StoreToOld(MORTAR::StrategyBase::pentrac);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoWriteRestart(
    std::map<std::string,Teuchos::RCP<Epetra_Vector> >& restart_vectors,
    bool forcedrestart)
{
  // initalize
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
  Teuchos::RCP<Epetra_Vector> sliptoggle = Teuchos::null;

  // write toggle
  restart_vectors["activetoggle"]=activetoggle;
  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    restart_vectors["sliptoggle"]=sliptoggle;
  }

  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      if (forcedrestart)
      {
        // set value active / inactive in toggle vector
        if (cnode->CoData().ActiveOld())
          (*activetoggle)[dof] = 1;
      }
      else
      {
        // set value active / inactive in toggle vector
        if (cnode->Active())
          (*activetoggle)[dof] = 1;
      }

      // set value slip / stick in the toggle vector
      if (friction_)
      {
        CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (forcedrestart)
        {
          if (frinode->FriData().SlipOld())
            (*sliptoggle)[dof] = 1;
        }
        else
        {
          if (frinode->FriData().Slip())
            (*sliptoggle)[dof] = 1;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::DoReadRestart(
    IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = DRT::INPUT::IntegralValue<int>(Params(),
      "RESTART_WITH_CONTACT");

  // set restart displacement state
  SetState("displacement", dis);
  SetState("olddisplacement", dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  InitMortar();
  InitEvalInterface();
  AssembleMortar();

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,
        false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp;
  }

  // read restart information on actice set and slip set (leave sets empty
  // if this is a restart with contact of a non-contact simulation run)
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(
      new Epetra_Vector(*gsnoderowmap_));
  if (!restartwithcontact)
    reader.ReadVector(activetoggle, "activetoggle");

  // friction
  Teuchos::RCP<Epetra_Vector> sliptoggle;
  Teuchos::RCP<Epetra_Vector> weightedwear;

  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    if (!restartwithcontact)
      reader.ReadVector(sliptoggle, "sliptoggle");
  }

  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < (interface_[i]->SlaveRowNodes())->NumMyElements(); ++j)
    {
      int gid = (interface_[i]->SlaveRowNodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof] == 1)
      {
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node)
          dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // set value active / inactive in cnode
        cnode->Active() = true;

        if (friction_)
        {
          // set value stick / slip in cnode
          if ((*sliptoggle)[dof] == 1)
            dynamic_cast<CONTACT::FriNode*>(cnode)->FriData().Slip() = true;
        }
      }
    }
  }

  // read restart information on Lagrange multipliers
  z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  if (!restartwithcontact)
    reader.ReadVector(LagrMult(), "lagrmultold");
  if (!restartwithcontact)
    reader.ReadVector(LagrMultOld(), "lagrmultold");

  // Lagrange multiplier increment is always zero (no restart value to be read)
  zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // store restart information on Lagrange multipliers into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // only for Uzawa augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  if (stype_ == INPAR::CONTACT::solution_uzawa)
  {
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (!restartwithcontact)
      reader.ReadVector(LagrMultUzawa(), "lagrmultold");
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
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,
        interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(),
        false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(),
        false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(),
        false);
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(),
          false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(),
          false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }
  }

  // update flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
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
  // check chosen output option
  INPAR::CONTACT::EmOutputType emtype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::EmOutputType>(Params(), "EMOUTPUT");

  // get out of here if no output wanted
  if (emtype == INPAR::CONTACT::output_none)
    return;

  // compute discrete slave and master interface forces
  Teuchos::RCP<Epetra_Vector> fcslavetemp = Teuchos::rcp(
      new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertemp = Teuchos::rcp(
      new Epetra_Vector(mmatrix_->DomainMap()));

  // for self contact, slave and master sets may have changed,
  // thus we have to export z to new D and M dimensions
  if (IsSelfContact())
  {
    Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(
        new Epetra_Vector(dmatrix_->RowMap()));
    if (dmatrix_->RowMap().NumGlobalElements())
      LINALG::Export(*z_, *zexp);
    dmatrix_->Multiply(true, *zexp, *fcslavetemp);
    mmatrix_->Multiply(true, *zexp, *fcmastertemp);
  }
  // if there is no self contact everything is ok
  else
  {
    dmatrix_->Multiply(true, *z_, *fcslavetemp);
    mmatrix_->Multiply(true, *z_, *fcmastertemp);
  }

  // export the interface forces to full dof layout
  Teuchos::RCP<Epetra_Vector> fcslave = Teuchos::rcp(
      new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> fcmaster = Teuchos::rcp(
      new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fcslavetemp, *fcslave);
  LINALG::Export(*fcmastertemp, *fcmaster);

  // contact forces and moments
  std::vector<double> gfcs(3);
  std::vector<double> ggfcs(3);
  std::vector<double> gfcm(3);
  std::vector<double> ggfcm(3);
  std::vector<double> gmcs(3);
  std::vector<double> ggmcs(3);
  std::vector<double> gmcm(3);
  std::vector<double> ggmcm(3);

  std::vector<double> gmcsnew(3);
  std::vector<double> ggmcsnew(3);
  std::vector<double> gmcmnew(3);
  std::vector<double> ggmcmnew(3);

  // weighted gap vector
  Teuchos::RCP<Epetra_Vector> gapslave = Teuchos::rcp(
      new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> gapmaster = Teuchos::rcp(
      new Epetra_Vector(mmatrix_->DomainMap()));

  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid < 0)
          dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d)
        gmcs[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      std::vector<int> lm(Dim());
      std::vector<int> lmowner(Dim());
      for (int d = 0; d < Dim(); ++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapslave, posnode, lm, lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j = 0; j < interface_[i]->MasterRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->MasterRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcmastertemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid < 0)
          dserror("ERROR: ContactForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d)
        gmcm[d] += nodemoment[d];

      // weighted gap
      Epetra_SerialDenseVector posnode(Dim());
      std::vector<int> lm(Dim());
      std::vector<int> lmowner(Dim());
      for (int d = 0; d < Dim(); ++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->Dofs()[d];
        lmowner[d] = cnode->Owner();
      }
      LINALG::Assemble(*gapmaster, posnode, lm, lmowner);
    }
  }

  // weighted gap
  Teuchos::RCP<Epetra_Vector> gapslavefinal = Teuchos::rcp(
      new Epetra_Vector(dmatrix_->RowMap()));
  Teuchos::RCP<Epetra_Vector> gapmasterfinal = Teuchos::rcp(
      new Epetra_Vector(mmatrix_->RowMap()));
  dmatrix_->Multiply(false, *gapslave, *gapslavefinal);
  mmatrix_->Multiply(false, *gapmaster, *gapmasterfinal);
  Teuchos::RCP<Epetra_Vector> gapfinal = Teuchos::rcp(
      new Epetra_Vector(dmatrix_->RowMap()));
  gapfinal->Update(1.0, *gapslavefinal, 0.0);
  gapfinal->Update(-1.0, *gapmasterfinal, 1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      std::vector<double> lm(3);
      std::vector<double> nodegaps(3);
      std::vector<double> nodegapm(3);

      // LMs and gaps
      for (int d = 0; d < Dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->Dofs()[d]);
        if (dofid < 0)
          dserror("ERROR: ContactForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = cnode->MoData().lm()[d];
      }

      // moments
      std::vector<double> nodemoments(3);
      std::vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1] * lm[2] - nodegaps[2] * lm[1];
      nodemoments[1] = nodegaps[2] * lm[0] - nodegaps[0] * lm[2];
      nodemoments[2] = nodegaps[0] * lm[1] - nodegaps[1] * lm[0];
      nodemomentm[0] = nodegapm[1] * lm[2] - nodegapm[2] * lm[1];
      nodemomentm[1] = nodegapm[2] * lm[0] - nodegapm[0] * lm[2];
      nodemomentm[2] = nodegapm[0] * lm[1] - nodegapm[1] * lm[0];
      for (int d = 0; d < 3; ++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }
    }
  }

  // summing up over all processors
  for (int i = 0; i < 3; ++i)
  {
    Comm().SumAll(&gfcs[i], &ggfcs[i], 1);
    Comm().SumAll(&gfcm[i], &ggfcm[i], 1);
    Comm().SumAll(&gmcs[i], &ggmcs[i], 1);
    Comm().SumAll(&gmcm[i], &ggmcm[i], 1);
    Comm().SumAll(&gmcsnew[i], &ggmcsnew[i], 1);
    Comm().SumAll(&gmcmnew[i], &ggmcmnew[i], 1);
  }

  // print interface results to file
  if (emtype == INPAR::CONTACT::output_file
      || emtype == INPAR::CONTACT::output_both)
  {
    // do this at end of time step only (output==true)!
    // processor 0 does all the work
    if (output && Comm().MyPID() == 0)
    {
      FILE* MyFile = NULL;
      std::ostringstream filename;
      filename << "interface.txt";
      MyFile = fopen(filename.str().c_str(), "at+");

      if (MyFile)
      {
        for (int i = 0; i < 3; i++)
          fprintf(MyFile, "%g\t", ggfcs[i]);
        for (int i = 0; i < 3; i++)
          fprintf(MyFile, "%g\t", ggfcm[i]);
        for (int i = 0; i < 3; i++)
          fprintf(MyFile, "%g\t", ggmcs[i]);
        for (int i = 0; i < 3; i++)
          fprintf(MyFile, "%g\t", ggmcm[i]);
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
  if (emtype == INPAR::CONTACT::output_screen
      || emtype == INPAR::CONTACT::output_both)
  {
    // do this during Newton steps only (output==false)!
    // processor 0 does all the work
    if (!output && Comm().MyPID() == 0)
    {
      double snorm = sqrt(
          ggfcs[0] * ggfcs[0] + ggfcs[1] * ggfcs[1] + ggfcs[2] * ggfcs[2]);
      double mnorm = sqrt(
          ggfcm[0] * ggfcm[0] + ggfcm[1] * ggfcm[1] + ggfcm[2] * ggfcm[2]);
      printf("Slave Contact Force:   % e  % e  % e \tNorm: % e\n", ggfcs[0],
          ggfcs[1], ggfcs[2], snorm);
      printf("Master Contact Force:  % e  % e  % e \tNorm: % e\n", ggfcm[0],
          ggfcm[1], ggfcm[2], mnorm);
      printf("Slave Contact Moment:  % e  % e  % e\n", ggmcs[0], ggmcs[1],
          ggmcs[2]);
      //printf("Slave Contact Moment:  % e  % e  % e\n",ggmcsnew[0],ggmcsnew[1],ggmcsnew[2]);
      printf("Master Contact Moment: % e  % e  % e\n", ggmcm[0], ggmcm[1],
          ggmcm[2]);
      //printf("Master Contact Moment: % e  % e  % e\n",ggmcmnew[0],ggmcmnew[1],ggmcmnew[2]);
      fflush(stdout);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    os << "--------------------------------- CONTACT::CoAbstractStrategy\n"
        << "Contact interfaces: " << (int) interface_.size() << std::endl
        << "-------------------------------------------------------------\n";
  }
  Comm().Barrier();
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    std::cout << *(interface_[i]);
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
  if (Comm().MyPID() == 0)
  {
    printf(
        "\nActive contact set--------------------------------------------------------------\n");
    fflush(stdout);
  }

  //**********************************************************************
  // detailed active set output
  //**********************************************************************

#ifdef CONTACTASOUTPUT
  // create storage for local and global data
  std::vector<int> lnid, gnid;
  std::vector<double> llmn, glmn;
  std::vector<double> lgap, ggap;

  std::vector<double> Xposl, Xposg;
  std::vector<double> Yposl, Yposg;
  std::vector<double> Zposl, Zposg;

  std::vector<double> xposl, xposg;
  std::vector<double> yposl, yposg;
  std::vector<double> zposl, zposg;

  // introduce integer variable status
  // (0=inactive, 1=active, 2=slip, 3=stick)
  // (this is necessary as all data will be written by proc 0, but
  // the knowledge of the above status ONLY exists on the owner
  // processor of the respective node. Thus this information must
  // also be communicated to proc 0 in addition to the actual data!)
  std::vector<int> lsta, gsta;

  // some more storage for local and global friction data
  std::vector<double> llmt, glmt;
  std::vector<double> ljtx, gjtx;
  std::vector<double> ljte, gjte;
  std::vector<double> lwear, gwear;

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
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // compute weighted gap
        double wgap = (*g_)[g_->Map().LID(gid)];

        double Xpos=cnode->X()[0];
        double Ypos=cnode->X()[1];
        double Zpos=cnode->X()[2];

        double xpos=cnode->xspatial()[0];
        double ypos=cnode->xspatial()[1];
        double zpos=cnode->xspatial()[2];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k=0;k<3;++k)
        nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);
        Xposl.push_back(Xpos);
        Yposl.push_back(Ypos);
        Zposl.push_back(Zpos);
        xposl.push_back(xpos);
        yposl.push_back(ypos);
        zposl.push_back(zpos);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active()) lsta.push_back(1);
        else lsta.push_back(0);
      }

      //--------------------------------------------------------------------
      // FRICTIONAL CASE
      //--------------------------------------------------------------------
      else
      {
        // cast to CoNode and FriNode
        CoNode* cnode = dynamic_cast<CoNode*>(node);
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

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
        if (weightedwear_) wear = frinode->WearData().WeightedWear();

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
          else lsta.push_back(3);
        }
        else
        {
          lsta.push_back(0);
        }
      }
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i=0; i<Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid,gnid,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(llmn,glmn,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(lgap,ggap,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<int>(lsta,gsta,(int)allproc.size(),&allproc[0],Comm());

  LINALG::Gather<double>(Xposl,Xposg,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(Yposl,Yposg,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(Zposl,Zposg,(int)allproc.size(),&allproc[0],Comm());

  LINALG::Gather<double>(xposl,xposg,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(yposl,yposg,(int)allproc.size(),&allproc[0],Comm());
  LINALG::Gather<double>(zposl,zposg,(int)allproc.size(),&allproc[0],Comm());

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
          printf("INACTIVE: %d \t wgap: % e \t lm: % e \t Xref: % e \t Yref: % e \t Zref: % e \n",gnid[k],ggap[k],glmn[k],Xposg[k],Yposg[k],Zposg[k]);
          fflush(stdout);
        }

        // print nodes of active set ***************************************
        else if (gsta[k]==1)
        {
          printf("ACTIVE:   %d \t wgap: % e \t lm: % e \t Xref: % e \t Yref: % e \t Zref: % e \n",gnid[k],ggap[k],glmn[k],Xposg[k],Yposg[k],Zposg[k]);
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
#ifdef CONTACTEXPORT
      // export variables
      double sum_jumpx=0.0;
      double sum_jumpe=0.0;
      double sum_jumpall=0.0;
      double iter_slip=0.0;
#endif

      // loop over all nodes
      for (int k=0;k<(int)gnid.size();++k)
      {
        // print nodes of slip set **************************************
        if (gsta[k]==2)
        {
          printf("SLIP:  %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",gnid[k],glmn[k],glmt[k],gjtx[k],gjte[k],gwear[k]);
          fflush(stdout);
#ifdef CONTACTEXPORT
          // preparation for output
          sum_jumpx+=gjtx[k];
          sum_jumpe+=gjte[k];
          sum_jumpall+=sqrt(gjtx[k]*gjtx[k] + gjte[k]*gjte[k]);
          iter_slip=iter_slip+1.0;
#endif
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

#ifdef CONTACTEXPORT
      // export averaged slip increments to xxx.jump
      double sum_jumpx_final=0.0;
      double sum_jumpe_final=0.0;
      double sum_jumpall_final=0.0;

      if (iter_slip>0.0)
      {
        sum_jumpx_final=sum_jumpx/iter_slip;
        sum_jumpe_final=sum_jumpe/iter_slip;
        sum_jumpall_final=sum_jumpall/iter_slip;
      }

      FILE* MyFile = NULL;
      std::ostringstream filename;
      const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
      filename << filebase << ".jump";
      MyFile = fopen(filename.str().c_str(), "at+");
      if (MyFile)
      {
        //fprintf(MyFile,valuename.c_str());
        fprintf(MyFile, "%g\t",sum_jumpx_final);
        fprintf(MyFile, "%g\t",sum_jumpe_final);
        fprintf(MyFile, "%g\n",sum_jumpall_final);
        fclose(MyFile);
      }
      else
      dserror("ERROR: File for Output could not be opened.");
#endif

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
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: PrintActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);

      // increase active counters
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      if (cnode->Active() or cnode->AugActive()) activenodes   += 1;
      else                                       inactivenodes += 1;

      // increase friction counters
      if (friction_)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);
        if (cnode->Active() && frinode->FriData().Slip())
          slipnodes += 1;
      }
    }
  }

  // sum among all processors
  Comm().SumAll(&activenodes, &gactivenodes, 1);
  Comm().SumAll(&inactivenodes, &ginactivenodes, 1);
  Comm().SumAll(&slipnodes, &gslipnodes, 1);

  // print active set information
  if (Comm().MyPID() == 0)
  {
    if (friction_)
    {
      std::cout << BLUE2_LIGHT << "Total     SLIP nodes:\t" << gslipnodes
          << END_COLOR << std::endl;
      std::cout << BLUE2_LIGHT << "Total    STICK nodes:\t"
          << gactivenodes - gslipnodes << END_COLOR << std::endl;
      std::cout << RED_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes
          << END_COLOR << std::endl;
    }
    else
    {
      std::cout << BLUE2_LIGHT << "Total   ACTIVE nodes:\t" << gactivenodes
          << END_COLOR << std::endl;
      std::cout << RED_LIGHT << "Total INACTIVE nodes:\t" << ginactivenodes
          << END_COLOR << std::endl;
    }
  }
#endif // #ifdef CONTACTASOUTPUT
  // output line
  Comm().Barrier();
  if (Comm().MyPID() == 0)
  {
    printf(
        "--------------------------------------------------------------------------------\n\n");
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
  for (int i = 0; i < (int) interface_.size(); ++i)
    interface_[i]->VisualizeGmsh(step, iter);
}

/*----------------------------------------------------------------------*
 | collect maps for preconditioner (AMG)                   wiesner 01/12|
 *----------------------------------------------------------------------*/
void CONTACT::CoAbstractStrategy::CollectMapsForPreconditioner(
    Teuchos::RCP<Epetra_Map>& MasterDofMap,
    Teuchos::RCP<Epetra_Map>& SlaveDofMap,
    Teuchos::RCP<Epetra_Map>& InnerDofMap,
    Teuchos::RCP<Epetra_Map>& ActiveDofMap)
{
  InnerDofMap = gndofrowmap_; // global internal dof row map
  ActiveDofMap = gactivedofs_; // global active slave dof row map

  // check if parallel redistribution is used
  // if parallel redistribution is activated, then use (original) maps before redistribution
  // otherwise we use just the standard master/slave maps
  if (pgsdofrowmap_ != Teuchos::null)
    SlaveDofMap = pgsdofrowmap_;
  else
    SlaveDofMap = gsdofrowmap_;
  if (pgmdofrowmap_ != Teuchos::null)
    MasterDofMap = pgmdofrowmap_;
  else
    MasterDofMap = gmdofrowmap_;
}
