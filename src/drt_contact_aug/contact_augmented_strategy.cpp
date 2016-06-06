/*----------------------------------------------------------------------*/
/*!
\file contact_augmented_strategy.cpp

\brief Augmented contact solving strategy with standard Lagrangian
       multipliers.

\level 3

\maintainer Michael Hiermeier

\date Apr 7, 2014

*/
/*----------------------------------------------------------------------*/
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_lagrange_strategy.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AugStratDataContainer::AugStratDataContainer()
    : gFdCheck_(false),
      wasincontactlastiter_(false),
      isactivesetconverged_(false),
      printlinearconservation_(false),
      printangularconservation_(false),
      dGLmSlLinMatrixPtr_(Teuchos::null),
      dGLmMaLinMatrixPtr_(Teuchos::null),
      dGGSlLinMatrixPtr_(Teuchos::null),
      dGGMaLinMatrixPtr_(Teuchos::null),
      dLmNWGapLinMatrixPtr_(Teuchos::null),
      dLmTLmTMatrixPtr_(Teuchos::null),
      dLmTLmTLinMatrixPtr_(Teuchos::null),
      inactiveMatrixPtr_(Teuchos::null),
      inactiveLinMatrixPtr_(Teuchos::null),
      aPtr_(Teuchos::null),
      kappaPtr_(Teuchos::null),
      lmNPtr_(Teuchos::null),
      dLmTLmTRhsPtr_(Teuchos::null),
      strContactRhsPtr_(Teuchos::null),
      cnPtr_(Teuchos::null),
      uCnPtr_(Teuchos::null),
      gsndofrowmapPtr_(Teuchos::null),
      gstdofrowmapPtr_(Teuchos::null),
      gOldActiveSlaveNodesPtr_(Teuchos::null)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AugmentedLagrangeStrategy::AugmentedLagrangeStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
    const int& dim,
    const Teuchos::RCP<const Epetra_Comm>& comm,
    const int& maxdof)
    : CoAbstractStrategy(data_ptr,DofRowMap,NodeRowMap,params,interfaces,dim,comm,0.0,maxdof)
{
  augDataPtr_ = Teuchos::rcp_dynamic_cast<CONTACT::AugStratDataContainer>(data_ptr);
  // store values of the parameter list
  Data().PrintLinearMomConservation() =
      DRT::INPUT::IntegralValue<bool>(params.sublist("AUGMENTED"),
          "PRINT_LINEAR_CONSERVATION");
  Data().PrintAngularMomConservation() =
        DRT::INPUT::IntegralValue<bool>(params.sublist("AUGMENTED"),
            "PRINT_ANGULAR_CONSERVATION");

  // cast to augmented interfaces
  for (int i=0; i<(int) interfaces.size(); ++i)
  {
    interface_.push_back(Teuchos::rcp_dynamic_cast<CONTACT::AugmentedInterface>(interfaces[i]));
    if (interface_[i]==Teuchos::null)
      dserror("AugmentedLagrangeStartegy: Interface-cast failed!");
  }

  // re-setup the global sndofrowmap_ and stdofrowmap_ for the augmented Lagrangian case
  AssembleGlobalSlNTDofRowMaps();

  // initialize cn
  InitializeCn(false,true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::PostSetup(
    const bool& redistributed,
    const bool& init)
{
  // just used for the redistributed case
  if (not redistributed) return;

  // redistribute the slave/master dirichlet boundary condition row map
  AssembleSlMaNoDbcDofRowMap(Teuchos::null);
  // reassemble the global slave normal/tangential dof row maps
  AssembleGlobalSlNTDofRowMaps();
  // redistribute the cn-vector
  InitializeCn(redistributed,false);
  // redistribute the global augmented old active slave nodes map
  if ((not Data().GOldActiveSlaveNodesPtr().is_null())and
      (Data().GOldActiveSlaveNodes().NumGlobalElements()>0))
    RedistributeRowMap(SlRowNodes(),Data().GOldActiveSlaveNodes());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AssembleGlobalSlNTDofRowMaps()
{
  Data().GSlNormalDofRowMapPtr()     = Teuchos::null;
  Data().GSlTangentialDofRowMapPtr() = Teuchos::null;

  for (int i=0; i<(int) interface_.size(); ++i)
  {
    Data().GSlNormalDofRowMapPtr()     =
        LINALG::MergeMap(Data().GSlNormalDofRowMapPtr(),interface_[i]->SlaveRowNDofs());
    Data().GSlTangentialDofRowMapPtr() =
        LINALG::MergeMap(Data().GSlTangentialDofRowMapPtr(),interface_[i]->SlaveRowTDofs());
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::InitializeCn(
    bool redistributed,
    bool init)
{
  if (init)
  {
    // get cn from the input file
    double cn = Params().get<double>("SEMI_SMOOTH_CN");

    if (Data().CnPtr().is_null() or Data().Cn().GlobalLength()==0)
      Data().CnPtr() = LINALG::CreateVector(SlRowNodes(),true);
    // set all nodal cn-values to the input value
    Data().Cn().PutScalar(cn);

    if (Data().UCnPtr().is_null() or Data().UCn().GlobalLength()==0)
      Data().UCnPtr() = LINALG::CreateVector(SlRowNodes(),true);
    // set all nodal entries of the uCn-vector
    Data().UCn().PutScalar(1.0e12);
  }

  if (redistributed)
  {
    // redistribute the cn-vector
    Teuchos::RCP<Epetra_Vector> newcn = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
    LINALG::Export(Data().Cn(),*newcn);
    Data().CnPtr() = newcn;
    // redistribute the cn-update-vector
    Teuchos::RCP<Epetra_Vector> newuCn = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
    LINALG::Export(Data().UCn(),*newuCn);
    Data().UCnPtr() = newuCn;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::InitMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing Dn, Mn etc.
  if (IsSelfContact()) UpdateMasterSlaveSetsGlobal();

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  Data().DMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(SlRowNodes(),100));
  Data().MMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(SlRowNodes(),100));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AssembleMortar()
{
  // for all interfaces
  for (std::size_t i = 0; i < interface_.size(); ++i)
    interface_[i]->AssembleAugDnMnMatrix(Data().DMatrix(),Data().MMatrix());

  Data().DMatrix().Complete(SlDoFRowMap(true),SlRowNodes());
  Data().MMatrix().Complete(MaDoFRowMap(true),SlRowNodes());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::SplitMortar()
{
  // temporal variables
  Teuchos::RCP<LINALG::SparseMatrix> DnActive, MnActive;

  if (Data().GActiveNodeRowMap().NumGlobalElements() == 0)
  {
    Data().DMatrixPtr() = Teuchos::rcp(
        new LINALG::SparseMatrix(Data().GActiveNDofRowMap(),100));
    Data().MMatrixPtr() = Teuchos::rcp(
        new LINALG::SparseMatrix(Data().GActiveNDofRowMap(),100));
    Data().DMatrix().Complete(SlDoFRowMap(true),Data().GActiveNDofRowMap());
    Data().MMatrix().Complete(MaDoFRowMap(true),Data().GActiveNDofRowMap());
  }
  else
  {
    // *** START - TIME MEASUREMENT ***
//    const double t_start = Teuchos::Time::wallTime();

    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx12, tempmtx21,tempmtx22;
    /*--------------------------------------------------------*
     |  Split Dn/Mn-matrices                                  |
     *--------------------------------------------------------*/
    // Split Dn-Matrix
    LINALG::SplitMatrix2x2(Data().DMatrixPtr(),Data().GActiveNodeRowMapPtr(),
        emptymap,Data().GSlDofRowMapPtr(),emptymap,DnActive,tempmtx12,tempmtx21,tempmtx22);
    // Split Mn-Matrix
    LINALG::SplitMatrix2x2(Data().MMatrixPtr(),Data().GActiveNodeRowMapPtr(),
        emptymap,Data().GMaDofRowMapPtr(),emptymap,MnActive,tempmtx12,tempmtx21,tempmtx22);

    // change row maps from nodal gids to ndof gids
    if (DnActive->RowMap().NumGlobalElements() == Data().GActiveNDofRowMap().NumGlobalElements())
      Data().DMatrixPtr() = MORTAR::MatrixRowTransformGIDs(DnActive,Data().GActiveNDofRowMapPtr());
    else
      dserror("The row-map can't be replaced! The number of entries does not fit.");

    if (MnActive->RowMap().NumGlobalElements() == Data().GActiveNDofRowMap().NumGlobalElements())
      Data().MMatrixPtr() = MORTAR::MatrixRowTransformGIDs(MnActive,Data().GActiveNDofRowMapPtr());
    else
      dserror("The row-map can't be replaced! The number of entries does not fit.");

    // *** END - TIME MEASUREMENT ***
//    const double t_end = Teuchos::Time::wallTime() - t_start;
//    std::cout << "=== Split Mortar Time: " << std::scientific << t_end << "[sec] === (PID: " << Comm().MyPID() << ")" << std::endl;
//    Comm().Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::Initialize()
{
  // *** (re)setup global matrices ***
  Data().DGLmSlLinMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      SlDoFRowMap(true),100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  Data().DGLmMaLinMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      MaDoFRowMap(true),100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  Data().DGGSlLinMatrixPtr()  = Teuchos::rcp(new LINALG::SparseMatrix(
      SlDoFRowMap(true),100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  Data().DGGMaLinMatrixPtr()  = Teuchos::rcp(new LINALG::SparseMatrix(
      MaDoFRowMap(true),100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  Data().DLmNWGapLinMatrixPtr() = Teuchos::rcp(
      new LINALG::SparseMatrix(Data().GActiveNDofRowMap(),100));
  Data().DLmTLmTMatrixPtr()     = Teuchos::rcp(
      new LINALG::SparseMatrix(Data().GActiveTDofRowMap(),100));
  Data().DLmTLmTLinMatrixPtr()  = Teuchos::rcp(
      new LINALG::SparseMatrix(Data().GActiveTDofRowMap(),100));

  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs =
      LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());

  Data().InactiveMatrixPtr()    =
      Teuchos::rcp(new LINALG::SparseMatrix(*gAugInactiveSlaveDofs,100));
  Data().InactiveLinMatrixPtr() =
      Teuchos::rcp(new LINALG::SparseMatrix(*gAugInactiveSlaveDofs,100));

  // *** (re)setup global augmented Epetra_Vectors ***
  Data().LmNPtr()         = LINALG::CreateVector(
      Data().GActiveNDofRowMap(),true);
  Data().AWGapPtr()       = LINALG::CreateVector(
      Data().GActiveNDofRowMap(),true);
  Data().DLmTLmTRhsPtr()  = LINALG::CreateVector(
      Data().GActiveTDofRowMap(),true);
  Data().InactiveRhsPtr() = LINALG::CreateVector(
      *gAugInactiveSlaveDofs,true);

  Data().AVecPtr()     = LINALG::CreateVector(SlRowNodes(),true);
  Data().KappaVecPtr() = LINALG::CreateVector(
      Data().GActiveNodeRowMap(),true);
  Data().WGapPtr()  = LINALG::CreateVector(
      Data().GActiveNDofRowMap(),true);

  if (Data().IsFriction())
    dserror("AugmentedLagrangeStrategy::Initialize: "
        "Frictional case is not yet considered!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalForceStiff(
    CONTACT::ParamsInterface& cparams)
{
  // call the evaluate force routine
  EvalForce(cparams);

  // --- Assemble stiffness matrix ---------------------------------------
  AssembleContactStiff();

  // --- DEBUGGING -------------------------------------------------------
  // Finite Difference check at Gauss-point level
  AugFDCheckGP(cparams);

  // Finite Difference check at global level
  AugFDCheckGlobal(cparams);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::InitEvalInterface(
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // get type of parallel strategy
  INPAR::MORTAR::ParallelStrategy strat = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::ParallelStrategy>(Params(), "PARALLEL_STRATEGY");

  // Evaluation for all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    // initialize / reset interfaces
    interface_[i]->Initialize();

    //store required integration time
    Data().IntTime() += interface_[i]->Inttime();
    switch (strat)
    {
      /*----------------------------------------------------------*
       |  Fully redundant ghosting of master side                 |
       *----------------------------------------------------------*/
      case INPAR::MORTAR::ghosting_redundant:
      {
        // evaluate averaged weighted gap
        interface_[i]->Evaluate(0,cparams_ptr);
        // Calculate weighted gap
        interface_[i]->WGap();
        // evaluate remaining entities and linearization
        interface_[i]->RedEvaluate(cparams_ptr);
        break;
      }
      default:
      {
        dserror("ERROR: Augmented Lagrange strategy supports only "
            "\"ghosting_redundant\" as \"PARALLEL_STRATEGY\".");
        break;
      }
    }
  } // end interface loop

  // check the parallel distribution
  CheckParallelDistribution(t_start);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::RecoverState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, the recovery is straight forward. */
  Teuchos::RCP<Epetra_Vector> znew_ptr =
      Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(true),true));
  LINALG::Export(xnew,*znew_ptr);
  Teuchos::RCP<Epetra_Vector> zdir_ptr =
      Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(true),true));
  LINALG::Export(dir,*zdir_ptr);
  // get the current step length
  const double stepLength = cparams.GetStepLength();
  // ---------------------------------------------------------------------
  // Update the current lagrange multiplier
  // ---------------------------------------------------------------------
  znew_ptr->ReplaceMap(Data().LmPtr()->Map());
  Data().LmPtr()->Scale(1.0,*znew_ptr);
  // ---------------------------------------------------------------------
  // store the SCALED Lagrange multiplier increment in the contact
  // strategy
  // ---------------------------------------------------------------------
  zdir_ptr->ReplaceMap(Data().LmIncrPtr()->Map());
  Data().LmIncrPtr()->Scale(stepLength,*zdir_ptr);
  // ---------------------------------------------------------------------
  // store the new Lagrange multiplier in the nodes
  // ---------------------------------------------------------------------
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalForce(
    CONTACT::ParamsInterface& cparams)
{
  /*--------------------------------------------------------------*
   | For self-contact the master/slave sets are updated within the|
   | contact search, see SelfBinaryTree.                          |
   | Therefore, we have to initialize the mortar matrices after   |
   | interface evaluations.                                       |
   *--------------------------------------------------------------*/
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
    Teuchos::rcp(&cparams,false);
  if (IsSelfContact())
  {
    // evaluate mortar terms (integrate...)
    InitEvalInterface(cparams_ptr);
    // initialize mortar matrices and vectors
    InitMortar();
    // assemble mortar terms into global matrices
    AssembleMortar();
  }
  else
  {
    // initialize mortar matrices and vectors
    InitMortar();
    // evaluate mortar terms (integrate...)
    InitEvalInterface(cparams_ptr);
    // assemble mortar terms into global matrices
    AssembleMortar();
  }

  if (cparams.IsPredictor())
  {
    // evaluate relative movement for friction
    EvaluateRelMovPredict();
  }
  else
    EvaluateRelMov();

  // update active set
  UpdateActiveSetSemiSmooth();

  /* Split the Dn/Mn matrices to get only the active rows
   * (only necessary for the augmented Lagrangian formulation) */
  SplitMortar();

  // initialize all rhs vectors and linearization matrices
  Initialize();

  // --- Assemble the ride hand side terms -------------------------------
  AssembleContactRHS();

  /* Evaluate structural and constraint rhs. This is also necessary, if the
   * rhs did not change during the predictor step, but a redistribution was
   * executed! */
  EvalStrContactRHS();   // update structural contact rhs
  EvalConstrRHS();       // update the constrRHS

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AssembleContactRHS()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;
  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE CONTACT RIGHT HAND SIDE                   |
   *--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*
   | calculate                                                          |
   *--------------------------------------------------------------------*
   | --> normal Lagrange multiplier                                     |
   | --> (averaged) weighted gap                                        |
   | --> tangential constraint right hand side for the frictionless case|
   | --> normal and tangential inactive rhs                             |
   *--------------------------------------------------------------------*/
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // --- augmented Lagrange formulation --------------------------------
    // --- Force Balance -------------------------------------------------
    interface_[i]->AssembleResidualVectors(Data().LmN(),
        Data().AWGap(), Data().WGap());
    // --- CONSTRAINTS ---------------------------------------------------
    // active - normal direction
    // --> wGapRhs_
    // active - tangential direction
    interface_[i]->AssembleDLmTLmTRhs(Data().DLmTLmTRhs());
    // inactive - all directions
    interface_[i]->AssembleAugInactiveRhs(Data().InactiveRhs(),Data().Cn());
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AssembleContactStiff()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;
  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE TANGENTIAL STIFFNESS MATRIX               |
   *--------------------------------------------------------------------*/
  // --- augmented Lagrange formulation ----------------------------------
  for (std::size_t i=0; i < interface_.size();++i)
  {
    // Calculate averaged weighted gap linearization for all active nodes
    interface_[i]->AWGapLin();
    // --- Force Balance ------------------------------------------------
    // linearization w.r.t. displ.
    interface_[i]->AssembleDGLmLinMatrix(Data().DGLmSlLinMatrix(),
                                         Data().DGLmMaLinMatrix());
    interface_[i]->AssembleDGGLinMatrix(Data().DGGSlLinMatrix(),
                                        Data().DGGMaLinMatrix(),
                                        Data().Cn());
    // --- Constraints --------------------------------------------------
    // linearization w.r.t. LM
    interface_[i]->AssembleDLmTLmTMatrix(Data().DLmTLmTMatrix());
    interface_[i]->AssembleAugInactiveMatrix(Data().InactiveMatrix(),Data().Cn());
    // linearization w.r.t. displ.
    // active - normal direction
    interface_[i]->AssembleDLmNWGapLinMatrix(Data().DLmNWGapLinMatrix());
    // active - tangential direction
    interface_[i]->AssembleDLmTLmTLinMatrix(Data().DLmTLmTLinMatrix());
    /*--------------------------------------------------------------------*
     | The linearization of the nodal area w.r.t. the displ. for inactive |
     | nodes can help to prevent oscillations of the active set, because  |
     | the reduction of the inactive lm-value is decelerated.             |
     |                                                                    |
     | In general, the following relation does NOT hold anymore:          |
     |                    Delta(z_{n,i}^{k+1}) = - z_{n,i}^{k}            |
     |                          z_{n,i}^{k+1}  =   0                      |
     *--------------------------------------------------------------------*/
    interface_[i]->AssembleAugInactiveLinMatrix(Data().InactiveLinMatrix(),Data().Cn());
  }

  // --- START - FillComplete matrices ----------------------------------
  // domainmap: columnmap | rangemap: rowmap
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = LINALG::SplitMap(SlDoFRowMap(true),
      Data().GActiveDofRowMap());

  // --- Force Balance --------------------------------------------------
  // linearization w.r.t. displ.
  Data().DGLmSlLinMatrix().Complete(SlMaDoFRowMap(true),SlDoFRowMap(true));
  Data().DGLmMaLinMatrix().Complete(SlMaDoFRowMap(true),MaDoFRowMap(true));
  Data().DGGSlLinMatrix().Complete(SlMaDoFRowMap(true),SlDoFRowMap(true));
  Data().DGGMaLinMatrix().Complete(SlMaDoFRowMap(true),MaDoFRowMap(true));
  // --- Constraints ----------------------------------------------------
  // linearization w.r.t. LM
  Data().DLmTLmTMatrix().Complete(Data().GActiveTDofRowMap(),Data().GActiveTDofRowMap());
  Data().InactiveMatrix().Complete(*gAugInactiveSlaveDofs,*gAugInactiveSlaveDofs);
  // linearization w.r.t. displ.
  Data().DLmNWGapLinMatrix().Complete(SlMaDoFRowMap(true),Data().GActiveNDofRowMap());
  Data().DLmTLmTLinMatrix().Complete(SlDoFRowMap(true),Data().GActiveTDofRowMap());
  Data().InactiveLinMatrix().Complete(SlDoFRowMap(true),*gAugInactiveSlaveDofs);
  // --- END - FillComplete matrices ------------------------------------

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::UpdateActiveSetSemiSmooth(
    const bool& correction)
{
  // get out of here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  bool semismooth = DRT::INPUT::IntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
  if (!semismooth)
  {
    // loop over all interfaces
    for (std::size_t i=0; i<interface_.size();++i)
    {
      // loop over all slave nodes on the current interface
      for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cnode = static_cast<CoNode*>(node);

        /* The nested active set strategy cannot deal with the case of
         * active nodes that have no integration segments/cells attached,
         * as this leads to zero rows in D and M and thus to singular systems.
         * However, this case might possibly happen when slave nodes slide
         * over the edge of a master body within one fixed active set step.
         * (Remark: Semi-smooth Newton has no problems in this case, as it
         * updates the active set after EACH Newton step, see below, and thus
         * would always set the corresponding nodes to INACTIVE.) */
        if (cnode->Active() && !cnode->HasSegment())
          dserror("ERROR: Active node %i without any segment/cell attached",cnode->Id());
      }
    }
    return;
  }

  // assume that active set has converged and check for opposite
  Data().IsActiveSetConverged() = true;

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave nodes of the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      CoNode* cnode = dynamic_cast<CoNode*>(interface_[i]->Discret().gNode(gid));
      if (!cnode) dserror("ERROR: AugmentedInterface::UpdateAugActiveSetSemiSmooth: "
          "Cannot find node with gid %",gid);

      /* read weighting factor cn
       * (this is necessary in semi-smooth Newton case, as the search for the
       * active set is now part of the Newton iteration. Thus, we do not know
       * the active / inactive status in advance and we can have a state in
       * which both the condition znormal = 0 and wgap = 0 are violated. Here
       * we have to weight the two violations via cn! */
      const double& cn = Data().Cn()[Data().Cn().Map().LID(gid)];

      // compute averaged weighted gap
      double kappa = cnode->CoData().GetKappa();
      double awgap = cnode->CoData().GetWGap();
      if (kappa != 1.0e12)
        awgap/= kappa;

      // get normal part of the Lagrange multiplier
      double nz = cnode->MoData().lm()[0];

      // check nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        // check for fulfillment of contact condition
        if (nz - cn*awgap > 0.0)
        {
          cnode->Active() = true;
        }
      }
      // check nodes of active set ***************************************
      else
      {
        if (nz-cn*awgap<=0.0)
        {
          cnode->Active() = false;
        }
      }
    }
  } // end loop over all interfaces

  // only if it's a full Newton step...
  if (Data().IterLS() <= 0)
  {
    // store the previous augmented active set
    if (Data().GActiveNodeRowMapPtr() != Teuchos::null)
      Data().GOldActiveSlaveNodesPtr() = Teuchos::rcp(new Epetra_Map(Data().GActiveNodeRowMap()));
    else
      Data().GOldActiveSlaveNodesPtr() = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
  }

  // (re)setup of the global Epetra_maps
  Data().GActiveNodeRowMapPtr() = Teuchos::null;
  Data().GActiveDofRowMapPtr()  = Teuchos::null;
  Data().GActiveNDofRowMapPtr() = Teuchos::null;
  Data().GActiveTDofRowMapPtr() = Teuchos::null;

  // loop over all interfaces
  for (int it=0;it<(int) interface_.size();++it)
  {
    // update active set Epetra_Maps
    interface_[it]->BuildActiveSet();

    // Update Active set
    Data().GActiveNodeRowMapPtr() = LINALG::MergeMap(
        Data().GActiveNodeRowMapPtr(),interface_[it]->ActiveNodes(),false);
    Data().GActiveDofRowMapPtr()  = LINALG::MergeMap(
        Data().GActiveDofRowMapPtr(),interface_[it]->ActiveDofs(),false);
    Data().GActiveNDofRowMapPtr() = LINALG::MergeMap(
        Data().GActiveNDofRowMapPtr(),interface_[it]->ActiveNDofs(),false);
    Data().GActiveTDofRowMapPtr() = LINALG::MergeMap(
        Data().GActiveTDofRowMapPtr(),interface_[it]->ActiveTDofs(),false);
  }

  // check the convergence of the active set
  Data().IsActiveSetConverged() =
      Data().GActiveNodeRowMap().SameAs(Data().GOldActiveSlaveNodes());

  // update the history information only if it's no correction step of the active set
  if (!correction)
  {
    // update flag for the contact status of the last iterate (history information)
    if (IsInContact())
      Data().WasInContactLastIter() = true;
    else
      Data().WasInContactLastIter() = false;
  }
  // update flag for global contact status
  if (Data().GActiveNodeRowMap().NumGlobalElements())
  {
    Data().IsInContact()  = true;
    Data().WasInContact() = true;
  }
  else
    Data().IsInContact() = false;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    Data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  Data().StrContactRhsPtr() =
      Teuchos::rcp(new Epetra_Vector(*ProblemDofs(),true));


  // For self contact, slave and master sets may have changed,
  // thus we have to export the products Dold^T * zold / D^T * z to fit
  // thus we have to export the products Mold^T * zold / M^T * z to fit
  if (IsSelfContact())
    dserror("ERROR: Augmented Lagrange Formulation: Self contact is not yet considered!");
  // if there is no self contact everything is ok
  else
  {
    Teuchos::RCP<Epetra_Vector> augLmN = LINALG::CreateVector(Data().GActiveNDofRowMap(),true);
    // build the augmented Lagrangian multiplier vector
    augLmN->Update(1.0,Data().AWGap(),0.0);
    // scale by cn
    MultiplyElementwise(Data().Cn(),Data().GActiveNodeRowMap(),*augLmN,false);
    augLmN->Update(1.0,Data().LmN(),-1.0);
    // --- add contact force terms
    // *** Slave side ***
    Teuchos::RCP<Epetra_Vector> augfs = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
    Data().DMatrix().Multiply(true,*augLmN,*augfs);
    Teuchos::RCP<Epetra_Vector> augfs_exp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*augfs,*augfs_exp);
    Data().StrContactRhs().Update(-1.0,*augfs_exp,0.0);

    // Master side
    Teuchos::RCP<Epetra_Vector> augfm = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)));
    Data().MMatrix().Multiply(true,*augLmN,*augfm);
    Teuchos::RCP<Epetra_Vector> augfm_exp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*augfm,*augfm_exp);
    Data().StrContactRhs().Update(-1.0,*augfm_exp,1.0);

    // Check linear and angular momentum conservation
    CheckConservationLaws(*augfs,*augfm);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalConstrRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    // (re)setup the vector
    Data().ConstrRhsPtr() = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> augConstrRhs = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true),true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!

  // ToDo Three export vectors seem unnecessary and can be replaced by scaling a general export vector
  // with zero, between the different steps! Check for possible performance gain!
  // Add active constraints in normal direction:
  Teuchos::RCP<Epetra_Vector> wGapRhs_exp = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  LINALG::Export(Data().WGap(),*wGapRhs_exp);
  augConstrRhs->Update(1.0,*wGapRhs_exp,0.0);

  // Add inactive constraints
  Teuchos::RCP<Epetra_Vector> augInactiveRhs_exp = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  LINALG::Export(Data().InactiveRhs(),*augInactiveRhs_exp);
  augConstrRhs->Update(1.0,*augInactiveRhs_exp,1.0);

  // Add tangential frictionless constraints
  Teuchos::RCP<Epetra_Vector> dLmTLmTRhs_exp = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  LINALG::Export(Data().DLmTLmTRhs(),*dLmTLmTRhs_exp);
  augConstrRhs->Update(1.0,*dLmTLmTRhs_exp,1.0);

  // replace row map
  augConstrRhs->ReplaceMap(LMDoFRowMap(true));

  // export and set constraint rhs vector
  if (ParRedist())
  {
    Data().ConstrRhsPtr() = Teuchos::rcp(new Epetra_Vector(LMDoFRowMap(false)));
    LINALG::Export(*augConstrRhs,*Data().ConstrRhsPtr());
  }
  else
    Data().ConstrRhsPtr() = augConstrRhs;


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::CheckConservationLaws(
    const Epetra_Vector& augfs,
    const Epetra_Vector& augfm)
{
  if (not Data().PrintLinearMomConservation() and
      not Data().PrintAngularMomConservation())
    return;
  /****************** SLAVE SIDE ********************************************************/
  // standard Lagrange multiplier fraction
  Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Teuchos::RCP<Epetra_Vector> z_exp = Teuchos::rcp(new Epetra_Vector(Data().GActiveNDofRowMap()));
  LINALG::Export(*z_,*z_exp);
  Data().DMatrix().Multiply(true,*z_exp,*augfs_lm);
  // regularization fraction
  Teuchos::RCP<Epetra_Vector> augfs_g = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Data().DMatrix().Multiply(true,Data().AWGap(),*augfs_g);
  // *** START: DEBUGGING OUTPUT  ***
//  // COMPARE two different assembly strategies
//  double cn= Params().get<double>("SEMI_SMOOTH_CN");
//  Teuchos::RCP<Epetra_Vector> augfs_check = Teuchos::rcp(new Epetra_Vector(*augfs_lm));
//  augfs_check->Update(-cn,*augfs_g,1.0);
//
//  std::cout << ">>> augfs <<<" << std::endl;
//  std::cout << augfs << std::endl;
//  std::cout << "======================================" << std::endl;
//  std::cout << ">>> augfs_lm <<<" << std::endl;
//  std::cout << *augfs_lm << std::endl;
//  std::cout << ">>> augfs_g <<<" << std::endl;
//  std::cout << *augfs_g << std::endl;
//  std::cout << "======================================" << std::endl;
//  std::cout << ">>> augfs_check<<<" << std::endl;
//  std::cout << *augfs_check << std::endl;
  // *** END: DEBUGGING OUTPUT  ***
  /****************** MASTER SIDE *******************************************************/
  // standard lagrange multiplier fraction
  Teuchos::RCP<Epetra_Vector> augfm_lm = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)));
  Data().MMatrix().Multiply(true,*z_exp,*augfm_lm);
  // regularization fraction
  Teuchos::RCP<Epetra_Vector> augfm_g = Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)));
  Data().MMatrix().Multiply(true,Data().AWGap(),*augfm_g);
  // *** START: DEBUGGING OUTPUT  ***
//  // COMPARE two different assembly strategies
//  Teuchos::RCP<Epetra_Vector> augfm_check = Teuchos::rcp(new Epetra_Vector(*augfm_lm));
//  augfm_check->Update(-cn,*augfm_g,1.0);
//
//    std::cout << ">>> augfm <<<" << std::endl;
//    std::cout << augfm << std::endl;
//    std::cout << "======================================" << std::endl;
//    std::cout << ">>> augfm_lm <<<" << std::endl;
//    std::cout << *augfm_lm << std::endl;
//    std::cout << ">>> augfm_g <<<" << std::endl;
//    std::cout << *augfm_g << std::endl;
//    std::cout << "======================================" << std::endl;
//    std::cout << ">>> augfm_check<<<" << std::endl;
//    std::cout << *augfm_check << std::endl;
  // *** END: DEBUGGING OUTPUT  ***
  /*-------------------------------*
   | LINEAR MOMENTUM CONSERVATION  |
   *-------------------------------*/
  if (Data().PrintLinearMomConservation())
  {
    double lssum = 0.0;   // local slave sum
    double gssum = 0.0;   // global slave sum
    double lmsum = 0.0;   // local master sum
    double gmsum = 0.0;   // global master sum
    double gcsum = 0.0;   // global complete sum
    // slave
    for (int i=0;i<augfs_lm->MyLength();++i) lssum+=(*augfs_lm)[i];
    Comm().SumAll(&lssum,&gssum,1);
    // master
    for (int i=0;i<augfm_lm->MyLength();++i) lmsum+=(*augfm_lm)[i];
    Comm().SumAll(&lmsum,&gmsum,1);
    // complete balance check
    gcsum = gssum+gmsum;
    if (abs(gcsum)>1.0e-11) dserror("Conservation of linear momentum is not fulfilled!");
    if (Comm().MyPID()==0)
    {
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << ">>      Linear Momentum Conservation      <<" << std::endl;
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      std::cout << "SLAVE:   " << std::setw(14) << gssum<< std::endl;
      std::cout << "MASTER:  " << std::setw(14) << gmsum << std::endl;
      std::cout << "Balance: " << std::setw(14) << gcsum << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
    }
    // slave
    lssum=0.0;
    for (int i=0;i<augfs_g->MyLength();++i) lssum+=(*augfs_g)[i];
    Comm().SumAll(&lssum,&gssum,1);
    // master
    lmsum=0.0;
    for (int i=0;i<augfm_g->MyLength();++i) lmsum+=(*augfm_g)[i];
    Comm().SumAll(&lmsum,&gmsum,1);
    // complete balance check
    gcsum = gssum+gmsum;
    if (abs(gcsum)>1.0e-11) dserror("Conservation of linear momentum is not fulfilled!");
    if (Comm().MyPID()==0)
    {
      std::cout << "SLAVE:   " << std::setw(14) << gssum<< std::endl;
      std::cout << "MASTER:  " << std::setw(14) << gmsum << std::endl;
      std::cout << "Balance: " << std::setw(14) << gcsum << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << ">>      Complete                          <<" << std::endl;
    }
    // slave
    lssum=0.0;
    for (int i=0;i<augfs.MyLength();++i) lssum+= augfs[i];
    Comm().SumAll(&lssum,&gssum,1);
    // master
    lmsum=0.0;
    for (int i=0;i<augfm.MyLength();++i) lmsum+= augfm[i];
    Comm().SumAll(&lmsum,&gmsum,1);
    // complete balance check
    gcsum = gssum+gmsum;
    if (abs(gcsum)>1.0e-11) dserror("Conservation of linear momentum is not fulfilled!");
    if (Comm().MyPID()==0)
    {
      std::cout << "SLAVE:   " << std::setw(14) << gssum<< std::endl;
      std::cout << "MASTER:  " << std::setw(14) << gmsum << std::endl;
      std::cout << "Balance: " << std::setw(14) << gcsum << std::endl;
    }
  }
  /*-------------------------------*
   | ANGULAR MOMENTUM CONSERVATION |
   *-------------------------------*/
  if (Data().PrintAngularMomConservation())
  {
    if (Comm().MyPID()==0)
    {
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << ">>      Angular Momentum Conservation     <<" << std::endl;
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
    }
    for (std::size_t i=0; i<interface_.size(); ++i)
    {
      if (Comm().MyPID()==0)
      {
        std::cout << ">>----- Interface " << std::setw(2) << i;
        std::cout << " ---------------------<<" << std::endl;
        std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      }
      interface_[i]->EvalResultantMoment(*augfs_lm,*augfm_lm);
      if (Comm().MyPID()==0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
      }
      interface_[i]->EvalResultantMoment(*augfs_g,*augfm_g);
      if (Comm().MyPID()==0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Complete                          <<" << std::endl;
      }
      interface_[i]->EvalResultantMoment(augfs,augfm);
    }
    if (Comm().MyPID()==0)
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AugForces(
    Epetra_Vector& gAugFs_lm,
    Epetra_Vector& gAugFs_g,
    Epetra_Vector& gAugFm_lm,
    Epetra_Vector& gAugFm_g) const
{
  if (!IsInContact()) return;

  /****************** SLAVE SIDE ****************************************/
  // *** standard Lagrange multiplier fraction ***
  Teuchos::RCP<Epetra_Vector> augfs_lm =
      Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Teuchos::RCP<Epetra_Vector> z_exp =
      Teuchos::rcp(new Epetra_Vector(*Data().GActiveNDofRowMapPtr()));
  LINALG::Export(*z_,*z_exp);
  /* Caution: We have to use the old D/M matrices, since the current ones
   *          are not redistributed! */
  Data().OldDMatrixPtr()->Multiply(true,*z_exp,*augfs_lm);

  // Export
  LINALG::Export(*augfs_lm,gAugFs_lm);

  // *** regularization fraction ***
  Teuchos::RCP<Epetra_Vector> augfs_g =
      Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Teuchos::RCP<Epetra_Vector> aWGap = Teuchos::rcp(
      new Epetra_Vector(*Data().GActiveNDofRowMapPtr(),true));
  LINALG::Export(*Data().AWGapPtr(),*aWGap);
  MultiplyElementwise(*Data().CnPtr(),*Data().GActiveNodeRowMapPtr(),
      *aWGap,false);
  aWGap->Scale(-1.0);

  Data().OldDMatrixPtr()->Multiply(true,*aWGap,*augfs_g);
  // Export
  LINALG::Export(*augfs_g,gAugFs_g);

  /****************** MASTER SIDE ***************************************/
  // *** standard lagrange multiplier fraction ***
  Teuchos::RCP<Epetra_Vector> augfm_lm =
      Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)));
  Data().OldMMatrixPtr()->Multiply(true,*z_exp,*augfm_lm);
  // Export
  LINALG::Export(*augfm_lm,gAugFm_lm);

  // *** regularization fraction ***
  Teuchos::RCP<Epetra_Vector> augfm_g =
      Teuchos::rcp(new Epetra_Vector(MaDoFRowMap(true)));
  Data().OldMMatrixPtr()->Multiply(true,*aWGap,*augfm_g);
  // Export
  LINALG::Export(*augfm_g,gAugFm_g);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::OutputStresses()
{
  // reset contact stress class variables
  Data().StressNormalPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Data().StressTangentialPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // get nodal normal and tangential directions
      double* nn  = cnode->MoData().n();
      double* nt1 = cnode->CoData().txi();
      double* nt2 = cnode->CoData().teta();
      double lmn  = cnode->MoData().lm()[0];
      double lmt1 = cnode->MoData().lm()[1];
      double lmt2 = cnode->MoData().lm()[2];

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(dim);

      // normal stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (Data().StressNormalPtr()->Map()).LID(cnode->Dofs()[dof]);
        if (DRT::INPUT::IntegralValue<int>(Params(),"LM_NODAL_SCALE")==false
            || cnode->MoData().GetScale()==0.)
          (*Data().StressNormalPtr())[locindex[dof]] = -lmn*nn[dof];
        else
          (*Data().StressNormalPtr())[locindex[dof]] = -lmn*nn[dof]/cnode->MoData().GetScale();
      }

      // tangential stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (Data().StressTangentialPtr()->Map()).LID(cnode->Dofs()[dof]);
        if (DRT::INPUT::IntegralValue<int>(Params(),"LM_NODAL_SCALE")==false
            || cnode->MoData().GetScale()==0.)
          (*Data().StressTangentialPtr())[locindex[dof]] = -lmt1*nt1[dof]-lmt2*nt2[dof];
        else
          (*Data().StressTangentialPtr())[locindex[dof]] = -lmt1*nt1[dof]-lmt2*nt2[dof]/cnode->MoData().GetScale();
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AugmentedLagrangeStrategy::
    GetRhsBlockPtr(const enum STR::VecBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  Teuchos::RCP<const Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case STR::block_displ:
    {
      vec_ptr = Data().StrContactRhsPtr();
      break;
    }
    case STR::block_constraint:
    {
      vec_ptr = Data().ConstrRhsPtr();
      break;
    }
    default:
    {
      dserror("Unknown STR::VecBlockType!");
      break;
    }
  }

  return vec_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> CONTACT::AugmentedLagrangeStrategy::
    GetMatrixBlockPtr(const enum STR::MatBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case STR::block_displ_displ:
    {
      mat_ptr = Teuchos::rcp(
          new LINALG::SparseMatrix(SlMaDoFRowMap(true),100,false,true));

      // build matrix kdd
      mat_ptr->Add(*Data().DGLmSlLinMatrixPtr(),false,-1.0,1.0);
      mat_ptr->Add(*Data().DGLmMaLinMatrixPtr(),false,-1.0,1.0);
      mat_ptr->Add(*Data().DGGSlLinMatrixPtr(),false,-1.0,1.0);
      mat_ptr->Add(*Data().DGGMaLinMatrixPtr(),false,-1.0,1.0);
      mat_ptr->Complete(SlMaDoFRowMap(true),SlMaDoFRowMap(true));

      // transform parallel row/column distribution of matrix kdd
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            SlMaDoFRowMapPtr(false),SlMaDoFRowMapPtr(false));

      break;
    }
    case STR::block_displ_lm:
    {
      // build constraint matrix kdz
      Teuchos::RCP<LINALG::SparseMatrix> kdz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(*Data().GDispDofRowMapPtr(),100,false,true));
      kdz_ptr->Add(*Data().DMatrixPtr(),true,-1.0,1.0);
      kdz_ptr->Add(*Data().MMatrixPtr(),true,-1.0,1.0);
      kdz_ptr->Complete(SlDoFRowMap(true),*Data().GDispDofRowMapPtr());

      // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
      mat_ptr = MORTAR::MatrixColTransformGIDs(kdz_ptr,LMDoFRowMapPtr(true));

      // transform parallel row/column distribution of matrix kdz
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            ProblemDofs(),LMDoFRowMapPtr(false));

      break;
    }
    case STR::block_lm_displ:
    {
      // build constraint matrix kzd
      Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs =
          LINALG::SplitMap(SlDoFRowMap(true),*Data().GActiveDofRowMapPtr());
      Teuchos::RCP<LINALG::SparseMatrix> kzd_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true),100,false,true));
      kzd_ptr->Add(*Data().DLmNWGapLinMatrixPtr(),false,1.0,1.0);
      kzd_ptr->Add(*Data().DLmTLmTLinMatrixPtr(),false,1.0,1.0);
      kzd_ptr->Add(*Data().InactiveLinMatrixPtr(),false,1.0,1.0);
      kzd_ptr->Complete(*Data().GDispDofRowMapPtr(),SlDoFRowMap(true));

      // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
      mat_ptr = MORTAR::MatrixRowTransformGIDs(kzd_ptr,LMDoFRowMapPtr(true));

      // transform parallel row/column distribution of matrix kzd
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            LMDoFRowMapPtr(false),ProblemDofs());

      break;
    }
    case STR::block_lm_lm:
    {
      // build constraint matrix kzz
      Teuchos::RCP<LINALG::SparseMatrix> kzz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true),100,false,true));
      kzz_ptr->Add(*Data().InactiveMatrixPtr(),false,1.0,1.0);
      kzz_ptr->Add(*Data().DLmTLmTMatrixPtr(),false,1.0,1.0);
      kzz_ptr->Complete(SlDoFRowMap(true),SlDoFRowMap(true));

      // transform constraint matrix kzz to lmdofmap
      mat_ptr = MORTAR::MatrixRowColTransformGIDs(kzz_ptr,
           LMDoFRowMapPtr(true),LMDoFRowMapPtr(true));

      // transform parallel row/column distribution of matrix kzz
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            LMDoFRowMapPtr(false),LMDoFRowMapPtr(false));

      break;
    }
    default:
    {
      dserror("Unknown STR::MatBlockType!");
      break;
    }
  }

  return mat_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CONTACT::AugmentedLagrangeStrategy::ConstraintNorm() const
{
  double nrm2 = 0.0;
  Data().ConstrRhsPtr()->Norm2(&nrm2);
  return nrm2;
};
