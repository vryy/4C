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
#include "contact_aug_potential.H"
#include "contact_aug_steepest_ascent_strategy.H"

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
CONTACT::AUG::DataContainer::DataContainer()
    : gFdCheck_( false ),
      wasincontactlastiter_( false ),
      isactivesetconverged_( false ),
      printlinearconservation_( false ),
      printangularconservation_( false ),
      is_semi_smooth_newton_( false ),
      matrix_maps_valid_( false ),
      vector_maps_valid_( false ),
      cn_( -1.0 ),
      parallel_strategy_( INPAR::MORTAR::ghosting_redundant ),
      potentialPtr_( Teuchos::null ),
      dGLmSlLinMatrixPtr_( Teuchos::null ),
      dGLmMaLinMatrixPtr_( Teuchos::null ),
      dGGSlLinMatrixPtr_( Teuchos::null ),
      dGGMaLinMatrixPtr_( Teuchos::null ),
      dLmNWGapLinMatrixPtr_( Teuchos::null ),
      dLmTLmTMatrixPtr_( Teuchos::null ),
      dLmTLmTLinMatrixPtr_( Teuchos::null ),
      inactiveLinMatrixPtr_( Teuchos::null ),
      inactiveDiagMatrixPtr_( Teuchos::null ),
      aPtr_( Teuchos::null ),
      kappaPtr_( Teuchos::null ),
      lmNPtr_( Teuchos::null ),
      dLmTLmTRhsPtr_( Teuchos::null ),
      slForceLmPtr_( Teuchos::null ),
      slForceGPtr_( Teuchos::null ),
      maForceLmPtr_( Teuchos::null ),
      maForceGPtr_( Teuchos::null ),
      cnPtr_( Teuchos::null ),
      gsndofrowmapPtr_( Teuchos::null ),
      gstdofrowmapPtr_( Teuchos::null ),
      gOldActiveSlaveNodesPtr_( Teuchos::null )
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::DataContainer::InitSubDataContainer(
    const INPAR::CONTACT::SolvingStrategy strat_type )
{
  switch ( strat_type )
  {
    case INPAR::CONTACT::solution_steepest_ascent:
      sa_data_ptr_ = Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::DataContainer>(
          new CONTACT::AUG::STEEPESTASCENT::DataContainer() );
      break;
    default:
      dserror( "There is no known sub-data container for the given "
          "strategy type! ( strat_type = %s )",
          INPAR::CONTACT::SolvingStrategy2String( strat_type ).c_str() );
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AUG::Strategy::Strategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params,
    const plain_interface_set& interfaces,
    int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm,
    int maxdof)
    : ::CONTACT::CoAbstractStrategy( data_ptr, DofRowMap, NodeRowMap, params,
        dim, comm, 0.0, maxdof ),
        augDataPtr_( Teuchos::rcp_dynamic_cast<CONTACT::AUG::DataContainer>(data_ptr, true) ),
        augData_( *augDataPtr_ )
{
  // store values of the parameter list
  const Teuchos::ParameterList& p_augmented = params.sublist( "AUGMENTED" );

  Data().PrintLinearMomConservation() =
      DRT::INPUT::IntegralValue<bool>( p_augmented,
          "PRINT_LINEAR_CONSERVATION");

  Data().PrintAngularMomConservation() =
      DRT::INPUT::IntegralValue<bool>( p_augmented,
          "PRINT_ANGULAR_CONSERVATION");

  Data().SetIsSemiSmoothNewton( DRT::INPUT::IntegralValue<bool>( params,
      "SEMI_SMOOTH_NEWTON") );

  Data().SetConstantCn( Params().get<double>( "SEMI_SMOOTH_CN" ) );

  Data().SetParallelStrategy( DRT::INPUT::IntegralValue<
      INPAR::MORTAR::ParallelStrategy>( Params(), "PARALLEL_STRATEGY") );

  Data().SetPotential( Teuchos::rcp( new Potential( *this, *augDataPtr_ ) ) );

  interface_.reserve( interfaces.size() );
  for ( plain_interface_set::const_iterator cit = interfaces.begin();
        cit != interfaces.end(); ++cit )
  {
    const Teuchos::RCP<CONTACT::CoInterface> & interface = *cit;
    // cast to augmented interfaces just as sanity check
    interface_.push_back( Teuchos::rcp_dynamic_cast<CONTACT::AUG::Interface>(
        interface, true ) );
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::Strategy::IsSaddlePointSystem() const
{
  return (IsInContact() or WasInContact() or WasInContactLastTimeStep());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AUG::plain_interface_set& CONTACT::AUG::Strategy::Interfaces()
{
  return interface_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const CONTACT::AUG::plain_interface_set&
CONTACT::AUG::Strategy::Interfaces() const
{
  return interface_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PostSetup(
    bool redistributed,
    bool init)
{
  if ( init )
  {
    // re-setup the global sndofrowmap_ and stdofrowmap_ for the augmented Lagrangian case
    AssembleGlobalSlNTDofRowMaps();

    // initialize cn
    InitializeCn( false, true );
  }

  // just used for the redistributed case
  if ( redistributed )
  {
    // reassemble the global slave normal/tangential dof row maps
    AssembleGlobalSlNTDofRowMaps();

    // redistribute the cn-vector
    InitializeCn( redistributed, false );

    // redistribute the global augmented old active slave nodes map
    if ( ( not Data().GOldActiveSlaveNodesPtr().is_null() ) and
         ( Data().GOldActiveSlaveNodes().NumGlobalElements() > 0 ) )
      RedistributeRowMap( SlRowNodes(), Data().GOldActiveSlaveNodes() );
  }

  // in both cases the maps change and we have to re-build all matrices
  Data().SetMatrixMapsValid( false );
  Data().SetVectorMapsValid( false );

  // setup the potential class with the current maps
  Data().Potential().Setup();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleGlobalSlNTDofRowMaps()
{
  Data().GSlNormalDofRowMapPtr()     = Teuchos::null;
  Data().GSlTangentialDofRowMapPtr() = Teuchos::null;

  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    const CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    Data().GSlNormalDofRowMapPtr()     =
        LINALG::MergeMap(Data().GSlNormalDofRowMapPtr(),interface.SlaveRowNDofs());
    Data().GSlTangentialDofRowMapPtr() =
        LINALG::MergeMap(Data().GSlTangentialDofRowMapPtr(),interface.SlaveRowTDofs());
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::InitializeCn(
    bool redistributed,
    bool init)
{
  if (init)
  {
    // get cn from the input file
    double cn = Data().ConstantdCn();

    if (Data().CnPtr().is_null() or Data().Cn().GlobalLength()==0)
      Data().CnPtr() = LINALG::CreateVector(SlRowNodes(),true);
    // set all nodal cn-values to the input value
    Data().Cn().PutScalar(cn);
  }

  if (redistributed)
  {
    // redistribute the cn-vector
    Teuchos::RCP<Epetra_Vector> newcn = Teuchos::rcp(new Epetra_Vector(SlRowNodes()));
    LINALG::Export(Data().Cn(),*newcn);
    Data().CnPtr() = newcn;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::DoReadRestart(
    IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis,
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr )
{
  CONTACT::CoAbstractStrategy::DoReadRestart( reader, dis, cparams_ptr );
  PostSetup( false, false );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::InitMortar()
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
void CONTACT::AUG::Strategy::AssembleMortar()
{
  // for all interfaces
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    const CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    interface.AssembleAugDnMnMatrix(Data().DMatrix(),Data().MMatrix());
  }

  Data().DMatrix().Complete(SlDoFRowMap(true),SlRowNodes());
  Data().MMatrix().Complete(MaDoFRowMap(true),SlRowNodes());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::SplitMortar()
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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ZeroizeMatrices()
{
  // *** zeroize existing global matrices ***
  Data().DGLmSlLinMatrix().Zero();
  Data().DGLmMaLinMatrix().Reset();
  Data().DGGSlLinMatrix().Zero();
  Data().DGGMaLinMatrix().Reset();

  Data().DLmNWGapLinMatrix().Reset();
  Data().DLmTLmTMatrix().Zero();
  Data().DLmTLmTLinMatrix().Zero();

  Data().InactiveLinMatrix().Zero();
  Data().InactiveDiagMatrix().PutScalar( 0.0 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CreateMatrices( const Epetra_Map &
    gAugInactiveSlaveDofs )
{
  // *** (re)setup global matrices ***
  Data().DGLmSlLinMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      SlDoFRowMap(true),100,false,true,LINALG::SparseMatrix::FE_MATRIX));
  Data().DGLmMaLinMatrixPtr() = Teuchos::rcp(new LINALG::SparseMatrix(
      MaDoFRowMap(true),100,false,false,LINALG::SparseMatrix::FE_MATRIX));
  Data().DGGSlLinMatrixPtr()  = Teuchos::rcp(new LINALG::SparseMatrix(
      SlDoFRowMap(true),100,false,true,LINALG::SparseMatrix::FE_MATRIX));
  Data().DGGMaLinMatrixPtr()  = Teuchos::rcp(new LINALG::SparseMatrix(
      MaDoFRowMap(true),100,false,false,LINALG::SparseMatrix::FE_MATRIX));

  Data().DLmNWGapLinMatrixPtr() = Teuchos::rcp(
      new LINALG::SparseMatrix(Data().GActiveNDofRowMap(),100,false,false));
  Data().DLmTLmTMatrixPtr()     = Teuchos::rcp(
      new LINALG::SparseMatrix(Data().GActiveTDofRowMap(),100,false,true));
  Data().DLmTLmTLinMatrixPtr()  = Teuchos::rcp(
      new LINALG::SparseMatrix(Data().GActiveTDofRowMap(),100,false,true));

  Data().InactiveDiagMatrixPtr() =
      LINALG::CreateVector( gAugInactiveSlaveDofs, true );
  Data().InactiveLinMatrixPtr() =
      Teuchos::rcp(new LINALG::SparseMatrix(gAugInactiveSlaveDofs,100,false,true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ZeroizeVectors()
{
  // *** zeroize existing global matrices ***
  Data().LmN().PutScalar( 0.0 );
  Data().AWGap().PutScalar( 0.0 );

  Data().DLmTLmTRhs().PutScalar( 0.0 );
  Data().InactiveRhs().PutScalar( 0.0 );

  Data().AVec().PutScalar( 0.0 );
  Data().KappaVec().PutScalar( 0.0 );
  Data().WGap().PutScalar( 0.0 );

  Data().SlForceLm().PutScalar( 0.0 );
  Data().SlForceG().PutScalar( 0.0 );
  Data().MaForceLm().PutScalar( 0.0 );
  Data().MaForceG().PutScalar( 0.0 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CreateVectors( const Epetra_Map &
    gAugInactiveSlaveDofs )
{
  // *** (re)setup global augmented Epetra_Vectors ***
  Data().LmNPtr()         = Teuchos::rcp( new Epetra_Vector(
      Data().GActiveNDofRowMap(), true ) );
  Data().AWGapPtr()       = Teuchos::rcp( new Epetra_Vector(
      Data().GActiveNDofRowMap(), true ) );
  Data().DLmTLmTRhsPtr()  = Teuchos::rcp( new Epetra_Vector(
      Data().GActiveTDofRowMap(), true ) );
  Data().InactiveRhsPtr() = Teuchos::rcp( new Epetra_Vector(
      gAugInactiveSlaveDofs, true ) );

  Data().AVecPtr()     = Teuchos::rcp( new Epetra_Vector(
      SlRowNodes(), true ) );
  Data().KappaVecPtr() = Teuchos::rcp( new Epetra_Vector(
      Data().GActiveNodeRowMap(), true ) );
  Data().WGapPtr()     = Teuchos::rcp( new Epetra_Vector(
      Data().GActiveNDofRowMap(), true ) );

  Data().SlForceLmPtr() = Teuchos::rcp( new Epetra_Vector(
      SlDoFRowMap(true), true ) );
  Data().SlForceGPtr()  = Teuchos::rcp( new Epetra_Vector(
      SlDoFRowMap(true) ), true );
  Data().MaForceLmPtr() = Teuchos::rcp( new Epetra_Vector(
      MaDoFRowMap(true), true ) );
  Data().MaForceGPtr()  = Teuchos::rcp( new Epetra_Vector(
      MaDoFRowMap(true) ),true );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::Initialize( enum MORTAR::ActionType actiontype )
{
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = Teuchos::null;
      LINALG::SplitMap(SlDoFRowMap(true), Data().GActiveDofRowMap());

  switch ( actiontype )
  {
    case MORTAR::eval_force_stiff:
    {
      if ( Data().MatrixMapsValid() )
      {
        ZeroizeMatrices();
      }
      else
      {
        gAugInactiveSlaveDofs = LINALG::SplitMap( SlDoFRowMap(true),
            Data().GActiveDofRowMap() );
        CreateMatrices( *gAugInactiveSlaveDofs );
        Data().SetMatrixMapsValid( true );
      }
    }
    case MORTAR::eval_force:
    {
      if ( Data().VectorMapsValid() )
      {
        ZeroizeVectors();
      }
      else
      {
        if ( gAugInactiveSlaveDofs.is_null() )
          gAugInactiveSlaveDofs = LINALG::SplitMap( SlDoFRowMap(true),
                      Data().GActiveDofRowMap() );
        CreateVectors( *gAugInactiveSlaveDofs );
        Data().SetVectorMapsValid( true );
      }

      break;
    }
    default:
    {
      dserror("Unsupported action type detected: %s",
          MORTAR::ActionType2String( actiontype ).c_str() );
      exit( EXIT_FAILURE );
    }
  }

  if (Data().IsFriction())
    dserror("AugmentedLagrangeStrategy::Initialize: "
        "Frictional case is not yet considered!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalForceStiff(
    CONTACT::ParamsInterface& cparams)
{
  // call the evaluate force routine
  EvalForce(cparams);

  // --- Assemble stiffness matrix ---------------------------------------
  AssembleContactStiff();

  PostEvalForceStiff( cparams );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PostEvalForceStiff(
    CONTACT::ParamsInterface& cparams )
{
  // --- DEBUGGING -------------------------------------------------------
  // Finite Difference check at Gauss-point level
  AugFDCheckGP(cparams);

  // Finite Difference check at global level
  AugFDCheckGlobal(cparams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::InitEvalInterface(
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // get type of parallel strategy
  INPAR::MORTAR::ParallelStrategy strat = Data().ParallelStrategy();

  // Evaluation for all interfaces
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    // initialize / reset interfaces
    interface.Initialize();

    //store required integration time
    Data().IntTime() += interface.Inttime();
    switch (strat)
    {
      /*----------------------------------------------------------*
       |  Fully redundant ghosting of master side                 |
       *----------------------------------------------------------*/
      case INPAR::MORTAR::ghosting_redundant:
      {
        // evaluate averaged weighted gap
        interface.Evaluate(0,cparams_ptr);
        // Calculate weighted gap
        interface.WGap();
        // evaluate remaining entities and linearization
        interface.RedEvaluate(cparams_ptr);
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
void CONTACT::AUG::Strategy::RecoverState(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  // set the new search direction in the potential object
  Data().Potential().SetDirection( dir );

  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we use this routine just to store the Lagrange
   * multiplier increment. */
  Epetra_Vector zincr_exp( LMDoFRowMap( true ) );
  LINALG::Export( dir, zincr_exp );

  // get the current step length
  const double stepLength = cparams.GetStepLength();
  // ---------------------------------------------------------------------
  /* store the SCALED Lagrange multiplier increment in the contact
   * strategy */
  // ---------------------------------------------------------------------
  Data().LmIncrPtr()->Update( stepLength, zincr_exp, 0.0 );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::ResetLagrangeMultipliers(
    const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we do not have to check if it is a saddle
   * point system. */
  Epetra_Vector& z = *Data().LmPtr();

  z.PutScalar( 0.0 );

  if ( z.ReplaceMap( LMDoFRowMap( true ) ) )
    dserror( "ReplaceMap failed!" );

  LINALG::Export( xnew, z );

  if ( z.ReplaceMap( SlDoFRowMap( true ) ) )
    dserror( "ReplaceMap failed!" );

  // ---------------------------------------------------------------------
  // store the new Lagrange multiplier in the nodes
  // ---------------------------------------------------------------------
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::PreEvalForce(
    CONTACT::ParamsInterface& cparams )
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
  Initialize( cparams.GetActionType() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalForce(
    CONTACT::ParamsInterface& cparams)
{
  // --- Prepare the evaluation and integrate all quantities ------------------
  PreEvalForce( cparams );

  // --- Assemble the gap vectors ---------------------------------------------
  AssembleGap();

  // --- compute the augmented forces -----------------------------------------
  EvalAugmentedForces();

  // --- Assemble the ride hand side terms ------------------------------------
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
void CONTACT::AUG::Strategy::AssembleGap()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;
  /*--------------------------------------------------------------------*
   | Assembly                                                          |
   *--------------------------------------------------------------------*
   | --> weighted gap                                                   |
   | --> averaged weighted gap                                          |
   *--------------------------------------------------------------------*/
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    interface.AssembleGapVectors( Data().AWGap(), Data().WGap() );
  }

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::Strategy::AssembleContactRHS()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return false;
  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE CONTACT RIGHT HAND SIDE                   |
   *--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*
   | Assembly                                                           |
   *--------------------------------------------------------------------*
   | --> normal Lagrange multiplier                                     |
   | --> (averaged) weighted gap                                        |
   | --> tangential constraint right hand side for the frictionless case|
   | --> normal and tangential inactive rhs                             |
   *--------------------------------------------------------------------*/
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    // --- augmented Lagrange formulation --------------------------------
    // --- FORCE BALANCE -------------------------------------------------
    interface.AssembleLmNVector(Data().LmN());

    // --- CONSTRAINTS ---------------------------------------------------
    // active - normal direction
    // --> wGapRhs_

    // active - tangential direction
    interface.AssembleDLmTLmTRhs(Data().DLmTLmTRhs());

    // inactive - all directions
    interface.AssembleAugInactiveRhs(Data().InactiveRhs(),Data().Cn());
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AssembleContactStiff()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;
  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE TANGENTIAL STIFFNESS MATRIX               |
   *--------------------------------------------------------------------*/
  // --- augmented Lagrange formulation ----------------------------------
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    const CONTACT::AUG::Interface& interface =
        dynamic_cast<CONTACT::AUG::Interface&>( **cit );

    // Calculate averaged weighted gap linearization for all active nodes
    interface.AWGapLin();
    // --- Force Balance ------------------------------------------------
    // linearization w.r.t. displ.
    interface.AssembleDGLmLinMatrix(Data().DGLmSlLinMatrix(),
                                    Data().DGLmMaLinMatrix());
    interface.AssembleDGGLinMatrix(Data().DGGSlLinMatrix(),
                                   Data().DGGMaLinMatrix(),
                                   Data().Cn());
    // --- Constraints --------------------------------------------------
    // linearization w.r.t. LM
    interface.AssembleDLmTLmTMatrix(Data().DLmTLmTMatrix());
    interface.AssembleAugInactiveDiagMatrix(Data().InactiveDiagMatrix(),Data().Cn());
    // linearization w.r.t. displ.
    // active - normal direction
    interface.AssembleDLmNWGapLinMatrix(Data().DLmNWGapLinMatrix());
    // active - tangential direction
    interface.AssembleDLmTLmTLinMatrix(Data().DLmTLmTLinMatrix());
    /*--------------------------------------------------------------------*
     | The linearization of the nodal area w.r.t. the displ. for inactive |
     | nodes can help to prevent oscillations of the active set, because  |
     | the reduction of the inactive lm-value is decelerated.             |
     |                                                                    |
     | In general, the following relation does NOT hold anymore:          |
     |                    Delta(z_{n,i}^{k+1}) = - z_{n,i}^{k}            |
     |                          z_{n,i}^{k+1}  =   0                      |
     *--------------------------------------------------------------------*/
    interface.AssembleAugInactiveLinMatrix(Data().InactiveLinMatrix(),Data().Cn());
  }

  // --- START - FillComplete matrices ----------------------------------

  // domainmap: columnmap | rangemap: rowmap
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs =
      LINALG::SplitMap( SlDoFRowMap(true), Data().GActiveDofRowMap() );

  // --- Force Balance --------------------------------------------------
  // linearization w.r.t. displ.
  Data().DGLmSlLinMatrix().Complete(SlMaDoFRowMap(true),SlDoFRowMap(true));
  Data().DGLmMaLinMatrix().Complete(SlMaDoFRowMap(true),MaDoFRowMap(true));
  Data().DGGSlLinMatrix().Complete(SlMaDoFRowMap(true),SlDoFRowMap(true));
  Data().DGGMaLinMatrix().Complete(SlMaDoFRowMap(true),MaDoFRowMap(true));

  // --- Constraints ----------------------------------------------------
  // linearization w.r.t. LM
  Data().DLmTLmTMatrix().Complete(Data().GActiveTDofRowMap(),Data().GActiveTDofRowMap());

  // linearization w.r.t. displ.
  Data().DLmNWGapLinMatrix().Complete(SlMaDoFRowMap(true),Data().GActiveNDofRowMap());
  Data().DLmTLmTLinMatrix().Complete(SlDoFRowMap(true),Data().GActiveTDofRowMap());
  Data().InactiveLinMatrix().Complete(SlDoFRowMap(true),*gAugInactiveSlaveDofs);

  // --- END - FillComplete matrices ------------------------------------

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::UpdateActiveSetSemiSmooth(
    const bool& correction)
{
  // get out of here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  if ( not Data().IsSemiSmoothNewton() )
  {
    // loop over all interfaces
    for ( plain_interface_set::const_iterator cit = interface_.begin();
          cit != interface_.end(); ++cit )
    {
      const CoInterface& interface = **cit;

      // loop over all slave nodes on the current interface
      for (int j=0;j<interface.SlaveRowNodes()->NumMyElements();++j)
      {
        int gid = interface.SlaveRowNodes()->GID(j);
        DRT::Node* node = interface.Discret().gNode(gid);
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
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    const CoInterface& interface = **cit;

    // loop over all slave nodes of the current interface
    const int num_my_slave_row_nodes = interface.SlaveRowNodes()->NumMyElements();
    int* my_slave_row_node_gids       = interface.SlaveRowNodes()->MyGlobalElements();
    for ( int j=0; j < num_my_slave_row_nodes; ++j )
    {
      const int gid = my_slave_row_node_gids[j];
      CoNode* cnode = dynamic_cast<CoNode*>(interface.Discret().gNode(gid));
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
      double kappa = cnode->AugData().GetKappa();
      double awgap = cnode->AugData().GetWGap();
      // TODO
//      std::cout << cnode->Id() << " | " << cnode->X()[0] << ", " << cnode->X()[1] <<
//          " | " << cnode->AugData().GetWGap() << std::endl;
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
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    CoInterface& interface = **cit;

    // update active set Epetra_Maps
    interface.BuildActiveSet();

    // Update Active set
    Data().GActiveNodeRowMapPtr() = LINALG::MergeMap(
        Data().GActiveNodeRowMapPtr(),interface.ActiveNodes(),false);
    Data().GActiveDofRowMapPtr()  = LINALG::MergeMap(
        Data().GActiveDofRowMapPtr(),interface.ActiveDofs(),false);
    Data().GActiveNDofRowMapPtr() = LINALG::MergeMap(
        Data().GActiveNDofRowMapPtr(),interface.ActiveNDofs(),false);
    Data().GActiveTDofRowMapPtr() = LINALG::MergeMap(
        Data().GActiveTDofRowMapPtr(),interface.ActiveTDofs(),false);
  }

  // check the convergence of the active set
  if ( Data().GActiveNodeRowMap().SameAs(Data().GOldActiveSlaveNodes()) )
  {
    Data().IsActiveSetConverged() = true;
  }
  else
  {
    Data().IsActiveSetConverged() = false;
    Data().SetVectorMapsValid( false );
    Data().SetMatrixMapsValid( false );

    // reset the active/inactive state vectors
    Data().Potential().Setup( true );
  }

  // set the new active/inactive state
  Data().Potential().SetActiveInactiveState();

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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalStrContactRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    Data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  Data().StrContactRhsPtr() =
      Teuchos::rcp(new Epetra_Vector(*ProblemDofs(),true));


  // For self contact, slave and master sets may have changed,
  if (IsSelfContact())
    dserror("ERROR: Augmented Lagrange Formulation: Self contact is not yet considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs( Data().SlForceLm() );
  augfs.Update( -1.0, Data().SlForceG(), 1.0 );

  Epetra_Vector augfs_exp( *ProblemDofs() );
  LINALG::Export( augfs, augfs_exp );
  Data().StrContactRhs().Scale( -1.0, augfs_exp );

  // Master side
  Epetra_Vector augfm( Data().MaForceLm() );
  augfm.Update( -1.0, Data().MaForceG(), 1.0 );

  Epetra_Vector augfm_exp( *ProblemDofs() );
  LINALG::Export( augfm, augfm_exp );
  Data().StrContactRhs().Update( -1.0, augfm_exp, 1.0 );

  // Check linear and angular momentum conservation
  CheckConservationLaws( augfs, augfm );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalConstrRHS()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
  {
    // (re)setup the vector
    Data().ConstrRhsPtr() = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> augConstrRhs =
      Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true),true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  AddContributionsToConstrRHS( *augConstrRhs );

  // replace row map
  augConstrRhs->ReplaceMap( LMDoFRowMap(true) );

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
void CONTACT::AUG::Strategy::AddContributionsToConstrRHS(
    Epetra_Vector& augConstrRhs ) const
{
  // Add active constraints in normal direction:
  LINALG::AssembleMyVector( 0.0, augConstrRhs, 1.0, *Data().WGapPtr() );

  // Add inactive constraints
  LINALG::AssembleMyVector( 1.0, augConstrRhs, 1.0, *Data().InactiveRhsPtr() );

  // Add tangential frictionless constraints
  LINALG::AssembleMyVector( 1.0, augConstrRhs, 1.0, *Data().DLmTLmTRhsPtr() );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::CheckConservationLaws(
    const Epetra_Vector& augfs,
    const Epetra_Vector& augfm)
{
  if (not Data().PrintLinearMomConservation() and
      not Data().PrintAngularMomConservation())
    return;
  /****************** SLAVE SIDE ********************************************************/
  // standard Lagrange multiplier fraction
  Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Teuchos::RCP<Epetra_Vector> z_exp = LINALG::ExtractMyVector( *z_, Data().GActiveNDofRowMap() );
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

    unsigned count = 0;
    for ( plain_interface_set::const_iterator cit = interface_.begin();
          cit != interface_.end(); ++cit )
    {
      const CoInterface& interface = **cit;

      if (Comm().MyPID()==0)
      {
        std::cout << ">>----- Interface " << std::setw(2) << count++;
        std::cout << " ---------------------<<" << std::endl;
        std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      }
      interface.EvalResultantMoment(*augfs_lm,*augfm_lm);
      if (Comm().MyPID()==0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
      }
      interface.EvalResultantMoment(*augfs_g,*augfm_g);
      if (Comm().MyPID()==0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Complete                          <<" << std::endl;
      }
      interface.EvalResultantMoment(augfs,augfm);
    }
    if (Comm().MyPID()==0)
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::EvalAugmentedForces()
{
  if (!IsInContact() and !WasInContact() and !WasInContactLastTimeStep())
    return;

  // augmented force terms
  enum ForceTerm { lm_term = 0, g_term = 1 };

  const Epetra_Vector& zn = Data().Potential().GetZnActive();
  Epetra_Vector awgapn( Data().GActiveNDofRowMap(), true );
  LINALG::ExtractMyVector( Data().AWGap(), awgapn );

  MultiplyElementwise( Data().Cn(), Data().GActiveNodeRowMap(),
      awgapn, false );

  double* values[2] = { NULL, NULL };

  values[ lm_term ] = zn.Values();
  values[ g_term ]  = awgapn.Values();

  const Epetra_MultiVector zn_awgapn( View, Data().GActiveNDofRowMap(), values, 2 );

  /****************** SLAVE SIDE ****************************************/
  double* f_values[2] = { NULL, NULL };
  f_values[ lm_term ] = Data().SlForceLm().Values();
  f_values[ g_term ]  = Data().SlForceG().Values();

  // interface forces on the slave side
  Epetra_MultiVector slForces( View, SlDoFRowMap( true ), f_values, 2 );

  Data().DMatrix().Multiply( true, zn_awgapn , slForces );

  /****************** MASTER SIDE ****************************************/
  f_values[ lm_term ] = Data().MaForceLm().Values();
  f_values[ g_term ]  = Data().MaForceG().Values();

  // interface forces on the slave side
  Epetra_MultiVector maForces( View, MaDoFRowMap( true ), f_values, 2 );

  Data().MMatrix().Multiply( true, zn_awgapn , maForces );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AugForces(
    Epetra_Vector& gAugFs_lm,
    Epetra_Vector& gAugFs_g,
    Epetra_Vector& gAugFm_lm,
    Epetra_Vector& gAugFm_g) const
{
  if (!IsInContact())
    return;

  /****************** SLAVE SIDE ****************************************/
  // *** standard Lagrange multiplier fraction ***
  // Export
  LINALG::Export( *Data().SlForceLmPtr(), gAugFs_lm );

  // *** regularization fraction ***
  // Export
  LINALG::Export( *Data().SlForceGPtr(), gAugFs_g );
  gAugFs_g.Scale( -1.0 );

  /****************** MASTER SIDE ***************************************/
  // *** standard lagrange multiplier fraction ***
  // Export
  LINALG::Export( *Data().MaForceLmPtr(), gAugFm_lm );

  // *** regularization fraction ***
  // Export
  LINALG::Export( *Data().MaForceGPtr(), gAugFm_g );
  gAugFm_g.Scale( -1.0 );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::OutputStresses()
{
  // reset contact stress class variables
  Data().StressNormalPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));
  Data().StressTangentialPtr() = Teuchos::rcp(new Epetra_Vector(SlDoFRowMap(true)));

  // loop over all interfaces
  for ( plain_interface_set::const_iterator cit = interface_.begin();
        cit != interface_.end(); ++cit )
  {
    const CoInterface& interface = **cit;

    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface.SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface.SlaveRowNodes()->GID(j);
      DRT::Node* node = interface.Discret().gNode(gid);
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
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::Strategy::
    GetRhsBlockPtr(const enum DRT::UTILS::VecBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  // get the desired vector and return it (read-only)
  Teuchos::RCP<const Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ:
    {
      vec_ptr = Data().StrContactRhsPtr();
      break;
    }
    case DRT::UTILS::block_constraint:
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
Teuchos::RCP<LINALG::SparseMatrix> CONTACT::AUG::Strategy::
    GetMatrixBlockPtr(const enum DRT::UTILS::MatBlockType& bt) const
{
  // if there are no active contact contributions
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case DRT::UTILS::block_displ_displ:
    {
      mat_ptr = Teuchos::rcp(
          new LINALG::SparseMatrix(SlMaDoFRowMap(true),100,false,true));

      // build matrix kdd
      AddContributionsToMatrixBlockDisplDispl( *mat_ptr );
      mat_ptr->Complete(SlMaDoFRowMap(true),SlMaDoFRowMap(true));

      // transform parallel row/column distribution of matrix kdd
      // (only necessary in the parallel redistribution case)
      if (ParRedist())
        mat_ptr = MORTAR::MatrixRowColTransform(mat_ptr,
            SlMaDoFRowMapPtr(false),SlMaDoFRowMapPtr(false));

      break;
    }
    case DRT::UTILS::block_displ_lm:
    {
      Teuchos::RCP<LINALG::SparseMatrix> kdz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(*Data().GDispDofRowMapPtr(),100,false,true));

      // build constraint matrix kdz
      AddContributionsToMatrixBlockDisplLm( *kdz_ptr );
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
    case DRT::UTILS::block_lm_displ:
    {
      Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs =
          LINALG::SplitMap(SlDoFRowMap(true),*Data().GActiveDofRowMapPtr());
      Teuchos::RCP<LINALG::SparseMatrix> kzd_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true),100,false,true));

      // build constraint matrix kzd
      AddContributionsToMatrixBlockLmDispl( *kzd_ptr );
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
    case DRT::UTILS::block_lm_lm:
    {
      Teuchos::RCP<LINALG::SparseMatrix> kzz_ptr =
          Teuchos::rcp(new LINALG::SparseMatrix(SlDoFRowMap(true),100,false,true));

      // build constraint matrix kzz
      AddContributionsToMatrixBlockLmLm( *kzz_ptr );
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
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockDisplDispl(
    LINALG::SparseMatrix& kdd ) const
{
  kdd.Add(*Data().DGLmSlLinMatrixPtr(),false,-1.0,1.0);
  kdd.Add(*Data().DGLmMaLinMatrixPtr(),false,-1.0,1.0);
  kdd.Add(*Data().DGGSlLinMatrixPtr(),false,-1.0,1.0);
  kdd.Add(*Data().DGGMaLinMatrixPtr(),false,-1.0,1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockDisplLm(
    LINALG::SparseMatrix& kdz ) const
{
  kdz.Add(*Data().DMatrixPtr(),true,-1.0,1.0);
  kdz.Add(*Data().MMatrixPtr(),true,-1.0,1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockLmDispl(
    LINALG::SparseMatrix& kzd ) const
{
  kzd.Add(*Data().DLmNWGapLinMatrixPtr(),false,1.0,1.0);
  kzd.Add(*Data().DLmTLmTLinMatrixPtr(),false,1.0,1.0);
  kzd.Add(*Data().InactiveLinMatrixPtr(),false,1.0,1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Strategy::AddContributionsToMatrixBlockLmLm(
    LINALG::SparseMatrix& kzz ) const
{

  if ( LINALG::InsertMyRowDiagonalIntoUnfilledMatrix( kzz,
      *Data().InactiveDiagMatrixPtr() ) )
  {
    Epetra_Vector kzz_diag = Epetra_Vector( kzz.RangeMap(), true );
    LINALG::AssembleMyVector( 0.0, kzz_diag, 1.0, *Data().InactiveDiagMatrixPtr() );

    // if the matrix is filled, we try to replace the diagonal
    if ( kzz.ReplaceDiagonalValues( kzz_diag ) )
      dserror( "ReplaceDiagonalValues failed!" );
  }

  kzz.Add( *Data().DLmTLmTMatrixPtr(), false, 1.0, 1.0 );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CONTACT::AUG::Strategy::ConstraintNorm() const
{
  double nrm2 = 0.0;
  Data().ConstrRhsPtr()->Norm2(&nrm2);
  return nrm2;
};
