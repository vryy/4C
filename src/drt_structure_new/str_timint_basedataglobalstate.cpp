/*-----------------------------------------------------------*/
/*!
\file str_timint_basedataglobalstate.cpp

\brief Global state data container for the structural (time)
       integration

\maintainer Michael Hiermeier

\date Jan 12, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_timint_basedataglobalstate.H"
#include "str_timint_basedatasdyn.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_generic.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_structure_new/str_utils.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"


#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_group_prepostoperator.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>

#include <NOX_Epetra_Vector.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState::BaseDataGlobalState()
    : isinit_(false),
      issetup_(false),
      datasdyn_(Teuchos::null),
      discret_(Teuchos::null),
      comm_(Teuchos::null),
      myRank_(-1),
      timenp_(0.0),
      timen_(Teuchos::null),
      dt_(Teuchos::null),
      stepn_(0),
      stepnp_(0),
      ispredict_(false),
      dis_(Teuchos::null),
      vel_(Teuchos::null),
      acc_(Teuchos::null),
      disnp_(Teuchos::null),
      velnp_(Teuchos::null),
      accnp_(Teuchos::null),
      fintnp_(Teuchos::null),
      fextnp_(Teuchos::null),
      freactnp_(Teuchos::null),
      finertialn_(Teuchos::null),
      finertialnp_(Teuchos::null),
      fviscon_(Teuchos::null),
      fvisconp_(Teuchos::null),
      fstructold_(Teuchos::null),
      jac_(Teuchos::null),
      mass_(Teuchos::null),
      damp_(Teuchos::null),
      timer_(Teuchos::null),
      dtsolve_(0.0),
      dtele_(0.0),
      max_block_num_(0),
      gproblem_map_ptr_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::Init(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::RCP<const BaseDataSDyn> datasdyn
    )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ----------------------------------------------------------
  // const pointer to the sDynData container
  // ----------------------------------------------------------
  datasdyn_ = datasdyn;

  // ----------------------------------------------------------
  // general purpose algorithm members
  // ----------------------------------------------------------
  {
    discret_ = discret;
    comm_ = Teuchos::rcpFromRef(discret_->Comm());
    myRank_  = comm_->MyPID();
  }

  // --------------------------------------
  // control parameters
  // --------------------------------------
  {
    timen_ = Teuchos::rcp(new ::TIMINT::TimIntMStep<double>(0, 0,
        sdynparams.get<double>("TIMEINIT")));
    dt_ = Teuchos::rcp(new ::TIMINT::TimIntMStep<double>(0, 0,
        sdynparams.get<double>("TIMESTEP")));
  }

  // end of initialization
  isinit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::Setup()
{
  // safety check
  CheckInit();

  // --------------------------------------
  // control parameters
  // --------------------------------------
  // set target time to initial time plus step size
  timenp_ = (*timen_)[0] + (*dt_)[0];
  stepnp_ = stepn_ + 1;

  timer_ = Teuchos::rcp(new Epetra_Time(*comm_));
  // --------------------------------------
  // vectors
  // --------------------------------------
  // displacements D_{n}
  dis_ = Teuchos::rcp(new ::TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  // velocities V_{n}
  vel_ = Teuchos::rcp(new ::TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  // accelerations A_{n}
  acc_ = Teuchos::rcp(new ::TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

  // displacements D_{n+1} at t_{n+1}
  disnp_ = LINALG::CreateVector(*DofRowMapView(), true);
  // velocities V_{n+1} at t_{n+1}
  velnp_ = LINALG::CreateVector(*DofRowMapView(), true);
  // accelerations A_{n+1} at t_{n+1}
  accnp_ = LINALG::CreateVector(*DofRowMapView(), true);

  fintn_ = LINALG::CreateVector(*DofRowMapView(), true);
  fintnp_ = LINALG::CreateVector(*DofRowMapView(), true);

  fextn_ = LINALG::CreateVector(*DofRowMapView(), true);
  fextnp_ = LINALG::CreateVector(*DofRowMapView(), true);

  freactnp_ = LINALG::CreateVector(*DofRowMapView(), true);

  finertialn_ = LINALG::CreateVector(*DofRowMapView(), true);
  finertialnp_ = LINALG::CreateVector(*DofRowMapView(), true);

  fviscon_ = LINALG::CreateVector(*DofRowMapView(), true);
  fvisconp_ = LINALG::CreateVector(*DofRowMapView(), true);

  fstructold_ = LINALG::CreateVector(*DofRowMapView(), true);

  // --------------------------------------
  // sparse operators
  // --------------------------------------
  mass_  = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));
  if (datasdyn_->GetDampingType() != INPAR::STR::damp_none)
  {
    if (datasdyn_->GetMassLinType() == INPAR::STR::ml_none)
    {
      damp_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));
    }
    else
    {
      /* Since our element evaluate routine is only designed for two input matrices
       * (stiffness and damping or stiffness and mass) its not possible, to have nonlinear
       * inertia forces AND material damping. */
      dserror("So far it is not possible to model nonlinear inertia forces and damping!");
    }
  }

  issetup_ = true;

  // Good bye
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Vector> STR::TIMINT::BaseDataGlobalState::
    CreateGlobalVector() const
{
  return CreateGlobalVector(vec_init_zero,Teuchos::null);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::BaseDataGlobalState::SetupBlockInformation(
    const STR::MODELEVALUATOR::Generic& me,
    const INPAR::STR::ModelType& mt)
{
  CheckInit();
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<const Epetra_Map> me_map_ptr = me.GetBlockDofRowMapPtr();

  model_maps_[mt] = me_map_ptr;

  switch (mt)
  {
    case INPAR::STR::model_structure:
    {
      // always called first, so we can use it to reset things
      gproblem_map_ptr_ = Teuchos::null;
      model_block_id_[mt] = 0;
      max_block_num_    = 1;
      break;
    }
    case INPAR::STR::model_contact:
    {
      enum INPAR::CONTACT::SystemType systype =
          DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(
              problem->ContactDynamicParams(),"SYSTEM");

      // --- saddle-point system
      if (systype == INPAR::CONTACT::system_saddlepoint)
      {
        model_block_id_[mt] = max_block_num_;
        ++max_block_num_;
      }
      // --- condensed system
      else
        model_block_id_[mt] = 0;
      break;
    }
    case INPAR::STR::model_cardiovascular0d:
    {
      // --- 2x2 block system
      model_block_id_[mt] = max_block_num_;
      ++max_block_num_;
      break;
    }
    case INPAR::STR::model_lag_pen_constraint:
    {
      // ---------------------------------------------------------------------------
      // check type of constraint conditions (Lagrange multiplier vs. penalty)
      // ---------------------------------------------------------------------------
      bool have_lag_constraint = false;
      std::vector<DRT::Condition*> lagcond_volconstr3d(0);
      std::vector<DRT::Condition*> lagcond_areaconstr3d(0);
      std::vector<DRT::Condition*> lagcond_areaconstr2d(0);
      std::vector<DRT::Condition*> lagcond_mpconline2d(0);
      std::vector<DRT::Condition*> lagcond_mpconplane3d(0);
      std::vector<DRT::Condition*> lagcond_mpcnormcomp3d(0);
      discret_->GetCondition("VolumeConstraint_3D",lagcond_volconstr3d);
      discret_->GetCondition("AreaConstraint_3D",lagcond_areaconstr3d);
      discret_->GetCondition("AreaConstraint_2D",lagcond_areaconstr2d);
      discret_->GetCondition("MPC_NodeOnLine_2D",lagcond_mpconline2d);
      discret_->GetCondition("MPC_NodeOnPlane_3D",lagcond_mpconplane3d);
      discret_->GetCondition("MPC_NormalComponent_3D",lagcond_mpcnormcomp3d);
      if (
             lagcond_volconstr3d.size()  or
             lagcond_areaconstr3d.size() or
             lagcond_areaconstr2d.size() or
             lagcond_mpconline2d.size()  or
             lagcond_mpconplane3d.size() or
             lagcond_mpcnormcomp3d.size()
          )
        have_lag_constraint = true;

      // --- 2x2 block system (saddle-point structure)
      if (have_lag_constraint)
      {
        model_block_id_[mt] = max_block_num_;
        ++max_block_num_;
      }
      // --- standard system
      else
      {
        model_block_id_[mt] = 0;
      }
      break;
    }
    case INPAR::STR::model_springdashpot:
    case INPAR::STR::model_beam_interaction_old:
    case INPAR::STR::model_browniandyn:
    case INPAR::STR::model_beaminteraction:
    {
      // structural block
      model_block_id_[mt] = 0;
      break;
    }
    case INPAR::STR::model_partitioned_coupling:
    {
      // do nothing
      break;
    }
    default:
    {
      // FixMe please
      dserror("Augment this function for your model type!");
      break;
    }
  }
  // create a global problem map
  gproblem_map_ptr_ = LINALG::MergeMap(gproblem_map_ptr_,me_map_ptr);

  return gproblem_map_ptr_->MaxAllGID();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::SetupMultiMapExtractor()
{
  CheckInit();
  /* copy the std::map into a std::vector and keep the numbering of the model-id
   * map */
  std::vector<Teuchos::RCP<const Epetra_Map> > maps_vec(MaxBlockNumber(),
      Teuchos::null);
  // Make sure, that the block ids and the vector entry ids coincide!
  std::map<INPAR::STR::ModelType,int>::const_iterator ci;
  for (ci=model_block_id_.begin();ci!=model_block_id_.end();++ci)
  {
    enum INPAR::STR::ModelType mt = ci->first;
    int bid = ci->second;
    maps_vec[bid] = model_maps_.at(mt);
  }
  blockextractor_.Setup(*gproblem_map_ptr_,maps_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::SetupRotVecMapExtractor()
{
  CheckInit();

  /* all additive DoFs, i.e. members of real value vector spaces
   * such as translational displacements, tangent vector displacements,
   * 1D rotation angles, ... */
  std::set<int> additdofset;
  /* DoFs which are non-additive and therefore e.g. can not be updated in usual
   * incremental manner, need special treatment in time integration ...
   * (currently only rotation pseudo-vector DoFs of beam elements) */
  std::set<int> rotvecdofset;

  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* nodeptr = discret_->lRowNode(i);

    const DRT::ELEMENTS::Beam3Base* beameleptr =
        dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(nodeptr->Elements()[0]);

    std::vector<int> nodaladditdofs;
    std::vector<int> nodalrotvecdofs;

    // so far we only expect DoFs of beam elements for the rotvecdofset
    if (beameleptr == NULL)
    {
      nodaladditdofs = discret_->Dof(0,nodeptr);
    }
    else
    {
      nodaladditdofs = beameleptr->GetAdditiveDofGIDs(*discret_,*nodeptr);
      nodalrotvecdofs = beameleptr->GetRotVecDofGIDs(*discret_,*nodeptr);

      if (nodaladditdofs.size() + nodalrotvecdofs.size() !=
          (unsigned) beameleptr->NumDofPerNode(*nodeptr))
        dserror("Expected %d DoFs for node with GID %d but collected %d DoFs",
            beameleptr->NumDofPerNode(*nodeptr),discret_->NodeRowMap()->GID(i),
            nodaladditdofs.size() + nodalrotvecdofs.size());
    }

    // add the DoFs of this node to the total set
    for (unsigned j=0; j<nodaladditdofs.size(); ++j)
      additdofset.insert(nodaladditdofs[j]);

    for (unsigned j=0; j<nodalrotvecdofs.size(); ++j)
      rotvecdofset.insert(nodalrotvecdofs[j]);

  } // loop over row nodes

  // create the required Epetra maps
  std::vector<int> additdofmapvec;
  additdofmapvec.reserve(additdofset.size());
  additdofmapvec.assign(additdofset.begin(), additdofset.end());
  additdofset.clear();
  Teuchos::RCP<Epetra_Map> additdofmap =
    Teuchos::rcp(new Epetra_Map(-1,additdofmapvec.size(),
                                &additdofmapvec[0],
                                0,
                                discret_->Comm()));
  additdofmapvec.clear();

  std::vector<int> rotvecdofmapvec;
  rotvecdofmapvec.reserve(rotvecdofset.size());
  rotvecdofmapvec.assign(rotvecdofset.begin(), rotvecdofset.end());
  rotvecdofset.clear();
  Teuchos::RCP<Epetra_Map> rotvecdofmap =
    Teuchos::rcp(new Epetra_Map(-1,rotvecdofmapvec.size(),
                                &rotvecdofmapvec[0],
                                0,
                                discret_->Comm()));
  rotvecdofmapvec.clear();

  std::vector<Teuchos::RCP<const Epetra_Map> > maps( 2 );
  maps[0] = additdofmap;
  maps[1] = rotvecdofmap;

  rotvecextractor_.Setup(*DofRowMapView(),maps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const LINALG::MultiMapExtractor& STR::TIMINT::BaseDataGlobalState::
    BlockExtractor() const
{
  // sanity check
  blockextractor_.CheckForValidMapExtractor();
  return blockextractor_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const LINALG::MultiMapExtractor& STR::TIMINT::BaseDataGlobalState::
    RotVecExtractor() const
{
  // sanity check
  rotvecextractor_.CheckForValidMapExtractor();
  return rotvecextractor_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Vector> STR::TIMINT::BaseDataGlobalState::
    CreateGlobalVector(const enum VecInitType& vecinittype,
    const Teuchos::RCP<const STR::ModelEvaluator>& modeleval_ptr) const
{
  CheckInit();
  Teuchos::RCP<Epetra_Vector> xvec_ptr =
      Teuchos::rcp(new Epetra_Vector(GlobalProblemMap(),true));

  // switch between the different vector initialization options
  switch (vecinittype)
  {
    /* use the last converged state to construct a new solution vector */
    case vec_init_last_time_step:
    {
      if (modeleval_ptr.is_null())
        dserror("We need access to the STR::ModelEvaluator object!");

      std::map<INPAR::STR::ModelType,int>::const_iterator ci;
      for (ci=model_block_id_.begin();ci!=model_block_id_.end();++ci)
      {
        // get the partial solution vector of the last time step
        Teuchos::RCP<const Epetra_Vector> model_sol_ptr =
            modeleval_ptr->Evaluator(ci->first).GetLastTimeStepSolutionPtr();
        // if there is a partial solution, we insert it into the full vector
        if (not model_sol_ptr.is_null())
          BlockExtractor().InsertVector(model_sol_ptr,ci->second,xvec_ptr);
        model_sol_ptr = Teuchos::null;
      }
      break;
    }
    /* use the current global state to construct a new solution vector */
    case vec_init_current_state:
    {
      if (modeleval_ptr.is_null())
        dserror("We need access to the STR::ModelEvaluator object!");

      std::map<INPAR::STR::ModelType,int>::const_iterator ci;
      for (ci=model_block_id_.begin();ci!=model_block_id_.end();++ci)
      {
        // get the partial solution vector of the current state
        Teuchos::RCP<const Epetra_Vector> model_sol_ptr =
            modeleval_ptr->Evaluator(ci->first).GetCurrentSolutionPtr();
        // if there is a partial solution, we insert it into the full vector
        if (not model_sol_ptr.is_null())
          BlockExtractor().InsertVector(model_sol_ptr,ci->second,xvec_ptr);
      }
      break;
    }
    /* construct a new solution vector filled with zeros */
    case vec_init_zero:
    default:
    {
      // nothing to do.
      break;
    }
  } // end of the switch-case statement

  //wrap and return
  return Teuchos::rcp(new NOX::Epetra::Vector(xvec_ptr));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> STR::TIMINT::BaseDataGlobalState::
    CreateJacobian()
{
  CheckInit();
  jac_ = Teuchos::null;

  if (max_block_num_>1)
  {
    jac_ = Teuchos::rcp(new LINALG::BlockSparseMatrix
         <LINALG::DefaultBlockMatrixStrategy>(BlockExtractor(),BlockExtractor(),
             81,true,true));
  }
  else
  {
    // pure structural case
    jac_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(),
        81, true, true));
  }

  return jac_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::TIMINT::BaseDataGlobalState::DofRowMap()
    const
{
  CheckInit();
  const Epetra_Map* dofrowmap_ptr = discret_->DofRowMap();
  // since it's const, we do not need to copy the map
  return Teuchos::rcp(dofrowmap_ptr,false);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::TIMINT::BaseDataGlobalState::DofRowMap(
    unsigned nds) const
{
  CheckInit();
  const Epetra_Map* dofrowmap_ptr = discret_->DofRowMap(nds);
  // since it's const, we do not need to copy the map
  return Teuchos::rcp(dofrowmap_ptr,false);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* STR::TIMINT::BaseDataGlobalState::DofRowMapView() const
{
  CheckInit();
  return discret_->DofRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* STR::TIMINT::BaseDataGlobalState::AdditiveDofRowMapView() const
{
  CheckInit();
  return RotVecExtractor().Map(0).get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* STR::TIMINT::BaseDataGlobalState::RotVecDofRowMapView() const
{
  CheckInit();
  return RotVecExtractor().Map(1).get();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::BaseDataGlobalState::ExtractDisplEntries(
    const Epetra_Vector& source) const
{
  return ExtractModelEntries(INPAR::STR::model_structure,source);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::BaseDataGlobalState::ExtractModelEntries(
    const INPAR::STR::ModelType& mt,
    const Epetra_Vector& source) const
{
  Teuchos::RCP<Epetra_Vector> model_ptr = Teuchos::null;
  // extract from the full state vector
  if (source.Map().NumGlobalElements() ==
      BlockExtractor().FullMap()->NumGlobalElements())
  {
    model_ptr =
        BlockExtractor().ExtractVector(source,model_block_id_.at(mt));
  }
  // copy the vector
  else if (source.Map().NumGlobalElements() ==
      model_maps_.at(mt)->NumGlobalElements())
  {
    model_ptr = Teuchos::rcp(new Epetra_Vector(source));
  }
  // otherwise do a standard export
  else
  {
    model_ptr = Teuchos::rcp(new Epetra_Vector(*model_maps_.at(mt)));
    LINALG::Export(source,*model_ptr);
  }


  return model_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::BaseDataGlobalState::ExtractAdditiveEntries(
    const Epetra_Vector& source) const
{
  Teuchos::RCP<Epetra_Vector> addit_ptr = RotVecExtractor().ExtractVector(source,0);

  return addit_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::BaseDataGlobalState::ExtractRotVecEntries(
    const Epetra_Vector& source) const
{
  Teuchos::RCP<Epetra_Vector> addit_ptr = RotVecExtractor().ExtractVector(source,1);

  return addit_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataGlobalState::AssignModelBlock(
    LINALG::SparseOperator& jac,
    const LINALG::SparseMatrix& matrix,
    const INPAR::STR::ModelType& mt,
    const MatBlockType& bt,
    const LINALG::DataAccess& access) const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* blockmat_ptr =
      dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* >(&jac);
  if (blockmat_ptr != NULL)
  {
    if (MaxBlockNumber()<2)
      dserror("The jacobian is a LINALG::BlockSparseMatrix but has less than"
          " two blocks! Seems wrong.");

    const int& b_id = model_block_id_.at(mt);
    switch (bt)
    {
      case block_displ_displ:
      {
        blockmat_ptr->Matrix(0,0).Assign(access,matrix);
        break;
      }
      case block_displ_lm:
      {
        blockmat_ptr->Matrix(0,b_id).Assign(access,matrix);
        break;
      }
      case block_lm_displ:
      {
        blockmat_ptr->Matrix(b_id,0).Assign(access,matrix);
        break;
      }
      case block_lm_lm:
      {
        blockmat_ptr->Matrix(b_id,b_id).Assign(access,matrix);
        break;
      }
    }
    return;
  }

  LINALG::SparseMatrix* stiff_ptr = dynamic_cast<LINALG::SparseMatrix*>(&jac);
  if (stiff_ptr!=NULL)
  {
    stiff_ptr->Assign(access,matrix);
    return;
  }

  dserror("The jacobian has the wrong type! (no LINALG::SparseMatrix "
      "and no LINALG::BlockSparseMatrix)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::BaseDataGlobalState::
    ExtractModelBlock(
    LINALG::SparseOperator& jac,
    const INPAR::STR::ModelType& mt,
    const MatBlockType& bt) const
{
  Teuchos::RCP<LINALG::SparseMatrix> block = Teuchos::null;
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* blockmat_ptr =
          dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* >(&jac);
  if (blockmat_ptr!=NULL)
  {
    if (MaxBlockNumber()<2)
      dserror("The jacobian is a LINALG::BlockSparseMatrix but has less than"
          " two blocks! Seems wrong.");
    const int& b_id = model_block_id_.at(mt);
    switch (bt)
    {
      case block_displ_displ:
      {
        block = Teuchos::rcp(&(blockmat_ptr->Matrix(0,0)),false);
        break;
      }
      case block_displ_lm:
      {
        block = Teuchos::rcp(&(blockmat_ptr->Matrix(0,b_id)),false);
        break;
      }
      case block_lm_displ:
      {
        block = Teuchos::rcp(&(blockmat_ptr->Matrix(b_id,0)),false);
        break;
      }
      case block_lm_lm:
      {
        block = Teuchos::rcp(&(blockmat_ptr->Matrix(b_id,b_id)),false);
        break;
      }
    }
    return block;
  }

  LINALG::SparseMatrix* stiff_ptr = dynamic_cast<LINALG::SparseMatrix*>(&jac);
  if (stiff_ptr!=NULL)
  {
    block = Teuchos::rcp(stiff_ptr,false);
    return block;
  }

  dserror("The jacobian has the wrong type! (no LINALG::SparseMatrix "
      "and no LINALG::BlockSparseMatrix)");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::BaseDataGlobalState::
    ExtractDisplBlock(LINALG::SparseOperator& jac) const
{
  return ExtractModelBlock(jac,INPAR::STR::model_structure,block_displ_displ);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> STR::TIMINT::BaseDataGlobalState::
    GetJacobianDisplBlock() const
{
  if (jac_.is_null())
    dserror("The jacobian is not initialized!");
  return ExtractDisplBlock(*jac_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::BaseDataGlobalState::
    JacobianDisplBlock()
{
  if (jac_.is_null())
    dserror("The jacobian is not initialized!");
  return ExtractDisplBlock(*jac_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOp::TIMINT::RotVecUpdater::RotVecUpdater(
    const Teuchos::RCP<const ::STR::TIMINT::BaseDataGlobalState>& gstate_ptr)
    : gstate_ptr_(gstate_ptr)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GROUP::PrePostOp::TIMINT::RotVecUpdater::runPreComputeX(
    const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir,
    const double& step,
    const NOX::NLN::Group& curr_grp)
{
  const Epetra_Vector& xold =
      dynamic_cast<const NOX::Epetra::Vector&>(input_grp.getX()).getEpetraVector();

  // cast the const away so that the new x vector can be set after the update
  NOX::NLN::Group& curr_grp_mutable = const_cast<NOX::NLN::Group&>(curr_grp);

  Teuchos::RCP<Epetra_Vector> xnew = Teuchos::rcp(new Epetra_Vector(xold.Map(),true));

  /* we do the multiplicative update only for those entries which belong to
   * rotation (pseudo-)vectors */
  Epetra_Vector x_rotvec = *gstate_ptr_->ExtractRotVecEntries(xold);
  Epetra_Vector dir_rotvec = *gstate_ptr_->ExtractRotVecEntries(dir);

  LINALG::Matrix<4,1> Qold;
  LINALG::Matrix<4,1> deltaQ;
  LINALG::Matrix<4,1> Qnew;

  /* since parallel distribution is node-wise, the three entries belonging to
   * a rotation vector should be stored on the same processor: safety-check */
  if (x_rotvec.Map().NumMyElements() %3 !=0 or dir_rotvec.Map().NumMyElements() %3 !=0)
    dserror("fatal error: apparently, the three DOFs of a nodal rotation vector are"
        " not stored on this processor. Can't apply multiplicative update!");

  // rotation vectors always consist of three consecutive DoFs
  for (int i=0; i<x_rotvec.Map().NumMyElements(); i=i+3)
  {
    // create a LINALG::Matrix from reference to three x vector entries
    LINALG::Matrix<3,1> theta(&x_rotvec[i],true);
    LARGEROTATIONS::angletoquaternion(theta,Qold);

    // same for relative rotation angle deltatheta
    LINALG::Matrix<3,1> deltatheta(&dir_rotvec[i],true);
    deltatheta.Scale(step);

    LARGEROTATIONS::angletoquaternion(deltatheta,deltaQ);
    LARGEROTATIONS::quaternionproduct(Qold,deltaQ,Qnew);
    LARGEROTATIONS::quaterniontoangle(Qnew,theta);
  }

  // first update entire x vector in an additive manner
  xnew->Update(1.0, xold, step, dir,0.0);

  // now replace the rotvec entries by the correct value computed before
  STR::AssembleVector(0.0,*xnew,1.0,x_rotvec);
  curr_grp_mutable.setX(xnew);

  /* tell the NOX::NLN::Group that the x vector has already been updated in
   * this preComputeX operator call */
  curr_grp_mutable.setSkipUpdateX(true);

  return;
}
