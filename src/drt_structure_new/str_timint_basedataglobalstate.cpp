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

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"

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
      fstructn_(Teuchos::null),
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
  disnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  // velocities V_{n+1} at t_{n+1}
  velnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  // accelerations A_{n+1} at t_{n+1}
  accnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  fintn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  fintnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  fextn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  fextnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  freactnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  finertialn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  finertialnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  fviscon_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
  fvisconp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  fstructn_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

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
  const int dof_off_set = me_map_ptr->MaxAllGID();

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
      // --- 2x2 block system (saddle-point structure)
      model_block_id_[mt] = max_block_num_;
      ++max_block_num_;
      break;
    }
    case INPAR::STR::model_springdashpot:
    {
      // structural block
      model_block_id_[mt] = 0;
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

  return dof_off_set;
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
const LINALG::MultiMapExtractor& STR::TIMINT::BaseDataGlobalState::
    BlockExtractor() const
{
  // sanity check
  blockextractor_.CheckForValidMapExtractor();
  return blockextractor_;
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
Teuchos::RCP<Epetra_Vector> STR::TIMINT::BaseDataGlobalState::ExportDisplEntries(
    const Epetra_Vector& source) const
{
  Teuchos::RCP<Epetra_Vector> displ_ptr = Teuchos::null;
  if (not source.Map().PointSameAs(*DofRowMapView()))
  {
      displ_ptr = Teuchos::rcp(new Epetra_Vector(*DofRowMapView()));
      LINALG::Export(source,*displ_ptr);
  }
  else
    displ_ptr = Teuchos::rcp(new Epetra_Vector(source));

  return displ_ptr;
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
  if (MaxBlockNumber()>1)
  {
    LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* blockmat =
        dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* >(&jac);
    if (blockmat == NULL)
      dserror("The jacobian has the wrong type! (no LINALG::BlockSparseMatrix)");

    const int& b_id = model_block_id_.at(mt);
    switch (bt)
    {
      case block_displ_displ:
      {
        blockmat->Matrix(0,0).Assign(access,matrix);
        break;
      }
      case block_displ_lm:
      {
        blockmat->Matrix(0,b_id).Assign(access,matrix);
        break;
      }
      case block_lm_displ:
      {
        blockmat->Matrix(b_id,0).Assign(access,matrix);
        break;
      }
      case block_lm_lm:
      {
        blockmat->Matrix(b_id,b_id).Assign(access,matrix);
        break;
      }
    }
  }
  else
  {
    LINALG::SparseMatrix* stiff_ptr = dynamic_cast<LINALG::SparseMatrix*>(&jac);
    if (stiff_ptr == NULL)
      dserror("The jacobian has the wrong type! (no LINALG::SparseMatrix)");
    stiff_ptr->Assign(access,matrix);
  }
  return;
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
  if (MaxBlockNumber()>1)
  {
    LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* blockmat =
        dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* >(&jac);

    if (blockmat == NULL)
      dserror("The jacobian has the wrong type! (no LINALG::BlockSparseMatrix)");
    const int& b_id = model_block_id_.at(mt);
    switch (bt)
    {
      case block_displ_displ:
      {
        block = Teuchos::rcp(&(blockmat->Matrix(0,0)),false);
        break;
      }
      case block_displ_lm:
      {
        block = Teuchos::rcp(&(blockmat->Matrix(0,b_id)),false);
        break;
      }
      case block_lm_displ:
      {
        block = Teuchos::rcp(&(blockmat->Matrix(b_id,0)),false);
        break;
      }
      case block_lm_lm:
      {
        block = Teuchos::rcp(&(blockmat->Matrix(b_id,b_id)),false);
        break;
      }
    }
  }
  else
  {
    LINALG::SparseMatrix* stiff_ptr = dynamic_cast<LINALG::SparseMatrix*>(&jac);
    if (stiff_ptr == NULL)
      dserror("The jacobian has the wrong type! (no LINALG::SparseMatrix)");
    block = Teuchos::rcp(stiff_ptr,false);
  }
  return block;
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
  return ExtractDisplBlock(*jac_);
}
