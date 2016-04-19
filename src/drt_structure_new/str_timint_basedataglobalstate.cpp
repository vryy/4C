/*-----------------------------------------------------------*/
/*!
\file str_timint_basedataglobalstate.cpp

\maintainer Michael Hiermeier

\date Jan 12, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_timint_basedataglobalstate.H"
#include "str_timint_basedatasdyn.H"

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
      dtele_(0.0)
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
  jac_ = CreateJacobian();
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
    CreateGlobalVector()
{
  return CreateGlobalVector(vec_init_zero);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Vector> STR::TIMINT::BaseDataGlobalState::
    CreateGlobalVector(const enum VecInitType& vecinittype)
{
  CheckInit();
  Teuchos::RCP<Epetra_Vector> xvec = Teuchos::null;
  DRT::Problem* problem = DRT::Problem::Instance();

  // get modeltypes
  const std::set<enum INPAR::STR::ModelType>& modeltypes =
      datasdyn_->GetModelTypes();

  // number of models
  std::size_t nummodels = modeltypes.size();
  int numblocks = 0;
  // ---------------------------------------------------------------------------
  // if there are more than one active model type
  // ---------------------------------------------------------------------------
  if(nummodels>1)
  {
    // count blocks
    std::set<enum INPAR::STR::ModelType>::const_iterator miter;
    for (miter=modeltypes.begin();miter!=modeltypes.end();++miter)
    {
      switch (*miter)
      {
      case INPAR::STR::model_structure:
      {
        ++numblocks;
        break;
      }
      case INPAR::STR::model_contact:
      {
        enum INPAR::CONTACT::SystemType systype =
            DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(
                problem->ContactDynamicParams(),"SYSTEM");

        if (systype == INPAR::CONTACT::system_saddlepoint)
          ++numblocks;
        break;
      }
      case INPAR::STR::model_springdashpot:
        break;
      default:
      {
        // ToDo
        dserror("Augment this function for your model type!");
        break;
      }
      }
    }
    // ToDo
    if(numblocks>1)
      dserror("The corresponding mnodel evaluators are necessary to create "
          "the global vector.");
  } // end of the case of more than one model type

  // ---------------------------------------------------------------------------
  // pure structural case
  // ---------------------------------------------------------------------------
  {
    // switch between the different vector initialization options
    switch (vecinittype)
    {
    /* use the last converged state to construct a new solution vector */
    case vec_init_last_time_step:
    {
      xvec = Teuchos::rcp(new Epetra_Vector(*DofRowMapView(),false));
      if (xvec->Scale(1.0,*GetDisN()))
        dserror("Scale operation failed!");
      break;
    }
    /* use the current global state to construct a new solution vector */
    case vec_init_current_state:
    {
      xvec = Teuchos::rcp(new Epetra_Vector(*DofRowMapView(),false));
      if (xvec->Scale(1.0,*GetDisNp()))
        dserror("Scale operation failed!");
      break;
    }
    /* construct a new solution vector filled with zeros */
    case vec_init_zero:
    default:
    {
      xvec = Teuchos::rcp(new Epetra_Vector(*DofRowMapView(),true));
      break;
    }
    } // end of the switch-case statement
  } // end of the pure structural problem case

  //wrap and return
  return Teuchos::rcp(new NOX::Epetra::Vector(xvec));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> STR::TIMINT::BaseDataGlobalState::
    CreateJacobian()
{
  CheckInit();
  Teuchos::RCP<LINALG::SparseOperator> jac = Teuchos::null;
  DRT::Problem* problem = DRT::Problem::Instance();

  // get modeltypes
  const std::set<enum INPAR::STR::ModelType>& modeltypes =
      datasdyn_->GetModelTypes();

  // number of models
  std::size_t nummodels = modeltypes.size();
  int numblocks = 0;
  if(nummodels>1)
  {
    std::set<enum INPAR::STR::ModelType>::const_iterator miter;
    for (miter=modeltypes.begin();miter!=modeltypes.end();++miter)
    {
      switch (*miter)
      {
        case INPAR::STR::model_structure:
        {
          ++numblocks;
          break;
        }
        case INPAR::STR::model_contact:
        {
          enum INPAR::CONTACT::SystemType systype =
              DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(
              problem->ContactDynamicParams(),"SYSTEM");

          if (systype == INPAR::CONTACT::system_saddlepoint)
            ++numblocks;
          break;
        }
        case INPAR::STR::model_springdashpot:
          break;
        default:
        {
          // FixMe please.
          dserror("Augment this function for your model type!");
          break;
        }
      }
    }
    // FixMe please.
    if(numblocks>1)
      dserror("The corresponding managers/model-evaluators are necessary to create "
          "the BlockMatrix.");
  }
  // pure structural case
  jac =
      Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMapView(), 81, true, true));
  return jac;
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
Teuchos::RCP<LINALG::SparseMatrix> STR::TIMINT::BaseDataGlobalState::
    ExtractDisplBlock(LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<LINALG::SparseMatrix> stiff = Teuchos::null;

  const LINALG::SparseOperator* testop = 0;
  testop = dynamic_cast<const LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* >(&jac);
  if (testop != NULL)
  {
    LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* blockmat =
        dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* >(&jac);
    // get block (0,0)
    stiff = Teuchos::rcp(&(blockmat->Matrix(0,0)),false);
  }

  testop = dynamic_cast<const LINALG::SparseMatrix*>(&jac);
  if (testop == NULL)
    dserror("Structural stiffness matrix has to be a LINALG::BlockSparseMatrix or "
        "a LINALG::SparseMatrix!");

  stiff = Teuchos::rcp(dynamic_cast<LINALG::SparseMatrix*>(&jac),false);

  return stiff;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> STR::TIMINT::BaseDataGlobalState::
    GetJacobianDisplBlock() const
{
  return ExtractDisplBlock(*jac_);
}
