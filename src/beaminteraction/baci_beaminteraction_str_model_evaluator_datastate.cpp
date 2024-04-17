/*-----------------------------------------------------------*/
/*! \file

\brief Global state data container for the structural (time)
       integration


\level 3

*/
/*-----------------------------------------------------------*/


#include "baci_beaminteraction_str_model_evaluator_datastate.hpp"

#include "baci_beaminteraction_periodic_boundingbox.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteractionDataState::BeamInteractionDataState()
    : isinit_(false),
      issetup_(false),
      myrank_(0),
      dis_(Teuchos::null),
      dis_restart_(Teuchos::null),
      dis_restart_col_(Teuchos::null),
      is_restart_coupling_(false),
      disnp_(Teuchos::null),
      discolnp_(Teuchos::null),
      forcen_(Teuchos::null),
      forcenp_(Teuchos::null),
      stiff_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteractionDataState::Init()
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // clear stl stuff
  bintorowelemap_.clear();
  exbintorowelemap_.clear();
  roweletobinmap_.clear();

  // end of initialization
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::BeamInteractionDataState::Setup(
    Teuchos::RCP<const DRT::Discretization> const& ia_discret)
{
  // safety check
  CheckInit();

  myrank_ = ia_discret->Comm().MyPID();

  // displacements
  dis_ = Teuchos::rcp(
      new TIMESTEPPING::TimIntMStep<Epetra_Vector>(0, 0, ia_discret->DofRowMap(), true));
  disnp_ = Teuchos::rcp(new Epetra_Vector(*ia_discret->DofColMap()));
  discolnp_ = Teuchos::rcp(new Epetra_Vector(*ia_discret->DofColMap()));

  // force
  forcen_ = Teuchos::rcp(new Epetra_FEVector(*ia_discret->DofRowMap()));
  forcenp_ = Teuchos::rcp(new Epetra_FEVector(*ia_discret->DofRowMap()));

  stiff_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *ia_discret->DofRowMap(), 81, true, true, CORE::LINALG::SparseMatrix::FE_MATRIX));

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
