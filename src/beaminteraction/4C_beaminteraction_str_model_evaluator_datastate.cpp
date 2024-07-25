/*-----------------------------------------------------------*/
/*! \file

\brief Global state data container for the structural (time)
       integration


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_geometry_periodic_boundingbox.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::ModelEvaluator::BeamInteractionDataState::BeamInteractionDataState()
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
void Solid::ModelEvaluator::BeamInteractionDataState::init()
{
  // We have to call setup() after init()
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
void Solid::ModelEvaluator::BeamInteractionDataState::setup(
    Teuchos::RCP<const Core::FE::Discretization> const& ia_discret)
{
  // safety check
  check_init();

  myrank_ = ia_discret->get_comm().MyPID();

  // displacements
  dis_ = Teuchos::rcp(
      new TimeStepping::TimIntMStep<Epetra_Vector>(0, 0, ia_discret->dof_row_map(), true));
  disnp_ = Teuchos::rcp(new Epetra_Vector(*ia_discret->dof_col_map()));
  discolnp_ = Teuchos::rcp(new Epetra_Vector(*ia_discret->dof_col_map()));

  // force
  forcen_ = Teuchos::rcp(new Epetra_FEVector(*ia_discret->dof_row_map()));
  forcenp_ = Teuchos::rcp(new Epetra_FEVector(*ia_discret->dof_row_map()));

  stiff_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *ia_discret->dof_row_map(), 81, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX));

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
