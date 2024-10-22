// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_geometry_periodic_boundingbox.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"

#include <Epetra_FEVector.h>

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
  dis_ = Teuchos::make_rcp<TimeStepping::TimIntMStep<Core::LinAlg::Vector<double>>>(
      0, 0, ia_discret->dof_row_map(), true);
  disnp_ = Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*ia_discret->dof_col_map());
  discolnp_ = Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*ia_discret->dof_col_map());

  // force
  forcen_ = Teuchos::make_rcp<Epetra_FEVector>(*ia_discret->dof_row_map());
  forcenp_ = Teuchos::make_rcp<Epetra_FEVector>(*ia_discret->dof_row_map());

  stiff_ = Teuchos::make_rcp<Core::LinAlg::SparseMatrix>(
      *ia_discret->dof_row_map(), 81, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
