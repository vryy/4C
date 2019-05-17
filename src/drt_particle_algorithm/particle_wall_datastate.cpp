/*---------------------------------------------------------------------------*/
/*!
\file particle_wall_datastate.cpp

\brief wall data state container for particle wall handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_wall_datastate.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::WallDataState::WallDataState()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init wall data state container                             sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallDataState::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup wall data state container                            sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallDataState::Setup(
    Teuchos::RCP<const DRT::Discretization> const& walldiscretization, bool ismoving, bool isloaded)
{
  // set current dof row and column map
  curr_dof_row_map_ = Teuchos::rcp(new Epetra_Map(*walldiscretization->DofRowMap()));

  // create states needed for moving walls
  if (ismoving)
  {
    disp_row_ = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    disp_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap()));
    disp_row_last_transfer_ = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    vel_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap()));
    acc_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap()));
  }

  // create states needed for loaded walls
  if (isloaded)
  {
    force_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap()));
  }
}

/*---------------------------------------------------------------------------*
 | check for correct maps                                     sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallDataState::CheckForCorrectMaps(
    Teuchos::RCP<const DRT::Discretization> const& walldiscretization)
{
  if (disp_row_ != Teuchos::null)
    if (not disp_row_->Map().SameAs(*walldiscretization->DofRowMap()))
      dserror("map of state 'disp_row_' corrupt!");

  if (disp_col_ != Teuchos::null)
    if (not disp_col_->Map().SameAs(*walldiscretization->DofColMap()))
      dserror("map of state 'disp_col_' corrupt!");

  if (disp_row_last_transfer_ != Teuchos::null)
    if (not disp_row_last_transfer_->Map().SameAs(*walldiscretization->DofRowMap()))
      dserror("map of state 'disp_row_last_transfer_' corrupt!");

  if (vel_col_ != Teuchos::null)
    if (not vel_col_->Map().SameAs(*walldiscretization->DofColMap()))
      dserror("map of state 'vel_col_' corrupt!");

  if (acc_col_ != Teuchos::null)
    if (not acc_col_->Map().SameAs(*walldiscretization->DofColMap()))
      dserror("map of state 'acc_col_' corrupt!");

  if (force_col_ != Teuchos::null)
    if (not force_col_->Map().SameAs(*walldiscretization->DofColMap()))
      dserror("map of state 'force_col_' corrupt!");
}

/*---------------------------------------------------------------------------*
 | update maps of state vectors                               sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::WallDataState::UpdateMapsOfStateVectors(
    Teuchos::RCP<const DRT::Discretization> const& walldiscretization)
{
  if (disp_row_ != Teuchos::null and disp_col_ != Teuchos::null)
  {
    // export row map based displacement vector
    Teuchos::RCP<Epetra_Vector> temp = disp_row_;
    disp_row_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofRowMap(), true));
    LINALG::Export(*temp, *disp_row_);

    // update column map based displacement vector
    disp_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap(), true));
    LINALG::Export(*disp_row_, *disp_col_);

    // store displacements after last transfer
    disp_row_last_transfer_ = Teuchos::rcp(new Epetra_Vector(*disp_row_));
  }

  if (vel_col_ != Teuchos::null)
  {
    // export old column to old row map based vector (no communication)
    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    LINALG::Export(*vel_col_, *temp);
    // export old row map based vector to new column map based vector
    vel_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap(), true));
    LINALG::Export(*temp, *vel_col_);
  }

  if (acc_col_ != Teuchos::null)
  {
    // export old column to old row map based vector (no communication)
    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    LINALG::Export(*acc_col_, *temp);
    // export old row map based vector to new column map based vector
    acc_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap(), true));
    LINALG::Export(*temp, *acc_col_);
  }

  if (force_col_ != Teuchos::null)
  {
    // export old column to old row map based vector (no communication)
    Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*curr_dof_row_map_));
    LINALG::Export(*force_col_, *temp);
    // export old row map based vector to new column map based vector
    force_col_ = Teuchos::rcp(new Epetra_Vector(*walldiscretization->DofColMap(), true));
    LINALG::Export(*temp, *force_col_);
  }

  // set new dof row map
  curr_dof_row_map_ = Teuchos::rcp(new Epetra_Map(*walldiscretization->DofRowMap()));
}
