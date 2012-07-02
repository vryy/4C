/*
 * ssi_partitioned.cpp
 *
 *  Created on: Jun 28, 2012
 *      Author: vuong
 */

#include "ssi_partitioned.H"
#include "../linalg/linalg_utils.H"

SSI::SSI_Part::SSI_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : SSI_Base(comm, timeparams)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetupSystem()
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetStructSolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetVelocityFields()
{

  scatra_->ScaTraField().SetVelocityField(
      zeros_, //convective vel.
      Teuchos::null, //acceleration
      structure_->ExtractVelnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      structure_->Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetMeshDisp()
{
  scatra_->ScaTraField().ApplyMeshMovement(
      structure_->Dispnp(),
      structure_->Discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void SSI::SSI_Part::TimeUpdateAndOutput()
{
  structure_->PrepareOutput();

  scatra_->ScaTraField().Update();

  structure_->Update();

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  structure_->Output();

  // output of fluid- and structure-based scalar transport
  scatra_->ScaTraField().Output();

  return;
}*/
