/*----------------------------------------------------------------------------*/
/*! \file

\brief Converter to use Adapter::Coupling type objects in both coupling directions

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */
#include "4C_coupling_adapter_converter.hpp"

#include "4C_coupling_adapter.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> Core::Adapter::CouplingMasterConverter::src_to_dst(
    Teuchos::RCP<const Epetra_Vector> source_vector) const
{
  return coup_.master_to_slave(source_vector);
}

Teuchos::RCP<Epetra_Vector> Core::Adapter::CouplingMasterConverter::dst_to_src(
    Teuchos::RCP<const Epetra_Vector> destination_vector) const
{
  return coup_.slave_to_master(destination_vector);
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingMasterConverter::src_map() const
{
  return coup_.master_dof_map();
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingMasterConverter::dst_map() const
{
  return coup_.slave_dof_map();
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingMasterConverter::perm_src_map() const
{
  return coup_.perm_master_dof_map();
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingMasterConverter::perm_dst_map() const
{
  return coup_.perm_slave_dof_map();
}

void Core::Adapter::CouplingMasterConverter::fill_src_to_dst_map(std::map<int, int>& rowmap) const
{
  coup_.fill_master_to_slave_map(rowmap);
}


Teuchos::RCP<Epetra_Vector> Core::Adapter::CouplingSlaveConverter::src_to_dst(
    Teuchos::RCP<const Epetra_Vector> source_vector) const
{
  return coup_.slave_to_master(source_vector);
}

Teuchos::RCP<Epetra_Vector> Core::Adapter::CouplingSlaveConverter::dst_to_src(
    Teuchos::RCP<const Epetra_Vector> destination_vector) const
{
  return coup_.master_to_slave(destination_vector);
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingSlaveConverter::src_map() const
{
  return coup_.slave_dof_map();
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingSlaveConverter::dst_map() const
{
  return coup_.master_dof_map();
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingSlaveConverter::perm_src_map() const
{
  return coup_.perm_slave_dof_map();
}

Teuchos::RCP<const Epetra_Map> Core::Adapter::CouplingSlaveConverter::perm_dst_map() const
{
  return coup_.perm_master_dof_map();
}

void Core::Adapter::CouplingSlaveConverter::fill_src_to_dst_map(std::map<int, int>& rowmap) const
{
  coup_.fill_slave_to_master_map(rowmap);
}

FOUR_C_NAMESPACE_CLOSE
