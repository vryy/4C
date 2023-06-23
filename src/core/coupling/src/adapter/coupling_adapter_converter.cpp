/*----------------------------------------------------------------------------*/
/*! \file

\brief Converter to use ADAPTER::Coupling type objects in both coupling directions

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */
#include "coupling_adapter_converter.H"
#include "coupling_adapter.H"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> CORE::ADAPTER::CouplingMasterConverter::SrcToDst(
    Teuchos::RCP<const Epetra_Vector> source_vector) const
{
  return coup_.MasterToSlave(source_vector);
}

Teuchos::RCP<Epetra_Vector> CORE::ADAPTER::CouplingMasterConverter::DstToSrc(
    Teuchos::RCP<const Epetra_Vector> destination_vector) const
{
  return coup_.SlaveToMaster(destination_vector);
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingMasterConverter::SrcMap() const
{
  return coup_.MasterDofMap();
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingMasterConverter::DstMap() const
{
  return coup_.SlaveDofMap();
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingMasterConverter::PermSrcMap() const
{
  return coup_.PermMasterDofMap();
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingMasterConverter::PermDstMap() const
{
  return coup_.PermSlaveDofMap();
}

void CORE::ADAPTER::CouplingMasterConverter::FillSrcToDstMap(std::map<int, int>& rowmap) const
{
  coup_.FillMasterToSlaveMap(rowmap);
}


Teuchos::RCP<Epetra_Vector> CORE::ADAPTER::CouplingSlaveConverter::SrcToDst(
    Teuchos::RCP<const Epetra_Vector> source_vector) const
{
  return coup_.SlaveToMaster(source_vector);
}

Teuchos::RCP<Epetra_Vector> CORE::ADAPTER::CouplingSlaveConverter::DstToSrc(
    Teuchos::RCP<const Epetra_Vector> destination_vector) const
{
  return coup_.MasterToSlave(destination_vector);
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingSlaveConverter::SrcMap() const
{
  return coup_.SlaveDofMap();
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingSlaveConverter::DstMap() const
{
  return coup_.MasterDofMap();
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingSlaveConverter::PermSrcMap() const
{
  return coup_.PermSlaveDofMap();
}

Teuchos::RCP<const Epetra_Map> CORE::ADAPTER::CouplingSlaveConverter::PermDstMap() const
{
  return coup_.PermMasterDofMap();
}

void CORE::ADAPTER::CouplingSlaveConverter::FillSrcToDstMap(std::map<int, int>& rowmap) const
{
  coup_.FillSlaveToMasterMap(rowmap);
}
