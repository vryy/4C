/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface


\level 3
*/

#include "baci_adapter_str_fpsiwrapper.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_structure_aux.H"

#include "baci_lib_prestress_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FPSIStructureWrapper::FPSIStructureWrapper(Teuchos::RCP<Structure> structure)
    : FSIStructureWrapper(structure)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FPSIStructureWrapper::ExtractInterfaceDispn(bool FPSI)
{
  if (!FPSI)
  {
    return ADAPTER::FSIStructureWrapper::ExtractInterfaceDispn();
  }
  else
  {
    // prestressing business
    if (UTILS::PRESTRESS::IsActive(TimeOld()))
    {
      return Teuchos::rcp(new Epetra_Vector(*interface_->FPSICondMap(), true));
    }
    else
    {
      return interface_->ExtractFPSICondVector(Dispn());
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FPSIStructureWrapper::ExtractInterfaceDispnp(bool FPSI)
{
  if (!FPSI)
  {
    return ADAPTER::FSIStructureWrapper::ExtractInterfaceDispnp();
  }
  else
  {
    // prestressing business
    if (UTILS::PRESTRESS::IsActive(Time()))
    {
      return Teuchos::rcp(new Epetra_Vector(*interface_->FPSICondMap(), true));
    }
    else
    {
      return interface_->ExtractFPSICondVector(Dispnp());
    }
  }
}
