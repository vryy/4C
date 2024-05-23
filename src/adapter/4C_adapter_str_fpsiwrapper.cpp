/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface


\level 3
*/

#include "4C_adapter_str_fpsiwrapper.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  bool PrestressIsActive(const double currentTime)
  {
    INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(
        GLOBAL::Problem::Instance()->structural_dynamic_params(), "PRESTRESS");
    const double pstime =
        GLOBAL::Problem::Instance()->structural_dynamic_params().get<double>("PRESTRESSTIME");
    return pstype != INPAR::STR::PreStress::none && currentTime <= pstime + 1.0e-15;
  }
}  // namespace


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FPSIStructureWrapper::FPSIStructureWrapper(Teuchos::RCP<Structure> structure)
    : FSIStructureWrapper(structure)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FPSIStructureWrapper::extract_interface_dispn(bool FPSI)
{
  if (!FPSI)
  {
    return ADAPTER::FSIStructureWrapper::extract_interface_dispn();
  }
  else
  {
    // prestressing business
    if (PrestressIsActive(TimeOld()))
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
Teuchos::RCP<Epetra_Vector> ADAPTER::FPSIStructureWrapper::extract_interface_dispnp(bool FPSI)
{
  if (!FPSI)
  {
    return ADAPTER::FSIStructureWrapper::extract_interface_dispnp();
  }
  else
  {
    // prestressing business
    if (PrestressIsActive(Time()))
    {
      return Teuchos::rcp(new Epetra_Vector(*interface_->FPSICondMap(), true));
    }
    else
    {
      return interface_->ExtractFPSICondVector(Dispnp());
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
