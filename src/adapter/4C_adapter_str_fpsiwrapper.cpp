/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface


\level 3
*/

#include "4C_adapter_str_fpsiwrapper.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  bool PrestressIsActive(const double currentTime)
  {
    Inpar::Solid::PreStress pstype = Teuchos::getIntegralValue<Inpar::Solid::PreStress>(
        Global::Problem::instance()->structural_dynamic_params(), "PRESTRESS");
    const double pstime =
        Global::Problem::instance()->structural_dynamic_params().get<double>("PRESTRESSTIME");
    return pstype != Inpar::Solid::PreStress::none && currentTime <= pstime + 1.0e-15;
  }
}  // namespace


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FPSIStructureWrapper::FPSIStructureWrapper(Teuchos::RCP<Structure> structure)
    : FSIStructureWrapper(structure)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FPSIStructureWrapper::extract_interface_dispn(bool FPSI)
{
  if (!FPSI)
  {
    return Adapter::FSIStructureWrapper::extract_interface_dispn();
  }
  else
  {
    // prestressing business
    if (PrestressIsActive(time_old()))
    {
      return Teuchos::rcp(new Epetra_Vector(*interface_->fpsi_cond_map(), true));
    }
    else
    {
      return interface_->extract_fpsi_cond_vector(dispn());
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FPSIStructureWrapper::extract_interface_dispnp(bool FPSI)
{
  if (!FPSI)
  {
    return Adapter::FSIStructureWrapper::extract_interface_dispnp();
  }
  else
  {
    // prestressing business
    if (PrestressIsActive(time()))
    {
      return Teuchos::rcp(new Epetra_Vector(*interface_->fpsi_cond_map(), true));
    }
    else
    {
      return interface_->extract_fpsi_cond_vector(dispnp());
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
