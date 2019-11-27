/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface

\maintainer Johannes Kremheller

\level 3
*/

#include "ad_str_fpsiwrapper.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../drt_structure/stru_aux.H"

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
    double time = 0.0;
    double pstime = -1.0;
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    INPAR::STR::PreStress pstype =
        DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
    if (pstype != INPAR::STR::prestress_none)
    {
      time = TimeOld();
      pstime = sdyn.get<double>("PRESTRESSTIME");
    }

    if (pstype != INPAR::STR::prestress_none && time <= pstime)
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
    double time = 0.0;
    double pstime = -1.0;
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    INPAR::STR::PreStress pstype =
        DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
    if (pstype != INPAR::STR::prestress_none)
    {
      time = Time();
      pstime = sdyn.get<double>("PRESTRESSTIME");
    }

    if (pstype != INPAR::STR::prestress_none && time <= pstime)
    {
      return Teuchos::rcp(new Epetra_Vector(*interface_->FPSICondMap(), true));
    }
    else
    {
      return interface_->ExtractFPSICondVector(Dispnp());
    }
  }
}
