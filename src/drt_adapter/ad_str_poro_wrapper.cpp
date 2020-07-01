/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for structure or poro time integration

\level 2


*/
/*----------------------------------------------------------------------*/

#include "ad_str_poro_wrapper.H"

#include "ad_str_fpsiwrapper.H"
#include "ad_fld_poro.H"
#include "../drt_poroelast/poroelast_monolithic.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure/stru_aux.H"

/// constructor
ADAPTER::StructurePoroWrapper::StructurePoroWrapper(
    Teuchos::RCP<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection)
    : FieldWrapper(field, type, NOXCorrection)
{
  switch (type_)
  {
    case FieldWrapper::type_StructureField:
      structure_ = Teuchos::rcp_dynamic_cast<FSIStructureWrapper>(field_);
      if (structure_ == Teuchos::null)
        dserror("StructurePoroWrapper: Cast from Field to FSIStructureWrapper failed!");
      poro_ = Teuchos::null;
      break;
    case FieldWrapper::type_PoroField:
      poro_ = Teuchos::rcp_dynamic_cast<POROELAST::Monolithic>(field_);
      if (poro_ == Teuchos::null)
        dserror("StructurePoroWrapper: Cast from Field to PoroBase failed!");
      structure_ = poro_->StructureField();
      break;
    default:
      dserror("StructurePoroWrapper - FieldWrapper::Fieldtype not available for this wrapper!");
      break;
  }
}

/// setup
void ADAPTER::StructurePoroWrapper::Setup()
{
  structure_->Setup();
  if (type_ == FieldWrapper::type_PoroField)
  {
    poro_->SetupSystem();
    poro_->SetupNewton();  // just to avoid modifications in poro (this sets iterinc_ there)
  }
}

//! unique map of all dofs that should be constrained with DBC
Teuchos::RCP<const Epetra_Map> ADAPTER::StructurePoroWrapper::CombinedDBCMap()
{
  switch (type_)
  {
    case FieldWrapper::type_StructureField:
      return structure_->GetDBCMapExtractor()->CondMap();
      break;
    case FieldWrapper::type_PoroField:
      return poro_->CombinedDBCMap();
      break;
    default:
      dserror("StructurePoroWrapper: type for this wrapper not considered!");
      return Teuchos::null;
      break;
  }
}

//   //! perform result test
void ADAPTER::StructurePoroWrapper::TestResults(DRT::Problem* problem)
{
  problem->AddFieldTest(structure_->CreateFieldTest());

  if (type_ == FieldWrapper::type_PoroField)
    problem->AddFieldTest(poro_->FluidField()->CreateFieldTest());
}

const Teuchos::RCP<POROELAST::Monolithic>& ADAPTER::StructurePoroWrapper::PoroField()
{
  if (type_ == ADAPTER::FieldWrapper::type_PoroField)
    return poro_;
  else
    dserror("StructurePoroWrapper - Field not a PoroField!");
  return poro_;  // do not remove dserror!!! - return just to make complier happy :-)
}

const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& ADAPTER::StructurePoroWrapper::StructureField()
{
  if (type_ == FieldWrapper::type_PoroField || type_ == FieldWrapper::type_StructureField)
    return structure_;
  else
    dserror("StructurePoroWrapper - Field not Structural- or Poro-Field!");
  return structure_;  // do not remove dserror!!! - return just to make complier happy :-)
}

//! return poro FluidField
const Teuchos::RCP<ADAPTER::FluidPoro>& ADAPTER::StructurePoroWrapper::FluidField()
{
  if (type_ == FieldWrapper::type_PoroField)
    return poro_->FluidField();
  else
    dserror("StructurePoroWrapper - Field not PoroField (no poro fluid field!");
  return poro_->FluidField();  // do not remove dserror!!! - return just to make complier happy :-)
}

//! Insert FSI Condition Vector
Teuchos::RCP<Epetra_Vector> ADAPTER::StructurePoroWrapper::InsertFSICondVector(
    Teuchos::RCP<const Epetra_Vector> cond)
{
  Teuchos::RCP<Epetra_Vector> tmpcond;
  switch (type_)
  {
    case FieldWrapper::type_StructureField:
      return Interface()->InsertFSICondVector(cond);
      break;
    case FieldWrapper::type_PoroField:
      tmpcond = Interface()->InsertFSICondVector(cond);
      return poro_->Extractor()->InsertVector(tmpcond, 0);  // into structural part = 0
      break;
    default:
      dserror("StructurePoroWrapper: type for this wrapper not considered!");
      return Teuchos::null;
      break;
  }
}

//! Recover Lagrange Multiplier during iteration (does nothing for structure)
void ADAPTER::StructurePoroWrapper::RecoverLagrangeMultiplierAfterNewtonStep(
    Teuchos::RCP<Epetra_Vector> iterinc)
{
  if (type_ == FieldWrapper::type_PoroField)
    poro_->RecoverLagrangeMultiplierAfterNewtonStep(iterinc);
}
