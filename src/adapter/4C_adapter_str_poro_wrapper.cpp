/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for structure or poro time integration

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_str_poro_wrapper.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/// constructor
Adapter::StructurePoroWrapper::StructurePoroWrapper(
    Teuchos::RCP<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection)
    : FieldWrapper(field, type, NOXCorrection)
{
  switch (type_)
  {
    case FieldWrapper::type_StructureField:
      structure_ = Teuchos::rcp_dynamic_cast<FSIStructureWrapper>(field_);
      if (structure_ == Teuchos::null)
        FOUR_C_THROW("StructurePoroWrapper: Cast from Field to FSIStructureWrapper failed!");
      poro_ = Teuchos::null;
      break;
    case FieldWrapper::type_PoroField:
      poro_ = Teuchos::rcp_dynamic_cast<PoroElast::Monolithic>(field_);
      if (poro_ == Teuchos::null)
        FOUR_C_THROW("StructurePoroWrapper: Cast from Field to PoroBase failed!");
      structure_ = poro_->structure_field();
      break;
    default:
      FOUR_C_THROW(
          "StructurePoroWrapper - FieldWrapper::Fieldtype not available for this wrapper!");
      break;
  }
}

/// setup
void Adapter::StructurePoroWrapper::setup()
{
  structure_->setup();
  if (type_ == FieldWrapper::type_PoroField)
  {
    poro_->SetupSystem();
    poro_->SetupNewton();  // just to avoid modifications in poro (this sets iterinc_ there)
  }
}

//! unique map of all dofs that should be constrained with DBC
Teuchos::RCP<const Epetra_Map> Adapter::StructurePoroWrapper::combined_dbc_map()
{
  switch (type_)
  {
    case FieldWrapper::type_StructureField:
      return structure_->GetDBCMapExtractor()->CondMap();
      break;
    case FieldWrapper::type_PoroField:
      return poro_->combined_dbc_map();
      break;
    default:
      FOUR_C_THROW("StructurePoroWrapper: type for this wrapper not considered!");
      return Teuchos::null;
      break;
  }
}

//   //! perform result test
void Adapter::StructurePoroWrapper::TestResults(Global::Problem* problem)
{
  problem->AddFieldTest(structure_->CreateFieldTest());

  if (type_ == FieldWrapper::type_PoroField)
    problem->AddFieldTest(poro_->fluid_field()->CreateFieldTest());
}

const Teuchos::RCP<PoroElast::Monolithic>& Adapter::StructurePoroWrapper::poro_field()
{
  if (type_ == Adapter::FieldWrapper::type_PoroField)
    return poro_;
  else
    FOUR_C_THROW("StructurePoroWrapper - Field not a poro_field!");
  return poro_;  // do not remove FOUR_C_THROW!!! - return just to make complier happy :-)
}

const Teuchos::RCP<Adapter::FSIStructureWrapper>& Adapter::StructurePoroWrapper::structure_field()
{
  if (type_ == FieldWrapper::type_PoroField || type_ == FieldWrapper::type_StructureField)
    return structure_;
  else
    FOUR_C_THROW("StructurePoroWrapper - Field not Structural- or Poro-Field!");
  return structure_;  // do not remove FOUR_C_THROW!!! - return just to make complier happy :-)
}

//! return poro fluid_field
const Teuchos::RCP<Adapter::FluidPoro>& Adapter::StructurePoroWrapper::fluid_field()
{
  if (type_ == FieldWrapper::type_PoroField)
    return poro_->fluid_field();
  else
    FOUR_C_THROW("StructurePoroWrapper - Field not poro_field (no poro fluid field!");
  return poro_
      ->fluid_field();  // do not remove FOUR_C_THROW!!! - return just to make complier happy :-)
}

//! Insert FSI Condition Vector
Teuchos::RCP<Epetra_Vector> Adapter::StructurePoroWrapper::InsertFSICondVector(
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
      FOUR_C_THROW("StructurePoroWrapper: type for this wrapper not considered!");
      return Teuchos::null;
      break;
  }
}

//! Recover Lagrange Multiplier during iteration (does nothing for structure)
void Adapter::StructurePoroWrapper::recover_lagrange_multiplier_after_newton_step(
    Teuchos::RCP<Epetra_Vector> iterinc)
{
  if (type_ == FieldWrapper::type_PoroField)
    poro_->recover_lagrange_multiplier_after_newton_step(iterinc);
}

FOUR_C_NAMESPACE_CLOSE
