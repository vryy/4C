/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed fluids

\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_immersed.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidImmersed::FluidImmersed(const Teuchos::ParameterList& prbdyn, std::string condname)
{
  fluid_ = Teuchos::rcp(new FluidBaseAlgorithm(
      prbdyn, GLOBAL::Problem::Instance()->FluidDynamicParams(), "fluid", false));
  fluidadapter_ = fluid_->fluid_field();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidImmersed::Discretization()
{
  return fluid_field()->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& ADAPTER::FluidImmersed::Interface() const
{
  return fluidadapter_->Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::prepare_time_step() { fluid_field()->prepare_time_step(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::Update() { fluid_field()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::Output() { fluid_field()->StatisticsAndOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidImmersed::read_restart(int step)
{
  fluid_field()->read_restart(step);
  return fluid_field()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::NonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  fluid_field()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::extract_interface_forces()
{
  return fluid_field()->extract_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::extract_interface_velnp()
{
  FOUR_C_THROW("Robin stuff");
  return fluid_field()->extract_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::extract_interface_veln()
{
  return fluid_field()->extract_interface_veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::integrate_interface_shape()
{
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CORE::UTILS::ResultTest> ADAPTER::FluidImmersed::CreateFieldTest()
{
  return fluid_field()->CreateFieldTest();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  fluid_field()->AddDirichCond(maptoadd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  fluid_field()->RemoveDirichCond(maptoremove);
}

FOUR_C_NAMESPACE_CLOSE
