/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fpsi. Can only be used in conjunction with #FluidImplicitTimeInt

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_fpsi.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fpsi_utils.hpp"

FOUR_C_NAMESPACE_OPEN


/* constructor */
Adapter::FluidFPSI::FluidFPSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Core::FE::Discretization> dis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond),
      fpsiinterface_(Teuchos::rcp(new FLD::UTILS::MapExtractor()))
{
  return;
}  // constructor


/* initialization */
void Adapter::FluidFPSI::init()
{
  // call base class init
  FluidFSI::init();

  fpsiinterface_->setup(*dis_, true, true);  // Always Create overlapping FPSI Interface

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFPSI::setup_interface(const int nds_master)
{
  // check nds_master
  if (nds_master != 0) FOUR_C_THROW("nds_master is supposed to be 0 here");

  interface_->setup(*dis_, false, true);  // create overlapping maps for fpsi problem
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFPSI::use_block_matrix(
    bool splitmatrix, Teuchos::RCP<FPSI::UTILS::MapExtractor> const& shapederivSplitter)
{
  Teuchos::RCP<std::set<int>> condelements =
      Interface()->conditioned_element_map(*discretization());
  Teuchos::RCP<std::set<int>> condelements_shape =
      shapederivSplitter->conditioned_element_map(*discretization());
  fluidimpl_->use_block_matrix(condelements, *Interface(), *Interface(), condelements_shape,
      *shapederivSplitter, *shapederivSplitter, splitmatrix);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFPSI::use_block_matrix(bool splitmatrix)
{
  Teuchos::RCP<std::set<int>> condelements =
      Interface()->conditioned_element_map(*discretization());
  fluidimpl_->use_block_matrix(condelements, *Interface(), *Interface(), condelements, *Interface(),
      *Interface(), splitmatrix);
}

FOUR_C_NAMESPACE_CLOSE
