/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fpsi.cpp

\brief Fluid field adapter for fpsi

Can only be used in conjunction with #FluidImplicitTimeInt

\level 2

\maintainer  Christoph Ager
*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fluid_fpsi.H"
#include "ad_fld_fluid_fsi.H"

#include "../drt_fpsi/fpsi_utils.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"


/* constructor */
ADAPTER::FluidFPSI::FluidFPSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond),
      fpsiinterface_(Teuchos::rcp(new FLD::UTILS::MapExtractor()))
{
  return;
}  // constructor


/* initialization */
void ADAPTER::FluidFPSI::Init()
{
  // call base class init
  FluidFSI::Init();

  fpsiinterface_->Setup(*dis_, true, true);  // Always Create overlapping FPSI Interface

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFPSI::SetupInterface()
{
  interface_->Setup(*dis_, false, true);  // create overlapping maps for fpsi problem
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFPSI::UseBlockMatrix(
    bool splitmatrix, Teuchos::RCP<FPSI::UTILS::MapExtractor> const& shapederivSplitter)
{
  Teuchos::RCP<std::set<int>> condelements = Interface()->ConditionedElementMap(*Discretization());
  Teuchos::RCP<std::set<int>> condelements_shape =
      shapederivSplitter->ConditionedElementMap(*Discretization());
  fluidimpl_->UseBlockMatrix(condelements, *Interface(), *Interface(), condelements_shape,
      *shapederivSplitter, *shapederivSplitter, splitmatrix);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFPSI::UseBlockMatrix(bool splitmatrix)
{
  Teuchos::RCP<std::set<int>> condelements = Interface()->ConditionedElementMap(*Discretization());
  fluidimpl_->UseBlockMatrix(condelements, *Interface(), *Interface(), condelements, *Interface(),
      *Interface(), splitmatrix);
}
