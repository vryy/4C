/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of monolithic XFSI algorithm that performs a coupling between the
       structural field equation and XFEM fluid field equations

\level 2

*/


#include "4C_fsi_xfem_algorithm.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_adapter_ale_fpsi.hpp"
#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_poro_wrapper.hpp"
#include "4C_fluid_xfluid.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_poroelast_monolithic.hpp"

FOUR_C_NAMESPACE_OPEN



//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                    schott 08/14 |
 *----------------------------------------------------------------------*/
FSI::AlgorithmXFEM::AlgorithmXFEM(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
    const ADAPTER::FieldWrapper::Fieldtype type)
    : AlgorithmBase(comm, timeparams),
      num_fields_(0),
      structp_block_(-1),
      fluid_block_(-1),
      fluidp_block_(-1),
      ale_i_block_(-1)
{
  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  const Teuchos::ParameterList& sdyn = GLOBAL::Problem::Instance()->structural_dynamic_params();
  const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();
  const Teuchos::ParameterList& xfdyn = GLOBAL::Problem::Instance()->XFluidDynamicParams();
  bool ale = CORE::UTILS::IntegralValue<bool>((xfdyn.sublist("GENERAL")), "ALE_XFluid");

  num_fields_ += 2;
  structp_block_ = 0;
  fluid_block_ = 1;

  if (type == ADAPTER::StructurePoroWrapper::type_StructureField)
  {
    // ask base algorithm for the structural time integrator
    // access the structural discretization
    Teuchos::RCP<DRT::Discretization> structdis = GLOBAL::Problem::Instance()->GetDis("structure");
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
            timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
    structureporo_ = Teuchos::rcp(new ADAPTER::StructurePoroWrapper(
        structure->StructureField(), ADAPTER::StructurePoroWrapper::type_StructureField, true));
  }
  else if (type == ADAPTER::StructurePoroWrapper::type_PoroField)
  {
    num_fields_ += 1;
    fluidp_block_ = 2;

    GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
    const Teuchos::ParameterList& poroelastdyn =
        problem->poroelast_dynamic_params();  // access the problem-specific parameter list
    Teuchos::RCP<POROELAST::Monolithic> poro = Teuchos::rcp_dynamic_cast<POROELAST::Monolithic>(
        POROELAST::UTILS::CreatePoroAlgorithm(poroelastdyn, comm, false));
    if (poro == Teuchos::null)  // safety check
      FOUR_C_THROW(
          "Couldn't cast poro to POROELAST::Monolithic --> check your COUPALGO in the "
          "POROELASTICITY DYNAMIC section!");
    if (CORE::UTILS::IntegralValue<INPAR::POROELAST::SolutionSchemeOverFields>(
            poroelastdyn, "COUPALGO") != INPAR::POROELAST::Monolithic)
      FOUR_C_THROW(
          "You created a different poroelast algorithm than monolithic (not combineable with xfpsi "
          "at the moment)--> check your COUPALGO in the POROELASTICITY DYNAMIC section!");
    structureporo_ = Teuchos::rcp(new ADAPTER::StructurePoroWrapper(
        poro, ADAPTER::StructurePoroWrapper::type_PoroField, true));
  }
  else
    FOUR_C_THROW("AlgorithmXFEM cannot handle this Fieldtype for structure!");

  if (ale)
  {
    num_fields_ += 1;
    ale_i_block_ = num_fields_ - 1;
    GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
    const Teuchos::ParameterList& fsidynparams = problem->FSIDynamicParams();
    // ask base algorithm for the ale time integrator
    Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
        new ADAPTER::AleBaseAlgorithm(fsidynparams, GLOBAL::Problem::Instance()->GetDis("ale")));
    ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFpsiWrapper>(ale->ale_field());
    if (ale_ == Teuchos::null)
      FOUR_C_THROW("Cast from ADAPTER::Ale to ADAPTER::AleFpsiWrapper failed");
  }
  else
  {
    ale_ = Teuchos::null;
  }


  //--------------------------------------------
  // ask base algorithm for the fluid time integrator
  // do not init in ale case!!! (will be done in MonolithicAFSI_XFEM::Setup System())
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams, fdyn, "fluid", ale, false));
  fluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluid>(fluid->fluid_field());
  if (fluid_ == Teuchos::null)
    FOUR_C_THROW("Cast of Fluid to XFluid failed! - Everything fine in setup_fluid()?");
  fluid_->Init(false);

  // Do setup of the fields here
  structureporo_->Setup();
  return;
}



/*----------------------------------------------------------------------*
 | setup (public)                                            ager 12/16 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::Setup()
{
  // Do setup of the fields here
  // structureporo_->Setup();
}

/*----------------------------------------------------------------------*
 | update (protected)                                      schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::Update()
{
  FOUR_C_THROW("currently unused");

  StructurePoro()->Update();
  fluid_field()->Update();
  if (HaveAle()) ale_field()->Update();

  return;
}



/*----------------------------------------------------------------------*
 | calculate stresses, strains, energies                   schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::prepare_output(bool force_prepare)
{
  StructurePoro()->prepare_output(force_prepare);
}

FOUR_C_NAMESPACE_CLOSE
