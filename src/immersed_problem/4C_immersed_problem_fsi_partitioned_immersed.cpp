/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi subclass

\level 1


*----------------------------------------------------------------------*/
#include "4C_immersed_problem_fsi_partitioned_immersed.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN


FSI::PartitionedImmersed::PartitionedImmersed(const Epetra_Comm& comm) : Partitioned(comm)
{
  // empty constructor
}


void FSI::PartitionedImmersed::Setup()
{
  // call setup of base class
  FSI::Partitioned::Setup();
}


void FSI::PartitionedImmersed::SetupCoupling(
    const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  if (Comm().MyPID() == 0)
    std::cout << "\n SetupCoupling in FSI::PartitionedImmersed ..." << std::endl;

  // for immersed fsi
  coupsfm_ = Teuchos::null;
  matchingnodes_ = false;

  // enable debugging
  if (CORE::UTILS::IntegralValue<int>(fsidyn, "DEBUGOUTPUT"))
    debugwriter_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField()->Discretization()));
}


void FSI::PartitionedImmersed::extract_previous_interface_solution()
{
  // not necessary in immersed fsi.
  // overrides version in fsi_paritioned with "do nothing".
}

FOUR_C_NAMESPACE_CLOSE
