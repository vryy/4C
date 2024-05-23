/*------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi
(fluidsplit) using XFEM and NOX

\level 3

*/
/*------------------------------------------------------*/

#include "4C_fsi_fluidfluidmonolithic_fluidsplit.hpp"

#include "4C_adapter_ale_xffsi.hpp"
#include "4C_adapter_fld_fluid_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidFluidMonolithicFluidSplit::FluidFluidMonolithicFluidSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicFluidSplit(comm, timeparams)
{
  // cast to problem-specific fluid-wrapper
  fluid_ = Teuchos::rcp_dynamic_cast<ADAPTER::FluidFluidFSI>(MonolithicFluidSplit::FluidField());

  // cast to problem-specific ALE-wrapper
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleXFFsiWrapper>(MonolithicFluidSplit::AleField());

  // XFFSI_Full_Newton is an invalid choice together with NOX,
  // because DOF-maps can change from one iteration step to the other (XFEM cut)
  if (FluidField()->monolithic_xffsi_approach() == INPAR::XFEM::XFFSI_Full_Newton)
    FOUR_C_THROW("NOX-based XFFSI Approach does not work with XFFSI_Full_Newton!");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplit::Update()
{
  // time to relax the ALE-mesh?
  if (FluidField()->IsAleRelaxationStep(Step()))
  {
    if (Comm().MyPID() == 0) IO::cout << "Relaxing Ale" << IO::endl;

    AleField()->Solve();
    FluidField()->apply_mesh_displacement(AleToFluid(AleField()->Dispnp()));
  }

  // update fields
  FSI::MonolithicFluidSplit::Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplit::PrepareTimeStep()
{
  // prepare time step on subsequent field & increment
  FSI::MonolithicFluidSplit::PrepareTimeStep();

  // when this is the first call or we haven't relaxed the ALE-mesh
  // previously, the DOF-maps have not
  // changed since system setup
  if (Step() == 0 || !FluidField()->IsAleRelaxationStep(Step() - 1)) return;

  // as the new xfem-cut may lead to a change in the fluid dof-map,
  // we have to refresh the block system matrix,
  // rebuild the merged DOF map & update map extractor for combined
  // Dirichlet maps
  FSI::MonolithicFluidSplit::create_combined_dof_row_map();
  setup_dbc_map_extractor();
  FSI::MonolithicFluidSplit::CreateSystemMatrix();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplit::setup_dbc_map_extractor()
{
  // merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;

  // structure DBC
  dbcmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  // fluid DBC (including background & embedded discretization)
  dbcmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  // ALE-DBC-maps, free of FSI DOF
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(AleField()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(AleField()->Interface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);
  dbcmaps.push_back(aleintersectionmap);

  Teuchos::RCP<const Epetra_Map> dbcmap = CORE::LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(*DofRowMap(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null)
  {
    FOUR_C_THROW("Creation of Dirichlet map extractor failed.");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplit::Output()
{
  StructureField()->Output();
  FluidField()->Output();

  // output Lagrange multiplier
  {
    // the Lagrange multiplier lives on the FSI interface
    // for output, we want to insert lambda into a full vector, defined on the embedded fluid field
    // 1. insert into vector containing all fluid DOF
    Teuchos::RCP<Epetra_Vector> lambdafull =
        FluidField()->Interface()->InsertFSICondVector(FSI::MonolithicFluidSplit::GetLambda());
    // 2. extract the embedded fluid part
    Teuchos::RCP<Epetra_Vector> lambdaemb =
        FluidField()->x_fluid_fluid_map_extractor()->ExtractFluidVector(lambdafull);

    const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
    const int upres = fsidyn.get<int>("RESULTSEVRY");
    if ((uprestart != 0 && FluidField()->Step() % uprestart == 0) ||
        FluidField()->Step() % upres == 0)
      FluidField()->DiscWriter()->WriteVector("fsilambda", lambdaemb);
  }
  AleField()->Output();

  if (StructureField()->get_constraint_manager()->HaveMonitor())
  {
    StructureField()->get_constraint_manager()->compute_monitor_values(StructureField()->Dispnp());
    if (Comm().MyPID() == 0) StructureField()->get_constraint_manager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplit::ReadRestart(int step)
{
  // Read Lagrange Multiplier (associated with embedded fluid)
  {
    Teuchos::RCP<Epetra_Vector> lambdaemb = Teuchos::rcp(
        new Epetra_Vector(*(FluidField()->x_fluid_fluid_map_extractor()->FluidMap()), true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(
        FluidField()->Discretization(), GLOBAL::Problem::Instance()->InputControlFile(), step);
    reader.ReadVector(lambdaemb, "fsilambda");
    // Insert into vector containing the whole merged fluid DOF
    Teuchos::RCP<Epetra_Vector> lambdafull =
        FluidField()->x_fluid_fluid_map_extractor()->InsertFluidVector(lambdaemb);
    FSI::MonolithicFluidSplit::SetLambda(
        FluidField()->Interface()->ExtractFSICondVector(lambdafull));
  }

  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);
  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), FluidField()->Step());
}

FOUR_C_NAMESPACE_CLOSE
