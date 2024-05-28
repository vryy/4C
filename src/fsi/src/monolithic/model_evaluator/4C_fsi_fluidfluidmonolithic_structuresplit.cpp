/*------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi
(structuresplit) using XFEM and NOX

\level 2

*/
/*------------------------------------------------------*/

#include "4C_fsi_fluidfluidmonolithic_structuresplit.hpp"

#include "4C_adapter_ale_xffsi.hpp"
#include "4C_adapter_fld_fluid_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_ale_utils_mapextractor.hpp"
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
FSI::FluidFluidMonolithicStructureSplit::FluidFluidMonolithicStructureSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicStructureSplit(comm, timeparams)
{
  // cast to problem-specific fluid-wrapper
  fluid_ =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidFluidFSI>(MonolithicStructureSplit::fluid_field());

  // cast to problem-specific ALE-wrapper
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleXFFsiWrapper>(MonolithicStructureSplit::ale_field());

  // XFFSI_Full_Newton is an invalid choice together with NOX,
  // because DOF-maps can change from one iteration step to the other (XFEM cut)
  if (fluid_field()->monolithic_xffsi_approach() == INPAR::XFEM::XFFSI_Full_Newton)
    FOUR_C_THROW("NOX-based XFFSI Approach does not work with XFFSI_Full_Newton!");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::update()
{
  // time to relax the ALE-mesh?
  if (fluid_field()->IsAleRelaxationStep(Step()))
  {
    if (Comm().MyPID() == 0) IO::cout << "Relaxing Ale" << IO::endl;

    ale_field()->Solve();
    fluid_field()->apply_mesh_displacement(ale_to_fluid(ale_field()->Dispnp()));
  }

  // update fields
  FSI::MonolithicStructureSplit::update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::prepare_time_step()
{
  // prepare time step on subsequent field & increment
  FSI::MonolithicStructureSplit::prepare_time_step();

  // when this is the first call or we haven't relaxed the ALE-mesh
  // previously, the DOF-maps have not
  // changed since system setup
  if (Step() == 0 || !fluid_field()->IsAleRelaxationStep(Step() - 1)) return;

  // REMARK:
  // as the new xfem-cut may lead to a change in the fluid dof-map,
  // we have to refresh the block system matrix,
  // rebuild the merged DOF map & update map extractor for combined
  // Dirichlet maps
  FSI::MonolithicStructureSplit::create_combined_dof_row_map();
  setup_dbc_map_extractor();
  FSI::MonolithicStructureSplit::create_system_matrix();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::setup_dbc_map_extractor()
{
  // merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;

  // structure DBC
  dbcmaps.push_back(structure_field()->GetDBCMapExtractor()->CondMap());
  // fluid DBC (including background & embedded discretization)
  dbcmaps.push_back(fluid_field()->GetDBCMapExtractor()->CondMap());
  // ALE-DBC-maps, free of FSI DOF
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(ale_field()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(ale_field()->Interface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);
  dbcmaps.push_back(aleintersectionmap);

  Teuchos::RCP<const Epetra_Map> dbcmap = CORE::LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(*dof_row_map(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null)
  {
    FOUR_C_THROW("Creation of Dirichlet map extractor failed.");
  }
}

FOUR_C_NAMESPACE_CLOSE
