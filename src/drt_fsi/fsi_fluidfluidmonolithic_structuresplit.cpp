/*------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi
(structuresplit) using XFEM and NOX

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*------------------------------------------------------*/

#include "fsi_fluidfluidmonolithic_structuresplit.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_fld_fluid_fluid_fsi.H"
#include "../drt_adapter/ad_ale_xffsi.H"

#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_ale.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidFluidMonolithicStructureSplit::FluidFluidMonolithicStructureSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicStructureSplit(comm, timeparams)
{
  // cast to problem-specific fluid-wrapper
  fluid_ =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidFluidFSI>(MonolithicStructureSplit::FluidField());

  // cast to problem-specific ALE-wrapper
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleXFFsiWrapper>(MonolithicStructureSplit::AleField());

  // XFFSI_Full_Newton is an invalid choice together with NOX,
  // because DOF-maps can change from one iteration step to the other (XFEM cut)
  if (FluidField()->MonolithicXffsiApproach() == INPAR::XFEM::XFFSI_Full_Newton)
    dserror("NOX-based XFFSI Approach does not work with XFFSI_Full_Newton!");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::Update()
{
  // time to relax the ALE-mesh?
  if (FluidField()->IsAleRelaxationStep(Step()))
  {
    if (Comm().MyPID() == 0) IO::cout << "Relaxing Ale" << IO::endl;

    AleField()->Solve();
    FluidField()->ApplyMeshDisplacement(AleToFluid(AleField()->Dispnp()));
  }

  // update fields
  FSI::MonolithicStructureSplit::Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::PrepareTimeStep()
{
  // prepare time step on subsequent field & increment
  FSI::MonolithicStructureSplit::PrepareTimeStep();

  // when this is the first call or we haven't relaxed the ALE-mesh
  // previously, the DOF-maps have not
  // changed since system setup
  if (Step() == 0 || !FluidField()->IsAleRelaxationStep(Step() - 1)) return;

  // REMARK:
  // as the new xfem-cut may lead to a change in the fluid dof-map,
  // we have to refresh the block system matrix,
  // rebuild the merged DOF map & update map extractor for combined
  // Dirichlet maps
  FSI::MonolithicStructureSplit::CreateCombinedDofRowMap();
  SetupDBCMapExtractor();
  FSI::MonolithicStructureSplit::CreateSystemMatrix();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupDBCMapExtractor()
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
      LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);
  dbcmaps.push_back(aleintersectionmap);

  Teuchos::RCP<const Epetra_Map> dbcmap = LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null)
  {
    dserror("Creation of Dirichlet map extractor failed.");
  }
}
