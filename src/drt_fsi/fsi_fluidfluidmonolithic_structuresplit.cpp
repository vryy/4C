/*------------------------------------------------------*/
/*!
\file fsi_fluidfluidmonolithic_structuresplit.cpp
\brief Control routine for monolithic fluid-fluid-fsi
(structuresplit) using XFEM and NOX

<pre>
Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            089 289 15249
</pre>
*/
/*------------------------------------------------------*/

#include "fsi_fluidfluidmonolithic_structuresplit.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale_new/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_ale_fsi.H"

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
FSI::FluidFluidMonolithicStructureSplit::FluidFluidMonolithicStructureSplit(const Epetra_Comm& comm,
                                                                            const Teuchos::ParameterList& timeparams)
  : MonolithicStructureSplit(comm,timeparams)
{
  // determine the type of monolithic approach
  const Teuchos::ParameterList& xfluiddyn  = DRT::Problem::Instance()->XFluidDynamicParams();
  enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>
               (xfluiddyn.sublist("GENERAL"),"MONOLITHIC_XFFSI_APPROACH");

  // XFFSI_Full_Newton is an invalid choice together with NOX,
  // because DOF-maps can change from one iteration step to the other (XFEM cut)
  if (monolithic_approach == INPAR::XFEM::XFFSI_Full_Newton)
    dserror("NOX-based XFFSI Approach does not work with XFFSI_Full_Newton!");

  // should ALE-relaxation be carried out?
  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(xfluiddyn.sublist("GENERAL"),"RELAXING_ALE");
  // get no. of timesteps, after which ALE-mesh should be relaxed
  relaxing_ale_every_ = xfluiddyn.sublist("GENERAL").get<int>("RELAXING_ALE_EVERY");

  if (! relaxing_ale_ && relaxing_ale_every_ != 0)
    dserror("You don't want to relax the ALE but provide a relaxation interval != 0 ?!");

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::Update()
{
  bool relaxing_ale = (relaxing_ale_ && relaxing_ale_every_ != 0) ? (Step() % relaxing_ale_every_ == 0) : false;

  if (relaxing_ale)
  {
    FluidField().ApplyEmbFixedMeshDisplacement(AleToFluid(AleField()->WriteAccessDispnp()));

    if (Comm().MyPID() == 0)
      IO::cout << "Relaxing Ale" << IO::endl;

    AleField()->SolveAleXFluidFluidFSI();
    FluidField().ApplyMeshDisplacement(AleToFluid(AleField()->WriteAccessDispnp()));
  }

  // update fields
  FSI::MonolithicStructureSplit::Update();

  // build ale system matrix for the next time step. Here first we
  // update the vectors then we set the fluid-fluid dirichlet values
  // in buildsystemmatrix
  if (relaxing_ale)
  {
    AleField()->CreateSystemMatrix(AleField()->Interface());
  }
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
  if (Step() == 0 || !relaxing_ale_)
    return;

  // rebuild maps and reset fluid matrix, if we relaxed the ALE-mesh in
  // the previous step
  if (relaxing_ale_every_ < 1)
    dserror("You want to relax the ALE-mesh, but provide a relaxation interval of %d?!", relaxing_ale_every_);

  // previous step was no relaxation step? leave!
  if ((Step()-1) % relaxing_ale_every_ != 0)
    return;

  // REMARK:
  // as the new xfem-cut may lead to a change in the fluid dof-map,
  // we have to refresh the block system matrix:
  // 1. rebuild the merged DOF map & update map extractor for combined
  // Dirichlet maps
  //FSI::MonolithicStructureSplit::CreateCombinedDofRowMap();
  //SetupDBCMapExtractor();
  // 2. refresh the fsi block system matrix
  // (structure and ALE-maps are unaffected)
  // Todo:
  // currently, there is no option in LINALG::BlockSparseMatrixBase, that allows
  // reallocation of single SparseMatrix blocks and change of the underlying MultiMapExtractor
  // something like: SystemMatrix()->RefreshBlocks(Extractor());
  // so as a temporary solution we reallocate the whole system matrix:
  // together with step 1, this is accomplished in SetupSystem(),
  // where also the costly setup of ADAPTER::Coupling takes place
  FSI::MonolithicStructureSplit::SetupSystem();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplit::SetupDBCMapExtractor()
{
  // merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map> > dbcmaps;

  // structure DBC
  dbcmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  // fluid DBC (including background & embedded discretization)
  dbcmaps.push_back(FluidField().GetDBCMapExtractor()->CondMap());
  // ALE-DBC-maps, free of FSI DOF
  std::vector<Teuchos::RCP<const Epetra_Map> > aleintersectionmaps;
  aleintersectionmaps.push_back(AleField()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(AleField()->Interface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap = LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);
  dbcmaps.push_back(aleintersectionmap);

  Teuchos::RCP<const Epetra_Map> dbcmap = LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(),dbcmap,true));
  if (dbcmaps_ == Teuchos::null) { dserror("Creation of Dirichlet map extractor failed."); }
}

