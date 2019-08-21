/*----------------------------------------------------------------------*/
/*! \file

 \brief porous medium algorithm with block matrices for splitting and condensation

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/

#include "poro_monolithicsplit.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_blocksparsematrix.H"

POROELAST::MonolithicSplit::MonolithicSplit(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : Monolithic(comm, timeparams)
{
  icoupfs_ = Teuchos::rcp(new ADAPTER::Coupling());

  evaluateinterface_ = false;

  fgcur_ = Teuchos::null;
  ddiinc_ = Teuchos::null;
  solipre_ = Teuchos::null;
  ddginc_ = Teuchos::null;
  solgpre_ = Teuchos::null;
  duiinc_ = Teuchos::null;
  solivelpre_ = Teuchos::null;
  duginc_ = Teuchos::null;
  solgvelpre_ = Teuchos::null;

  fsibcmap_ = Teuchos::null;
  fsibcextractor_ = Teuchos::null;

  ddi_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplit::PrepareTimeStep()
{
  // call base class
  POROELAST::Monolithic::PrepareTimeStep();

  // counter and print header
  // IncrementTimeAndStep();
  // PrintHeader();

  //  // call the predictor
  //  StructureField()->PrepareTimeStep();
  //  FluidField()->PrepareTimeStep();

  if (evaluateinterface_)
  {
    // here we account for DBCs and preconditioning on the FSI-Interface. In both cases the
    // structure field decides, what to do (I don't think this is the best solution, but at least
    // the
    // easiest one)

    double timescale = FluidField()->TimeScaling();

    Teuchos::RCP<const Epetra_Vector> idispnp =
        StructureField()->Interface()->ExtractFSICondVector(StructureField()->Dispnp());
    Teuchos::RCP<const Epetra_Vector> idispn =
        StructureField()->Interface()->ExtractFSICondVector(StructureField()->Dispn());
    Teuchos::RCP<const Epetra_Vector> ivelnp =
        StructureField()->Interface()->ExtractFSICondVector(StructureField()->Velnp());
    Teuchos::RCP<Epetra_Vector> ifvelnp = FluidField()->ExtractInterfaceVelnp();
    Teuchos::RCP<Epetra_Vector> ifveln = FluidField()->ExtractInterfaceVeln();

    ddi_->Update(1.0, *idispnp, -1.0, *idispn, 0.0);
    ddi_->Update(-1.0, *ifveln, timescale);

    if (fsibcmap_->NumGlobalElements())
    {
      // if there are DBCs on FSI conditioned nodes, they have to be treated seperately

      Teuchos::RCP<Epetra_Vector> ibcveln =
          fsibcextractor_->ExtractCondVector(StructureToFluidAtInterface(ivelnp));
      Teuchos::RCP<Epetra_Vector> inobcveln =
          fsibcextractor_->ExtractOtherVector(StructureToFluidAtInterface(ddi_));

      // DBCs at FSI-Interface
      fsibcextractor_->InsertCondVector(ibcveln, ifvelnp);
      // any preconditioned values at the FSI-Interface
      fsibcextractor_->InsertOtherVector(inobcveln, ifvelnp);
    }
    else
      // no DBCs on FSI interface -> just make preconditioners consistent (structure decides)
      ifvelnp = StructureToFluidAtInterface(ddi_);

    FluidField()->ApplyInterfaceVelocities(ifvelnp);
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::MonolithicSplit::StructureToFluidAtInterface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfs_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::MonolithicSplit::FluidToStructureAtInterface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfs_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC and FSI Coupling Condition on structure side
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> POROELAST::MonolithicSplit::FSIDBCMap()
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicSplit::FSIDBCMap");

  // get interface map and DBC map of fluid
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidmaps;
  fluidmaps.push_back(FluidField()->Interface()->FSICondMap());
  fluidmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());

  // build vector of dbc and fsi coupling of fluid field
  Teuchos::RCP<Epetra_Map> fluidfsibcmap = LINALG::MultiMapExtractor::IntersectMaps(fluidmaps);

  if (fluidfsibcmap->NumMyElements())
    dserror("Dirichlet boundary conditions on fluid and FSI interface not supported!!");

  // get interface map and DBC map of structure
  std::vector<Teuchos::RCP<const Epetra_Map>> structmaps;
  structmaps.push_back(StructureField()->Interface()->FSICondMap());
  structmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());

  // vector of dbc and fsi coupling of structure field
  Teuchos::RCP<Epetra_Map> structfsibcmap = LINALG::MultiMapExtractor::IntersectMaps(structmaps);

  Teuchos::RCP<Epetra_Vector> gidmarker_struct =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(), true));

  // Todo this is ugly, fix it!
  const int mylength = structfsibcmap->NumMyElements();  // on each processor (lids)
  const int* mygids = structfsibcmap->MyGlobalElements();

  // mark gids with fsi and DBC Condition
  for (int i = 0; i < mylength; ++i)
  {
    int gid = mygids[i];
    // dsassert(slavemastermap.count(gid),"master gid not found on slave side");
    int err = gidmarker_struct->ReplaceGlobalValue(gid, 0, 1.0);
    if (err) dserror("ReplaceMyValue failed for gid %i error code %d", gid, err);
  }

  // transfer to fluid side
  Teuchos::RCP<Epetra_Vector> gidmarker_fluid = StructureToFluidAtInterface(gidmarker_struct);

  std::vector<int> structfsidbcvector;
  const int numgids = gidmarker_fluid->MyLength();  // on each processor (lids)
  double* mygids_fluid = gidmarker_fluid->Values();
  const int* fluidmap = gidmarker_fluid->Map().MyGlobalElements();
  for (int i = 0; i < numgids; ++i)
  {
    double val = mygids_fluid[i];
    if (val == 1.0) structfsidbcvector.push_back(fluidmap[i]);
  }

  Teuchos::RCP<Epetra_Map> structfsidbcmap = Teuchos::rcp(
      new Epetra_Map(-1, structfsidbcvector.size(), &structfsidbcvector[0], 0, Comm()));
  // dsassert(fluidfsidbcmap->UniqueGIDs(),"fsidbcmap is not unique!");

  return structfsidbcmap;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 11/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplit::SetupCouplingAndMatrices()
{
  const int ndim = DRT::Problem::Instance()->NDim();
  icoupfs_->SetupConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->FSICondMap(), *FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  fsibcmap_ = FSIDBCMap();

  evaluateinterface_ = StructureField()->Interface()->FSICondRelevant();

  if (evaluateinterface_)
  {
    if (fsibcmap_->NumGlobalElements())
    {
      const Teuchos::RCP<ADAPTER::FluidPoro>& fluidfield =
          Teuchos::rcp_dynamic_cast<ADAPTER::FluidPoro>(FluidField());
      fluidfield->AddDirichCond(fsibcmap_);

      fsibcextractor_ = Teuchos::rcp(
          new LINALG::MapExtractor(*FluidField()->Interface()->FSICondMap(), fsibcmap_));
    }

    Teuchos::RCP<const Epetra_Vector> idispnp =
        StructureField()->Interface()->ExtractFSICondVector(StructureField()->Dispnp());
    ddi_ = Teuchos::rcp(new Epetra_Vector(idispnp->Map(), true));
  }

  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *Extractor(), *Extractor(), 81, false, true));

  // initialize coupling matrices
  k_fs_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *(StructureField()->Interface()), *(FluidField()->Interface()), 81, false, true));

  k_sf_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *(FluidField()->Interface()), *(StructureField()->Interface()), 81, false, true));
}

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC            vuong 11/12  |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplit::BuildCombinedDBCMap()
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicSplit::CombinedDBCMap");

  // first, get DBC-maps from structure and fluid field and merge them to one map
  const Teuchos::RCP<const Epetra_Map> scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap = FluidField()->GetDBCMapExtractor()->CondMap();

  std::vector<Teuchos::RCP<const Epetra_Map>> vectoroverallfsimaps;
  vectoroverallfsimaps.push_back(scondmap);
  vectoroverallfsimaps.push_back(fcondmap);

  Teuchos::RCP<Epetra_Map> overallfsidbcmaps =
      LINALG::MultiMapExtractor::MergeMaps(vectoroverallfsimaps);

  // now we intersect the global dof map with the DBC map to get all dofs with DBS applied, which
  // are in the global
  // system, i.e. are not condensed
  std::vector<Teuchos::RCP<const Epetra_Map>> vectordbcmaps;
  vectordbcmaps.push_back(overallfsidbcmaps);
  vectordbcmaps.push_back(fullmap_);

  combinedDBCMap_ = LINALG::MultiMapExtractor::IntersectMaps(vectordbcmaps);

  return;
}

/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration            vuong 11/12    |
 *----------------------------------------------------------------------*/
void POROELAST::MonolithicSplit::Solve()
{
  // solve monolithic system by newton iteration
  Monolithic::Solve();

  // recover Lagrange multiplier \lambda_\Gamma at the interface at the end of each time step
  // (i.e. condensed forces onto the structure) needed for rhs in next time step
  RecoverLagrangeMultiplierAfterTimeStep();
}
