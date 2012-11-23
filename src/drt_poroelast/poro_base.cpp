/*----------------------------------------------------------------------*/
/*!
 \file poro_base.cpp

 \brief  Basis of all porous media algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  vuong 01/12 |
 *----------------------------------------------------------------------*/
#include "poro_base.H"
#include "poroelast_defines.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"


#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_solver.H"

#include <Teuchos_TimeMonitor.hpp>
// needed for PrintNewton
#include <sstream>

// include this header for coupling stiffness terms
#include "../drt_lib/drt_assemblestrategy.H"

#include "../drt_lib/drt_utils.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_structure/stru_aux.H"

/*----------------------------------------------------------------------*
 | constructor (public)                                    vuong 01/12  |
 *----------------------------------------------------------------------*/
POROELAST::PoroBase::PoroBase(const Epetra_Comm& comm,
                                          const Teuchos::ParameterList& timeparams) :
      AlgorithmBase(comm, timeparams),
      nopencond_(0),
      condIDs_(Teuchos::null),
      subnodemap_(Teuchos::null),
      subelemap_(Teuchos::null)
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // ask base algorithm for the structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams, structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams,true));
  fluid_ = Teuchos::rcp_dynamic_cast<ADAPTER::FluidPoro>(fluid->FluidFieldrcp());

  // access the problem-specific parameter lists
  const Teuchos::ParameterList& sdyn
  = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& fdyn
  = DRT::Problem::Instance()->FluidDynamicParams();

  // check time integration algo -> currently only one-step-theta scheme supported
  INPAR::STR::DynamicType structtimealgo
  = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo
  = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  if ( structtimealgo != INPAR::STR::dyna_onesteptheta or
      fluidtimealgo != INPAR::FLUID::timeint_one_step_theta )
  dserror("porous media problem is limited in functionality (only one-step-theta scheme possible)");

  //check for submeshes and build of subnode and subelement map if necessary
  submeshes_ = BuildSubMaps();

  // the problem is two way coupled, thus each discretization must know the other discretization
  Teuchos::RCP<DRT::DofSet> structdofset = Teuchos::null;
  Teuchos::RCP<DRT::DofSet> fluiddofset = Teuchos::null;

  /* When coupling porous media with a pure structure or fluid we will have two discretizations
   * of different size. In this case we need a special proxy, which can handle submeshes.
   */
  if(submeshes_)
  {
    // build a proxy of the structure discretization for the fluid field
    structdofset = StructureField()->Discretization()->GetDofSetProxy(subnodemap_, subelemap_);
    // build a proxy of the fluid discretization for the structure field
    fluiddofset = FluidField().Discretization()->GetDofSetProxy(subnodemap_, subelemap_);
  }
  else
  {
    // build a proxy of the structure discretization for the fluid field
    structdofset = StructureField()->Discretization()->GetDofSetProxy();
    // build a proxy of the fluid discretization for the structure field
    fluiddofset = FluidField().Discretization()->GetDofSetProxy();
  }

  // check if FluidField has 2 discretizations, so that coupling is possible
  if (FluidField().Discretization()->AddDofSet(structdofset) != 1)
    dserror("unexpected dof sets in fluid field");
  if (StructureField()->Discretization()->AddDofSet(fluiddofset)!=1)
    dserror("unexpected dof sets in structure field");

  // the fluid-ale coupling not always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* structurenodemap = StructureField()->Discretization()->NodeRowMap();

  coupfs_ = Teuchos::rcp(new ADAPTER::Coupling());
  const int ndim = DRT::Problem::Instance()->NDim();

  coupfs_->SetupCoupling(*StructureField()->Discretization(),
                         *FluidField().Discretization(),
                         *structurenodemap,
                         *fluidnodemap,
                          ndim,
                          not submeshes_);

  if(submeshes_)
    psiextractor_ = Teuchos::rcp(new LINALG::MapExtractor(*StructureField()->DofRowMap(), coupfs_->MasterDofMap()));

  //FluidField().SetMeshMap(coupfs_->MasterDofMap());
  FluidField().SetMeshMap(coupfs_->SlaveDofMap());

  //extractor for constraints on structure phase
  //
  // when using constraints applied via Lagrange-Multipliers there is a
  // difference between StructureField()->DofRowMap() and StructureField()->DofRowMap(0).
  // StructureField()->DofRowMap(0) returns the DofRowMap
  // known to the discretization (without lagrange multipliers)
  // while StructureField()->DofRowMap() returns the DofRowMap known to
  // the constraint manager (with lagrange multipliers)
  consplitter_ = Teuchos::rcp(new LINALG::MapExtractor(*StructureField()->DofRowMap(),
      StructureField()->DofRowMap(0)));

  FluidField().Discretization()->GetCondition("NoPenetration", nopencond_);

  //do some checks
  {
    std::vector<DRT::Condition*> porocoupl;
    FluidField().Discretization()->GetCondition("PoroCoupling", porocoupl);
    if ( porocoupl.size() == 0 )
      dserror("no Poro Coupling Condition defined for porous media problem. Fix your input file!");
  }
}

/*----------------------------------------------------------------------*
 | destructor (public)                                    vuong 01/12   |
 *----------------------------------------------------------------------*/
POROELAST::PoroBase::~PoroBase()
{
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::ReadRestart( int restart)
{
  if (restart)
  {
    FluidField().ReadRestart(restart);
    StructureField()->ReadRestart(restart);

    // apply current velocity and pressures to structure
    StructureField()->ApplyVelAndPress(FluidField().Velnp());

    Teuchos::RCP<Epetra_Vector> dispn;
    if (StructureField()->HaveConstraint())
    {
      //displacment vector without lagrange-multipliers
      dispn = consplitter_->ExtractCondVector(StructureField()->Dispnp());
    }
    else
      dispn = StructureField()->ExtractDispnp();

    // transfer the current structure displacement to the fluid field
    Teuchos::RCP<Epetra_Vector> structdisp = StructureToFluidField(dispn);
    FluidField().ApplyMeshDisplacement(structdisp);

    // transfer the current structure velocity to the fluid field
    Teuchos::RCP<Epetra_Vector> structvel = StructureToFluidField(
        StructureField()->ExtractVelnp());
    FluidField().ApplyMeshVelocity(structvel);

    // second ReadRestart needed due to the coupling variables
    FluidField().ReadRestart(restart);
    StructureField()->ReadRestart(restart);

    SetTimeStep(FluidField().Time(), restart);
  }

  return;
}

/*----------------------------------------------------------------------*
 | prepare time step (protected)                         vuong 01/12       |
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::PrepareTimeStep()
{
// counter and print header
  IncrementTimeAndStep();
  PrintHeader();

// call the predictor
  StructureField()-> PrepareTimeStep();
  FluidField().PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | update (protected)                                     vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::Update()
{
  StructureField()->Update();
  FluidField().Update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::PrepareOutput()
{
  StructureField()->PrepareOutput();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(StructureField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(FluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::PoroBase::StructureToFluidField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  if(submeshes_)
    return coupfs_->MasterToSlave(psiextractor_->ExtractCondVector(iv));
  else
    return coupfs_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::PoroBase::FluidToStructureField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  if(submeshes_)
    return coupfs_->SlaveToMaster(psiextractor_->ExtractCondVector(iv));
  else
    return coupfs_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PoroBase::BuidNoPenetrationMap()
{
  std::vector<int> condIDs;
  std::set<int>::iterator it;
  for(it=condIDs_->begin();it!=condIDs_->end();it++)
  {
    condIDs.push_back(*it);
  }
  Teuchos::RCP<Epetra_Map> nopendofmap = Teuchos::rcp(new Epetra_Map(-1, condIDs.size(), &condIDs[0], 0, FluidField().Discretization()->Comm()));

  nopenetration_ = Teuchos::rcp(new LINALG::MapExtractor(*FluidField().DofRowMap(), nopendofmap));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PoroBase::SetStructSolution()
{
  Teuchos::RCP<Epetra_Vector> dispn;
    // apply current displacements and velocities to the fluid field
  if (StructureField()->HaveConstraint())
    //displacement vector without lagrange-multipliers
    dispn = consplitter_->ExtractCondVector(StructureField()->Dispnp());
  else
    dispn = StructureField()->ExtractDispnp();

  Teuchos::RCP<Epetra_Vector> veln = StructureField()->ExtractVelnp();

  // transfer the current structure displacement to the fluid field
  Teuchos::RCP<Epetra_Vector> structdisp = StructureToFluidField(dispn);
  FluidField().ApplyMeshDisplacement(structdisp);

  // transfer the current structure velocity to the fluid field
  Teuchos::RCP<Epetra_Vector> structvel = StructureToFluidField(veln);
  FluidField().ApplyMeshVelocity(structvel);

  CalculateSurfPoro("PoroPartInt");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PoroBase::SetFluidSolution()
{
  StructureField()->ApplyVelAndPress(FluidField().Velnp());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PoroBase::TimeLoop()
{

  while (NotFinished())
  {
    //solve one time step
    Solve();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  StructureField()->Output();
} // Monolithic::Output()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool POROELAST::PoroBase::BuildSubMaps()
{
  std::set<int> subnodegids;
  const Epetra_Map* structnodemap = StructureField()->Discretization()->NodeColMap();
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeColMap();
  const int numstructnodes = structnodemap->NumMyElements();
  const int numfluidnodes = fluidnodemap->NumMyElements();

  if(numstructnodes != numfluidnodes)
  {

    for(int i=0; i<(numfluidnodes);i++)
    {
      const int fluidgid = fluidnodemap->GID(i);
      if(fluidgid == -1) dserror("Node with gobal id gid=%d not stored on this proc",fluidgid);
      const int structlid = structnodemap->LID(fluidgid);
      if(structlid != -1)
        subnodegids.insert(fluidgid);
    }

    subnodemap_ = LINALG::CreateMap(subnodegids, Comm());

    std::set<int> subelegids;
    const Epetra_Map* structelemap = StructureField()->Discretization()->ElementColMap();
    const Epetra_Map* fluidelemap = FluidField().Discretization()->ElementColMap();

    for(int i=0; i<(fluidelemap->NumMyElements());i++)
    {
      const int fluidgid = fluidelemap->GID(i);
      if(fluidgid == -1) dserror("Element with gobal id gid=%d not stored on this proc",fluidgid);
      const int structlid = structelemap->LID(fluidgid);
      if(structlid != -1)
        subelegids.insert(fluidgid);
    }

    subelemap_ = LINALG::CreateMap(subelegids, Comm());

    return true;
  }
  else
    return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::CalculateSurfPoro(const string& condstring)
{
  //-------------------------------
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  p.set("action", "calc_struct_area_poro");
  // other parameters that might be needed by the elements
  p.set("total time", Time());
  p.set("delta time", Dt());

  const int numdim = DRT::Problem::Instance()->NDim();

  Teuchos::RCP<DRT::Discretization> structdis = StructureField()->Discretization();
  Teuchos::RCP<DRT::Discretization> fluiddis = FluidField().Discretization();

  // set vector values needed by elements
  structdis->ClearState();
  // extended SetState(0,...) in case of multiple dofsets (e.g. TSI)
  structdis->SetState(0,"displacement", StructureField()->Dispnp());

  structdis->SetState(1,"fluidvel",FluidField().Velnp());

  // velocities (always three velocity components per node)
  // (get noderowmap of discretization for creating this multivector)
  const Epetra_Map* noderowmap = structdis->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> convel = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim+1,true));
  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<structdis->NumMyRowNodes();lnodeid++)
  {
    // get the processor local node
    DRT::Node* lnode = structdis->lRowNode(lnodeid);

    vector<int> nodedofs = fluiddis->Dof(0,lnode);
    for(int index=0;index<numdim+1;++index)
    {
      // get global and local ID
      const int gid = nodedofs[index];
      // const int lid = dofrowmap->LID(gid);
      const int lid = FluidField().Velnp()->Map().LID(gid);
      if (lid < 0) dserror("Local ID not found in map for given global ID!");

      double convelocity = (*(FluidField().Velnp()))[lid];

      // insert velocity value into node-based vector
      int err = convel->ReplaceMyValue(lnodeid,index,convelocity);
      if (err != 0) dserror("Error while inserting value into vector convel_!");
    }
  }

 AddMultiVectorToParameterList(p,"convective velocity field",convel,fluiddis);

  StructureField()->Discretization()->EvaluateCondition(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condstring);
  StructureField()->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 | export multivector to column map & add it to parameter list  vuong 11/12|
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::AddMultiVectorToParameterList
(Teuchos::ParameterList& p,
 const std::string name,
 Teuchos::RCP<Epetra_MultiVector> vec,
 Teuchos::RCP<DRT::Discretization> discret
)
{
  if (vec != Teuchos::null)
  {
    //provide data in node-based multi-vector for usage on element level
    // -> export to column map is necessary for parallel evaluation
    //SetState cannot be used since this multi-vector is nodebased and not dofbased!
    const Epetra_Map* nodecolmap = discret->NodeColMap();
    int numcol = vec->NumVectors();
    RCP<Epetra_MultiVector> tmp = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,numcol));
    LINALG::Export(*vec,*tmp);
    p.set(name,tmp);
  }
  else
    p.set(name,Teuchos::null);

  return;
}
