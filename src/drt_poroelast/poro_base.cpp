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
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/inpar_solver.H"

#include <Teuchos_TimeMonitor.hpp>
// needed for PrintNewton
#include <sstream>

// include this header for coupling stiffness terms
#include "../drt_lib/drt_assemblestrategy.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_structure/stru_aux.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.

/*----------------------------------------------------------------------*
 | constructor (public)                                    vuong 01/12  |
 *----------------------------------------------------------------------*/
POROELAST::PoroBase::PoroBase(const Epetra_Comm& comm,
                                          const Teuchos::ParameterList& timeparams) :
      AlgorithmBase(comm, timeparams)
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // ask base algorithm for the structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams, structdis));
  structure_ = rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams,true));
  fluid_ = fluid->FluidFieldrcp();
  StructureField()->Discretization();

  // monolithic Poroelasticity must know the other discretization
  // build a proxy of the structure discretization for the fluid field
  Teuchos::RCP<DRT::DofSet> structdofset =
      StructureField()->Discretization()->GetDofSetProxy();
  // build a proxy of the fluid discretization for the structure field
  Teuchos::RCP<DRT::DofSet> fluiddofset =
      FluidField().Discretization()->GetDofSetProxy();

  // check if FluidField has 2 discretizations, so that coupling is possible
  if (FluidField().Discretization()->AddDofSet(structdofset) != 1)
    dserror("unexpected dof sets in fluid field");
  if (StructureField()->Discretization()->AddDofSet(fluiddofset)!=1)
    dserror("unexpected dof sets in structure field");

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
  dserror("monolithic Poroelasticity is limited in functionality (only one-step-theta scheme possible)");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* structurenodemap = StructureField()->Discretization()->NodeRowMap();

  coupfs_ = Teuchos::rcp(new ADAPTER::Coupling());
  const int ndim = DRT::Problem::Instance()->NDim();
  coupfs_->SetupCoupling(*FluidField().Discretization(),
                         *StructureField()->Discretization(),
                         *fluidnodemap,
                         *structurenodemap,
                          ndim);

  FluidField().SetMeshMap(coupfs_->MasterDofMap());

  //extractor for constraints on structure phase
  //
  // when using constraints applied via Lagrange-Multipliers there is a
  // difference between StructureField()->DofRowMap() and StructureField()->DofRowMap(0).
  // StructureField()->DofRowMap(0) returns the DofRowMap
  // known to the discretization (without lagrange multipliers)
  // while StructureField()->DofRowMap() returns the DofRowMap known to
  // the constraint manager (with lagrange multipliers)
  consplitter_=LINALG::MapExtractor(*StructureField()->DofRowMap(),
      StructureField()->DofRowMap(0));

  FluidField().Discretization()->GetCondition("NoPenetration", nopencond_);
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
      dispn = consplitter_.ExtractCondVector(StructureField()->Dispnp());
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
 | prepare time step (public)                         vuong 01/12       |
 *----------------------------------------------------------------------*/
void POROELAST::PoroBase::PrepareTimeStep()
{
  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

  // call the predictor
  StructureField()->PrepareTimeStep();
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
  return coupfs_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::PoroBase::FluidToStructureField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfs_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::PoroBase::StructureToFluidAtInterface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfs_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> POROELAST::PoroBase::FluidToStructureAtInterface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfs_->MasterToSlave(iv);
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
  Teuchos::RCP<Epetra_Map> nopendofmap = rcp(new Epetra_Map(-1, condIDs.size(), &condIDs[0], 0, FluidField().Discretization()->Comm()));

  nopenetration_ = LINALG::MapExtractor(*FluidField().DofRowMap(), nopendofmap);

  return;
}
