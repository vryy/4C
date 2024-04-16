/*----------------------------------------------------------------------*/
/*! \file

\level 1


\brief Entry routines for FSI problems and some other problem types as well
*/

/*----------------------------------------------------------------------*/

#include "baci_fsi_dyn.hpp"

#include "baci_adapter_ale_fsi.hpp"
#include "baci_adapter_fld_fluid_fsi.hpp"
#include "baci_adapter_fld_fluid_xfsi.hpp"
#include "baci_adapter_fld_moving_boundary.hpp"
#include "baci_adapter_str_fpsiwrapper.hpp"
#include "baci_adapter_str_fsiwrapper.hpp"
#include "baci_adapter_str_poro_wrapper.hpp"
#include "baci_adapter_str_structure.hpp"
#include "baci_ale_utils_clonestrategy.hpp"
#include "baci_binstrategy.hpp"
#include "baci_coupling_adapter.hpp"
#include "baci_coupling_adapter_mortar.hpp"
#include "baci_fluid_xfluid.hpp"
#include "baci_fluid_xfluid_fluid.hpp"
#include "baci_fsi_constrmonolithic_fluidsplit.hpp"
#include "baci_fsi_constrmonolithic_structuresplit.hpp"
#include "baci_fsi_dirichletneumann_disp.hpp"
#include "baci_fsi_dirichletneumann_factory.hpp"
#include "baci_fsi_dirichletneumann_vel.hpp"
#include "baci_fsi_dirichletneumann_volcoupl.hpp"
#include "baci_fsi_dirichletneumannslideale.hpp"
#include "baci_fsi_fluid_ale.hpp"
#include "baci_fsi_fluidfluidmonolithic_fluidsplit.hpp"
#include "baci_fsi_fluidfluidmonolithic_fluidsplit_nonox.hpp"
#include "baci_fsi_fluidfluidmonolithic_structuresplit.hpp"
#include "baci_fsi_fluidfluidmonolithic_structuresplit_nonox.hpp"
#include "baci_fsi_free_surface_monolithic.hpp"
#include "baci_fsi_lungmonolithic.hpp"
#include "baci_fsi_lungmonolithic_fluidsplit.hpp"
#include "baci_fsi_lungmonolithic_structuresplit.hpp"
#include "baci_fsi_monolithicfluidsplit.hpp"
#include "baci_fsi_monolithicstructuresplit.hpp"
#include "baci_fsi_mortarmonolithic_fluidsplit.hpp"
#include "baci_fsi_mortarmonolithic_fluidsplit_sp.hpp"
#include "baci_fsi_mortarmonolithic_structuresplit.hpp"
#include "baci_fsi_resulttest.hpp"
#include "baci_fsi_slidingmonolithic_fluidsplit.hpp"
#include "baci_fsi_slidingmonolithic_structuresplit.hpp"
#include "baci_fsi_utils.hpp"
#include "baci_fsi_xfem_fluid.hpp"
#include "baci_fsi_xfem_monolithic.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_fbi.hpp"
#include "baci_inpar_fsi.hpp"
#include "baci_lib_condition_selector.hpp"
#include "baci_lib_condition_utils.hpp"
#include "baci_lib_discret_xfem.hpp"
#include "baci_lib_dofset_fixed_size.hpp"
#include "baci_lib_resulttest.hpp"
#include "baci_lib_utils_createdis.hpp"
#include "baci_lib_utils_parallel.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_poroelast_utils_clonestrategy.hpp"
#include "baci_poroelast_utils_setup.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <functional>
#include <set>
#include <string>
#include <vector>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
// entry point for Fluid on Ale in DRT
/*----------------------------------------------------------------------*/
void fluid_ale_drt()
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  const Epetra_Comm& comm = problem->GetDis("fluid")->Comm();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  // check for xfem discretization
  if (CORE::UTILS::IntegralValue<bool>(
          (problem->XFluidDynamicParams().sublist("GENERAL")), "XFLUIDFLUID"))
  {
    FLD::XFluid::SetupFluidDiscretization();
  }
  else
  {
    fluiddis->FillComplete();
  }

  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
  aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes() == 0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);
    aledis->FillComplete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else  // filled ale discretization
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
      dserror(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));
  const int restart = problem->Restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    fluid->ReadRestart(restart);
  }
  fluid->Timeloop();

  GLOBAL::Problem::Instance()->AddFieldTest(fluid->MBFluidField()->CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(comm);
}

/*----------------------------------------------------------------------*/
// entry point for Fluid on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_xfem_drt()
{
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("structure")->Comm();

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  soliddis->FillComplete();

  FLD::XFluid::SetupFluidDiscretization();

  const Teuchos::ParameterList xfluid = problem->XFluidDynamicParams();
  bool alefluid = CORE::UTILS::IntegralValue<bool>((xfluid.sublist("GENERAL")), "ALE_XFluid");

  if (alefluid)  // in ale case
  {
    Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
    aledis->FillComplete();

    // create ale elements if the ale discretization is empty
    if (aledis->NumGlobalNodes() == 0)
    {
      DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(
          problem->GetDis("fluid"), aledis);
      aledis->FillComplete();
      // setup material in every ALE element
      Teuchos::ParameterList params;
      params.set<std::string>("action", "setup_material");
      aledis->Evaluate(params);
    }
    else  // filled ale discretization
    {
      if (!FSI::UTILS::FluidAleNodesDisjoint(problem->GetDis("fluid"), aledis))
        dserror(
            "Fluid and ALE nodes have the same node numbers. "
            "This it not allowed since it causes problems with Dirichlet BCs. "
            "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
    }
  }

  if (alefluid)
  {
    // create instance of fluid xfem algorithm, for moving interfaces
    Teuchos::RCP<FSI::FluidXFEMAlgorithm> fluidalgo =
        Teuchos::rcp(new FSI::FluidXFEMAlgorithm(comm));

    const int restart = GLOBAL::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fluidalgo->ReadRestart(restart);
    }

    // run the simulation
    fluidalgo->Timeloop();

    // perform result tests if required
    problem->AddFieldTest(fluidalgo->MBFluidField()->CreateFieldTest());
    problem->TestAll(comm);
  }
  else
  {
    //--------------------------------------------------------------
    // create instance of fluid basis algorithm
    const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();

    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo =
        Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

    //--------------------------------------------------------------
    // restart the simulation
    const int restart = GLOBAL::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fluidalgo->FluidField()->ReadRestart(restart);
    }

    //--------------------------------------------------------------
    // run the simulation
    fluidalgo->FluidField()->Integrate();

    //--------------------------------------------------------------
    // perform result tests if required
    problem->AddFieldTest(fluidalgo->FluidField()->CreateFieldTest());
    problem->TestAll(comm);
  }
}

/*----------------------------------------------------------------------*/
// entry point for (pure) free surface in DRT
/*----------------------------------------------------------------------*/
void fluid_freesurf_drt()
{
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("fluid")->Comm();

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->GetDis("fluid")->FillComplete();
  problem->GetDis("ale")->FillComplete();

  // get discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes() == 0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);
    aledis->FillComplete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else  // filled ale discretization
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
      dserror(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  int coupling = CORE::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
    {
      Teuchos::RCP<FSI::MonolithicMainFS> fsi;

      // Monolithic Free Surface Algorithm

      fsi = Teuchos::rcp(new FSI::MonolithicFS(comm, fsidyn));

      const int restart = GLOBAL::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        fsi->ReadRestart(restart);
      }

      fsi->Timeloop(fsi);

      GLOBAL::Problem::Instance()->AddFieldTest(fsi->FluidField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->TestAll(comm);
      break;
    }
    default:
    {
      Teuchos::RCP<FSI::FluidAleAlgorithm> fluid;

      // Partitioned FS Algorithm
      fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));

      fluid->Timeloop();

      GLOBAL::Problem::Instance()->AddFieldTest(fluid->MBFluidField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->TestAll(comm);
      break;
    }
  }
}
/*----------------------------------------------------------------------*/
// entry point for FSI using multidimensional immersed method (FBI)
/*----------------------------------------------------------------------*/
void fsi_immersed_drt()
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  const Epetra_Comm& comm = structdis->Comm();

  // Redistribute beams in the case of point coupling conditions

  if (structdis->GetCondition("PointCoupling") != nullptr)
  {
    structdis->FillComplete(false, false, false);
    DRT::UTILS::RedistributeDiscretizationsByBinning({structdis}, true);
  }
  else if (not structdis->Filled() || not structdis->HaveDofs())
  {
    structdis->FillComplete();
  }

  problem->GetDis("fluid")->FillComplete();

  // get discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");

  // create vector of discr.
  std::vector<Teuchos::RCP<DRT::Discretization>> dis;
  dis.push_back(fluiddis);
  dis.push_back(structdis);

  // binning strategy is created
  Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
      Teuchos::rcp(new BINSTRATEGY::BinningStrategy());

  const Teuchos::ParameterList& fbidyn = problem->FBIParams();

  INPAR::FBI::BeamToFluidPreSortStrategy presort_strategy =
      Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidPreSortStrategy>(fbidyn, "PRESORT_STRATEGY");

  // redistribute discr. with help of binning strategy
  if (presort_strategy == INPAR::FBI::BeamToFluidPreSortStrategy::binning)
  {
    std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
    std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;
    binningstrategy->Init(dis);
    Teuchos::RCP<Epetra_Map> rowbins =
        binningstrategy->DoWeightedPartitioningOfBinsAndExtendGhostingOfDiscretToOneBinLayer(
            dis, stdelecolmap, stdnodecolmap);
    binningstrategy->FillBinsIntoBinDiscretization(rowbins);
    binningstrategy->FillBinsIntoBinDiscretization(rowbins);
  }

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  // Any partitioned algorithm.
  Teuchos::RCP<FSI::Partitioned> fsi;

  INPAR::FSI::PartitionedCouplingMethod method =
      CORE::UTILS::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(
          fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");
  if (method == INPAR::FSI::DirichletNeumann)
  {
    fsi = FSI::DirichletNeumannFactory::CreateAlgorithm(comm, fsidyn);
    Teuchos::rcp_dynamic_cast<FSI::DirichletNeumann>(fsi, true)->Setup();
  }
  else
    dserror("unsupported partitioned FSI scheme");

  if (presort_strategy == INPAR::FBI::BeamToFluidPreSortStrategy::binning)
  {
    Teuchos::rcp_dynamic_cast<FSI::DirichletNeumannVel>(fsi, true)->SetBinning(binningstrategy);
  }

  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    fsi->ReadRestart(restart);
  }

  fsi->Timeloop(fsi);

  // create result tests for single fields
  GLOBAL::Problem::Instance()->AddFieldTest(fsi->MBFluidField()->CreateFieldTest());
  GLOBAL::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

  // do the actual testing
  GLOBAL::Problem::Instance()->TestAll(comm);
  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}
/*----------------------------------------------------------------------*/
// entry point for FSI using ALE in DRT
/*----------------------------------------------------------------------*/
void fsi_ale_drt()
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  const Epetra_Comm& comm = structdis->Comm();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  if (structdis->GetCondition("PointCoupling") != nullptr)
  {
    structdis->FillComplete(false, false, false);
    DRT::UTILS::RedistributeDiscretizationsByBinning({structdis}, true);
  }
  else if (not structdis->Filled() || not structdis->HaveDofs())
  {
    structdis->FillComplete();
  }

  if (CORE::UTILS::IntegralValue<bool>(
          (problem->XFluidDynamicParams().sublist("GENERAL")), "XFLUIDFLUID"))
  {
    FLD::XFluid::SetupFluidDiscretization();
  }
  else
    problem->GetDis("fluid")->FillComplete();

  problem->GetDis("ale")->FillComplete();

  // get discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes() == 0)  // empty ale discretization
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);
    aledis->FillComplete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else  // filled ale discretization (i.e. read from input file)
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
      dserror(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");

    if ((not CORE::UTILS::IntegralValue<bool>(problem->FSIDynamicParams(), "MATCHGRID_FLUIDALE")) or
        (not CORE::UTILS::IntegralValue<bool>(problem->FSIDynamicParams(), "MATCHGRID_STRUCTALE")))
    {
      // create vector of discr.
      std::vector<Teuchos::RCP<DRT::Discretization>> dis;
      dis.push_back(structdis);
      dis.push_back(fluiddis);
      dis.push_back(aledis);

      std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
      std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;

      // redistribute discr. with help of binning strategy
      if (structdis->Comm().NumProc() > 1)
      {
        // binning strategy is created and parallel redistribution is performed
        Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
            Teuchos::rcp(new BINSTRATEGY::BinningStrategy());
        binningstrategy->Init(dis);
        binningstrategy->DoWeightedPartitioningOfBinsAndExtendGhostingOfDiscretToOneBinLayer(
            dis, stdelecolmap, stdnodecolmap);
      }
    }
  }

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  FSI_COUPLING coupling = CORE::UTILS::IntegralValue<FSI_COUPLING>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
    case fsi_iter_lung_monolithicstructuresplit:
    case fsi_iter_lung_monolithicfluidsplit:
    case fsi_iter_constr_monolithicfluidsplit:
    case fsi_iter_constr_monolithicstructuresplit:
    case fsi_iter_mortar_monolithicstructuresplit:
    case fsi_iter_mortar_monolithicfluidsplit:
    case fsi_iter_mortar_monolithicfluidsplit_saddlepoint:
    case fsi_iter_fluidfluid_monolithicfluidsplit:
    case fsi_iter_fluidfluid_monolithicstructuresplit:
    case fsi_iter_sliding_monolithicfluidsplit:
    case fsi_iter_sliding_monolithicstructuresplit:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

      Teuchos::RCP<FSI::Monolithic> fsi;

      INPAR::FSI::LinearBlockSolver linearsolverstrategy =
          CORE::UTILS::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      // call constructor to initialize the base class
      if (coupling == fsi_iter_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_lung_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::LungMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_lung_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::LungMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_constr_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::ConstrMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_constr_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::ConstrMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_mortar_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::MortarMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_mortar_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::MortarMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_mortar_monolithicfluidsplit_saddlepoint)
      {
        fsi = Teuchos::rcp(new FSI::MortarMonolithicFluidSplitSaddlePoint(comm, fsidyn));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_sliding_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::SlidingMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::SlidingMonolithicStructureSplit(comm, fsidyn));
      }
      else
      {
        dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",
            coupling, linearsolverstrategy);
      }

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a Redistribute call
      const int restart = GLOBAL::Problem::Instance()->Restart();
      if (restart)
      {
        fsi->ReadRestart(restart);
      }

      // now do the coupling setup and create the combined dofmap
      fsi->SetupSystem();

      // possibly redistribute domain decomposition
      {
        const INPAR::FSI::Redistribute redistribute =
            CORE::UTILS::IntegralValue<INPAR::FSI::Redistribute>(fsimono, "REDISTRIBUTE");

        const double weight1 = fsimono.get<double>("REDIST_WEIGHT1");
        const double weight2 = fsimono.get<double>("REDIST_WEIGHT2");
        if (redistribute == INPAR::FSI::Redistribute_structure or
            redistribute == INPAR::FSI::Redistribute_fluid)
        {
          // redistribute either structure or fluid domain
          fsi->RedistributeDomainDecomposition(redistribute, coupling, weight1, weight2, comm, 0);

          // do setup again after redistribution
          fsi->SetupSystem();

          const double secweight1 =
              fsidyn.sublist("MONOLITHIC SOLVER").get<double>("REDIST_SECONDWEIGHT1");
          const double secweight2 =
              fsidyn.sublist("MONOLITHIC SOLVER").get<double>("REDIST_SECONDWEIGHT2");
          if (secweight1 != -1.0)
          {
            // redistribute either structure or fluid domain
            fsi->RedistributeDomainDecomposition(
                redistribute, coupling, secweight1, secweight2, comm, 0);

            // do setup again after redistribution
            fsi->SetupSystem();
          }
        }
        else if (redistribute == INPAR::FSI::Redistribute_both)
        {
          int numproc = comm.NumProc();

          // redistribute both structure and fluid domain
          fsi->RedistributeDomainDecomposition(
              INPAR::FSI::Redistribute_structure, coupling, weight1, weight2, comm, numproc / 2);

          // do setup again after redistribution (do this again here in between because the P matrix
          // changed!)
          fsi->SetupSystem();
          fsi->RedistributeDomainDecomposition(
              INPAR::FSI::Redistribute_fluid, coupling, weight1, weight2, comm, numproc / 2);

          // do setup again after redistribution
          fsi->SetupSystem();
        }
        else if (redistribute == INPAR::FSI::Redistribute_monolithic)
        {
          fsi->RedistributeMonolithicGraph(coupling, comm);

          // do setup again after redistribution
          fsi->SetupSystem();
        }
      }

      // here we go...
      fsi->Timeloop(fsi);

      // calculate errors in comparison to analytical solution
      fsi->FluidField()->CalculateError();

      // create result tests for single fields
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->AleField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->FluidField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi, fsidyn));
      GLOBAL::Problem::Instance()->AddFieldTest(fsitest);

      // do the actual testing
      GLOBAL::Problem::Instance()->TestAll(comm);

      break;
    }
    case fsi_iter_fluidfluid_monolithicfluidsplit_nonox:
    case fsi_iter_fluidfluid_monolithicstructuresplit_nonox:
    {
      Teuchos::RCP<FSI::MonolithicNoNOX> fsi;
      if (coupling == fsi_iter_fluidfluid_monolithicfluidsplit_nonox)
      {
        fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicFluidSplitNoNOX(comm, fsidyn));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nonox)
      {
        fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicStructureSplitNoNOX(comm, fsidyn));
      }
      else
        dserror("Unsupported monolithic XFFSI scheme");

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a Redistribute call
      const int restart = GLOBAL::Problem::Instance()->Restart();
      if (restart)
      {
        fsi->ReadRestart(restart);
      }

      // now do the coupling setup and create the combined dofmap
      fsi->SetupSystem();

      // here we go...
      fsi->Timeloop();

      // calculate errors in comparison to analytical solution
      fsi->FluidField()->CalculateError();

      // create result tests for single fields
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->FluidField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi, fsidyn));
      GLOBAL::Problem::Instance()->AddFieldTest(fsitest);

      // do the actual testing
      GLOBAL::Problem::Instance()->TestAll(comm);

      break;
    }
    default:
    {
      // Any partitioned algorithm.

      Teuchos::RCP<FSI::Partitioned> fsi;

      INPAR::FSI::PartitionedCouplingMethod method =
          CORE::UTILS::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(
              fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");

      switch (method)
      {
        case INPAR::FSI::DirichletNeumann:
        case INPAR::FSI::DirichletNeumannSlideale:
        case INPAR::FSI::DirichletNeumannVolCoupl:
          fsi = FSI::DirichletNeumannFactory::CreateAlgorithm(comm, fsidyn);
          Teuchos::rcp_dynamic_cast<FSI::DirichletNeumann>(fsi, true)->Setup();
          break;
        default:
          dserror("unsupported partitioned FSI scheme");
          break;
      }
      const int restart = GLOBAL::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        fsi->ReadRestart(restart);
      }

      fsi->Timeloop(fsi);

      // create result tests for single fields
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->MBFluidField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

      // do the actual testing
      GLOBAL::Problem::Instance()->TestAll(comm);

      break;
    }
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}

/*----------------------------------------------------------------------*/
// entry point for FSI using XFEM in DRT (also for ale case)
/*----------------------------------------------------------------------*/
void xfsi_drt()
{
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("structure")->Comm();

  if (comm.MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "       @..@    " << std::endl;
    std::cout << "      (----)      " << std::endl;
    std::cout << "     ( >__< )   " << std::endl;
    std::cout << "     ^^ ~~ ^^  " << std::endl;
    std::cout << "     _     _ _______ _______ _____" << std::endl;
    std::cout << "      \\\\__/  |______ |______   |  " << std::endl;
    std::cout << "     _/  \\\\_ |       ______| __|__" << std::endl;
    std::cout << std::endl << std::endl;
  }

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  Teuchos::RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  soliddis->FillComplete();

  FLD::XFluid::SetupFluidDiscretization();

  Teuchos::RCP<DRT::Discretization> fluiddis = GLOBAL::Problem::Instance()->GetDis(
      "fluid");  // at the moment, 'fluid'-discretization is used for ale!!!

  // CREATE ALE
  const Teuchos::ParameterList& xfdyn = problem->XFluidDynamicParams();
  bool ale = CORE::UTILS::IntegralValue<bool>((xfdyn.sublist("GENERAL")), "ALE_XFluid");
  Teuchos::RCP<DRT::Discretization> aledis;
  if (ale)
  {
    aledis = problem->GetDis("ale");
    if (aledis == Teuchos::null) dserror("XFSI DYNAMIC: ALE Discretization empty!!!");

    aledis->FillComplete(true, true, true);

    // Create ALE elements if the ale discretization is empty
    if (aledis->NumGlobalNodes() == 0)  // ALE discretization still empty
    {
      DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);
      aledis->FillComplete();
      // setup material in every ALE element
      Teuchos::ParameterList params;
      params.set<std::string>("action", "setup_material");
      aledis->Evaluate(params);
    }
    else  // ALE discretization already filled
    {
      if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
        dserror(
            "Fluid and ALE nodes have the same node numbers. "
            "This it not allowed since it causes problems with Dirichlet BCs. "
            "Use the ALE cloning functionality or ensure non-overlapping node numbering!");
    }
  }

  int coupling = CORE::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_xfem_monolithic:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

      INPAR::FSI::LinearBlockSolver linearsolverstrategy =
          CORE::UTILS::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      if (linearsolverstrategy != INPAR::FSI::PreconditionedKrylov)
        dserror("Only Newton-Krylov scheme with XFEM fluid");

      // create the MonolithicXFEM object that does the whole work
      Teuchos::RCP<FSI::AlgorithmXFEM> fsi = Teuchos::rcp(new FSI::MonolithicXFEM(comm, fsidyn));

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a Redistribute call
      const int restart = GLOBAL::Problem::Instance()->Restart();
      if (restart)
      {
        fsi->ReadRestart(restart);
      }

      // setup the system (block-DOF-row maps, systemmatrix etc.) for the monolithic XFEM system
      fsi->SetupSystem();

      // here we go...
      fsi->Timeloop();

      GLOBAL::Problem::Instance()->AddFieldTest(fsi->FluidField()->CreateFieldTest());
      fsi->StructurePoro()->TestResults(GLOBAL::Problem::Instance());

      //    // create FSI specific result test
      //    Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new
      //    FSI::FSIResultTest(fsi,fsidyn)); GLOBAL::Problem::Instance()->AddFieldTest(fsitest);

      // do the actual testing
      GLOBAL::Problem::Instance()->TestAll(comm);

      break;
    }
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
      dserror("Unreasonable choice");
      break;
    default:
    {
      // Any partitioned algorithm. Stable of working horses.

      Teuchos::RCP<FSI::Partitioned> fsi;

      INPAR::FSI::PartitionedCouplingMethod method =
          CORE::UTILS::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(
              fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");

      switch (method)
      {
        case INPAR::FSI::DirichletNeumann:
          fsi = FSI::DirichletNeumannFactory::CreateAlgorithm(comm, fsidyn);
          Teuchos::rcp_dynamic_cast<FSI::DirichletNeumann>(fsi, true)->Setup();
          break;
        default:
          dserror("only Dirichlet-Neumann partitioned schemes with XFEM");
          break;
      }

      const int restart = GLOBAL::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        fsi->ReadRestart(restart);
      }

      fsi->Timeloop(fsi);

      GLOBAL::Problem::Instance()->AddFieldTest(fsi->MBFluidField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
      GLOBAL::Problem::Instance()->TestAll(comm);

      break;
    }
  }

  Teuchos::TimeMonitor::summarize();
}

/*----------------------------------------------------------------------*/
// entry point for FPSI using XFEM in DRT
/*----------------------------------------------------------------------*/
void xfpsi_drt()
{
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("structure")->Comm();

  if (comm.MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "       @..@    " << std::endl;
    std::cout << "      (----)      " << std::endl;
    std::cout << "     ( >__< )   " << std::endl;
    std::cout << "     ^^ ~~ ^^  " << std::endl;
    std::cout << "     _     _ _______  ______  _______ _____" << std::endl;
    std::cout << "      \\\\__/  |______ ||____|  |______   |  " << std::endl;
    std::cout << "     _/  \\\\_ |       ||       ______| __|__" << std::endl;
    std::cout << std::endl << std::endl;
  }
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // 1.-Initialization.
  // setup of the discretizations, including clone strategy
  POROELAST::UTILS::SetupPoro<POROELAST::UTILS::PoroelastCloneStrategy>();

  // setup of discretization for xfluid
  FLD::XFluid::SetupFluidDiscretization();
  Teuchos::RCP<DRT::Discretization> fluiddis = GLOBAL::Problem::Instance()->GetDis(
      "fluid");  // at the moment, 'fluid'-discretization is used for ale!!!

  Teuchos::RCP<DRT::Discretization> aledis;
  const Teuchos::ParameterList& xfdyn = problem->XFluidDynamicParams();
  bool ale = CORE::UTILS::IntegralValue<bool>((xfdyn.sublist("GENERAL")), "ALE_XFluid");
  if (ale)
  {
    aledis = problem->GetDis("ale");
    if (aledis == Teuchos::null) dserror("Ale Discretization empty!");

    aledis->FillComplete(true, true, true);

    // 3.- Create ALE elements if the ale discretization is empty
    if (aledis->NumGlobalNodes() == 0)  // ALE discretization still empty
    {
      DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);
      aledis->FillComplete();
      // setup material in every ALE element
      Teuchos::ParameterList params;
      params.set<std::string>("action", "setup_material");
      aledis->Evaluate(params);
    }
    else  // ALE discretization already filled
    {
      if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
        dserror(
            "Fluid and ALE nodes have the same node numbers. "
            "This it not allowed since it causes problems with Dirichlet BCs. "
            "Use the ALE cloning functionality or ensure non-overlapping node numbering!");
    }
  }

  // print all dofsets
  fluiddis->GetDofSetProxy()->PrintAllDofsets(fluiddis->Comm());

  // 2.- Parameter reading
  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  int coupling = CORE::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");

  switch (coupling)
  {
    case fsi_iter_xfem_monolithic:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
      INPAR::FSI::LinearBlockSolver linearsolverstrategy =
          CORE::UTILS::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      if (linearsolverstrategy != INPAR::FSI::PreconditionedKrylov)
        dserror("Only Newton-Krylov scheme with XFEM fluid");

      Teuchos::RCP<FSI::AlgorithmXFEM> fsi = Teuchos::rcp(
          new FSI::MonolithicXFEM(comm, fsidyn, ADAPTER::FieldWrapper::type_PoroField));

      // read the restart information, set vectors and variables ---

      // be careful, dofmaps might be changed here in a Redistribute call
      const int restart =
          GLOBAL::Problem::Instance()
              ->Restart();  // not adapated at the moment .... Todo check it .. ChrAg
      if (restart)
      {
        fsi->ReadRestart(restart);
      }

      fsi->SetupSystem();


      // 3.2.- redistribute the FPSI interface
      // Todo .... fsi->RedistributeInterface(); // this is required for paralles fpi-condition (not
      // included in this commit)

      // here we go...
      fsi->Timeloop();

      GLOBAL::Problem::Instance()->AddFieldTest(fsi->FluidField()->CreateFieldTest());
      fsi->StructurePoro()->TestResults(GLOBAL::Problem::Instance());

      // do the actual testing
      GLOBAL::Problem::Instance()->TestAll(comm);
      break;
    }
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
      dserror("Unreasonable choice");
      break;
    default:
    {
      dserror("FPSI_XFEM: No Partitioned Algorithms implemented !!!");
      break;
    }
  }
  Teuchos::TimeMonitor::summarize();
}

BACI_NAMESPACE_CLOSE
