/*!----------------------------------------------------------------------
\file fps3i_partitioned.cpp
\brief General algorithmic routines for partitioned solution approaches
       to fluid-porous-structure-scalar-scalar interaction (FPS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively.

 <pre>
   Maintainer: Moritz Thon & Andre Hemmler
               thon@mhpc.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>

*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_io/io_control.H"
#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_fpsi/fpsi_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_poroelast/poroelast_monolithic.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_poroelast/poro_utils_clonestrategy.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_inpar/inpar_scatra.H"

#include "fps3i_partitioned.H"


/*----------------------------------------------------------------------*
 |  Constructor                                           hemmler 07/14 |
 *----------------------------------------------------------------------*/
FS3I::PartFPS3I::PartFPS3I(const Epetra_Comm& comm)
  : FS3I_Base(),
    comm_(comm)
{
  if(comm.MyPID()==0)
  {
    //##################       0.- Warning          //#########################
    std::cout<<std::endl;
    std::cout<<"##############################################################################"<<std::endl;
    std::cout<<"################################# WARNING!!! #################################"<<std::endl;
    std::cout<<"##############################################################################"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"This version of Fluid-porous-structure-scatra-scatra interaction (FPS3I) does NOT"<<std::endl;
    std::cout<<"account for the convective scalar transport at the fluid-poro interface!"<<std::endl;
    std::cout<<"The conservation of mass at the interface is only guaranteed for purely diffusive transport"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"##############################################################################"<<std::endl;
    std::cout<<"################################# WARNING!!! #################################"<<std::endl;
    std::cout<<"##############################################################################"<<std::endl;
    std::cout<<std::endl;
  }
  //##################       1.- Parameter reading          //#########################
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fs3icontrol = problem->FS3IControlParams();
  const Teuchos::ParameterList& fpsidynparams       = problem->FPSIDynamicParams();
  const Teuchos::ParameterList& poroelastdynparams  = problem->PoroelastDynamicParams();

  double dt_fpsi = fpsidynparams.get<double>("TIMESTEP");
  double dt_poroelast = poroelastdynparams.get<double>("TIMESTEP");
  if(dt_fpsi!=dt_poroelast)
  {
    dserror("Please set \"TIMESTEP\" in \"POROELASTICITY DYNAMIC\" to the same value as in \"FPSI DYNAMIC\"!");
  }

  Teuchos::RCP<FPSI::UTILS> FPSI_UTILS = FPSI::UTILS::Instance();

  //##################    2.- Creation of Poroelastic + Fluid problem. (Discretization called inside)     //##################
  Teuchos::RCP<FPSI::FPSI_Base> fpsi_algo = Teuchos::null;

  fpsi_algo = FPSI_UTILS->SetupDiscretizations(comm, fpsidynparams,poroelastdynparams);

  //only monolithic coupling of fpsi problem is supported!
  int coupling = DRT::INPUT::IntegralValue<int>(fpsidynparams,"COUPALGO");
  if (coupling == fpsi_monolithic_plain)
  {
    //Cast needed because functions such as PoroField() and FluidField() are just a member-functions of the derived class Monolithic_Plain, but not of the base class FPSI_Base
    fpsi_= Teuchos::rcp_dynamic_cast<FPSI::Monolithic_Plain>(fpsi_algo);
  }
  else
  {
    dserror("Partitioned solution scheme not implemented for FPSI, yet. "
            "Make sure that the parameter COUPALGO is set to 'fpsi_monolithic_plain', "
            "and the parameter PARITIONED is set to 'monolithic'. ");
  }

  //##################      3. Discretization of Scatra problem       //##################
  problem->GetDis("scatra1")->FillComplete();
  problem->GetDis("scatra2")->FillComplete();

  //---------------------------------------------------------------------
  // access discretizations for poro (structure) and fluid as well as fluid-
  // and poro-based scalar transport and get material map for fluid
  // and scalar transport elements
  //---------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluidscatradis = problem->GetDis("scatra1");
  Teuchos::RCP<DRT::Discretization> structscatradis = problem->GetDis("scatra2");

  //---------------------------------------------------------------------
  // create discretization for fluid-based scalar transport from and
  // according to fluid discretization
  //---------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");


  //std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap = problem->CloningMaterialMap();
  //if (clonefieldmatmap.size() < 2)
  //  dserror("At least two material lists required for partitioned FS3I!");

  // create fluid-based scalar transport elements if fluid-based scalar
  // transport discretization is empty
  if (fluidscatradis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,fluidscatradis);
  }
  else
    dserror("Fluid AND ScaTra discretization present. This is not supported.");

  //---------------------------------------------------------------------
  // create discretization for poro-based scalar transport from and
  // according to poro (structure) discretization
  //--------------------------------------------------------------------

  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  if(!structscatradis->Filled())
    structscatradis->FillComplete();
  if (structscatradis->NumGlobalNodes()==0)
  {
    // create the poro scatra discretization
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroScatraCloneStrategy>(structdis,structscatradis);

    // redistribute FPSI interface here, since if done before the PoroScatra cloning does not work
    //fpsi_->RedistributeInterface();
    // after redistributing the interface we have to fix the material pointers of the structure-scatra discretisation
    //POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,structscatradis);
  }
  else
  dserror("Structure AND ScaTra discretization present. This is not supported.");

  //##################      End of discretization       //##################

  //---------------------------------------------------------------------
  // create instances for fluid- and poro (structure)-based scalar transport
  // solver and arrange them in combined vector
  //---------------------------------------------------------------------
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3icontrol.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3icontrol.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    dserror("no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");
  if (linsolver2number == (-1))
    dserror("no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(fs3icontrol,true,"scatra1",problem->SolverParams(linsolver1number)));
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(fs3icontrol,true,"scatra2",problem->SolverParams(linsolver2number)));

  scatravec_.push_back(fluidscatra);
  scatravec_.push_back(structscatra);

  //Manipulates the scatratype for the Poro-domain to scatratype_poro
  INPAR::SCATRA::ScaTraType new_ScaTraType = INPAR::SCATRA::scatratype_pororeac;
  scatravec_[1]->ScaTraField()->ManipulateScaTraType(new_ScaTraType);

  //---------------------------------------------------------------------
  // check various input parameters
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& structdyn = problem->StructuralDynamicParams();
  const Teuchos::ParameterList& fluiddyn  = problem->FluidDynamicParams();
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();
  // check consistency of time-integration schemes in input file
  // (including parameter theta itself in case of one-step-theta scheme)
  // and rule out unsupported versions of generalized-alpha time-integration
  // scheme (as well as other inappropriate schemes) for fluid subproblem
  INPAR::SCATRA::TimeIntegrationScheme scatratimealgo = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fluiddyn,"TIMEINTEGR");

  INPAR::STR::DynamicType structtimealgo = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdyn,"DYNAMICTYP");

  if (fluidtimealgo  == INPAR::FLUID::timeint_one_step_theta)
  {
    if (scatratimealgo != INPAR::SCATRA::timeint_one_step_theta or
        structtimealgo != INPAR::STR::dyna_onesteptheta)
      dserror("Partitioned FS3I computations should feature consistent time-integration schemes for the subproblems; in this case, a one-step-theta scheme is intended to be used for the fluid subproblem, and different schemes are intended to be used for the structure and/or scalar transport subproblems!");

    if (scatradyn.get<double>("THETA") != fluiddyn.get<double>("THETA") or
        scatradyn.get<double>("THETA") != structdyn.sublist("ONESTEPTHETA").get<double>("THETA"))
    dserror("Parameter(s) theta for one-step-theta time-integration scheme defined in one or more of the individual fields do(es) not match for partitioned FS3I computation.");
  }
  else if (fluidtimealgo  == INPAR::FLUID::timeint_afgenalpha)
  {
    if (scatratimealgo != INPAR::SCATRA::timeint_gen_alpha or
        structtimealgo != INPAR::STR::dyna_genalpha)
      dserror("Partitioned FS3I computations should feature consistent time-integration schemes for the subproblems; in this case, a (alpha_f-based) generalized-alpha scheme is intended to be used for the fluid subproblem, and different schemes are intended to be used for the structure and/or scalar transport subproblems!");
  }
  else if (fluidtimealgo  == INPAR::FLUID::timeint_npgenalpha)
  {
      dserror("Partitioned FS3I computations do not support n+1-based generalized-alpha time-integration schemes for the fluid subproblem!");
  }
  else if (fluidtimealgo  == INPAR::FLUID::timeint_bdf2 or
           fluidtimealgo  == INPAR::FLUID::timeint_stationary)
  {
      dserror("Partitioned FS3I computations do not support stationary of BDF2 time-integration schemes for the fluid subproblem!");
  }

  // check that incremental formulation is used for scalar transport field,
  // according to structure and fluid field
  if (scatravec_[0]->ScaTraField()->IsIncremental() == false)
    dserror("Incremental formulation required for partitioned FS3I computations!");

  //---------------------------------------------------------------------
  // check existence of scatra coupling conditions for both
  // discretizations and definition of the permeability coefficient
  //---------------------------------------------------------------------
  CheckFS3IInputs();

  //in case of FPS3I we have to handle the conductivity, too
  Teuchos::RCP<DRT::Discretization> dis = scatravec_[0]->ScaTraField()->Discretization();
  std::vector<DRT::Condition*> coupcond;
  dis->GetCondition("ScaTraCoupling",coupcond);
  double myconduct = coupcond[0]->GetDouble("hydraulic conductivity"); //here we assume the conductivity to be the same in every BC

  //conductivity is not only needed in scatracoupling but also in FPSI coupling!
  if(myconduct==0.0)
  {
    dserror("conductivity of 0.0 is not allowed!!! Should be set in \"DESIGN SCATRA COUPLING SURF CONDITIONS\"");
  }
  fpsi_->SetConductivity(myconduct);

}


/*----------------------------------------------------------------------*
 |  Restart                                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ReadRestart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  DRT::Problem* problem = DRT::Problem::Instance();
  const int restart = problem->Restart();

  if (restart)
  {
   //restart of FPSI problem
    fpsi_->ReadRestart(restart);

    //restart of scatra problem
    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
      currscatra->ScaTraField()->ReadRestart(restart);
    }

    time_ = fpsi_->FluidField()->Time();
    step_ = fpsi_->FluidField()->Step();
  }
}

/*----------------------------------------------------------------------*
 | redistribute the FPSI interface                           thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::RedistributeInterface()
{
  fpsi_->RedistributeInterface();

  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> structscatradis = problem->GetDis("scatra2");

  // after redistributing the interface we have to fix the material pointers of the structure-scatra discretisation
  POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,structscatradis);
}

/*----------------------------------------------------------------------*
 |  System Setup                                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetupSystem()
{
  //do the coupling setup and create the combined dofmap

  // Setup FPSI system
  fpsi_->SetupSystem();

  /*----------------------------------------------------------------------*/
  /*                  General set up for scalar fields                    */
  /*----------------------------------------------------------------------*/

  // create map extractors needed for scatra condition coupling
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
    const int numscal = currscatra->ScaTraField()->NumScal();
    Teuchos::RCP<DRT::Discretization> currdis = currscatra->ScaTraField()->Discretization();
    Teuchos::RCP<LINALG::MultiMapExtractor> mapex = Teuchos::rcp(new LINALG::MultiMapExtractor());
    DRT::UTILS::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::NDimConditionSelector(*currdis,"ScaTraCoupling",0,numscal)));
    mcs.SetupExtractor(*currdis,*currdis->DofRowMap(),*mapex);
    scatrafieldexvec_.push_back(mapex);
  }

  scatracoup_->SetupConditionCoupling(*(scatravec_[0]->ScaTraField()->Discretization()),
                                     scatrafieldexvec_[0]->Map(1),
                                     *(scatravec_[1]->ScaTraField()->Discretization()),
                                     scatrafieldexvec_[1]->Map(1),
                                     "ScaTraCoupling",
                                     scatravec_[0]->ScaTraField()->NumScal()); //we assume here that both discretisation have the same number of scalars

  // create map extractor for coupled scatra fields
  // the second field is always split
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;

  // In the limiting case of an infinite permeability of the interface between
  // different scatra fields, the concentrations on both sides of the interface are
  // constrained to be equal. In this case, we keep the fluid scatra dofs at the
  // interface as unknowns in the overall system, whereas the poro (structure) scatra
  // dofs are condensed (cf. "structuresplit" in a monolithic FPSI
  // system). Otherwise, both concentrations are kept in the overall system
  // and the equality of fluxes is considered explicitly.
  if (infperm_)
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->Map(0));
  }
  else
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->FullMap());
  }
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_->Setup(*fullmap,maps);

  // create coupling vectors and matrices (only needed for finite surface permeabilities)
  if (not infperm_)
  {
    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<Epetra_Vector> scatracoupforce =
      Teuchos::rcp(new Epetra_Vector(*(scatraglobalex_->Map(i)),true));
      scatracoupforce_.push_back(scatracoupforce);

      Teuchos::RCP<LINALG::SparseMatrix> scatracoupmat =
        Teuchos::rcp(new LINALG::SparseMatrix(*(scatraglobalex_->Map(i)),27,false,true));
      scatracoupmat_.push_back(scatracoupmat);


      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
      const Epetra_Map* dofrowmap = scatra->ScaTraField()->Discretization()->DofRowMap();
      Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);
      scatrazeros_.push_back(zeros);
    }
  }
  // create scatra block matrix
  scatrasystemmatrix_ =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*scatraglobalex_,
                                                                                   *scatraglobalex_,
                                                                                   27,
                                                                                   false,
                                                                                   true));
  // create scatra rhs vector
  scatrarhs_ = Teuchos::rcp(new Epetra_Vector(*scatraglobalex_->FullMap(),true));
  // create scatra increment vector
  scatraincrement_ = Teuchos::rcp(new Epetra_Vector(*scatraglobalex_->FullMap(),true));
  // check whether potential Dirichlet conditions at scatra interface are
  // defined for both discretizations
  CheckInterfaceDirichletBC();
  // scatra solver
  Teuchos::RCP<DRT::Discretization> firstscatradis = (scatravec_[0])->ScaTraField()->Discretization();
#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<Teuchos::ParameterList> scatrasolvparams = Teuchos::rcp(new Teuchos::ParameterList);
  scatrasolvparams->set("solver","umfpack");
  scatrasolver_ = Teuchos::rcp(new LINALG::Solver(scatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));
#else
  const Teuchos::ParameterList& fs3icontrol = DRT::Problem::Instance()->FS3IControlParams();
  // get solver number used for fs3i
  const int linsolvernumber = fs3icontrol.get<int>("COUPLED_LINEAR_SOLVER");
  // check if LOMA solvers has a valid number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for FS3I problems. Please set COUPLED_LINEAR_SOLVER in FS3I CONTROL to a valid number!");
  const Teuchos::ParameterList& coupledscatrasolvparams =
    DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype =
    DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(coupledscatrasolvparams,"SOLVER");
  if (solvertype != INPAR::SOLVER::aztec_msr)
    dserror("aztec solver expected");
  const int azprectype =
    DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(coupledscatrasolvparams,"AZPREC");
  if (azprectype != INPAR::SOLVER::azprec_BGS2x2)
    dserror("Block Gauss-Seidel preconditioner expected");

  // use coupled scatra solver object
  scatrasolver_ = Teuchos::rcp(new LINALG::Solver(coupledscatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3icontrol.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3icontrol.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    dserror("no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");
  if (linsolver2number == (-1))
    dserror("no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");
  scatrasolver_->PutSolverParamsToSubParams("Inverse1",DRT::Problem::Instance()->SolverParams(linsolver1number));
  scatrasolver_->PutSolverParamsToSubParams("Inverse2",DRT::Problem::Instance()->SolverParams(linsolver2number));
  (scatravec_[0])->ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])->ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse2"));

#endif
}


/*----------------------------------------------------------------------*
 |  Test results                                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fpsi_->FluidField()->CreateFieldTest());

  fpsi_->PoroField()->StructureField()->CreateFieldTest();
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    DRT::Problem::Instance()->AddFieldTest(scatra->CreateScaTraFieldTest());
  }
  DRT::Problem::Instance()->TestAll(comm);

}


/*----------------------------------------------------------------------*
 |  Transfer FPSI solution                                hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetFPSISolution()
{

  SetMeshDisp();
  SetVelocityFields();
  SetWallShearStresses();
  SetPressureFields();


}

/*----------------------------------------------------------------------*
 |  Transfer scatra solution                              hemmler 07/14 |
 *----------------------------------------------------------------------*/
//only needed for two-way coupling; at the moment function is not used
void FS3I::PartFPS3I::SetStructScatraSolution()
{
  fpsi_->PoroField()->StructureField()->Discretization()->SetState(1,"temperature",(scatravec_[1])->ScaTraField()->Phinp());
}


/*----------------------------------------------------------------------*
 |  Set displacements                                     hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetMeshDisp()
{
  // fluid field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra = scatravec_[0];
  ADAPTER::Fluid& fluidadapter = *(fpsi_->FluidField());
  fluidscatra->ScaTraField()->ApplyMeshMovement(fluidadapter.Dispnp(),
                                               fluidadapter.Discretization());

  // Poro field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra = scatravec_[1];
  const Teuchos::RCP<ADAPTER::Fluid>& structadapter = fpsi_->PoroField()->FluidField();
  structscatra->ScaTraField()->ApplyMeshMovement(structadapter->Dispnp(),
                                                structadapter->Discretization());
}


/*----------------------------------------------------------------------*
 |  Set velocities                                        hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetVelocityFields()
{
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& scatradyn       = problem->ScalarTransportDynamicParams();
  int cdvel = DRT::INPUT::IntegralValue<int>(scatradyn,"VELOCITYFIELD");
  switch(cdvel)
  {
    case INPAR::SCATRA::velocity_zero:
    case INPAR::SCATRA::velocity_function:
    case INPAR::SCATRA::velocity_function_and_curve:
    {
      for (unsigned i=0; i<scatravec_.size(); ++i)
      {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
      scatra->ScaTraField()->SetVelocityField();
      }
      break;
    }
    case INPAR::SCATRA::velocity_Navier_Stokes:
    {
      std::vector<Teuchos::RCP<const Epetra_Vector> > convel;
      std::vector<Teuchos::RCP<const Epetra_Vector> > vel;
      ExtractVel(convel, vel);


      std::vector<Teuchos::RCP<DRT::Discretization> > discret;

      discret.push_back(fpsi_->FluidField()->Discretization());
      discret.push_back(fpsi_->PoroField()->FluidField()->Discretization());

      for (unsigned i=0; i<scatravec_.size(); ++i)
      {
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
        scatra->ScaTraField()->SetVelocityField(convel[i],
                                               Teuchos::null,
                                               vel[i],
                                               Teuchos::null,
                                               Teuchos::null,
                                               discret[i]);
      }
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  Set wall shear stresses                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetWallShearStresses()
{
  std::vector<Teuchos::RCP<const Epetra_Vector> > wss;
  ExtractWSS(wss);

  std::vector<Teuchos::RCP<DRT::Discretization> > discret;

  discret.push_back(fpsi_->FluidField()->Discretization());
  discret.push_back(fpsi_->PoroField()->FluidField()->Discretization());

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->SetWallShearStresses(wss[i],Teuchos::null,discret[i]);
  }
}

/*----------------------------------------------------------------------*
 |  Set presures                                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetPressureFields()
{
  std::vector<Teuchos::RCP<const Epetra_Vector> > pressure;
  ExtractPressure(pressure);

  std::vector<Teuchos::RCP<DRT::Discretization> > discret;

  discret.push_back(fpsi_->FluidField()->Discretization());
  discret.push_back(fpsi_->PoroField()->FluidField()->Discretization());

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->SetPressureFields(pressure[i],Teuchos::null,discret[i]);
  }
}

/*----------------------------------------------------------------------*
 |  Set mean concentration                                hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetMeanConcentration()
{
    std::vector<Teuchos::RCP<Epetra_Vector> > MeanConc;
    ExtractMeanConcentration(MeanConc);

    std::vector<Teuchos::RCP<DRT::Discretization> > discret;

    discret.push_back(scatravec_[0]->ScaTraField()->Discretization());
    discret.push_back(scatravec_[1]->ScaTraField()->Discretization());

    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
      scatra->ScaTraField()->SetMeanConcentration(MeanConc[i]);
    }
}

/*----------------------------------------------------------------------*
 |  Evaluate scatra fields                                hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::EvaluateScatraFields()
{
  //mean concentration at the interface needed for membrane equation of Kedem and Katchalsky.
  //needs to be set here, since it depends on the scalar interface values on both discretisations
  SetMeanConcentration();

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_adap = scatravec_[i];
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatra_adap->ScaTraField();

    scatra->PrepareLinearSolve();

    // add contributions due to finite interface permeability
    if (!infperm_)
    {
      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<LINALG::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      // evaluate interface; second Kedem-Katchalsky equation for coupling of solute flux
      scatra->KedemKatchalsky(coupmat,coupforce);

      // apply Dirichlet boundary conditions to coupling matrix and vector
      Teuchos::RCP<Epetra_Vector> zeros = scatrazeros_[i];
      const Teuchos::RCP<const LINALG::MapExtractor> dbcmapex = scatra->DirichMaps();
      const Teuchos::RCP< const Epetra_Map > dbcmap = dbcmapex->CondMap();
      coupmat->ApplyDirichlet(*dbcmap,false);
      LINALG::ApplyDirichlettoSystem(coupforce,zeros,*dbcmap);
    }
  }
}


/*----------------------------------------------------------------------*
 |  Extract velocities                                    hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& convel,
                      std::vector<Teuchos::RCP<const Epetra_Vector> >& vel)
{
  //############ Fluid Field ###############
  convel.push_back(fpsi_->FluidField()->ConvectiveVel());
  vel.push_back(fpsi_->FluidField()->Velnp());

  //############ Poro Field ###############
  convel.push_back(fpsi_->PoroField()->FluidField()->ConvectiveVel());
  vel.push_back(fpsi_->PoroField()->FluidField()->Velnp());
}


/*----------------------------------------------------------------------*
 |  Extract wall shear stresses                           hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss)
{
  //############ Fluid Field ###############
  Teuchos::RCP<Epetra_Vector> WallShearStress = LINALG::CreateVector(*(fpsi_->FluidField()->Discretization()->DofRowMap()),true);
  WallShearStress = CalcWallShearStress();
  wss.push_back(WallShearStress);

  //############ Poro Field ###############

  //Hint: The Wall shear stresses in the fluid field at the Interface are equal to the ones of the poro structure
  //      Therefore, we map the results of the wss of the fluid field to the dofs of the poro field without computing
  //      them explicitly in the poro field

  // extract FPSI-Interface from fluid field
  WallShearStress = fpsi_->FluidField()->FPSIInterface()->ExtractFPSICondVector(WallShearStress);

  // replace global fluid interface dofs through porofluid interface dofs
  WallShearStress = fpsi_->FluidToPorofluid_FPSI(WallShearStress);

  // insert porofluid interface entries into vector with full porofluid length
  Teuchos::RCP<Epetra_Vector> porofluid = LINALG::CreateVector(*(fpsi_->PoroField()->Interface().FullMap()),true);

  //Parameter int block of function InsertVector: (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
  fpsi_->PoroField()->Interface().InsertVector(WallShearStress,3,porofluid);
  wss.push_back(porofluid);
}

/*----------------------------------------------------------------------*
 |  Extract pressures                                     hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractPressure(std::vector<Teuchos::RCP<const Epetra_Vector> >& pressure)
{
  //############ Fluid Field ###############
  pressure.push_back(fpsi_->FluidField()->Velnp()); //this actually extracts the velocities as well. We sort them out later.

  //############ Poro Field ###############
  pressure.push_back(fpsi_->PoroField()->FluidField()->Velnp()); //this actually extracts the velocities as well. We sort them out later.
}

/*----------------------------------------------------------------------*
 |  Extract mean concentration                            hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractMeanConcentration(std::vector<Teuchos::RCP<Epetra_Vector> >& MeanConc)
{
  //############ Fluid Field ###############
  Teuchos::RCP<Epetra_Vector> MeanConcentration1 = CalcMeanConcentration();
  MeanConc.push_back(MeanConcentration1);

  //############ Poro Field ###############
  //Hint: The mean concentration is not calculated again; we just map the values from the Fluid Field into the Poro Field

  //extract interface values
  Teuchos::RCP<Epetra_Vector> interface_phin = scatrafieldexvec_[0]->ExtractVector(MeanConcentration1,1);

  //insert interface values from Fluid Field into Poro Field;
  Teuchos::RCP<Epetra_Vector> MeanConcentration2 = scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(interface_phin),1);
  MeanConc.push_back(MeanConcentration2);
}


/*----------------------------------------------------------------------*
 |  Calculate wall shear stresses                         hemmler 07/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::PartFPS3I::CalcWallShearStress()
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action",FLD::ba_calc_node_normal);
  std::string condstring("ScaTraCoupling");

  //get rowmap
  Teuchos::RCP<DRT::Discretization> discret = fpsi_->FluidField()->Discretization();
  const Epetra_Map* dofrowmap = discret->DofRowMap();

  //vector ndnorm0 with pressure-entries is needed for EvaluateCondition
  Teuchos::RCP<Epetra_Vector> ndnorm0 = LINALG::CreateVector(*dofrowmap,true);

  // evaluate the normals of the surface, note: normal vectors do not yet have length = 1.0
  discret->ClearState();
  Teuchos::RCP<const Epetra_Vector> dispnp = LINALG::CreateVector(*dofrowmap,true);
  dispnp=fpsi_->FluidField()->Dispnp();
  discret->SetState("dispnp", dispnp);
  discret->EvaluateCondition(eleparams,ndnorm0,condstring);
  discret->ClearState();

  //dimensions
  const int ndim = DRT::Problem::Instance()->NDim();

  // normalize the normal vectors
  for (int i = 0; i < ndnorm0->MyLength();i+=ndim+1)
  {
    // calculate the length of the normal
    double L = 0.0;
    for (int j = 0; j<ndim; j++)
    {
      L += ((*ndnorm0)[i+j])*((*ndnorm0)[i+j]);
    }
    L = sqrt(L);

    // normalise the normal vector (if present for the current node)
    if (L > EPS15)
    {
      for (int j = 0; j < ndim; j++)
      {
        (*ndnorm0)[i+j] /=  L;
      }
    }
  }

  // evaluate the wall shear stress from the traction by removing
  // the normal stresses

  // set action for elements
  eleparams.set<int>("action",FLD::integrate_Shapefunction);
  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = LINALG::CreateVector(*dofrowmap,true);

  //integrate shape functions
  discret->ClearState();
  discret->SetState("dispnp", dispnp);
  discret->EvaluateCondition(eleparams,integratedshapefunc,condstring);
  discret->ClearState();

  //true (rescaled) residual vector without zeros at dirichlet positions (Unit: Newton [N])
  Teuchos::RCP<const Epetra_Vector> trueresidual  =fpsi_->FluidField()->TrueResidual();

  // compute traction values at specified nodes; otherwise do not touch the zero values
  for (int i=0;i<integratedshapefunc->MyLength();i++)
  {
    if ((*integratedshapefunc)[i] != 0.0)
    {
      // overwrite integratedshapefunc values with the calculated traction coefficients,
      // which are reconstructed out of the nodal forces (trueresidual_) using the
      // same shape functions on the boundary as for velocity and pressure.
      (*integratedshapefunc)[i] = (*trueresidual)[i]/(*integratedshapefunc)[i];
    }
  }

  //this is just a renaming
  Teuchos::RCP<Epetra_Vector> wss = integratedshapefunc;

  // loop over all entities within the traction vector
  for (int i = 0; i < ndnorm0->MyLength();i+=ndim+1)
  {
    // evaluate the normal stress = < traction . normal >
    double normal_stress = 0.0;
    for (int j = 0; j<ndim; j++)
    {
      normal_stress += (*wss)[i+j] * (*ndnorm0)[i+j];
    }

    // subtract the normal stresses from traction
    for (int j = 0; j<ndim; j++)
    {
      (*wss)[i+j] -= normal_stress * (*ndnorm0)[i+j];
    }
  }

  // return the wall_shear_stress vector
  return wss;
}

/*----------------------------------------------------------------------*
 |  Calculate mean concentration                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::PartFPS3I::CalcMeanConcentration()
{
  //Get concentration phi2 in scatrafield2
  //Hint: in the following we talk of phi1 and phi2, but they mean the same concentration just on different scatrafields
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra2 = scatravec_[1];
  Teuchos::RCP<Epetra_Vector> scatrafield2_phi2n= scatra2->ScaTraField()->Phin();

  //extract interface values from phi2 but we are still on scatrafield2
  Teuchos::RCP<Epetra_Vector> interface2_phi2n = scatrafieldexvec_[1]->ExtractVector(scatrafield2_phi2n,1);

  //insert interface values from scatrafield2 into scatrafield1; scatrafield1_phi2n is again of full length, i.e. of size of scatrafield1; all values that do not belong to the interface are zero
  Teuchos::RCP<Epetra_Vector> scatrafield1_phi2n = scatrafieldexvec_[0]->InsertVector(Scatra2ToScatra1(interface2_phi2n),1);

  //Get concentration phi1 in scatrafield1
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra1 = scatravec_[0];
  Teuchos::RCP<Epetra_Vector> scatrafield1_phi1n= scatra1->ScaTraField()->Phin();

  //extract interface values from phi1 but we are still on scatrafield1
  Teuchos::RCP<Epetra_Vector> interface1_phi1n = scatrafieldexvec_[0]->ExtractVector(scatrafield1_phi1n,1);

  //insert interface values interface1_phi1n from scatrafield1 into the full scatrafield1 again; this is just to obtain a vector whose entries are zero except for the nodes of the interface
  Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0]->InsertVector(interface1_phi1n,1);

  //nodewise calculation of mean concentration in the interface

  for(int i=0; i<temp->MyLength();i++)
  {
    //here the unweighted average is uses. One could also use a logarithmic average...
    (*temp)[i]=0.5*((*temp)[i]+(*scatrafield1_phi2n)[i]); //log. average: ((*temp)[i]-(*scatrafield1_phi2n)[i])/log(((*temp)[i])/((*scatrafield1_phi2n)[i]));
                                                          //linear approach: 0.5*((*temp)[i]+(*scatrafield1_phi2n)[i]);
  }

  //return mean concentration in the interface
  //this vector now belongs to scatrafield1!!!
  return temp;
}

