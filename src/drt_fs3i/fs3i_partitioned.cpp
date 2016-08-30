/*!----------------------------------------------------------------------
\file fs3i_partitioned.cpp
\brief General algorithmic routines for partitioned solution approaches
       to fluid-structure-scalar-scalar interaction (FS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively

\level 2

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364


*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>
//#include <initializer_list>

#include "fs3i_partitioned.H"

//IO
#include "../drt_io/io_control.H"
//FSI
#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_condition_utils.H"
//LINALG
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
//INPAR
#include "../drt_inpar/drt_validparameters.H"
//ALE
#include "../drt_ale/ale_utils_clonestrategy.H"
//ADAPTER
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_adapter/adapter_coupling_volmortar.H"
//SCATRA
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"
//FOR WSS CALCULATIONS
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"


//#define SCATRABLOCKMATRIXMERGE

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I::PartFS3I(const Epetra_Comm& comm)
  : FS3I_Base(),
    comm_(comm)
{
  volume_fieldcouplings_.push_back(DRT::INPUT::IntegralValue<INPAR::FS3I::VolumeCoupling>(DRT::Problem::Instance()->FS3IDynamicParams(),"FLUIDSCAL_FIELDCOUPLING"));
  volume_fieldcouplings_.push_back(DRT::INPUT::IntegralValue<INPAR::FS3I::VolumeCoupling>(DRT::Problem::Instance()->FS3IDynamicParams(),"STRUCTSCAL_FIELDCOUPLING"));

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fs3idyn = problem->FS3IDynamicParams();

  //---------------------------------------------------------------------
  // ensure correct order of three discretizations, with dof-numbering
  // such that structure dof < fluid dof < ale dofs
  // (ordering required at certain non-intuitive points)
  //---------------------------------------------------------------------
  problem->GetDis("structure")->FillComplete();
  problem->GetDis("fluid")->FillComplete();
  problem->GetDis("ale")->FillComplete();
  problem->GetDis("scatra1")->FillComplete();
  problem->GetDis("scatra2")->FillComplete();

  //---------------------------------------------------------------------
  // access discretizations for structure, fluid, ale as well as fluid-
  // and structure-based scalar transport and get material map for fluid
  // and scalar transport elements
  //---------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluidscatradis = problem->GetDis("scatra1");
  Teuchos::RCP<DRT::Discretization> structscatradis = problem->GetDis("scatra2");
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  //---------------------------------------------------------------------
  // create ale discretization as a clone from fluid discretization
  //---------------------------------------------------------------------
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else
    dserror("Providing an ALE mesh is not supported for FS3I problems.");

  //std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap = problem->CloningMaterialMap();
  //if (clonefieldmatmap.size() < 2)
  //  dserror("At least two material lists required for partitioned FS3I!");

  // determine type of scalar transport
  const INPAR::SCATRA::ImplType impltype_fluid =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(DRT::Problem::Instance()->FS3IDynamicParams(),"FLUIDSCAL_SCATRATYPE");

  //---------------------------------------------------------------------
  // create discretization for fluid-based scalar transport from and
  // according to fluid discretization
  //---------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  // create fluid-based scalar transport elements if fluid-based scalar
  // transport discretization is empty
  if (fluidscatradis->NumGlobalNodes()==0)
  {
    if ( not (volume_fieldcouplings_[0]==INPAR::FS3I::coupling_match) )
      dserror("If you clone your fluid-scatra mesh from the fluid use FLUIDSCAL_FIELDCOUPLING 'volume_matching'!");

    // fill fluid-based scatra discretization by cloning fluid discretization
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,fluidscatradis);

    // set implementation type of cloned scatra elements to advanced reactions
    for(int i=0; i<fluidscatradis->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(fluidscatradis->lColElement(i));
      if(element == NULL)
        dserror("Invalid element type!");
      else
        element->SetImplType(impltype_fluid);
    }

    volume_coupling_objects_.push_back(Teuchos::null);

    // care for secondary dof sets:
    // add proxy of fluid degrees of freedom to scatra discretization
    if(fluidscatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
      dserror("Fluid scatra discretization has illegal number of dofsets!");
  }
  else
  {
    if ( not (volume_fieldcouplings_[0]==INPAR::FS3I::coupling_nonmatch) )
      dserror("If you have specified the fluid-scalar by TRANSPORT ELEMENTS use FLUIDSCAL_FIELDCOUPLING 'volume_nonmatching'!");

    if ( not (impltype_fluid == INPAR::SCATRA::impltype_undefined) )
      dserror("Be aware that your FLUIDSCAL_SCATRATYPE will be ignored and the impltype from the TRANSPORT ELMENTS section will be utilized. Use FLUIDSCAL_SCATRATYPE 'Undefined'!");

    volume_coupling_objects_.push_back( CreateVolMortarObject(fluiddis,fluidscatradis) );

    dserror("Mortar volume coupling for the fluid-scalar is yet not tested. So be careful!");
  }


  // determine type of scalar transport
    const INPAR::SCATRA::ImplType impltype_struct =
        DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(DRT::Problem::Instance()->FS3IDynamicParams(),"STRUCTSCAL_SCATRATYPE");

  //---------------------------------------------------------------------
  // create discretization for structure-based scalar transport from and
  // according to structure discretization
  //---------------------------------------------------------------------
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create structure-based scalar transport elements if structure-based
  // scalar transport discretization is empty
  if (structscatradis->NumGlobalNodes()==0)
  {
    if ( not (volume_fieldcouplings_[1]==INPAR::FS3I::coupling_match) )
      dserror("If you clone your structure-scatra mesh from the structure use STRUCTSCAL_FIELDCOUPLING 'volume_matching'!");

    // fill structure-based scatra discretization by cloning structure discretization
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,structscatradis);

    // set implementation type of cloned scatra elements to advanced reactions
    for(int i=0; i<structscatradis->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(structscatradis->lColElement(i));
      if(element == NULL)
        dserror("Invalid element type!");
      else
        element->SetImplType(impltype_struct);
    }

    volume_coupling_objects_.push_back(Teuchos::null);

    // care for secondary dof sets:
    // add proxy of structure scatra degrees of freedom to structure discretization
    if (structdis->AddDofSet( structscatradis->GetDofSetProxy() )!=1)
      dserror("Structure discretization has illegal number of dofsets!");

    // add proxy of structure degrees of freedom to scatra discretization
    if(structscatradis->AddDofSet(structdis->GetDofSetProxy()) != 1)
      dserror("Structure scatra discretization has illegal number of dofsets!");
  }
  else
  {
    if ( not (volume_fieldcouplings_[1]==INPAR::FS3I::coupling_nonmatch) )
      dserror("If you have specified the structure-scalar by TRANSPORT2 ELEMENTS use STRUCTSCAL_FIELDCOUPLING 'volume_nonmatching'!");

    if ( not (impltype_struct == INPAR::SCATRA::impltype_undefined) )
          dserror("Be aware that your STRUCTSCAL_SCATRATYPE will be ignored and the impltype from the TRANSPORT2 ELMENTS section will be utilized. Use STRUCTSCAL_SCATRATYPE 'Undefined'!");

    volume_coupling_objects_.push_back( CreateVolMortarObject(structdis,structscatradis) );
  }

  //safety check
  if ( not (volume_coupling_objects_.size() == 2) )
    dserror("Unexpected size of volmortar object vector!");

  //Note: in the scatra fields we have now the following dof-sets:
  // structure dofset 0: structure dofset
  // structure dofset 1: structscatra dofset
  // fluidscatra dofset 0: fluidscatra dofset
  // fluidscatra dofset 1: fluid dofset
  // structscatra dofset 0: structscatra dofset
  // structscatra dofset 1: structure dofset

  //---------------------------------------------------------------------
  // get FSI coupling algorithm
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();
  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");

  const Teuchos::ParameterList& fsitimeparams = ManipulateFsiTimeParams(fs3idyn);

  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    {
      fsi_ = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm,fsitimeparams));
      break;
    }
    case fsi_iter_monolithicstructuresplit:
    {
      fsi_ = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm,fsitimeparams));
      break;
    }
    default:
    {
      dserror("Unknown FSI coupling algorithm");
      break;
    }
  }

  //---------------------------------------------------------------------
  // create instances for fluid- and structure-based scalar transport
  // solver and arrange them in combined vector
  //---------------------------------------------------------------------
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3idyn.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3idyn.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    dserror("no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I DYNAMIC to a valid number!");
  if (linsolver2number == (-1))
    dserror("no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I DYNAMIC to a valid number!");

  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(fs3idyn,problem->ScalarTransportDynamicParams(),problem->SolverParams(linsolver1number),"scatra1",true));
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(fs3idyn,problem->ScalarTransportDynamicParams(),problem->SolverParams(linsolver2number),"scatra2",true));

  scatravec_.push_back(fluidscatra);
  scatravec_.push_back(structscatra);

  //---------------------------------------------------------------------
  // check existence of scatra coupling conditions for both
  // discretizations and definition of the permeability coefficient
  //---------------------------------------------------------------------
  CheckFS3IInputs();

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP< ::ADAPTER::MortarVolCoupl> FS3I::PartFS3I::CreateVolMortarObject(
    Teuchos::RCP<DRT::Discretization> masterdis,
    Teuchos::RCP<DRT::Discretization> slavedis )
{
  // copy conditions
  // this is actually only needed for copying TRANSPORT DIRICHLET/NEUMANN CONDITIONS
  // as standard DIRICHLET/NEUMANN CONDITIONS
  std::map<std::string,std::string> conditions_to_copy;
  SCATRA::ScatraFluidCloneStrategy clonestrategy;
  conditions_to_copy = clonestrategy.ConditionsToCopy();
  DRT::UTILS::DiscretizationCreatorBase creator;
  creator.CopyConditions(*slavedis,*slavedis,conditions_to_copy);

  //first call FillComplete for single discretizations.
  //This way the physical dofs are numbered successively
  masterdis->FillComplete();
  slavedis->FillComplete();

  //build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = slavedis->NumDof(0,slavedis->lRowNode(0));
  const int ndofperelement_scatra  = 0;
  const int ndofpernode_struct = masterdis->NumDof(0,masterdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  if (masterdis->BuildDofSetAuxProxy(ndofpernode_scatra, ndofperelement_scatra, 0, true ) != 1)
    dserror("unexpected dof sets in structure field");
  if (slavedis->BuildDofSetAuxProxy(ndofpernode_struct, ndofperelement_struct, 0, true) != 1)
    dserror("unexpected dof sets in scatra field");

  //call AssignDegreesOfFreedom also for auxiliary dofsets
  //note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  masterdis->FillComplete(true, false,false);
  slavedis->FillComplete(true, false,false);


  // Scheme: non matching meshes --> volumetric mortar coupling...
  Teuchos::RCP< ::ADAPTER::MortarVolCoupl> volume_coupling_object = Teuchos::rcp(new ADAPTER::MortarVolCoupl() );

  //setup projection matrices (use default material strategy)
  volume_coupling_object->Setup(masterdis,slavedis);

  return volume_coupling_object;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList& FS3I::PartFS3I::ManipulateFsiTimeParams(const Teuchos::ParameterList& fs3idyn) const
{
  // NOTE: we can not do this in the AC-fs3i class were it would belong,
  // since overloading a function inside the constructor does not work :(

  Teuchos::ParameterList& timeparams= *( new Teuchos::ParameterList(fs3idyn));

  const int fsisubcycles = fs3idyn.sublist("AC").get<int>("FSI_STEPS_PER_SCATRA_STEP");

  if (fsisubcycles != 1) //if we have subcycling for ac_fsi
  {
    timeparams.set<double>("TIMESTEP", fs3idyn.get<double>("TIMESTEP")/(double)fsisubcycles);
  }
  return timeparams;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ReadRestart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  const int restart = DRT::Problem::Instance()->Restart();

  if (restart)
  {
    const Teuchos::ParameterList& fs3idynac = DRT::Problem::Instance()->FS3IDynamicParams();
    const bool restartfrompartfsi = DRT::INPUT::IntegralValue<int>(fs3idynac,"RESTART_FROM_PART_FSI");

    if (not restartfrompartfsi) //standard restart
    {
      fsi_->ReadRestart(restart);

      for (unsigned i=0; i<scatravec_.size(); ++i)
      {
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
        currscatra->ScaTraField()->ReadRestart(restart);
      }
    }
    else //we do not want to read the scatras values and the lagrange multiplyer, since we start from a partitioned FSI
    {
      fsi_->ReadRestart(restart);

      // We may want to calculate the initial time deriv, since the scatra field is new.
      // This has to be done before setting the step in the scatra fields, because step_ it needs to be zero.
      const Teuchos::ParameterList& scatradynparams = DRT::Problem::Instance()->ScalarTransportDynamicParams();
      const bool skiptimederiv = DRT::INPUT::IntegralValue<int>(scatradynparams,"SKIPINITDER");

      if (not skiptimederiv)
      {
        SetMeshDisp();
        SetVelocityFields();
        //SetWallShearStresses(); //Note: We can not do this since we don't have trueresidual_ in the fluid discretisation

        for (unsigned i=0; i<scatravec_.size(); ++i)
        {
          Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
          currscatra->ScaTraField()->PrepareFirstTimeStep();
        }
      }

      //we need to set time and step in scatra to have it matching with the fsi ones
      for (unsigned i=0; i<scatravec_.size(); ++i)
      {
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
        currscatra->ScaTraField()->SetTimeStep(fsi_->FluidField()->Time(),fsi_->FluidField()->Step());
      }
    }

    time_ = fsi_->FluidField()->Time();
    step_ = fsi_->FluidField()->Step();
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetupSystem()
{
  // now do the coupling setup and create the combined dofmap
  fsi_->SetupSystem();

  /*----------------------------------------------------------------------*/
  /*                  General set up for scalar fields                    */
  /*----------------------------------------------------------------------*/

  // create map extractors needed for scatra condition coupling
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
    Teuchos::RCP<DRT::Discretization> currdis = currscatra->ScaTraField()->Discretization();
    const int numscal = currscatra->ScaTraField()->NumScal();
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
  // the second field (currently structure) is always split
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;

  // In the limiting case of an infinite permeability of the interface between
  // different scatra fields, the concentrations on both sides of the interface are
  // constrained to be equal. In this case, we keep the fluid scatra dofs at the
  // interface as unknowns in the overall system, whereas the structure scatra
  // dofs are condensed (cf. "structuresplit" in a monolithic FSI
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
  if (!infperm_)
  {
    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<Epetra_Vector> scatracoupforce =
      Teuchos::rcp(new Epetra_Vector(*(scatraglobalex_->Map(i)),true));
      scatracoupforce_.push_back(scatracoupforce);

      Teuchos::RCP<LINALG::SparseMatrix> scatracoupmat =
        Teuchos::rcp(new LINALG::SparseMatrix(*(scatraglobalex_->Map(i)),27,false,true));
      scatracoupmat_.push_back(scatracoupmat);

      const Epetra_Map* dofrowmap = scatravec_[i]->ScaTraField()->Discretization()->DofRowMap();
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
  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  // get solver number used for fs3i
  const int linsolvernumber = fs3idyn.get<int>("COUPLED_LINEAR_SOLVER");
  // check if LOMA solvers has a valid number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for FS3I problems. Please set COUPLED_LINEAR_SOLVER in FS3I DYNAMIC to a valid number!");

  const Teuchos::ParameterList& coupledscatrasolvparams = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(coupledscatrasolvparams,"SOLVER");
  if (solvertype != INPAR::SOLVER::aztec_msr)
    dserror("aztec solver expected");

  const int azprectype = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(coupledscatrasolvparams,"AZPREC");
  if (azprectype != INPAR::SOLVER::azprec_BGS2x2)
    dserror("Block Gauss-Seidel preconditioner expected");

  // use coupled scatra solver object
  scatrasolver_ = Teuchos::rcp(new LINALG::Solver(coupledscatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));

  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3idyn.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3idyn.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    dserror("no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER1 in FS3I DYNAMIC to a valid number!");
  if (linsolver2number == (-1))
    dserror("no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I DYNAMIC to a valid number!");

  scatrasolver_->PutSolverParamsToSubParams("Inverse1",DRT::Problem::Instance()->SolverParams(linsolver1number));
  scatrasolver_->PutSolverParamsToSubParams("Inverse2",DRT::Problem::Instance()->SolverParams(linsolver2number));

  (scatravec_[0])->ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])->ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse2"));
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fsi_->FluidField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(fsi_->AleField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(fsi_->StructureField()->CreateFieldTest());

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    DRT::Problem::Instance()->AddFieldTest(scatra->CreateScaTraFieldTest());
  }
  DRT::Problem::Instance()->TestAll(comm);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetFSISolution()
{
  //we clear every state, including the states of the secondary dof sets
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    scatravec_[i]->ScaTraField()->Discretization()->ClearState(true);
    // we have to manually clear this since this can not be saved directly in the
    // primary dof set (because it is cleared in between)
    scatravec_[i]->ScaTraField()->ClearExternalConcentrations();
  }

  SetMeshDisp();
  SetVelocityFields();
  SetWallShearStresses();
  SetMembraneConcentration();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetStructScatraSolution() const
{
  fsi_->StructureField()->Discretization()->SetState( 1,"temperature",StructureScalarToStructure(scatravec_[1]->ScaTraField()->Phinp()) );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetMeshDisp() const
{
  // fluid field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra = scatravec_[0];
  fluidscatra->ScaTraField()->ApplyMeshMovement( FluidToFluidScalar(fsi_->FluidField()->Dispnp()),1 );

  // structure field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra = scatravec_[1];
  structscatra->ScaTraField()->ApplyMeshMovement( StructureToStructureScalar(fsi_->StructureField()->Dispnp()),1 );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetVelocityFields() const
{
  std::vector<Teuchos::RCP<const Epetra_Vector> > convel;
  std::vector<Teuchos::RCP<const Epetra_Vector> > vel;
  ExtractVel(convel, vel);

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->SetVelocityField(VolMortarMasterToSlavei(i,convel[i]),
                                           Teuchos::null,
                                           VolMortarMasterToSlavei(i,vel[i]),
                                           Teuchos::null,
                                           1);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& convel,
                      std::vector<Teuchos::RCP<const Epetra_Vector> >& vel) const
{
  // extract fluid velocities

  convel.push_back(fsi_->FluidField()->ConvectiveVel());
  vel.push_back(fsi_->FluidField()->Velnp());

  // extract structure velocities and accelerations

  Teuchos::RCP<Epetra_Vector> velocity = Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Velnp())));
  vel.push_back(velocity);
  // structure ScaTra: velocity and grid velocity are identical!
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(velocity->Map(),true));
  convel.push_back(zeros);
}

/*----------------------------------------------------------------------*
 |  Set wall shear stresses                                  Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetWallShearStresses() const
{
  std::vector<Teuchos::RCP<const Epetra_Vector> > wss;
  ExtractWSS(wss);

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->SetWallShearStresses( VolMortarMasterToSlavei(i,wss[i]),1 );
  }
}

/*----------------------------------------------------------------------*
 |  Extract wall shear stresses                              Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFS3I::ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector> >& wss) const
{
  //############ Fluid Field ###############

  Teuchos::RCP<ADAPTER::FluidFSI> fluid = Teuchos::rcp_dynamic_cast<ADAPTER::FluidFSI>( fsi_->FluidField() );
  if (fluid == Teuchos::null)
    dserror("Dynamic cast to ADAPTER::FluidFSI failed!");

  Teuchos::RCP<Epetra_Vector> WallShearStress = fluid->CalculateWallShearStresses();

  if ( DRT::INPUT::IntegralValue<INPAR::FLUID::WSSType>(DRT::Problem::Instance()->FluidDynamicParams() ,"WSS_TYPE") != INPAR::FLUID::wss_standard)
    dserror("WSS_TYPE not supported for FS3I!");

  wss.push_back(WallShearStress);

  //############ Structure Field ###############

  // extract FSI-Interface from fluid field
  WallShearStress = fsi_->FluidField()->Interface()->ExtractFSICondVector(WallShearStress);

  // replace global fluid interface dofs through structure interface dofs
  WallShearStress = fsi_->FluidToStruct(WallShearStress);

  // insert structure interface entries into vector with full structure length
  Teuchos::RCP<Epetra_Vector> structure = LINALG::CreateVector(*(fsi_->StructureField()->Interface()->FullMap()),true);

  //Parameter int block of function InsertVector: (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
  fsi_->StructureField()->Interface()->InsertVector(WallShearStress,1,structure);
  wss.push_back(structure);
}

/*----------------------------------------------------------------------*
 |  transport quantity from fluid to fluid-scalar            Thon 08/16 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::FluidToFluidScalar(const Teuchos::RCP<const Epetra_Vector> fluidvector) const
{
  return VolMortarMasterToSlavei(0,fluidvector);
}

/*----------------------------------------------------------------------*
 |  transport quantity from fluid-scalar to fluid            Thon 08/16 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::FluidScalarToFluid(const Teuchos::RCP<const Epetra_Vector> fluidscalarvector) const
{
  return VolMortarSlaveToMasteri(0,fluidscalarvector);
}

/*----------------------------------------------------------------------*
 |  transport quantity from structure to structure-scalar    Thon 08/16 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::StructureToStructureScalar(const Teuchos::RCP<const Epetra_Vector> structurevector) const
{
  return VolMortarMasterToSlavei(1,structurevector);
}

/*----------------------------------------------------------------------*
 |  transport quantity from structure-scalar to structure    Thon 08/16 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::StructureScalarToStructure(const Teuchos::RCP<const Epetra_Vector> structurescalavector) const
{
  return VolMortarSlaveToMasteri(1,structurescalavector);
}

/*-------------------------------------------------------------------------------------*
 |  transport quantity from i-th volmortar master to i-th volmortar slave   Thon 08/16 |
 *-------------------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::VolMortarMasterToSlavei(const int i, const Teuchos::RCP<const Epetra_Vector> mastervector) const
{
  switch(volume_fieldcouplings_[i])
  {
  case INPAR::FS3I::coupling_match:
    return mastervector;
    break;
  case INPAR::FS3I::coupling_nonmatch:
    return volume_coupling_objects_[i]->ApplyVectorMapping21(mastervector);
    break;
  default:
    dserror("unknown field coupling type");
    return Teuchos::null;
    break;
  }
}

/*-------------------------------------------------------------------------------------*
 |  transport quantity from i-th volmortar slave to i-th volmortar master   Thon 08/16 |
 *-------------------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> FS3I::PartFS3I::VolMortarSlaveToMasteri(const int i, const Teuchos::RCP<const Epetra_Vector> slavevector) const
{
  switch(volume_fieldcouplings_[i])
  {
  case INPAR::FS3I::coupling_match:
    return slavevector;
    break;
  case INPAR::FS3I::coupling_nonmatch:
    return volume_coupling_objects_[i]->ApplyVectorMapping12(slavevector);
    break;
  default:
    dserror("unknown field coupling type");
    return Teuchos::null;
    break;
  }
}

