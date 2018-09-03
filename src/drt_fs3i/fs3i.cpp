/*!----------------------------------------------------------------------
\file fs3i.cpp
\brief cpp-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I) and fluid-porous-structure-scalar-scalar
       interaction (FPS3I).

\level 2

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364

*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fsi/fsi_dyn.H"
#include "../drt_fsi/fs_monolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fsi/fsi_utils.H"

#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_fs3i.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/fluidresulttest.H"

#include "../drt_ssi/ssi_clonestrategy.H"

#include "fs3i.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::FS3I_Base::FS3I_Base()
    : infperm_(DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->FS3IDynamicParams(), "INF_PERM")),
      timemax_(DRT::Problem::Instance()->FS3IDynamicParams().get<double>("MAXTIME")),
      numstep_(DRT::Problem::Instance()->FS3IDynamicParams().get<int>("NUMSTEP")),
      dt_(DRT::Problem::Instance()->FS3IDynamicParams().get<double>("TIMESTEP")),
      time_(0.0),
      step_(0),
      issetup_(false),
      isinit_(false)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::Init()
{
  SetIsSetup(false);

  scatracoup_ = Teuchos::rcp(new ADAPTER::Coupling());
  scatraglobalex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  sbbtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform());
  sbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());
  sibtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform());
  fbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());

  SetIsInit(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::Setup()
{
  CheckIsInit();

  SetIsSetup(true);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::CheckInterfaceDirichletBC()
{
  Teuchos::RCP<DRT::Discretization> masterdis = scatravec_[0]->ScaTraField()->Discretization();
  Teuchos::RCP<DRT::Discretization> slavedis = scatravec_[1]->ScaTraField()->Discretization();

  Teuchos::RCP<const Epetra_Map> mastermap = scatracoup_->MasterDofMap();
  Teuchos::RCP<const Epetra_Map> permmastermap = scatracoup_->PermMasterDofMap();
  Teuchos::RCP<const Epetra_Map> slavemap = scatracoup_->SlaveDofMap();
  Teuchos::RCP<const Epetra_Map> permslavemap = scatracoup_->PermSlaveDofMap();

  const Teuchos::RCP<const LINALG::MapExtractor> masterdirichmapex =
      scatravec_[0]->ScaTraField()->DirichMaps();
  const Teuchos::RCP<const Epetra_Map> masterdirichmap = masterdirichmapex->CondMap();

  // filter out master dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> masterifdirich = Teuchos::rcp(new Epetra_Vector(*mastermap, true));
  for (int i = 0; i < mastermap->NumMyElements(); ++i)
  {
    int gid = mastermap->GID(i);
    if (masterdirichmap->MyGID(gid))
    {
      (*masterifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_slaveifdirich = scatracoup_->MasterToSlave(masterifdirich);

  const Teuchos::RCP<const LINALG::MapExtractor> slavedirichmapex =
      scatravec_[1]->ScaTraField()->DirichMaps();
  const Teuchos::RCP<const Epetra_Map> slavedirichmap = slavedirichmapex->CondMap();

  // filter out slave dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> slaveifdirich = Teuchos::rcp(new Epetra_Vector(*slavemap, true));
  for (int i = 0; i < slavemap->NumMyElements(); ++i)
  {
    int gid = slavemap->GID(i);
    if (slavedirichmap->MyGID(gid))
    {
      (*slaveifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_masterifdirich = scatracoup_->SlaveToMaster(slaveifdirich);

  // check if the locations of non-zero entries do not match
  for (int i = 0; i < slavedis->DofRowMap()->NumMyElements(); ++i)
  {
    int gid = slavedis->DofRowMap()->GID(i);
    if (slavemap->MyGID(gid))  // in this case, the current dof is part of the interface
    {
      if ((*test_slaveifdirich)[slavemap->LID(gid)] == 1.0 and
          (*slaveifdirich)[slavemap->LID(gid)] != 1.0)
      {
        dserror("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }

  for (int i = 0; i < masterdis->DofRowMap()->NumMyElements(); ++i)
  {
    int gid = masterdis->DofRowMap()->GID(i);
    if (mastermap->MyGID(gid))  // in this case, the current dof is part of the interface
    {
      if ((*test_masterifdirich)[mastermap->LID(gid)] == 1.0 and
          (*masterifdirich)[mastermap->LID(gid)] != 1.0)
      {
        dserror("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | Check FS3I specific inputs                                Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::FS3I_Base::CheckFS3IInputs()
{
  // Check FS3I dynamic parameters
  DRT::Problem* problem = DRT::Problem::Instance();
  // const Teuchos::ParameterList& ioparams = problem->IOParams();
  const Teuchos::ParameterList& fs3idyn = problem->FS3IDynamicParams();
  const Teuchos::ParameterList& structdynparams = problem->StructuralDynamicParams();
  const Teuchos::ParameterList& scatradynparams = problem->ScalarTransportDynamicParams();
  // const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  const Teuchos::ParameterList& fluiddynparams = problem->FluidDynamicParams();

  // check consistency of time-integration schemes in input file
  // (including parameter theta itself in case of one-step-theta scheme)
  // and rule out unsupported versions of generalized-alpha time-integration
  // scheme (as well as other inappropriate schemes) for fluid subproblem
  INPAR::SCATRA::TimeIntegrationScheme scatratimealgo =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(
          scatradynparams, "TIMEINTEGR");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo =
      DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fluiddynparams, "TIMEINTEGR");
  INPAR::STR::DynamicType structtimealgo =
      DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdynparams, "DYNAMICTYP");

  if (fluidtimealgo == INPAR::FLUID::timeint_one_step_theta)
  {
    if (scatratimealgo != INPAR::SCATRA::timeint_one_step_theta or
        structtimealgo != INPAR::STR::dyna_onesteptheta)
      dserror(
          "Partitioned FS3I computations should feature consistent time-integration schemes for "
          "the subproblems; in this case, a one-step-theta scheme is intended to be used for the "
          "fluid subproblem, and different schemes are intended to be used for the structure "
          "and/or scalar transport subproblems!");

    if (scatradynparams.get<double>("THETA") != fluiddynparams.get<double>("THETA") or
        scatradynparams.get<double>("THETA") !=
            structdynparams.sublist("ONESTEPTHETA").get<double>("THETA"))
      dserror(
          "Parameter(s) theta for one-step-theta time-integration scheme defined in one or more of "
          "the individual fields do(es) not match for partitioned FS3I computation.");
  }
  else if (fluidtimealgo == INPAR::FLUID::timeint_afgenalpha)
  {
    if (scatratimealgo != INPAR::SCATRA::timeint_gen_alpha or
        structtimealgo != INPAR::STR::dyna_genalpha)
      dserror(
          "Partitioned FS3I computations should feature consistent time-integration schemes for "
          "the subproblems; in this case, a (alpha_f-based) generalized-alpha scheme is intended "
          "to be used for the fluid subproblem, and different schemes are intended to be used for "
          "the structure and/or scalar transport subproblems!");
  }
  else if (fluidtimealgo == INPAR::FLUID::timeint_npgenalpha)
  {
    dserror(
        "Partitioned FS3I computations do not support n+1-based generalized-alpha time-integration "
        "schemes for the fluid subproblem!");
  }
  else if (fluidtimealgo == INPAR::FLUID::timeint_bdf2 or
           fluidtimealgo == INPAR::FLUID::timeint_stationary)
  {
    dserror(
        "Partitioned FS3I computations do not support stationary of BDF2 time-integration schemes "
        "for the fluid subproblem!");
  }

  // check that incremental formulation is used for scalar transport field,
  // according to structure and fluid field
  if (scatravec_[0]->ScaTraField()->IsIncremental() == false)
    dserror("Incremental formulation required for partitioned FS3I computations!");


  // is scatra calculated conservative?
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(fs3idyn, "STRUCTSCAL_CONVFORM") ==
          INPAR::SCATRA::convform_convective and
      DRT::INPUT::IntegralValue<INPAR::FS3I::VolumeCoupling>(fs3idyn, "STRUCTSCAL_FIELDCOUPLING") ==
          INPAR::FS3I::coupling_match)
  {
    // get structure discretization
    Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
    SSI::ScatraStructureCloneStrategy clonestrategy;

    for (int i = 0; i < structdis->NumMyColElements(); ++i)
    {
      if (clonestrategy.GetImplType(structdis->lColElement(i)) !=
          INPAR::SCATRA::impltype_refconcreac)
        dserror(
            "Your scalar fields have to be calculated in conservative form, "
            "since the velocity field in the structure is NOT divergence free!");
    }
  }

  // is structure calculated dynamic when not prestressing?
  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdynparams, "DYNAMICTYP") ==
          INPAR::STR::dyna_statics and
      DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(structdynparams, "DYNAMICTYP") !=
          INPAR::STR::prestress_mulf)
    dserror(
        "Since we need a velocity field in the structure domain for the scalar field you need do "
        "calculate the structure dynamically! Exception: when prestressing..");


  // Check DESIGN SCATRA COUPLING SURF CONDITIONS
  std::vector<std::set<int>> condIDs;
  std::set<int> fluidIDs;
  std::set<int> structIDs;
  condIDs.push_back(fluidIDs);
  condIDs.push_back(structIDs);
  std::vector<std::map<int, std::vector<double>*>> PermCoeffs;
  std::map<int, std::vector<double>*> fluidcoeff;
  std::map<int, std::vector<double>*> structcoeff;
  PermCoeffs.push_back(fluidcoeff);
  PermCoeffs.push_back(structcoeff);
  const int numscal = scatravec_[0]->ScaTraField()->NumScal();

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<DRT::Discretization> disscatra = (scatravec_[i])->ScaTraField()->Discretization();
    std::vector<DRT::Condition*> coupcond;
    disscatra->GetCondition("ScaTraCoupling", coupcond);

    for (unsigned iter = 0; iter < coupcond.size(); ++iter)
    {
      int myID = (coupcond[iter])->GetInt("coupling id");
      condIDs[i].insert(myID);

      if (!infperm_)  // get all FS3I interface condition parameters from the input file
      {
        // initialize a large enough vector
        std::vector<double>* params = new std::vector<double>(7, true);
        params->at(0) = (coupcond[iter])->GetDouble("permeability coefficient");
        params->at(1) = (coupcond[iter])->GetDouble("hydraulic conductivity");
        params->at(2) = (coupcond[iter])->GetDouble("filtration coefficient");
        params->at(3) = (double)(coupcond[iter])->GetInt("wss onoff");
        const std::vector<double>* mywsscoeffs =
            (coupcond[iter])->Get<std::vector<double>>("wss coeffs");
        params->at(4) = mywsscoeffs->at(0);
        params->at(5) = mywsscoeffs->at(1);
        params->at(6) = (double)((coupcond[iter])->GetInt("numscal"));
        const std::vector<int>* onoffs = (coupcond[iter])->Get<std::vector<int>>("onoff");
        for (int k = 0; k < numscal; k++)
        {
          params->push_back((double)(onoffs->at(k)));
        }

        if (scatravec_[i]->ScaTraField()->NumScal() != params->at(6))
          dserror(
              "Number of scalars NUMSCAL in ScaTra coupling conditions with COUPID %i does not "
              "equal the number of scalars your scalar field has!",
              myID);

        if ((bool)params->at(3))  // if we have WSS depended interface permeabiliy
        {
          std::vector<DRT::Condition*> FSCCond;
          problem->GetDis("fluid")->GetCondition("FluidStressCalc", FSCCond);

          if (FSCCond.size() == 0)
            dserror(
                "If you have a WSS dependent interface permeablity you need at least one FLUID "
                "STRESS CALC CONDITION to specify the region you want to evaluate the WSS. "
                "Typically this region is equal to the SSI interface...");
        }

        PermCoeffs[i].insert(std::pair<int, std::vector<double>*>(myID, params));
        ;
      }
    }
  }

  if (condIDs[0].size() != condIDs[1].size())
    dserror("ScaTra coupling conditions need to be defined on both discretizations");

  if (!infperm_)  // now do the testing
  {
    std::map<int, std::vector<double>*> fluid_PermCoeffs = PermCoeffs[0];
    std::map<int, std::vector<double>*> struct_PermCoeffs = PermCoeffs[1];

    std::vector<int>* onoff_sum = new std::vector<int>(numscal, 0);

    for (std::map<int, std::vector<double>*>::iterator fit = fluid_PermCoeffs.begin();
         fit != fluid_PermCoeffs.end(); ++fit)  // loop over all fluid-scatra COUPIDs
    {
      const int ID = (*fit).first;
      std::vector<double>* fluid_permcoeffs =
          (*fit).second;  // get the pointer to the fluid-scatra params

      std::map<int, std::vector<double>*>::iterator sit = struct_PermCoeffs.find(
          ID);  // get corresponding structure-scatra condition with same COUPID
      std::vector<double>* structure_permcoeffs =
          (*sit).second;  // get the pointer to the structure-scatra params

      // no the actual testing
      if (fluid_permcoeffs->at(0) != structure_permcoeffs->at(0))
        dserror(
            "Permeability coefficient PERMCOEF of ScaTra couplings with COUPID %i needs to be the "
            "same!",
            ID);
      if (fluid_permcoeffs->at(1) != structure_permcoeffs->at(1))
        dserror(
            "Hydraulic conductivity coefficient CONDUCT of ScaTra couplings with COUPID %i needs "
            "to be the same!",
            ID);
      if (fluid_permcoeffs->at(2) != structure_permcoeffs->at(2))
        dserror(
            "Filtration coefficient coefficient FILTR of ScaTra couplings with COUPID %i needs to "
            "be the same!",
            ID);
      if (fluid_permcoeffs->at(2) < 0 or fluid_permcoeffs->at(2) > 1)
        dserror(
            "The filtration coefficient FILTR of ScaTra couplings with COUPID %i must be in [0;1], "
            "since it is the ratio of average pore size per area!",
            ID);
      if (fluid_permcoeffs->at(3) != structure_permcoeffs->at(3))
        dserror(
            "WSS onoff flag WSSONOFF of ScaTra couplings with COUPID %i needs to be the same!", ID);
      if (fluid_permcoeffs->at(4) != structure_permcoeffs->at(4))
        dserror(
            "First WSS coefficient WSSCOEFFS of ScaTra couplings with COUPID %i needs to be the "
            "same!",
            ID);
      if (fluid_permcoeffs->at(5) != structure_permcoeffs->at(5))
        dserror(
            "Second WSS coefficient WSSCOEFFS of ScaTra couplings with COUPID %i needs to be the "
            "same!",
            ID);
      if (fluid_permcoeffs->at(6) != structure_permcoeffs->at(6))
        dserror(
            "Number of scalars NUMSCAL of ScaTra couplings with COUPID %i needs to be the same!",
            ID);

      for (int k = 0; k < numscal; k++)
      {
        if (fluid_permcoeffs->at(7 + k) != structure_permcoeffs->at(7 + k))
          dserror("ONOFF vector of ScaTra couplings with COUPID %i needs to be the same!", ID);

        onoff_sum->at(k) += fluid_permcoeffs->at(7 + k);
      }
    }

    for (int j = 0; j < numscal; j++)
    {
      if (onoff_sum->at(j) > 1)
        dserror(
            "In the ONOFF vector the %i-th scalar has been switched on multiple times. The ON is "
            "allowed only once per scalar!",
            j);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ScatraOutput()
{
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->Output(i);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::UpdateScatraFields()
{
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->Update(i);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ScatraEvaluateSolveIterUpdate()
{
  EvaluateScatraFields();
  SetupCoupledScatraSystem();
  LinearSolveScatra();
  ScatraIterUpdate();

  // generalized-alpha time integration: compute intermediate values
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->ComputeIntermediateValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::EvaluateScatraFields()
{
  // membrane concentration at the interface needed for simplified membrane equation of Kedem and
  // Katchalsky. NOTE: needs to be set here, since it depends on the scalar interface values on both
  // discretisations changing with each Newton iteration
  SetMembraneConcentration();

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatravec_[i]->ScaTraField();

    // evaluate scatra field
    scatra->PrepareLinearSolve();
    // add contributions due to finite interface permeability
    if (!infperm_)
    {
      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<LINALG::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      scatra->SurfacePermeability(coupmat, coupforce);

      // apply Dirichlet boundary conditions to coupling matrix and vector
      Teuchos::RCP<Epetra_Vector> zeros = scatrazeros_[i];
      const Teuchos::RCP<const LINALG::MapExtractor> dbcmapex = scatra->DirichMaps();
      const Teuchos::RCP<const Epetra_Map> dbcmap = dbcmapex->CondMap();
      coupmat->ApplyDirichlet(*dbcmap, false);
      LINALG::ApplyDirichlettoSystem(coupforce, zeros, *dbcmap);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Set Membrane concentration in scatra fields              Thon 08/16 |
 *----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetMembraneConcentration() const
{
  std::vector<Teuchos::RCP<Epetra_Vector>> MembraneConc;
  ExtractMembraneConcentration(MembraneConc);

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->SetMembraneConcentration(MembraneConc[i]);
  }
}

/*----------------------------------------------------------------------*
 |  Extract membrane concentration                           thon 08/16 |
 *----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ExtractMembraneConcentration(
    std::vector<Teuchos::RCP<Epetra_Vector>>& MembraneConcentration) const
{
  //############ Fluid Field ###############
  Teuchos::RCP<Epetra_Vector> MembraneConcentration1 = CalcMembraneConcentration();
  MembraneConcentration.push_back(MembraneConcentration1);

  //############ Poro Field ###############
  // Hint: The mean concentration is not calculated again; we just map the values from the
  // Fluid-Scatra Field into the Structure-Scatra Field

  // extract interface values
  Teuchos::RCP<Epetra_Vector> interface_phin =
      scatrafieldexvec_[0]->ExtractVector(MembraneConcentration1, 1);

  // insert interface values from Fluid Field into Poro Field;
  Teuchos::RCP<Epetra_Vector> MembraneConcentration2 =
      scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(interface_phin), 1);
  MembraneConcentration.push_back(MembraneConcentration2);
}

/*----------------------------------------------------------------------*
 |  Calculate membrane concentration                         thon 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3I_Base::CalcMembraneConcentration() const
{
  // Get concentration phi2 in scatrafield2
  // Hint: in the following we talk of phi1 and phi2, but they mean the same concentration just on
  // different scatrafields
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra2 = scatravec_[1];
  Teuchos::RCP<Epetra_Vector> scatrafield2_phi2np = scatra2->ScaTraField()->Phinp();

  // extract interface values from phi2 but we are still on scatrafield2
  Teuchos::RCP<Epetra_Vector> interface2_phi2np =
      scatrafieldexvec_[1]->ExtractVector(scatrafield2_phi2np, 1);

  // insert interface values from scatrafield2 into scatrafield1; scatrafield1_phi2n is again of
  // full length, i.e. of size of scatrafield1; all values that do not belong to the interface are
  // zero
  Teuchos::RCP<Epetra_Vector> scatrafield1_phi2np =
      scatrafieldexvec_[0]->InsertVector(Scatra2ToScatra1(interface2_phi2np), 1);

  // Get concentration phi1 in scatrafield1
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra1 = scatravec_[0];
  Teuchos::RCP<Epetra_Vector> scatrafield1_phi1np = scatra1->ScaTraField()->Phinp();

  // extract interface values from phi1 but we are still on scatrafield1
  Teuchos::RCP<Epetra_Vector> interface1_phi1np =
      scatrafieldexvec_[0]->ExtractVector(scatrafield1_phi1np, 1);

  // insert interface values interface1_phi1n from scatrafield1 into the full scatrafield1 again;
  // this is just to obtain a vector whose entries are zero except for the nodes of the interface
  Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0]->InsertVector(interface1_phi1np, 1);

  // nodewise calculation of mean concentration in the interface

  for (int i = 0; i < temp->MyLength(); i++)
  {
    // here the unweighted average is uses. One could also use a logarithmic average...
    (*temp)[i] =
        0.5 *
        ((*temp)[i] +
            (*scatrafield1_phi2np)
                [i]);  // log. average:
                       // ((*temp)[i]-(*scatrafield1_phi2n)[i])/log(((*temp)[i])/((*scatrafield1_phi2n)[i]));
                       // linear approach: 0.5*((*temp)[i]+(*scatrafield1_phi2n)[i]);
  }

  // return mean concentration in the interface
  // this vector now belongs to scatrafield1!!!
  return temp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetupCoupledScatraSystem()
{
  // set up scatra rhs
  SetupCoupledScatraRHS();

  // set up scatra system matrix
  SetupCoupledScatraMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetupCoupledScatraRHS()
{
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Residual();
  SetupCoupledScatraVector(scatrarhs_, scatra1, scatra2);

  // additional contributions in case of finite interface permeability
  if (!infperm_)
  {
    Teuchos::RCP<Epetra_Vector> coup1 = scatracoupforce_[0];
    Teuchos::RCP<Epetra_Vector> coup2 = scatracoupforce_[1];

    // contribution of the same field
    scatraglobalex_->AddVector(*coup1, 0, *scatrarhs_, 1.0);
    scatraglobalex_->AddVector(*coup2, 1, *scatrarhs_, 1.0);

    // contribution of the respective other field
    Teuchos::RCP<Epetra_Vector> coup1_boundary = scatrafieldexvec_[0]->ExtractVector(coup1, 1);
    Teuchos::RCP<Epetra_Vector> temp =
        scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(coup1_boundary), 1);
    temp->Scale(-1.0);
    scatraglobalex_->AddVector(*temp, 1, *scatrarhs_);

    Teuchos::RCP<Epetra_Vector> coup2_boundary = scatrafieldexvec_[1]->ExtractVector(coup2, 1);
    temp = scatrafieldexvec_[0]->InsertVector(Scatra2ToScatra1(coup2_boundary), 1);
    temp->Scale(-1.0);
    scatraglobalex_->AddVector(*temp, 0, *scatrarhs_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetupCoupledScatraVector(Teuchos::RCP<Epetra_Vector> globalvec,
    Teuchos::RCP<const Epetra_Vector>& vec1, Teuchos::RCP<const Epetra_Vector>& vec2)
{
  if (infperm_)
  {
    // concentrations are assumed to be equal at the interface
    // extract the inner (uncoupled) dofs from second field
    Teuchos::RCP<Epetra_Vector> vec2_other = scatrafieldexvec_[1]->ExtractVector(vec2, 0);

    Teuchos::RCP<Epetra_Vector> vec2_boundary = scatrafieldexvec_[1]->ExtractVector(vec2, 1);
    Teuchos::RCP<Epetra_Vector> temp =
        scatrafieldexvec_[0]->InsertVector(Scatra2ToScatra1(vec2_boundary), 1);
    temp->Update(1.0, *vec1, 1.0);

    scatraglobalex_->InsertVector(*temp, 0, *globalvec);
    scatraglobalex_->InsertVector(*vec2_other, 1, *globalvec);
  }
  else
  {
    scatraglobalex_->InsertVector(*vec1, 0, *globalvec);
    scatraglobalex_->InsertVector(*vec2, 1, *globalvec);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetupCoupledScatraMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField()->SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField()->SystemMatrix();

  if (scatra1 == Teuchos::null) dserror("expect fluid scatra block matrix");
  if (scatra2 == Teuchos::null) dserror("expect structure scatra block matrix");

  if (infperm_)
  {
    // Uncomplete system matrix to be able to deal with slightly defective
    // interface meshes.
    scatra1->UnComplete();

    // structure scatra
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockscatra2 =
        scatra2->Split<LINALG::DefaultBlockMatrixStrategy>(
            *(scatrafieldexvec_[1]), *(scatrafieldexvec_[1]));
    blockscatra2->Complete();

    scatrasystemmatrix_->Assign(1, 1, LINALG::View, blockscatra2->Matrix(0, 0));

    (*sibtransform_)(blockscatra2->FullRowMap(), blockscatra2->FullColMap(),
        blockscatra2->Matrix(0, 1), 1.0, ADAPTER::CouplingSlaveConverter(*scatracoup_),
        scatrasystemmatrix_->Matrix(1, 0));
    (*sbitransform_)(blockscatra2->Matrix(1, 0), 1.0, ADAPTER::CouplingSlaveConverter(*scatracoup_),
        scatrasystemmatrix_->Matrix(0, 1));
    (*sbbtransform_)(blockscatra2->Matrix(1, 1), 1.0, ADAPTER::CouplingSlaveConverter(*scatracoup_),
        ADAPTER::CouplingSlaveConverter(*scatracoup_), *scatra1, true, true);

    // fluid scatra
    scatrasystemmatrix_->Assign(0, 0, LINALG::View, *scatra1);
  }
  else
  {
    // conventional contributions
    scatrasystemmatrix_->Assign(0, 0, LINALG::View, *scatra1);
    scatrasystemmatrix_->Assign(1, 1, LINALG::View, *scatra2);

    // additional contributions due to interface permeability (-> coupling terms)
    // contribution of the same field
    Teuchos::RCP<LINALG::SparseMatrix> coup1 = scatracoupmat_[0];
    Teuchos::RCP<LINALG::SparseMatrix> coup2 = scatracoupmat_[1];

    scatrasystemmatrix_->Matrix(0, 0).Add(*coup1, false, 1.0, 1.0);
    scatrasystemmatrix_->Matrix(1, 1).Add(*coup2, false, 1.0, 1.0);

    // contribution of the respective other field
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> coupblock1 =
        coup1->Split<LINALG::DefaultBlockMatrixStrategy>(
            *(scatrafieldexvec_[0]), *(scatrafieldexvec_[0]));
    coupblock1->Complete();
    (*fbitransform_)(coupblock1->Matrix(1, 1), -1.0, ADAPTER::CouplingMasterConverter(*scatracoup_),
        scatrasystemmatrix_->Matrix(1, 0));

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> coupblock2 =
        coup2->Split<LINALG::DefaultBlockMatrixStrategy>(
            *(scatrafieldexvec_[1]), *(scatrafieldexvec_[1]));
    coupblock2->Complete();
    (*sbitransform_)(coupblock2->Matrix(1, 1), -1.0, ADAPTER::CouplingSlaveConverter(*scatracoup_),
        scatrasystemmatrix_->Matrix(0, 1));
  }

  scatrasystemmatrix_->Complete();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3I_Base::Scatra2ToScatra1(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return scatracoup_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3I_Base::Scatra1ToScatra2(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return scatracoup_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::LinearSolveScatra()
{
  scatraincrement_->PutScalar(0.0);

#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<LINALG::SparseMatrix> sparse = scatrasystemmatrix_->Merge();

  scatrasolver_->Solve(sparse->EpetraMatrix(), scatraincrement_, scatrarhs_, true);
#else
  scatrasolver_->Solve(
      scatrasystemmatrix_->EpetraOperator(), scatraincrement_, scatrarhs_, true, true);
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ScatraIterUpdate()
{
  // define incremental vectors for fluid- and structure-based scatra
  // fields and extract respective vectors
  Teuchos::RCP<const Epetra_Vector> inc1;
  Teuchos::RCP<const Epetra_Vector> inc2;
  ExtractScatraFieldVectors(scatraincrement_, inc1, inc2);

  // update both fluid- and structure-based solution vectors
  scatravec_[0]->ScaTraField()->UpdateIter(inc1);
  scatravec_[1]->ScaTraField()->UpdateIter(inc2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ExtractScatraFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
    Teuchos::RCP<const Epetra_Vector>& vec1, Teuchos::RCP<const Epetra_Vector>& vec2)
{
  if (infperm_)
  {
    // process fluid scatra unknowns
    vec1 = scatraglobalex_->ExtractVector(globalvec, 0);

    // process structure scatra unknowns at the boundary
    Teuchos::RCP<Epetra_Vector> vec1_boundary = scatrafieldexvec_[0]->ExtractVector(vec1, 1);
    Teuchos::RCP<const Epetra_Vector> vec2_inner = scatraglobalex_->ExtractVector(globalvec, 1);
    Teuchos::RCP<Epetra_Vector> vec2_boundary = Scatra1ToScatra2(vec1_boundary);

    Teuchos::RCP<Epetra_Vector> vec2_temp = scatrafieldexvec_[1]->InsertVector(vec2_inner, 0);
    scatrafieldexvec_[1]->InsertVector(vec2_boundary, 1, vec2_temp);
    vec2 = vec2_temp;
  }
  else
  {
    vec1 = scatraglobalex_->ExtractVector(globalvec, 0);
    vec2 = scatraglobalex_->ExtractVector(globalvec, 1);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::CheckIsSetup()
{
  if (not IsSetup()) dserror("Setup() was not called.");
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::CheckIsInit()
{
  if (not IsInit()) dserror("Init(...) was not called.");
};
