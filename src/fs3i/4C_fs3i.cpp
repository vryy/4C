/*----------------------------------------------------------------------*/
/*! \file
\brief cpp-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I) and fluid-porous-structure-scalar-scalar
       interaction (FPS3I).

\level 2


*----------------------------------------------------------------------*/


#include "4C_fs3i.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_structure_scatra_ele.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_result_test.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_fsi_dyn.hpp"
#include "4C_fsi_free_surface_monolithic.hpp"
#include "4C_fsi_monolithicfluidsplit.hpp"
#include "4C_fsi_monolithicstructuresplit.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fs3i.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_ssi_clonestrategy.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::FS3IBase::FS3IBase()
    : infperm_(Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->FS3IDynamicParams(), "INF_PERM")),
      timemax_(Global::Problem::Instance()->FS3IDynamicParams().get<double>("MAXTIME")),
      numstep_(Global::Problem::Instance()->FS3IDynamicParams().get<int>("NUMSTEP")),
      dt_(Global::Problem::Instance()->FS3IDynamicParams().get<double>("TIMESTEP")),
      time_(0.0),
      step_(0),
      issetup_(false),
      isinit_(false)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::init()
{
  set_is_setup(false);

  scatracoup_ = Teuchos::rcp(new Core::Adapter::Coupling());
  scatraglobalex_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor());
  sbbtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform());
  sbitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform());
  sibtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform());
  fbitransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform());

  set_is_init(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::setup()
{
  check_is_init();

  set_is_setup(true);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::check_interface_dirichlet_bc()
{
  Teuchos::RCP<Core::FE::Discretization> masterdis = scatravec_[0]->ScaTraField()->discretization();
  Teuchos::RCP<Core::FE::Discretization> slavedis = scatravec_[1]->ScaTraField()->discretization();

  Teuchos::RCP<const Epetra_Map> mastermap = scatracoup_->MasterDofMap();
  Teuchos::RCP<const Epetra_Map> permmastermap = scatracoup_->PermMasterDofMap();
  Teuchos::RCP<const Epetra_Map> slavemap = scatracoup_->SlaveDofMap();
  Teuchos::RCP<const Epetra_Map> permslavemap = scatracoup_->PermSlaveDofMap();

  const Teuchos::RCP<const Core::LinAlg::MapExtractor> masterdirichmapex =
      scatravec_[0]->ScaTraField()->DirichMaps();
  const Teuchos::RCP<const Epetra_Map> masterdirichmap = masterdirichmapex->cond_map();

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

  const Teuchos::RCP<const Core::LinAlg::MapExtractor> slavedirichmapex =
      scatravec_[1]->ScaTraField()->DirichMaps();
  const Teuchos::RCP<const Epetra_Map> slavedirichmap = slavedirichmapex->cond_map();

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
  for (int i = 0; i < slavedis->dof_row_map()->NumMyElements(); ++i)
  {
    int gid = slavedis->dof_row_map()->GID(i);
    if (slavemap->MyGID(gid))  // in this case, the current dof is part of the interface
    {
      if ((*test_slaveifdirich)[slavemap->LID(gid)] == 1.0 and
          (*slaveifdirich)[slavemap->LID(gid)] != 1.0)
      {
        FOUR_C_THROW("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }

  for (int i = 0; i < masterdis->dof_row_map()->NumMyElements(); ++i)
  {
    int gid = masterdis->dof_row_map()->GID(i);
    if (mastermap->MyGID(gid))  // in this case, the current dof is part of the interface
    {
      if ((*test_masterifdirich)[mastermap->LID(gid)] == 1.0 and
          (*masterifdirich)[mastermap->LID(gid)] != 1.0)
      {
        FOUR_C_THROW("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | Check FS3I specific inputs                                Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::FS3IBase::CheckFS3IInputs()
{
  // Check FS3I dynamic parameters
  Global::Problem* problem = Global::Problem::Instance();
  // const Teuchos::ParameterList& ioparams = problem->IOParams();
  const Teuchos::ParameterList& fs3idyn = problem->FS3IDynamicParams();
  const Teuchos::ParameterList& structdynparams = problem->structural_dynamic_params();
  const Teuchos::ParameterList& scatradynparams = problem->scalar_transport_dynamic_params();
  // const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  const Teuchos::ParameterList& fluiddynparams = problem->FluidDynamicParams();

  // check consistency of time-integration schemes in input file
  // (including parameter theta itself in case of one-step-theta scheme)
  // and rule out unsupported versions of generalized-alpha time-integration
  // scheme (as well as other inappropriate schemes) for fluid subproblem
  Inpar::ScaTra::TimeIntegrationScheme scatratimealgo =
      Core::UTILS::IntegralValue<Inpar::ScaTra::TimeIntegrationScheme>(
          scatradynparams, "TIMEINTEGR");
  Inpar::FLUID::TimeIntegrationScheme fluidtimealgo =
      Core::UTILS::IntegralValue<Inpar::FLUID::TimeIntegrationScheme>(fluiddynparams, "TIMEINTEGR");
  Inpar::STR::DynamicType structtimealgo =
      Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(structdynparams, "DYNAMICTYP");

  if (fluidtimealgo == Inpar::FLUID::timeint_one_step_theta)
  {
    if (scatratimealgo != Inpar::ScaTra::timeint_one_step_theta or
        structtimealgo != Inpar::STR::dyna_onesteptheta)
      FOUR_C_THROW(
          "Partitioned FS3I computations should feature consistent time-integration schemes for "
          "the subproblems; in this case, a one-step-theta scheme is intended to be used for the "
          "fluid subproblem, and different schemes are intended to be used for the structure "
          "and/or scalar transport subproblems!");

    if (scatradynparams.get<double>("THETA") != fluiddynparams.get<double>("THETA") or
        scatradynparams.get<double>("THETA") !=
            structdynparams.sublist("ONESTEPTHETA").get<double>("THETA"))
      FOUR_C_THROW(
          "Parameter(s) theta for one-step-theta time-integration scheme defined in one or more of "
          "the individual fields do(es) not match for partitioned FS3I computation.");
  }
  else if (fluidtimealgo == Inpar::FLUID::timeint_afgenalpha)
  {
    if (scatratimealgo != Inpar::ScaTra::timeint_gen_alpha or
        structtimealgo != Inpar::STR::dyna_genalpha)
      FOUR_C_THROW(
          "Partitioned FS3I computations should feature consistent time-integration schemes for "
          "the subproblems; in this case, a (alpha_f-based) generalized-alpha scheme is intended "
          "to be used for the fluid subproblem, and different schemes are intended to be used for "
          "the structure and/or scalar transport subproblems!");
  }
  else if (fluidtimealgo == Inpar::FLUID::timeint_npgenalpha)
  {
    FOUR_C_THROW(
        "Partitioned FS3I computations do not support n+1-based generalized-alpha time-integration "
        "schemes for the fluid subproblem!");
  }
  else if (fluidtimealgo == Inpar::FLUID::timeint_bdf2 or
           fluidtimealgo == Inpar::FLUID::timeint_stationary)
  {
    FOUR_C_THROW(
        "Partitioned FS3I computations do not support stationary of BDF2 time-integration schemes "
        "for the fluid subproblem!");
  }

  // check that incremental formulation is used for scalar transport field,
  // according to structure and fluid field
  if (scatravec_[0]->ScaTraField()->IsIncremental() == false)
    FOUR_C_THROW("Incremental formulation required for partitioned FS3I computations!");


  // is scatra calculated conservative?
  if (Core::UTILS::IntegralValue<Inpar::ScaTra::ConvForm>(fs3idyn, "STRUCTSCAL_CONVFORM") ==
          Inpar::ScaTra::convform_convective and
      Core::UTILS::IntegralValue<Inpar::FS3I::VolumeCoupling>(
          fs3idyn, "STRUCTSCAL_FIELDCOUPLING") == Inpar::FS3I::coupling_match)
  {
    // get structure discretization
    Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");

    for (int i = 0; i < structdis->NumMyColElements(); ++i)
    {
      if (Adapter::GetScaTraImplType(structdis->lColElement(i)) !=
          Inpar::ScaTra::impltype_refconcreac)
        FOUR_C_THROW(
            "Your scalar fields have to be calculated in conservative form, "
            "since the velocity field in the structure is NOT divergence free!");
    }
  }
  Inpar::STR::PreStress pstype = Teuchos::getIntegralValue<Inpar::STR::PreStress>(
      Global::Problem::Instance()->structural_dynamic_params(), "PRESTRESS");
  // is structure calculated dynamic when not prestressing?
  if (Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(structdynparams, "DYNAMICTYP") ==
          Inpar::STR::dyna_statics and
      pstype != Inpar::STR::PreStress::mulf)
    FOUR_C_THROW(
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
    Teuchos::RCP<Core::FE::Discretization> disscatra =
        (scatravec_[i])->ScaTraField()->discretization();
    std::vector<Core::Conditions::Condition*> coupcond;
    disscatra->GetCondition("ScaTraCoupling", coupcond);

    for (auto& iter : coupcond)
    {
      int myID = iter->parameters().get<int>("coupling id");
      condIDs[i].insert(myID);

      if (!infperm_)  // get all FS3I interface condition parameters from the input file
      {
        // initialize a large enough vector
        auto* params = new std::vector<double>(7, true);
        params->at(0) = iter->parameters().get<double>("permeability coefficient");
        params->at(1) = iter->parameters().get<double>("hydraulic conductivity");
        params->at(2) = iter->parameters().get<double>("filtration coefficient");
        params->at(3) = (double)iter->parameters().get<int>("wss onoff");
        const auto& mywsscoeffs = iter->parameters().get<std::vector<double>>("wss coeffs");
        params->at(4) = mywsscoeffs.at(0);
        params->at(5) = mywsscoeffs.at(1);
        params->at(6) = (double)(iter->parameters().get<int>("numscal"));
        const auto& onoffs = iter->parameters().get<std::vector<int>>("onoff");
        for (int k = 0; k < numscal; k++)
        {
          params->push_back((double)(onoffs.at(k)));
        }

        if (scatravec_[i]->ScaTraField()->NumScal() != params->at(6))
          FOUR_C_THROW(
              "Number of scalars NUMSCAL in ScaTra coupling conditions with COUPID %i does not "
              "equal the number of scalars your scalar field has!",
              myID);

        if ((bool)params->at(3))  // if we have WSS depended interface permeabiliy
        {
          std::vector<Core::Conditions::Condition*> FSCCond;
          problem->GetDis("fluid")->GetCondition("FluidStressCalc", FSCCond);

          if (FSCCond.size() == 0)
            FOUR_C_THROW(
                "If you have a WSS dependent interface permeablity you need at least one FLUID "
                "STRESS CALC CONDITION to specify the region you want to evaluate the WSS. "
                "Typically this region is equal to the SSI interface...");
        }

        PermCoeffs[i].insert(std::pair<int, std::vector<double>*>(myID, params));
      }
    }
  }

  if (condIDs[0].size() != condIDs[1].size())
    FOUR_C_THROW("ScaTra coupling conditions need to be defined on both discretizations");

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
        FOUR_C_THROW(
            "Permeability coefficient PERMCOEF of ScaTra couplings with COUPID %i needs to be the "
            "same!",
            ID);
      if (fluid_permcoeffs->at(1) != structure_permcoeffs->at(1))
        FOUR_C_THROW(
            "Hydraulic conductivity coefficient CONDUCT of ScaTra couplings with COUPID %i needs "
            "to be the same!",
            ID);
      if (fluid_permcoeffs->at(2) != structure_permcoeffs->at(2))
        FOUR_C_THROW(
            "Filtration coefficient coefficient FILTR of ScaTra couplings with COUPID %i needs to "
            "be the same!",
            ID);
      if (fluid_permcoeffs->at(2) < 0 or fluid_permcoeffs->at(2) > 1)
        FOUR_C_THROW(
            "The filtration coefficient FILTR of ScaTra couplings with COUPID %i must be in [0;1], "
            "since it is the ratio of average pore size per area!",
            ID);
      if (fluid_permcoeffs->at(3) != structure_permcoeffs->at(3))
        FOUR_C_THROW(
            "WSS onoff flag WSSONOFF of ScaTra couplings with COUPID %i needs to be the same!", ID);
      if (fluid_permcoeffs->at(4) != structure_permcoeffs->at(4))
        FOUR_C_THROW(
            "First WSS coefficient WSSCOEFFS of ScaTra couplings with COUPID %i needs to be the "
            "same!",
            ID);
      if (fluid_permcoeffs->at(5) != structure_permcoeffs->at(5))
        FOUR_C_THROW(
            "Second WSS coefficient WSSCOEFFS of ScaTra couplings with COUPID %i needs to be the "
            "same!",
            ID);
      if (fluid_permcoeffs->at(6) != structure_permcoeffs->at(6))
        FOUR_C_THROW(
            "Number of scalars NUMSCAL of ScaTra couplings with COUPID %i needs to be the same!",
            ID);

      for (int k = 0; k < numscal; k++)
      {
        if (fluid_permcoeffs->at(7 + k) != structure_permcoeffs->at(7 + k))
          FOUR_C_THROW("ONOFF vector of ScaTra couplings with COUPID %i needs to be the same!", ID);

        onoff_sum->at(k) += fluid_permcoeffs->at(7 + k);
      }
    }

    for (int j = 0; j < numscal; j++)
    {
      if (onoff_sum->at(j) > 1)
        FOUR_C_THROW(
            "In the ONOFF vector the %i-th scalar has been switched on multiple times. The ON is "
            "allowed only once per scalar!",
            j);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::ScatraOutput()
{
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->check_and_write_output_and_restart();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::increment_time_and_step()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::UpdateScatraFields()
{
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->update();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::scatra_evaluate_solve_iter_update()
{
  evaluate_scatra_fields();
  setup_coupled_scatra_system();
  LinearSolveScatra();
  ScatraIterUpdate();

  // generalized-alpha time integration: compute intermediate values
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->compute_intermediate_values();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::evaluate_scatra_fields()
{
  // membrane concentration at the interface needed for simplified membrane equation of Kedem and
  // Katchalsky. NOTE: needs to be set here, since it depends on the scalar interface values on both
  // discretisations changing with each Newton iteration
  set_membrane_concentration();

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra = scatravec_[i]->ScaTraField();

    // evaluate scatra field
    scatra->PrepareLinearSolve();
    // add contributions due to finite interface permeability
    if (!infperm_)
    {
      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<Core::LinAlg::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      scatra->SurfacePermeability(coupmat, coupforce);

      // apply Dirichlet boundary conditions to coupling matrix and vector
      Teuchos::RCP<Epetra_Vector> zeros = scatrazeros_[i];
      const Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmapex = scatra->DirichMaps();
      const Teuchos::RCP<const Epetra_Map> dbcmap = dbcmapex->cond_map();
      coupmat->ApplyDirichlet(*dbcmap, false);
      Core::LinAlg::apply_dirichlet_to_system(*coupforce, *zeros, *dbcmap);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Set Membrane concentration in scatra fields              Thon 08/16 |
 *----------------------------------------------------------------------*/
void FS3I::FS3IBase::set_membrane_concentration() const
{
  std::vector<Teuchos::RCP<Epetra_Vector>> MembraneConc;
  extract_membrane_concentration(MembraneConc);

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->set_membrane_concentration(MembraneConc[i]);
  }
}

/*----------------------------------------------------------------------*
 |  Extract membrane concentration                           thon 08/16 |
 *----------------------------------------------------------------------*/
void FS3I::FS3IBase::extract_membrane_concentration(
    std::vector<Teuchos::RCP<Epetra_Vector>>& MembraneConcentration) const
{
  // ############ Fluid Field ###############
  Teuchos::RCP<Epetra_Vector> MembraneConcentration1 = calc_membrane_concentration();
  MembraneConcentration.push_back(MembraneConcentration1);

  // ############ Poro Field ###############
  //  Hint: The mean concentration is not calculated again; we just map the values from the
  //  Fluid-Scatra Field into the Structure-Scatra Field

  // extract interface values
  Teuchos::RCP<Epetra_Vector> interface_phin =
      scatrafieldexvec_[0]->extract_vector(MembraneConcentration1, 1);

  // insert interface values from Fluid Field into Poro Field;
  Teuchos::RCP<Epetra_Vector> MembraneConcentration2 =
      scatrafieldexvec_[1]->insert_vector(Scatra1ToScatra2(interface_phin), 1);
  MembraneConcentration.push_back(MembraneConcentration2);
}

/*----------------------------------------------------------------------*
 |  Calculate membrane concentration                         thon 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3IBase::calc_membrane_concentration() const
{
  // Get concentration phi2 in scatrafield2
  // Hint: in the following we talk of phi1 and phi2, but they mean the same concentration just on
  // different scatrafields
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra2 = scatravec_[1];
  Teuchos::RCP<Epetra_Vector> scatrafield2_phi2np = scatra2->ScaTraField()->Phinp();

  // extract interface values from phi2 but we are still on scatrafield2
  Teuchos::RCP<Epetra_Vector> interface2_phi2np =
      scatrafieldexvec_[1]->extract_vector(scatrafield2_phi2np, 1);

  // insert interface values from scatrafield2 into scatrafield1; scatrafield1_phi2n is again of
  // full length, i.e. of size of scatrafield1; all values that do not belong to the interface are
  // zero
  Teuchos::RCP<Epetra_Vector> scatrafield1_phi2np =
      scatrafieldexvec_[0]->insert_vector(Scatra2ToScatra1(interface2_phi2np), 1);

  // Get concentration phi1 in scatrafield1
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra1 = scatravec_[0];
  Teuchos::RCP<Epetra_Vector> scatrafield1_phi1np = scatra1->ScaTraField()->Phinp();

  // extract interface values from phi1 but we are still on scatrafield1
  Teuchos::RCP<Epetra_Vector> interface1_phi1np =
      scatrafieldexvec_[0]->extract_vector(scatrafield1_phi1np, 1);

  // insert interface values interface1_phi1n from scatrafield1 into the full scatrafield1 again;
  // this is just to obtain a vector whose entries are zero except for the nodes of the interface
  Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0]->insert_vector(interface1_phi1np, 1);

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
void FS3I::FS3IBase::setup_coupled_scatra_system()
{
  // set up scatra rhs
  setup_coupled_scatra_rhs();

  // set up scatra system matrix
  setup_coupled_scatra_matrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::setup_coupled_scatra_rhs()
{
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Residual();
  setup_coupled_scatra_vector(scatrarhs_, scatra1, scatra2);

  // additional contributions in case of finite interface permeability
  if (!infperm_)
  {
    Teuchos::RCP<Epetra_Vector> coup1 = scatracoupforce_[0];
    Teuchos::RCP<Epetra_Vector> coup2 = scatracoupforce_[1];

    // contribution of the same field
    scatraglobalex_->add_vector(*coup1, 0, *scatrarhs_, 1.0);
    scatraglobalex_->add_vector(*coup2, 1, *scatrarhs_, 1.0);

    // contribution of the respective other field
    Teuchos::RCP<Epetra_Vector> coup1_boundary = scatrafieldexvec_[0]->extract_vector(coup1, 1);
    Teuchos::RCP<Epetra_Vector> temp =
        scatrafieldexvec_[1]->insert_vector(Scatra1ToScatra2(coup1_boundary), 1);
    temp->Scale(-1.0);
    scatraglobalex_->add_vector(*temp, 1, *scatrarhs_);

    Teuchos::RCP<Epetra_Vector> coup2_boundary = scatrafieldexvec_[1]->extract_vector(coup2, 1);
    temp = scatrafieldexvec_[0]->insert_vector(Scatra2ToScatra1(coup2_boundary), 1);
    temp->Scale(-1.0);
    scatraglobalex_->add_vector(*temp, 0, *scatrarhs_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::setup_coupled_scatra_vector(Teuchos::RCP<Epetra_Vector> globalvec,
    Teuchos::RCP<const Epetra_Vector>& vec1, Teuchos::RCP<const Epetra_Vector>& vec2)
{
  if (infperm_)
  {
    // concentrations are assumed to be equal at the interface
    // extract the inner (uncoupled) dofs from second field
    Teuchos::RCP<Epetra_Vector> vec2_other = scatrafieldexvec_[1]->extract_vector(vec2, 0);

    Teuchos::RCP<Epetra_Vector> vec2_boundary = scatrafieldexvec_[1]->extract_vector(vec2, 1);
    Teuchos::RCP<Epetra_Vector> temp =
        scatrafieldexvec_[0]->insert_vector(Scatra2ToScatra1(vec2_boundary), 1);
    temp->Update(1.0, *vec1, 1.0);

    scatraglobalex_->insert_vector(*temp, 0, *globalvec);
    scatraglobalex_->insert_vector(*vec2_other, 1, *globalvec);
  }
  else
  {
    scatraglobalex_->insert_vector(*vec1, 0, *globalvec);
    scatraglobalex_->insert_vector(*vec2, 1, *globalvec);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::setup_coupled_scatra_matrix()
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField()->SystemMatrix();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField()->SystemMatrix();

  if (scatra1 == Teuchos::null) FOUR_C_THROW("expect fluid scatra block matrix");
  if (scatra2 == Teuchos::null) FOUR_C_THROW("expect structure scatra block matrix");

  if (infperm_)
  {
    // Uncomplete system matrix to be able to deal with slightly defective
    // interface meshes.
    scatra1->UnComplete();

    // structure scatra
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockscatra2 =
        scatra2->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
            *(scatrafieldexvec_[1]), *(scatrafieldexvec_[1]));
    blockscatra2->Complete();

    scatrasystemmatrix_->Assign(1, 1, Core::LinAlg::View, blockscatra2->Matrix(0, 0));

    (*sibtransform_)(blockscatra2->FullRowMap(), blockscatra2->FullColMap(),
        blockscatra2->Matrix(0, 1), 1.0, Core::Adapter::CouplingSlaveConverter(*scatracoup_),
        scatrasystemmatrix_->Matrix(1, 0));
    (*sbitransform_)(blockscatra2->Matrix(1, 0), 1.0,
        Core::Adapter::CouplingSlaveConverter(*scatracoup_), scatrasystemmatrix_->Matrix(0, 1));
    (*sbbtransform_)(blockscatra2->Matrix(1, 1), 1.0,
        Core::Adapter::CouplingSlaveConverter(*scatracoup_),
        Core::Adapter::CouplingSlaveConverter(*scatracoup_), *scatra1, true, true);

    // fluid scatra
    scatrasystemmatrix_->Assign(0, 0, Core::LinAlg::View, *scatra1);
  }
  else
  {
    // conventional contributions
    scatrasystemmatrix_->Assign(0, 0, Core::LinAlg::View, *scatra1);
    scatrasystemmatrix_->Assign(1, 1, Core::LinAlg::View, *scatra2);

    // additional contributions due to interface permeability (-> coupling terms)
    // contribution of the same field
    Teuchos::RCP<Core::LinAlg::SparseMatrix> coup1 = scatracoupmat_[0];
    Teuchos::RCP<Core::LinAlg::SparseMatrix> coup2 = scatracoupmat_[1];

    scatrasystemmatrix_->Matrix(0, 0).Add(*coup1, false, 1.0, 1.0);
    scatrasystemmatrix_->Matrix(1, 1).Add(*coup2, false, 1.0, 1.0);

    // contribution of the respective other field
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> coupblock1 =
        coup1->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
            *(scatrafieldexvec_[0]), *(scatrafieldexvec_[0]));
    coupblock1->Complete();
    (*fbitransform_)(coupblock1->Matrix(1, 1), -1.0,
        Core::Adapter::CouplingMasterConverter(*scatracoup_), scatrasystemmatrix_->Matrix(1, 0));

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> coupblock2 =
        coup2->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
            *(scatrafieldexvec_[1]), *(scatrafieldexvec_[1]));
    coupblock2->Complete();
    (*sbitransform_)(coupblock2->Matrix(1, 1), -1.0,
        Core::Adapter::CouplingSlaveConverter(*scatracoup_), scatrasystemmatrix_->Matrix(0, 1));
  }

  scatrasystemmatrix_->Complete();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3IBase::Scatra2ToScatra1(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return scatracoup_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3IBase::Scatra1ToScatra2(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return scatracoup_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::LinearSolveScatra()
{
  scatraincrement_->PutScalar(0.0);

#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse = scatrasystemmatrix_->Merge();

  scatrasolver_->Solve(sparse->EpetraMatrix(), scatraincrement_, scatrarhs_, true);
#else
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  scatrasolver_->Solve(
      scatrasystemmatrix_->EpetraOperator(), scatraincrement_, scatrarhs_, solver_params);
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::ScatraIterUpdate()
{
  // define incremental vectors for fluid- and structure-based scatra
  // fields and extract respective vectors
  Teuchos::RCP<const Epetra_Vector> inc1;
  Teuchos::RCP<const Epetra_Vector> inc2;
  extract_scatra_field_vectors(scatraincrement_, inc1, inc2);

  // update both fluid- and structure-based solution vectors
  scatravec_[0]->ScaTraField()->UpdateIter(inc1);
  scatravec_[1]->ScaTraField()->UpdateIter(inc2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::extract_scatra_field_vectors(Teuchos::RCP<const Epetra_Vector> globalvec,
    Teuchos::RCP<const Epetra_Vector>& vec1, Teuchos::RCP<const Epetra_Vector>& vec2)
{
  if (infperm_)
  {
    // process fluid scatra unknowns
    vec1 = scatraglobalex_->extract_vector(globalvec, 0);

    // process structure scatra unknowns at the boundary
    Teuchos::RCP<Epetra_Vector> vec1_boundary = scatrafieldexvec_[0]->extract_vector(vec1, 1);
    Teuchos::RCP<const Epetra_Vector> vec2_inner = scatraglobalex_->extract_vector(globalvec, 1);
    Teuchos::RCP<Epetra_Vector> vec2_boundary = Scatra1ToScatra2(vec1_boundary);

    Teuchos::RCP<Epetra_Vector> vec2_temp = scatrafieldexvec_[1]->insert_vector(vec2_inner, 0);
    scatrafieldexvec_[1]->insert_vector(vec2_boundary, 1, vec2_temp);
    vec2 = vec2_temp;
  }
  else
  {
    vec1 = scatraglobalex_->extract_vector(globalvec, 0);
    vec2 = scatraglobalex_->extract_vector(globalvec, 1);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::check_is_setup()
{
  if (not is_setup()) FOUR_C_THROW("setup() was not called.");
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3IBase::check_is_init()
{
  if (not is_init()) FOUR_C_THROW("init(...) was not called.");
};

FOUR_C_NAMESPACE_CLOSE
