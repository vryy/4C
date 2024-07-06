/*----------------------------------------------------------------------*/
/*! \file

\brief FPSI Coupling Object: Holds all objects on the Fluid-Poro-Interface and evaluates the
Fluid-Poro-Coupling Matrixes!


\level 3
*/
///*----------------------------------------------------------------------*/
// GENERAL includes
#include "4C_fpsi_coupling.hpp"

#include "4C_adapter_ale_fpsi.hpp"
#include "4C_adapter_fld_fluid.hpp"
#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_structure_aux.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::FpsiCoupling::FpsiCoupling(Teuchos::RCP<PoroElast::Monolithic> poro,
    Teuchos::RCP<Adapter::Fluid> fluid, Teuchos::RCP<Adapter::AleFpsiWrapper> ale,
    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap,
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap)
    : poro_(poro),
      fluid_(fluid),
      ale_(ale),
      fluidvelpres_extractor_(Teuchos::null),
      fluidvel_extractor_(Teuchos::null),
      porofluid_extractor_(Teuchos::null),
      porostruct_extractor_(Teuchos::null),
      poro_extractor_(Teuchos::null),
      fluid_fsifpsi_extractor_(Teuchos::null),
      isfirstcall_(true),
      fluid_poro_fluid_interface_map_(Fluid_PoroFluid_InterfaceMap),  // to be removed later
      poro_fluid_fluid_interface_map_(PoroFluid_Fluid_InterfaceMap),
      conductivity_(0.0)
{
  setup_interface_coupling();
  init_coupling_matrixes_rhs();
  re_init_coupling_matrix_transform();
  return;
}
/*----------------------------------------------------------------------/
| Initialize Coupling Matrixes and Coupling RHS              ager 12/14 |
/----------------------------------------------------------------------*/
void FPSI::FpsiCoupling::init_coupling_matrixes_rhs()
{
  // fluid extractor
  Core::LinAlg::MapExtractor fluidextractor(*fluid_->dof_row_map(), fluid_->dof_row_map(), false);
  // ale extractor
  Core::LinAlg::MapExtractor aleextractor(
      *(ale_->interface()->other_map()), ale_->interface()->other_map(), false);

  c_pp_ = Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>(
      new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *poro_->extractor(), *poro_->extractor(), 81, true, true));
  c_ff_ = Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>(
      new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *fluid_->interface(), *fluid_->interface(), 81, true, true));
  c_pf_ = Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>(
      new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          fluidextractor, *poro_->extractor(), 81, true, true));
  c_fp_ = Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>(
      new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *poro_->extractor(), fluidextractor, 81, true, true));
  c_pa_ = Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>(
      new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          aleextractor, *poro_->extractor(), 81, true, true));
  c_fa_ = Teuchos::RCP<Core::LinAlg::SparseMatrix>(
      new Core::LinAlg::SparseMatrix(*fluid_->dof_row_map(), 81, true, true));

  c_rhs_s_ = Teuchos::RCP<Epetra_Vector>(
      new Epetra_Vector(*poro_->structure_field()->dof_row_map(), true));
  c_rhs_pf_ =
      Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*poro_->fluid_field()->dof_row_map(), true));
  c_rhs_f_ = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*fluid_->dof_row_map(), true));

  return;
}

/*----------------------------------------------------------------------/
| Setup the Coupling Object                                  ager 12/14 |
/----------------------------------------------------------------------*/
void FPSI::FpsiCoupling::setup_interface_coupling()
{
  const int ndim = Global::Problem::instance()->n_dim();

  Teuchos::RCP<Core::FE::Discretization> fluiddis = fluid_field()->discretization();
  Teuchos::RCP<Core::FE::Discretization> porofluiddis =
      poro_field()->fluid_field()->discretization();
  Teuchos::RCP<Core::FE::Discretization> porostructdis =
      poro_field()->structure_field()->discretization();

  {
    porofluid_extractor_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.add_selector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*porofluiddis, "fpsi_coupling", 0, ndim + 1)));
    mcs.setup_extractor(*porofluiddis, *(porofluiddis->dof_row_map()), *porofluid_extractor_);
  }

  {
    porostruct_extractor_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.add_selector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*porostructdis, "fpsi_coupling", 0, ndim)));
    mcs.setup_extractor(*porostructdis, *(porostructdis->dof_row_map()), *porostruct_extractor_);
  }

  {
    fluidvelpres_extractor_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.add_selector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*fluiddis, "fpsi_coupling", 0, ndim + 1)));
    mcs.setup_extractor(*fluiddis, *(fluiddis->dof_row_map()), *fluidvelpres_extractor_);
  }

  {
    fluidvel_extractor_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.add_selector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*fluiddis, "fpsi_coupling", 0, ndim)));
    mcs.setup_extractor(*fluiddis, *(fluiddis->dof_row_map()), *fluidvel_extractor_);
  }

  {
    fluid_fsifpsi_extractor_ = Teuchos::rcp(new FPSI::UTILS::MapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.add_selector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*fluiddis, "FSICoupling", 0, ndim)));
    mcs.add_selector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*fluiddis, "fpsi_coupling", 0, ndim)));
    mcs.setup_extractor(*fluiddis, *(fluiddis->dof_row_map()), *fluid_fsifpsi_extractor_);
  }

  {
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    // Split poro_field into:
    //                      --> Structure (inside + FSI-Interface)
    //                      --> Structure FPSI-Interface
    //                      --> PoroFluid (inside + FSI-Interface)
    //                      --> PoroFluid FPSI-Interface

    Teuchos::RCP<const Epetra_Map> s_other_map = Core::LinAlg::MergeMap(
        poro_field()->structure_field()->interface()->Map(Solid::MapExtractor::cond_other),
        poro_field()->structure_field()->interface()->Map(Solid::MapExtractor::cond_fsi));
    vecSpaces.push_back(s_other_map);  // other map
    vecSpaces.push_back(poro_field()->structure_field()->interface()->Map(
        Solid::MapExtractor::cond_fpsi));                    // fpsi_coupling
    vecSpaces.push_back(porofluid_extractor_->other_map());  // other map
    vecSpaces.push_back(porofluid_extractor_->cond_map());   // fpsi_coupling

    Teuchos::RCP<Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::merge_maps(vecSpaces);
    // full Poroelasticity-blockmap
    poro_extractor_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor());
    poro_extractor_->setup(*fullmap, vecSpaces);
  }

  // porous fluid to fluid
  icoup_pf_f_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoup_pf_f_->setup_condition_coupling(*porofluiddis, porofluid_extractor_->cond_map(), *fluiddis,
      fluidvelpres_extractor_->cond_map(), "fpsi_coupling", ndim + 1, false);

  // porous structure to fluid
  icoup_ps_f_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoup_ps_f_->setup_condition_coupling(*porostructdis, porostruct_extractor_->cond_map(),
      *fluiddis, fluidvel_extractor_->cond_map(), "fpsi_coupling", ndim, false);

  // porous structure to ale
  icoup_ps_a_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoup_ps_a_->setup_condition_coupling(*porostructdis, porostruct_extractor_->cond_map(),
      *ale_field()->discretization(), ale_field()->interface()->fpsi_cond_map(), "fpsi_coupling",
      ndim, false);

  return;
}

/*-----------------------------------------------------------------------/
| Method reinitializes the matrix transformation objects      ager 12/14 |
/-----------------------------------------------------------------------*/
void FPSI::FpsiCoupling::re_init_coupling_matrix_transform()
{
  // create transformation objects for coupling terms
  couplingrowtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  couplingrowtransform2_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  couplingrowtransform3_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  couplingrowtransform4_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  couplingrowtransform5_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  couplingcoltransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  couplingcoltransform2_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  couplingrowcoltransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
  couplingrowcoltransform2_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
}

/*-------------------------------------------------------------------------------/
| Evaluate Coupling Matrixes and Coupling RHS    orig. rauch / modif. ager 12/14 |
/-------------------------------------------------------------------------------*/
void FPSI::FpsiCoupling::evaluate_coupling_matrixes_rhs()
{
  // Evaluates all Coupling Matrixes ...
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::fpsi_coupling::evaluate_coupling_matrixes_rhs");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_fp_porofluid =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *(poro_field()->fluid_field()->dof_row_map()), 81, true, true));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_pf_porofluid =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(fluid_field()->dof_row_map()), 81, true, true));

  // set all coupling matrixes to zero!!
  c_pp_->zero();
  c_ff_->zero();
  c_pf_->zero();
  c_fp_->zero();
  c_fa_->zero();
  c_pa_->zero();

  c_rhs_s_->PutScalar(0.0);
  c_rhs_pf_->PutScalar(0.0);
  c_rhs_f_->PutScalar(0.0);

  k_pf_porofluid->zero();

  const Core::Adapter::Coupling& couppff_fpsi = *icoup_pf_f_;
  const Core::Adapter::Coupling& coupsf_fpsi = *icoup_ps_f_;
  const Core::Adapter::Coupling& coup_ps_a_fpsi = *icoup_ps_a_;

  Global::Problem* problem = Global::Problem::instance();
  const Teuchos::ParameterList& fpsidynparams = problem->fpsi_dynamic_params();
  Inpar::FPSI::PartitionedCouplingMethod method =
      Core::UTILS::IntegralValue<Inpar::FPSI::PartitionedCouplingMethod>(
          fpsidynparams, "PARTITIONED");

  if (method != Inpar::FPSI::nocoupling)
  {
    // set general vector values needed by elements

    poro_field()->fluid_field()->discretization()->clear_state();

    poro_field()->fluid_field()->discretization()->set_state(
        0, "dispnp", poro_field()->fluid_field()->dispnp());

    poro_field()->fluid_field()->discretization()->set_state(
        0, "gridv", poro_field()->fluid_field()->grid_vel());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "dispn", poro_field()->fluid_field()->dispn());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "veln", poro_field()->fluid_field()->veln());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "velaf", poro_field()->fluid_field()->velnp());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "velnp", poro_field()->fluid_field()->velnp());

    fluid_field()->discretization()->clear_state();

    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->dispnp());
    fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->grid_vel());
    fluid_field()->discretization()->set_state(0, "dispn", fluid_field()->dispn());
    fluid_field()->discretization()->set_state(0, "veln", fluid_field()->veln());
    fluid_field()->discretization()->set_state(0, "velaf", fluid_field()->velnp());
    fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->velnp());

    // create the parameters for the discretization
    Teuchos::ParameterList fparams;

    // action for elements
    fparams.set<int>("action", FLD::fpsi_coupling);
    fparams.set("timescale", poro_field()->fluid_field()->residual_scaling());

    fparams.set("dt", fpsidynparams.get<double>("TIMESTEP"));
    fparams.set<int>("Physical Type", poro_field()->fluid_field()->physical_type());

    if (method == Inpar::FPSI::monolithic)
    {
      fparams.set<std::string>("fillblock", "Porofluid_Freefluid");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
          0,                                       // fluiddofset for column
          k_pf_porofluid,                          // coupling matrix with fluid rowmap
          Teuchos::null,                           // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      // what's the current problem type? Is it a fps3i problem?
      Core::ProblemType probtype = Global::Problem::instance()->get_problem_type();

      if (probtype == Core::ProblemType::fps3i)
      {
        if (conductivity_ == 0.0)
        {
          FOUR_C_THROW(
              "In the case of FPS3I, a positive conductivity must be set in DESIGN SCATRA COUPLING "
              "SURF CONDITIONS");
        }
        fparams.set("membrane conductivity", conductivity_);
      }

      fluid_field()->discretization()->evaluate_condition(fparams, fluidstrategy, "fpsi_coupling");
      k_pf_porofluid->complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform_)(*k_pf_porofluid, 1.0,
            Core::Adapter::CouplingSlaveConverter(couppff_fpsi), c_pf_->matrix(1, 0), true);
      }

      fparams.set<std::string>("fillblock", "Porofluid_Structure");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->zero();
      Core::FE::AssembleStrategy fluidstrategy21(0,  // porofluiddofset for row
          0,                                         // porofluiddofset for column
          k_pf_porofluid,                            // coupling matrix with fluid rowmap
          Teuchos::null,                             // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, fluidstrategy21, "fpsi_coupling");

      k_pf_porofluid->complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());
      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowcoltransform2_)(*k_pf_porofluid, 1.0,
            Core::Adapter::CouplingSlaveConverter(
                couppff_fpsi),  // row converter: important to use slave converter
            Core::Adapter::CouplingSlaveConverter(
                coupsf_fpsi),  //  col converter: important to use slave converter
            c_pp_->matrix(1, 0),
            false,  // bool exactmatch = true (default)
            true);
      }

      fparams.set<std::string>("fillblock", "Fluid_Porofluid");
      fparams.set("InterfaceFacingElementMap", poro_fluid_fluid_interface_map_);
      k_fp_porofluid->zero();
      Core::FE::AssembleStrategy porofluidstrategy(0,  // porofluiddofset for row
          0,                                           // porofluiddofset for column
          k_fp_porofluid,                              // porofluid-structure coupling matrix
          Teuchos::null,                               // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      poro_field()->fluid_field()->discretization()->evaluate_condition(
          fparams, porofluidstrategy, "fpsi_coupling");
      k_fp_porofluid->complete(poro_field()->fluid_domain_map(), poro_field()->fluid_range_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform2_)(*k_fp_porofluid, 1.0,
            Core::Adapter::CouplingMasterConverter(couppff_fpsi), c_fp_->matrix(0, 1),
            true);  // add
      }

      fparams.set<std::string>("fillblock", "Fluid_Structure");
      fparams.set("InterfaceFacingElementMap", poro_fluid_fluid_interface_map_);

      // move me somewhere else
      Teuchos::RCP<Core::LinAlg::SparseMatrix> k_pfs_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *(poro_field()->fluid_field()->discretization()->dof_row_map()), 81, true, true));

      k_pfs_->un_complete();

      Core::FE::AssembleStrategy structurestrategy(0,  // porofluiddofset for row
          1,                                           // structuredofset for column
          k_pfs_,                                      // coupling matrix with porofluid rowmap
          Teuchos::null,                               // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      poro_field()->fluid_field()->discretization()->evaluate_condition(
          fparams, structurestrategy, "fpsi_coupling");

      k_pfs_->complete(poro_field()->structure_domain_map(), poro_field()->fluid_range_map());
      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform3_)(*k_pfs_, 1.0,
            Core::Adapter::CouplingMasterConverter(couppff_fpsi), c_fp_->matrix(0, 0),
            true);  // add
      }

      ///// Fluid_Structure (fluid part / linearization of tangentials with respect to
      /// displacements)
      k_pf_porofluid->reset();
      k_pf_porofluid->un_complete();

      Core::FE::AssembleStrategy structurestrategy2(0,  // porofluiddofset for row
          0,                                            // fluiddofset for column
          k_pf_porofluid,                               // coupling matrix with porofluid rowmap
          Teuchos::null,                                // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      fluid_field()->discretization()->evaluate_condition(
          fparams, structurestrategy2, "fpsi_coupling");

      k_pf_porofluid->complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingcoltransform_)(fluid_field()->block_system_matrix()->full_row_map(),
            fluid_field()->block_system_matrix()->full_col_map(), *k_pf_porofluid, 1.0,
            Core::Adapter::CouplingSlaveConverter(
                coupsf_fpsi),  // row converter: important to use slave converter
            c_fp_->matrix(0, 0),
            false,  // bool exactmatch = true (default)
            true);
      }

      fparams.set<std::string>("fillblock", "Fluid_Fluid");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->zero();
      k_pf_porofluid->un_complete();

      Core::FE::AssembleStrategy fluidfluidstrategy(0,  // fluiddofset for row
          0,                                            // fluiddofset for column
          c_ff_,                                        // porofluid-structure coupling matrix
          Teuchos::null,                                // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, fluidfluidstrategy, "fpsi_coupling");

      fparams.set<std::string>("fillblock", "Structure_Fluid");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->zero();
      k_pf_porofluid->un_complete();

      Core::FE::AssembleStrategy structurefluidstrategy(0,  // fluid dofset for row
          0,                                                // fluid dofset for column
          k_pf_porofluid,                                   // coupling matrix with fluid rowmap
          Teuchos::null,                                    // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, structurefluidstrategy, "fpsi_coupling");

      k_pf_porofluid->complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());
      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform4_)(
            *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(k_pf_porofluid), 1.0,
            Core::Adapter::CouplingSlaveConverter(coupsf_fpsi),  // important to use slave converter
            c_pf_->matrix(0, 0),
            true);  // add
      }

      fparams.set<std::string>("fillblock", "Structure_Structure");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->zero();
      k_pf_porofluid->un_complete();

      Core::FE::AssembleStrategy structurestructurestrategy(0,  // fluid dofset for row
          0,                                                    // fluid dofset for column
          k_pf_porofluid,                                       // coupling matrix with fluid rowmap
          Teuchos::null,                                        // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, structurestructurestrategy, "fpsi_coupling");
      // condense linearization with respect to the ale mesh motion (interface structural
      // displacements = interface ale displacements)
      k_pf_porofluid->complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowcoltransform_)(*k_pf_porofluid, 1.0,
            Core::Adapter::CouplingSlaveConverter(
                coupsf_fpsi),  // row converter: important to use slave converter
            Core::Adapter::CouplingSlaveConverter(
                coupsf_fpsi),  // col converter: important to use slave converter
            c_pp_->matrix(0, 0),
            false,  // bool exactmatch = true (default)
            true);  // add
      }

      // Process inner ale dofs
      fparams.set<std::string>("fillblock", "Structure_Ale");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);

      // temporal matrix
      // todo (initialization should be avoided in every iteration...)
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>
          temp6 = Teuchos::rcp(
              new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                  *ale_field()->interface(), *fluid_field()->interface(), 81, false,
                  false));  // Use fluid_field()->Interface =

      // assemble into fluid row and column dofs -> need to transform rows to structure dofs and
      // cols to ale dofs
      Core::FE::AssembleStrategy structurealestrategy(0,  // fluid dofset for row
          1,                                              // ale dofset for column
          temp6,                                          // coupling matrix with fluid rowmap
          Teuchos::null,                                  // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      // evaluate coupling terms
      fluid_field()->discretization()->evaluate_condition(
          fparams, structurealestrategy, "fpsi_coupling");
      temp6->complete();  // for row transform!

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform5_)(temp6->matrix(FLD::UTILS::MapExtractor::cond_other,
                                      ALE::UTILS::MapExtractor::cond_other),
            1.0,
            Core::Adapter::CouplingSlaveConverter(coupsf_fpsi),  // important to use slave converter
            c_pa_->matrix(0, 0), false);
      }
    }  // if monolithic
    else
    {
      if (isfirstcall_)
        std::cout << "LINEARIZATION OF FPSI INTERFACE IS TURNED OFF \n"
                     "IF YOU THINK THAT SUCKS. SET 'PARITIONED' to 'monolithic' \n"
                     "IN THE FPSI SECTION OF YOUR INPUT FILE !!! "
                  << std::endl;
      isfirstcall_ = false;
    }

    //////////////////////////////////////
    ///////      __       ___       //////
    ///////      |_|  |_| |__       //////
    ///////      | \  | |  __|      //////
    //////                          //////
    //////////////////////////////////////

    Teuchos::RCP<Epetra_Vector> temprhs = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> temprhs2 = Teuchos::null;

    fparams.set<std::string>("fillblock", "conti");
    fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
    temprhs2 = Teuchos::rcp(new Epetra_Vector(*poro_field()->dof_row_map(), true));
    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    Core::FE::AssembleStrategy rhscontistrategy(0,  // fluid dofset for row
        0,                                          // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    fluid_field()->discretization()->evaluate_condition(fparams, rhscontistrategy, "fpsi_coupling");

    // extract FPSI part of the fluid field
    temprhs = fluidvelpres_extractor_->extract_cond_vector(temprhs);

    // replace global fluid interface dofs through porofluid interface dofs
    temprhs = i_fluid_to_porofluid(temprhs);

    // insert porofluid interface entries into vector with full porofield length (0: inner dofs of
    // structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of
    // porofluid )
    porofluid_extractor_->insert_cond_vector(temprhs, c_rhs_pf_);

    // add vector with full porofield length to global rhs

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    fparams.set<std::string>("fillblock", "structure");
    fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
    temprhs2 = Teuchos::rcp(new Epetra_Vector(*poro_field()->dof_row_map(), true));

    Core::FE::AssembleStrategy rhsstructurestrategy(0,  // fluid dofset for row
        0,                                              // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    fluid_field()->discretization()->evaluate_condition(
        fparams, rhsstructurestrategy, "fpsi_coupling");

    // extract FPSI part of the fluid field
    temprhs = fluidvel_extractor_->extract_cond_vector(temprhs);  //
    // replace global fluid interface dofs through porofluid interface dofs
    temprhs = i_fluid_to_porostruct(temprhs);  //
    // insert porofluid interface
    porostruct_extractor_->add_cond_vector(temprhs, c_rhs_s_);

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    fparams.set<std::string>("fillblock", "fluid");
    fparams.set("InterfaceFacingElementMap", poro_fluid_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*poro_field()->fluid_field()->dof_row_map(0), true));
    temprhs2 = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(0), true));

    Core::FE::AssembleStrategy rhsfluidstrategy(0,  // fluid dofset for row
        0,                                          // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    poro_field()->fluid_field()->discretization()->evaluate_condition(
        fparams, rhsfluidstrategy, "fpsi_coupling");
    // extract FPSI part of the poro fluid field
    temprhs = porofluid_extractor_->extract_cond_vector(temprhs);  //

    // replace global fluid interface dofs through porofluid interface dofs
    temprhs = i_porofluid_to_fluid(temprhs);
    // insert porofluid interface entries into vector with full fluidfield length
    fluidvelpres_extractor_->insert_cond_vector(temprhs, temprhs2);
    // add vector with full porofield length to global rhs
    c_rhs_f_->Update(1.0, *temprhs2, 0.0);

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    fparams.set<std::string>("fillblock", "fluidfluid");  // (wot,tangentialfac*uot) part
    fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(0), true));

    Core::FE::AssembleStrategy rhsfluidfluidstrategy(0,  // fluid dofset for row
        0,                                               // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    fluid_field()->discretization()->evaluate_condition(
        fparams, rhsfluidfluidstrategy, "fpsi_coupling");

    c_rhs_f_->Update(1.0, *temprhs, 1.0);

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    //////////////////////////////////////
    //////                          //////
    //////   NEUMANN INTEGRATION    //////
    //////                          //////
    //////////////////////////////////////
    if (method == Inpar::FPSI::monolithic or method == Inpar::FPSI::RobinNeumann)
    {
      fparams.set<std::string>("fillblock", "NeumannIntegration");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);

      Core::FE::AssembleStrategy rhsfluidfluidstrategy2(0,  // fluid dofset for row
          0,                                                // fluid dofset for column
          c_ff_,                                            // coupling matrix with fluid rowmap
          Teuchos::null,
          temprhs,  // rhs vector
          Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, rhsfluidfluidstrategy2, "NeumannIntegration");

      c_rhs_f_->Update(1.0, *temprhs, 1.0);
      temprhs->PutScalar(0.0);

      {
        fparams.set<std::string>("fillblock", "NeumannIntegration_Ale");
        fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);

        Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> tmp_c_fa = Teuchos::rcp(
            new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                *ale_field()->interface(), *fluid_field()->fpsi_interface(), 81, false, true));

        Core::FE::AssembleStrategy rhsfluidfluidstrategy3(0,  // fluid dofset for row
            1,                                                // ale dofset for column
            tmp_c_fa,                                         // coupling matrix with fluid rowmap
            Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

        fluid_field()->discretization()->evaluate_condition(
            fparams, rhsfluidfluidstrategy3, "NeumannIntegration");
        tmp_c_fa->complete();

        // Add all inner ale parts to c_fa_ directly
        c_fa_->add(tmp_c_fa->matrix(
                       FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_other),
            false, 1.0, 0.0);
        c_fa_->add(tmp_c_fa->matrix(
                       FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_other),
            false, 1.0, 1.0);

        //-->now transform ale fpsi block to structure (is condensed)!!!
        // still in ale domain map!!!
        Core::LinAlg::SparseMatrix tmp_c_fp =
            Core::LinAlg::SparseMatrix(*fluid_field()->dof_row_map(), 81);

        // Add all condensed parts to tmp_c_fa...
        tmp_c_fp.add(tmp_c_fa->matrix(
                         FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_fpsi),
            false, 1.0, 0.0);
        tmp_c_fp.add(tmp_c_fa->matrix(
                         FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_fpsi),
            false, 1.0, 1.0);
        tmp_c_fp.add(tmp_c_fa->matrix(
                         FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_fsi),
            false, 1.0, 1.0);
        tmp_c_fp.add(tmp_c_fa->matrix(
                         FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_fsi),
            false, 1.0, 1.0);
        tmp_c_fp.complete(
            *ale_field()->interface()->fpsi_cond_map(), *fluid_field()->dof_row_map());

        // For Ale Condensation ==> AleColumns to StructuralColumns
        (*couplingcoltransform2_)(ale_field()->block_system_matrix()->full_row_map(),
            ale_field()->block_system_matrix()->full_col_map(), tmp_c_fp, 1.0,
            Core::Adapter::CouplingSlaveConverter(
                coup_ps_a_fpsi),  // row converter: important to use slave converter
            c_fp_->matrix(0, 0),
            false,  // bool exactmatch = true (default)
            true);  // bool add = false (default)
      }
    }
    else  // only fill rhs
    {
      fparams.set<std::string>("fillblock", "NeumannIntegration");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      temprhs->PutScalar(0.0);
      temprhs2->PutScalar(0.0);
      fluid_field()->discretization()->evaluate_condition(
          fparams, rhsfluidfluidstrategy, "NeumannIntegration");

      c_rhs_f_->Update(1.0, *temprhs, 1.0);
    }

    ////////////////////////////
    // DONE -> CLEAR STATES
    poro_field()->fluid_field()->discretization()->clear_state();
    fluid_field()->discretization()->clear_state();

  }  // if not nocoupling

  ////////////////////////////
  // DONE -> Complete Coupling Matrixes
  c_ff_->complete();
  c_pp_->complete();
}

/*----------------------------------------------------------------------/
| set hydraulic conductivity
/----------------------------------------------------------------------*/
void FPSI::FpsiCoupling::set_conductivity(double conduct) { conductivity_ = conduct; }

FOUR_C_NAMESPACE_CLOSE
