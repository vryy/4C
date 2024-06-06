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
#include "4C_discretization_condition_selector.hpp"
#include "4C_discretization_fem_general_assemblestrategy.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_lib_discret.hpp"
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
FPSI::FpsiCoupling::FpsiCoupling(Teuchos::RCP<POROELAST::Monolithic> poro,
    Teuchos::RCP<ADAPTER::Fluid> fluid, Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale,
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
  CORE::LINALG::MapExtractor fluidextractor(*fluid_->dof_row_map(), fluid_->dof_row_map(), false);
  // ale extractor
  CORE::LINALG::MapExtractor aleextractor(
      *(ale_->Interface()->OtherMap()), ale_->Interface()->OtherMap(), false);

  c_pp_ = Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>(
      new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *poro_->Extractor(), *poro_->Extractor(), 81, true, true));
  c_ff_ = Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>(
      new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *fluid_->Interface(), *fluid_->Interface(), 81, true, true));
  c_pf_ = Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>(
      new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          fluidextractor, *poro_->Extractor(), 81, true, true));
  c_fp_ = Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>(
      new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *poro_->Extractor(), fluidextractor, 81, true, true));
  c_pa_ = Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>(
      new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          aleextractor, *poro_->Extractor(), 81, true, true));
  c_fa_ = Teuchos::RCP<CORE::LINALG::SparseMatrix>(
      new CORE::LINALG::SparseMatrix(*fluid_->dof_row_map(), 81, true, true));

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
  const int ndim = GLOBAL::Problem::Instance()->NDim();

  Teuchos::RCP<DRT::Discretization> fluiddis = fluid_field()->discretization();
  Teuchos::RCP<DRT::Discretization> porofluiddis = poro_field()->fluid_field()->discretization();
  Teuchos::RCP<DRT::Discretization> porostructdis =
      poro_field()->structure_field()->discretization();

  {
    porofluid_extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor());
    CORE::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*porofluiddis, "fpsi_coupling", 0, ndim + 1)));
    mcs.SetupExtractor(*porofluiddis, *(porofluiddis->dof_row_map()), *porofluid_extractor_);
  }

  {
    porostruct_extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor());
    CORE::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*porostructdis, "fpsi_coupling", 0, ndim)));
    mcs.SetupExtractor(*porostructdis, *(porostructdis->dof_row_map()), *porostruct_extractor_);
  }

  {
    fluidvelpres_extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor());
    CORE::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*fluiddis, "fpsi_coupling", 0, ndim + 1)));
    mcs.SetupExtractor(*fluiddis, *(fluiddis->dof_row_map()), *fluidvelpres_extractor_);
  }

  {
    fluidvel_extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor());
    CORE::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*fluiddis, "fpsi_coupling", 0, ndim)));
    mcs.SetupExtractor(*fluiddis, *(fluiddis->dof_row_map()), *fluidvel_extractor_);
  }

  {
    fluid_fsifpsi_extractor_ = Teuchos::rcp(new FPSI::UTILS::MapExtractor());
    CORE::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*fluiddis, "FSICoupling", 0, ndim)));
    mcs.AddSelector(Teuchos::rcp(
        new CORE::Conditions::NDimConditionSelector(*fluiddis, "fpsi_coupling", 0, ndim)));
    mcs.SetupExtractor(*fluiddis, *(fluiddis->dof_row_map()), *fluid_fsifpsi_extractor_);
  }

  {
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    // Split poro_field into:
    //                      --> Structure (inside + FSI-Interface)
    //                      --> Structure FPSI-Interface
    //                      --> PoroFluid (inside + FSI-Interface)
    //                      --> PoroFluid FPSI-Interface

    Teuchos::RCP<const Epetra_Map> s_other_map = CORE::LINALG::MergeMap(
        poro_field()->structure_field()->Interface()->Map(STR::MapExtractor::cond_other),
        poro_field()->structure_field()->Interface()->Map(STR::MapExtractor::cond_fsi));
    vecSpaces.push_back(s_other_map);  // other map
    vecSpaces.push_back(poro_field()->structure_field()->Interface()->Map(
        STR::MapExtractor::cond_fpsi));                     // fpsi_coupling
    vecSpaces.push_back(porofluid_extractor_->OtherMap());  // other map
    vecSpaces.push_back(porofluid_extractor_->CondMap());   // fpsi_coupling

    Teuchos::RCP<Epetra_Map> fullmap = CORE::LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    poro_extractor_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
    poro_extractor_->Setup(*fullmap, vecSpaces);
  }

  // porous fluid to fluid
  icoup_pf_f_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  icoup_pf_f_->setup_condition_coupling(*porofluiddis, porofluid_extractor_->CondMap(), *fluiddis,
      fluidvelpres_extractor_->CondMap(), "fpsi_coupling", ndim + 1, false);

  // porous structure to fluid
  icoup_ps_f_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  icoup_ps_f_->setup_condition_coupling(*porostructdis, porostruct_extractor_->CondMap(), *fluiddis,
      fluidvel_extractor_->CondMap(), "fpsi_coupling", ndim, false);

  // porous structure to ale
  icoup_ps_a_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  icoup_ps_a_->setup_condition_coupling(*porostructdis, porostruct_extractor_->CondMap(),
      *ale_field()->discretization(), ale_field()->Interface()->FPSICondMap(), "fpsi_coupling",
      ndim, false);

  return;
}

/*-----------------------------------------------------------------------/
| Method reinitializes the matrix transformation objects      ager 12/14 |
/-----------------------------------------------------------------------*/
void FPSI::FpsiCoupling::re_init_coupling_matrix_transform()
{
  // create transformation objects for coupling terms
  couplingrowtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  couplingrowtransform2_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  couplingrowtransform3_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  couplingrowtransform4_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  couplingrowtransform5_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  couplingcoltransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  couplingcoltransform2_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  couplingrowcoltransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  couplingrowcoltransform2_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
}

/*-------------------------------------------------------------------------------/
| Evaluate Coupling Matrixes and Coupling RHS    orig. rauch / modif. ager 12/14 |
/-------------------------------------------------------------------------------*/
void FPSI::FpsiCoupling::evaluate_coupling_matrixes_rhs()
{
  // Evaluates all Coupling Matrixes ...
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::fpsi_coupling::evaluate_coupling_matrixes_rhs");

  Teuchos::RCP<CORE::LINALG::SparseMatrix> k_fp_porofluid =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *(poro_field()->fluid_field()->dof_row_map()), 81, true, true));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> k_pf_porofluid =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*(fluid_field()->dof_row_map()), 81, true, true));

  // set all coupling matrixes to zero!!
  c_pp_->Zero();
  c_ff_->Zero();
  c_pf_->Zero();
  c_fp_->Zero();
  c_fa_->Zero();
  c_pa_->Zero();

  c_rhs_s_->PutScalar(0.0);
  c_rhs_pf_->PutScalar(0.0);
  c_rhs_f_->PutScalar(0.0);

  k_pf_porofluid->Zero();

  const CORE::ADAPTER::Coupling& couppff_fpsi = *icoup_pf_f_;
  const CORE::ADAPTER::Coupling& coupsf_fpsi = *icoup_ps_f_;
  const CORE::ADAPTER::Coupling& coup_ps_a_fpsi = *icoup_ps_a_;

  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  INPAR::FPSI::PartitionedCouplingMethod method =
      CORE::UTILS::IntegralValue<INPAR::FPSI::PartitionedCouplingMethod>(
          fpsidynparams, "PARTITIONED");

  if (method != INPAR::FPSI::nocoupling)
  {
    // set general vector values needed by elements

    poro_field()->fluid_field()->discretization()->ClearState();

    poro_field()->fluid_field()->discretization()->set_state(
        0, "dispnp", poro_field()->fluid_field()->Dispnp());

    poro_field()->fluid_field()->discretization()->set_state(
        0, "gridv", poro_field()->fluid_field()->GridVel());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "dispn", poro_field()->fluid_field()->Dispn());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "veln", poro_field()->fluid_field()->Veln());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "velaf", poro_field()->fluid_field()->Velnp());
    poro_field()->fluid_field()->discretization()->set_state(
        0, "velnp", poro_field()->fluid_field()->Velnp());

    fluid_field()->discretization()->ClearState();

    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->Dispnp());
    fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->GridVel());
    fluid_field()->discretization()->set_state(0, "dispn", fluid_field()->Dispn());
    fluid_field()->discretization()->set_state(0, "veln", fluid_field()->Veln());
    fluid_field()->discretization()->set_state(0, "velaf", fluid_field()->Velnp());
    fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->Velnp());

    // create the parameters for the discretization
    Teuchos::ParameterList fparams;

    // action for elements
    fparams.set<int>("action", FLD::fpsi_coupling);
    fparams.set("timescale", poro_field()->fluid_field()->residual_scaling());

    fparams.set("dt", fpsidynparams.get<double>("TIMESTEP"));
    fparams.set<int>("Physical Type", poro_field()->fluid_field()->PhysicalType());

    if (method == INPAR::FPSI::monolithic)
    {
      fparams.set<std::string>("fillblock", "Porofluid_Freefluid");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      CORE::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
          0,                                       // fluiddofset for column
          k_pf_porofluid,                          // coupling matrix with fluid rowmap
          Teuchos::null,                           // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      // what's the current problem type? Is it a fps3i problem?
      CORE::ProblemType probtype = GLOBAL::Problem::Instance()->GetProblemType();

      if (probtype == CORE::ProblemType::fps3i)
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
      k_pf_porofluid->Complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform_)(*k_pf_porofluid, 1.0,
            CORE::ADAPTER::CouplingSlaveConverter(couppff_fpsi), c_pf_->Matrix(1, 0), true);
      }

      fparams.set<std::string>("fillblock", "Porofluid_Structure");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->Zero();
      CORE::FE::AssembleStrategy fluidstrategy21(0,  // porofluiddofset for row
          0,                                         // porofluiddofset for column
          k_pf_porofluid,                            // coupling matrix with fluid rowmap
          Teuchos::null,                             // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, fluidstrategy21, "fpsi_coupling");

      k_pf_porofluid->Complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());
      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowcoltransform2_)(*k_pf_porofluid, 1.0,
            CORE::ADAPTER::CouplingSlaveConverter(
                couppff_fpsi),  // row converter: important to use slave converter
            CORE::ADAPTER::CouplingSlaveConverter(
                coupsf_fpsi),  //  col converter: important to use slave converter
            c_pp_->Matrix(1, 0),
            false,  // bool exactmatch = true (default)
            true);
      }

      fparams.set<std::string>("fillblock", "Fluid_Porofluid");
      fparams.set("InterfaceFacingElementMap", poro_fluid_fluid_interface_map_);
      k_fp_porofluid->Zero();
      CORE::FE::AssembleStrategy porofluidstrategy(0,  // porofluiddofset for row
          0,                                           // porofluiddofset for column
          k_fp_porofluid,                              // porofluid-structure coupling matrix
          Teuchos::null,                               // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      poro_field()->fluid_field()->discretization()->evaluate_condition(
          fparams, porofluidstrategy, "fpsi_coupling");
      k_fp_porofluid->Complete(poro_field()->FluidDomainMap(), poro_field()->FluidRangeMap());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform2_)(*k_fp_porofluid, 1.0,
            CORE::ADAPTER::CouplingMasterConverter(couppff_fpsi), c_fp_->Matrix(0, 1),
            true);  // add
      }

      fparams.set<std::string>("fillblock", "Fluid_Structure");
      fparams.set("InterfaceFacingElementMap", poro_fluid_fluid_interface_map_);

      // move me somewhere else
      Teuchos::RCP<CORE::LINALG::SparseMatrix> k_pfs_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *(poro_field()->fluid_field()->discretization()->dof_row_map()), 81, true, true));

      k_pfs_->UnComplete();

      CORE::FE::AssembleStrategy structurestrategy(0,  // porofluiddofset for row
          1,                                           // structuredofset for column
          k_pfs_,                                      // coupling matrix with porofluid rowmap
          Teuchos::null,                               // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      poro_field()->fluid_field()->discretization()->evaluate_condition(
          fparams, structurestrategy, "fpsi_coupling");

      k_pfs_->Complete(poro_field()->StructureDomainMap(), poro_field()->FluidRangeMap());
      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform3_)(*k_pfs_, 1.0,
            CORE::ADAPTER::CouplingMasterConverter(couppff_fpsi), c_fp_->Matrix(0, 0),
            true);  // add
      }

      ///// Fluid_Structure (fluid part / linearization of tangentials with respect to
      /// displacements)
      k_pf_porofluid->Reset();
      k_pf_porofluid->UnComplete();

      CORE::FE::AssembleStrategy structurestrategy2(0,  // porofluiddofset for row
          0,                                            // fluiddofset for column
          k_pf_porofluid,                               // coupling matrix with porofluid rowmap
          Teuchos::null,                                // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      fluid_field()->discretization()->evaluate_condition(
          fparams, structurestrategy2, "fpsi_coupling");

      k_pf_porofluid->Complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingcoltransform_)(fluid_field()->BlockSystemMatrix()->FullRowMap(),
            fluid_field()->BlockSystemMatrix()->FullColMap(), *k_pf_porofluid, 1.0,
            CORE::ADAPTER::CouplingSlaveConverter(
                coupsf_fpsi),  // row converter: important to use slave converter
            c_fp_->Matrix(0, 0),
            false,  // bool exactmatch = true (default)
            true);
      }

      fparams.set<std::string>("fillblock", "Fluid_Fluid");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->Zero();
      k_pf_porofluid->UnComplete();

      CORE::FE::AssembleStrategy fluidfluidstrategy(0,  // fluiddofset for row
          0,                                            // fluiddofset for column
          c_ff_,                                        // porofluid-structure coupling matrix
          Teuchos::null,                                // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, fluidfluidstrategy, "fpsi_coupling");

      fparams.set<std::string>("fillblock", "Structure_Fluid");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->Zero();
      k_pf_porofluid->UnComplete();

      CORE::FE::AssembleStrategy structurefluidstrategy(0,  // fluid dofset for row
          0,                                                // fluid dofset for column
          k_pf_porofluid,                                   // coupling matrix with fluid rowmap
          Teuchos::null,                                    // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, structurefluidstrategy, "fpsi_coupling");

      k_pf_porofluid->Complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());
      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform4_)(
            *Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(k_pf_porofluid), 1.0,
            CORE::ADAPTER::CouplingSlaveConverter(coupsf_fpsi),  // important to use slave converter
            c_pf_->Matrix(0, 0),
            true);  // add
      }

      fparams.set<std::string>("fillblock", "Structure_Structure");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
      k_pf_porofluid->Zero();
      k_pf_porofluid->UnComplete();

      CORE::FE::AssembleStrategy structurestructurestrategy(0,  // fluid dofset for row
          0,                                                    // fluid dofset for column
          k_pf_porofluid,                                       // coupling matrix with fluid rowmap
          Teuchos::null,                                        // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      fluid_field()->discretization()->evaluate_condition(
          fparams, structurestructurestrategy, "fpsi_coupling");
      // condense linearization with respect to the ale mesh motion (interface structural
      // displacements = interface ale displacements)
      k_pf_porofluid->Complete(*fluid_field()->dof_row_map(), *fluid_field()->dof_row_map());

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowcoltransform_)(*k_pf_porofluid, 1.0,
            CORE::ADAPTER::CouplingSlaveConverter(
                coupsf_fpsi),  // row converter: important to use slave converter
            CORE::ADAPTER::CouplingSlaveConverter(
                coupsf_fpsi),  // col converter: important to use slave converter
            c_pp_->Matrix(0, 0),
            false,  // bool exactmatch = true (default)
            true);  // add
      }

      // Process inner ale dofs
      fparams.set<std::string>("fillblock", "Structure_Ale");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);

      // temporal matrix
      // todo (initialization should be avoided in every iteration...)
      Teuchos::RCP<CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>>
          temp6 = Teuchos::rcp(
              new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
                  *ale_field()->Interface(), *fluid_field()->Interface(), 81, false,
                  false));  // Use fluid_field()->Interface =

      // assemble into fluid row and column dofs -> need to transform rows to structure dofs and
      // cols to ale dofs
      CORE::FE::AssembleStrategy structurealestrategy(0,  // fluid dofset for row
          1,                                              // ale dofset for column
          temp6,                                          // coupling matrix with fluid rowmap
          Teuchos::null,                                  // no other matrix or vectors
          Teuchos::null, Teuchos::null, Teuchos::null);

      // evaluate coupling terms
      fluid_field()->discretization()->evaluate_condition(
          fparams, structurealestrategy, "fpsi_coupling");
      temp6->Complete();  // for row transform!

      {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowtransform5_)(temp6->Matrix(FLD::UTILS::MapExtractor::cond_other,
                                      ALE::UTILS::MapExtractor::cond_other),
            1.0,
            CORE::ADAPTER::CouplingSlaveConverter(coupsf_fpsi),  // important to use slave converter
            c_pa_->Matrix(0, 0), false);
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

    CORE::FE::AssembleStrategy rhscontistrategy(0,  // fluid dofset for row
        0,                                          // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    fluid_field()->discretization()->evaluate_condition(fparams, rhscontistrategy, "fpsi_coupling");

    // extract FPSI part of the fluid field
    temprhs = fluidvelpres_extractor_->ExtractCondVector(temprhs);

    // replace global fluid interface dofs through porofluid interface dofs
    temprhs = iFluidToPorofluid(temprhs);

    // insert porofluid interface entries into vector with full porofield length (0: inner dofs of
    // structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of
    // porofluid )
    porofluid_extractor_->InsertCondVector(temprhs, c_rhs_pf_);

    // add vector with full porofield length to global rhs

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    fparams.set<std::string>("fillblock", "structure");
    fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
    temprhs2 = Teuchos::rcp(new Epetra_Vector(*poro_field()->dof_row_map(), true));

    CORE::FE::AssembleStrategy rhsstructurestrategy(0,  // fluid dofset for row
        0,                                              // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    fluid_field()->discretization()->evaluate_condition(
        fparams, rhsstructurestrategy, "fpsi_coupling");

    // extract FPSI part of the fluid field
    temprhs = fluidvel_extractor_->ExtractCondVector(temprhs);  //
    // replace global fluid interface dofs through porofluid interface dofs
    temprhs = iFluidToPorostruct(temprhs);  //
    // insert porofluid interface
    porostruct_extractor_->AddCondVector(temprhs, c_rhs_s_);

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    fparams.set<std::string>("fillblock", "fluid");
    fparams.set("InterfaceFacingElementMap", poro_fluid_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*poro_field()->fluid_field()->dof_row_map(0), true));
    temprhs2 = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(0), true));

    CORE::FE::AssembleStrategy rhsfluidstrategy(0,  // fluid dofset for row
        0,                                          // fluid dofset for column
        Teuchos::null, Teuchos::null,
        temprhs,  // rhs vector
        Teuchos::null, Teuchos::null);

    poro_field()->fluid_field()->discretization()->evaluate_condition(
        fparams, rhsfluidstrategy, "fpsi_coupling");
    // extract FPSI part of the poro fluid field
    temprhs = porofluid_extractor_->ExtractCondVector(temprhs);  //

    // replace global fluid interface dofs through porofluid interface dofs
    temprhs = iPorofluidToFluid(temprhs);
    // insert porofluid interface entries into vector with full fluidfield length
    fluidvelpres_extractor_->InsertCondVector(temprhs, temprhs2);
    // add vector with full porofield length to global rhs
    c_rhs_f_->Update(1.0, *temprhs2, 0.0);

    temprhs->PutScalar(0.0);
    temprhs2->PutScalar(0.0);

    fparams.set<std::string>("fillblock", "fluidfluid");  // (wot,tangentialfac*uot) part
    fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);
    temprhs = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(0), true));

    CORE::FE::AssembleStrategy rhsfluidfluidstrategy(0,  // fluid dofset for row
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
    if (method == INPAR::FPSI::monolithic or method == INPAR::FPSI::RobinNeumann)
    {
      fparams.set<std::string>("fillblock", "NeumannIntegration");
      fparams.set("InterfaceFacingElementMap", fluid_poro_fluid_interface_map_);

      CORE::FE::AssembleStrategy rhsfluidfluidstrategy2(0,  // fluid dofset for row
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

        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> tmp_c_fa = Teuchos::rcp(
            new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
                *ale_field()->Interface(), *fluid_field()->FPSIInterface(), 81, false, true));

        CORE::FE::AssembleStrategy rhsfluidfluidstrategy3(0,  // fluid dofset for row
            1,                                                // ale dofset for column
            tmp_c_fa,                                         // coupling matrix with fluid rowmap
            Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

        fluid_field()->discretization()->evaluate_condition(
            fparams, rhsfluidfluidstrategy3, "NeumannIntegration");
        tmp_c_fa->Complete();

        // Add all inner ale parts to c_fa_ directly
        c_fa_->Add(tmp_c_fa->Matrix(
                       FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_other),
            false, 1.0, 0.0);
        c_fa_->Add(tmp_c_fa->Matrix(
                       FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_other),
            false, 1.0, 1.0);

        //-->now transform ale fpsi block to structure (is condensed)!!!
        // still in ale domain map!!!
        CORE::LINALG::SparseMatrix tmp_c_fp =
            CORE::LINALG::SparseMatrix(*fluid_field()->dof_row_map(), 81);

        // Add all condensed parts to tmp_c_fa...
        tmp_c_fp.Add(tmp_c_fa->Matrix(
                         FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_fpsi),
            false, 1.0, 0.0);
        tmp_c_fp.Add(tmp_c_fa->Matrix(
                         FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_fpsi),
            false, 1.0, 1.0);
        tmp_c_fp.Add(tmp_c_fa->Matrix(
                         FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_fsi),
            false, 1.0, 1.0);
        tmp_c_fp.Add(tmp_c_fa->Matrix(
                         FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_fsi),
            false, 1.0, 1.0);
        tmp_c_fp.Complete(*ale_field()->Interface()->FPSICondMap(), *fluid_field()->dof_row_map());

        // For Ale Condensation ==> AleColumns to StructuralColumns
        (*couplingcoltransform2_)(ale_field()->BlockSystemMatrix()->FullRowMap(),
            ale_field()->BlockSystemMatrix()->FullColMap(), tmp_c_fp, 1.0,
            CORE::ADAPTER::CouplingSlaveConverter(
                coup_ps_a_fpsi),  // row converter: important to use slave converter
            c_fp_->Matrix(0, 0),
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
    poro_field()->fluid_field()->discretization()->ClearState();
    fluid_field()->discretization()->ClearState();

  }  // if not nocoupling

  ////////////////////////////
  // DONE -> Complete Coupling Matrixes
  c_ff_->Complete();
  c_pp_->Complete();
}

/*----------------------------------------------------------------------/
| set hydraulic conductivity
/----------------------------------------------------------------------*/
void FPSI::FpsiCoupling::SetConductivity(double conduct) { conductivity_ = conduct; }

FOUR_C_NAMESPACE_CLOSE
