/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for XFEM

\level 2


*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_xfem.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_cut.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::XFEM::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL", false, "");

  // OUTPUT options
  Core::UTILS::BoolParameter("GMSH_DEBUG_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", &xfem_general);
  Core::UTILS::BoolParameter("GMSH_DEBUG_OUT_SCREEN", "No",
      "Do you want to be informed, if Gmsh output is written?", &xfem_general);
  Core::UTILS::BoolParameter("GMSH_SOL_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", &xfem_general);
  Core::UTILS::BoolParameter("GMSH_TIMINT_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", &xfem_general);
  Core::UTILS::BoolParameter("GMSH_EOS_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", &xfem_general);
  Core::UTILS::BoolParameter("GMSH_DISCRET_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", &xfem_general);
  Core::UTILS::BoolParameter("GMSH_CUT_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", &xfem_general);

  Core::UTILS::IntParameter(
      "MAX_NUM_DOFSETS", 3, "Maximum number of volumecells in the XFEM element", &xfem_general);

  setStringToIntegralParameter<int>("NODAL_DOFSET_STRATEGY", "full",
      "Strategy used for the nodal dofset management per node",
      tuple<std::string>(
          "OneDofset_PerNodeAndPosition", "ConnectGhostDofsets_PerNodeAndPosition", "full"),
      tuple<int>(Inpar::Cut::NDS_Strategy_OneDofset_PerNodeAndPosition,
          Inpar::Cut::NDS_Strategy_ConnectGhostDofsets_PerNodeAndPosition,
          Inpar::Cut::NDS_Strategy_full),
      &xfem_general);

  // Integration options
  setStringToIntegralParameter<int>("VOLUME_GAUSS_POINTS_BY", "Tessellation",
      "Method for finding Gauss Points for the cut volumes",
      tuple<std::string>("Tessellation", "MomentFitting", "DirectDivergence"),
      tuple<int>(Inpar::Cut::VCellGaussPts_Tessellation, Inpar::Cut::VCellGaussPts_MomentFitting,
          Inpar::Cut::VCellGaussPts_DirectDivergence),
      &xfem_general);

  setStringToIntegralParameter<int>("BOUNDARY_GAUSS_POINTS_BY", "Tessellation",
      "Method for finding Gauss Points for the boundary cells",
      tuple<std::string>("Tessellation", "MomentFitting", "DirectDivergence"),
      tuple<int>(Inpar::Cut::BCellGaussPts_Tessellation, Inpar::Cut::BCellGaussPts_MomentFitting,
          Inpar::Cut::BCellGaussPts_DirectDivergence),
      &xfem_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_dyn = list->sublist("XFLUID DYNAMIC", false, "");

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_general = xfluid_dyn.sublist("GENERAL", false, "");

  // Do we use more than one fluid discretization?
  Core::UTILS::BoolParameter("XFLUIDFLUID", "no", "Use an embedded fluid patch.", &xfluid_general);

  // How many monolithic steps we keep the fluidfluid-boundary fixed
  Core::UTILS::IntParameter(
      "RELAXING_ALE_EVERY", 1, "Relaxing Ale after how many monolithic steps", &xfluid_general);

  Core::UTILS::BoolParameter("RELAXING_ALE", "yes",
      "switch on/off for relaxing Ale in monolithic fluid-fluid-fsi", &xfluid_general);

  Core::UTILS::DoubleParameter(
      "XFLUIDFLUID_SEARCHRADIUS", 1.0, "Radius of the search tree", &xfluid_general);

  // xfluidfluid-fsi-monolithic approach
  setStringToIntegralParameter<int>("MONOLITHIC_XFFSI_APPROACH", "xffsi_fixedALE_partitioned",
      "The monolithic approach for xfluidfluid-fsi",
      tuple<std::string>(
          "xffsi_full_newton", "xffsi_fixedALE_interpolation", "xffsi_fixedALE_partitioned"),
      tuple<int>(Inpar::XFEM::XFFSI_Full_Newton,      // xffsi with no fixed xfem-coupling
          Inpar::XFEM::XFFSI_FixedALE_Interpolation,  // xffsi with fixed xfem-coupling in every
                                                      // newtonstep and interpolations for
                                                      // embedded-dis afterwards
          Inpar::XFEM::XFFSI_FixedALE_Partitioned     // xffsi with fixed xfem-coupling in every
                                                      // newtonstep and solving fluid-field again
          ),
      &xfluid_general);

  // xfluidfluid time integration approach
  setStringToIntegralParameter<int>("XFLUIDFLUID_TIMEINT", "Xff_TimeInt_FullProj",
      "The xfluidfluid-timeintegration approach",
      tuple<std::string>("Xff_TimeInt_FullProj", "Xff_TimeInt_ProjIfMoved",
          "Xff_TimeInt_KeepGhostValues", "Xff_TimeInt_IncompProj"),
      tuple<int>(Inpar::XFEM::Xff_TimeInt_FullProj,  // always project nodes from embedded to
                                                     // background nodes
          Inpar::XFEM::Xff_TimeInt_ProjIfMoved,  // project nodes just if the status of background
                                                 // nodes changed
          Inpar::XFEM::Xff_TimeInt_KeepGhostValues,  // always keep the ghost values of the
                                                     // background discretization
          Inpar::XFEM::Xff_TimeInt_IncompProj  // after projecting nodes do a incompressibility
                                               // projection
          ),
      &xfluid_general);

  setStringToIntegralParameter<int>("XFLUID_TIMEINT", "STD=COPY/SL_and_GHOST=COPY/GP",
      "The xfluid time integration approach",
      tuple<std::string>("STD=COPY_and_GHOST=COPY/GP", "STD=COPY/SL_and_GHOST=COPY/GP",
          "STD=SL(boundary-zone)_and_GHOST=GP", "STD=COPY/PROJ_and_GHOST=COPY/PROJ/GP"),
      tuple<int>(
          Inpar::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP,  // STD= only copy,
                                                                              // GHOST= copy or
                                                                              // ghost penalty
                                                                              // reconstruction
          Inpar::XFEM::
              Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP,  // STD= copy or SL,
                                                                           // GHOST= copy or ghost
                                                                           // penalty reconstruction
          Inpar::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP,  // STD= only SL on
                                                                             // whole boundary zone,
                                                                             // GHOST= ghost penalty
                                                                             // reconstruction
          Inpar::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP),
      &xfluid_general);

  Core::UTILS::BoolParameter("ALE_XFluid", "no", "XFluid is Ale Fluid?", &xfluid_general);

  // for new OST-implementation: which interface terms to be evaluated for previous time step
  setStringToIntegralParameter<int>("INTERFACE_TERMS_PREVIOUS_STATE",
      "PreviousState_only_consistency",
      "how to treat interface terms from previous time step (new OST)",
      tuple<std::string>("PreviousState_only_consistency", "PreviousState_full"),
      tuple<int>(Inpar::XFEM::PreviousState_only_consistency,  /// evaluate only consistency terms
                                                               /// for previous time step
          Inpar::XFEM::PreviousState_full  /// evaluate consistency, adjoint consistency and penalty
                                           /// terms or previous time step
          ),
      &xfluid_general);

  Core::UTILS::BoolParameter("XFLUID_TIMEINT_CHECK_INTERFACETIPS", "Yes",
      "Xfluid TimeIntegration Special Check if node has changed the side!", &xfluid_general);

  Core::UTILS::BoolParameter("XFLUID_TIMEINT_CHECK_SLIDINGONSURFACE", "Yes",
      "Xfluid TimeIntegration Special Check if node is sliding on surface!", &xfluid_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_stab = xfluid_dyn.sublist("STABILIZATION", false, "");

  // Boundary-Coupling options
  setStringToIntegralParameter<int>("COUPLING_METHOD", "Nitsche",
      "method how to enforce embedded boundary/coupling conditions at the interface",
      tuple<std::string>("Hybrid_LM_Cauchy_stress", "Hybrid_LM_viscous_stress", "Nitsche"),
      tuple<int>(
          Inpar::XFEM::Hybrid_LM_Cauchy_stress,   // Cauchy stress-based mixed/hybrid formulation
          Inpar::XFEM::Hybrid_LM_viscous_stress,  // viscous stress-based mixed/hybrid formulation
          Inpar::XFEM::Nitsche                    // Nitsche's formulation
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("HYBRID_LM_L2_PROJ", "part_ele_proj",
      "perform the L2 projection between stress fields on whole element or on fluid part?",
      tuple<std::string>("full_ele_proj", "part_ele_proj"),
      tuple<int>(
          Inpar::XFEM::Hybrid_LM_L2_Proj_full,  // L2 stress projection on whole fluid element
          Inpar::XFEM::Hybrid_LM_L2_Proj_part   // L2 stress projection on partial fluid element
                                                // volume
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_ADJOINT_SYMMETRY", "yes",
      "viscous and adjoint viscous interface terms with matching sign?",
      tuple<std::string>("yes", "no", "sym", "skew", "none"),
      tuple<int>(Inpar::XFEM::adj_sym, Inpar::XFEM::adj_skew, Inpar::XFEM::adj_sym,
          Inpar::XFEM::adj_skew, Inpar::XFEM::adj_none),
      &xfluid_stab);

  // viscous and convective Nitsche/MSH stabilization parameter
  Core::UTILS::DoubleParameter(
      "NIT_STAB_FAC", 35.0, " ( stabilization parameter for Nitsche's penalty term", &xfluid_stab);
  Core::UTILS::DoubleParameter("NIT_STAB_FAC_TANG", 35.0,
      " ( stabilization parameter for Nitsche's penalty tangential term", &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_TRACE_ESTIMATE", "CT_div_by_hk",
      "how to estimate the scaling from the trace inequality in Nitsche's method",
      tuple<std::string>("CT_div_by_hk", "eigenvalue"),
      tuple<int>(
          Inpar::XFEM::ViscStab_TraceEstimate_CT_div_by_hk,  // estimate the trace inequality by a
                                                             // trace-constant CT and hk: CT/hk
          Inpar::XFEM::ViscStab_TraceEstimate_eigenvalue     // estimate the trace inequality by
                                                          // solving a local element-wise eigenvalue
                                                          // problem
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("UPDATE_EIGENVALUE_TRACE_ESTIMATE", "every_iter",
      "how often should the local eigenvalue problem be updated",
      tuple<std::string>("every_iter", "every_timestep", "once"),
      tuple<int>(Inpar::XFEM::Eigenvalue_update_every_iter,
          Inpar::XFEM::Eigenvalue_update_every_timestep, Inpar::XFEM::Eigenvalue_update_once),
      &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_HK", "ele_vol_div_by_max_ele_surf",
      "how to define the characteristic element length in cut elements",
      tuple<std::string>("vol_equivalent", "cut_vol_div_by_cut_surf", "ele_vol_div_by_cut_surf",
          "ele_vol_div_by_ele_surf", "ele_vol_div_by_max_ele_surf"),
      tuple<int>(Inpar::XFEM::ViscStab_hk_vol_equivalent,    /// volume equivalent element diameter
          Inpar::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf,  /// physical partial/cut volume divided
                                                             /// by physical partial/cut surface
                                                             /// measure ( used to estimate the
                                                             /// cut-dependent inverse estimate on
                                                             /// cut elements, not useful for sliver
                                                             /// and/or dotted cut situations)
          Inpar::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf,  /// full element volume divided by
                                                             /// physical partial/cut surface
                                                             /// measure ( used to estimate the
                                                             /// cut-dependent inverse estimate on
                                                             /// cut elements, however, avoids
                                                             /// problems with sliver cuts, not
                                                             /// useful for dotted cuts)
          Inpar::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf,  /// full element volume divided by
                                                             /// surface measure ( used for uncut
                                                             /// situations, standard weak Dirichlet
                                                             /// boundary/coupling conditions)
          Inpar::XFEM::
              ViscStab_hk_ele_vol_div_by_max_ele_surf  /// default: full element volume divided by
                                                       /// maximal element surface measure ( used to
                                                       /// estimate the trace inequality for
                                                       /// stretched elements in combination with
                                                       /// ghost-penalties)
          ),
      &xfluid_stab);


  setStringToIntegralParameter<int>("CONV_STAB_SCALING", "none",
      "scaling factor for viscous interface stabilization (Nitsche, MSH)",
      tuple<std::string>("inflow", "abs_inflow", "none"),
      tuple<int>(Inpar::XFEM::ConvStabScaling_inflow,  // scaling with max(0,-u*n)
          Inpar::XFEM::ConvStabScaling_abs_inflow,     // scaling with |u*n|
          Inpar::XFEM::ConvStabScaling_none            // no convective stabilization
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("XFF_CONV_STAB_SCALING", "none",
      "scaling factor for convective interface stabilization of fluid-fluid Coupling",
      tuple<std::string>("inflow", "averaged", "none"),
      tuple<int>(Inpar::XFEM::XFF_ConvStabScaling_upwinding,  // one-sided inflow stabilization
          Inpar::XFEM::XFF_ConvStabScaling_only_averaged,     // averaged inflow stabilization
          Inpar::XFEM::XFF_ConvStabScaling_none               // no convective stabilization
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_COMBO", "max",
      "choose the maximum from viscous and convective contributions or just sum both up",
      tuple<std::string>("max", "sum"),
      tuple<int>(Inpar::XFEM::MassConservationCombination_max,  /// use the maximum contribution
          Inpar::XFEM::MassConservationCombination_sum  /// sum viscous and convective contributions
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_SCALING", "only_visc",
      "apply additional scaling of penalty term to enforce mass conservation for "
      "convection-dominated flow",
      tuple<std::string>("full", "only_visc"),
      tuple<int>(Inpar::XFEM::MassConservationScaling_full,  /// apply mass-conserving convective
                                                             /// scaling additionally
          Inpar::XFEM::MassConservationScaling_only_visc     /// use only the viscous scaling
          ),
      &xfluid_stab);

  Core::UTILS::BoolParameter("GHOST_PENALTY_STAB", "no",
      "switch on/off ghost penalty interface stabilization", &xfluid_stab);

  Core::UTILS::BoolParameter("GHOST_PENALTY_TRANSIENT_STAB", "no",
      "switch on/off ghost penalty transient interface stabilization", &xfluid_stab);

  Core::UTILS::BoolParameter("GHOST_PENALTY_2nd_STAB", "no",
      "switch on/off ghost penalty interface stabilization for 2nd order derivatives",
      &xfluid_stab);
  Core::UTILS::BoolParameter("GHOST_PENALTY_2nd_STAB_NORMAL", "no",
      "switch between ghost penalty interface stabilization for 2nd order derivatives in normal or "
      "all spatial directions",
      &xfluid_stab);


  Core::UTILS::DoubleParameter("GHOST_PENALTY_FAC", 0.1,
      "define stabilization parameter ghost penalty interface stabilization", &xfluid_stab);

  Core::UTILS::DoubleParameter("GHOST_PENALTY_TRANSIENT_FAC", 0.001,
      "define stabilization parameter ghost penalty transient interface stabilization",
      &xfluid_stab);

  //  // NOT chosen optimally!
  //  Core::UTILS::DoubleParameter("GHOST_PENALTY_2nd_FAC", 1.0,"define stabilization parameter
  //  ghost penalty 2nd order viscous interface stabilization",&xfluid_stab);
  //  Core::UTILS::DoubleParameter("GHOST_PENALTY_PRESSURE_2nd_FAC", 1.0,"define stabilization
  //  parameter ghost penalty 2nd order pressure interface stabilization",&xfluid_stab);

  Core::UTILS::DoubleParameter("GHOST_PENALTY_2nd_FAC", 0.05,
      "define stabilization parameter ghost penalty 2nd order viscous interface stabilization",
      &xfluid_stab);
  Core::UTILS::DoubleParameter("GHOST_PENALTY_PRESSURE_2nd_FAC", 0.05,
      "define stabilization parameter ghost penalty 2nd order pressure interface stabilization",
      &xfluid_stab);


  Core::UTILS::BoolParameter("XFF_EOS_PRES_EMB_LAYER", "no",
      "switch on/off edge-based pressure stabilization on interface-contributing elements of the "
      "embedded fluid",
      &xfluid_stab);

  Core::UTILS::BoolParameter("IS_PSEUDO_2D", "no",
      "modify viscous interface stabilization due to the vanishing polynomial in third dimension "
      "when using strong Dirichlet conditions to block polynomials in one spatial dimension",
      &xfluid_stab);

  Core::UTILS::BoolParameter("GHOST_PENALTY_ADD_INNER_FACES", "no",
      "Apply ghost penalty stabilization also for inner faces if this is possible due to the "
      "dofsets",
      &xfluid_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfsi_monolithic = xfluid_dyn.sublist("XFPSI MONOLITHIC", false, "");

  Core::UTILS::IntParameter(
      "ITEMIN", 1, "How many iterations are performed minimal", &xfsi_monolithic);
  Core::UTILS::IntParameter(
      "ITEMAX_OUTER", 5, "How many outer iterations are performed maximal", &xfsi_monolithic);
  Core::UTILS::BoolParameter("ND_NEWTON_DAMPING", "no",
      "Activate Newton damping based on residual and increment", &xfsi_monolithic);
  Core::UTILS::DoubleParameter("ND_MAX_DISP_ITERINC", -1.0,
      "Maximal displacement increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("ND_MAX_VEL_ITERINC", -1.0,
      "Maximal fluid velocity increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("ND_MAX_PRES_ITERINC", -1.0,
      "Maximal fluid pressure increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("ND_MAX_PVEL_ITERINC", -1.0,
      "Maximal porofluid velocity increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("ND_MAX_PPRES_ITERINC", -1.0,
      "Maximal porofluid pressure increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("CUT_EVALUATE_MINTOL", 0.0,
      "Minimal value of the maximal strucutral displacement for which the CUT is evaluate in this "
      "iteration!",
      &xfsi_monolithic);
  Core::UTILS::IntParameter("CUT_EVALUATE_MINITER", 0,
      "Minimal number of nonlinear iterations, before the CUT is potentially not evaluated",
      &xfsi_monolithic);
  Core::UTILS::BoolParameter("EXTRAPOLATE_TO_ZERO", "no",
      "the extrapolation of the fluid stress in the contact zone is relaxed to zero after a "
      "certain distance",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("POROCONTACTFPSI_HFRACTION", 1.0,
      "factor of element size, when transition between FPSI and PSCI is started!",
      &xfsi_monolithic);
  Core::UTILS::DoubleParameter("POROCONTACTFPSI_FULLPCFRACTION", 0.0,
      "ration of gap/(POROCONTACTFPSI_HFRACTION*h) when full PSCI is started!", &xfsi_monolithic);
  Core::UTILS::BoolParameter("USE_PORO_PRESSURE", "yes",
      "the extrapolation of the fluid stress in the contact zone is relaxed to zero after a "
      "certtain distance",
      &xfsi_monolithic);
}


void Inpar::XFEM::SetValidConditions(
    const std::vector<Teuchos::RCP<Input::LineComponent>>& dirichletbundcomponents,
    const std::vector<Teuchos::RCP<Input::LineComponent>>& neumanncomponents,
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  std::vector<Teuchos::RCP<Input::LineComponent>> xfemcomponents;

  xfemcomponents.push_back(Teuchos::rcp(new Input::IntComponent("label")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> movingfluid = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN FLUID MESH VOL CONDITIONS", "FluidMesh",
          "Fluid Mesh", Core::Conditions::FluidMesh, true, Core::Conditions::geometry_type_volume));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> fluidfluidcoupling = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN FLUID FLUID COUPLING SURF CONDITIONS",
          "FluidFluidCoupling", "FLUID FLUID Coupling", Core::Conditions::FluidFluidCoupling, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> ALEfluidcoupling =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN ALE FLUID COUPLING SURF CONDITIONS", "ALEFluidCoupling", "ALE FLUID Coupling",
          Core::Conditions::ALEFluidCoupling, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < xfemcomponents.size(); ++i)
  {
    movingfluid->add_component(xfemcomponents[i]);
    fluidfluidcoupling->add_component(xfemcomponents[i]);
    ALEfluidcoupling->add_component(xfemcomponents[i]);
  }

  condlist.push_back(fluidfluidcoupling);
  condlist.push_back(movingfluid);
  condlist.push_back(ALEfluidcoupling);


  /*--------------------------------------------------------------------*/
  // XFEM coupling conditions

  //*----------------*/
  // Displacement surface condition for XFEM WDBC and Neumann boundary conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_displacement = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN XFEM DISPLACEMENT SURF CONDITIONS",
          "XFEMSurfDisplacement", "XFEM Surf Displacement",
          Core::Conditions::XFEM_Surf_Displacement, true, Core::Conditions::geometry_type_surface));

  xfem_surf_displacement->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_displacement->add_component(Teuchos::rcp(new Input::IntComponent("label")));
  xfem_surf_displacement->add_component(Teuchos::rcp(new Input::SeparatorComponent("EVALTYPE")));
  xfem_surf_displacement->add_component(Teuchos::rcp(new Input::SelectionComponent("evaltype",
      "funct", Teuchos::tuple<std::string>("zero", "funct", "implementation"),
      Teuchos::tuple<std::string>("zero", "funct", "implementation"), true)));


  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_surf_displacement->add_component(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_surf_displacement);



  //*----------------*/
  // Levelset field condition components

  std::vector<Teuchos::RCP<Input::LineComponent>> levelsetfield_components;

  levelsetfield_components.push_back(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  levelsetfield_components.push_back(Teuchos::rcp(new Input::IntComponent("label")));

  levelsetfield_components.push_back(
      Teuchos::rcp(new Input::SeparatorComponent("LEVELSETFIELDNO")));
  levelsetfield_components.push_back(Teuchos::rcp(new Input::IntComponent("levelsetfieldno")));

  // define which boolean operator is used for combining this level-set field with the previous one
  // with smaller coupling id
  levelsetfield_components.push_back(Teuchos::rcp(new Input::SeparatorComponent("BOOLEANTYPE")));
  levelsetfield_components.push_back(Teuchos::rcp(new Input::SelectionComponent("booleantype",
      "none", Teuchos::tuple<std::string>("none", "cut", "union", "difference", "sym_difference"),
      Teuchos::tuple<std::string>("none", "cut", "union", "difference", "sym_difference"), false)));

  // define which complementary operator is applied after combining the level-set field with a
  // boolean operator with the previous one
  levelsetfield_components.push_back(Teuchos::rcp(new Input::SeparatorComponent("COMPLEMENTARY")));
  levelsetfield_components.push_back(Teuchos::rcp(new Input::IntComponent("complementary")));


  //*----------------*/
  // Levelset based Weak Dirichlet conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_levelset_wdbc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM LEVELSET WEAK DIRICHLET VOL CONDITIONS", "XFEMLevelsetWeakDirichlet",
          "XFEM Levelset Weak Dirichlet", Core::Conditions::XFEM_Levelset_Weak_Dirichlet, true,
          Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_wdbc->add_component(levelsetfield_components[i]);
  }

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_levelset_wdbc->add_component(dirichletbundcomponents[i]);
  }

  // optional: allow for random noise, set percentage used in uniform random distribution
  xfem_levelset_wdbc->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("RANDNOISE", "", true)));
  xfem_levelset_wdbc->add_component(
      Teuchos::rcp(new Input::RealComponent("randnoise", {0., true})));

  condlist.push_back(xfem_levelset_wdbc);

  //*----------------*/
  // Levelset based Neumann conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_levelset_neumann = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN XFEM LEVELSET NEUMANN VOL CONDITIONS",
          "XFEMLevelsetNeumann", "XFEM Levelset Neumann", Core::Conditions::XFEM_Levelset_Neumann,
          true, Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_neumann->add_component(levelsetfield_components[i]);
  }

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_levelset_neumann->add_component(neumanncomponents[i]);
  }

  // define if we use inflow stabilization on the xfem neumann surf condition
  xfem_levelset_neumann->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("INFLOW_STAB", "", true)));
  xfem_levelset_neumann->add_component(
      Teuchos::rcp(new Input::BoolComponent("InflowStab", false, true)));

  condlist.push_back(xfem_levelset_neumann);

  //*----------------*/
  // Levelset based Navier Slip conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_levelset_navier_slip =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM LEVELSET NAVIER SLIP VOL CONDITIONS", "XFEMLevelsetNavierSlip",
          "XFEM Levelset Navier Slip", Core::Conditions::XFEM_Levelset_Navier_Slip, true,
          Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_navier_slip->add_component(levelsetfield_components[i]);
  }

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SURFACE_PROJECTION")));
  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SelectionComponent("SURFACE_PROJECTION", "proj_normal",
          Teuchos::tuple<std::string>(
              "proj_normal", "proj_smoothed", "proj_normal_smoothed_comb", "proj_normal_phi"),
          Teuchos::tuple<int>(Inpar::XFEM::Proj_normal, Inpar::XFEM::Proj_smoothed,
              Inpar::XFEM::Proj_normal_smoothed_comb, Inpar::XFEM::Proj_normal_phi),
          true)));

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("L2_PROJECTION_SOLVER")));
  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("l2projsolv", {0, true, true, false})));

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_DIRICHLET_ID")));
  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id_dirch", {0, true, true, false})));

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_NEUMANN_ID")));
  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id_neumann", {0, true, true, false})));

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIPCOEFFICIENT")));
  xfem_levelset_navier_slip->add_component(Teuchos::rcp(new Input::RealComponent("slipcoeff")));

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIP_FUNCT", "", true)));
  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("funct", {0, false, false, true})));

  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("FORCE_ONLY_TANG_VEL", "", true)));
  xfem_levelset_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("force_tang_vel", {0, false, false, true})));

  condlist.push_back(xfem_levelset_navier_slip);

  // Add condition XFEM DIRICHLET/NEUMANN?

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_navier_slip_robin_dirch =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM ROBIN DIRICHLET VOL CONDITIONS", "XFEMRobinDirichletVol",
          "XFEM Robin Dirichlet Volume", Core::Conditions::XFEM_Robin_Dirichlet_Volume, true,
          Core::Conditions::geometry_type_volume));

  xfem_navier_slip_robin_dirch->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_DIRICHLET_ID")));
  xfem_navier_slip_robin_dirch->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id", {0, true, true, false})));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_navier_slip_robin_dirch->add_component(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_dirch);

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_navier_slip_robin_neumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM ROBIN NEUMANN VOL CONDITIONS", "XFEMRobinNeumannVol",
          "XFEM Robin Neumann Volume", Core::Conditions::XFEM_Robin_Neumann_Volume, true,
          Core::Conditions::geometry_type_volume));

  xfem_navier_slip_robin_neumann->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_NEUMANN_ID")));
  xfem_navier_slip_robin_neumann->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id", {0, true, true, false})));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_navier_slip_robin_neumann->add_component(neumanncomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_neumann);


  //*----------------*/
  // Levelset based Twophase conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_levelset_twophase = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS",
          "XFEMLevelsetTwophase", "XFEM Levelset Twophase",
          Core::Conditions::XFEM_Levelset_Twophase, true, Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_twophase->add_component(levelsetfield_components[i]);
  }

  condlist.push_back(xfem_levelset_twophase);

  //*----------------*/
  // Levelset based Combustion conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_levelset_combustion =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM LEVELSET COMBUSTION VOL CONDITIONS", "XFEMLevelsetCombustion",
          "XFEM Levelset Combustion", Core::Conditions::XFEM_Levelset_Combustion, true,
          Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_combustion->add_component(levelsetfield_components[i]);
  }

  // "The laminar flamespeed incorporates all chemical kinetics into the problem for now"
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("LAMINAR_FLAMESPEED", "", false)));
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::RealComponent("laminar_flamespeed")));

  // "Molecular diffusivity"
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("MOL_DIFFUSIVITY", "", false)));
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::RealComponent("mol_diffusivity")));

  // "The Markstein length takes flame curvature into account"
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("MARKSTEIN_LENGTH", "", false)));
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::RealComponent("markstein_length")));

  // interface transport in all directions or just in normal direction?
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("TRANSPORT_DIRECTIONS")));


  // define if the curvature shall be accounted for in computing the transport velocity
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("TRANSPORT_CURVATURE")));
  xfem_levelset_combustion->add_component(
      Teuchos::rcp(new Input::IntComponent("transport_curvature")));


  condlist.push_back(xfem_levelset_combustion);

  //*----------------*/
  // Surface Fluid-Fluid coupling conditions


  std::vector<Teuchos::RCP<Input::LineComponent>> xfluidfluidsurfcomponents;

  xfluidfluidsurfcomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfluidfluidsurfcomponents.push_back(Teuchos::rcp(new Input::IntComponent("label")));

  xfluidfluidsurfcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("COUPSTRATEGY",
      "xfluid", Teuchos::tuple<std::string>("xfluid", "embedded", "mean"),
      Teuchos::tuple<int>(
          Inpar::XFEM::Xfluid_Sided, Inpar::XFEM::Embedded_Sided, Inpar::XFEM::Mean))));


  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_fluidfluid =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM FLUIDFLUID SURF CONDITIONS", "XFEMSurfFluidFluid", "XFEM Surf FluidFluid",
          Core::Conditions::XFEM_Surf_FluidFluid, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < xfluidfluidsurfcomponents.size(); ++i)
    xfem_surf_fluidfluid->add_component(xfluidfluidsurfcomponents[i]);

  condlist.push_back(xfem_surf_fluidfluid);

  //*----------------*/
  // Surface partitioned XFSI boundary conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_fsi_part =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM FSI PARTITIONED SURF CONDITIONS", "XFEMSurfFSIPart", "XFEM Surf FSI Part",
          Core::Conditions::XFEM_Surf_FSIPart, true, Core::Conditions::geometry_type_surface));

  xfem_surf_fsi_part->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_fsi_part->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  // COUPSTRATEGY IS FLUID SIDED

  xfem_surf_fsi_part->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("INTLAW", "", true)));
  xfem_surf_fsi_part->add_component(Teuchos::rcp(new Input::SelectionComponent("INTLAW", "noslip",
      Teuchos::tuple<std::string>("noslip", "noslip_splitpen", "slip", "navslip"),
      Teuchos::tuple<int>(Inpar::XFEM::noslip, Inpar::XFEM::noslip_splitpen, Inpar::XFEM::slip,
          Inpar::XFEM::navierslip),
      true)));

  xfem_surf_fsi_part->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIPCOEFFICIENT", "", true)));
  xfem_surf_fsi_part->add_component(
      Teuchos::rcp(new Input::RealComponent("slipcoeff", {0., true})));

  xfem_surf_fsi_part->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIP_FUNCT", "", true)));
  xfem_surf_fsi_part->add_component(
      Teuchos::rcp(new Input::IntComponent("funct", {0, false, false, true})));

  condlist.push_back(xfem_surf_fsi_part);

  //*----------------*/
  // Surface monolithic XFSI coupling conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_fsi_mono =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM FSI MONOLITHIC SURF CONDITIONS", "XFEMSurfFSIMono", "XFEM Surf FSI Mono",
          Core::Conditions::XFEM_Surf_FSIMono, true, Core::Conditions::geometry_type_surface));

  xfem_surf_fsi_mono->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_fsi_mono->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  xfem_surf_fsi_mono->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("COUPSTRATEGY", "", true)));
  xfem_surf_fsi_mono->add_component(Teuchos::rcp(new Input::SelectionComponent("COUPSTRATEGY",
      "xfluid", Teuchos::tuple<std::string>("xfluid", "solid", "mean", "harmonic"),
      Teuchos::tuple<int>(Inpar::XFEM::Xfluid_Sided, Inpar::XFEM::Embedded_Sided, Inpar::XFEM::Mean,
          Inpar::XFEM::Harmonic),
      true)));

  xfem_surf_fsi_mono->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("INTLAW", "", true)));
  xfem_surf_fsi_mono->add_component(Teuchos::rcp(new Input::SelectionComponent("INTLAW", "noslip",
      Teuchos::tuple<std::string>(
          "noslip", "noslip_splitpen", "slip", "navslip", "navslip_contact"),
      Teuchos::tuple<int>(Inpar::XFEM::noslip, Inpar::XFEM::noslip_splitpen, Inpar::XFEM::slip,
          Inpar::XFEM::navierslip, Inpar::XFEM::navierslip_contact),
      true)));

  xfem_surf_fsi_mono->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIPCOEFFICIENT", "", true)));
  xfem_surf_fsi_mono->add_component(
      Teuchos::rcp(new Input::RealComponent("slipcoeff", {0., true})));

  xfem_surf_fsi_mono->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIP_FUNCT", "", true)));
  xfem_surf_fsi_mono->add_component(
      Teuchos::rcp(new Input::IntComponent("funct", {0, false, false, true})));

  condlist.push_back(xfem_surf_fsi_mono);

  //*----------------*/
  // Surface monolithic XFPI coupling conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_fpi_mono =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM FPI MONOLITHIC SURF CONDITIONS", "XFEMSurfFPIMono", "XFEM Surf FPI Mono",
          Core::Conditions::XFEM_Surf_FPIMono, true, Core::Conditions::geometry_type_surface));

  xfem_surf_fpi_mono->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_fpi_mono->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  xfem_surf_fpi_mono->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("BJ_COEFF", "", true)));
  xfem_surf_fpi_mono->add_component(Teuchos::rcp(new Input::RealComponent("bj_coeff", {0., true})));

  xfem_surf_fpi_mono->add_component(Teuchos::rcp(new Input::SelectionComponent("Variant", "BJ",
      Teuchos::tuple<std::string>("BJ", "BJS"), Teuchos::tuple<std::string>("BJ", "BJS"), true)));

  xfem_surf_fpi_mono->add_component(Teuchos::rcp(new Input::SelectionComponent("Method", "NIT",
      Teuchos::tuple<std::string>("NIT", "SUB"), Teuchos::tuple<std::string>("NIT", "SUB"), true)));

  xfem_surf_fpi_mono->add_component(Teuchos::rcp(new Input::SelectionComponent("Contact",
      "contact_no", Teuchos::tuple<std::string>("contact_no", "contact_yes"),
      Teuchos::tuple<std::string>("contact_no", "contact_yes"), true)));

  condlist.push_back(xfem_surf_fpi_mono);


  //*----------------*/
  // Surface Weak Dirichlet conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_wdbc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM WEAK DIRICHLET SURF CONDITIONS", "XFEMSurfWeakDirichlet",
          "XFEM Surf Weak Dirichlet", Core::Conditions::XFEM_Surf_Weak_Dirichlet, true,
          Core::Conditions::geometry_type_surface));

  xfem_surf_wdbc->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_wdbc->add_component(Teuchos::rcp(new Input::IntComponent("label")));
  xfem_surf_wdbc->add_component(Teuchos::rcp(new Input::SeparatorComponent("EVALTYPE")));
  xfem_surf_wdbc->add_component(
      Teuchos::rcp(new Input::SelectionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_surf_wdbc->add_component(dirichletbundcomponents[i]);
  }

  // optional: allow for random noise, set percentage used in uniform random distribution
  xfem_surf_wdbc->add_component(Teuchos::rcp(new Input::SeparatorComponent("RANDNOISE", "", true)));
  xfem_surf_wdbc->add_component(Teuchos::rcp(new Input::RealComponent("randnoise", {0., true})));

  condlist.push_back(xfem_surf_wdbc);


  //*----------------*/
  // Surface Neumann conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_neumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN XFEM NEUMANN SURF CONDITIONS",
          "XFEMSurfNeumann", "XFEM Surf Neumann", Core::Conditions::XFEM_Surf_Neumann, true,
          Core::Conditions::geometry_type_surface));

  xfem_surf_neumann->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_neumann->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_surf_neumann->add_component(neumanncomponents[i]);
  }

  // define if we use inflow stabilization on the xfem neumann surf condition
  xfem_surf_neumann->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("INFLOW_STAB", "", true)));
  xfem_surf_neumann->add_component(
      Teuchos::rcp(new Input::BoolComponent("InflowStab", false, true)));

  condlist.push_back(xfem_surf_neumann);

  //*----------------*/
  // Surface Navier Slip conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_navier_slip =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM NAVIER SLIP SURF CONDITIONS", "XFEMSurfNavierSlip", "XFEM Surf Navier Slip",
          Core::Conditions::XFEM_Surf_Navier_Slip, true, Core::Conditions::geometry_type_surface));

  xfem_surf_navier_slip->add_component(Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_navier_slip->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  xfem_surf_navier_slip->add_component(Teuchos::rcp(new Input::SeparatorComponent("EVALTYPE")));
  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::SelectionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));


  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_DIRICHLET_ID")));
  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id_dirch", {0, true, true, false})));

  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_NEUMANN_ID")));
  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id_neumann", {0, true, true, false})));

  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIPCOEFFICIENT")));
  xfem_surf_navier_slip->add_component(Teuchos::rcp(new Input::RealComponent("slipcoeff")));

  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIP_FUNCT", "", true)));
  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("funct", {0, false, false, true})));

  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("FORCE_ONLY_TANG_VEL", "", true)));
  xfem_surf_navier_slip->add_component(
      Teuchos::rcp(new Input::IntComponent("force_tang_vel", {0, false, false, true})));

  condlist.push_back(xfem_surf_navier_slip);

  //*----------------*/
  // Surface Navier Slip conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_surf_navier_slip_tpf =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM NAVIER SLIP TWO PHASE SURF CONDITIONS", "XFEMSurfNavierSlipTwoPhase",
          "XFEM Surf Navier Slip", Core::Conditions::XFEM_Surf_Navier_Slip_Twophase, true,
          Core::Conditions::geometry_type_surface));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_surf_navier_slip_tpf->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  xfem_surf_navier_slip_tpf->add_component(Teuchos::rcp(new Input::SeparatorComponent("EVALTYPE")));
  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SelectionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_DIRICHLET_ID")));
  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id_dirch", {0, true, true, false})));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_NEUMANN_ID")));
  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id_neumann", {0, true, true, false})));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIP_SMEAR", "", true)));
  xfem_surf_navier_slip_tpf->add_component(Teuchos::rcp(new Input::RealComponent("slipsmear")));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("NORMAL_PENALTY_SCALING", "", true)));
  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::RealComponent("normalpen_scaling")));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIPCOEFFICIENT")));
  xfem_surf_navier_slip_tpf->add_component(Teuchos::rcp(new Input::RealComponent("slipcoeff")));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("SLIP_FUNCT", "", true)));
  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::IntComponent("funct", {0, false, false, true})));

  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("FORCE_ONLY_TANG_VEL", "", true)));
  xfem_surf_navier_slip_tpf->add_component(
      Teuchos::rcp(new Input::IntComponent("force_tang_vel", {0, false, false, true})));

  condlist.push_back(xfem_surf_navier_slip_tpf);

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_navier_slip_robin_dirch_surf =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM ROBIN DIRICHLET SURF CONDITIONS", "XFEMRobinDirichletSurf",
          "XFEM Robin Dirichlet Volume", Core::Conditions::XFEM_Robin_Dirichlet_Surf, true,
          Core::Conditions::geometry_type_surface));

  // this implementation should be reviewed at some point as it requires these conditions
  //  to have a couplingID. In theory this should not be necessary.
  xfem_navier_slip_robin_dirch_surf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_navier_slip_robin_dirch_surf->add_component(Teuchos::rcp(new Input::IntComponent("label")));

  xfem_navier_slip_robin_dirch_surf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_DIRICHLET_ID")));
  xfem_navier_slip_robin_dirch_surf->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id", {0, true, true, false})));

  // Likely, not necessary. But needed for the current structure.
  xfem_navier_slip_robin_dirch_surf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("EVALTYPE")));
  xfem_navier_slip_robin_dirch_surf->add_component(
      Teuchos::rcp(new Input::SelectionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_navier_slip_robin_dirch_surf->add_component(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_dirch_surf);

  Teuchos::RCP<Core::Conditions::ConditionDefinition> xfem_navier_slip_robin_neumann_surf =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN XFEM ROBIN NEUMANN SURF CONDITIONS", "XFEMRobinNeumannSurf",
          "XFEM Robin Neumann Volume", Core::Conditions::XFEM_Robin_Neumann_Surf, true,
          Core::Conditions::geometry_type_surface));

  // this implementation should be reviewed at some point as it requires these conditions
  //  to have a couplingID. In theory this should not be necessary.
  xfem_navier_slip_robin_neumann_surf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("COUPLINGID")));
  xfem_navier_slip_robin_neumann_surf->add_component(
      Teuchos::rcp(new Input::IntComponent("label")));

  xfem_navier_slip_robin_neumann_surf->add_component(
      Teuchos::rcp(new Input::SeparatorComponent("ROBIN_NEUMANN_ID")));
  xfem_navier_slip_robin_neumann_surf->add_component(
      Teuchos::rcp(new Input::IntComponent("robin_id", {0, true, true, false})));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_navier_slip_robin_neumann_surf->add_component(neumanncomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_neumann_surf);
}

FOUR_C_NAMESPACE_CLOSE
