/*----------------------------------------------------------------------*/
/*!
\file inpar_xfem.cpp

\brief Input parameters for XFEM

\level 2

<pre>
\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_xfem.H"
#include "inpar_twophase.H"
#include "inpar_cut.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::XFEM::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL", false, "");

  // OUTPUT options
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", yesnotuple, yesnovalue,
      &xfem_general);
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT_SCREEN", "No",
      "Do you want to be informed, if Gmsh output is written?", yesnotuple, yesnovalue,
      &xfem_general);
  setStringToIntegralParameter<int>("GMSH_SOL_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", yesnotuple, yesnovalue,
      &xfem_general);
  setStringToIntegralParameter<int>("GMSH_TIMINT_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", yesnotuple, yesnovalue,
      &xfem_general);
  setStringToIntegralParameter<int>("GMSH_EOS_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", yesnotuple, yesnovalue,
      &xfem_general);
  setStringToIntegralParameter<int>("GMSH_DISCRET_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", yesnotuple, yesnovalue,
      &xfem_general);
  setStringToIntegralParameter<int>("GMSH_CUT_OUT", "No",
      "Do you want to write extended Gmsh output for each timestep?", yesnotuple, yesnovalue,
      &xfem_general);


  IntParameter(
      "MAX_NUM_DOFSETS", 3, "Maximum number of volumecells in the XFEM element", &xfem_general);


  setStringToIntegralParameter<int>("NODAL_DOFSET_STRATEGY", "full",
      "which strategy used for the nodal dofset management per node?",
      tuple<std::string>(
          "OneDofset_PerNodeAndPosition", "ConnectGhostDofsets_PerNodeAndPosition", "full"),
      tuple<int>(INPAR::CUT::NDS_Strategy_OneDofset_PerNodeAndPosition,
          INPAR::CUT::NDS_Strategy_ConnectGhostDofsets_PerNodeAndPosition,
          INPAR::CUT::NDS_Strategy_full),
      &xfem_general);


  // Integration options
  setStringToIntegralParameter<int>("VOLUME_GAUSS_POINTS_BY", "Tessellation",
      "how to find Gauss Points for the cut volumes",
      tuple<std::string>("Tessellation", "MomentFitting", "DirectDivergence"),
      tuple<int>(INPAR::CUT::VCellGaussPts_Tessellation, INPAR::CUT::VCellGaussPts_MomentFitting,
          INPAR::CUT::VCellGaussPts_DirectDivergence),
      &xfem_general);

  setStringToIntegralParameter<int>("BOUNDARY_GAUSS_POINTS_BY", "Tessellation",
      "how to find Gauss Points for the boundary cells",
      tuple<std::string>("Tessellation", "MomentFitting", "DirectDivergence"),
      tuple<int>(INPAR::CUT::BCellGaussPts_Tessellation, INPAR::CUT::BCellGaussPts_MomentFitting,
          INPAR::CUT::BCellGaussPts_DirectDivergence),
      &xfem_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_dyn = list->sublist("XFLUID DYNAMIC", false, "");

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_general = xfluid_dyn.sublist("GENERAL", false, "");

  // Do we use more than one fluid discretization?
  BoolParameter("XFLUIDFLUID", "no", "Use an embedded fluid patch.", &xfluid_general);

  // How many monolithic steps we keep the fluidfluid-boundary fixed
  IntParameter(
      "RELAXING_ALE_EVERY", 1, "Relaxing Ale after how many monolithic steps", &xfluid_general);

  BoolParameter("RELAXING_ALE", "yes",
      "switch on/off for relaxing Ale in monolithic fluid-fluid-fsi", &xfluid_general);

  DoubleParameter("XFLUIDFLUID_SEARCHRADIUS", 1.0, "Radius of the search tree", &xfluid_general);

  // xfluidfluid-fsi-monolithic approach
  setStringToIntegralParameter<int>("MONOLITHIC_XFFSI_APPROACH", "xffsi_fixedALE_partitioned",
      "The monolithic apporach for xfluidfluid-fsi",
      tuple<std::string>(
          "xffsi_full_newton", "xffsi_fixedALE_interpolation", "xffsi_fixedALE_partitioned"),
      tuple<int>(INPAR::XFEM::XFFSI_Full_Newton,      // xffsi with no fixed xfem-coupling
          INPAR::XFEM::XFFSI_FixedALE_Interpolation,  // xffsi with fixed xfem-coupling in every
                                                      // newtonstep and interpolations for
                                                      // embedded-dis afterwards
          INPAR::XFEM::XFFSI_FixedALE_Partitioned     // xffsi with fixed xfem-coupling in every
                                                      // newtonstep and solving fluid-field again
          ),
      &xfluid_general);

  // xfluidfluid time integration approach
  setStringToIntegralParameter<int>("XFLUIDFLUID_TIMEINT", "Xff_TimeInt_FullProj",
      "The xfluidfluid-timeintegration approach",
      tuple<std::string>("Xff_TimeInt_FullProj", "Xff_TimeInt_ProjIfMoved",
          "Xff_TimeInt_KeepGhostValues", "Xff_TimeInt_IncompProj"),
      tuple<int>(INPAR::XFEM::Xff_TimeInt_FullProj,  // always project nodes from embedded to
                                                     // background nodes
          INPAR::XFEM::Xff_TimeInt_ProjIfMoved,  // project nodes just if the status of background
                                                 // nodes changed
          INPAR::XFEM::Xff_TimeInt_KeepGhostValues,  // always keep the ghost values of the
                                                     // background discretization
          INPAR::XFEM::Xff_TimeInt_IncompProj  // after projecting nodes do a incompressibility
                                               // projection
          ),
      &xfluid_general);

  setStringToIntegralParameter<int>("XFLUID_TIMEINT", "STD=COPY/SL_and_GHOST=COPY/GP",
      "The xfluid time integration approach",
      tuple<std::string>("STD=COPY_and_GHOST=COPY/GP", "STD=COPY/SL_and_GHOST=COPY/GP",
          "STD=SL(boundary-zone)_and_GHOST=GP", "STD=COPY/PROJ_and_GHOST=COPY/PROJ/GP"),
      tuple<int>(
          INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP,  // STD= only copy,
                                                                              // GHOST= copy or
                                                                              // ghost penalty
                                                                              // reconstruction
          INPAR::XFEM::
              Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP,  // STD= copy or SL,
                                                                           // GHOST= copy or ghost
                                                                           // penalty reconstruction
          INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP,  // STD= only SL on
                                                                             // whole boundary zone,
                                                                             // GHOST= ghost penalty
                                                                             // reconstruction
          INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP),
      &xfluid_general);

  BoolParameter("ALE_XFluid", "no", "XFluid is Ale Fluid?", &xfluid_general);

  // for new OST-implementation: which interface terms to be evaluated for previous time step
  setStringToIntegralParameter<int>("INTERFACE_TERMS_PREVIOUS_STATE",
      "PreviousState_only_consistency",
      "how to treat interface terms from previous time step (new OST)",
      tuple<std::string>("PreviousState_only_consistency", "PreviousState_full"),
      tuple<int>(INPAR::XFEM::PreviousState_only_consistency,  /// evaluate only consistency terms
                                                               /// for previous time step
          INPAR::XFEM::PreviousState_full  /// evaluate consistency, adjoint consistency and penalty
                                           /// terms or previous time step
          ),
      &xfluid_general);

  setStringToIntegralParameter<int>("XFLUID_TIMEINT_CHECK_INTERFACETIPS", "Yes",
      "Xfluid TimeIntegration Special Check if node has changed the side!", yesnotuple, yesnovalue,
      &xfluid_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_stab = xfluid_dyn.sublist("STABILIZATION", false, "");

  // Boundary-Coupling options
  setStringToIntegralParameter<int>("COUPLING_METHOD", "Nitsche",
      "method how to enforce embedded boundary/coupling conditions at the interface",
      tuple<std::string>("Hybrid_LM_Cauchy_stress", "Hybrid_LM_viscous_stress", "Nitsche"),
      tuple<int>(
          INPAR::XFEM::Hybrid_LM_Cauchy_stress,   // Cauchy stress-based mixed/hybrid formulation
          INPAR::XFEM::Hybrid_LM_viscous_stress,  // viscous stress-based mixed/hybrid formulation
          INPAR::XFEM::Nitsche                    // Nitsche's formulation
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("HYBRID_LM_L2_PROJ", "part_ele_proj",
      "perform the L2 projection between stress fields on whole element or on fluid part?",
      tuple<std::string>("full_ele_proj", "part_ele_proj"),
      tuple<int>(
          INPAR::XFEM::Hybrid_LM_L2_Proj_full,  // L2 stress projection on whole fluid element
          INPAR::XFEM::Hybrid_LM_L2_Proj_part   // L2 stress projection on partial fluid element
                                                // volume
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_ADJOINT_SYMMETRY", "yes",
      "viscous and adjoint viscous interface terms with matching sign?",
      tuple<std::string>("yes", "no", "sym", "skew", "none"),
      tuple<int>(INPAR::XFEM::adj_sym, INPAR::XFEM::adj_skew, INPAR::XFEM::adj_sym,
          INPAR::XFEM::adj_skew, INPAR::XFEM::adj_none),
      &xfluid_stab);

  // viscous and convective Nitsche/MSH stabilization parameter
  DoubleParameter(
      "NIT_STAB_FAC", 35.0, " ( stabilization parameter for Nitsche's penalty term", &xfluid_stab);
  DoubleParameter("NIT_STAB_FAC_TANG", 35.0,
      " ( stabilization parameter for Nitsche's penalty tangential term", &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_TRACE_ESTIMATE", "CT_div_by_hk",
      "how to estimate the scaling from the trace inequality in Nitsche's method",
      tuple<std::string>("CT_div_by_hk", "eigenvalue"),
      tuple<int>(
          INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk,  // estimate the trace inequality by a
                                                             // trace-constant CT and hk: CT/hk
          INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue     // estimate the trace inequality by
                                                          // solving a local element-wise eigenvalue
                                                          // problem
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("UPDATE_EIGENVALUE_TRACE_ESTIMATE", "every_iter",
      "how often should the local eigenvalue problem be updated",
      tuple<std::string>("every_iter", "every_timestep", "once"),
      tuple<int>(INPAR::XFEM::Eigenvalue_update_every_iter,
          INPAR::XFEM::Eigenvalue_update_every_timestep, INPAR::XFEM::Eigenvalue_update_once),
      &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_HK", "ele_vol_div_by_max_ele_surf",
      "how to define the characteristic element length in cut elements",
      tuple<std::string>("vol_equivalent", "cut_vol_div_by_cut_surf", "ele_vol_div_by_cut_surf",
          "ele_vol_div_by_ele_surf", "ele_vol_div_by_max_ele_surf"),
      tuple<int>(INPAR::XFEM::ViscStab_hk_vol_equivalent,    /// volume equivalent element diameter
          INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf,  /// physical partial/cut volume divided
                                                             /// by physical partial/cut surface
                                                             /// measure ( used to estimate the
                                                             /// cut-dependent inverse estimate on
                                                             /// cut elements, not useful for sliver
                                                             /// and/or dotted cut situations)
          INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf,  /// full element volume divided by
                                                             /// physical partial/cut surface
                                                             /// measure ( used to estimate the
                                                             /// cut-dependent inverse estimate on
                                                             /// cut elements, however, avoids
                                                             /// problems with sliver cuts, not
                                                             /// useful for dotted cuts)
          INPAR::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf,  /// full element volume divided by
                                                             /// surface measure ( used for uncut
                                                             /// situations, standard weak Dirichlet
                                                             /// boundary/coupling conditions)
          INPAR::XFEM::
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
      tuple<int>(INPAR::XFEM::ConvStabScaling_inflow,  // scaling with max(0,-u*n)
          INPAR::XFEM::ConvStabScaling_abs_inflow,     // scaling with |u*n|
          INPAR::XFEM::ConvStabScaling_none            // no convective stabilization
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("XFF_CONV_STAB_SCALING", "none",
      "scaling factor for convective interface stabilization of fluid-fluid Coupling",
      tuple<std::string>("inflow", "averaged", "none"),
      tuple<int>(INPAR::XFEM::XFF_ConvStabScaling_upwinding,  // one-sided inflow stabilization
          INPAR::XFEM::XFF_ConvStabScaling_only_averaged,     // averaged inflow stabilization
          INPAR::XFEM::XFF_ConvStabScaling_none               // no convective stabilization
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_COMBO", "max",
      "choose the maximum from viscous and convective contributions or just sum both up",
      tuple<std::string>("max", "sum"),
      tuple<int>(INPAR::XFEM::MassConservationCombination_max,  /// use the maximum contribution
          INPAR::XFEM::MassConservationCombination_sum  /// sum viscous and convective contributions
          ),
      &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_SCALING", "only_visc",
      "apply additional scaling of penalty term to enforce mass conservation for "
      "convection-dominated flow",
      tuple<std::string>("full", "only_visc"),
      tuple<int>(INPAR::XFEM::MassConservationScaling_full,  /// apply mass-conserving convective
                                                             /// scaling additionally
          INPAR::XFEM::MassConservationScaling_only_visc     /// use only the viscous scaling
          ),
      &xfluid_stab);

  BoolParameter("GHOST_PENALTY_STAB", "no", "switch on/off ghost penalty interface stabilization",
      &xfluid_stab);

  BoolParameter("GHOST_PENALTY_TRANSIENT_STAB", "no",
      "switch on/off ghost penalty transient interface stabilization", &xfluid_stab);

  BoolParameter("GHOST_PENALTY_2nd_STAB", "no",
      "switch on/off ghost penalty interface stabilization for 2nd order derivatives",
      &xfluid_stab);
  BoolParameter("GHOST_PENALTY_2nd_STAB_NORMAL", "no",
      "switch between ghost penalty interface stabilization for 2nd order derivatives in normal or "
      "all spatial directions",
      &xfluid_stab);


  DoubleParameter("GHOST_PENALTY_FAC", 0.1,
      "define stabilization parameter ghost penalty interface stabilization", &xfluid_stab);

  DoubleParameter("GHOST_PENALTY_TRANSIENT_FAC", 0.001,
      "define stabilization parameter ghost penalty transient interface stabilization",
      &xfluid_stab);

  //  // NOT chosen optimally!
  //  DoubleParameter("GHOST_PENALTY_2nd_FAC", 1.0,"define stabilization parameter ghost penalty 2nd
  //  order viscous interface stabilization",&xfluid_stab);
  //  DoubleParameter("GHOST_PENALTY_PRESSURE_2nd_FAC", 1.0,"define stabilization parameter ghost
  //  penalty 2nd order pressure interface stabilization",&xfluid_stab);

  DoubleParameter("GHOST_PENALTY_2nd_FAC", 0.05,
      "define stabilization parameter ghost penalty 2nd order viscous interface stabilization",
      &xfluid_stab);
  DoubleParameter("GHOST_PENALTY_PRESSURE_2nd_FAC", 0.05,
      "define stabilization parameter ghost penalty 2nd order pressure interface stabilization",
      &xfluid_stab);


  BoolParameter("XFF_EOS_PRES_EMB_LAYER", "no",
      "switch on/off edge-based pressure stabilization on interface-contributing elements of the "
      "embedded fluid",
      &xfluid_stab);

  BoolParameter("IS_PSEUDO_2D", "no",
      "modify viscous interface stabilization due to the vanishing polynomial in third dimension "
      "when using strong Dirichlet conditions to block polynomials in one spatial dimension",
      &xfluid_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfsi_monolithic = xfluid_dyn.sublist("XFPSI MONOLITHIC", false, "");

  IntParameter("ITEMIN", 1, "How many iterations are performed minimal", &xfsi_monolithic);
  IntParameter(
      "ITEMAX_OUTER", 5, "How many outer iterations are performed maximal", &xfsi_monolithic);
  BoolParameter("ND_NEWTON_DAMPING", "no",
      "Activate Newton damping based on residual and increment", &xfsi_monolithic);
  DoubleParameter("ND_MAX_DISP_ITERINC", -1.0,
      "Maximal displacement increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  DoubleParameter("ND_MAX_VEL_ITERINC", -1.0,
      "Maximal fluid velocity increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  DoubleParameter("ND_MAX_PRES_ITERINC", -1.0,
      "Maximal fluid pressure increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  DoubleParameter("ND_MAX_PVEL_ITERINC", -1.0,
      "Maximal porofluid velocity increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
  DoubleParameter("ND_MAX_PPRES_ITERINC", -1.0,
      "Maximal porofluid pressure increment to apply full newton --> otherwise damp newton",
      &xfsi_monolithic);
}



void INPAR::XFEM::SetValidConditions(
    const std::vector<Teuchos::RCP<DRT::INPUT::ConditionComponent>>& dirichletbundcomponents,
    const std::vector<Teuchos::RCP<DRT::INPUT::ConditionComponent>>& neumanncomponents,
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  std::vector<Teuchos::RCP<ConditionComponent>> xfemcomponents;

  xfemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("label")));

  Teuchos::RCP<ConditionDefinition> movingfluid =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID MESH VOL CONDITIONS", "FluidMesh",
          "Fluid Mesh", DRT::Condition::FluidMesh, true, DRT::Condition::Volume));
  Teuchos::RCP<ConditionDefinition> fluidfluidcoupling = Teuchos::rcp(new ConditionDefinition(
      "DESIGN FLUID FLUID COUPLING SURF CONDITIONS", "FluidFluidCoupling", "FLUID FLUID Coupling",
      DRT::Condition::FluidFluidCoupling, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> ALEfluidcoupling = Teuchos::rcp(
      new ConditionDefinition("DESIGN ALE FLUID COUPLING SURF CONDITIONS", "ALEFluidCoupling",
          "ALE FLUID Coupling", DRT::Condition::ALEFluidCoupling, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < xfemcomponents.size(); ++i)
  {
    movingfluid->AddComponent(xfemcomponents[i]);
    fluidfluidcoupling->AddComponent(xfemcomponents[i]);
    ALEfluidcoupling->AddComponent(xfemcomponents[i]);
  }

  condlist.push_back(fluidfluidcoupling);
  condlist.push_back(movingfluid);
  condlist.push_back(ALEfluidcoupling);


  /*--------------------------------------------------------------------*/
  // XFEM coupling conditions

  //*----------------*/
  // Displacement surface condition for XFEM WDBC and Neumann boundary conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_displacement = Teuchos::rcp(new ConditionDefinition(
      "DESIGN XFEM DISPLACEMENT SURF CONDITIONS", "XFEMSurfDisplacement", "XFEM Surf Displacement",
      DRT::Condition::XFEM_Surf_Displacement, true, DRT::Condition::Surface));

  xfem_surf_displacement->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_displacement->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));
  xfem_surf_displacement->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_surf_displacement->AddComponent(Teuchos::rcp(new StringConditionComponent("evaltype",
      "funct", Teuchos::tuple<std::string>("zero", "funct", "implementation"),
      Teuchos::tuple<std::string>("zero", "funct", "implementation"), true)));


  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_surf_displacement->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_surf_displacement);



  //*----------------*/
  // Levelset field condition components

  std::vector<Teuchos::RCP<ConditionComponent>> levelsetfield_components;

  levelsetfield_components.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  levelsetfield_components.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  levelsetfield_components.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("LEVELSETFIELDNO")));
  levelsetfield_components.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("levelsetfieldno", 1, false, false, false)));

  // define which boolean operator is used for combining this level-set field with the previous one
  // with smaller coupling id
  levelsetfield_components.push_back(Teuchos::rcp(new SeparatorConditionComponent("BOOLEANTYPE")));
  levelsetfield_components.push_back(Teuchos::rcp(new StringConditionComponent("booleantype",
      "none", Teuchos::tuple<std::string>("none", "cut", "union", "difference", "sym_difference"),
      Teuchos::tuple<std::string>("none", "cut", "union", "difference", "sym_difference"), false)));

  // define which complementary operator is applied after combining the level-set field with a
  // boolean operator with the previous one
  levelsetfield_components.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("COMPLEMENTARY")));
  levelsetfield_components.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("complementary", 1, false, false, false)));


  //*----------------*/
  // Levelset based Weak Dirichlet conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_wdbc =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM LEVELSET WEAK DIRICHLET VOL CONDITIONS",
          "XFEMLevelsetWeakDirichlet", "XFEM Levelset Weak Dirichlet",
          DRT::Condition::XFEM_Levelset_Weak_Dirichlet, true, DRT::Condition::Volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_wdbc->AddComponent(levelsetfield_components[i]);
  }

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_levelset_wdbc->AddComponent(dirichletbundcomponents[i]);
  }

  // optional: allow for random noise, set percentage used in uniform random distribution
  xfem_levelset_wdbc->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("RANDNOISE", true)));
  xfem_levelset_wdbc->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("randnoise", 1, true)));

  condlist.push_back(xfem_levelset_wdbc);

  //*----------------*/
  // Levelset based Neumann conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_neumann = Teuchos::rcp(new ConditionDefinition(
      "DESIGN XFEM LEVELSET NEUMANN VOL CONDITIONS", "XFEMLevelsetNeumann", "XFEM Levelset Neumann",
      DRT::Condition::XFEM_Levelset_Neumann, true, DRT::Condition::Volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_neumann->AddComponent(levelsetfield_components[i]);
  }

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_levelset_neumann->AddComponent(neumanncomponents[i]);
  }

  // define if we use inflow stabilization on the xfem neumann surf condition
  xfem_levelset_neumann->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("INFLOW_STAB", true)));
  xfem_levelset_neumann->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("InflowStab", 1, true, true, true)));

  condlist.push_back(xfem_levelset_neumann);

  //*----------------*/
  // Levelset based Navier Slip conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_navier_slip =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM LEVELSET NAVIER SLIP VOL CONDITIONS",
          "XFEMLevelsetNavierSlip", "XFEM Levelset Navier Slip",
          DRT::Condition::XFEM_Levelset_Navier_Slip, true, DRT::Condition::Volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_navier_slip->AddComponent(levelsetfield_components[i]);
  }

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SURFACE_PROJECTION")));
  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new StringConditionComponent("SURFACE_PROJECTION", "proj_normal",
          Teuchos::tuple<std::string>(
              "proj_normal", "proj_smoothed", "proj_normal_smoothed_comb", "proj_normal_phi"),
          Teuchos::tuple<int>(INPAR::XFEM::Proj_normal, INPAR::XFEM::Proj_smoothed,
              INPAR::XFEM::Proj_normal_smoothed_comb, INPAR::XFEM::Proj_normal_phi),
          true)));

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("L2_PROJECTION_SOLVER")));
  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("l2projsolv", 1, true, true, false)));

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_DIRICHLET_ID")));
  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id_dirch", 1, true, true, false)));

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_NEUMANN_ID")));
  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id_neumann", 1, true, true, false)));

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT")));
  xfem_levelset_navier_slip->AddComponent(Teuchos::rcp(new RealConditionComponent("slipcoeff")));

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIP_FUNCT", true)));
  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("FORCE_ONLY_TANG_VEL", true)));
  xfem_levelset_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("force_tang_vel", 1, false, false, true)));

  condlist.push_back(xfem_levelset_navier_slip);

  // Add condition XFEM DIRICHLET/NEUMANN?

  Teuchos::RCP<ConditionDefinition> xfem_navier_slip_robin_dirch =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM ROBIN DIRICHLET VOL CONDITIONS",
          "XFEMRobinDirichletVol", "XFEM Robin Dirichlet Volume",
          DRT::Condition::XFEM_Robin_Dirichlet_Volume, true, DRT::Condition::Volume));

  xfem_navier_slip_robin_dirch->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_DIRICHLET_ID")));
  xfem_navier_slip_robin_dirch->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id", 1, true, true, false)));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_navier_slip_robin_dirch->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_dirch);

  Teuchos::RCP<ConditionDefinition> xfem_navier_slip_robin_neumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM ROBIN NEUMANN VOL CONDITIONS",
          "XFEMRobinNeumannVol", "XFEM Robin Neumann Volume",
          DRT::Condition::XFEM_Robin_Neumann_Volume, true, DRT::Condition::Volume));

  xfem_navier_slip_robin_neumann->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_NEUMANN_ID")));
  xfem_navier_slip_robin_neumann->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id", 1, true, true, false)));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_navier_slip_robin_neumann->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_neumann);


  //*----------------*/
  // Levelset based Twophase conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_twophase =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS",
          "XFEMLevelsetTwophase", "XFEM Levelset Twophase", DRT::Condition::XFEM_Levelset_Twophase,
          true, DRT::Condition::Volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_twophase->AddComponent(levelsetfield_components[i]);
  }

  condlist.push_back(xfem_levelset_twophase);

  //*----------------*/
  // Levelset based Combustion conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_combustion =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM LEVELSET COMBUSTION VOL CONDITIONS",
          "XFEMLevelsetCombustion", "XFEM Levelset Combustion",
          DRT::Condition::XFEM_Levelset_Combustion, true, DRT::Condition::Volume));

  for (unsigned i = 0; i < levelsetfield_components.size(); ++i)
  {
    xfem_levelset_combustion->AddComponent(levelsetfield_components[i]);
  }

  // "The laminar flamespeed incorporates all chemical kinetics into the problem for now"
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("LAMINAR_FLAMESPEED", false)));
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("laminar_flamespeed", 1, false)));

  // "Molecular diffusivity"
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("MOL_DIFFUSIVITY", false)));
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("mol_diffusivity", 1, false)));

  // "The Markstein length takes flame curvature into account"
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("MARKSTEIN_LENGTH", false)));
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("markstein_length", 1, false)));

  // interface transport in all directions or just in normal direction?
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("TRANSPORT_DIRECTIONS")));
  xfem_levelset_combustion->AddComponent(Teuchos::rcp(new StringConditionComponent(
      "TRANSPORT_DIRECTIONS", "all", Teuchos::tuple<std::string>("all", "normal"),
      Teuchos::tuple<int>(
          INPAR::TWOPHASE::transport_dir_all, INPAR::TWOPHASE::transport_dir_normal),
      true)));


  // define if the curvature shall be accounted for in computing the transport velocity
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("TRANSPORT_CURVATURE")));
  xfem_levelset_combustion->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("transport_curvature", 1, false, false, false)));


  condlist.push_back(xfem_levelset_combustion);

  //*----------------*/
  // Surface Fluid-Fluid coupling conditions


  std::vector<Teuchos::RCP<ConditionComponent>> xfluidfluidsurfcomponents;

  xfluidfluidsurfcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfluidfluidsurfcomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfluidfluidsurfcomponents.push_back(Teuchos::rcp(new StringConditionComponent("COUPSTRATEGY",
      "xfluid", Teuchos::tuple<std::string>("xfluid", "embedded", "mean"),
      Teuchos::tuple<int>(
          INPAR::XFEM::Xfluid_Sided, INPAR::XFEM::Embedded_Sided, INPAR::XFEM::Mean))));


  Teuchos::RCP<ConditionDefinition> xfem_surf_fluidfluid = Teuchos::rcp(new ConditionDefinition(
      "DESIGN XFEM FLUIDFLUID SURF CONDITIONS", "XFEMSurfFluidFluid", "XFEM Surf FluidFluid",
      DRT::Condition::XFEM_Surf_FluidFluid, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < xfluidfluidsurfcomponents.size(); ++i)
    xfem_surf_fluidfluid->AddComponent(xfluidfluidsurfcomponents[i]);

  condlist.push_back(xfem_surf_fluidfluid);

  //*----------------*/
  // Surface partitioned XFSI boundary conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_fsi_part = Teuchos::rcp(
      new ConditionDefinition("DESIGN XFEM FSI PARTITIONED SURF CONDITIONS", "XFEMSurfFSIPart",
          "XFEM Surf FSI Part", DRT::Condition::XFEM_Surf_FSIPart, true, DRT::Condition::Surface));

  xfem_surf_fsi_part->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_fsi_part->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  // COUPSTRATEGY IS FLUID SIDED

  xfem_surf_fsi_part->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("INTLAW", true)));
  xfem_surf_fsi_part->AddComponent(Teuchos::rcp(new StringConditionComponent("INTLAW", "noslip",
      Teuchos::tuple<std::string>("noslip", "noslip_splitpen", "slip", "navslip"),
      Teuchos::tuple<int>(INPAR::XFEM::noslip, INPAR::XFEM::noslip_splitpen, INPAR::XFEM::slip,
          INPAR::XFEM::navierslip),
      true)));

  xfem_surf_fsi_part->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT", true)));
  xfem_surf_fsi_part->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("slipcoeff", 1, true)));

  xfem_surf_fsi_part->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIP_FUNCT", true)));
  xfem_surf_fsi_part->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  condlist.push_back(xfem_surf_fsi_part);

  //*----------------*/
  // Surface monolithic XFSI coupling conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_fsi_mono = Teuchos::rcp(
      new ConditionDefinition("DESIGN XFEM FSI MONOLITHIC SURF CONDITIONS", "XFEMSurfFSIMono",
          "XFEM Surf FSI Mono", DRT::Condition::XFEM_Surf_FSIMono, true, DRT::Condition::Surface));

  xfem_surf_fsi_mono->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_fsi_mono->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfem_surf_fsi_mono->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("COUPSTRATEGY", true)));
  xfem_surf_fsi_mono->AddComponent(Teuchos::rcp(new StringConditionComponent("COUPSTRATEGY",
      "xfluid", Teuchos::tuple<std::string>("xfluid", "solid", "mean", "harmonic"),
      Teuchos::tuple<int>(INPAR::XFEM::Xfluid_Sided, INPAR::XFEM::Embedded_Sided, INPAR::XFEM::Mean,
          INPAR::XFEM::Harmonic),
      true)));

  xfem_surf_fsi_mono->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("INTLAW", true)));
  xfem_surf_fsi_mono->AddComponent(Teuchos::rcp(new StringConditionComponent("INTLAW", "noslip",
      Teuchos::tuple<std::string>("noslip", "noslip_splitpen", "slip", "navslip"),
      Teuchos::tuple<int>(INPAR::XFEM::noslip, INPAR::XFEM::noslip_splitpen, INPAR::XFEM::slip,
          INPAR::XFEM::navierslip),
      true)));

  xfem_surf_fsi_mono->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT", true)));
  xfem_surf_fsi_mono->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("slipcoeff", 1, true)));

  xfem_surf_fsi_mono->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIP_FUNCT", true)));
  xfem_surf_fsi_mono->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  condlist.push_back(xfem_surf_fsi_mono);

  //*----------------*/
  // Surface monolithic XFPI coupling conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_fpi_mono = Teuchos::rcp(
      new ConditionDefinition("DESIGN XFEM FPI MONOLITHIC SURF CONDITIONS", "XFEMSurfFPIMono",
          "XFEM Surf FPI Mono", DRT::Condition::XFEM_Surf_FPIMono, true, DRT::Condition::Surface));

  xfem_surf_fpi_mono->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_fpi_mono->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfem_surf_fpi_mono->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("BJ_COEFF", true)));
  xfem_surf_fpi_mono->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("bj_coeff", 1, true)));

  xfem_surf_fpi_mono->AddComponent(Teuchos::rcp(new StringConditionComponent("Variant", "BJ",
      Teuchos::tuple<std::string>("BJ", "BJS"), Teuchos::tuple<std::string>("BJ", "BJS"), true)));

  xfem_surf_fpi_mono->AddComponent(Teuchos::rcp(new StringConditionComponent("Method", "NIT",
      Teuchos::tuple<std::string>("NIT", "SUB"), Teuchos::tuple<std::string>("NIT", "SUB"), true)));

  condlist.push_back(xfem_surf_fpi_mono);


  //*----------------*/
  // Surface Weak Dirichlet conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_wdbc =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM WEAK DIRICHLET SURF CONDITIONS",
          "XFEMSurfWeakDirichlet", "XFEM Surf Weak Dirichlet",
          DRT::Condition::XFEM_Surf_Weak_Dirichlet, true, DRT::Condition::Surface));

  xfem_surf_wdbc->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_wdbc->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));
  xfem_surf_wdbc->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_surf_wdbc->AddComponent(
      Teuchos::rcp(new StringConditionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_surf_wdbc->AddComponent(dirichletbundcomponents[i]);
  }

  // optional: allow for random noise, set percentage used in uniform random distribution
  xfem_surf_wdbc->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("RANDNOISE", true)));
  xfem_surf_wdbc->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("randnoise", 1, true)));

  condlist.push_back(xfem_surf_wdbc);


  //*----------------*/
  // Surface Neumann conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_neumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM NEUMANN SURF CONDITIONS", "XFEMSurfNeumann",
          "XFEM Surf Neumann", DRT::Condition::XFEM_Surf_Neumann, true, DRT::Condition::Surface));

  xfem_surf_neumann->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_neumann->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_surf_neumann->AddComponent(neumanncomponents[i]);
  }

  // define if we use inflow stabilization on the xfem neumann surf condition
  xfem_surf_neumann->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("INFLOW_STAB", true)));
  xfem_surf_neumann->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("InflowStab", 1, true, true, true)));

  condlist.push_back(xfem_surf_neumann);

  //*----------------*/
  // Surface Navier Slip conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_navier_slip = Teuchos::rcp(new ConditionDefinition(
      "DESIGN XFEM NAVIER SLIP SURF CONDITIONS", "XFEMSurfNavierSlip", "XFEM Surf Navier Slip",
      DRT::Condition::XFEM_Surf_Navier_Slip, true, DRT::Condition::Surface));

  xfem_surf_navier_slip->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfem_surf_navier_slip->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new StringConditionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));


  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_DIRICHLET_ID")));
  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id_dirch", 1, true, true, false)));

  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_NEUMANN_ID")));
  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id_neumann", 1, true, true, false)));

  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT")));
  xfem_surf_navier_slip->AddComponent(Teuchos::rcp(new RealConditionComponent("slipcoeff")));

  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIP_FUNCT", true)));
  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("FORCE_ONLY_TANG_VEL", true)));
  xfem_surf_navier_slip->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("force_tang_vel", 1, false, false, true)));

  condlist.push_back(xfem_surf_navier_slip);

  //*----------------*/
  // Surface Navier Slip conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_navier_slip_tpf =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM NAVIER SLIP TWO PHASE SURF CONDITIONS",
          "XFEMSurfNavierSlipTwoPhase", "XFEM Surf Navier Slip",
          DRT::Condition::XFEM_Surf_Navier_Slip_Twophase, true, DRT::Condition::Surface));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new StringConditionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_DIRICHLET_ID")));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id_dirch", 1, true, true, false)));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_NEUMANN_ID")));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id_neumann", 1, true, true, false)));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIP_SMEAR", true)));
  xfem_surf_navier_slip_tpf->AddComponent(Teuchos::rcp(new RealConditionComponent("slipsmear")));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("NORMAL_PENALTY_SCALING", true)));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new RealConditionComponent("normalpen_scaling")));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT")));
  xfem_surf_navier_slip_tpf->AddComponent(Teuchos::rcp(new RealConditionComponent("slipcoeff")));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("SLIP_FUNCT", true)));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("FORCE_ONLY_TANG_VEL", true)));
  xfem_surf_navier_slip_tpf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("force_tang_vel", 1, false, false, true)));

  condlist.push_back(xfem_surf_navier_slip_tpf);

  Teuchos::RCP<ConditionDefinition> xfem_navier_slip_robin_dirch_surf =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM ROBIN DIRICHLET SURF CONDITIONS",
          "XFEMRobinDirichletSurf", "XFEM Robin Dirichlet Volume",
          DRT::Condition::XFEM_Robin_Dirichlet_Surf, true, DRT::Condition::Surface));

  // this implementation should be reviewed at some point as it requires these conditions
  //  to have a couplingID. In theory this should not be necessary.
  xfem_navier_slip_robin_dirch_surf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_navier_slip_robin_dirch_surf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfem_navier_slip_robin_dirch_surf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_DIRICHLET_ID")));
  xfem_navier_slip_robin_dirch_surf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id", 1, true, true, false)));

  // Likely, not necessary. But needed for the current structure.
  xfem_navier_slip_robin_dirch_surf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_navier_slip_robin_dirch_surf->AddComponent(
      Teuchos::rcp(new StringConditionComponent("evaltype", "funct_interpolated",
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          Teuchos::tuple<std::string>("zero", "funct_interpolated", "funct_gausspoint",
              "displacement_1storder_wo_initfunct", "displacement_2ndorder_wo_initfunct",
              "displacement_1storder_with_initfunct", "displacement_2ndorder_with_initfunct"),
          true)));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    xfem_navier_slip_robin_dirch_surf->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_dirch_surf);

  Teuchos::RCP<ConditionDefinition> xfem_navier_slip_robin_neumann_surf =
      Teuchos::rcp(new ConditionDefinition("DESIGN XFEM ROBIN NEUMANN SURF CONDITIONS",
          "XFEMRobinNeumannSurf", "XFEM Robin Neumann Volume",
          DRT::Condition::XFEM_Robin_Neumann_Surf, true, DRT::Condition::Surface));

  // this implementation should be reviewed at some point as it requires these conditions
  //  to have a couplingID. In theory this should not be necessary.
  xfem_navier_slip_robin_neumann_surf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_navier_slip_robin_neumann_surf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("label", 1, false, false, false)));

  xfem_navier_slip_robin_neumann_surf->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("ROBIN_NEUMANN_ID")));
  xfem_navier_slip_robin_neumann_surf->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("robin_id", 1, true, true, false)));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    xfem_navier_slip_robin_neumann_surf->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(xfem_navier_slip_robin_neumann_surf);
}
