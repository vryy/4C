/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid conditions for input

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_inpar_validconditions.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_inpar_beampotential.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_inpar_cardiac_monodomain.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_inpar_ehl.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_inpar_elemag.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_immersed.hpp"
#include "4C_inpar_levelset.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_inpar_mpc_rve.hpp"
#include "4C_inpar_particle.hpp"
#include "4C_inpar_poromultiphase_scatra.hpp"
#include "4C_inpar_s2i.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_inpar_ssti.hpp"
#include "4C_inpar_sti.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_inpar_thermo.hpp"
#include "4C_inpar_xfem.hpp"
#include "4C_io_linecomponent.hpp"

FOUR_C_NAMESPACE_OPEN


void Input::print_empty_condition_definitions(std::ostream& stream,
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  for (unsigned i = 0; i < condlist.size(); ++i)
  {
    condlist[i]->print(stream, nullptr);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void print_condition_dat_header()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> condlist =
      Input::valid_conditions();
  Input::print_empty_condition_definitions(std::cout, *condlist);
}


namespace Input
{
  // collect some problem-specific conditions that do not fit in the generic sections
  void set_miscellaneous_conditions(
      std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
  {
    /*--------------------------------------------------------------------*/
    // microscale boundary

    Teuchos::RCP<Core::Conditions::ConditionDefinition> microscale =
        Teuchos::RCP(new Core::Conditions::ConditionDefinition("MICROSCALE CONDITIONS",
            "MicroBoundary", "Microscale Boundary", Core::Conditions::MicroBoundary, true,
            Core::Conditions::geometry_type_surface));

    condlist.push_back(microscale);

    /*--------------------------------------------------------------------*/
    // stc layer condition

    Teuchos::RCP<Core::Conditions::ConditionDefinition> stclayer =
        Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL STC LAYER", "STC Layer",
            "Layer for Multilayered STC", Core::Conditions::VolSTCLayer, true,
            Core::Conditions::geometry_type_volume));

    stclayer->add_component(Teuchos::RCP(new IntComponent("ConditionID")));

    condlist.push_back(stclayer);
  }
}  // namespace Input


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>>
Input::valid_conditions()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> vc =
      Teuchos::RCP(new std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>());

  std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist = *vc;

  /*--------------------------------------------------------------------*/
  // Neumann conditions
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT NEUMANN CONDITIONS",
          "PointNeumann", "Point Neumann", Core::Conditions::PointNeumann, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointneumanneb =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT MOMENT EB CONDITIONS",
          "PointNeumannEB", "Point Neumann Moment for an Euler-Bernoulli beam",
          Core::Conditions::PointNeumannEB, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE NEUMANN CONDITIONS",
          "LineNeumann", "Line Neumann", Core::Conditions::LineNeumann, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF NEUMANN CONDITIONS",
          "SurfaceNeumann", "Surface Neumann", Core::Conditions::SurfaceNeumann, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL NEUMANN CONDITIONS",
          "VolumeNeumann", "Volume Neumann", Core::Conditions::VolumeNeumann, true,
          Core::Conditions::geometry_type_volume));

  // Neumann conditions for transport problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointtransportneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT TRANSPORT NEUMANN CONDITIONS", "TransportPointNeumann", "Point Neumann",
          Core::Conditions::PointNeumann, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linetransportneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE TRANSPORT NEUMANN CONDITIONS", "TransportLineNeumann", "Line Neumann",
          Core::Conditions::LineNeumann, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surftransportneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF TRANSPORT NEUMANN CONDITIONS", "TransportSurfaceNeumann", "Surface Neumann",
          Core::Conditions::SurfaceNeumann, true, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> voltransportneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL TRANSPORT NEUMANN CONDITIONS", "TransportVolumeNeumann", "Volume Neumann",
          Core::Conditions::VolumeNeumann, true, Core::Conditions::geometry_type_volume));

  // Neumann conditions for thermo problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO NEUMANN CONDITIONS", "ThermoPointNeumann", "Point Neumann",
          Core::Conditions::PointNeumann, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linethermoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE THERMO NEUMANN CONDITIONS", "ThermoLineNeumann", "Line Neumann",
          Core::Conditions::LineNeumann, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfthermoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF THERMO NEUMANN CONDITIONS", "ThermoSurfaceNeumann", "Surface Neumann",
          Core::Conditions::SurfaceNeumann, true, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volthermoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL THERMO NEUMANN CONDITIONS",
          "ThermoVolumeNeumann", "Volume Neumann", Core::Conditions::VolumeNeumann, true,
          Core::Conditions::geometry_type_volume));

  // Neumann conditions for poroelasticity problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointporoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT PORO NEUMANN CONDITIONS",
          "PoroPointNeumann", "Point Neumann", Core::Conditions::PointNeumann, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineporoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE PORO NEUMANN CONDITIONS",
          "PoroLineNeumann", "Line Neumann", Core::Conditions::LineNeumann, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfporoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF PORO NEUMANN CONDITIONS",
          "PoroSurfaceNeumann", "Surface Neumann", Core::Conditions::SurfaceNeumann, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volporoneumann =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL PORO NEUMANN CONDITIONS",
          "PoroVolumeNeumann", "Volume Neumann", Core::Conditions::VolumeNeumann, true,
          Core::Conditions::geometry_type_volume));

  auto all_neumann_conditions = {pointneumann, pointneumanneb, lineneumann, surfneumann, volneumann,
      pointtransportneumann, linetransportneumann, surftransportneumann, voltransportneumann,
      pointthermoneumann, linethermoneumann, surfthermoneumann, volthermoneumann, pointporoneumann,
      lineporoneumann, surfporoneumann, volporoneumann};

  for (const auto& cond : all_neumann_conditions)
  {
    add_named_int(cond, "NUMDOF");
    add_named_int_vector(cond, "ONOFF", "", "NUMDOF");
    add_named_real_vector(cond, "VAL", "", "NUMDOF");
    add_named_int_vector(cond, "FUNCT", "", "NUMDOF", 0, false, true);

    // optional
    cond->add_component(Teuchos::RCP(new SelectionComponent("TYPE", "Live",
        Teuchos::tuple<std::string>("Live", "Dead", "PrescribedDomainLoad", "constHydro_z",
            "increaseHydro_z", "pseudo_orthopressure", "orthopressure", "LAS", "PressureGrad",
            "Torque"),
        Teuchos::tuple<std::string>("neum_live", "neum_dead", "pres_domain_load",
            "neum_consthydro_z", "neum_increhydro_z", "neum_pseudo_orthopressure",
            "neum_orthopressure", "neum_LAS", "neum_pgrad", "neum_torque"),
        true)));
    cond->add_component(Teuchos::RCP(
        new SelectionComponent("surface", "Mid", Teuchos::tuple<std::string>("Mid", "Top", "Bot"),
            Teuchos::tuple<std::string>("mid", "top", "bot"), true)));

    condlist.emplace_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // Dirichlet conditions
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT DIRICH CONDITIONS",
          "Dirichlet", "Point Dirichlet", Core::Conditions::PointDirichlet, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linedirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE DIRICH CONDITIONS",
          "Dirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF DIRICH CONDITIONS",
          "Dirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> voldirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL DIRICH CONDITIONS",
          "Dirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointaledirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Point Dirichlet", Core::Conditions::PointDirichlet, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealedirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfaledirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volaledirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  // Dirichlet conditions for transport problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointtransportdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Point Dirichlet",
          Core::Conditions::PointDirichlet, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linetransportdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Line Dirichlet",
          Core::Conditions::LineDirichlet, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surftransportdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Surface Dirichlet",
          Core::Conditions::SurfaceDirichlet, false, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> voltransportdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Volume Dirichlet",
          Core::Conditions::VolumeDirichlet, false, Core::Conditions::geometry_type_volume));

  // Dirichlet conditions for thermo problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO DIRICH CONDITIONS", "ThermoDirichlet", "Point Dirichlet",
          Core::Conditions::PointDirichlet, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linethermodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE THERMO DIRICH CONDITIONS",
          "ThermoDirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfthermodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF THERMO DIRICH CONDITIONS",
          "ThermoDirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volthermodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL THERMO DIRICH CONDITIONS",
          "ThermoDirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  // Dirichlet conditions for poroelasticity problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointporodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Point Dirichlet", Core::Conditions::PointDirichlet, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineporodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfporodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volporodirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointnurbslsdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Point Dirichlet",
          Core::Conditions::PointDirichlet, true, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linenurbslsdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Line Dirichlet",
          Core::Conditions::LineDirichlet, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfnurbslsdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Surface Dirichlet",
          Core::Conditions::SurfaceDirichlet, true, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volnurbslsdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Volume Dirichlet",
          Core::Conditions::VolumeDirichlet, true, Core::Conditions::geometry_type_volume));

  auto all_dirichlet_conditions = {pointdirichlet, linedirichlet, surfdirichlet, voldirichlet,
      pointaledirichlet, linealedirichlet, surfaledirichlet, volaledirichlet,
      pointtransportdirichlet, linetransportdirichlet, surftransportdirichlet,
      voltransportdirichlet, pointthermodirichlet, linethermodirichlet, surfthermodirichlet,
      volthermodirichlet, pointporodirichlet, lineporodirichlet, surfporodirichlet,
      volporodirichlet, pointnurbslsdirichlet, linenurbslsdirichlet, surfnurbslsdirichlet,
      volnurbslsdirichlet};

  for (const auto& cond : all_dirichlet_conditions)
  {
    add_named_int(cond, "NUMDOF");
    add_named_int_vector(cond, "ONOFF", "", "NUMDOF");
    add_named_real_vector(cond, "VAL", "", "NUMDOF");
    add_named_int_vector(cond, "FUNCT", "", "NUMDOF", 0, false, true);

    // optional
    add_named_selection_component(cond, "TAG", "", "none",
        Teuchos::tuple<std::string>("none", "monitor_reaction"),
        Teuchos::tuple<std::string>("none", "monitor_reaction"), true);

    condlist.emplace_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // Point coupling (e.g. joints - couple X out of Y nodal DoFs)

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointcoupling =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT COUPLING CONDITIONS",
          "PointCoupling", "Point Coupling", Core::Conditions::PointCoupling, false,
          Core::Conditions::geometry_type_point));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermocoupling =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO COUPLING CONDITIONS", "PointThermoCoupling", "Point Coupling",
          Core::Conditions::PointCoupling, false, Core::Conditions::geometry_type_point));

  for (const auto& cond : {pointcoupling, pointthermocoupling})
  {
    add_named_int(cond, "NUMDOF");
    add_named_int_vector(cond, "ONOFF", "", "NUMDOF");

    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // Initial fields

  // general initial field conditions
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT INITIAL FIELD CONDITIONS", "Initfield", "Point Initfield",
          Core::Conditions::PointInitfield, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE INITIAL FIELD CONDITIONS",
          "Initfield", "Line Initfield", Core::Conditions::LineInitfield, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF INITIAL FIELD CONDITIONS",
          "Initfield", "Surface Initfield", Core::Conditions::SurfaceInitfield, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL INITIAL FIELD CONDITIONS",
          "Initfield", "Volume Initfield", Core::Conditions::VolumeInitfield, false,
          Core::Conditions::geometry_type_volume));

  for (const auto& cond : {pointinitfields, lineinitfields, surfinitfields, volinitfields})
  {
    cond->add_component(Teuchos::RCP(new SelectionComponent("Field", "Undefined",
        Teuchos::tuple<std::string>("Undefined", "Velocity", "Pressure", "Temperature", "ScaTra",
            "Porosity", "PoroMultiFluid", "Artery"),
        Teuchos::tuple<std::string>("Undefined", "Velocity", "Pressure", "Temperature", "ScaTra",
            "Porosity", "PoroMultiFluid", "Artery"))));

    // give function id - always one single integer
    // (for initial vector fields, use the COMPONENT option of our functions)
    cond->add_component(Teuchos::RCP(new IntComponent("funct")));

    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // define initial field that can be set on thermo simulations that use the ScaTra
  // discretization e.g. STI, SSTI

  // initial field conditions for temperature on ScaTra discretizations
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermoinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on points",
          Core::Conditions::PointInitfield, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linethermoinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on lines",
          Core::Conditions::LineInitfield, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfthermoinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on surfaces",
          Core::Conditions::SurfaceInitfield, false, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volthermoinitfields =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on volumes",
          Core::Conditions::VolumeInitfield, false, Core::Conditions::geometry_type_volume));

  for (const auto& cond :
      {pointthermoinitfields, linethermoinitfields, surfthermoinitfields, volthermoinitfields})
  {
    cond->add_component(Teuchos::RCP(new SelectionComponent("Field", "Undefined",
        Teuchos::tuple<std::string>("Undefined", "ScaTra"),
        Teuchos::tuple<std::string>("Undefined", "ScaTra"))));

    // give function id - always one single integer
    cond->add_component(Teuchos::RCP(new IntComponent("funct")));

    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // compute domain integrals, i.e., cumulative volumes of 3D domain elements or cumulative
  // surface areas of 2D domain elements

  // definition of surface and volume conditions for domain integral computation
  Teuchos::RCP<Core::Conditions::ConditionDefinition> domainintegralsurf = Teuchos::RCP(
      new Core::Conditions::ConditionDefinition("DESIGN DOMAIN INTEGRAL SURF CONDITIONS",
          "DomainIntegral", "compute cumulative surface areas of 2D domain elements",
          Core::Conditions::DomainIntegral, true, Core::Conditions::geometry_type_surface));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> domainintegralvol = Teuchos::RCP(
      new Core::Conditions::ConditionDefinition("DESIGN DOMAIN INTEGRAL VOL CONDITIONS",
          "DomainIntegral", "compute cumulative volumes of 3D domain elements",
          Core::Conditions::DomainIntegral, true, Core::Conditions::geometry_type_volume));

  for (const auto& cond : {domainintegralsurf, domainintegralvol})
  {
    // add input file line components to condition definitions
    cond->add_component(Teuchos::RCP(new SeparatorComponent("ID")));
    cond->add_component(Teuchos::RCP(new IntComponent("ConditionID")));

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // compute boundary integrals, i.e., cumulative surface areas of 2D boundary elements

  Teuchos::RCP<Core::Conditions::ConditionDefinition> boundaryintegralsurf = Teuchos::RCP(
      new Core::Conditions::ConditionDefinition("DESIGN BOUNDARY INTEGRAL SURF CONDITIONS",
          "BoundaryIntegral", "compute cumulative surface areas of 2D boundary elements",
          Core::Conditions::BoundaryIntegral, true, Core::Conditions::geometry_type_surface));

  // add input file line components to condition definition
  boundaryintegralsurf->add_component(Teuchos::RCP(new SeparatorComponent("ID")));
  boundaryintegralsurf->add_component(Teuchos::RCP(new IntComponent("ConditionID")));

  // insert condition definition into global list of valid condition definitions
  condlist.push_back(boundaryintegralsurf);


  /*--------------------------------------------------------------------*/
  // wear in ALE description

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealewear = Teuchos::RCP(
      new Core::Conditions::ConditionDefinition("DESIGN LINE ALE WEAR CONDITIONS 2D", "AleWear",
          "Line Ale Wear", Core::Conditions::AleWear, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfalewear =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURFACE WEAR CONDITIONS 3D",
          "AleWear", "Surface Ale Wear", Core::Conditions::AleWear, true,
          Core::Conditions::geometry_type_surface));

  condlist.push_back(linealewear);
  condlist.push_back(surfalewear);

  /*--------------------------------------------------------------------*/
  // local coordinate systems

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointlocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT LOCSYS CONDITIONS",
          "Locsys", "Point local coordinate system", Core::Conditions::PointLocsys, true,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linelocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE LOCSYS CONDITIONS",
          "Locsys", "Line local coordinate system", Core::Conditions::LineLocsys, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surflocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF LOCSYS CONDITIONS",
          "Locsys", "Surface local coordinate system", Core::Conditions::SurfaceLocsys, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> vollocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL LOCSYS CONDITIONS",
          "Locsys", "Volume local coordinate system", Core::Conditions::VolumeLocsys, true,
          Core::Conditions::geometry_type_volume));

  // Ale
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointalelocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN POINT ALE LOCSYS CONDITIONS",
          "AleLocsys", "Point local coordinate system", Core::Conditions::PointLocsys, true,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealelocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE ALE LOCSYS CONDITIONS",
          "AleLocsys", "Line local coordinate system", Core::Conditions::LineLocsys, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfalelocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURF ALE LOCSYS CONDITIONS",
          "AleLocsys", "Surface local coordinate system", Core::Conditions::SurfaceLocsys, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volalelocsys =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN VOL ALE LOCSYS CONDITIONS",
          "AleLocsys", "Volume local coordinate system", Core::Conditions::VolumeLocsys, true,
          Core::Conditions::geometry_type_volume));

  // add components to condition definitions
  for (const auto& cond : {pointlocsys, linelocsys, surflocsys, vollocsys, pointalelocsys,
           linealelocsys, surfalelocsys, volalelocsys})
  {
    add_named_real_vector(cond, "ROTANGLE", "", 3);
    add_named_int_vector(cond, "FUNCT", "", 3);
    add_named_int(cond, "USEUPDATEDNODEPOS");
  }

  // add node normal system option only for lines and surfaces
  for (const auto& cond : {linelocsys, surflocsys, linealelocsys, surfalelocsys})
  {
    add_named_int(cond, "USECONSISTENTNODENORMAL");
  }

  // add conditions to global list of conditions
  for (const auto& cond : {pointlocsys, linelocsys, surflocsys, vollocsys, pointalelocsys,
           linealelocsys, surfalelocsys, volalelocsys})
  {
    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // periodic boundary

  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineperiodic =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE PERIODIC BOUNDARY CONDITIONS", "LinePeriodic", "Line Periodic",
          Core::Conditions::LinePeriodic, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfperiodic =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF PERIODIC BOUNDARY CONDITIONS", "SurfacePeriodic", "Surface Periodic",
          Core::Conditions::SurfacePeriodic, false, Core::Conditions::geometry_type_surface));

  for (const auto& cond : {lineperiodic, surfperiodic})
  {
    cond->add_component(Teuchos::RCP(
        new IntComponent("Id of periodic boundary condition", {0, true, false, false})));
    cond->add_component(Teuchos::RCP(new SelectionComponent("Is slave periodic boundary condition",
        "Master", Teuchos::tuple<std::string>("Master", "Slave"),
        Teuchos::tuple<std::string>("Master", "Slave"))));
    add_named_selection_component(cond, "PLANE", "degrees of freedom for the pbc plane", "xy",
        Teuchos::tuple<std::string>("xy", "yx", "yz", "zy", "xz", "zx", "xyz"),
        Teuchos::tuple<std::string>("xy", "xy", "yz", "yz", "xz", "xz", "xyz"));
    add_named_int(cond, "LAYER", "layer of periodic boundary condition", 0, false, false, true);
    add_named_real(cond, "ANGLE", "angle of rotation");
    add_named_real(cond, "ABSTREETOL", "tolerance for nodematching in octree");

    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // weak Dirichlet conditions
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineweakdirichlet =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE WEAK DIRICHLET CONDITIONS", "LineWeakDirichlet", "LineWeakDirichlet",
          Core::Conditions::LineWeakDirichlet, true, Core::Conditions::geometry_type_line));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfweakdirichlet = Teuchos::RCP(
      new Core::Conditions::ConditionDefinition("DESIGN SURFACE WEAK DIRICHLET CONDITIONS",
          "SurfaceWeakDirichlet", "SurfaceWeakDirichlet", Core::Conditions::SurfaceWeakDirichlet,
          true, Core::Conditions::geometry_type_surface));

  // attach all components to those condition
  for (const auto& cond : {lineweakdirichlet, surfweakdirichlet})
  {
    // weak DBCs can be imposed adjoint consistent or adjoint inconsistent
    cond->add_component(
        Teuchos::RCP(new SelectionComponent("Choice of gamma parameter", "adjoint-consistent",
            Teuchos::tuple<std::string>("adjoint-consistent", "diffusive-optimal"),
            Teuchos::tuple<std::string>("adjoint-consistent", "diffusive-optimal"))));

    // weak DBCs can be imposed in all directions or only in normal direction
    // (SCATRA: not checked, only in all_directions so far)
    cond->add_component(Teuchos::RCP(new SelectionComponent("Directions to apply weak dbc",
        "all_directions", Teuchos::tuple<std::string>("all_directions", "only_in_normal_direction"),
        Teuchos::tuple<std::string>("all_directions", "only_in_normal_direction"))));

    // FLUID: penalty parameter either computed dynamically (using Spaldings law of
    // the wall) or by a fixed value; SCATRA: not checked, only constant value so far
    cond->add_component(Teuchos::RCP(new SelectionComponent("Definition of penalty parameter",
        "constant", Teuchos::tuple<std::string>("constant", "Spalding"),
        Teuchos::tuple<std::string>("constant", "Spalding"))));

    // scaling factor for penalty parameter tauB or
    // stabilization parameter alpha for Nitsche term
    // (SCATRA: if stabilization parameter negative -> mixed-hybrid formulation)
    cond->add_component(Teuchos::RCP(new RealComponent("TauBscaling")));

    // linearisation strategies --- the linearisation (i.e. the matrix
    // contribution) of the convective term on the inflow could be
    // suppressed, since the flux is a kink function and including this one
    // might result in even worse convergence behaviour
    // (SCATRA: not checked)
    cond->add_component(Teuchos::RCP(new SelectionComponent("Linearisation", "lin_all",
        Teuchos::tuple<std::string>("lin_all", "no_lin_conv_inflow"),
        Teuchos::tuple<std::string>("lin_all", "no_lin_conv_inflow"))));

    // we provide a vector of 3 values for velocities
    cond->add_component(Teuchos::RCP(new RealVectorComponent("VAL", 3)));

    // and optional spatial functions
    cond->add_component(Teuchos::RCP(new IntVectorComponent("FUNCT", 3, {0, false, false, true})));

    // append the condition to the list of all conditions
    condlist.push_back(cond);
  }

  /*--------------------------------------------------------------------*/
  // boundary for superconvergent patch recovery (SPR)

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linespr =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN PATCH RECOVERY BOUNDARY LINE CONDITIONS", "SPRboundary", "Boundary for SPR",
          Core::Conditions::SPRboundary, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfspr =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN PATCH RECOVERY BOUNDARY SURF CONDITIONS", "SPRboundary", "Boundary for SPR",
          Core::Conditions::SPRboundary, false, Core::Conditions::geometry_type_surface));

  condlist.push_back(linespr);
  condlist.push_back(surfspr);

  /*--------------------------------------------------------------------*/
  // volume constraint

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volumeconstraint =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURFACE VOLUME CONSTRAINT 3D",
          "VolumeConstraint_3D", "Surface Volume Constraint", Core::Conditions::VolumeConstraint_3D,
          true, Core::Conditions::geometry_type_surface));

  volumeconstraint->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  volumeconstraint->add_component(Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  volumeconstraint->add_component(Teuchos::RCP(new RealComponent("activTime")));
  volumeconstraint->add_component(Teuchos::RCP(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(volumeconstraint);

  /*--------------------------------------------------------------------*/
  // volume constraint penalty

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volumeconstraintpen =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE VOLUME CONSTRAINT 3D PEN", "VolumeConstraint_3D_Pen",
          "Surface Volume Constraint Penalty", Core::Conditions::VolumeConstraint_3D_pen, true,
          Core::Conditions::geometry_type_surface));

  volumeconstraintpen->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  volumeconstraintpen->add_component(
      Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  volumeconstraintpen->add_component(Teuchos::RCP(new RealComponent("activTime")));
  volumeconstraintpen->add_component(Teuchos::RCP(new RealComponent("penalty")));
  volumeconstraintpen->add_component(Teuchos::RCP(new RealComponent("rho")));
  volumeconstraintpen->add_component(Teuchos::RCP(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(volumeconstraintpen);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areaconstraint =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D",
          "AreaConstraint_3D", "Surface Area Constraint", Core::Conditions::AreaConstraint_3D, true,
          Core::Conditions::geometry_type_surface));

  areaconstraint->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  areaconstraint->add_component(Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  areaconstraint->add_component(Teuchos::RCP(new RealComponent("activTime")));

  condlist.push_back(areaconstraint);

  /*--------------------------------------------------------------------*/
  // area constraint penalty

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areaconstraintpen = Teuchos::RCP(
      new Core::Conditions::ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D PEN",
          "AreaConstraint_3D_Pen", "Surface Area Constraint Penalty",
          Core::Conditions::AreaConstraint_3D_pen, true, Core::Conditions::geometry_type_surface));

  areaconstraintpen->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  areaconstraintpen->add_component(Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  areaconstraintpen->add_component(Teuchos::RCP(new RealComponent("activTime")));
  areaconstraintpen->add_component(Teuchos::RCP(new RealComponent("penalty")));
  areaconstraintpen->add_component(Teuchos::RCP(new RealComponent("rho")));
  areaconstraintpen->add_component(Teuchos::RCP(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(areaconstraintpen);


  /*--------------------------------------------------------------------*/
  // volume monitor

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volumemonitor =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURFACE VOLUME MONITOR 3D",
          "VolumeMonitor_3D", "Surface Volume Monitor", Core::Conditions::VolumeMonitor_3D, true,
          Core::Conditions::geometry_type_surface));

  volumemonitor->add_component(Teuchos::RCP(new IntComponent("ConditionID")));

  condlist.push_back(volumemonitor);

  /*--------------------------------------------------------------------*/
  // area monitor 3D

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areamonitor =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN SURFACE AREA MONITOR 3D",
          "AreaMonitor_3D", "Surface Area Monitor", Core::Conditions::AreaMonitor_3D, true,
          Core::Conditions::geometry_type_surface));

  areamonitor->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  areamonitor->add_component(Teuchos::RCP(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(areamonitor);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areaconstraint2D =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE AREA CONSTRAINT 2D",
          "AreaConstraint_2D", "Line Area Constraint", Core::Conditions::AreaConstraint_2D, true,
          Core::Conditions::geometry_type_line));

  areaconstraint2D->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  areaconstraint2D->add_component(Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  areaconstraint2D->add_component(Teuchos::RCP(new RealComponent("activTime")));
  condlist.push_back(areaconstraint2D);

  /*--------------------------------------------------------------------*/
  // area monitor 2D

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areamonitor2D =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE AREA MONITOR 2D",
          "AreaMonitor_2D", "Line Area Monitor", Core::Conditions::AreaMonitor_2D, true,
          Core::Conditions::geometry_type_line));

  areamonitor2D->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  condlist.push_back(areamonitor2D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D for a node over a plane

  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodeonplaneconst3D =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE MULTIPNT CONSTRAINT 3D", "MPC_NodeOnPlane_3D", "Node on Plane Constraint",
          Core::Conditions::MPC_NodeOnPlane_3D, false, Core::Conditions::geometry_type_surface));

  nodeonplaneconst3D->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  nodeonplaneconst3D->add_component(Teuchos::RCP(new RealComponent("amplitude")));
  nodeonplaneconst3D->add_component(
      Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  nodeonplaneconst3D->add_component(Teuchos::RCP(new RealComponent("activTime")));
  nodeonplaneconst3D->add_component(Teuchos::RCP(new IntVectorComponent("planeNodes", 3)));
  nodeonplaneconst3D->add_component(Teuchos::RCP(new SelectionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodeonplaneconst3D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously

  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodemasterconst3D =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D", "MPC_NormalComponent_3D",
          "Node on Plane Constraint", Core::Conditions::MPC_NormalComponent_3D, false,
          Core::Conditions::geometry_type_surface));

  nodemasterconst3D->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  nodemasterconst3D->add_component(Teuchos::RCP(new RealComponent("amplitude")));
  nodemasterconst3D->add_component(Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  nodemasterconst3D->add_component(Teuchos::RCP(new RealComponent("activTime")));
  nodemasterconst3D->add_component(Teuchos::RCP(new IntComponent("masterNode")));
  nodemasterconst3D->add_component(Teuchos::RCP(new RealVectorComponent("direction", 3)));
  nodemasterconst3D->add_component(Teuchos::RCP(new SelectionComponent("value", "disp",
      Teuchos::tuple<std::string>("disp", "x"), Teuchos::tuple<std::string>("disp", "x"), true)));
  nodemasterconst3D->add_component(Teuchos::RCP(new SelectionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodemasterconst3D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously, penalty based

  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodemasterconst3Dpen =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D PEN", "MPC_NormalComponent_3D_Pen",
          "Node on Plane Constraint Penalty", Core::Conditions::MPC_NormalComponent_3D_pen, false,
          Core::Conditions::geometry_type_surface));

  nodemasterconst3Dpen->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new RealComponent("amplitude")));
  nodemasterconst3Dpen->add_component(
      Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new RealComponent("activTime")));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new RealComponent("penalty")));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new IntComponent("masterNode")));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new RealVectorComponent("direction", 3)));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new SelectionComponent("value", "disp",
      Teuchos::tuple<std::string>("disp", "x"), Teuchos::tuple<std::string>("disp", "x"), true)));
  nodemasterconst3Dpen->add_component(Teuchos::RCP(new SelectionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodemasterconst3Dpen);
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 2D for a node on a line
  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodeonlineconst2D =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition("DESIGN LINE MULTIPNT CONSTRAINT 2D",
          "MPC_NodeOnLine_2D", "Node on Line Constraint", Core::Conditions::MPC_NodeOnLine_2D,
          false, Core::Conditions::geometry_type_line));

  nodeonlineconst2D->add_component(Teuchos::RCP(new IntComponent("ConditionID")));
  nodeonlineconst2D->add_component(Teuchos::RCP(new RealComponent("amplitude")));
  nodeonlineconst2D->add_component(Teuchos::RCP(new IntComponent("curve", {0, true, true, false})));
  nodeonlineconst2D->add_component(Teuchos::RCP(new IntComponent("constrNode 1")));
  nodeonlineconst2D->add_component(Teuchos::RCP(new IntComponent("constrNode 2")));
  nodeonlineconst2D->add_component(Teuchos::RCP(new IntComponent("constrNode 3")));
  nodeonlineconst2D->add_component(Teuchos::RCP(
      new SelectionComponent("control value", "dist", Teuchos::tuple<std::string>("dist", "angle"),
          Teuchos::tuple<std::string>("dist", "angle"), true)));
  nodeonlineconst2D->add_component(Teuchos::RCP(new RealComponent("activTime")));
  condlist.push_back(nodeonlineconst2D);

  /*--------------------------------------------------------------------*/
  // Krylov Space Projection:
  // ========================
  // specify an unsupported (i.e. rigid body or zero energy) mode orthogonal
  // that is not excited by the body force
  //
  // examples are:
  // fluid:  pure Dirichlet (i.e. velocity) boundary conditions - pressure
  //         level is undetermined
  // scatra: pure Neumann boundary condition - level of scalar quantity is
  //         undetermined
  // solid:  insufficient support of translational or rotational rigid body
  //         modes
  //
  // for fluid and scatra, NUMDOF needs to be the number of dofs per node. the
  // following ONOFF values trigger the fixation of the quantity level (for
  // fluid, only pressure is allowed).
  // for solid NUMDOF is the number of potential rigid body modes (e.g. 6 for
  // 3-D solid, 3 for 2-D solid), where ONOFF triggers first the
  // translational followed by the rotational modes, each in/around x to z

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfrigidbodymode =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF MODE FOR KRYLOV SPACE PROJECTION", "KrylovSpaceProjection",
          "Surface mode for Krylov space projection", Core::Conditions::SurfaceModeKrylovProjection,
          true, Core::Conditions::geometry_type_surface));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volrigidbodymode =
      Teuchos::RCP(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL MODE FOR KRYLOV SPACE PROJECTION", "KrylovSpaceProjection",
          "Volume mode for Krylov space projection", Core::Conditions::VolumeModeKrylovProjection,
          true, Core::Conditions::geometry_type_volume));

  for (const auto& cond : {surfrigidbodymode, volrigidbodymode})
  {
    cond->add_component(Teuchos::RCP(new SelectionComponent("discretization", "fluid",
        Teuchos::tuple<std::string>("fluid", "scatra", "solid"),
        Teuchos::tuple<std::string>("fluid", "scatra", "solid"))));
    add_named_int(cond, "NUMMODES");
    add_named_int_vector(cond, "ONOFF", "", "NUMMODES");
    cond->add_component(Teuchos::RCP(new SelectionComponent("weight vector definition",
        "integration", Teuchos::tuple<std::string>("integration", "pointvalues"),
        Teuchos::tuple<std::string>("integration", "pointvalues"))));

    condlist.push_back(cond);
  }


  // Finally, add the problem-specific conditions from the various modules
  Inpar::Mortar::set_valid_conditions(condlist);
  Inpar::S2I::set_valid_conditions(condlist);
  Inpar::ScaTra::set_valid_conditions(condlist);
  Inpar::STI::set_valid_conditions(condlist);
  Inpar::ElCh::set_valid_conditions(condlist);
  Inpar::ElectroPhysiology::set_valid_conditions(condlist);
  Inpar::FLUID::set_valid_conditions(condlist);
  Inpar::ALE::set_valid_conditions(condlist);
  Inpar::FSI::set_valid_conditions(condlist);
  Inpar::FPSI::set_valid_conditions(condlist);
  Inpar::Immersed::set_valid_conditions(condlist);
  Inpar::XFEM::set_valid_conditions(condlist);
  Inpar::BioFilm::set_valid_conditions(condlist);
  Inpar::ArteryNetwork::set_valid_conditions(condlist);
  Inpar::ReducedLung::set_valid_conditions(condlist);
  Inpar::Cardiovascular0D::set_valid_conditions(condlist);
  Inpar::Solid::set_valid_conditions(condlist);
  Inpar::Thermo::set_valid_conditions(condlist);
  Inpar::SSI::set_valid_conditions(condlist);
  Inpar::SSTI::set_valid_conditions(condlist);
  Inpar::PARTICLE::set_valid_conditions(condlist);
  Inpar::LevelSet::set_valid_conditions(condlist);
  Inpar::EleMag::set_valid_conditions(condlist);
  Inpar::BEAMPOTENTIAL::set_valid_conditions(condlist);
  Inpar::RveMpc::set_valid_conditions(condlist);
  Inpar::BEAMINTERACTION::set_valid_conditions(condlist);
  Inpar::EHL::set_valid_conditions(condlist);
  Inpar::PoroMultiPhaseScaTra::set_valid_conditions(condlist);

  // finally some conditions that do not have their own files yet are problem-specific
  set_miscellaneous_conditions(condlist);

  return vc;
}

FOUR_C_NAMESPACE_CLOSE
