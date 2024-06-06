/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid conditions for input

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_inpar_validconditions.hpp"

#include "4C_discretization_condition_definition.hpp"
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

FOUR_C_NAMESPACE_OPEN


void Input::PrintEmptyConditionDefinitions(std::ostream& stream,
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  for (unsigned i = 0; i < condlist.size(); ++i)
  {
    condlist[i]->Print(stream, nullptr);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintConditionDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> condlist =
      Input::ValidConditions();
  Input::PrintEmptyConditionDefinitions(std::cout, *condlist);
}


namespace Input
{
  // collect some problem-specific conditions that do not fit in the generic sections
  void SetMiscellaneousConditions(
      std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
  {
    /*--------------------------------------------------------------------*/
    // microscale boundary

    Teuchos::RCP<Core::Conditions::ConditionDefinition> microscale =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition("MICROSCALE CONDITIONS",
            "MicroBoundary", "Microscale Boundary", Core::Conditions::MicroBoundary, true,
            Core::Conditions::geometry_type_surface));

    condlist.push_back(microscale);

    /*--------------------------------------------------------------------*/
    // stc layer condition

    Teuchos::RCP<Core::Conditions::ConditionDefinition> stclayer =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL STC LAYER", "STC Layer",
            "Layer for Multilayered STC", Core::Conditions::VolSTCLayer, true,
            Core::Conditions::geometry_type_volume));

    stclayer->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));

    condlist.push_back(stclayer);
  }
}  // namespace Input


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>>
Input::ValidConditions()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> vc =
      Teuchos::rcp(new std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>());

  std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist = *vc;

  /*--------------------------------------------------------------------*/
  // Neumann
  std::vector<Teuchos::RCP<LineComponent>> neumanncomponents;

  neumanncomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("NUMDOF")));
  neumanncomponents.emplace_back(Teuchos::rcp(new IntComponent("numdof")));

  neumanncomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("ONOFF")));
  neumanncomponents.emplace_back(
      Teuchos::rcp(new IntVectorComponent("onoff", LengthFromInt("numdof"))));
  neumanncomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("VAL")));
  neumanncomponents.emplace_back(
      Teuchos::rcp(new RealVectorComponent("val", LengthFromInt("numdof"))));
  neumanncomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("FUNCT")));
  neumanncomponents.emplace_back(Teuchos::rcp(
      new IntVectorComponent("funct", LengthFromInt("numdof"), {0, false, true, false})));

  // optional
  neumanncomponents.emplace_back(Teuchos::rcp(new SelectionComponent("type", "Live",
      Teuchos::tuple<std::string>("Live", "Dead", "PrescribedDomainLoad", "constHydro_z",
          "increaseHydro_z", "pseudo_orthopressure", "orthopressure", "LAS", "PressureGrad",
          "Torque"),
      Teuchos::tuple<std::string>("neum_live", "neum_dead", "pres_domain_load", "neum_consthydro_z",
          "neum_increhydro_z", "neum_pseudo_orthopressure", "neum_orthopressure", "neum_LAS",
          "neum_pgrad", "neum_torque"),
      true)));
  neumanncomponents.emplace_back(Teuchos::rcp(
      new SelectionComponent("surface", "Mid", Teuchos::tuple<std::string>("Mid", "Top", "Bot"),
          Teuchos::tuple<std::string>("mid", "top", "bot"), true)));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT NEUMANN CONDITIONS",
          "PointNeumann", "Point Neumann", Core::Conditions::PointNeumann, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointneumanneb =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT MOMENT EB CONDITIONS",
          "PointNeumannEB", "Point Neumann Moment for an Euler-Bernoulli beam",
          Core::Conditions::PointNeumannEB, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE NEUMANN CONDITIONS",
          "LineNeumann", "Line Neumann", Core::Conditions::LineNeumann, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF NEUMANN CONDITIONS",
          "SurfaceNeumann", "Surface Neumann", Core::Conditions::SurfaceNeumann, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL NEUMANN CONDITIONS",
          "VolumeNeumann", "Volume Neumann", Core::Conditions::VolumeNeumann, true,
          Core::Conditions::geometry_type_volume));

  // Neumann conditions for transport problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointtransportneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT TRANSPORT NEUMANN CONDITIONS", "TransportPointNeumann", "Point Neumann",
          Core::Conditions::PointNeumann, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linetransportneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE TRANSPORT NEUMANN CONDITIONS", "TransportLineNeumann", "Line Neumann",
          Core::Conditions::LineNeumann, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surftransportneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF TRANSPORT NEUMANN CONDITIONS", "TransportSurfaceNeumann", "Surface Neumann",
          Core::Conditions::SurfaceNeumann, true, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> voltransportneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL TRANSPORT NEUMANN CONDITIONS", "TransportVolumeNeumann", "Volume Neumann",
          Core::Conditions::VolumeNeumann, true, Core::Conditions::geometry_type_volume));

  // Neumann conditions for thermo problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO NEUMANN CONDITIONS", "ThermoPointNeumann", "Point Neumann",
          Core::Conditions::PointNeumann, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linethermoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE THERMO NEUMANN CONDITIONS", "ThermoLineNeumann", "Line Neumann",
          Core::Conditions::LineNeumann, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfthermoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF THERMO NEUMANN CONDITIONS", "ThermoSurfaceNeumann", "Surface Neumann",
          Core::Conditions::SurfaceNeumann, true, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volthermoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL THERMO NEUMANN CONDITIONS",
          "ThermoVolumeNeumann", "Volume Neumann", Core::Conditions::VolumeNeumann, true,
          Core::Conditions::geometry_type_volume));

  // Neumann conditions for poroelasticity problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointporoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT PORO NEUMANN CONDITIONS",
          "PoroPointNeumann", "Point Neumann", Core::Conditions::PointNeumann, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineporoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE PORO NEUMANN CONDITIONS",
          "PoroLineNeumann", "Line Neumann", Core::Conditions::LineNeumann, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfporoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF PORO NEUMANN CONDITIONS",
          "PoroSurfaceNeumann", "Surface Neumann", Core::Conditions::SurfaceNeumann, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volporoneumann =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL PORO NEUMANN CONDITIONS",
          "PoroVolumeNeumann", "Volume Neumann", Core::Conditions::VolumeNeumann, true,
          Core::Conditions::geometry_type_volume));

  for (const auto& neumanncomponent : neumanncomponents)
  {
    pointneumann->AddComponent(neumanncomponent);
    pointneumanneb->AddComponent(neumanncomponent);
    lineneumann->AddComponent(neumanncomponent);
    surfneumann->AddComponent(neumanncomponent);
    volneumann->AddComponent(neumanncomponent);

    pointtransportneumann->AddComponent(neumanncomponent);
    linetransportneumann->AddComponent(neumanncomponent);
    surftransportneumann->AddComponent(neumanncomponent);
    voltransportneumann->AddComponent(neumanncomponent);

    pointthermoneumann->AddComponent(neumanncomponent);
    linethermoneumann->AddComponent(neumanncomponent);
    surfthermoneumann->AddComponent(neumanncomponent);
    volthermoneumann->AddComponent(neumanncomponent);

    pointporoneumann->AddComponent(neumanncomponent);
    lineporoneumann->AddComponent(neumanncomponent);
    surfporoneumann->AddComponent(neumanncomponent);
    volporoneumann->AddComponent(neumanncomponent);
  }

  condlist.push_back(pointneumann);
  condlist.push_back(pointneumanneb);
  condlist.push_back(lineneumann);
  condlist.push_back(surfneumann);
  condlist.push_back(volneumann);

  condlist.push_back(pointtransportneumann);
  condlist.push_back(linetransportneumann);
  condlist.push_back(surftransportneumann);
  condlist.push_back(voltransportneumann);

  condlist.push_back(pointthermoneumann);
  condlist.push_back(linethermoneumann);
  condlist.push_back(surfthermoneumann);
  condlist.push_back(volthermoneumann);

  condlist.push_back(pointporoneumann);
  condlist.push_back(lineporoneumann);
  condlist.push_back(surfporoneumann);
  condlist.push_back(volporoneumann);

  /*--------------------------------------------------------------------*/
  // Dirichlet

  std::vector<Teuchos::RCP<LineComponent>> dirichletbundcomponents;

  dirichletbundcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("NUMDOF")));
  dirichletbundcomponents.emplace_back(Teuchos::rcp(new IntComponent("numdof")));

  dirichletbundcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("ONOFF")));
  dirichletbundcomponents.emplace_back(
      Teuchos::rcp(new IntVectorComponent("onoff", LengthFromInt("numdof"))));
  dirichletbundcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("VAL")));
  dirichletbundcomponents.emplace_back(
      Teuchos::rcp(new RealVectorComponent("val", LengthFromInt("numdof"))));
  dirichletbundcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("FUNCT")));
  dirichletbundcomponents.emplace_back(Teuchos::rcp(
      new IntVectorComponent("funct", LengthFromInt("numdof"), {0, false, true, false})));

  // optional
  dirichletbundcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("TAG", "", true)));
  dirichletbundcomponents.emplace_back(Teuchos::rcp(
      new SelectionComponent("tag", "none", Teuchos::tuple<std::string>("none", "monitor_reaction"),
          Teuchos::tuple<std::string>("none", "monitor_reaction"), true)));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT DIRICH CONDITIONS",
          "Dirichlet", "Point Dirichlet", Core::Conditions::PointDirichlet, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linedirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE DIRICH CONDITIONS",
          "Dirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF DIRICH CONDITIONS",
          "Dirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> voldirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL DIRICH CONDITIONS",
          "Dirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointaledirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Point Dirichlet", Core::Conditions::PointDirichlet, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealedirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfaledirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volaledirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL ALE DIRICH CONDITIONS",
          "ALEDirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  // Dirichlet conditions for transport problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointtransportdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Point Dirichlet",
          Core::Conditions::PointDirichlet, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linetransportdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Line Dirichlet",
          Core::Conditions::LineDirichlet, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surftransportdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Surface Dirichlet",
          Core::Conditions::SurfaceDirichlet, false, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> voltransportdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL TRANSPORT DIRICH CONDITIONS", "TransportDirichlet", "Volume Dirichlet",
          Core::Conditions::VolumeDirichlet, false, Core::Conditions::geometry_type_volume));

  // Dirichlet conditions for thermo problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO DIRICH CONDITIONS", "ThermoDirichlet", "Point Dirichlet",
          Core::Conditions::PointDirichlet, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linethermodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE THERMO DIRICH CONDITIONS",
          "ThermoDirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfthermodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF THERMO DIRICH CONDITIONS",
          "ThermoDirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volthermodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL THERMO DIRICH CONDITIONS",
          "ThermoDirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  // Dirichlet conditions for poroelasticity problems
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointporodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Point Dirichlet", Core::Conditions::PointDirichlet, false,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineporodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Line Dirichlet", Core::Conditions::LineDirichlet, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfporodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Surface Dirichlet", Core::Conditions::SurfaceDirichlet, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volporodirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL PORO DIRICH CONDITIONS",
          "PoroDirichlet", "Volume Dirichlet", Core::Conditions::VolumeDirichlet, false,
          Core::Conditions::geometry_type_volume));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointnurbslsdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Point Dirichlet",
          Core::Conditions::PointDirichlet, true, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linenurbslsdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Line Dirichlet",
          Core::Conditions::LineDirichlet, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfnurbslsdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Surface Dirichlet",
          Core::Conditions::SurfaceDirichlet, true, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volnurbslsdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet", "Volume Dirichlet",
          Core::Conditions::VolumeDirichlet, true, Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    pointdirichlet->AddComponent(dirichletbundcomponents[i]);
    linedirichlet->AddComponent(dirichletbundcomponents[i]);
    surfdirichlet->AddComponent(dirichletbundcomponents[i]);
    voldirichlet->AddComponent(dirichletbundcomponents[i]);

    pointaledirichlet->AddComponent(dirichletbundcomponents[i]);
    linealedirichlet->AddComponent(dirichletbundcomponents[i]);
    surfaledirichlet->AddComponent(dirichletbundcomponents[i]);
    volaledirichlet->AddComponent(dirichletbundcomponents[i]);

    pointtransportdirichlet->AddComponent(dirichletbundcomponents[i]);
    linetransportdirichlet->AddComponent(dirichletbundcomponents[i]);
    surftransportdirichlet->AddComponent(dirichletbundcomponents[i]);
    voltransportdirichlet->AddComponent(dirichletbundcomponents[i]);

    pointthermodirichlet->AddComponent(dirichletbundcomponents[i]);
    linethermodirichlet->AddComponent(dirichletbundcomponents[i]);
    surfthermodirichlet->AddComponent(dirichletbundcomponents[i]);
    volthermodirichlet->AddComponent(dirichletbundcomponents[i]);

    pointporodirichlet->AddComponent(dirichletbundcomponents[i]);
    lineporodirichlet->AddComponent(dirichletbundcomponents[i]);
    surfporodirichlet->AddComponent(dirichletbundcomponents[i]);
    volporodirichlet->AddComponent(dirichletbundcomponents[i]);

    pointnurbslsdirichlet->AddComponent(dirichletbundcomponents[i]);
    linenurbslsdirichlet->AddComponent(dirichletbundcomponents[i]);
    surfnurbslsdirichlet->AddComponent(dirichletbundcomponents[i]);
    volnurbslsdirichlet->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(pointdirichlet);
  condlist.push_back(linedirichlet);
  condlist.push_back(surfdirichlet);
  condlist.push_back(voldirichlet);

  condlist.push_back(pointaledirichlet);
  condlist.push_back(linealedirichlet);
  condlist.push_back(surfaledirichlet);
  condlist.push_back(volaledirichlet);

  condlist.push_back(pointtransportdirichlet);
  condlist.push_back(linetransportdirichlet);
  condlist.push_back(surftransportdirichlet);
  condlist.push_back(voltransportdirichlet);

  condlist.push_back(pointthermodirichlet);
  condlist.push_back(linethermodirichlet);
  condlist.push_back(surfthermodirichlet);
  condlist.push_back(volthermodirichlet);

  condlist.push_back(pointporodirichlet);
  condlist.push_back(lineporodirichlet);
  condlist.push_back(surfporodirichlet);
  condlist.push_back(volporodirichlet);

  condlist.push_back(pointnurbslsdirichlet);
  condlist.push_back(linenurbslsdirichlet);
  condlist.push_back(surfnurbslsdirichlet);
  condlist.push_back(volnurbslsdirichlet);

  /*--------------------------------------------------------------------*/
  // Point coupling (e.g. joints - couple X out of Y nodal DoFs)

  std::vector<Teuchos::RCP<LineComponent>> couplingcomponents;

  couplingcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("NUMDOF")));
  couplingcomponents.emplace_back(Teuchos::rcp(new IntComponent("numdof")));
  couplingcomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("ONOFF")));
  couplingcomponents.emplace_back(
      Teuchos::rcp(new IntVectorComponent("onoff", LengthFromInt("numdof"))));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointcoupling =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT COUPLING CONDITIONS",
          "PointCoupling", "Point Coupling", Core::Conditions::PointCoupling, false,
          Core::Conditions::geometry_type_point));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermocoupling =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO COUPLING CONDITIONS", "PointThermoCoupling", "Point Coupling",
          Core::Conditions::PointCoupling, false, Core::Conditions::geometry_type_point));

  for (const auto& couplingcomponent : couplingcomponents)
  {
    pointcoupling->AddComponent(couplingcomponent);
    pointthermocoupling->AddComponent(couplingcomponent);
  }

  condlist.push_back(pointcoupling);
  condlist.push_back(pointthermocoupling);

  /*--------------------------------------------------------------------*/
  // Initial fields

  // define initial fields that can be set
  std::vector<Teuchos::RCP<LineComponent>> initial_field_components;
  initial_field_components.emplace_back(Teuchos::rcp(new SelectionComponent("Field", "Undefined",
      Teuchos::tuple<std::string>("Undefined", "Velocity", "Pressure", "Temperature", "ScaTra",
          "Porosity", "PoroMultiFluid", "Artery"),
      Teuchos::tuple<std::string>("Undefined", "Velocity", "Pressure", "Temperature", "ScaTra",
          "Porosity", "PoroMultiFluid", "Artery"))));

  // give function id - always one single integer
  // (for initial vector fields, use the COMPONENT option of our functions)
  initial_field_components.emplace_back(Teuchos::rcp(new IntComponent("funct")));

  // general initial field conditions
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT INITIAL FIELD CONDITIONS", "Initfield", "Point Initfield",
          Core::Conditions::PointInitfield, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE INITIAL FIELD CONDITIONS",
          "Initfield", "Line Initfield", Core::Conditions::LineInitfield, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF INITIAL FIELD CONDITIONS",
          "Initfield", "Surface Initfield", Core::Conditions::SurfaceInitfield, false,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL INITIAL FIELD CONDITIONS",
          "Initfield", "Volume Initfield", Core::Conditions::VolumeInitfield, false,
          Core::Conditions::geometry_type_volume));

  for (const auto& initial_field_component : initial_field_components)
  {
    pointinitfields->AddComponent(initial_field_component);
    lineinitfields->AddComponent(initial_field_component);
    surfinitfields->AddComponent(initial_field_component);
    volinitfields->AddComponent(initial_field_component);
  }

  condlist.push_back(pointinitfields);
  condlist.push_back(lineinitfields);
  condlist.push_back(surfinitfields);
  condlist.push_back(volinitfields);

  /*--------------------------------------------------------------------*/
  // define initial field that can be set on thermo simulations that use the ScaTra
  // discretization e.g. STI, SSTI
  std::vector<Teuchos::RCP<LineComponent>> initial_field_components_thermo_on_scatra_dis;
  initial_field_components_thermo_on_scatra_dis.emplace_back(Teuchos::rcp(new SelectionComponent(
      "Field", "Undefined", Teuchos::tuple<std::string>("Undefined", "ScaTra"),
      Teuchos::tuple<std::string>("Undefined", "ScaTra"))));

  // give function id - always one single integer
  initial_field_components_thermo_on_scatra_dis.emplace_back(
      Teuchos::rcp(new IntComponent("funct")));

  // initial field conditions for temperature on ScaTra discretizations
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointthermoinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on points",
          Core::Conditions::PointInitfield, false, Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linethermoinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on lines",
          Core::Conditions::LineInitfield, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfthermoinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on surfaces",
          Core::Conditions::SurfaceInitfield, false, Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volthermoinitfields =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN VOL THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Set the initial temperature field if the thermo field is solved using a ScaTra "
          "discretization (e.g. STI, SSTI) on volumes",
          Core::Conditions::VolumeInitfield, false, Core::Conditions::geometry_type_volume));

  for (const auto& component : initial_field_components_thermo_on_scatra_dis)
  {
    pointthermoinitfields->AddComponent(component);
    linethermoinitfields->AddComponent(component);
    surfthermoinitfields->AddComponent(component);
    volthermoinitfields->AddComponent(component);
  }

  condlist.push_back(pointthermoinitfields);
  condlist.push_back(linethermoinitfields);
  condlist.push_back(surfthermoinitfields);
  condlist.push_back(volthermoinitfields);

  /*--------------------------------------------------------------------*/
  // compute domain integrals, i.e., cumulative volumes of 3D domain elements or cumulative
  // surface areas of 2D domain elements
  {
    // definition of surface and volume conditions for domain integral computation
    Teuchos::RCP<Core::Conditions::ConditionDefinition> domainintegralsurf = Teuchos::rcp(
        new Core::Conditions::ConditionDefinition("DESIGN DOMAIN INTEGRAL SURF CONDITIONS",
            "DomainIntegral", "compute cumulative surface areas of 2D domain elements",
            Core::Conditions::DomainIntegral, true, Core::Conditions::geometry_type_surface));

    Teuchos::RCP<Core::Conditions::ConditionDefinition> domainintegralvol = Teuchos::rcp(
        new Core::Conditions::ConditionDefinition("DESIGN DOMAIN INTEGRAL VOL CONDITIONS",
            "DomainIntegral", "compute cumulative volumes of 3D domain elements",
            Core::Conditions::DomainIntegral, true, Core::Conditions::geometry_type_volume));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<LineComponent>> domainintegralcomponents;

    {
      domainintegralcomponents.push_back(Teuchos::rcp(new SeparatorComponent("ID")));
      domainintegralcomponents.push_back(Teuchos::rcp(new IntComponent("ConditionID")));
    }

    // insert input file line components into condition definitions
    for (unsigned i = 0; i < domainintegralcomponents.size(); ++i)
    {
      domainintegralsurf->AddComponent(domainintegralcomponents[i]);
      domainintegralvol->AddComponent(domainintegralcomponents[i]);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(domainintegralsurf);
    condlist.push_back(domainintegralvol);
  }

  /*--------------------------------------------------------------------*/
  // compute boundary integrals, i.e., cumulative surface areas of 2D boundary elements
  {
    Teuchos::RCP<Core::Conditions::ConditionDefinition> boundaryintegralsurf = Teuchos::rcp(
        new Core::Conditions::ConditionDefinition("DESIGN BOUNDARY INTEGRAL SURF CONDITIONS",
            "BoundaryIntegral", "compute cumulative surface areas of 2D boundary elements",
            Core::Conditions::BoundaryIntegral, true, Core::Conditions::geometry_type_surface));

    // equip condition definition with input file line components
    std::vector<Teuchos::RCP<LineComponent>> boundaryintegralcomponents;

    {
      boundaryintegralcomponents.push_back(Teuchos::rcp(new SeparatorComponent("ID")));
      boundaryintegralcomponents.push_back(Teuchos::rcp(new IntComponent("ConditionID")));
    }

    // insert input file line components into condition definition
    for (unsigned i = 0; i < boundaryintegralcomponents.size(); ++i)
      boundaryintegralsurf->AddComponent(boundaryintegralcomponents[i]);

    // insert condition definition into global list of valid condition definitions
    condlist.push_back(boundaryintegralsurf);
  }

  /*--------------------------------------------------------------------*/
  // wear in ALE description

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealewear = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE ALE WEAR CONDITIONS 2D", "AleWear",
          "Line Ale Wear", Core::Conditions::AleWear, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfalewear =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURFACE WEAR CONDITIONS 3D",
          "AleWear", "Surface Ale Wear", Core::Conditions::AleWear, true,
          Core::Conditions::geometry_type_surface));

  condlist.push_back(linealewear);
  condlist.push_back(surfalewear);

  /*--------------------------------------------------------------------*/
  // local coordinate systems

  std::vector<Teuchos::RCP<LineComponent>> locsyscomponents;

  locsyscomponents.push_back(Teuchos::rcp(new SeparatorComponent("ROTANGLE")));
  locsyscomponents.push_back(Teuchos::rcp(new RealVectorComponent("rotangle", 3)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorComponent("FUNCT")));
  locsyscomponents.push_back(Teuchos::rcp(new IntVectorComponent("funct", 3)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorComponent("USEUPDATEDNODEPOS")));
  locsyscomponents.push_back(Teuchos::rcp(new IntComponent("useupdatednodepos")));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorComponent("USECONSISTENTNODENORMAL")));
  locsyscomponents.push_back(Teuchos::rcp(new IntComponent("useconsistentnodenormal")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointlocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT LOCSYS CONDITIONS",
          "Locsys", "Point local coordinate system", Core::Conditions::PointLocsys, true,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linelocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE LOCSYS CONDITIONS",
          "Locsys", "Line local coordinate system", Core::Conditions::LineLocsys, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surflocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF LOCSYS CONDITIONS",
          "Locsys", "Surface local coordinate system", Core::Conditions::SurfaceLocsys, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> vollocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL LOCSYS CONDITIONS",
          "Locsys", "Volume local coordinate system", Core::Conditions::VolumeLocsys, true,
          Core::Conditions::geometry_type_volume));

  // Ale
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointalelocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN POINT ALE LOCSYS CONDITIONS",
          "AleLocsys", "Point local coordinate system", Core::Conditions::PointLocsys, true,
          Core::Conditions::geometry_type_point));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealelocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE ALE LOCSYS CONDITIONS",
          "AleLocsys", "Line local coordinate system", Core::Conditions::LineLocsys, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfalelocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF ALE LOCSYS CONDITIONS",
          "AleLocsys", "Surface local coordinate system", Core::Conditions::SurfaceLocsys, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volalelocsys =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN VOL ALE LOCSYS CONDITIONS",
          "AleLocsys", "Volume local coordinate system", Core::Conditions::VolumeLocsys, true,
          Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < locsyscomponents.size() - 2; ++i)
  {
    pointlocsys->AddComponent(locsyscomponents[i]);
    linelocsys->AddComponent(locsyscomponents[i]);
    surflocsys->AddComponent(locsyscomponents[i]);
    vollocsys->AddComponent(locsyscomponents[i]);

    // Ale
    pointalelocsys->AddComponent(locsyscomponents[i]);
    linealelocsys->AddComponent(locsyscomponents[i]);
    surfalelocsys->AddComponent(locsyscomponents[i]);
    volalelocsys->AddComponent(locsyscomponents[i]);
  }

  // Add node normal system option only for lines and surfaces
  linelocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 2]);
  linelocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 1]);
  surflocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 2]);
  surflocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 1]);
  linealelocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 2]);
  linealelocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 1]);
  surfalelocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 2]);
  surfalelocsys->AddComponent(locsyscomponents[locsyscomponents.size() - 1]);

  condlist.push_back(pointlocsys);
  condlist.push_back(linelocsys);
  condlist.push_back(surflocsys);
  condlist.push_back(vollocsys);

  // Ale
  condlist.push_back(pointalelocsys);
  condlist.push_back(linealelocsys);
  condlist.push_back(surfalelocsys);
  condlist.push_back(volalelocsys);

  /*--------------------------------------------------------------------*/
  // periodic boundary

  std::vector<Teuchos::RCP<LineComponent>> pbccomponents;

  pbccomponents.push_back(
      Teuchos::rcp(new IntComponent("Id of periodic boundary condition", {0, true})));
  pbccomponents.push_back(
      Teuchos::rcp(new SelectionComponent("Is slave periodic boundary condition", "Master",
          Teuchos::tuple<std::string>("Master", "Slave"),
          Teuchos::tuple<std::string>("Master", "Slave"))));
  pbccomponents.push_back(Teuchos::rcp(new SeparatorComponent("PLANE")));
  pbccomponents.push_back(
      Teuchos::rcp(new SelectionComponent("degrees of freedom for the pbc plane", "xy",
          Teuchos::tuple<std::string>("xy", "yx", "yz", "zy", "xz", "zx", "xyz"),
          Teuchos::tuple<std::string>("xy", "xy", "yz", "yz", "xz", "xz", "xyz"))));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorComponent("LAYER")));
  pbccomponents.push_back(
      Teuchos::rcp(new IntComponent("Layer of periodic boundary condition", {0, true})));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorComponent("ANGLE")));
  pbccomponents.push_back(Teuchos::rcp(new RealComponent("Angle of rotation")));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorComponent("ABSTREETOL")));
  pbccomponents.push_back(Teuchos::rcp(new RealComponent("Tolerance for nodematching in octree")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineperiodic =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE PERIODIC BOUNDARY CONDITIONS", "LinePeriodic", "Line Periodic",
          Core::Conditions::LinePeriodic, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfperiodic =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF PERIODIC BOUNDARY CONDITIONS", "SurfacePeriodic", "Surface Periodic",
          Core::Conditions::SurfacePeriodic, false, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < pbccomponents.size(); ++i)
  {
    lineperiodic->AddComponent(pbccomponents[i]);
    surfperiodic->AddComponent(pbccomponents[i]);
  }

  condlist.push_back(lineperiodic);
  condlist.push_back(surfperiodic);

  /*--------------------------------------------------------------------*/
  // weak Dirichlet conditions

  std::vector<Teuchos::RCP<LineComponent>> weakDirichletcomponents;

  // weak DBCs can be imposed adjoint consistent or adjoint inconsistent
  weakDirichletcomponents.push_back(Teuchos::rcp(new SelectionComponent("Choice of gamma parameter",
      "adjoint-consistent", Teuchos::tuple<std::string>("adjoint-consistent", "diffusive-optimal"),
      Teuchos::tuple<std::string>("adjoint-consistent", "diffusive-optimal"))));

  // weak DBCs can be imposed in all directions or only in normal direction
  // (SCATRA: not checked, only in all_directions so far)
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new SelectionComponent("Directions to apply weak dbc", "all_directions",
          Teuchos::tuple<std::string>("all_directions", "only_in_normal_direction"),
          Teuchos::tuple<std::string>("all_directions", "only_in_normal_direction"))));

  // FLUID: penalty parameter either computed dynamically (using Spaldings law of
  // the wall) or by a fixed value; SCATRA: not checked, only constant value so far
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new SelectionComponent("Definition of penalty parameter", "constant",
          Teuchos::tuple<std::string>("constant", "Spalding"),
          Teuchos::tuple<std::string>("constant", "Spalding"))));

  // scaling factor for penalty parameter tauB or
  // stabilization parameter alpha for Nitsche term
  // (SCATRA: if stabilization parameter negative -> mixed-hybrid formulation)
  weakDirichletcomponents.push_back(Teuchos::rcp(new RealComponent("TauBscaling")));

  // linearisation strategies --- the linearisation (i.e. the matrix
  // contribution) of the convective term on the inflow could be
  // suppressed, since the flux is a kink function and including this one
  // might result in even worse convergence behaviour
  // (SCATRA: not checked)
  weakDirichletcomponents.push_back(Teuchos::rcp(new SelectionComponent("Linearisation", "lin_all",
      Teuchos::tuple<std::string>("lin_all", "no_lin_conv_inflow"),
      Teuchos::tuple<std::string>("lin_all", "no_lin_conv_inflow"))));

  // we provide a vector of 3 values for velocities
  weakDirichletcomponents.push_back(Teuchos::rcp(new RealVectorComponent("val", 3)));

  // and optional spatial functions
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new IntVectorComponent("funct", 3, {0, false, false, true})));


  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineweakdirichlet =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE WEAK DIRICHLET CONDITIONS", "LineWeakDirichlet", "LineWeakDirichlet",
          Core::Conditions::LineWeakDirichlet, true, Core::Conditions::geometry_type_line));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfweakdirichlet = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN SURFACE WEAK DIRICHLET CONDITIONS",
          "SurfaceWeakDirichlet", "SurfaceWeakDirichlet", Core::Conditions::SurfaceWeakDirichlet,
          true, Core::Conditions::geometry_type_surface));

  // we attach all the components of this condition to this weak line DBC
  for (unsigned i = 0; i < weakDirichletcomponents.size(); ++i)
  {
    lineweakdirichlet->AddComponent(weakDirichletcomponents[i]);
    surfweakdirichlet->AddComponent(weakDirichletcomponents[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(lineweakdirichlet);
  condlist.push_back(surfweakdirichlet);

  /*--------------------------------------------------------------------*/
  // boundary for superconvergent patch recovery (SPR)

  std::vector<Teuchos::RCP<LineComponent>> sprcomponents;

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linespr =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN PATCH RECOVERY BOUNDARY LINE CONDITIONS", "SPRboundary", "Boundary for SPR",
          Core::Conditions::SPRboundary, false, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfspr =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN PATCH RECOVERY BOUNDARY SURF CONDITIONS", "SPRboundary", "Boundary for SPR",
          Core::Conditions::SPRboundary, false, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < sprcomponents.size(); ++i)
  {
    linespr->AddComponent(sprcomponents[i]);
    surfspr->AddComponent(sprcomponents[i]);
  }

  condlist.push_back(linespr);
  condlist.push_back(surfspr);

  /*--------------------------------------------------------------------*/
  // volume constraint

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volumeconstraint =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURFACE VOLUME CONSTRAINT 3D",
          "VolumeConstraint_3D", "Surface Volume Constraint", Core::Conditions::VolumeConstraint_3D,
          true, Core::Conditions::geometry_type_surface));

  volumeconstraint->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  volumeconstraint->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  volumeconstraint->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  volumeconstraint->AddComponent(Teuchos::rcp(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(volumeconstraint);

  /*--------------------------------------------------------------------*/
  // volume constraint penalty

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volumeconstraintpen =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE VOLUME CONSTRAINT 3D PEN", "VolumeConstraint_3D_Pen",
          "Surface Volume Constraint Penalty", Core::Conditions::VolumeConstraint_3D_pen, true,
          Core::Conditions::geometry_type_surface));

  volumeconstraintpen->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealComponent("penalty")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealComponent("rho")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(volumeconstraintpen);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areaconstraint =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D",
          "AreaConstraint_3D", "Surface Area Constraint", Core::Conditions::AreaConstraint_3D, true,
          Core::Conditions::geometry_type_surface));

  areaconstraint->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  areaconstraint->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  areaconstraint->AddComponent(Teuchos::rcp(new RealComponent("activTime")));

  condlist.push_back(areaconstraint);

  /*--------------------------------------------------------------------*/
  // area constraint penalty

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areaconstraintpen = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D PEN",
          "AreaConstraint_3D_Pen", "Surface Area Constraint Penalty",
          Core::Conditions::AreaConstraint_3D_pen, true, Core::Conditions::geometry_type_surface));

  areaconstraintpen->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealComponent("penalty")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealComponent("rho")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(areaconstraintpen);


  /*--------------------------------------------------------------------*/
  // volume monitor

  Teuchos::RCP<Core::Conditions::ConditionDefinition> volumemonitor =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURFACE VOLUME MONITOR 3D",
          "VolumeMonitor_3D", "Surface Volume Monitor", Core::Conditions::VolumeMonitor_3D, true,
          Core::Conditions::geometry_type_surface));

  volumemonitor->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));

  condlist.push_back(volumemonitor);

  /*--------------------------------------------------------------------*/
  // area monitor 3D

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areamonitor =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURFACE AREA MONITOR 3D",
          "AreaMonitor_3D", "Surface Area Monitor", Core::Conditions::AreaMonitor_3D, true,
          Core::Conditions::geometry_type_surface));

  areamonitor->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  areamonitor->AddComponent(Teuchos::rcp(new SelectionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(areamonitor);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areaconstraint2D =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE AREA CONSTRAINT 2D",
          "AreaConstraint_2D", "Line Area Constraint", Core::Conditions::AreaConstraint_2D, true,
          Core::Conditions::geometry_type_line));

  areaconstraint2D->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  areaconstraint2D->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  areaconstraint2D->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  condlist.push_back(areaconstraint2D);

  /*--------------------------------------------------------------------*/
  // area monitor 2D

  Teuchos::RCP<Core::Conditions::ConditionDefinition> areamonitor2D =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE AREA MONITOR 2D",
          "AreaMonitor_2D", "Line Area Monitor", Core::Conditions::AreaMonitor_2D, true,
          Core::Conditions::geometry_type_line));

  areamonitor2D->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  condlist.push_back(areamonitor2D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D for a node over a plane

  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodeonplaneconst3D =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE MULTIPNT CONSTRAINT 3D", "MPC_NodeOnPlane_3D", "Node on Plane Constraint",
          Core::Conditions::MPC_NodeOnPlane_3D, false, Core::Conditions::geometry_type_surface));

  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new RealComponent("amplitude")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntVectorComponent("planeNodes", 3)));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new SelectionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodeonplaneconst3D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously

  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodemasterconst3D =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D", "MPC_NormalComponent_3D",
          "Node on Plane Constraint", Core::Conditions::MPC_NormalComponent_3D, false,
          Core::Conditions::geometry_type_surface));

  nodemasterconst3D->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new RealComponent("amplitude")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new IntComponent("masterNode")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new RealVectorComponent("direction", 3)));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new SelectionComponent("value", "disp",
      Teuchos::tuple<std::string>("disp", "x"), Teuchos::tuple<std::string>("disp", "x"), true)));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new SelectionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodemasterconst3D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously, penalty based

  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodemasterconst3Dpen =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D PEN", "MPC_NormalComponent_3D_Pen",
          "Node on Plane Constraint Penalty", Core::Conditions::MPC_NormalComponent_3D_pen, false,
          Core::Conditions::geometry_type_surface));

  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealComponent("amplitude")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealComponent("penalty")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntComponent("masterNode")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealVectorComponent("direction", 3)));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new SelectionComponent("value", "disp",
      Teuchos::tuple<std::string>("disp", "x"), Teuchos::tuple<std::string>("disp", "x"), true)));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new SelectionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodemasterconst3Dpen);
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 2D for a node on a line
  Teuchos::RCP<Core::Conditions::ConditionDefinition> nodeonlineconst2D =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE MULTIPNT CONSTRAINT 2D",
          "MPC_NodeOnLine_2D", "Node on Line Constraint", Core::Conditions::MPC_NodeOnLine_2D,
          false, Core::Conditions::geometry_type_line));

  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntComponent("ConditionID")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new RealComponent("amplitude")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntComponent("curve", {0, true, true})));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntComponent("constrNode 1")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntComponent("constrNode 2")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntComponent("constrNode 3")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(
      new SelectionComponent("control value", "dist", Teuchos::tuple<std::string>("dist", "angle"),
          Teuchos::tuple<std::string>("dist", "angle"), true)));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new RealComponent("activTime")));
  condlist.push_back(nodeonlineconst2D);


  {
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


    std::vector<Teuchos::RCP<LineComponent>> rigidbodymodecomponents;

    rigidbodymodecomponents.emplace_back(Teuchos::rcp(new SelectionComponent("discretization",
        "fluid", Teuchos::tuple<std::string>("fluid", "scatra", "solid"),
        Teuchos::tuple<std::string>("fluid", "scatra", "solid"))));

    rigidbodymodecomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("NUMMODES")));
    rigidbodymodecomponents.emplace_back(Teuchos::rcp(new IntComponent("NUMMODES")));
    rigidbodymodecomponents.emplace_back(Teuchos::rcp(new SeparatorComponent("ONOFF")));
    rigidbodymodecomponents.emplace_back(
        Teuchos::rcp(new IntVectorComponent("ONOFF", LengthFromInt("NUMMODES"))));

    rigidbodymodecomponents.emplace_back(
        Teuchos::rcp(new SelectionComponent("weight vector definition", "integration",
            Teuchos::tuple<std::string>("integration", "pointvalues"),
            Teuchos::tuple<std::string>("integration", "pointvalues"))));

    Teuchos::RCP<Core::Conditions::ConditionDefinition> surfrigidbodymode = Teuchos::rcp(
        new Core::Conditions::ConditionDefinition("DESIGN SURF MODE FOR KRYLOV SPACE PROJECTION",
            "KrylovSpaceProjection", "Surface mode for Krylov space projection",
            Core::Conditions::SurfaceModeKrylovProjection, true,
            Core::Conditions::geometry_type_surface));

    Teuchos::RCP<Core::Conditions::ConditionDefinition> volrigidbodymode =
        Teuchos::rcp(new Core::Conditions::ConditionDefinition(
            "DESIGN VOL MODE FOR KRYLOV SPACE PROJECTION", "KrylovSpaceProjection",
            "Volume mode for Krylov space projection", Core::Conditions::VolumeModeKrylovProjection,
            true, Core::Conditions::geometry_type_volume));

    for (unsigned i = 0; i < rigidbodymodecomponents.size(); ++i)
    {
      surfrigidbodymode->AddComponent(rigidbodymodecomponents[i]);
      volrigidbodymode->AddComponent(rigidbodymodecomponents[i]);
    }
    condlist.push_back(surfrigidbodymode);
    condlist.push_back(volrigidbodymode);
  }

  // Finally, add the problem-specific conditions from the various modules

  Inpar::Mortar::SetValidConditions(condlist);

  Inpar::S2I::SetValidConditions(condlist);

  Inpar::ScaTra::SetValidConditions(condlist);

  Inpar::STI::SetValidConditions(condlist);

  Inpar::ElCh::SetValidConditions(condlist);

  Inpar::ElectroPhysiology::SetValidConditions(condlist);

  Inpar::FLUID::SetValidConditions(condlist);

  Inpar::ALE::SetValidConditions(condlist);

  Inpar::FSI::SetValidConditions(condlist);

  Inpar::FPSI::SetValidConditions(condlist);

  Inpar::Immersed::SetValidConditions(condlist);

  Inpar::XFEM::SetValidConditions(dirichletbundcomponents, neumanncomponents, condlist);

  Inpar::BioFilm::SetValidConditions(condlist);

  Inpar::ArteryNetwork::SetValidConditions(condlist);

  Inpar::ReducedLung::SetValidConditions(condlist);

  Inpar::CARDIOVASCULAR0D::SetValidConditions(condlist);

  Inpar::STR::SetValidConditions(condlist);

  Inpar::THR::SetValidConditions(condlist);

  Inpar::SSI::SetValidConditions(condlist);

  Inpar::SSTI::SetValidConditions(condlist);

  Inpar::PARTICLE::SetValidConditions(condlist);

  Inpar::LevelSet::SetValidConditions(condlist);

  Inpar::EleMag::SetValidConditions(condlist);

  Inpar::BEAMPOTENTIAL::SetValidConditions(condlist);

  Inpar::RveMpc::SetValidConditions(condlist);

  Inpar::BEAMINTERACTION::SetValidConditions(condlist);

  Inpar::EHL::SetValidConditions(condlist);

  Inpar::PoroMultiPhaseScaTra::SetValidConditions(condlist);

  // finally some conditions that do not have their own files yet are problem-specific
  SetMiscellaneousConditions(condlist);

  return vc;
}

FOUR_C_NAMESPACE_CLOSE
