/*----------------------------------------------------------------------*/
/*!
\file drt_validconditions.cpp

\brief Setup of the list of valid conditions for input

\level 1

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/


#include "drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "inpar_s2i.H"
#include "inpar_scatra.H"
#include "inpar_sti.H"
#include "inpar_fluid.H"
#include "inpar_structure.H"
#include "inpar_thermo.H"
#include "inpar_ale.H"
#include "inpar_fsi.H"
#include "inpar_immersed.H"
#include "inpar_cell.H"
#include "inpar_fpsi.H"
#include "inpar_xfem.H"
#include "inpar_mortar.H"
#include "inpar_elch.H"
#include "inpar_cardiac_monodomain.H"
#include "inpar_acou.H"
#include "inpar_invanalysis.H"
#include "inpar_levelset.H"
#include "inpar_bio.H"
#include "inpar_cardiovascular0d.H"
#include "inpar_ssi.H"
#include "inpar_ehl.H"
#include "inpar_particle_old.H"
#include "inpar_cavitation.H"
#include "inpar_beampotential.H"
#include "inpar_beaminteraction.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyConditionDefinitions(std::ostream& stream,
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist, bool color)
{
  for (unsigned i = 0; i < condlist.size(); ++i)
  {
    condlist[i]->Print(stream, NULL, color);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintConditionDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>> condlist =
      DRT::INPUT::ValidConditions();
  DRT::INPUT::PrintEmptyConditionDefinitions(std::cout, *condlist);
}



namespace DRT
{
  namespace INPUT
  {
    // collect some problem-specific conditions that do not fit in the generic sections
    void SetMiscellaneousConditions(
        std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
    {
      /*--------------------------------------------------------------------*/
      // Force Sensor

      // declaration of a variable which contains all the components of the condition
      std::vector<Teuchos::RCP<ConditionComponent>> forcesensorcomponents;

      /*the condition consists of one component which has to read an integer value (the number of
       * the dof of the node with respect to which the force is to be measured; note: numbering
       * starts at zero):*/
      forcesensorcomponents.push_back(Teuchos::rcp(new IntConditionComponent("DOF Number")));

      // the condition itself hast to be defined so that it is clear how to read or write such a
      // condition in the dat file
      Teuchos::RCP<ConditionDefinition> forcesensor =
          Teuchos::rcp(new ConditionDefinition("FORCE SENSORS",  // name of input file section
              "ForceSensor",                // name to get the condition from a discretization
              "Force Sensor",               // description of condition
              DRT::Condition::ForceSensor,  // type of condition in DRT (cf. drt_condition.H)
              false,  // should explicit elements be generated (e.g. line elements for line contact
                      // condition in 2D)?
              DRT::Condition::Point));  // type of geometry the condition lives on

      // after definition of the condition all its components have to be added:
      for (unsigned i = 0; i < forcesensorcomponents.size(); ++i)
        forcesensor->AddComponent(forcesensorcomponents[i]);

      // the condition itself has to be added to the condition list
      condlist.push_back(forcesensor);

      /*--------------------------------------------------------------------*/
      // microscale boundary

      Teuchos::RCP<ConditionDefinition> microscale =
          Teuchos::rcp(new ConditionDefinition("MICROSCALE CONDITIONS", "MicroBoundary",
              "Microscale Boundary", DRT::Condition::MicroBoundary, true, DRT::Condition::Surface));

      condlist.push_back(microscale);

      /*--------------------------------------------------------------------*/
      // stc layer condition

      Teuchos::RCP<ConditionDefinition> stclayer = Teuchos::rcp(
          new ConditionDefinition("DESIGN VOL STC LAYER", "STC Layer", "Layer for Multilayered STC",
              DRT::Condition::VolSTCLayer, true, DRT::Condition::Volume));

      stclayer->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

      condlist.push_back(stclayer);
    }
  }  // namespace INPUT
}  // namespace DRT


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>>
DRT::INPUT::ValidConditions()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>> vc =
      Teuchos::rcp(new std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>());

  std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist = *vc;

  /*--------------------------------------------------------------------*/
  // Neumann

  std::vector<Teuchos::RCP<SeparatorConditionComponent>> neumannintsepveccomponents;
  std::vector<Teuchos::RCP<IntVectorConditionComponent>> neumannintveccomponents;
  std::vector<Teuchos::RCP<SeparatorConditionComponent>> neumannrealsepveccomponents;
  std::vector<Teuchos::RCP<RealVectorConditionComponent>> neumannrealveccomponents;
  std::vector<Teuchos::RCP<ConditionComponent>> neumanncomponents;

  neumannintsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  neumannintveccomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff", 1)));
  neumannrealsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  neumannrealveccomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  neumannintsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  neumannintveccomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, true, false)));

  neumanncomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
  neumanncomponents.push_back(Teuchos::rcp(new DirichletNeumannBundle("neumannbund",
      Teuchos::rcp(new IntConditionComponent("numdof")), neumannintsepveccomponents,
      neumannintveccomponents, neumannrealsepveccomponents, neumannrealveccomponents)));

  // optional
  neumanncomponents.push_back(Teuchos::rcp(new StringConditionComponent("type", "Live",
      Teuchos::tuple<std::string>("Live", "Dead", "PrescribedDomainLoad", "constHydro_z",
          "increaseHydro_z", "pseudo_orthopressure", "orthopressure", "LAS", "PressureGrad",
          "Torque"),
      Teuchos::tuple<std::string>("neum_live", "neum_dead", "pres_domain_load", "neum_consthydro_z",
          "neum_increhydro_z", "neum_pseudo_orthopressure", "neum_orthopressure", "neum_LAS",
          "neum_pgrad", "neum_torque"),
      true)));
  neumanncomponents.push_back(Teuchos::rcp(new StringConditionComponent("surface", "Mid",
      Teuchos::tuple<std::string>("Mid", "Top", "Bot"),
      Teuchos::tuple<std::string>("mid", "top", "bot"), true)));

  Teuchos::RCP<ConditionDefinition> pointneumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT NEUMANN CONDITIONS", "PointNeumann",
          "Point Neumann", DRT::Condition::PointNeumann, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> pointneumanneb =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT MOMENT EB CONDITIONS", "PointNeumannEB",
          "Point Neumann Moment auf Euler-Bernoulli Balken", DRT::Condition::PointNeumannEB, false,
          DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineneumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE NEUMANN CONDITIONS", "LineNeumann",
          "Line Neumann", DRT::Condition::LineNeumann, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfneumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF NEUMANN CONDITIONS", "SurfaceNeumann",
          "Surface Neumann", DRT::Condition::SurfaceNeumann, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volneumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL NEUMANN CONDITIONS", "VolumeNeumann",
          "Volume Neumann", DRT::Condition::VolumeNeumann, true, DRT::Condition::Volume));

  // Neumann conditions for transport problems
  Teuchos::RCP<ConditionDefinition> pointtransportneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT TRANSPORT NEUMANN CONDITIONS", "TransportPointNeumann",
          "Point Neumann", DRT::Condition::PointNeumann, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linetransportneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE TRANSPORT NEUMANN CONDITIONS", "TransportLineNeumann",
          "Line Neumann", DRT::Condition::LineNeumann, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF TRANSPORT NEUMANN CONDITIONS", "TransportSurfaceNeumann",
          "Surface Neumann", DRT::Condition::SurfaceNeumann, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> voltransportneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN VOL TRANSPORT NEUMANN CONDITIONS", "TransportVolumeNeumann",
          "Volume Neumann", DRT::Condition::VolumeNeumann, true, DRT::Condition::Volume));

  // Neumann conditions for thermo problems
  Teuchos::RCP<ConditionDefinition> pointthermoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT THERMO NEUMANN CONDITIONS", "ThermoPointNeumann",
          "Point Neumann", DRT::Condition::PointNeumann, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linethermoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE THERMO NEUMANN CONDITIONS", "ThermoLineNeumann",
          "Line Neumann", DRT::Condition::LineNeumann, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfthermoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF THERMO NEUMANN CONDITIONS", "ThermoSurfaceNeumann",
          "Surface Neumann", DRT::Condition::SurfaceNeumann, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volthermoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN VOL THERMO NEUMANN CONDITIONS", "ThermoVolumeNeumann",
          "Volume Neumann", DRT::Condition::VolumeNeumann, true, DRT::Condition::Volume));

  // Neumann conditions for poroelasticity problems
  Teuchos::RCP<ConditionDefinition> pointporoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT PORO NEUMANN CONDITIONS", "PoroPointNeumann",
          "Point Neumann", DRT::Condition::PointNeumann, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineporoneumann =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO NEUMANN CONDITIONS", "PoroLineNeumann",
          "Line Neumann", DRT::Condition::LineNeumann, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfporoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF PORO NEUMANN CONDITIONS", "PoroSurfaceNeumann",
          "Surface Neumann", DRT::Condition::SurfaceNeumann, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volporoneumann = Teuchos::rcp(
      new ConditionDefinition("DESIGN VOL PORO NEUMANN CONDITIONS", "PoroVolumeNeumann",
          "Volume Neumann", DRT::Condition::VolumeNeumann, true, DRT::Condition::Volume));

  for (unsigned i = 0; i < neumanncomponents.size(); ++i)
  {
    pointneumann->AddComponent(neumanncomponents[i]);
    pointneumanneb->AddComponent(neumanncomponents[i]);
    lineneumann->AddComponent(neumanncomponents[i]);
    surfneumann->AddComponent(neumanncomponents[i]);
    volneumann->AddComponent(neumanncomponents[i]);

    pointtransportneumann->AddComponent(neumanncomponents[i]);
    linetransportneumann->AddComponent(neumanncomponents[i]);
    surftransportneumann->AddComponent(neumanncomponents[i]);
    voltransportneumann->AddComponent(neumanncomponents[i]);

    pointthermoneumann->AddComponent(neumanncomponents[i]);
    linethermoneumann->AddComponent(neumanncomponents[i]);
    surfthermoneumann->AddComponent(neumanncomponents[i]);
    volthermoneumann->AddComponent(neumanncomponents[i]);

    pointporoneumann->AddComponent(neumanncomponents[i]);
    lineporoneumann->AddComponent(neumanncomponents[i]);
    surfporoneumann->AddComponent(neumanncomponents[i]);
    volporoneumann->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(pointneumann);
  condlist.push_back(pointneumanneb);
  condlist.push_back(lineneumann);
  condlist.push_back(surfneumann);
  condlist.push_back(volneumann);

  condlist.push_back(pointtransportneumann);
  condlist.push_back(linetransportneumann);
  condlist.push_back(surftransportneumann);

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

  std::vector<Teuchos::RCP<SeparatorConditionComponent>> dirichletintsepveccomponents;
  std::vector<Teuchos::RCP<IntVectorConditionComponent>> dirichletintveccomponents;
  std::vector<Teuchos::RCP<SeparatorConditionComponent>> dirichletrealsepveccomponents;
  std::vector<Teuchos::RCP<RealVectorConditionComponent>> dirichletrealveccomponents;
  std::vector<Teuchos::RCP<ConditionComponent>> dirichletbundcomponents;

  dirichletintsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  dirichletintveccomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff", 1)));
  dirichletrealsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  dirichletrealveccomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  dirichletintsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  dirichletintveccomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, true, false)));

  dirichletbundcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
  dirichletbundcomponents.push_back(Teuchos::rcp(new DirichletNeumannBundle("dirichbund",
      Teuchos::rcp(new IntConditionComponent("numdof")), dirichletintsepveccomponents,
      dirichletintveccomponents, dirichletrealsepveccomponents, dirichletrealveccomponents)));

  // optional
  dirichletbundcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("TAG", true)));
  dirichletbundcomponents.push_back(Teuchos::rcp(new StringConditionComponent("tag", "none",
      Teuchos::tuple<std::string>("none", "monitor_reaction"),
      Teuchos::tuple<std::string>("none", "monitor_reaction"), true)));

  Teuchos::RCP<ConditionDefinition> pointdirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT DIRICH CONDITIONS", "Dirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linedirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE DIRICH CONDITIONS", "Dirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfdirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF DIRICH CONDITIONS", "Dirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> voldirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL DIRICH CONDITIONS", "Dirichlet",
          "Volume Dirichlet", DRT::Condition::VolumeDirichlet, false, DRT::Condition::Volume));

  Teuchos::RCP<ConditionDefinition> particledirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE DIRICH CONDITIONS", "Dirichlet",
          "Particle Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Particle));

  Teuchos::RCP<ConditionDefinition> pointaledirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT ALE DIRICH CONDITIONS", "ALEDirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linealedirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE ALE DIRICH CONDITIONS", "ALEDirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfaledirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF ALE DIRICH CONDITIONS", "ALEDirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volaledirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL ALE DIRICH CONDITIONS", "ALEDirichlet",
          "Volume Dirichlet", DRT::Condition::VolumeDirichlet, false, DRT::Condition::Volume));

  // Dirichlet conditions for transport problems
  Teuchos::RCP<ConditionDefinition> pointtransportdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT TRANSPORT DIRICH CONDITIONS", "TransportDirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linetransportdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE TRANSPORT DIRICH CONDITIONS", "TransportDirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF TRANSPORT DIRICH CONDITIONS", "TransportDirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> voltransportdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN VOL TRANSPORT DIRICH CONDITIONS", "TransportDirichlet",
          "Volume Dirichlet", DRT::Condition::VolumeDirichlet, false, DRT::Condition::Volume));

  // Dirichlet conditions for thermo problems
  Teuchos::RCP<ConditionDefinition> pointthermodirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT THERMO DIRICH CONDITIONS", "ThermoDirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linethermodirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE THERMO DIRICH CONDITIONS", "ThermoDirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfthermodirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF THERMO DIRICH CONDITIONS", "ThermoDirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volthermodirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL THERMO DIRICH CONDITIONS", "ThermoDirichlet",
          "Volume Dirichlet", DRT::Condition::VolumeDirichlet, false, DRT::Condition::Volume));

  // Dirichlet conditions for poroelasticity problems
  Teuchos::RCP<ConditionDefinition> pointporodirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT PORO DIRICH CONDITIONS", "PoroDirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineporodirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO DIRICH CONDITIONS", "PoroDirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfporodirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF PORO DIRICH CONDITIONS", "PoroDirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volporodirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL PORO DIRICH CONDITIONS", "PoroDirichlet",
          "Volume Dirichlet", DRT::Condition::VolumeDirichlet, false, DRT::Condition::Volume));

  Teuchos::RCP<ConditionDefinition> pointnurbslsdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, true, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linenurbslsdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfnurbslsdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volnurbslsdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN VOL NURBS LS DIRICH CONDITIONS", "NurbsLSDirichlet",
          "Volume Dirichlet", DRT::Condition::VolumeDirichlet, true, DRT::Condition::Volume));

  for (unsigned i = 0; i < dirichletbundcomponents.size(); ++i)
  {
    pointdirichlet->AddComponent(dirichletbundcomponents[i]);
    linedirichlet->AddComponent(dirichletbundcomponents[i]);
    surfdirichlet->AddComponent(dirichletbundcomponents[i]);
    voldirichlet->AddComponent(dirichletbundcomponents[i]);

    particledirichlet->AddComponent(dirichletbundcomponents[i]);

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

  condlist.push_back(particledirichlet);

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

  std::vector<Teuchos::RCP<SeparatorConditionComponent>> couplingintsepveccomponents;
  std::vector<Teuchos::RCP<IntVectorConditionComponent>> couplingintveccomponents;
  std::vector<Teuchos::RCP<ConditionComponent>> couplingcomponents;
  std::vector<Teuchos::RCP<SeparatorConditionComponent>> couplingrealsepveccomponents;
  std::vector<Teuchos::RCP<RealVectorConditionComponent>> couplingrealveccomponents;

  couplingintsepveccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  couplingintveccomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff", 1)));

  couplingcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
  couplingcomponents.push_back(Teuchos::rcp(new IntRealBundle("couplingbund",
      Teuchos::rcp(new IntConditionComponent("numdof")), couplingintsepveccomponents,
      couplingintveccomponents, couplingrealsepveccomponents, couplingrealveccomponents)));

  Teuchos::RCP<ConditionDefinition> pointcoupling =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT COUPLING CONDITIONS", "PointCoupling",
          "Point Coupling", DRT::Condition::PointCoupling, false, DRT::Condition::Point));

  Teuchos::RCP<ConditionDefinition> pointthermocoupling = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT THERMO COUPLING CONDITIONS", "PointThermoCoupling",
          "Point Coupling", DRT::Condition::PointCoupling, false, DRT::Condition::Point));

  for (unsigned i = 0; i < couplingcomponents.size(); ++i)
  {
    pointcoupling->AddComponent(couplingcomponents[i]);
    pointthermocoupling->AddComponent(couplingcomponents[i]);
  }

  condlist.push_back(pointcoupling);
  condlist.push_back(pointthermocoupling);

  /*--------------------------------------------------------------------*/
  // Initial fields

  std::vector<Teuchos::RCP<ConditionComponent>> initfieldscomponents;

  initfieldscomponents.push_back(Teuchos::rcp(new StringConditionComponent("Field", "Undefined",
      Teuchos::tuple<std::string>("Undefined", "Velocity", "Pressure", "Temperature", "ScaTra",
          "Porosity", "PoroMultiFluid", "Artery"),
      Teuchos::tuple<std::string>("Undefined", "Velocity", "Pressure", "Temperature", "ScaTra",
          "Porosity", "PoroMultiFluid", "Artery"))));

  // give function id - always one single integer
  // (for initial vector fields, use the COMPONENT option of our functions)
  initfieldscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct", 1)));

  // general initial field conditions
  Teuchos::RCP<ConditionDefinition> pointinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT INITIAL FIELD CONDITIONS", "Initfield",
          "Point Initfield", DRT::Condition::PointInitfield, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE INITIAL FIELD CONDITIONS", "Initfield",
          "Line Initfield", DRT::Condition::LineInitfield, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF INITIAL FIELD CONDITIONS", "Initfield",
          "Surface Initfield", DRT::Condition::SurfaceInitfield, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL INITIAL FIELD CONDITIONS", "Initfield",
          "Volume Initfield", DRT::Condition::VolumeInitfield, false, DRT::Condition::Volume));

  // initial field conditions for particle problems
  Teuchos::RCP<ConditionDefinition> particleinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INITIAL FIELD CONDITIONS", "Initfield",
          "Particle Initfield", DRT::Condition::PointInitfield, false, DRT::Condition::Particle));

  // initial field conditions for temperature
  Teuchos::RCP<ConditionDefinition> pointthermoinitfields = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Point Thermo Initfield", DRT::Condition::PointInitfield, false, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linethermoinitfields = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield",
          "Line Thermo Initfield", DRT::Condition::LineInitfield, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfthermoinitfields = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURF THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield", "Surface Thermo Initfield",
      DRT::Condition::SurfaceInitfield, false, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volthermoinitfields = Teuchos::rcp(new ConditionDefinition(
      "DESIGN VOL THERMO INITIAL FIELD CONDITIONS", "ThermoInitfield", "Volume Thermo Initfield",
      DRT::Condition::VolumeInitfield, false, DRT::Condition::Volume));

  for (unsigned i = 0; i < initfieldscomponents.size(); ++i)
  {
    pointinitfields->AddComponent(initfieldscomponents[i]);
    lineinitfields->AddComponent(initfieldscomponents[i]);
    surfinitfields->AddComponent(initfieldscomponents[i]);
    volinitfields->AddComponent(initfieldscomponents[i]);

    particleinitfields->AddComponent(initfieldscomponents[i]);

    pointthermoinitfields->AddComponent(initfieldscomponents[i]);
    linethermoinitfields->AddComponent(initfieldscomponents[i]);
    surfthermoinitfields->AddComponent(initfieldscomponents[i]);
    volthermoinitfields->AddComponent(initfieldscomponents[i]);
  }

  condlist.push_back(pointinitfields);
  condlist.push_back(lineinitfields);
  condlist.push_back(surfinitfields);
  condlist.push_back(volinitfields);

  condlist.push_back(particleinitfields);

  condlist.push_back(pointthermoinitfields);
  condlist.push_back(linethermoinitfields);
  condlist.push_back(surfthermoinitfields);
  condlist.push_back(volthermoinitfields);

  /*--------------------------------------------------------------------*/
  // compute domain integrals, i.e., cumulative volumes of 3D domain elements or cumulative surface
  // areas of 2D domain elements
  {
    // definition of surface and volume conditions for domain integral computation
    Teuchos::RCP<ConditionDefinition> domainintegralsurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN DOMAIN INTEGRAL SURF CONDITIONS",
            "DomainIntegral", "compute cumulative surface areas of 2D domain elements",
            DRT::Condition::DomainIntegral, true, DRT::Condition::Surface));

    Teuchos::RCP<ConditionDefinition> domainintegralvol =
        Teuchos::rcp(new ConditionDefinition("DESIGN DOMAIN INTEGRAL VOL CONDITIONS",
            "DomainIntegral", "compute cumulative volumes of 3D domain elements",
            DRT::Condition::DomainIntegral, true, DRT::Condition::Volume));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> domainintegralcomponents;

    {
      domainintegralcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
      domainintegralcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
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
    Teuchos::RCP<ConditionDefinition> boundaryintegralsurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN BOUNDARY INTEGRAL SURF CONDITIONS",
            "BoundaryIntegral", "compute cumulative surface areas of 2D boundary elements",
            DRT::Condition::BoundaryIntegral, true, DRT::Condition::Surface));

    // equip condition definition with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> boundaryintegralcomponents;

    {
      boundaryintegralcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
      boundaryintegralcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
    }

    // insert input file line components into condition definition
    for (unsigned i = 0; i < boundaryintegralcomponents.size(); ++i)
      boundaryintegralsurf->AddComponent(boundaryintegralcomponents[i]);

    // insert condition definition into global list of valid condition definitions
    condlist.push_back(boundaryintegralsurf);
  }

  /*--------------------------------------------------------------------*/
  // wear in ALE description

  Teuchos::RCP<ConditionDefinition> linealewear =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE ALE WEAR CONDITIONS 2D", "AleWear",
          "Line Ale Wear", DRT::Condition::AleWear, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfalewear =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE WEAR CONDITIONS 3D", "AleWear",
          "Surface Ale Wear", DRT::Condition::AleWear, true, DRT::Condition::Surface));

  condlist.push_back(linealewear);
  condlist.push_back(surfalewear);

  /*--------------------------------------------------------------------*/
  // local coordinate systems

  std::vector<Teuchos::RCP<ConditionComponent>> locsyscomponents;

  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ROTANGLE")));
  locsyscomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("rotangle", 3)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  locsyscomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 3, false, false)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("USEUPDATEDNODEPOS")));
  locsyscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("useupdatednodepos", 1)));
  locsyscomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("USECONSISTENTNODENORMAL")));
  locsyscomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("useconsistentnodenormal", 1)));

  Teuchos::RCP<ConditionDefinition> pointlocsys = Teuchos::rcp(new ConditionDefinition(
      "DESIGN POINT LOCSYS CONDITIONS", "Locsys", "Point local coordinate system",
      DRT::Condition::PointLocsys, true, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linelocsys =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE LOCSYS CONDITIONS", "Locsys",
          "Line local coordinate system", DRT::Condition::LineLocsys, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surflocsys = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURF LOCSYS CONDITIONS", "Locsys", "Surface local coordinate system",
      DRT::Condition::SurfaceLocsys, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> vollocsys = Teuchos::rcp(new ConditionDefinition(
      "DESIGN VOL LOCSYS CONDITIONS", "Locsys", "Volume local coordinate system",
      DRT::Condition::VolumeLocsys, true, DRT::Condition::Volume));

  // Ale
  Teuchos::RCP<ConditionDefinition> pointalelocsys = Teuchos::rcp(new ConditionDefinition(
      "DESIGN POINT ALE LOCSYS CONDITIONS", "AleLocsys", "Point local coordinate system",
      DRT::Condition::PointLocsys, true, DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linealelocsys =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE ALE LOCSYS CONDITIONS", "AleLocsys",
          "Line local coordinate system", DRT::Condition::LineLocsys, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfalelocsys = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURF ALE LOCSYS CONDITIONS", "AleLocsys", "Surface local coordinate system",
      DRT::Condition::SurfaceLocsys, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volalelocsys = Teuchos::rcp(new ConditionDefinition(
      "DESIGN VOL ALE LOCSYS CONDITIONS", "AleLocsys", "Volume local coordinate system",
      DRT::Condition::VolumeLocsys, true, DRT::Condition::Volume));

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

  std::vector<Teuchos::RCP<ConditionComponent>> pbccomponents;

  pbccomponents.push_back(
      Teuchos::rcp(new IntConditionComponent("Id of periodic boundary condition", true)));
  pbccomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("Is slave periodic boundary condition", "Master",
          Teuchos::tuple<std::string>("Master", "Slave"),
          Teuchos::tuple<std::string>("Master", "Slave"))));
  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("PLANE")));
  pbccomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("degrees of freedom for the pbc plane", "xy",
          Teuchos::tuple<std::string>("xy", "yx", "yz", "zy", "xz", "zx", "xyz"),
          Teuchos::tuple<std::string>("xy", "xy", "yz", "yz", "xz", "xz", "xyz"))));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("LAYER")));
  pbccomponents.push_back(
      Teuchos::rcp(new IntConditionComponent("Layer of periodic boundary condition", true)));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ANGLE")));
  pbccomponents.push_back(Teuchos::rcp(new RealConditionComponent("Angle of rotation")));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ABSTREETOL")));
  pbccomponents.push_back(
      Teuchos::rcp(new RealConditionComponent("Tolerance for nodematching in octree")));

  Teuchos::RCP<ConditionDefinition> lineperiodic = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE PERIODIC BOUNDARY CONDITIONS", "LinePeriodic",
          "Line Periodic", DRT::Condition::LinePeriodic, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfperiodic = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF PERIODIC BOUNDARY CONDITIONS", "SurfacePeriodic",
          "Surface Periodic", DRT::Condition::SurfacePeriodic, false, DRT::Condition::Surface));

  for (unsigned i = 0; i < pbccomponents.size(); ++i)
  {
    lineperiodic->AddComponent(pbccomponents[i]);
    surfperiodic->AddComponent(pbccomponents[i]);
  }

  condlist.push_back(lineperiodic);
  condlist.push_back(surfperiodic);

  /*--------------------------------------------------------------------*/
  // weak Dirichlet conditions

  std::vector<Teuchos::RCP<ConditionComponent>> weakDirichletcomponents;

  // weak DBCs can be imposed adjoint consistent or adjoint inconsistent
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("Choice of gamma parameter", "adjoint-consistent",
          Teuchos::tuple<std::string>("adjoint-consistent", "diffusive-optimal"),
          Teuchos::tuple<std::string>("adjoint-consistent", "diffusive-optimal"))));

  // weak DBCs can be imposed in all directions or only in normal direction
  // (SCATRA: not checked, only in all_directions so far)
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("Directions to apply weak dbc", "all_directions",
          Teuchos::tuple<std::string>("all_directions", "only_in_normal_direction"),
          Teuchos::tuple<std::string>("all_directions", "only_in_normal_direction"))));

  // FLUID: penalty parameter either computed dynamically (using Spaldings law of
  // the wall) or by a fixed value; SCATRA: not checked, only constant value so far
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("Definition of penalty parameter", "constant",
          Teuchos::tuple<std::string>("constant", "Spalding"),
          Teuchos::tuple<std::string>("constant", "Spalding"))));

  // scaling factor for penalty parameter tauB or
  // stabilization parameter alpha for Nitsche term
  // (SCATRA: if stabilization parameter negative -> mixed-hybrid formulation)
  weakDirichletcomponents.push_back(Teuchos::rcp(new RealConditionComponent("TauBscaling")));

  // linearisation strategies --- the linearisation (i.e. the matrix
  // contribution) of the convective term on the inflow could be
  // suppressed, since the flux is a kink function and including this one
  // might result in even worse convergence behaviour
  // (SCATRA: not checked)
  weakDirichletcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Linearisation",
      "lin_all", Teuchos::tuple<std::string>("lin_all", "no_lin_conv_inflow"),
      Teuchos::tuple<std::string>("lin_all", "no_lin_conv_inflow"))));

  // we provide a vector of 3 values for velocities
  weakDirichletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val", 3)));

  // and optional spatial functions
  weakDirichletcomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 3, false, false, true)));


  Teuchos::RCP<ConditionDefinition> lineweakdirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE WEAK DIRICHLET CONDITIONS", "LineWeakDirichlet",
          "LineWeakDirichlet", DRT::Condition::LineWeakDirichlet, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfweakdirichlet = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE WEAK DIRICHLET CONDITIONS", "SurfaceWeakDirichlet", "SurfaceWeakDirichlet",
      DRT::Condition::SurfaceWeakDirichlet, true, DRT::Condition::Surface));

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

  std::vector<Teuchos::RCP<ConditionComponent>> sprcomponents;

  Teuchos::RCP<ConditionDefinition> linespr = Teuchos::rcp(
      new ConditionDefinition("DESIGN PATCH RECOVERY BOUNDARY LINE CONDITIONS", "SPRboundary",
          "Boundary for SPR", DRT::Condition::SPRboundary, false, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfspr = Teuchos::rcp(
      new ConditionDefinition("DESIGN PATCH RECOVERY BOUNDARY SURF CONDITIONS", "SPRboundary",
          "Boundary for SPR", DRT::Condition::SPRboundary, false, DRT::Condition::Surface));

  for (unsigned i = 0; i < sprcomponents.size(); ++i)
  {
    linespr->AddComponent(sprcomponents[i]);
    surfspr->AddComponent(sprcomponents[i]);
  }

  condlist.push_back(linespr);
  condlist.push_back(surfspr);

  /*--------------------------------------------------------------------*/
  // volume constraint

  Teuchos::RCP<ConditionDefinition> volumeconstraint = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE VOLUME CONSTRAINT 3D", "VolumeConstraint_3D", "Surface Volume Constraint",
      DRT::Condition::VolumeConstraint_3D, true, DRT::Condition::Surface));

  volumeconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  volumeconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  volumeconstraint->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  volumeconstraint->AddComponent(Teuchos::rcp(new StringConditionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(volumeconstraint);

  /*--------------------------------------------------------------------*/
  // volume constraint penalty

  Teuchos::RCP<ConditionDefinition> volumeconstraintpen =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE VOLUME CONSTRAINT 3D PEN",
          "VolumeConstraint_3D_Pen", "Surface Volume Constraint Penalty",
          DRT::Condition::VolumeConstraint_3D_pen, true, DRT::Condition::Surface));

  volumeconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("penalty")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("rho")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new StringConditionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(volumeconstraintpen);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<ConditionDefinition> areaconstraint = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE AREA CONSTRAINT 3D", "AreaConstraint_3D", "Surface Area Constraint",
      DRT::Condition::AreaConstraint_3D, true, DRT::Condition::Surface));

  areaconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areaconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  areaconstraint->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));

  condlist.push_back(areaconstraint);

  /*--------------------------------------------------------------------*/
  // area constraint penalty

  Teuchos::RCP<ConditionDefinition> areaconstraintpen =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D PEN",
          "AreaConstraint_3D_Pen", "Surface Area Constraint Penalty",
          DRT::Condition::AreaConstraint_3D_pen, true, DRT::Condition::Surface));

  areaconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("penalty")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("rho")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new StringConditionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(areaconstraintpen);


  /*--------------------------------------------------------------------*/
  // volume monitor

  Teuchos::RCP<ConditionDefinition> volumemonitor = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE VOLUME MONITOR 3D", "VolumeMonitor_3D", "Surface Volume Monitor",
      DRT::Condition::VolumeMonitor_3D, true, DRT::Condition::Surface));

  volumemonitor->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(volumemonitor);

  /*--------------------------------------------------------------------*/
  // area monitor 3D

  Teuchos::RCP<ConditionDefinition> areamonitor =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE AREA MONITOR 3D", "AreaMonitor_3D",
          "Surface Area Monitor", DRT::Condition::AreaMonitor_3D, true, DRT::Condition::Surface));

  areamonitor->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areamonitor->AddComponent(Teuchos::rcp(new StringConditionComponent("projection", "none",
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"),
      Teuchos::tuple<std::string>("none", "xy", "yz", "xz"), true)));

  condlist.push_back(areamonitor);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<ConditionDefinition> areaconstraint2D =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE AREA CONSTRAINT 2D", "AreaConstraint_2D",
          "Line Area Constraint", DRT::Condition::AreaConstraint_2D, true, DRT::Condition::Line));

  areaconstraint2D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areaconstraint2D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  areaconstraint2D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  condlist.push_back(areaconstraint2D);

  /*--------------------------------------------------------------------*/
  // area monitor 2D

  Teuchos::RCP<ConditionDefinition> areamonitor2D =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE AREA MONITOR 2D", "AreaMonitor_2D",
          "Line Area Monitor", DRT::Condition::AreaMonitor_2D, true, DRT::Condition::Line));

  areamonitor2D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  condlist.push_back(areamonitor2D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D for a node over a plane

  Teuchos::RCP<ConditionDefinition> nodeonplaneconst3D = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE MULTIPNT CONSTRAINT 3D", "MPC_NodeOnPlane_3D", "Node on Plane Constraint",
      DRT::Condition::MPC_NodeOnPlane_3D, false, DRT::Condition::Surface));

  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("planeNodes", 3)));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new StringConditionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodeonplaneconst3D);

  /*--------------------------------------------------------------------*/
  //  // Multi point constraint in 3D for a node over a plane
  //
  //  Teuchos::RCP<ConditionDefinition> nodeonlineconst3D =
  //    Teuchos::rcp(new ConditionDefinition("DESIGN LINE MULTIPNT CONSTRAINT 3D",
  //                                         "MPC_NodeOnLine_3D",
  //                                         "Node on Line Constraint",
  //                                         DRT::Condition::MPC_NodeOnLine_3D,
  //                                         false,
  //                                         DRT::Condition::Line));
  //
  //  nodeonlineconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  //  nodeonlineconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  //  nodeonlineconst3D->AddComponent(Teuchos::rcp(new
  //  IntVectorConditionComponent("curveNodes",2))); condlist.push_back(nodeonlineconst3D);
  //
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously

  Teuchos::RCP<ConditionDefinition> nodemasterconst3D =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D",
          "MPC_NormalComponent_3D", "Node on Plane Constraint",
          DRT::Condition::MPC_NormalComponent_3D, false, DRT::Condition::Surface));

  nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("masterNode")));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("direction", 3)));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new StringConditionComponent("value", "disp",
      Teuchos::tuple<std::string>("disp", "x"), Teuchos::tuple<std::string>("disp", "x"), true)));
  nodemasterconst3D->AddComponent(Teuchos::rcp(new StringConditionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodemasterconst3D);

  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously, penalty based

  Teuchos::RCP<ConditionDefinition> nodemasterconst3Dpen =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D PEN",
          "MPC_NormalComponent_3D_Pen", "Node on Plane Constraint Penalty",
          DRT::Condition::MPC_NormalComponent_3D_pen, false, DRT::Condition::Surface));

  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealConditionComponent("penalty")));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntConditionComponent("masterNode")));
  nodemasterconst3Dpen->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("direction", 3)));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new StringConditionComponent("value", "disp",
      Teuchos::tuple<std::string>("disp", "x"), Teuchos::tuple<std::string>("disp", "x"), true)));
  nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new StringConditionComponent("control", "rel",
      Teuchos::tuple<std::string>("rel", "abs"), Teuchos::tuple<std::string>("rel", "abs"), true)));
  condlist.push_back(nodemasterconst3Dpen);
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 2D for a node on a line
  Teuchos::RCP<ConditionDefinition> nodeonlineconst2D = Teuchos::rcp(new ConditionDefinition(
      "DESIGN LINE MULTIPNT CONSTRAINT 2D", "MPC_NodeOnLine_2D", "Node on Line Constraint",
      DRT::Condition::MPC_NodeOnLine_2D, false, DRT::Condition::Line));

  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("constrNode 1")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("constrNode 2")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("constrNode 3")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new StringConditionComponent("control value", "dist",
      Teuchos::tuple<std::string>("dist", "angle"), Teuchos::tuple<std::string>("dist", "angle"),
      true)));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
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


    std::vector<Teuchos::RCP<ConditionComponent>> rigidbodymodecomponents;

    rigidbodymodecomponents.push_back(Teuchos::rcp(new StringConditionComponent("discretization",
        "fluid", Teuchos::tuple<std::string>("fluid", "scatra", "solid"),
        Teuchos::tuple<std::string>("fluid", "scatra", "solid"))));

    // input: trigger as int-vector using an IntRealBundle
    // definition separator for int vectors
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepveccomp;
    intsepveccomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
    // definition int vectors
    std::vector<Teuchos::RCP<IntVectorConditionComponent>> intveccomp;
    intveccomp.push_back(Teuchos::rcp(new IntVectorConditionComponent("ONOFF", 1)));
    // definition separator for real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepveccomp;
    // definition real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<RealVectorConditionComponent>> realveccomp;

    rigidbodymodecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMMODES")));
    rigidbodymodecomponents.push_back(Teuchos::rcp(
        new IntRealBundle("intreal bundle", Teuchos::rcp(new IntConditionComponent("NUMMODES")),
            intsepveccomp, intveccomp, realsepveccomp, realveccomp)));

    rigidbodymodecomponents.push_back(
        Teuchos::rcp(new StringConditionComponent("weight vector definition", "integration",
            Teuchos::tuple<std::string>("integration", "pointvalues"),
            Teuchos::tuple<std::string>("integration", "pointvalues"))));

    Teuchos::RCP<ConditionDefinition> surfrigidbodymode =
        Teuchos::rcp(new ConditionDefinition("DESIGN SURF MODE FOR KRYLOV SPACE PROJECTION",
            "KrylovSpaceProjection", "Surface mode for Krylov space projection",
            DRT::Condition::SurfaceModeKrylovProjection, true, DRT::Condition::Surface));

    Teuchos::RCP<ConditionDefinition> volrigidbodymode =
        Teuchos::rcp(new ConditionDefinition("DESIGN VOL MODE FOR KRYLOV SPACE PROJECTION",
            "KrylovSpaceProjection", "Volume mode for Krylov space projection",
            DRT::Condition::VolumeModeKrylovProjection, true, DRT::Condition::Volume));

    for (unsigned i = 0; i < rigidbodymodecomponents.size(); ++i)
    {
      surfrigidbodymode->AddComponent(rigidbodymodecomponents[i]);
      volrigidbodymode->AddComponent(rigidbodymodecomponents[i]);
    }
    condlist.push_back(surfrigidbodymode);
    condlist.push_back(volrigidbodymode);
  }

  // Finally, add the problem-specific conditions from the various modules

  INPAR::MORTAR::SetValidConditions(condlist);

  INPAR::S2I::SetValidConditions(condlist);

  INPAR::SCATRA::SetValidConditions(condlist);

  INPAR::STI::SetValidConditions(condlist);

  INPAR::ELCH::SetValidConditions(condlist);

  INPAR::EP::SetValidConditions(condlist);

  INPAR::FLUID::SetValidConditions(condlist);

  INPAR::ALE::SetValidConditions(condlist);

  INPAR::FSI::SetValidConditions(condlist);

  INPAR::FPSI::SetValidConditions(condlist);

  INPAR::IMMERSED::SetValidConditions(condlist);

  INPAR::CELL::SetValidConditions(condlist);

  INPAR::XFEM::SetValidConditions(dirichletbundcomponents, neumanncomponents, condlist);

  INPAR::INVANA::SetValidConditions(condlist);

  INPAR::BIOFILM::SetValidConditions(condlist);

  INPAR::ARTNET::SetValidConditions(condlist);

  INPAR::REDAIRWAYS::SetValidConditions(condlist);

  INPAR::CARDIOVASCULAR0D::SetValidConditions(condlist);

  INPAR::STR::SetValidConditions(condlist);

  INPAR::THR::SetValidConditions(condlist);

  INPAR::SSI::SetValidConditions(condlist);

  INPAR::PARTICLEOLD::SetValidConditions(condlist);

  INPAR::PATSPEC::SetValidConditions(condlist);

  INPAR::LEVELSET::SetValidConditions(condlist);

  INPAR::ACOU::SetValidConditions(condlist);

  INPAR::BEAMPOTENTIAL::SetValidConditions(condlist);

  INPAR::BEAMINTERACTION::SetValidConditions(condlist);

  INPAR::EHL::SetValidConditions(condlist);

  // finally some conditions that do not have their own files yet are problem-specific
  SetMiscellaneousConditions(condlist);

  return vc;
}
