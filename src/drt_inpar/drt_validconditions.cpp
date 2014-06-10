/*----------------------------------------------------------------------*/
/*!
\file drt_validconditions.cpp

\brief Setup of the list of valid conditions for input

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/


#include "drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "../drt_inpar/inpar_scatra.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyConditionDefinitions(std::ostream& stream,
                                                std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist,
                                                bool color)
{
  for (unsigned i=0; i<condlist.size(); ++i)
  {
    condlist[i]->Print(stream,NULL,color);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintConditionDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > condlist = DRT::INPUT::ValidConditions();
  DRT::INPUT::PrintEmptyConditionDefinitions(std::cout,*condlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > DRT::INPUT::ValidConditions()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > vc =
    Teuchos::rcp(new std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >());

  std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist = *vc;

  /*--------------------------------------------------------------------*/
  // Neumann

  std::vector<Teuchos::RCP<ConditionComponent> > neumanncomponents;
  neumanncomponents.push_back(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  neumanncomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",6)));
  neumanncomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",6)));

  // optional
  neumanncomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "type","Live",
        Teuchos::tuple<std::string>("Live","Dead","PrescribedDomainLoad","constHydro_z","increaseHydro_z","pseudo_orthopressure","orthopressure","LAS","PressureGrad"),
        Teuchos::tuple<std::string>("neum_live","neum_dead","pres_domain_load","neum_consthydro_z","neum_increhydro_z","neum_pseudo_orthopressure","neum_orthopressure","neum_LAS","neum_pgrad"),
        true)));
  neumanncomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "surface","Mid",
        Teuchos::tuple<std::string>("Mid","Top","Bot"),
        Teuchos::tuple<std::string>("mid","top","bot"),
        true)));

  // optional
  neumanncomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",6,false,false,true)));

  Teuchos::RCP<ConditionDefinition> pointneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT NEUMANN CONDITIONS",
                                         "PointNeumann",
                                         "Point Neumann",
                                         DRT::Condition::PointNeumann,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> pointneumanneb =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT MOMENT EB CONDITIONS",
                                         "PointNeumannEB",
                                         "Point Neumann Moment auf Euler-Bernoulli Balken",
                                         DRT::Condition::PointNeumannEB,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE NEUMANN CONDITIONS",
                                         "LineNeumann",
                                         "Line Neumann",
                                         DRT::Condition::LineNeumann,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF NEUMANN CONDITIONS",
                                         "SurfaceNeumann",
                                         "Surface Neumann",
                                         DRT::Condition::SurfaceNeumann,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL NEUMANN CONDITIONS",
                                         "VolumeNeumann",
                                         "Volume Neumann",
                                         DRT::Condition::VolumeNeumann,
                                         true,
                                         DRT::Condition::Volume));

  // Neumann conditions for transport problems
  Teuchos::RCP<ConditionDefinition> pointtransportneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT TRANSPORT NEUMANN CONDITIONS",
                                         "TransportPointNeumann",
                                         "Point Neumann",
                                         DRT::Condition::PointNeumann,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linetransportneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE TRANSPORT NEUMANN CONDITIONS",
                                         "TransportLineNeumann",
                                         "Line Neumann",
                                         DRT::Condition::LineNeumann,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF TRANSPORT NEUMANN CONDITIONS",
                                         "TransportSurfaceNeumann",
                                         "Surface Neumann",
                                         DRT::Condition::SurfaceNeumann,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> voltransportneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL TRANSPORT NEUMANN CONDITIONS",
                                         "TransportVolumeNeumann",
                                         "Volume Neumann",
                                         DRT::Condition::VolumeNeumann,
                                         true,
                                         DRT::Condition::Volume));

  // Neumann conditions for thermo problems
  Teuchos::RCP<ConditionDefinition> pointthermoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT THERMO NEUMANN CONDITIONS",
                                         "ThermoPointNeumann",
                                         "Point Neumann",
                                         DRT::Condition::PointNeumann,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linethermoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE THERMO NEUMANN CONDITIONS",
                                         "ThermoLineNeumann",
                                         "Line Neumann",
                                         DRT::Condition::LineNeumann,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfthermoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF THERMO NEUMANN CONDITIONS",
                                         "ThermoSurfaceNeumann",
                                         "Surface Neumann",
                                         DRT::Condition::SurfaceNeumann,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volthermoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL THERMO NEUMANN CONDITIONS",
                                         "ThermoVolumeNeumann",
                                         "Volume Neumann",
                                         DRT::Condition::VolumeNeumann,
                                         true,
                                         DRT::Condition::Volume));

  // Neumann conditions for xfem fluid problems
  Teuchos::RCP<ConditionDefinition> pointXFEMneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN XFEM POINT NEUMANN CONDITIONS",
                                         "PointXFEMNeumann",
                                         "Point XFEM Neumann",
                                         DRT::Condition::PointNeumann,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineXFEMneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN XFEM LINE NEUMANN CONDITIONS",
                                         "LineXFEMNeumann",
                                         "Line XFEM Neumann",
                                         DRT::Condition::LineNeumann,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfXFEMneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN XFEM SURF NEUMANN CONDITIONS",
                                         "SurfaceXFEMNeumann",
                                         "Surface XFEM Neumann",
                                         DRT::Condition::SurfaceNeumann,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volXFEMneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN XFEM VOL NEUMANN CONDITIONS",
                                         "VolumeXFEMNeumann",
                                         "Volume XFEM Neumann",
                                         DRT::Condition::VolumeNeumann,
                                         true,
                                         DRT::Condition::Volume));

  // Neumann conditions for poroelasticity problems
  Teuchos::RCP<ConditionDefinition> pointporoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT PORO NEUMANN CONDITIONS",
                                         "PoroPointNeumann",
                                         "Point Neumann",
                                         DRT::Condition::PointNeumann,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineporoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO NEUMANN CONDITIONS",
                                         "PoroLineNeumann",
                                         "Line Neumann",
                                         DRT::Condition::LineNeumann,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfporoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF PORO NEUMANN CONDITIONS",
                                         "PoroSurfaceNeumann",
                                         "Surface Neumann",
                                         DRT::Condition::SurfaceNeumann,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volporoneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL PORO NEUMANN CONDITIONS",
                                         "PoroVolumeNeumann",
                                         "Volume Neumann",
                                         DRT::Condition::VolumeNeumann,
                                         true,
                                         DRT::Condition::Volume));

  for (unsigned i=0; i<neumanncomponents.size(); ++i)
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

    pointXFEMneumann->AddComponent(neumanncomponents[i]);
    lineXFEMneumann->AddComponent(neumanncomponents[i]);
    surfXFEMneumann->AddComponent(neumanncomponents[i]);
    volXFEMneumann->AddComponent(neumanncomponents[i]);

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

  condlist.push_back(pointXFEMneumann);
  condlist.push_back(lineXFEMneumann);
  condlist.push_back(surfXFEMneumann);
  condlist.push_back(volXFEMneumann);

  condlist.push_back(pointporoneumann);
  condlist.push_back(lineporoneumann);
  condlist.push_back(surfporoneumann);
  condlist.push_back(volporoneumann);

  /*--------------------------------------------------------------------*/
  // Dirichlet

  std::vector<Teuchos::RCP<SeparatorConditionComponent> > dirichletintsepveccomponents;
  std::vector<Teuchos::RCP<IntVectorConditionComponent> > dirichletintveccomponents;
  std::vector<Teuchos::RCP<SeparatorConditionComponent> > dirichletrealsepveccomponents;
  std::vector<Teuchos::RCP<RealVectorConditionComponent> > dirichletrealveccomponents;
  std::vector<Teuchos::RCP<ConditionComponent> > dirichletbundcomponents;

  dirichletintsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  dirichletintveccomponents.push_back(
    Teuchos::rcp(new IntVectorConditionComponent("onoff", 1)));
  dirichletrealsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  dirichletrealveccomponents.push_back(
    Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  dirichletintsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("CURVE")));
  dirichletintveccomponents.push_back(
    Teuchos::rcp(new IntVectorConditionComponent("curve", 1, true, true)));
  dirichletintsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("FUNCT",true)));
  dirichletintveccomponents.push_back(
    Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  dirichletbundcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
  dirichletbundcomponents.push_back(
      Teuchos::rcp(
        new DirichletNeumannBundle(
            "dirichbund",
            Teuchos::rcp(new IntConditionComponent("numdof")),
            dirichletintsepveccomponents,
            dirichletintveccomponents,
            dirichletrealsepveccomponents,
            dirichletrealveccomponents)));

  Teuchos::RCP<ConditionDefinition> pointdirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT DIRICH CONDITIONS",
                                         "Dirichlet",
                                         "Point Dirichlet",
                                         DRT::Condition::PointDirichlet,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linedirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE DIRICH CONDITIONS",
                                         "Dirichlet",
                                         "Line Dirichlet",
                                         DRT::Condition::LineDirichlet,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfdirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF DIRICH CONDITIONS",
                                         "Dirichlet",
                                         "Surface Dirichlet",
                                         DRT::Condition::SurfaceDirichlet,
                                         false,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> voldirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL DIRICH CONDITIONS",
                                         "Dirichlet",
                                         "Volume Dirichlet",
                                         DRT::Condition::VolumeDirichlet,
                                         false,
                                         DRT::Condition::Volume));

  Teuchos::RCP<ConditionDefinition> particledirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE DIRICH CONDITIONS",
                                         "Dirichlet",
                                         "Particle Dirichlet",
                                         DRT::Condition::PointDirichlet,
                                         false,
                                         DRT::Condition::Particle));

  Teuchos::RCP<ConditionDefinition> pointaledirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT ALE DIRICH CONDITIONS",
                                         "ALEDirichlet",
                                         "Point Dirichlet",
                                         DRT::Condition::PointDirichlet,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linealedirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE ALE DIRICH CONDITIONS",
                                         "ALEDirichlet",
                                         "Line Dirichlet",
                                         DRT::Condition::LineDirichlet,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfaledirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF ALE DIRICH CONDITIONS",
                                         "ALEDirichlet",
                                         "Surface Dirichlet",
                                         DRT::Condition::SurfaceDirichlet,
                                         false,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volaledirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL ALE DIRICH CONDITIONS",
                                         "ALEDirichlet",
                                         "Volume Dirichlet",
                                         DRT::Condition::VolumeDirichlet,
                                         false,
                                         DRT::Condition::Volume));

  // Dirichlet conditions for transport problems
  Teuchos::RCP<ConditionDefinition> pointtransportdirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT TRANSPORT DIRICH CONDITIONS",
                                         "TransportDirichlet",
                                         "Point Dirichlet",
                                         DRT::Condition::PointDirichlet,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linetransportdirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE TRANSPORT DIRICH CONDITIONS",
                                         "TransportDirichlet",
                                         "Line Dirichlet",
                                         DRT::Condition::LineDirichlet,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportdirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF TRANSPORT DIRICH CONDITIONS",
                                         "TransportDirichlet",
                                         "Surface Dirichlet",
                                         DRT::Condition::SurfaceDirichlet,
                                         false,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> voltransportdirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL TRANSPORT DIRICH CONDITIONS",
                                         "TransportDirichlet",
                                         "Volume Dirichlet",
                                         DRT::Condition::VolumeDirichlet,
                                         false,
                                         DRT::Condition::Volume));

  // Dirichlet conditions for thermo problems
  Teuchos::RCP<ConditionDefinition> pointthermodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT THERMO DIRICH CONDITIONS",
                                         "ThermoDirichlet",
                                         "Point Dirichlet",
                                         DRT::Condition::PointDirichlet,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linethermodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE THERMO DIRICH CONDITIONS",
                                         "ThermoDirichlet",
                                         "Line Dirichlet",
                                         DRT::Condition::LineDirichlet,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfthermodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF THERMO DIRICH CONDITIONS",
                                         "ThermoDirichlet",
                                         "Surface Dirichlet",
                                         DRT::Condition::SurfaceDirichlet,
                                         false,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volthermodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL THERMO DIRICH CONDITIONS",
                                         "ThermoDirichlet",
                                         "Volume Dirichlet",
                                         DRT::Condition::VolumeDirichlet,
                                         false,
                                         DRT::Condition::Volume));

  // Dirichlet conditions for poroelasticity problems
  Teuchos::RCP<ConditionDefinition> pointporodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT PORO DIRICH CONDITIONS",
                                         "PoroDirichlet",
                                         "Point Dirichlet",
                                         DRT::Condition::PointDirichlet,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineporodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO DIRICH CONDITIONS",
                                         "PoroDirichlet",
                                         "Line Dirichlet",
                                         DRT::Condition::LineDirichlet,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfporodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF PORO DIRICH CONDITIONS",
                                         "PoroDirichlet",
                                         "Surface Dirichlet",
                                         DRT::Condition::SurfaceDirichlet,
                                         false,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volporodirichlet =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL PORO DIRICH CONDITIONS",
                                         "PoroDirichlet",
                                         "Volume Dirichlet",
                                         DRT::Condition::VolumeDirichlet,
                                         false,
                                         DRT::Condition::Volume));


  for (unsigned i=0; i<dirichletbundcomponents.size(); ++i)
  {
    pointdirichlet->AddComponent(dirichletbundcomponents[i]);
    linedirichlet->AddComponent(dirichletbundcomponents[i]);
    surfdirichlet->AddComponent(dirichletbundcomponents[i]);
    voldirichlet->AddComponent(dirichletbundcomponents[i]);

    particledirichlet->AddComponent(dirichletbundcomponents[i]);

    pointaledirichlet->AddComponent(dirichletbundcomponents[i]);
    linealedirichlet ->AddComponent(dirichletbundcomponents[i]);
    surfaledirichlet ->AddComponent(dirichletbundcomponents[i]);
    volaledirichlet  ->AddComponent(dirichletbundcomponents[i]);

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


  /*--------------------------------------------------------------------*/
  // Biofilm growth Dirichlet

  std::vector<Teuchos::RCP<SeparatorConditionComponent> > biodirichletintsepveccomponents;
  std::vector<Teuchos::RCP<IntVectorConditionComponent> > biodirichletintveccomponents;
  std::vector<Teuchos::RCP<SeparatorConditionComponent> > biodirichletrealsepveccomponents;
  std::vector<Teuchos::RCP<RealVectorConditionComponent> > biodirichletrealveccomponents;
  std::vector<Teuchos::RCP<ConditionComponent> > biodirichletbundcomponents;

  biodirichletintsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  biodirichletintveccomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("onoff", 1)));
  biodirichletrealsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  biodirichletrealveccomponents.push_back(
      Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  biodirichletintsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("CURVE")));
  biodirichletintveccomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("curve", 1, true, true)));
  biodirichletintsepveccomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("FUNCT",true)));
  biodirichletintveccomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, false, true)));

  biodirichletbundcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
  biodirichletbundcomponents.push_back(
      Teuchos::rcp(
          new DirichletNeumannBundle(
              "dirichbund",
              Teuchos::rcp(new IntConditionComponent("numdof")),
              biodirichletintsepveccomponents,
              biodirichletintveccomponents,
              biodirichletrealsepveccomponents,
              biodirichletrealveccomponents)));

  // Dirichlet conditions for biofilm growth problems
  Teuchos::RCP<ConditionDefinition> pointbiofilmdirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT BIOFILM DIRICH CONDITIONS",
          "BioDirichlet",
          "Point BioDirichlet",
          DRT::Condition::PointDirichlet,
          false,
          DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linebiofilmdirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE BIOFILM DIRICH CONDITIONS",
          "BioDirichlet",
          "Line BioDirichlet",
          DRT::Condition::LineDirichlet,
          false,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfbiofilmdirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF BIOFILM DIRICH CONDITIONS",
          "BioDirichlet",
          "Surface BioDirichlet",
          DRT::Condition::SurfaceDirichlet,
          false,
          DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volbiofilmdirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL BIOFILM DIRICH CONDITIONS",
          "BioDirichlet",
          "Volume BioDirichlet",
          DRT::Condition::VolumeDirichlet,
          false,
          DRT::Condition::Volume));
  for (unsigned i=0; i<dirichletbundcomponents.size(); ++i)
  {
    pointbiofilmdirichlet->AddComponent(biodirichletbundcomponents[i]);
    linebiofilmdirichlet->AddComponent(biodirichletbundcomponents[i]);
    surfbiofilmdirichlet->AddComponent(biodirichletbundcomponents[i]);
    volbiofilmdirichlet->AddComponent(biodirichletbundcomponents[i]);
  }


  condlist.push_back(pointbiofilmdirichlet);
  condlist.push_back(linebiofilmdirichlet);
  condlist.push_back(surfbiofilmdirichlet);
  condlist.push_back(volbiofilmdirichlet);


  /*--------------------------------------------------------------------*/
  // Initial fields

  std::vector<Teuchos::RCP<ConditionComponent> > initfieldscomponents;

  initfieldscomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Field","Undefined",
        Teuchos::tuple<std::string>("Undefined","Velocity","Pressure","Temperature","ScaTra","Porosity"),
        Teuchos::tuple<std::string>("Undefined","Velocity","Pressure","Temperature","ScaTra","Porosity"))));

  // give function id - always one single integer
  // (for initial vector fields, use the COMPONENT option of our functions)
  initfieldscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",1)));

  Teuchos::RCP<ConditionDefinition> pointinitfields =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT INITIAL FIELD CONDITIONS",
                                         "Initfield",
                                         "Point Initfield",
                                         DRT::Condition::PointInitfield,
                                         false,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> lineinitfields =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE INITIAL FIELD CONDITIONS",
                                         "Initfield",
                                         "Line Initfield",
                                         DRT::Condition::LineInitfield,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfinitfields =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF INITIAL FIELD CONDITIONS",
                                         "Initfield",
                                         "Surface Initfield",
                                         DRT::Condition::SurfaceInitfield,
                                         false,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volinitfields =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL INITIAL FIELD CONDITIONS",
                                         "Initfield",
                                         "Volume Initfield",
                                         DRT::Condition::VolumeInitfield,
                                         false,
                                         DRT::Condition::Volume));

  Teuchos::RCP<ConditionDefinition> particleinitfields =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INITIAL FIELD CONDITIONS",
                                         "Initfield",
                                         "Particle Initfield",
                                         DRT::Condition::PointInitfield,
                                         false,
                                         DRT::Condition::Particle));

  for (unsigned i=0; i<initfieldscomponents.size(); ++i)
  {
    pointinitfields->AddComponent(initfieldscomponents[i]);
    lineinitfields->AddComponent(initfieldscomponents[i]);
    surfinitfields->AddComponent(initfieldscomponents[i]);
    volinitfields->AddComponent(initfieldscomponents[i]);

    particleinitfields->AddComponent(initfieldscomponents[i]);
  }

  condlist.push_back(pointinitfields);
  condlist.push_back(lineinitfields);
  condlist.push_back(surfinitfields);
  condlist.push_back(volinitfields);

  condlist.push_back(particleinitfields);

  /*--------------------------------------------------------------------*/
  // mortar coupling (for ALL kinds of interface problems)

  std::vector<Teuchos::RCP<ConditionComponent> > mortarcomponents;

  mortarcomponents.push_back(Teuchos::rcp(new IntConditionComponent("Interface ID")));
  mortarcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Side","Master",
        Teuchos::tuple<std::string>("Master","Slave","Selfcontact"),
        Teuchos::tuple<std::string>("Master","Slave","Selfcontact"))));
  mortarcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Initialization","Inactive",
        Teuchos::tuple<std::string>("Inactive","Active"),
        Teuchos::tuple<std::string>("Inactive","Active"),true)));

  mortarcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FrCoeffOrBound",true)));
  mortarcomponents.push_back(Teuchos::rcp(new RealConditionComponent("FrCoeffOrBound")));

  Teuchos::RCP<ConditionDefinition> linemortar =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE MORTAR COUPLING CONDITIONS 2D",
                                         "Mortar",
                                         "Line Mortar Coupling",
                                         DRT::Condition::Mortar,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfmortar =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF MORTAR COUPLING CONDITIONS 3D",
                                         "Mortar",
                                         "Surface Mortar Coupling",
                                         DRT::Condition::Mortar,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<mortarcomponents.size(); ++i)
  {
    linemortar->AddComponent(mortarcomponents[i]);
    surfmortar->AddComponent(mortarcomponents[i]);
  }

  condlist.push_back(linemortar);
  condlist.push_back(surfmortar);

  /*--------------------------------------------------------------------*/
  // mortar coupling symmetry condition

  std::vector<Teuchos::RCP<ConditionComponent> > mrtrsymcomponents;
  mrtrsymcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  mrtrsymcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",3)));

  Teuchos::RCP<ConditionDefinition> linemrtrsym =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE MORTAR SYMMETRY CONDITIONS 3D",
                                         "mrtrsym",
                                         "Symmetry plane normal for 3D contact",
                                         DRT::Condition::LineMrtrSym,
                                         true,
                                         DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> pointmrtrsym =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT MORTAR SYMMETRY CONDITIONS 2D/3D",
                                         "mrtrsym",
                                         "Symmetry plane normal for 2D/3D contact",
                                         DRT::Condition::PointMrtrSym,
                                         true,
                                         DRT::Condition::Point));

  for (unsigned i=0; i<mrtrsymcomponents.size(); ++i)
  {
    linemrtrsym->AddComponent(mrtrsymcomponents[i]);
    pointmrtrsym->AddComponent(mrtrsymcomponents[i]);
  }

  condlist.push_back(linemrtrsym);
  condlist.push_back(pointmrtrsym);




  /*--------------------------------------------------------------------*/
  // wear in ALE description

  Teuchos::RCP<ConditionDefinition> linealewear =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE ALE WEAR CONDITIONS 2D",
                                         "AleWear",
                                         "Line Ale Wear",
                                         DRT::Condition::AleWear,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfalewear =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE WEAR CONDITIONS 3D",
                                         "AleWear",
                                         "Surface Ale Wear",
                                         DRT::Condition::AleWear,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(linealewear);
  condlist.push_back(surfalewear);

  /*--------------------------------------------------------------------*/
  // local coordinate systems

  std::vector<Teuchos::RCP<ConditionComponent> > locsyscomponents;

  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ROTANGLE")));
  locsyscomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("rotangle",3)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("CURVE")));
  locsyscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve", 3, true, true)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  locsyscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",3,false,false)));
  locsyscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("USEUPDATEDNODEPOS")));
  locsyscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("useupdatednodepos", 1)));

  Teuchos::RCP<ConditionDefinition> pointlocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT LOCSYS CONDITIONS",
                                         "Locsys",
                                         "Point local coordinate system",
                                         DRT::Condition::PointLocsys,
                                         true,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linelocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE LOCSYS CONDITIONS",
                                         "Locsys",
                                         "Line local coordinate system",
                                         DRT::Condition::LineLocsys,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surflocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF LOCSYS CONDITIONS",
                                         "Locsys",
                                         "Surface local coordinate system",
                                         DRT::Condition::SurfaceLocsys,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> vollocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL LOCSYS CONDITIONS",
                                         "Locsys",
                                         "Volume local coordinate system",
                                         DRT::Condition::VolumeLocsys,
                                         true,
                                         DRT::Condition::Volume));

  // Ale
  Teuchos::RCP<ConditionDefinition> pointalelocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT ALE LOCSYS CONDITIONS",
                                         "AleLocsys",
                                         "Point local coordinate system",
                                         DRT::Condition::PointLocsys,
                                         true,
                                         DRT::Condition::Point));
  Teuchos::RCP<ConditionDefinition> linealelocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE ALE LOCSYS CONDITIONS",
                                         "AleLocsys",
                                         "Line local coordinate system",
                                         DRT::Condition::LineLocsys,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfalelocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF ALE LOCSYS CONDITIONS",
                                         "AleLocsys",
                                         "Surface local coordinate system",
                                         DRT::Condition::SurfaceLocsys,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volalelocsys =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL ALE LOCSYS CONDITIONS",
                                         "AleLocsys",
                                         "Volume local coordinate system",
                                         DRT::Condition::VolumeLocsys,
                                         true,
                                         DRT::Condition::Volume));

  for (unsigned i=0; i<locsyscomponents.size(); ++i)
  {
    pointlocsys->AddComponent(locsyscomponents[i]);
    linelocsys->AddComponent(locsyscomponents[i]);
    surflocsys->AddComponent(locsyscomponents[i]);
    vollocsys->AddComponent(locsyscomponents[i]);

    //Ale
    pointalelocsys->AddComponent(locsyscomponents[i]);
    linealelocsys->AddComponent(locsyscomponents[i]);
    surfalelocsys->AddComponent(locsyscomponents[i]);
    volalelocsys->AddComponent(locsyscomponents[i]);
  }

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

  std::vector<Teuchos::RCP<ConditionComponent> > pbccomponents;

  pbccomponents.push_back(Teuchos::rcp(new IntConditionComponent("Id of periodic boundary condition",true)));
  pbccomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Is slave periodic boundary condition","Master",
        Teuchos::tuple<std::string>("Master","Slave"),
        Teuchos::tuple<std::string>("Master","Slave"))));
  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("PLANE")));
  pbccomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "degrees of freedom for the pbc plane","xy",
        Teuchos::tuple<std::string>("xy","yx","yz","zy","xz","zx","xyz"),
        Teuchos::tuple<std::string>("xy","xy","yz","yz","xz","xz","xyz"))));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("LAYER")));
  pbccomponents.push_back(Teuchos::rcp(new IntConditionComponent("Layer of periodic boundary condition",true)));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ANGLE")));
  pbccomponents.push_back(Teuchos::rcp(new RealConditionComponent("Angle of rotation")));

  pbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ABSTREETOL")));
  pbccomponents.push_back(Teuchos::rcp(new RealConditionComponent("Tolerance for nodematching in octree")));

  Teuchos::RCP<ConditionDefinition> lineperiodic =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE PERIODIC BOUNDARY CONDITIONS",
                                         "LinePeriodic",
                                         "Line Periodic",
                                         DRT::Condition::LinePeriodic,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfperiodic =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF PERIODIC BOUNDARY CONDITIONS",
                                         "SurfacePeriodic",
                                         "Surface Periodic",
                                         DRT::Condition::SurfacePeriodic,
                                         false,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<pbccomponents.size(); ++i)
  {
    lineperiodic->AddComponent(pbccomponents[i]);
    surfperiodic->AddComponent(pbccomponents[i]);
  }

  condlist.push_back(lineperiodic);
  condlist.push_back(surfperiodic);


  /*--------------------------------------------------------------------*/
  // transfer boundary condition for turbulent inflow

  std::vector<Teuchos::RCP<ConditionComponent> > tbc_turb_inflow_components;

  tbc_turb_inflow_components.push_back(
    Teuchos::rcp(new SeparatorConditionComponent("ID")));
  tbc_turb_inflow_components.push_back(
    Teuchos::rcp(new IntConditionComponent("id",true)));
  tbc_turb_inflow_components.push_back(
    Teuchos::rcp(new StringConditionComponent(
                   "toggle","master",
                   Teuchos::tuple<std::string>("master","slave"),
                   Teuchos::tuple<std::string>("master","slave"))));
  tbc_turb_inflow_components.push_back(
    Teuchos::rcp(new SeparatorConditionComponent("DIRECTION")));
  tbc_turb_inflow_components.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "transfer direction","x",
        Teuchos::tuple<std::string>("x","y","z"),
        Teuchos::tuple<std::string>("x","y","z"))));
  tbc_turb_inflow_components.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));

  Teuchos::RCP<ConditionDefinition> tbc_turb_inflow =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF TURBULENT INFLOW TRANSFER",
                                         "TransferTurbulentInflow",
                                         "TransferTurbulentInflow",
                                         DRT::Condition::TransferTurbulentInflow,
                                         true,
                                         DRT::Condition::Surface));

  // we attach all the components of this condition to this weak line DBC
  for (unsigned i=0; i<tbc_turb_inflow_components.size(); ++i)
  {
    tbc_turb_inflow->AddComponent(tbc_turb_inflow_components[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(tbc_turb_inflow);

  /*--------------------------------------------------------------------*/
  // separate domain for turbulent inflow generation

  Teuchos::RCP<ConditionDefinition> turbulentinflowgeneration =
    Teuchos::rcp(new ConditionDefinition("FLUID TURBULENT INFLOW VOLUME",
                                         "TurbulentInflowSection",
                                         "TurbulentInflowSection",
                                         DRT::Condition::TurbulentInflowSection,
                                         true,
                                         DRT::Condition::Volume));

   condlist.push_back(turbulentinflowgeneration);

  /*--------------------------------------------------------------------*/
  // volume condition to blend material parameters smoothly for turbulent combustion and two-phase flow

  std::vector<Teuchos::RCP<ConditionComponent> > blendmaterialcomponents;

  blendmaterialcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",1)));
  blendmaterialcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  blendmaterialcomponents.push_back(Teuchos::rcp(new StringConditionComponent("domain","plus",
              Teuchos::tuple<std::string>("plus","minus"),
              Teuchos::tuple<std::string>("plus","minus"))));

  Teuchos::RCP<ConditionDefinition> blendmaterial =
    Teuchos::rcp(new ConditionDefinition("COMBUST BLEND MATERIAL VOLUME",
                                         "BlendMaterial",
                                         "BlendMaterial",
                                         DRT::Condition::BlendMaterial,
                                         true,
                                         DRT::Condition::Volume));

  for (unsigned i=0; i<blendmaterialcomponents.size(); ++i)
    blendmaterial->AddComponent(blendmaterialcomponents[i]);

  condlist.push_back(blendmaterial);

   /*--------------------------------------------------------------------*/
  // flow-dependent pressure conditions

  std::vector<Teuchos::RCP<ConditionComponent> > flowdeppressurecomponents;

  // flow-dependent pressure conditions can be imposed either based on
  // (out)flow rate or (out)flow volume (e.g., for air-cushion condition)
  flowdeppressurecomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "type of flow dependence","flow_rate",
        Teuchos::tuple<std::string>("flow_rate","flow_volume","fixed_pressure"),
        Teuchos::tuple<std::string>("flow_rate","flow_volume","fixed_pressure"))));

  // constant coefficient for (linear) flow-rate-based condition
  // and constant fixed pressure
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("ConstCoeff")));

  // linear coefficient for (linear) flow-rate-based condition
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("LinCoeff")));

  // initial (air-cushion) volume outside of boundary
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("InitialVolume")));

  // reference pressure outside of boundary
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("ReferencePressure")));

  // adiabatic exponent
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("AdiabaticExponent")));

  // values for time curve
  flowdeppressurecomponents.push_back(Teuchos::rcp(new IntConditionComponent("curve",true,true)));


  Teuchos::RCP<ConditionDefinition> lineflowdeppressure
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE FLOW-DEPENDENT PRESSURE CONDITIONS",
                                         "LineFlowDepPressure",
                                         "LineFlowDepPressure",
                                         DRT::Condition::LineFlowDepPressure,
                                         true,
                                         DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfflowdeppressure
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE FLOW-DEPENDENT PRESSURE CONDITIONS",
                                         "SurfaceFlowDepPressure",
                                         "SurfaceFlowDepPressure",
                                         DRT::Condition::SurfaceFlowDepPressure,
                                         true,
                                         DRT::Condition::Surface));

  // we attach all the components of this condition to this weak line DBC
  for (unsigned i=0; i<flowdeppressurecomponents.size(); ++i)
  {
    lineflowdeppressure->AddComponent(flowdeppressurecomponents[i]);
    surfflowdeppressure->AddComponent(flowdeppressurecomponents[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(lineflowdeppressure);
  condlist.push_back(surfflowdeppressure);

  /*--------------------------------------------------------------------*/
  // BC-free boundary conditions

  std::vector<Teuchos::RCP<ConditionComponent> > bcfreecomponents;

  bcfreecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCTFORNODENORMAL")));
  bcfreecomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("functfornodenormal",3,false,false)));
  bcfreecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("USEUPDATEDNODEPOS")));
  bcfreecomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("useupdatednodepos", 1)));

  Teuchos::RCP<ConditionDefinition> linebcfree =
   Teuchos::rcp(new ConditionDefinition("DESIGN LINE BC-FREE BOUNDARY CONDITIONS",
                                        "LineBCFree",
                                        "LineBCFree",
                                        DRT::Condition::LineBCFree,
                                        true,
                                        DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfbcfree =
   Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE BC-FREE BOUNDARY CONDITIONS",
                                        "SurfaceBCFree",
                                        "SurfaceBCFree",
                                        DRT::Condition::SurfaceBCFree,
                                        true,
                                        DRT::Condition::Surface));

  for (unsigned i=0; i<bcfreecomponents.size(); ++i)
   {
     linebcfree->AddComponent(bcfreecomponents[i]);
     surfbcfree->AddComponent(bcfreecomponents[i]);
   }

  condlist.push_back(linebcfree);
  condlist.push_back(surfbcfree);

  /*--------------------------------------------------------------------*/
  // Navier-slip boundary conditions

  std::vector<Teuchos::RCP<ConditionComponent> > navierslipcomponents;

  navierslipcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT")));
  navierslipcomponents.push_back(Teuchos::rcp(new RealConditionComponent("slipcoefficient")));

  Teuchos::RCP<ConditionDefinition> linenavierslip =
   Teuchos::rcp(new ConditionDefinition("DESIGN LINE NAVIER-SLIP BOUNDARY CONDITIONS",
                                        "LineNavierSlip",
                                        "LineNavierSlip",
                                        DRT::Condition::LineNavierSlip,
                                        true,
                                        DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfnavierslip =
   Teuchos::rcp(new ConditionDefinition("DESIGN SURF NAVIER-SLIP BOUNDARY CONDITIONS",
                                        "SurfNavierSlip",
                                        "SurfNavierSlip",
                                        DRT::Condition::SurfNavierSlip,
                                        true,
                                        DRT::Condition::Surface));

  for (unsigned i=0; i<navierslipcomponents.size(); ++i)
   {
     linenavierslip->AddComponent(navierslipcomponents[i]);
     surfnavierslip->AddComponent(navierslipcomponents[i]);
   }

  condlist.push_back(linenavierslip);
  condlist.push_back(surfnavierslip);

   /*--------------------------------------------------------------------*/
  // weak Dirichlet conditions

  std::vector<Teuchos::RCP<ConditionComponent> > weakDirichletcomponents;

  // weak DBCs can be imposed adjoint consistent or adjoint inconsistent
  weakDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Choice of gamma parameter","adjoint-consistent",
        Teuchos::tuple<std::string>("adjoint-consistent","diffusive-optimal"),
        Teuchos::tuple<std::string>("adjoint-consistent","diffusive-optimal"))));

  // weak DBCs can be imposed in all directions or only in normal direction
  // (SCATRA: not checked, only in all_directions so far)
  weakDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Directions to apply weak dbc","all_directions",
        Teuchos::tuple<std::string>("all_directions","only_in_normal_direction"),
        Teuchos::tuple<std::string>("all_directions","only_in_normal_direction"))));

  // FLUID: penalty parameter either computed dynamically (using Spaldings law of
  // the wall) or by a fixed value; SCATRA: not checked, only constant value so far
  weakDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Definition of penalty parameter","constant",
        Teuchos::tuple<std::string>("constant","Spalding"),
        Teuchos::tuple<std::string>("constant","Spalding"))));

  // scaling factor for penalty parameter tauB or
  // stabilization parameter alpha for Nitsche term
  // (SCATRA: if stabilization parameter negative -> mixed-hybrid formulation)
  weakDirichletcomponents.push_back(Teuchos::rcp(new RealConditionComponent("TauBscaling")));

  // linearisation strategies --- the linearisation (i.e. the matrix
  // contribution) of the convective term on the inflow could be
  // suppressed, since the flux is a kink function and including this one
  // might result in even worse convergence behaviour
  // (SCATRA: not checked)
  weakDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Linearisation","lin_all",
        Teuchos::tuple<std::string>("lin_all","no_lin_conv_inflow"),
        Teuchos::tuple<std::string>("lin_all","no_lin_conv_inflow"))));

  // we provide a vector of 3 values for velocities
  weakDirichletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",3)));

  // values for curves
  weakDirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",3,true,true)));

  // and optional spatial functions
  weakDirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",3,false,false,true)));


  Teuchos::RCP<ConditionDefinition> lineweakdirichlet
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE WEAK DIRICHLET CONDITIONS",
                                         "LineWeakDirichlet",
                                         "LineWeakDirichlet",
                                         DRT::Condition::LineWeakDirichlet,
                                         true,
                                         DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfweakdirichlet
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE WEAK DIRICHLET CONDITIONS",
                                         "SurfaceWeakDirichlet",
                                         "SurfaceWeakDirichlet",
                                         DRT::Condition::SurfaceWeakDirichlet,
                                         true,
                                         DRT::Condition::Surface));

  // we attach all the components of this condition to this weak line DBC
  for (unsigned i=0; i<weakDirichletcomponents.size(); ++i)
  {
    lineweakdirichlet->AddComponent(weakDirichletcomponents[i]);
    surfweakdirichlet->AddComponent(weakDirichletcomponents[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(lineweakdirichlet);
  condlist.push_back(surfweakdirichlet);

  /*--------------------------------------------------------------------*/
  // mixed/hybrid Dirichlet conditions

  std::vector<Teuchos::RCP<ConditionComponent> > mixhybDirichletcomponents;

  // we provide a vector of 3 values for velocities
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",3)));

  // values for curves
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",3,true,true)));

  // and optional spatial functions
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",3,false,false,true)));

  // characteristic velocity
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new RealConditionComponent("u_C")));

  // the penalty parameter could be computed dynamically (using Spaldings
  // law of the wall) or using a fixed value (1)
  mixhybDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Definition of penalty parameter","constant",
        Teuchos::tuple<std::string>("constant","Spalding"),
        Teuchos::tuple<std::string>("constant","Spalding"))));

  // scaling factor for penalty parameter tauB
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new RealConditionComponent("hB_divided_by")));

  // if Spaldings law is used, this defines the way how the traction at y is computed from utau
  mixhybDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "utau_computation","at_wall",
        Teuchos::tuple<std::string>("at_wall","viscous_tangent"),
        Teuchos::tuple<std::string>("at_wall","viscous_tangent"))));


  Teuchos::RCP<ConditionDefinition> linemixhybDirichlet
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE MIXED/HYBRID DIRICHLET CONDITIONS",
                                         "LineMixHybDirichlet",
                                         "LineMixHybDirichlet",
                                         DRT::Condition::LineMixHybDirichlet,
                                         true,
                                         DRT::Condition::Line));


  Teuchos::RCP<ConditionDefinition> surfmixhybDirichlet
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE MIXED/HYBRID DIRICHLET CONDITIONS",
                                         "SurfaceMixHybDirichlet",
                                         "SurfaceMixHybDirichlet",
                                         DRT::Condition::SurfaceMixHybDirichlet,
                                         true,
                                         DRT::Condition::Surface));

  // we attach all the components of this condition to this condition
  for (unsigned i=0; i<mixhybDirichletcomponents.size(); ++i)
  {
    linemixhybDirichlet->AddComponent(mixhybDirichletcomponents[i]);
    surfmixhybDirichlet->AddComponent(mixhybDirichletcomponents[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(linemixhybDirichlet);
  condlist.push_back(surfmixhybDirichlet);

  /*--------------------------------------------------------------------*/
  // consistent outflow bcs for conservative element formulations

  Teuchos::RCP<ConditionDefinition> surfconsistentoutflowconsistency
    =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE CONSERVATIVE OUTFLOW CONSISTENCY",
                                         "SurfaceConservativeOutflowConsistency",
                                         "SurfaceConservativeOutflowConsistency",
                                         DRT::Condition::SurfaceConservativeOutflowConsistency,
                                         true,
                                         DRT::Condition::Surface));

   condlist.push_back(surfconsistentoutflowconsistency);

  /*--------------------------------------------------------------------*/
  // Neumann inflow for FLUID

  Teuchos::RCP<ConditionDefinition> linefluidneumanninflow =
    Teuchos::rcp(new ConditionDefinition("FLUID NEUMANN INFLOW LINE CONDITIONS",
                                         "FluidNeumannInflow",
                                         "Line Fluid Neumann Inflow",
                                         DRT::Condition::FluidNeumannInflow,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffluidneumanninflow =
    Teuchos::rcp(new ConditionDefinition("FLUID NEUMANN INFLOW SURF CONDITIONS",
                                         "FluidNeumannInflow",
                                         "Surface Fluid Neumann Inflow",
                                         DRT::Condition::FluidNeumannInflow,
                                         true,
                                         DRT::Condition::Surface));

   condlist.push_back(linefluidneumanninflow);
   condlist.push_back(surffluidneumanninflow);

  /*--------------------------------------------------------------------*/
  // Neumann inflow for SCATRA

  Teuchos::RCP<ConditionDefinition> linetransportneumanninflow =
    Teuchos::rcp(new ConditionDefinition("TRANSPORT NEUMANN INFLOW LINE CONDITIONS",
                                         "TransportNeumannInflow",
                                         "Line Transport Neumann Inflow",
                                         DRT::Condition::TransportNeumannInflow,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportneumanninflow =
    Teuchos::rcp(new ConditionDefinition("TRANSPORT NEUMANN INFLOW SURF CONDITIONS",
                                         "TransportNeumannInflow",
                                         "Surface Transport Neumann Inflow",
                                         DRT::Condition::TransportNeumannInflow,
                                         true,
                                         DRT::Condition::Surface));

   condlist.push_back(linetransportneumanninflow);
   condlist.push_back(surftransportneumanninflow);

   /*--------------------------------------------------------------------*/
   // Taylor Galerkin outflow Boundaries for level set transport equation

   Teuchos::RCP<ConditionDefinition> surfOutflowTaylorGalerkin =
     Teuchos::rcp(new ConditionDefinition("TAYLOR GALERKIN OUTFLOW SURF CONDITIONS",
                                          "TaylorGalerkinOutflow",
                                          "Surface Taylor Galerkin Outflow",
                                          DRT::Condition::TaylorGalerkinOutflow,
                                          true,
                                          DRT::Condition::Surface));

    condlist.push_back(surfOutflowTaylorGalerkin);

  /*--------------------------------------------------------------------*/

    Teuchos::RCP<ConditionDefinition> surfneumanninflowTaylorGalerkin =
      Teuchos::rcp(new ConditionDefinition("TAYLOR GALERKIN NEUMANN INFLOW SURF CONDITIONS",
                                           "TaylorGalerkinNeumannInflow",
                                           "Surface Taylor Galerkin Neumann Inflow",
                                           DRT::Condition::TaylorGalerkinNeumannInflow,
                                           true,
                                           DRT::Condition::Surface));

     condlist.push_back(surfneumanninflowTaylorGalerkin);


  /*--------------------------------------------------------------------*/
  // Characteristic Galerkin Boundaries for LevelSet-Reinitialization

  Teuchos::RCP<ConditionDefinition> surfreinitializationtaylorgalerkin =
    Teuchos::rcp(new ConditionDefinition("REINITIALIZATION TAYLOR GALERKIN SURF CONDITIONS",
                                         "ReinitializationTaylorGalerkin",
                                         "Surface Reinitialization Taylor Galerkin",
                                         DRT::Condition::ReinitializationTaylorGalerkin,
                                         true,
                                         DRT::Condition::Surface));

   condlist.push_back(surfreinitializationtaylorgalerkin);

  // FSI

  std::vector<Teuchos::RCP<ConditionComponent> > fsicomponents;

  fsicomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linefsi =
    Teuchos::rcp(new ConditionDefinition("DESIGN FSI COUPLING LINE CONDITIONS",
                                         "FSICoupling",
                                         "FSI Coupling",
                                         DRT::Condition::FSICoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffsi =
    Teuchos::rcp(new ConditionDefinition("DESIGN FSI COUPLING SURF CONDITIONS",
                                         "FSICoupling",
                                         "FSI Coupling",
                                         DRT::Condition::FSICoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<fsicomponents.size(); ++i)
  {
    linefsi->AddComponent(fsicomponents[i]);
    surffsi->AddComponent(fsicomponents[i]);
  }

  condlist.push_back(linefsi);
  condlist.push_back(surffsi);

  /*--------------------------------------------------------------------*/
  // FSI without sliding

  Teuchos::RCP<ConditionDefinition> linefsins =
    Teuchos::rcp(new ConditionDefinition("DESIGN FSI COUPLING NO SLIDE LINE CONDITIONS",
                                         "FSICouplingNoSlide",
                                         "FSI Coupling No Slide",
                                         DRT::Condition::FSICouplingNoSlide,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffsins =
    Teuchos::rcp(new ConditionDefinition("DESIGN FSI COUPLING NO SLIDE SURF CONDITIONS",
                                         "FSICouplingNoSlide",
                                         "FSI Coupling No Slide",
                                         DRT::Condition::FSICouplingNoSlide,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(linefsins);
  condlist.push_back(surffsins);

  /*--------------------------------------------------------------------*/
  // FSI define centerdisp for sliding interfaces

  Teuchos::RCP<ConditionDefinition> linefsicd =
    Teuchos::rcp(new ConditionDefinition("DESIGN FSI COUPLING CENTER DISP LINE CONDITIONS",
                                         "FSICouplingCenterDisp",
                                         "FSI Coupling Center Disp",
                                         DRT::Condition::FSICouplingCenterDisp,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffsicd =
    Teuchos::rcp(new ConditionDefinition("DESIGN FSI COUPLING CENTER DISP SURF CONDITIONS",
                                         "FSICouplingCenterDisp",
                                         "FSI Coupling Center Disp",
                                         DRT::Condition::FSICouplingCenterDisp,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(linefsicd);
  condlist.push_back(surffsicd);

  /*--------------------------------------------------------------------*/
  // FPSI

  std::vector<Teuchos::RCP<ConditionComponent> > fpsicomponents;

  fpsicomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linefpsi =
    Teuchos::rcp(new ConditionDefinition("DESIGN FPSI COUPLING LINE CONDITIONS",
                                         "FPSICoupling",
                                         "FPSI Coupling",
                                         DRT::Condition::FPSICoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffpsi =
    Teuchos::rcp(new ConditionDefinition("DESIGN FPSI COUPLING SURF CONDITIONS",
                                         "FPSICoupling",
                                         "FPSI Coupling",
                                         DRT::Condition::FPSICoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<fsicomponents.size(); ++i)
  {
    linefpsi->AddComponent(fpsicomponents[i]);
    surffpsi->AddComponent(fpsicomponents[i]);
  }

  condlist.push_back(linefpsi);
  condlist.push_back(surffpsi);

/*--------------------------------------------------------------------*/
  // IMMERSED FSI

  Teuchos::RCP<ConditionDefinition> immersedsearchbox =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOLUME IMMERSED SEARCHBOX",
                                         "ImmersedSearchbox",
                                         "Immersed Searchbox",
                                         DRT::Condition::ImmersedSearchbox,
                                         true,
                                         DRT::Condition::Volume));

  condlist.push_back(immersedsearchbox);

  /*--------------------------------------------------------------------*/
  // FREESURF

  std::vector<Teuchos::RCP<ConditionComponent> > freesurfcomponents;

  freesurfcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "field","fluid",
        Teuchos::tuple<std::string>("fluid","ale"),
        Teuchos::tuple<std::string>("fluid","ale"))));

  freesurfcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "coupling","lagrange",
        Teuchos::tuple<std::string>("lagrange","heightfunction"),
        Teuchos::tuple<std::string>("lagrange","heightfunction"),
        true)));

  Teuchos::RCP<ConditionDefinition> linefreesurf =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUID FREE SURFACE LINE CONDITIONS",
                                         "FREESURFCoupling",
                                         "FREESURF Coupling",
                                         DRT::Condition::FREESURFCoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffreesurf =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUID FREE SURFACE SURF CONDITIONS",
                                         "FREESURFCoupling",
                                         "FREESURF Coupling",
                                         DRT::Condition::FREESURFCoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<freesurfcomponents.size(); ++i)
  {
    linefreesurf->AddComponent(freesurfcomponents[i]);
    surffreesurf->AddComponent(freesurfcomponents[i]);
  }

  condlist.push_back(linefreesurf);
  condlist.push_back(surffreesurf);

  /*--------------------------------------------------------------------*/
  // Local Lagrange boundary condition

  Teuchos::RCP<ConditionDefinition> linelocallagrange =
    Teuchos::rcp(new ConditionDefinition("DESIGN LOCAL LAGRANGE LINE CONDITIONS",
                                         "LOCALLAGRANGECoupling",
                                         "LOCALLAGRANGE Coupling",
                                         DRT::Condition::LOCALLAGRANGECoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surflocallagrange =
    Teuchos::rcp(new ConditionDefinition("DESIGN LOCAL LAGRANGE SURF CONDITIONS",
                                         "LOCALLAGRANGECoupling",
                                         "LOCALLAGRANGE Coupling",
                                         DRT::Condition::LOCALLAGRANGECoupling,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(linelocallagrange);
  condlist.push_back(surflocallagrange);

  /*--------------------------------------------------------------------*/
  // Additional coupling of structure and ale fields (for lung fsi)

  std::vector<Teuchos::RCP<ConditionComponent> > saccomponents;

  saccomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  saccomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "field","structure",
        Teuchos::tuple<std::string>("structure","fluid"),
        Teuchos::tuple<std::string>("structure","fluid"))));

  Teuchos::RCP<ConditionDefinition> surfsac =
    Teuchos::rcp(new ConditionDefinition("DESIGN STRUCTURE ALE COUPLING SURF CONDITIONS",
                                         "StructAleCoupling",
                                         "StructAleCoupling",
                                         DRT::Condition::StructAleCoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<saccomponents.size(); ++i)
    surfsac->AddComponent(saccomponents[i]);

  condlist.push_back(surfsac);

  /*--------------------------------------------------------------------*/
  // Additional coupling of structure and fluid volumes (for lung fsi)

  std::vector<Teuchos::RCP<ConditionComponent> > sfvcomponents;

  sfvcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  sfvcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "field","structure",
        Teuchos::tuple<std::string>("structure","fluid"),
        Teuchos::tuple<std::string>("structure","fluid"))));

  Teuchos::RCP<ConditionDefinition> surfsfv =
    Teuchos::rcp(new ConditionDefinition("DESIGN STRUCTURE FLUID VOLUME COUPLING SURF CONDITIONS",
                                         "StructFluidSurfCoupling",
                                         "StructFluidSurfCoupling",
                                         DRT::Condition::StructFluidSurfCoupling,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volsfv =
    Teuchos::rcp(new ConditionDefinition("DESIGN STRUCTURE FLUID VOLUME COUPLING VOL CONDITIONS",
                                         "StructFluidVolCoupling",
                                         "StructFluidVolCoupling",
                                         DRT::Condition::StructFluidVolCoupling,
                                         false,
                                         DRT::Condition::Volume));

  for (unsigned i=0; i<sfvcomponents.size(); ++i)
  {
    surfsfv->AddComponent(sfvcomponents[i]);
    volsfv->AddComponent(sfvcomponents[i]);
  }

  condlist.push_back(surfsfv);
  condlist.push_back(volsfv);

  /*--------------------------------------------------------------------*/
  // Additional coupling for biofilm growth

  std::vector<Teuchos::RCP<ConditionComponent> > biogrcomponents;

  biogrcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linebiogr =
    Teuchos::rcp(new ConditionDefinition("DESIGN BIOFILM GROWTH COUPLING LINE CONDITIONS",
                                         "BioGrCoupling",
                                         "BioGrCoupling",
                                         DRT::Condition::BioGrCoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfbiogr =
    Teuchos::rcp(new ConditionDefinition("DESIGN BIOFILM GROWTH COUPLING SURF CONDITIONS",
                                         "BioGrCoupling",
                                         "BioGrCoupling",
                                         DRT::Condition::BioGrCoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<biogrcomponents.size(); ++i)
  {
    linebiogr->AddComponent(biogrcomponents[i]);
    surfbiogr->AddComponent(biogrcomponents[i]);
  }

  condlist.push_back(linebiogr);
  condlist.push_back(surfbiogr);

  /*--------------------------------------------------------------------*/
  // xfem

  std::vector<Teuchos::RCP<ConditionComponent> > xfemcomponents;

  xfemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("label")));

  Teuchos::RCP<ConditionDefinition> linexfem =
    Teuchos::rcp(new ConditionDefinition("DESIGN XFEM COUPLING LINE CONDITIONS",
                                         "XFEMCoupling",
                                         "XFEM Coupling",
                                         DRT::Condition::XFEMCoupling,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfxfem =
    Teuchos::rcp(new ConditionDefinition("DESIGN XFEM COUPLING SURF CONDITIONS",
                                         "XFEMCoupling",
                                         "XFEM Coupling",
                                         DRT::Condition::XFEMCoupling,
                                         true,
                                         DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> movingfluid =
        Teuchos::rcp(new ConditionDefinition("DESIGN MOVING FLUID VOL CONDITIONS",
                                             "MovingFluid",
                                             "Moving Fluid",
                                             DRT::Condition::MovingFluid,
                                             true,
                                             DRT::Condition::Volume));
  Teuchos::RCP<ConditionDefinition> fluidfluidcoupling =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID FLUID COUPLING SURF CONDITIONS",
                                           "FluidFluidCoupling",
                                           "FLUID FLUID Coupling",
                                           DRT::Condition::FluidFluidCoupling,
                                           true,
                                           DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> ALEfluidcoupling =
      Teuchos::rcp(new ConditionDefinition("DESIGN ALE FLUID COUPLING SURF CONDITIONS",
                                           "ALEFluidCoupling",
                                           "ALE FLUID Coupling",
                                           DRT::Condition::ALEFluidCoupling,
                                           true,
                                           DRT::Condition::Surface));

  for (unsigned i=0; i<xfemcomponents.size(); ++i)
  {
    linexfem->AddComponent(xfemcomponents[i]);
    surfxfem->AddComponent(xfemcomponents[i]);
    movingfluid->AddComponent(xfemcomponents[i]);
    fluidfluidcoupling->AddComponent(xfemcomponents[i]);
    ALEfluidcoupling->AddComponent(xfemcomponents[i]);
  }

  condlist.push_back(linexfem);
  condlist.push_back(surfxfem);
  condlist.push_back(fluidfluidcoupling);
  condlist.push_back(movingfluid);
  condlist.push_back(ALEfluidcoupling);

  /*--------------------------------------------------------------------*/
  // crack

  std::vector<Teuchos::RCP<ConditionComponent> > crackdef;
  Teuchos::RCP<ConditionDefinition> crmaster =
      Teuchos::rcp(new ConditionDefinition("DESIGN CRACK MASTER SURFACE",
                                           "masterCrackSurface",
                                           "master crack surface",
                                           DRT::Condition::CrackMastersurface,
                                           true,
                                           DRT::Condition::Surface));

  Teuchos::RCP<ConditionDefinition> crslave =
        Teuchos::rcp(new ConditionDefinition("DESIGN CRACK SLAVE SURFACE",
                                             "slaveCrackSurface",
                                             "slave crack surface",
                                             DRT::Condition::CrackSlavesurface,
                                             true,
                                             DRT::Condition::Surface));

  Teuchos::RCP<ConditionDefinition> initialCrack =
          Teuchos::rcp(new ConditionDefinition("DESIGN CRACK INITIATION POINTS",
                                               "CrackInitiationPoints",
                                               "Points at which crack initiates",
                                               DRT::Condition::CrackInitPoints,
                                               true,
                                               DRT::Condition::Point));

  Teuchos::RCP<ConditionDefinition> crackBoundary =
          Teuchos::rcp(new ConditionDefinition("DESIGN CRACK BOUNDARY POINTS",
                                               "CrackBoundaryPoints",
                                               "All nodes on the boundary of the domain",
                                               DRT::Condition::CrackBoundaryPoints,
                                               true,
                                               DRT::Condition::Point));

  for (unsigned i=0; i<crackdef.size(); ++i)
  {
    crmaster->AddComponent(crackdef[i]);
    crslave->AddComponent(crackdef[i]);
    initialCrack->AddComponent(crackdef[i]);
    crackBoundary->AddComponent(crackdef[i]);
  }

  condlist.push_back(crmaster);
  condlist.push_back(crslave);
  condlist.push_back(initialCrack);
  condlist.push_back(crackBoundary);

  /*--------------------------------------------------------------------*/

  /*--------------------------------------------------------------------*/
  // surface tension

  Teuchos::RCP<ConditionDefinition> surftension =
    Teuchos::rcp(new ConditionDefinition("SURFACE TENSION CONDITIONS",
                                         "SurfaceStress",
                                         "Surface Stress (ideal water)",
                                         DRT::Condition::SurfaceTension,
                                         true,
                                         DRT::Condition::Surface));

  surftension->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedReal(surftension,"gamma");

  condlist.push_back(surftension);

  /*--------------------------------------------------------------------*/
  // surfactant

  Teuchos::RCP<ConditionDefinition> surfactant =
    Teuchos::rcp(new ConditionDefinition("SURFACTANT CONDITIONS",
                                         "SurfaceStress",
                                         "Surface Stress (surfactant)",
                                         DRT::Condition::Surfactant,
                                         true,
                                         DRT::Condition::Surface));

  surfactant->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedReal(surfactant,"k1xCbulk");
  AddNamedReal(surfactant,"k2");
  AddNamedReal(surfactant,"m1");
  AddNamedReal(surfactant,"m2");
  AddNamedReal(surfactant,"gamma_0");
  AddNamedReal(surfactant,"gamma_min");

  condlist.push_back(surfactant);

  /*--------------------------------------------------------------------*/
  // Lennard Jones potential volume
  Teuchos::RCP<ConditionDefinition> lj_potential_volume =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential_Volume",
                                         DRT::Condition::LJ_Potential_Volume,
                                         true,
                                         DRT::Condition::Volume));

  lj_potential_volume->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential_volume,"label");
  AddNamedReal(lj_potential_volume,"depth");
  AddNamedReal(lj_potential_volume,"rootDist");
  AddNamedReal(lj_potential_volume,"cutOff");
  AddNamedReal(lj_potential_volume,"exvollength");
  AddNamedReal(lj_potential_volume,"beta");

  condlist.push_back(lj_potential_volume);

  /*--------------------------------------------------------------------*/
  // Lennard Jones potential surface

  Teuchos::RCP<ConditionDefinition> lj_potential_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential_Surface",
                                         DRT::Condition::LJ_Potential_Surface,
                                         true,
                                         DRT::Condition::Surface));

  lj_potential_surface->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential_surface,"label");
  AddNamedReal(lj_potential_surface,"depth");
  AddNamedReal(lj_potential_surface,"rootDist");
  AddNamedReal(lj_potential_surface,"cutOff");
  AddNamedReal(lj_potential_surface,"exvollength");
  AddNamedReal(lj_potential_surface,"beta");

  condlist.push_back(lj_potential_surface);


  /*--------------------------------------------------------------------*/
  // Lennard Jones potential line

  Teuchos::RCP<ConditionDefinition> lj_potential_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential_Line",
                                         DRT::Condition::LJ_Potential_Line,
                                         true,
                                         DRT::Condition::Line));

  lj_potential_line->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential_line,"label");
  AddNamedReal(lj_potential_line,"depth");
  AddNamedReal(lj_potential_line,"rootDist");
  AddNamedReal(lj_potential_line,"cutOff");
  AddNamedReal(lj_potential_line,"exvollength");
  AddNamedReal(lj_potential_line,"beta");

  condlist.push_back(lj_potential_line);


  /*--------------------------------------------------------------------*/
  // Van der Waals potential volume
  Teuchos::RCP<ConditionDefinition> vanderwaals_potential_volume =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL VAN DER WAALS POTENTIAL CONDITIONS",
                                         "Potential",
                                         "VanDerWaals_Potential_Volume",
                                         DRT::Condition::VanDerWaals_Potential_Volume,
                                         true,
                                         DRT::Condition::Volume));

  vanderwaals_potential_volume->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(vanderwaals_potential_volume,"label");
  AddNamedReal(vanderwaals_potential_volume,"lambda");
  AddNamedReal(vanderwaals_potential_volume,"cutOff");
  AddNamedReal(vanderwaals_potential_volume,"beta");
  AddNamedReal(vanderwaals_potential_volume,"exvollength");

  condlist.push_back(vanderwaals_potential_volume);

  /*--------------------------------------------------------------------*/
  // Van der Waals potential surface

  Teuchos::RCP<ConditionDefinition> vanderwaals_potential_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF VAN DER WAALS POTENTIAL CONDITIONS",
                                         "Potential",
                                         "VanDerWaals_Potential_Surface",
                                         DRT::Condition::VanDerWaals_Potential_Surface,
                                         true,
                                         DRT::Condition::Surface));

  vanderwaals_potential_surface->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(vanderwaals_potential_surface,"label");
  AddNamedReal(vanderwaals_potential_surface,"lambda");
  AddNamedReal(vanderwaals_potential_surface,"cutOff");
  AddNamedReal(vanderwaals_potential_surface,"beta");
  AddNamedReal(vanderwaals_potential_surface,"exvollength");

  condlist.push_back(vanderwaals_potential_surface);


  /*--------------------------------------------------------------------*/
  // Van der Waals line

  Teuchos::RCP<ConditionDefinition> vanderwaals_potential_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE VAN DER WAALS POTENTIAL CONDITIONS",
                                         "Potential",
                                         "VanDerWaals_Potential_Line",
                                         DRT::Condition::VanDerWaals_Potential_Line,
                                         true,
                                         DRT::Condition::Line));

  vanderwaals_potential_line->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(vanderwaals_potential_line,"label");
  AddNamedReal(vanderwaals_potential_line,"lambda");
  AddNamedReal(vanderwaals_potential_line,"cutOff");
  AddNamedReal(vanderwaals_potential_line,"beta");
  AddNamedReal(vanderwaals_potential_line,"exvollength");

  condlist.push_back(vanderwaals_potential_line);


  /*-------------------------------------------------------------------*/
  // Electrostatic Repulsion Surface
  Teuchos::RCP<ConditionDefinition> electro_repulsion_potential_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF ELECTRO REPULSION CONDITIONS",
                                       "Potential",
                                       "ElectroRepulsion_Potential_Surface",
                                       DRT::Condition::ElectroRepulsion_Potential_Surface,
                                       true,
                                       DRT::Condition::Surface));

  electro_repulsion_potential_surface->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(electro_repulsion_potential_surface,"label");
  AddNamedReal(electro_repulsion_potential_surface,"zeta_param_1");
  AddNamedReal(electro_repulsion_potential_surface,"zeta_param_2");
  AddNamedReal(electro_repulsion_potential_surface,"cutOff");
  AddNamedReal(electro_repulsion_potential_surface,"beta");

  condlist.push_back(electro_repulsion_potential_surface);


  /*-------------------------------------------------------------------*/
  // Electrostatic Repulsion Line
  Teuchos::RCP<ConditionDefinition> electro_repulsion_potential_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE ELECTRO REPULSION CONDITIONS",
                                       "Potential",
                                       "ElectroRepulsion_Potential_Line",
                                       DRT::Condition::ElectroRepulsion_Potential_Line,
                                       true,
                                       DRT::Condition::Line));

  electro_repulsion_potential_line->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(electro_repulsion_potential_line,"label");
  AddNamedReal(electro_repulsion_potential_line,"zeta_param_1");
  AddNamedReal(electro_repulsion_potential_line,"zeta_param_2");
  AddNamedReal(electro_repulsion_potential_line,"cutOff");
  AddNamedReal(electro_repulsion_potential_line,"beta");

  condlist.push_back(electro_repulsion_potential_line);


  /*--------------------------------------------------------------------*/
  // Fluctuating Hydrodynamics Statistics on a surface

  std::vector<Teuchos::RCP<ConditionComponent> > flucthydrostatsurfcomponents;
  flucthydrostatsurfcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  flucthydrostatsurfcomponents.push_back(
      Teuchos::rcp(
           new StringConditionComponent(
             "evaluation type",
             "nodalbased",
             Teuchos::tuple<std::string>("elebased","nodalbased", "ele_and_nodalbased"),
             Teuchos::tuple<std::string>("elebased","nodalbased", "ele_and_nodalbased"))));


  Teuchos::RCP<ConditionDefinition> fluctHydro_statisticsSurf =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUCTHYDRO STATISTICS SURF CONDITIONS",
                                         "FluctHydroStatisticsSurf",
                                         "FluctHydro_StatisticsSurf",
                                         DRT::Condition::FluctHydro_StatisticsSurf,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<flucthydrostatsurfcomponents.size(); ++i)
    fluctHydro_statisticsSurf->AddComponent(flucthydrostatsurfcomponents[i]);

  condlist.push_back(fluctHydro_statisticsSurf);

  /*--------------------------------------------------------------------*/
  // Fluctuating Hydrodynamics Statistics on a line

  std::vector<Teuchos::RCP<ConditionComponent> > flucthydrostatlinecomponents;
  flucthydrostatlinecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  flucthydrostatlinecomponents.push_back(
        Teuchos::rcp(
          new StringConditionComponent(
            "evaluation type",
            "nodalbased",
            Teuchos::tuple<std::string>("elebased","nodalbased", "ele_and_nodalbased"),
            Teuchos::tuple<std::string>("elebased","nodalbased", "ele_and_nodalbased"))));

  Teuchos::RCP<ConditionDefinition> fluctHydro_statisticsLine =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUCTHYDRO STATISTICS LINE CONDITIONS",
                                         "FluctHydroStatisticsLine",
                                         "FluctHydro_StatisticsLine",
                                         DRT::Condition::FluctHydro_StatisticsLine,
                                         true,
                                         DRT::Condition::Line));

  for (unsigned i=0; i<flucthydrostatlinecomponents.size(); ++i)
      fluctHydro_statisticsLine->AddComponent(flucthydrostatlinecomponents[i]);

  condlist.push_back(fluctHydro_statisticsLine);

  /*--------------------------------------------------------------------*/
  // Brownian Motion

  Teuchos::RCP<ConditionDefinition> brownian_motion =
    Teuchos::rcp(new ConditionDefinition("DESIGN BROWNIAN MOTION SURF CONDITIONS",
                                         "BrownianMotion",
                                         "Brownian_Motion",
                                         DRT::Condition::Brownian_Motion,
                                         true,
                                         DRT::Condition::Surface));

  brownian_motion->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(brownian_motion,"label");
  AddNamedReal(brownian_motion,"boltz_const");
  AddNamedReal(brownian_motion,"temperatur");
  AddNamedReal(brownian_motion,"frict_coeff");

  condlist.push_back(brownian_motion);

  /*--------------------------------------------------------------------*/
  // Filament Number

  //declaration of a variable which contains all the components of the condition
  std::vector<Teuchos::RCP<ConditionComponent> > filamentnumbercomponents;

  //the condition consists of one component which has to read an integer value (the so called filament number):
  filamentnumbercomponents.push_back(Teuchos::rcp(new IntConditionComponent("Filament Number")));

  //the condition itself hast to be defined so that it is clear how to read or write such a condition in the dat file
  Teuchos::RCP<ConditionDefinition> filamentnumber =
    Teuchos::rcp(new ConditionDefinition("FILAMENT NUMBERS",              //name of input file section
                                         "FilamentNumber",                //name to get the condition from a discretization
                                         "Filament Number",               //description of condition
                                         DRT::Condition::FilamentNumber,  //type of condition in DRT (cf. drt_condition.H)
                                         false,                           //should explicit elements be generated (e.g. line elements for line contact condition in 2D)?
                                         DRT::Condition::Line));          //type of geometry the condition lives on

  //after definition of the condition all its components have to be added:
  for (unsigned i =0 ; i < filamentnumbercomponents.size(); ++i)
    filamentnumber->AddComponent(filamentnumbercomponents[i]);

  //the condition itself has to be added to the condition list
  condlist.push_back(filamentnumber);

  /*--------------------------------------------------------------------*/
  // Force Sensor

  //declaration of a variable which contains all the components of the condition
  std::vector<Teuchos::RCP<ConditionComponent> > forcesensorcomponents;

  /*the condition consists of one component which has to read an integer value (the number of the dof of the node with respect
   * to which the force is to be measured; note: numbering starts at zero):*/
  forcesensorcomponents.push_back(Teuchos::rcp(new IntConditionComponent("DOF Number")));

  //the condition itself hast to be defined so that it is clear how to read or write such a condition in the dat file
  Teuchos::RCP<ConditionDefinition> forcesensor =
    Teuchos::rcp(new ConditionDefinition("FORCE SENSORS",                 //name of input file section
                                         "ForceSensor",                   //name to get the condition from a discretization
                                         "Force Sensor",                  //description of condition
                                         DRT::Condition::ForceSensor   ,  //type of condition in DRT (cf. drt_condition.H)
                                         false,                           //should explicit elements be generated (e.g. line elements for line contact condition in 2D)?
                                         DRT::Condition::Point));         //type of geometry the condition lives on

  //after definition of the condition all its components have to be added:
  for (unsigned i =0 ; i < forcesensorcomponents.size(); ++i)
    forcesensor->AddComponent(forcesensorcomponents[i]);

  //the condition itself has to be added to the condition list
  condlist.push_back(forcesensor);

  /*--------------------------------------------------------------------*/
  // microscale boundary

  Teuchos::RCP<ConditionDefinition> microscale =
    Teuchos::rcp(new ConditionDefinition("MICROSCALE CONDITIONS",
                                         "MicroBoundary",
                                         "Microscale Boundary",
                                         DRT::Condition::MicroBoundary,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(microscale);

  /*--------------------------------------------------------------------*/
  // fluid stress

  std::vector<Teuchos::RCP<ConditionComponent> > fluidstresscomponents;

  Teuchos::RCP<ConditionDefinition> linefluidstress =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUID STRESS CALC LINE CONDITIONS",
                                         "FluidStressCalc",
                                         "Line Fluid Stress Calculation",
                                         DRT::Condition::FluidStressCalc,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffluidstress =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUID STRESS CALC SURF CONDITIONS",
                                         "FluidStressCalc",
                                         "Surf Fluid Stress Calculation",
                                         DRT::Condition::FluidStressCalc,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<fluidstresscomponents.size(); ++i)
  {
    linefluidstress->AddComponent(fluidstresscomponents[i]);
    surffluidstress->AddComponent(fluidstresscomponents[i]);
  }

  condlist.push_back(linefluidstress);
  condlist.push_back(surffluidstress);

  /*--------------------------------------------------------------------*/
  // lift & drag

  std::vector<Teuchos::RCP<ConditionComponent> > liftdragcomponents;

  liftdragcomponents.push_back(Teuchos::rcp(new IntConditionComponent("label")));
  liftdragcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("CENTER")));
  liftdragcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("centerCoord",3)));
  // optional
  liftdragcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("AXIS",true)));
  liftdragcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("axis",3,true)));

  Teuchos::RCP<ConditionDefinition> lineliftdrag =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUID LINE LIFT&DRAG",
                                         "LIFTDRAG",
                                         "Line LIFTDRAG",
                                         DRT::Condition::LineLIFTDRAG,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfliftdrag =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUID SURF LIFT&DRAG",
                                         "LIFTDRAG",
                                         "Surface LIFTDRAG",
                                         DRT::Condition::SurfLIFTDRAG,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<liftdragcomponents.size(); ++i)
  {
    lineliftdrag->AddComponent(liftdragcomponents[i]);
    surfliftdrag->AddComponent(liftdragcomponents[i]);
  }

  condlist.push_back(lineliftdrag);
  condlist.push_back(surfliftdrag);

  /*--------------------------------------------------------------------*/
  // stc layer condition

  Teuchos::RCP<ConditionDefinition> stclayer =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL STC LAYER",
                                         "STC Layer",
                                         "Layer for Multilayered STC",
                                         DRT::Condition::VolSTCLayer,
                                         true,
                                         DRT::Condition::Volume));

  stclayer->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(stclayer);

  /*--------------------------------------------------------------------*/
  // volume constraint

  Teuchos::RCP<ConditionDefinition> volumeconstraint =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE VOLUME CONSTRAINT 3D",
                                         "VolumeConstraint_3D",
                                         "Surface Volume Constraint",
                                         DRT::Condition::VolumeConstraint_3D,
                                         true,
                                         DRT::Condition::Surface));

  volumeconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  volumeconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  volumeconstraint->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  volumeconstraint->AddComponent(Teuchos::rcp(new StringConditionComponent("projection","none",
      Teuchos::tuple<std::string>("none","xy","yz","xz"),
      Teuchos::tuple<std::string>("none","xy","yz","xz"),
      true)));

  condlist.push_back(volumeconstraint);

  /*--------------------------------------------------------------------*/
  // volume constraint penalty

  Teuchos::RCP<ConditionDefinition> volumeconstraintpen =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE VOLUME CONSTRAINT 3D PEN",
                                         "VolumeConstraint_3D_Pen",
                                         "Surface Volume Constraint Penalty",
                                         DRT::Condition::VolumeConstraint_3D_pen,
                                         true,
                                         DRT::Condition::Surface));

  volumeconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("penalty")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("rho")));
  volumeconstraintpen->AddComponent(Teuchos::rcp(new StringConditionComponent("projection","none",
      Teuchos::tuple<std::string>("none","xy","yz","xz"),
      Teuchos::tuple<std::string>("none","xy","yz","xz"),
      true)));

  condlist.push_back(volumeconstraintpen);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<ConditionDefinition> areaconstraint =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D",
                                         "AreaConstraint_3D",
                                         "Surface Area Constraint",
                                         DRT::Condition::AreaConstraint_3D,
                                         true,
                                         DRT::Condition::Surface));

  areaconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areaconstraint->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  areaconstraint->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));

  condlist.push_back(areaconstraint);

  /*--------------------------------------------------------------------*/
  // area constraint penalty

  Teuchos::RCP<ConditionDefinition> areaconstraintpen =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE AREA CONSTRAINT 3D PEN",
                                         "AreaConstraint_3D_Pen",
                                         "Surface Area Constraint Penalty",
                                         DRT::Condition::AreaConstraint_3D_pen,
                                         true,
                                         DRT::Condition::Surface));

  areaconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("penalty")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new RealConditionComponent("rho")));
  areaconstraintpen->AddComponent(Teuchos::rcp(new StringConditionComponent("projection","none",
      Teuchos::tuple<std::string>("none","xy","yz","xz"),
      Teuchos::tuple<std::string>("none","xy","yz","xz"),
      true)));

  condlist.push_back(areaconstraintpen);


  /*--------------------------------------------------------------------*/
  // volume monitor

  Teuchos::RCP<ConditionDefinition> volumemonitor =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE VOLUME MONITOR 3D",
                                         "VolumeMonitor_3D",
                                         "Surface Volume Monitor",
                                         DRT::Condition::VolumeMonitor_3D,
                                         true,
                                         DRT::Condition::Surface));

  volumemonitor->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(volumemonitor);

  /*--------------------------------------------------------------------*/
  // area monitor 3D

  Teuchos::RCP<ConditionDefinition> areamonitor =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE AREA MONITOR 3D",
                                         "AreaMonitor_3D",
                                         "Surface Area Monitor",
                                         DRT::Condition::AreaMonitor_3D,
                                         true,
                                         DRT::Condition::Surface));

  areamonitor->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areamonitor->AddComponent(Teuchos::rcp(new StringConditionComponent("projection","none",
    Teuchos::tuple<std::string>("none","xy","yz","xz"),
    Teuchos::tuple<std::string>("none","xy","yz","xz"),
    true)));

  condlist.push_back(areamonitor);

  /*--------------------------------------------------------------------*/
  // area constraint

  Teuchos::RCP<ConditionDefinition> areaconstraint2D =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE AREA CONSTRAINT 2D",
                                         "AreaConstraint_2D",
                                         "Line Area Constraint",
                                         DRT::Condition::AreaConstraint_2D,
                                         true,
                                         DRT::Condition::Line));

  areaconstraint2D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  areaconstraint2D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  areaconstraint2D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  condlist.push_back(areaconstraint2D);

  /*--------------------------------------------------------------------*/
  // area monitor 2D

  Teuchos::RCP<ConditionDefinition> areamonitor2D =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE AREA MONITOR 2D",
                                         "AreaMonitor_2D",
                                         "Line Area Monitor",
                                         DRT::Condition::AreaMonitor_2D,
                                         true,
                                         DRT::Condition::Line));

  areamonitor2D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  condlist.push_back(areamonitor2D);

  /*--------------------------------------------------------------------*/
  // Impedance condition

  Teuchos::RCP<ConditionDefinition> impedancebc =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF IMPEDANCE CONDITIONS",
                                         "ImpedanceCond",
                                         "Impedance boundary condition",
                                         DRT::Condition::ImpedanceCond,
                                         true,
                                         DRT::Condition::Surface));

  impedancebc->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  AddNamedReal(impedancebc,"timeperiod");
  impedancebc->AddComponent(Teuchos::rcp(new StringConditionComponent("tree", "lung",
                                                                      Teuchos::tuple<std::string>("lung","artery","windkessel","windkessel_freq_indp","resistive"),
                                                                      Teuchos::tuple<std::string>("lung","artery","windkessel","windkessel_freq_indp","resistive"),
                                                                      true)));
  AddNamedReal(impedancebc,"termradius");
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("k1")));
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("k2")));
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("k3")));

  AddNamedReal(impedancebc,"stiffness");
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("h1")));
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("h2")));
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("h3")));

  condlist.push_back(impedancebc);

  /*--------------------------------------------------------------------*/
  // Frequency indipendent precalibrated impedance condition

  Teuchos::RCP<ConditionDefinition> impedance_calb_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF IMPEDANCE CALIBRATION CONDITIONS",
                                         "ImpedanceCalbCond",
                                         "Impedance calibration boundary condition",
                                         DRT::Condition::Impedance_Calb_Cond,
                                         true,
                                         DRT::Condition::Surface));

  impedance_calb_bc->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  impedance_calb_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("tree", "windkessel_freq_indp",
                                                                            Teuchos::tuple<std::string>("windkessel_freq_indp"),
                                                                            Teuchos::tuple<std::string>("windkessel_freq_indp"),
                                                                            true)));
  AddNamedReal(impedance_calb_bc,"Pin_n");
  AddNamedReal(impedance_calb_bc,"Pin_np");
  AddNamedReal(impedance_calb_bc,"Pc_n");
  AddNamedReal(impedance_calb_bc,"Pc_np");

  condlist.push_back(impedance_calb_bc);
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D for a node over a plane

  Teuchos::RCP<ConditionDefinition> nodeonplaneconst3D =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE MULTIPNT CONSTRAINT 3D",
                                         "MPC_NodeOnPlane_3D",
                                         "Node on Plane Constraint",
                                         DRT::Condition::MPC_NodeOnPlane_3D,
                                         false,
                                         DRT::Condition::Surface));

  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("planeNodes",3)));
  nodeonplaneconst3D->AddComponent(Teuchos::rcp(new StringConditionComponent("control","rel",
      Teuchos::tuple<std::string>("rel","abs"),
      Teuchos::tuple<std::string>("rel","abs"),
      true)));
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
//  nodeonlineconst3D->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curveNodes",2)));
//  condlist.push_back(nodeonlineconst3D);
//
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 3D, moving all constraint nodes synchronously

   Teuchos::RCP<ConditionDefinition> nodemasterconst3D =
     Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D",
                                          "MPC_NormalComponent_3D",
                                          "Node on Plane Constraint",
                                          DRT::Condition::MPC_NormalComponent_3D,
                                          false,
                                          DRT::Condition::Surface));

   nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("masterNode")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("direction",3)));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new StringConditionComponent("value","disp",
       Teuchos::tuple<std::string>("disp","x"),
       Teuchos::tuple<std::string>("disp","x"),
       true)));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new StringConditionComponent("control","rel",
       Teuchos::tuple<std::string>("rel","abs"),
       Teuchos::tuple<std::string>("rel","abs"),
       true)));
   condlist.push_back(nodemasterconst3D);

   /*--------------------------------------------------------------------*/
   // Multi point constraint in 3D, moving all constraint nodes synchronously, penalty based

    Teuchos::RCP<ConditionDefinition> nodemasterconst3Dpen =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D PEN",
                                           "MPC_NormalComponent_3D_Pen",
                                           "Node on Plane Constraint Penalty",
                                           DRT::Condition::MPC_NormalComponent_3D_pen,
                                           false,
                                           DRT::Condition::Surface));

    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealConditionComponent("penalty")));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new IntConditionComponent("masterNode")));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("direction",3)));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new StringConditionComponent("value","disp",
        Teuchos::tuple<std::string>("disp","x"),
        Teuchos::tuple<std::string>("disp","x"),
        true)));
    nodemasterconst3Dpen->AddComponent(Teuchos::rcp(new StringConditionComponent("control","rel",
        Teuchos::tuple<std::string>("rel","abs"),
        Teuchos::tuple<std::string>("rel","abs"),
        true)));
    condlist.push_back(nodemasterconst3Dpen);
  /*--------------------------------------------------------------------*/
  // Multi point constraint in 2D for a node on a line
  Teuchos::RCP<ConditionDefinition> nodeonlineconst2D =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE MULTIPNT CONSTRAINT 2D",
                                           "MPC_NodeOnLine_2D",
                                           "Node on Line Constraint",
                                           DRT::Condition::MPC_NodeOnLine_2D,
                                           false,
                                           DRT::Condition::Line));

  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new RealConditionComponent("amplitude")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("constrNode 1")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("constrNode 2")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new IntConditionComponent("constrNode 3")));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new StringConditionComponent("control value","dist",
        Teuchos::tuple<std::string>("dist","angle"),
        Teuchos::tuple<std::string>("dist","angle"),
        true)));
  nodeonlineconst2D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  condlist.push_back(nodeonlineconst2D);


  /*--------------------------------------------------------------------*/
   // Electrode kinetics (Electrochemistry)
  {
    std::vector<Teuchos::RCP<CondCompBundle> > reactionmodel;

    // Butler-Volmer
    std::vector<Teuchos::RCP<ConditionComponent> > butlervolmer;
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
    butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
    butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("i0")));
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer",
                                                             butlervolmer,
                                                             INPAR::SCATRA::butler_volmer)));

    // Butler-Volmer Yang
    // parameter are identical to Butler-Volmer
    std::vector<Teuchos::RCP<ConditionComponent> > butlervolmeryang;
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
    butlervolmeryang.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
    butlervolmeryang.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    butlervolmeryang.push_back(Teuchos::rcp(new RealConditionComponent("i0")));
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    butlervolmeryang.push_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    butlervolmeryang.push_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    butlervolmeryang.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    butlervolmeryang.push_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Yang1997",
                                                             butlervolmeryang,
                                                             INPAR::SCATRA::butler_volmer_yang1997)));

    // Tafel kinetics
    std::vector<Teuchos::RCP<ConditionComponent> > tafel;
    tafel.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha")));
    tafel.push_back(Teuchos::rcp(new RealConditionComponent("alpha")));
    tafel.push_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    tafel.push_back(Teuchos::rcp(new RealConditionComponent("i0")));
    tafel.push_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    tafel.push_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    tafel.push_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    tafel.push_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    tafel.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    tafel.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("Tafel",
                                                             tafel,
                                                             INPAR::SCATRA::tafel)));

    // linear kinetics
    std::vector<Teuchos::RCP<ConditionComponent> > linear;
    linear.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha")));
    linear.push_back(Teuchos::rcp(new RealConditionComponent("alpha")));
    linear.push_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    linear.push_back(Teuchos::rcp(new RealConditionComponent("i0")));
    linear.push_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    linear.push_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    linear.push_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    linear.push_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    linear.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    linear.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    linear.push_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("linear",
                                                                   linear,
                                                                   INPAR::SCATRA::linear)));

    // Butler-Volmer-Newman: "Newman (book), 2004, p. 213, eq. 8.26"
    //                       "Wittmann (Bachelor thesis), 2011, p. 15, eq. 2.30"
    std::vector<Teuchos::RCP<ConditionComponent> > bvnewman;
    bvnewman.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_a")));
    bvnewman.push_back(Teuchos::rcp(new RealConditionComponent("k_a")));
    bvnewman.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_c")));
    bvnewman.push_back(Teuchos::rcp(new RealConditionComponent("k_c")));
    bvnewman.push_back(Teuchos::rcp(new SeparatorConditionComponent("beta")));
    bvnewman.push_back(Teuchos::rcp(new RealConditionComponent("beta")));
    bvnewman.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    bvnewman.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    bvnewman.push_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Newman",
                                                              bvnewman,
                                                              INPAR::SCATRA::butler_volmer_newman)));

    // Butler-Volmer-Newman: "Bard (book), 2001, p. 99, eq. 3.4.10"
    //                       "Wittmann (Bachelor thesis), 2011, p. 16, eq. 2.32"
    std::vector<Teuchos::RCP<ConditionComponent> > bvbard;
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("e0")));
    bvbard.push_back(Teuchos::rcp(new RealConditionComponent("e0")));
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("k0")));
    bvbard.push_back(Teuchos::rcp(new RealConditionComponent("k0")));
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("beta")));
    bvbard.push_back(Teuchos::rcp(new RealConditionComponent("beta")));
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("c_c0")));
    bvbard.push_back(Teuchos::rcp(new RealConditionComponent("c_c0")));
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("c_a0")));
    bvbard.push_back(Teuchos::rcp(new RealConditionComponent("c_a0")));
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    bvbard.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    bvbard.push_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Bard",
                                                              bvbard,
                                                              INPAR::SCATRA::butler_volmer_bard)));

    // Nernst  equation:
    std::vector<Teuchos::RCP<ConditionComponent> > nernst;
    nernst.push_back(Teuchos::rcp(new SeparatorConditionComponent("e0")));
    nernst.push_back(Teuchos::rcp(new RealConditionComponent("e0")));
    nernst.push_back(Teuchos::rcp(new SeparatorConditionComponent("c0")));
    nernst.push_back(Teuchos::rcp(new RealConditionComponent("c0")));
    nernst.push_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    nernst.push_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    reactionmodel.push_back(Teuchos::rcp(new CondCompBundle("Nernst",
                                                              nernst,
                                                              INPAR::SCATRA::nernst)));

    // input: stoichiometry for reaction mechanism (IntRealBundle)
    // definition separator for int vectors
    std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepveccomp;
    intsepveccomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("stoich")));

    // definition int vectors
    std::vector<Teuchos::RCP<IntVectorConditionComponent> > intveccomp;
    intveccomp.push_back(Teuchos::rcp(new IntVectorConditionComponent("stoich",2)));

    // definition separator for real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepveccomp;

    // definition real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<RealVectorConditionComponent> > realveccomp;


    std::vector<Teuchos::RCP<ConditionComponent> > elechemcomponents;
    elechemcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("pot")));
    elechemcomponents.push_back(Teuchos::rcp(new RealConditionComponent("pot")));
    elechemcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("curve")));
    elechemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
    elechemcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));
    elechemcomponents.push_back(Teuchos::rcp(new IntRealBundle(
        "intreal bundle",
        Teuchos::rcp(new IntConditionComponent("numscal")),
        intsepveccomp,
        intveccomp,
        realsepveccomp,
        realveccomp)));
    elechemcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
    elechemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
    elechemcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("zero_cur")));
    elechemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("zero_cur")));
    elechemcomponents.push_back(Teuchos::rcp(new CondCompBundleSelector(
        "kinetic model bundle",
        Teuchos::rcp(new StringConditionComponent(
           "kinetic model",
           "Butler-Volmer",
           Teuchos::tuple<std::string>("Butler-Volmer","Butler-Volmer-Yang1997","Tafel","linear",
                                       "Butler-Volmer-Newman","Butler-Volmer-Bard","Nernst","zero"),
           Teuchos::tuple<int>(INPAR::SCATRA::butler_volmer,
                               INPAR::SCATRA::butler_volmer_yang1997,
                               INPAR::SCATRA::tafel,
                               INPAR::SCATRA::linear,
                               INPAR::SCATRA::butler_volmer_newman,
                               INPAR::SCATRA::butler_volmer_bard,
                               INPAR::SCATRA::nernst,
                               INPAR::SCATRA::zero))),
        reactionmodel)));


    Teuchos::RCP<ConditionDefinition> lineelec =
      Teuchos::rcp(new ConditionDefinition("ELECTRODE KINETICS LINE CONDITIONS",
                                           "ElectrodeKinetics",
                                           "Line Electrode Kinetics",
                                           DRT::Condition::ElectrodeKinetics,
                                           true,
                                           DRT::Condition::Line));
    Teuchos::RCP<ConditionDefinition> surfelec =
      Teuchos::rcp(new ConditionDefinition("ELECTRODE KINETICS SURF CONDITIONS",
                                           "ElectrodeKinetics",
                                           "Surface Electrode Kinetics",
                                           DRT::Condition::ElectrodeKinetics,
                                           true,
                                           DRT::Condition::Surface));

    for (unsigned i=0; i<elechemcomponents.size(); ++i)
    {
      lineelec->AddComponent(elechemcomponents[i]);
      surfelec->AddComponent(elechemcomponents[i]);
    }

    condlist.push_back(lineelec);
    condlist.push_back(surfelec);
  }

  /*-------------------------------------------------------------------*/
  // Scalar transport reaction terms

  //Separators for "numscal xx stoich xx xx xx":
  // definition separator for int vectors
  std::vector<Teuchos::RCP<SeparatorConditionComponent> > HSTRintsepveccompstoich;
  HSTRintsepveccompstoich.push_back(Teuchos::rcp(new SeparatorConditionComponent("stoich")));
  // definition int vectors
  std::vector<Teuchos::RCP<IntVectorConditionComponent> > HSTRintveccompstoich;
  HSTRintveccompstoich.push_back(Teuchos::rcp(new IntVectorConditionComponent("stoich",2)));
  // definition separator for real vectors: length of the real vector is zero -> nothing is read
  std::vector<Teuchos::RCP<SeparatorConditionComponent> > HSTRrealsepveccompstoich;
  // definition real vectors: length of the real vector is zero -> nothing is read
  std::vector<Teuchos::RCP<RealVectorConditionComponent> > HSTRrealveccompstoich;


  std::vector<Teuchos::RCP<ConditionComponent> > homvolscatracoupcomp;

  homvolscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")) );
  homvolscatracoupcomp.push_back(Teuchos::rcp(new IntRealBundle(
                                 "intreal bundle numscal",
                                Teuchos::rcp(new IntConditionComponent("numscal")),
                                HSTRintsepveccompstoich,
                                HSTRintveccompstoich,
                                HSTRrealsepveccompstoich,
                                HSTRrealveccompstoich)) );
  homvolscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("reaccoeff")) );
  homvolscatracoupcomp.push_back(Teuchos::rcp(new RealConditionComponent("reaccoeff")) );
  homvolscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("coupling")) );
  homvolscatracoupcomp.push_back(Teuchos::rcp(new StringConditionComponent("coupling",
                                  "simple_multiplicative",
                                  Teuchos::tuple<std::string>("simple_multiplicative","other"),
                                  /*Teuchos::tuple<int>(INPAR::SCATRA::simple_multiplicative,
                                            INPAR::SCATRA::other)*/
                                  Teuchos::tuple<std::string>("simple_multiplicative","other"))) );
  homvolscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("reacstart",true)) ); //optional input
  homvolscatracoupcomp.push_back(Teuchos::rcp(new RealConditionComponent("reacstart")) );

  Teuchos::RCP<ConditionDefinition> homvolscatracoup =
    Teuchos::rcp(new ConditionDefinition("DESIGN HOMOGENEOUS SCATRA COUPLING VOLUME CONDITIONS",
                 "HomoScaTraCoupling",
                 "Homogenous ScaTra Volume Coupling",
                 DRT::Condition::HomoScaTraCoupling,
                 true,
                 DRT::Condition::Volume) );

  for (unsigned i=0; i<homvolscatracoupcomp.size(); ++i)
  {
	  homvolscatracoup->AddComponent(homvolscatracoupcomp[i]);
  }
  condlist.push_back(homvolscatracoup);


  /*--------------------------------------------------------------------*/
  // Boundary flux evaluation condition for scalar transport

  std::vector<Teuchos::RCP<ConditionComponent> > fluxeval;
  //OUTPUT:
  // - default: scalar flux is only evaluated for standard degrees of freedom (NUMSCAL)
  //            this was the only option until 05/13
  // - all: scalar flux is evaluated for all degrees of freedom (NUMDOF - including potential dof)
  //        e.g. diffusion conduction formulation with div i -> potential dof gives current flow
  fluxeval.push_back(Teuchos::rcp(new SeparatorConditionComponent("OUTPUT")));
  fluxeval.push_back(Teuchos::rcp(new StringConditionComponent("output",
                        "standard",
                        Teuchos::tuple<std::string>("standard","alldof"),
                        Teuchos::tuple<int>(INPAR::SCATRA::fluxeval_standard,
                                            INPAR::SCATRA::fluxeval_alldof))));

  Teuchos::RCP<ConditionDefinition> linebndryfluxeval =
    Teuchos::rcp(new ConditionDefinition("SCATRA FLUX CALC LINE CONDITIONS",
                                         "ScaTraFluxCalc",
                                         "Scalar Transport Boundary Flux Calculation",
                                         DRT::Condition::ScaTraFluxCalc,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfbndryfluxeval =
    Teuchos::rcp(new ConditionDefinition("SCATRA FLUX CALC SURF CONDITIONS",
                                         "ScaTraFluxCalc",
                                         "Scalar Transport Boundary Flux Calculation",
                                         DRT::Condition::ScaTraFluxCalc,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<fluxeval.size(); ++i)
  {
    linebndryfluxeval->AddComponent(fluxeval[i]);
    surfbndryfluxeval->AddComponent(fluxeval[i]);
  }

  condlist.push_back(linebndryfluxeval);
  condlist.push_back(surfbndryfluxeval);

  /*--------------------------------------------------------------------*/
  // Coupling of different scalar transport fields

  std::vector<Teuchos::RCP<ConditionComponent> > scatracoupcomponents;

  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPID")));
  scatracoupcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FIELDNUM")));
  scatracoupcomponents.push_back(Teuchos::rcp(new IntConditionComponent("field number")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("PERMCOEF")));
  scatracoupcomponents.push_back(Teuchos::rcp(new RealConditionComponent("permeability coefficient")));

  Teuchos::RCP<ConditionDefinition> surfscatracoup =
    Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA COUPLING SURF CONDITIONS",
                                         "ScaTraCoupling",
                                         "ScaTra Coupling",
                                         DRT::Condition::ScaTraCoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<scatracoupcomponents.size(); ++i)
  {
    surfscatracoup->AddComponent(scatracoupcomponents[i]);
  }

  condlist.push_back(surfscatracoup);

  std::vector<Teuchos::RCP<ConditionComponent> > linescatracoupcomp;

  linescatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPID")));
  linescatracoupcomp.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  linescatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("FIELDNUM")));
  linescatracoupcomp.push_back(Teuchos::rcp(new IntConditionComponent("field number")));
  linescatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("PERMCOEF")));
  linescatracoupcomp.push_back(Teuchos::rcp(new RealConditionComponent("permeability coefficient")));

  Teuchos::RCP<ConditionDefinition> linescatracoup =
    Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA COUPLING LINE CONDITIONS",
                                         "ScaTraCoupling",
                                         "ScaTra Coupling",
                                         DRT::Condition::ScaTraCoupling,
                                         true,
                                         DRT::Condition::Line));

  for (unsigned i=0; i<linescatracoupcomp.size(); ++i)
  {
    linescatracoup->AddComponent(linescatracoupcomp[i]);
  }

  condlist.push_back(linescatracoup);

  std::vector<Teuchos::RCP<ConditionComponent> > pointscatracoupcomp;

  pointscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPID")));
  pointscatracoupcomp.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  pointscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("FIELDNUM")));
  pointscatracoupcomp.push_back(Teuchos::rcp(new IntConditionComponent("field number")));
  pointscatracoupcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("PERMCOEF")));
  pointscatracoupcomp.push_back(Teuchos::rcp(new RealConditionComponent("permeability coefficient")));

  Teuchos::RCP<ConditionDefinition> pointscatracoup =
    Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA COUPLING POINT CONDITIONS",
                                         "ScaTraCoupling",
                                         "ScaTra Coupling",
                                         DRT::Condition::ScaTraCoupling,
                                         true,
                                         DRT::Condition::Point));

  for (unsigned i=0; i<pointscatracoupcomp.size(); ++i)
  {
    pointscatracoup->AddComponent(pointscatracoupcomp[i]);
  }

  condlist.push_back(pointscatracoup);

  /*--------------------------------------------------------------------*/
  // flow rate through line

  std::vector<Teuchos::RCP<ConditionComponent> > lineflowratecomponents;
  lineflowratecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> lineflowrate =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLOW RATE LINE CONDITIONS",
                                         "LineFlowRate",
                                         "Line Flow Rate",
                                         DRT::Condition::FlowRateThroughLine_2D,
                                         true,
                                         DRT::Condition::Line));

  for (unsigned i=0; i<lineflowratecomponents.size(); ++i)
  {
    lineflowrate->AddComponent(lineflowratecomponents[i]);
  }
  condlist.push_back(lineflowrate);

  /*--------------------------------------------------------------------*/
  // flow rate through surface

  std::vector<Teuchos::RCP<ConditionComponent> > flowratecomponents;
  flowratecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfflowrate =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLOW RATE SURF CONDITIONS",
                                         "SurfFlowRate",
                                         "Surface Flow Rate",
                                         DRT::Condition::FlowRateThroughSurface_3D,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<flowratecomponents.size(); ++i)
  {
    surfflowrate->AddComponent(flowratecomponents[i]);
  }
  condlist.push_back(surfflowrate);

  /*--------------------------------------------------------------------*/
  // impuls rate through surface

  std::vector<Teuchos::RCP<ConditionComponent> > impulsratecomponents;
  impulsratecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfimpulsrate =
    Teuchos::rcp(new ConditionDefinition("DESIGN IMPULS RATE SURF CONDITIONS",
                                         "SurfImpulsRate",
                                         "Surface Impuls Rate",
                                         DRT::Condition::ImpulsRateThroughSurface_3D,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<impulsratecomponents.size(); ++i)
  {
    surfimpulsrate->AddComponent(impulsratecomponents[i]);
  }
  condlist.push_back(surfimpulsrate);

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


    std::vector<Teuchos::RCP<ConditionComponent> > rigidbodymodecomponents;

    rigidbodymodecomponents.push_back(
      Teuchos::rcp(
        new StringConditionComponent(
          "discretization",
          "fluid",
          Teuchos::tuple<std::string>("fluid","scatra","solid"),
          Teuchos::tuple<std::string>("fluid","scatra","solid"))));

    // input: trigger as int-vector using an IntRealBundle
    // definition separator for int vectors
    std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepveccomp;
    intsepveccomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
    // definition int vectors
    std::vector<Teuchos::RCP<IntVectorConditionComponent> > intveccomp;
    intveccomp.push_back(Teuchos::rcp(new IntVectorConditionComponent("ONOFF",1)));
    // definition separator for real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepveccomp;
    // definition real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<RealVectorConditionComponent> > realveccomp;

    rigidbodymodecomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMMODES")));
    rigidbodymodecomponents.push_back(Teuchos::rcp(new IntRealBundle(
        "intreal bundle",
        Teuchos::rcp(new IntConditionComponent("NUMMODES")),
        intsepveccomp,
        intveccomp,
        realsepveccomp,
        realveccomp)));

    rigidbodymodecomponents.push_back(
      Teuchos::rcp(
        new StringConditionComponent(
          "weight vector definition",
          "integration",
          Teuchos::tuple<std::string>("integration","pointvalues"),
          Teuchos::tuple<std::string>("integration","pointvalues"))));

    Teuchos::RCP<ConditionDefinition> surfrigidbodymode =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF MODE FOR KRYLOV SPACE PROJECTION",
                                           "KrylovSpaceProjection",
                                           "Surface mode for Krylov space projection",
                                           DRT::Condition::SurfaceModeKrylovProjection,
                                           true,
                                           DRT::Condition::Surface));

    Teuchos::RCP<ConditionDefinition> volrigidbodymode =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOL MODE FOR KRYLOV SPACE PROJECTION",
                                           "KrylovSpaceProjection",
                                           "Volume mode for Krylov space projection",
                                           DRT::Condition::VolumeModeKrylovProjection,
                                           true,
                                           DRT::Condition::Volume));

    for (unsigned i=0; i<rigidbodymodecomponents.size(); ++i)
    {
      surfrigidbodymode->AddComponent(rigidbodymodecomponents[i]);
      volrigidbodymode ->AddComponent(rigidbodymodecomponents[i]);
    }
    condlist.push_back(surfrigidbodymode);
    condlist.push_back(volrigidbodymode);
  }

  /*--------------------------------------------------------------------*/
  // inverse analysis fitted surface

  std::vector<Teuchos::RCP<ConditionComponent> > invanacomponents;
  invanacomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfinvana =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE INV ANALYSIS",
                                         "SurfInvAna",
                                         "Inverse Analysis Surface",
                                         DRT::Condition::InvAnaSurface,
                                         true,
                                         DRT::Condition::Surface));

  surfinvana->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(surfinvana);

  /*--------------------------------------------------------------------*/
  // 1D-Artery connector condition

  Teuchos::RCP<ConditionDefinition> art_connection_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE 1D ARTERY JUNCTION CONDITIONS",
                                         "ArtJunctionCond",
                                         "Artery junction boundary condition",
                                         DRT::Condition::ArtJunctionCond,
                                         true,
                                         DRT::Condition::Point));

  art_connection_bc->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  art_connection_bc->AddComponent(Teuchos::rcp(new RealConditionComponent("Kr")));

  condlist.push_back(art_connection_bc);

  /*--------------------------------------------------------------------*/
  // Export 1D-Arterial network in gnuplot format

  Teuchos::RCP<ConditionDefinition> art_write_gnuplot_c =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE EXPORT 1D-ARTERIAL NETWORK GNUPLOT FORMAT",
                                         "ArtWriteGnuplotCond",
                                         "Artery write gnuplot format condition",
                                         DRT::Condition::ArtWriteGnuplotCond,
                                         false,
                                         DRT::Condition::Line));

  art_write_gnuplot_c->AddComponent(Teuchos::rcp(new IntConditionComponent("ArteryNumber")));

  condlist.push_back(art_write_gnuplot_c);

  /*--------------------------------------------------------------------*/
  // 1D artery prescribed BC

  Teuchos::RCP<ConditionDefinition> art_in_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE 1D ARTERY PRESCRIBED CONDITIONS",
                                         "ArtPrescribedCond",
                                         "Artery prescribed boundary condition",
                                         DRT::Condition::ArtPrescribedCond,
                                         true,
                                         DRT::Condition::Point));

  art_in_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("boundarycond", "flow",
    Teuchos::tuple<std::string>("flow","pressure","velocity","area","characteristicWave"),
    Teuchos::tuple<std::string>("flow","pressure","velocity","area","characteristicWave"),
    true)));
  art_in_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("type", "forced",
    Teuchos::tuple<std::string>("forced","absorbing"),
    Teuchos::tuple<std::string>("forced","absorbing"),
    true)));

  std::vector<Teuchos::RCP<ConditionComponent> > artinletcomponents;
  artinletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",2)));
  artinletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",2,true,true)));
  for (unsigned i=0; i<artinletcomponents.size(); ++i)
    art_in_bc->AddComponent(artinletcomponents[i]);

  condlist.push_back(art_in_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery reflective BC
  Teuchos::RCP<ConditionDefinition> art_rf_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE 1D ARTERY REFLECTIVE CONDITIONS",
                                         "ArtRfCond",
                                         "Artery reflection condition",
                                         DRT::Condition::ArtRfCond,
                                         true,
                                         DRT::Condition::Point));

  std::vector<Teuchos::RCP<ConditionComponent> > artrfcomponents;
  artrfcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  artrfcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  for (unsigned i=0; i<artrfcomponents.size(); ++i)
    art_rf_bc->AddComponent(artrfcomponents[i]);

  condlist.push_back(art_rf_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery windkessel BC
  Teuchos::RCP<ConditionDefinition> art_wk_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE 1D ARTERY WINDKESSEL CONDITIONS",
                                         "ArtWkCond",
                                         "Artery windkessel condition",
                                         DRT::Condition::ArtWkCond,
                                         true,
                                         DRT::Condition::Point));

  std::vector<Teuchos::RCP<ConditionComponent> > artwkcomponents;

  art_wk_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("intigrationType", "ExplicitWindkessel",
    Teuchos::tuple<std::string>("ExplicitWindkessel", "ImpedaceWindkessel"),
    Teuchos::tuple<std::string>("ExplicitWindkessel", "ImpedaceWindkessel"),
    true)));

  art_wk_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("windkesselType", "RCR",
    Teuchos::tuple<std::string>("R","RC", "RCR", "RCRL", "none"),
    Teuchos::tuple<std::string>("R","RC", "RCR", "RCRL", "none"),
    true)));

  artwkcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",5)));
  artwkcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",5,true,true)));
  for (unsigned i=0; i<artwkcomponents.size(); ++i)
    art_wk_bc->AddComponent(artwkcomponents[i]);

  condlist.push_back(art_wk_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery in/out condition

  Teuchos::RCP<ConditionDefinition> art_in_outlet_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE 1D ARTERY IN_OUTLET CONDITIONS",
                                         "ArtInOutCond",
                                         "Artery terminal in_outlet condition",
                                         DRT::Condition::ArtInOutletCond,
                                         true,
                                         DRT::Condition::Point));

  art_in_outlet_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("terminaltype", "inlet",
    Teuchos::tuple<std::string>("inlet","outlet"),
    Teuchos::tuple<std::string>("inlet","outlet"),
    true)));

  condlist.push_back(art_in_outlet_bc);
  /*--------------------------------------------------------------------*/
  // 1D artery scalar transport condition
  Teuchos::RCP<ConditionDefinition> art_scatra_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE 1D ARTERY SCATRA PRESCRIBED CONDITIONS",
                                         "ArtPrescribedScatraCond",
                                         "Artery prescribed scatra boundary condition",
                                         DRT::Condition::ArtPrescribedScatraCond,
                                         true,
                                         DRT::Condition::Point));

  std::vector<Teuchos::RCP<ConditionComponent> > artscatracomponents;
  artscatracomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  artscatracomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  for (unsigned i=0; i<artscatracomponents.size(); ++i)
    art_scatra_bc->AddComponent(artscatracomponents[i]);

  condlist.push_back(art_scatra_bc);


  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<ConditionDefinition> art_red_to_3d_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE REDUCED D To 3D FLOW COUPLING CONDITIONS",
                                         "Art_redD_3D_CouplingCond",
                                         "Artery reduced D 3D coupling condition",
                                         DRT::Condition::ArtRedTo3DCouplingCond,
                                         true,
                                         DRT::Condition::Point));

  art_red_to_3d_bc->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  art_red_to_3d_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("CouplingType", "forced",
                                                                           Teuchos::tuple<std::string>("forced","absorbing"),
                                                                           Teuchos::tuple<std::string>("forced","absorbing"),
                                                                           true)));

  art_red_to_3d_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("ReturnedVariable", "pressure",
                                                                           Teuchos::tuple<std::string>("pressure","flow"),
                                                                           Teuchos::tuple<std::string>("pressure","flow"),
                                                                           true)));
  AddNamedReal(art_red_to_3d_bc,"Tolerance");
  AddNamedInt (art_red_to_3d_bc,"MaximumIterations");

  condlist.push_back(art_red_to_3d_bc);

  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<ConditionDefinition> art_3d_to_red_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF 3D To REDUCED D FLOW COUPLING CONDITIONS",
                                         "Art_3D_redD_CouplingCond",
                                         "Artery 3D reduced D coupling condition",
                                         DRT::Condition::Art3DToRedCouplingCond,
                                         true,
                                         DRT::Condition::Surface));

  art_3d_to_red_bc->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  art_3d_to_red_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("ReturnedVariable", "flow",
                                                                           Teuchos::tuple<std::string>("pressure","flow"),
                                                                           Teuchos::tuple<std::string>("pressure","flow"),
                                                                           true)));
  AddNamedReal(art_3d_to_red_bc,"Tolerance");
  AddNamedInt (art_3d_to_red_bc,"MaximumIterations");

  condlist.push_back(art_3d_to_red_bc);


  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<ConditionDefinition> windkessel_optim_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF WINDKESSEL OPTIMIZATION CONDITIONS",
                                         "Windkessel_Optimization_Cond",
                                         "windkessel optimization condition",
                                         DRT::Condition::WindkesselOptimCond,
                                         true,
                                         DRT::Condition::Surface));

  windkessel_optim_bc->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));


  windkessel_optim_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("ObjectiveFunction", "Psys_Pdia",
                                                                              Teuchos::tuple<std::string>("Psys_Pdia","Psys_Pdia_Pavg"),
                                                                              Teuchos::tuple<std::string>("Psys_Pdia","Psys_Pdia_Pavg"),
                                                                              true)));
  windkessel_optim_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("DesignVariables", "R_C",
                                                                              Teuchos::tuple<std::string>("R_C","R1_R2_C"),
                                                                              Teuchos::tuple<std::string>("R_C","R1_R2_C"),
                                                                              true)));
  AddNamedReal(windkessel_optim_bc,"Psystolic");
  AddNamedReal(windkessel_optim_bc,"Pdiastolic");
  AddNamedReal(windkessel_optim_bc,"R1R2_ratio");
  AddNamedReal(windkessel_optim_bc,"Tolerance");

  condlist.push_back(windkessel_optim_bc);


  /*--------------------------------------------------------------------*/
  // Coupling of 3D tissue models and reduced-D airway tree

  std::vector<Teuchos::RCP<ConditionComponent> > redairtiscomponents;

  redairtiscomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> surfredairtis =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF TISSUE REDAIRWAY CONDITIONS",
                                         "SurfaceNeumann",
                                         "tissue RedAirway coupling surface condition",
                                         DRT::Condition::RedAirwayTissue,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<redairtiscomponents.size(); ++i)
  {
    surfredairtis->AddComponent(redairtiscomponents[i]);
  }

  condlist.push_back(surfredairtis);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a four-element Windkessel - mhv 11/13

  Teuchos::RCP<ConditionDefinition> windkesselcondition =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF WINDKESSEL CONDITIONS",
                                         "WindkesselStdStructureCond",
                                         "Surface Windkessel",
                                         DRT::Condition::WindkesselStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselcondition,"id");
  AddNamedReal(windkesselcondition,"C");
  AddNamedReal(windkesselcondition,"R_p");
  AddNamedReal(windkesselcondition,"Z_c");
  AddNamedReal(windkesselcondition,"L");
  AddNamedReal(windkesselcondition,"p_ref");
  AddNamedReal(windkesselcondition,"p_0");

  condlist.push_back(windkesselcondition);


  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a special heart valve arterial four-element Windkessel - mhv 02/14

  Teuchos::RCP<ConditionDefinition> windkesselheartvalvearterialcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF HEART VALVE ARTERIAL WINDKESSEL CONDITIONS",
                                         "WindkesselHeartValveArterialStructureCond",
                                         "Surface heart valve arterial Windkessel",
                                         DRT::Condition::WindkesselHeartValveArterialStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselheartvalvearterialcond,"id");
  AddNamedReal(windkesselheartvalvearterialcond,"R_av_max");
  AddNamedReal(windkesselheartvalvearterialcond,"R_av_min");
  AddNamedReal(windkesselheartvalvearterialcond,"R_mv_max");
  AddNamedReal(windkesselheartvalvearterialcond,"R_mv_min");
  AddNamedReal(windkesselheartvalvearterialcond,"k_p");
  AddNamedReal(windkesselheartvalvearterialcond,"C");
  AddNamedReal(windkesselheartvalvearterialcond,"R_p");
  AddNamedReal(windkesselheartvalvearterialcond,"Z_c");
  AddNamedReal(windkesselheartvalvearterialcond,"L");
  AddNamedReal(windkesselheartvalvearterialcond,"p_ref");
  AddNamedReal(windkesselheartvalvearterialcond,"p_ar_0");
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("P_AT")));
  AddNamedReal(windkesselheartvalvearterialcond,"fac");
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("crv")));
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("VALVE")));
  windkesselheartvalvearterialcond->AddComponent(Teuchos::rcp(new StringConditionComponent("valvelaw", "smooth",
                                                                                       Teuchos::tuple<std::string>("smooth","pwlin"),
                                                                                       Teuchos::tuple<std::string>("smooth","pwlin"),
                                                                                       true)));

  condlist.push_back(windkesselheartvalvearterialcond);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and a heart valve arterial Windkessel accounting for proximal and distal arterial pressure
  // formulation proposed by Cristobal Bertoglio - mhv 03/14

  Teuchos::RCP<ConditionDefinition> windkesselheartvalvearterialproxdistcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF HEART VALVE ARTERIAL PROX DIST WINDKESSEL CONDITIONS",
                                         "WindkesselHeartValveArterialProxDistStructureCond",
                                         "Surface heart valve arterial proximal and distal Windkessel",
                                         DRT::Condition::WindkesselHeartValveArterialProxDistStructure,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselheartvalvearterialproxdistcond,"id");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_av_max");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_av_min");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_mv_max");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_mv_min");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"k_p");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"L_arp");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"C_arp");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_arp");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"C_ard");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"R_ard");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"p_ref");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"p_arp_0");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"y_arp_0");
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"p_ard_0");
  windkesselheartvalvearterialproxdistcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("P_AT")));
  AddNamedReal(windkesselheartvalvearterialproxdistcond,"fac");
  windkesselheartvalvearterialproxdistcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("crv")));
  windkesselheartvalvearterialproxdistcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));

  condlist.push_back(windkesselheartvalvearterialproxdistcond);

  /*--------------------------------------------------------------------*/
  // Monolithic coupling of structure and Windkessel: Neumann coupling surface - mhv 11/13

  Teuchos::RCP<ConditionDefinition> windkesselstructurecouplingcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF WINDKESSEL STRUCTURE COUPLING CONDITIONS",
                                         "SurfaceNeumann",
                                         "structure windkessel coupling surface condition",
                                         DRT::Condition::WindkesselStructureCoupling,
                                         true,
                                         DRT::Condition::Surface));

  AddNamedInt(windkesselstructurecouplingcond,"coupling_id");

  condlist.push_back(windkesselstructurecouplingcond);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways

  std::vector<Teuchos::RCP<ConditionComponent> > noderedairtiscomponents;

  noderedairtiscomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> noderedairtis =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE TISSUE REDAIRWAY CONDITIONS",
                                         "RedAirwayPrescribedCond",
                                         "tissue RedAirway coupling node condition",
                                         DRT::Condition::RedAirwayNodeTissue,
                                         true,
                                         DRT::Condition::Point));


  for (unsigned i=0; i<noderedairtiscomponents.size(); ++i)
  {
    noderedairtis->AddComponent(noderedairtiscomponents[i]);
  }

  condlist.push_back(noderedairtis);



  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways

  Teuchos::RCP<ConditionDefinition> raw_in_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE Reduced D AIRWAYS PRESCRIBED CONDITIONS",
                                         "RedAirwayPrescribedCond",
                                         "Reduced d airway prescribed boundary condition",
                                         DRT::Condition::RedAirwayPrescribedCond,
                                         true,
                                         DRT::Condition::Point));

  raw_in_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("boundarycond", "flow",
    Teuchos::tuple<std::string>("flow","pressure", "VolumeDependentPleuralPressure"),
    Teuchos::tuple<std::string>("flow","pressure", "VolumeDependentPleuralPressure"),
    true)));

  std::vector<Teuchos::RCP<ConditionComponent> > redairwayinletcomponents;
  redairwayinletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  redairwayinletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",2,true,true)));
  redairwayinletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",1, false, false, true)));
  for (unsigned i=0; i<redairwayinletcomponents.size(); ++i)
    raw_in_bc->AddComponent(redairwayinletcomponents[i]);

  condlist.push_back(raw_in_bc);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways external pressure

  Teuchos::RCP<ConditionDefinition> raw_pext_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE Reduced D AIRWAYS PRESCRIBED EXTERNAL PRESSURE CONDITIONS",
                                         "RedAirwayPrescribedExternalPressure",
                                         "Reduced d airway prescribed external pressure boundary condition",
                                         DRT::Condition::RedAirwayPrescribedExternalPressure,
                                         true,
                                         DRT::Condition::Line));

  raw_pext_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("boundarycond","ExternalPressure",
    Teuchos::tuple<std::string>("ExternalPressure"),
    Teuchos::tuple<std::string>("ExternalPressure"),
    true)));

  std::vector<Teuchos::RCP<ConditionComponent> > redairwaypextcomponents;
  redairwaypextcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  redairwaypextcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",2,true,true)));
  redairwaypextcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",1, false, false, true)));
  for (unsigned i=0; i<redairwaypextcomponents.size(); ++i)
    raw_pext_bc->AddComponent(redairwaypextcomponents[i]);

  condlist.push_back(raw_pext_bc);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional scalar transpirt in airways

  Teuchos::RCP<ConditionDefinition> raw_in_scatra_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE Reduced D AIRWAYS PRESCRIBED SCATRA CONDITIONS",
                                         "RedAirwayPrescribedScatraCond",
                                         "Reduced d airway prescribed scatra boundary condition",
                                         DRT::Condition::RedAirwayPrescribedScatraCond,
                                         true,
                                         DRT::Condition::Point));

  std::vector<Teuchos::RCP<ConditionComponent> > redairwayinletscatracomponents;
  redairwayinletscatracomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  redairwayinletscatracomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  redairwayinletscatracomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",1, false, false, true)));
  for (unsigned i=0; i<redairwayinletscatracomponents.size(); ++i)
    raw_in_scatra_bc->AddComponent(redairwayinletscatracomponents[i]);

  condlist.push_back(raw_in_scatra_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for intial values of the scalar transport in reduced dimensional airways

  Teuchos::RCP<ConditionDefinition> raw_int_scatra_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE Reduced D AIRWAYS INITIAL SCATRA CONDITIONS",
                                         "RedAirwayInitialScatraCond",
                                         "Reduced d airway initial scatra boundary condition",
                                         DRT::Condition::RedAirwayInitialScatraCond,
                                         true,
                                         DRT::Condition::Line));

  raw_int_scatra_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("scalar", "O2",
    Teuchos::tuple<std::string>("O2","CO2"),
    Teuchos::tuple<std::string>("O2","CO2"),
    true)));

  AddNamedReal(raw_int_scatra_bc,"CONCENTRATION");
  condlist.push_back(raw_int_scatra_bc);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions of scatra exchange
  Teuchos::RCP<ConditionDefinition> scatra_exchange_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE Reduced D AIRWAYS SCATRA EXCHANGE CONDITIONS",
                                         "RedAirwayScatraExchangeCond",
                                         "scatra exchange condition",
                                         DRT::Condition::RedAirwayScatraExchangeCond,
                                         true,
                                         DRT::Condition::Line));

  scatra_exchange_cond->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(scatra_exchange_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<ConditionDefinition> scatra_hemoglobin_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE Reduced D AIRWAYS HEMOGLOBIN CONDITIONS",
                                         "RedAirwayScatraHemoglobinCond",
                                         "scatra hemoglobin condition",
                                         DRT::Condition::RedAirwayScatraHemoglobinCond,
                                         false,
                                         DRT::Condition::Line));

  AddNamedReal(scatra_hemoglobin_cond,"INITIAL_CONCENTRATION");
  condlist.push_back(scatra_hemoglobin_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<ConditionDefinition> scatra_air_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE Reduced D AIRWAYS AIR CONDITIONS",
                                         "RedAirwayScatraAirCond",
                                         "scatra air condition",
                                         DRT::Condition::RedAirwayScatraAirCond,
                                         false,
                                         DRT::Condition::Line));

  AddNamedReal(scatra_air_cond,"INITIAL_CONCENTRATION");
  condlist.push_back(scatra_air_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<ConditionDefinition> scatra_capillary_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE Reduced D AIRWAYS CAPILLARY CONDITIONS",
                                         "RedAirwayScatraCapillaryCond",
                                         "scatra capillary condition",
                                         DRT::Condition::RedAirwayScatraCapillaryCond,
                                         false,
                                         DRT::Condition::Line));

  condlist.push_back(scatra_capillary_cond);

  /*--------------------------------------------------------------------*/
  // Prescribed Ventilator BC for reduced dimensional airways

  Teuchos::RCP<ConditionDefinition> raw_vent_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN NODE Reduced D AIRWAYS VENTILATOR CONDITIONS",
                                         "RedAirwayVentilatorCond",
                                         "Reduced d airway prescribed ventilator condition",
                                         DRT::Condition::RedAirwayVentilatorCond,
                                         true,
                                         DRT::Condition::Point));

  raw_vent_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("phase1", "flow",
                                                                      Teuchos::tuple<std::string>("flow","pressure"),
                                                                      Teuchos::tuple<std::string>("flow","pressure"),
                                                                      true)));

  raw_vent_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("Phase1Smoothness", "smooth",
                                                                      Teuchos::tuple<std::string>("smooth","discontinous"),
                                                                      Teuchos::tuple<std::string>("smooth","discontinous"),
                                                                      true)));

  raw_vent_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("phase2", "pressure",
                                                                      Teuchos::tuple<std::string>("pressure","flow"),
                                                                      Teuchos::tuple<std::string>("pressure","flow"),
                                                                      true)));

  raw_vent_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("Phase2Smoothness", "smooth",
                                                                      Teuchos::tuple<std::string>("smooth","discontinous"),
                                                                      Teuchos::tuple<std::string>("smooth","discontinous"),
                                                                      true)));

  AddNamedReal(raw_vent_bc,"period");
  AddNamedReal(raw_vent_bc,"phase1_period");
  AddNamedReal(raw_vent_bc,"smoothness_period1");
  AddNamedReal(raw_vent_bc,"smoothness_period2");

  std::vector<Teuchos::RCP<ConditionComponent> > redairwayventcomponents;
  redairwayventcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",2)));
  redairwayventcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",2,true,true)));
  for (unsigned i=0; i<redairwayventcomponents.size(); ++i)
    raw_vent_bc->AddComponent(redairwayventcomponents[i]);

  condlist.push_back(raw_vent_bc);




  /*--------------------------------------------------------------------*/
  // Prescribed volume dependent pleural pressure for reduced dimensional airways

  Teuchos::RCP<ConditionDefinition> raw_volPpl_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE REDUCED D AIRWAYS VOL DEPENDENT PLEURAL PRESSURE CONDITIONS",
                                         "RedAirwayVolDependentPleuralPressureCond",
                                         "Reduced D airways volume-dependent peural pressure condition",
                                         DRT::Condition::RedAirwayVolDependentPleuralPressureCond,
                                         true,
                                         DRT::Condition::Line));

  raw_volPpl_bc->AddComponent(Teuchos::rcp(new StringConditionComponent("TYPE", "Exponential",
                                                                      Teuchos::tuple<std::string>("Exponential","Polynomial"),
                                                                      Teuchos::tuple<std::string>("Exponential","Polynomial"),
                                                                      true)));

  AddNamedReal(raw_volPpl_bc,"TLC");
  AddNamedReal(raw_volPpl_bc,"VFR");

  AddNamedReal(raw_volPpl_bc,"P_PLEURAL_0");
  AddNamedReal(raw_volPpl_bc,"P_PLEURAL_LIN");
  AddNamedReal(raw_volPpl_bc,"P_PLEURAL_NONLIN");
  AddNamedReal(raw_volPpl_bc,"TAU");


  std::vector<Teuchos::RCP<ConditionComponent> > raw_volPpl_bc_components;
  raw_volPpl_bc_components.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  raw_volPpl_bc_components.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  for (unsigned i=0; i<raw_volPpl_bc_components.size(); ++i)
    raw_volPpl_bc->AddComponent(raw_volPpl_bc_components[i]);

  condlist.push_back(raw_volPpl_bc);

  /*--------------------------------------------------------------------*/
  // Evaluate lung volume condition for reduced dimensional airways

  Teuchos::RCP<ConditionDefinition> raw_eval_lungV_bc =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE REDUCED D AIRWAYS EVALUATE LUNG VOLUME CONDITIONS",
                                         "RedAirwayEvalLungVolCond",
                                         "Reduced D airways evaluate lung volume condition",
                                         DRT::Condition::RedAirwayEvalLungVolCond,
                                         true,
                                         DRT::Condition::Line));


  condlist.push_back(raw_eval_lungV_bc);

  /*--------------------------------------------------------------------*/
  // Volumetric surface flow profile condition
  Teuchos::RCP<ConditionDefinition> volumetric_surface_flow_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF VOLUMETRIC FLOW CONDITIONS",
                                         "VolumetricSurfaceFlowCond",
                                         "volumetric surface flow condition",
                                         DRT::Condition::VolumetricSurfaceFlowCond,
                                         true,
                                         DRT::Condition::Surface));

  std::vector<Teuchos::RCP<ConditionComponent> > inflownormalcomponents;

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("ConditionType", "WOMERSLEY",
                                                                                       Teuchos::tuple<std::string>("WOMERSLEY","POLYNOMIAL"),
                                                                                       Teuchos::tuple<std::string>("WOMERSLEY","POLYNOMIAL"),
                                                                                       true)));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("prebiased", "NOTPREBIASED",
                                                                                       Teuchos::tuple<std::string>("NOTPREBIASED","PREBIASED","FORCED"),
                                                                                       Teuchos::tuple<std::string>("NOTPREBIASED","PREBIASED","FORCED"),
                                                                                       true)));


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("FlowType", "InFlow",
                                                                                       Teuchos::tuple<std::string>("InFlow","OutFlow"),
                                                                                       Teuchos::tuple<std::string>("InFlow","OutFlow"),
                                                                                       true)));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("CorrectionFlag", "WithOutCorrection",
                                                                                       Teuchos::tuple<std::string>("WithOutCorrection","WithCorrection"),
                                                                                       Teuchos::tuple<std::string>("WithOutCorrection","WithCorrection"),
                                                                                       true)));
  AddNamedReal(volumetric_surface_flow_cond,"Period");
  AddNamedInt (volumetric_surface_flow_cond,"Order");
  AddNamedInt (volumetric_surface_flow_cond,"Harmonics");

  std::vector<Teuchos::RCP<ConditionComponent> > inflowprofilecomponents;
  inflowprofilecomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  inflowprofilecomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  for (unsigned i=0; i<inflowprofilecomponents.size(); ++i)
    volumetric_surface_flow_cond->AddComponent(inflowprofilecomponents[i]);


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("NORMAL", "SelfEvaluateNormal",
                                                                                       Teuchos::tuple<std::string>("SelfEvaluateNormal","UsePrescribedNormal"),
                                                                                       Teuchos::tuple<std::string>("SelfEvaluateNormal","UsePrescribedNormal"),
                                                                                       true)));


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n1")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n2")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n3")));


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("CenterOfMass", "SelfEvaluateCenterOfMass",
                                                                                       Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass","UsePrescribedCenterOfMass"),
                                                                                       Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass","UsePrescribedCenterOfMass"),
                                                                                       true)));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c1")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c2")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c3")));

  condlist.push_back(volumetric_surface_flow_cond);



  /*--------------------------------------------------------------------*/
  // Volumetric flow border nodes condition

  Teuchos::RCP<ConditionDefinition> volumetric_border_nodes_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE VOLUMETRIC FLOW BORDER NODES",
                                         "VolumetricFlowBorderNodesCond",
                                         "volumetric flow border nodes condition",
                                         DRT::Condition::VolumetricFlowBorderNodes,
                                         true,
                                         DRT::Condition::Line));

  volumetric_border_nodes_cond->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));


  condlist.push_back(volumetric_border_nodes_cond);

  /*--------------------------------------------------------------------*/
  // Volumetric surface total traction corrector
  Teuchos::RCP<ConditionDefinition> total_traction_correction_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF TOTAL TRACTION CORRECTION CONDITIONS",
                                         "TotalTractionCorrectionCond",
                                         "total traction correction condition",
                                         DRT::Condition::TotalTractionCorrectionCond,
                                         true,
                                         DRT::Condition::Surface));


  total_traction_correction_cond->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("ConditionType", "POLYNOMIAL",
                                                                                         Teuchos::tuple<std::string>("POLYNOMIAL","WOMERSLEY"),
                                                                                         Teuchos::tuple<std::string>("POLYNOMIAL","WOMERSLEY"),
                                                                                         true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("prebiased", "NOTPREBIASED",
                                                                                         Teuchos::tuple<std::string>("NOTPREBIASED","PREBIASED","FORCED"),
                                                                                         Teuchos::tuple<std::string>("NOTPREBIASED","PREBIASED","FORCED"),
                                                                                         true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("FlowType", "InFlow",
                                                                                         Teuchos::tuple<std::string>("InFlow","OutFlow"),
                                                                                         Teuchos::tuple<std::string>("InFlow","OutFlow"),
                                                                                         true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("CorrectionFlag", "WithOutCorrection",
                                                                                         Teuchos::tuple<std::string>("WithOutCorrection","WithCorrection"),
                                                                                         Teuchos::tuple<std::string>("WithOutCorrection","WithCorrection"),
                                                                                         true)));
  AddNamedReal(total_traction_correction_cond,"Period");
  AddNamedInt (total_traction_correction_cond,"Order");
  AddNamedInt (total_traction_correction_cond,"Harmonics");

  std::vector<Teuchos::RCP<ConditionComponent> > flowbiasingcomponents;
  flowbiasingcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",1)));
  flowbiasingcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",1,true,true)));
  for (unsigned i=0; i<flowbiasingcomponents.size(); ++i)
  total_traction_correction_cond->AddComponent(flowbiasingcomponents[i]);

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("NORMAL", "SelfEvaluateNormal",
                                                                                         Teuchos::tuple<std::string>("SelfEvaluateNormal","UsePrescribedNormal"),
                                                                                         Teuchos::tuple<std::string>("SelfEvaluateNormal","UsePrescribedNormal"),
                                                                                         true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n1")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n2")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n3")));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("CenterOfMass", "SelfEvaluateCenterOfMass",
                                                                                         Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass","UsePrescribedCenterOfMass"),
                                                                                         Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass","UsePrescribedCenterOfMass"),
                                                                                         true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c1")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c2")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c3")));

  condlist.push_back(total_traction_correction_cond);

  /*--------------------------------------------------------------------*/
  // Volumetric flow traction correction border nodes condition

  Teuchos::RCP<ConditionDefinition> traction_corrector_border_nodes_cond =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE TOTAL TRACTION CORRECTION BORDER NODES",
                                         "TotalTractionCorrectionBorderNodesCond",
                                         "total traction correction border nodes condition",
                                         DRT::Condition::TotalTractionCorrectionBorderNodes,
                                         true,
                                         DRT::Condition::Line));

  traction_corrector_border_nodes_cond->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));


  condlist.push_back(traction_corrector_border_nodes_cond);


  /*--------------------------------------------------------------------*/
  // Convective heat transfer (Newton's law of heat transfer)

  std::vector<Teuchos::RCP<ConditionComponent> > thermoconvectcomponents;

  // decide here if approximation is sufficient
  // --> Tempn (old temperature T_n)
  // or if the exact solution is needed
  // --> Tempnp (current temperature solution T_n+1) with linearisation
  thermoconvectcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "temperature state","Tempnp",
        Teuchos::tuple<std::string>("Tempnp","Tempn"),
        Teuchos::tuple<std::string>("Tempnp","Tempn"))));
  // heat transfer coefficient h
  thermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("coeff")));
  thermoconvectcomponents.push_back(Teuchos::rcp(new RealConditionComponent("coeff")));
  // surrounding (fluid) temperature T_oo
  thermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("surtemp")));
  thermoconvectcomponents.push_back(Teuchos::rcp(new RealConditionComponent("surtemp")));
  // time curve to increase the surrounding (fluid) temperature T_oo in time
  thermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("surtempcurve")));
  thermoconvectcomponents.push_back(Teuchos::rcp(new IntConditionComponent("surtempcurve",true,true)));
  // time curve to increase the complete boundary condition, i.e., the heat flux
  thermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("curve")));
  thermoconvectcomponents.push_back(Teuchos::rcp(new IntConditionComponent("curve",true,true)));

  Teuchos::RCP<ConditionDefinition> linethermoconvect =
    Teuchos::rcp(new ConditionDefinition("DESIGN THERMO CONVECTION LINE CONDITIONS",
                                         "ThermoConvections",
                                         "Line Thermo Convections",
                                         DRT::Condition::ThermoConvections,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfthermoconvect =
    Teuchos::rcp(new ConditionDefinition("DESIGN THERMO CONVECTION SURF CONDITIONS",
                                         "ThermoConvections",
                                         "Surface Thermo Convections",
                                         DRT::Condition::ThermoConvections,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<thermoconvectcomponents.size(); ++i)
  {
    linethermoconvect->AddComponent(thermoconvectcomponents[i]);
    surfthermoconvect->AddComponent(thermoconvectcomponents[i]);
  }

  condlist.push_back(linethermoconvect);
  condlist.push_back(surfthermoconvect);


  /*--------------------------------------------------------------------*/

  // structural spring dashpot boundary condition (spring and dashpot in parallel) - mhv 01/14

  Teuchos::RCP<ConditionDefinition> springdashpotcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF SPRING DASHPOT CONDITIONS",
                                         "SpringDashpot",
                                         "Spring Dashpot",
                                         DRT::Condition::SpringDashpot,
                                         true,
                                         DRT::Condition::Surface));


  AddNamedReal(springdashpotcond,"SPRING_STIFF_TENS");
  AddNamedReal(springdashpotcond,"SPRING_STIFF_COMP");
  AddNamedReal(springdashpotcond,"SPRING_OFFSET");
  AddNamedReal(springdashpotcond,"DASHPOT_VISCOSITY");
  springdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("DIRECTION")));
  springdashpotcond->AddComponent(Teuchos::rcp(new StringConditionComponent("direction", "all",
                                                                                       Teuchos::tuple<std::string>("all","refsurfnormal"),
                                                                                       Teuchos::tuple<std::string>("all","refsurfnormal"),
                                                                                       true)));

  condlist.push_back(springdashpotcond);

  /*--------------------------------------------------------------------*/
  // no penetration for darcy flow in porous media

  Teuchos::RCP<ConditionDefinition> nopenetration_surf =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NORMAL NO PENETRATION CONDITION",
                                         "NoPenetration",
                                         "No Penetration",
                                         DRT::Condition::NoPenetration,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(nopenetration_surf);

  /*--------------------------------------------------------------------*/
  // no penetration for darcy flow in porous media

  Teuchos::RCP<ConditionDefinition> nopenetration_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE NORMAL NO PENETRATION CONDITION",
                                         "NoPenetration",
                                         "No Penetration",
                                         DRT::Condition::NoPenetration,
                                         true,
                                         DRT::Condition::Line));

  condlist.push_back(nopenetration_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of coupling terms in porous media

  Teuchos::RCP<ConditionDefinition> porocoupling_vol =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOLUME POROCOUPLING CONDITION",
                                         "PoroCoupling",
                                         "Poro Coupling",
                                         DRT::Condition::PoroCoupling,
                                         true,
                                         DRT::Condition::Volume));

  condlist.push_back(porocoupling_vol);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of coupling terms in porous media

  Teuchos::RCP<ConditionDefinition> porocoupling_surf =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE POROCOUPLING CONDITION",
                                         "PoroCoupling",
                                         "Poro Coupling",
                                         DRT::Condition::PoroCoupling,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(porocoupling_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropartint_surf =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PORO PARTIAL INTEGRATION",
                                         "PoroPartInt",
                                         "Poro Partial Integration",
                                         DRT::Condition::PoroPartInt,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(poropartint_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropartint_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO PARTIAL INTEGRATION",
                                         "PoroPartInt",
                                         "Poro Partial Integration",
                                         DRT::Condition::PoroPartInt,
                                         true,
                                         DRT::Condition::Line));

  condlist.push_back(poropartint_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropresint_surf =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PORO PRESSURE INTEGRATION",
                                         "PoroPresInt",
                                         "Poro Pressure Integration",
                                         DRT::Condition::PoroPresInt,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(poropresint_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropresint_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO PRESSURE INTEGRATION",
                                         "PoroPresInt",
                                         "Poro Pressure Integration",
                                         DRT::Condition::PoroPresInt,
                                         true,
                                         DRT::Condition::Line));

  condlist.push_back(poropresint_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems
  // necessary where neumann term needs to be integrated in interface
  // elements which share a node with the fpsi interface. Tangential
  // Beaver-Joseph-Condition must not be overwritten by prescribed value!

  Teuchos::RCP<ConditionDefinition> neumannintegration_surf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NEUMANN INTEGRATION",
                                           "NeumannIntegration",
                                           "Neumann Integration",
                                           DRT::Condition::NeumannIntegration,
                                           true,
                                           DRT::Condition::Surface));

  condlist.push_back(neumannintegration_surf);

   /*--------------------------------------------------------------------*/
   // condition for evaluation of boundary terms in fpsi problems

   Teuchos::RCP<ConditionDefinition> neumannintegration_line =
     Teuchos::rcp(new ConditionDefinition("DESIGN LINE NEUMANN INTEGRATION",
                                          "NeumannIntegration",
                                          "Neumann Integration",
                                          DRT::Condition::NeumannIntegration,
                                          true,
                                          DRT::Condition::Line));

   condlist.push_back(neumannintegration_line);


  /*--------------------------------------------------------------------*/
  // particle inflow condition

  std::vector<Teuchos::RCP<ConditionComponent> > particleinflowcomponents;
  // two vertices describing the bounding box for the inflow
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex1")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex1",3)));
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex2")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex2",3)));
  // number of particles to inflow in each direction in bounding box
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("num_per_dir")));
  particleinflowcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("num_per_dir",3)));
  // particle inflow velocity
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_vel")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("inflow_vel",3)));
  // particle inflow velocity can be superposed with a time curve
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_vel_curve")));
  particleinflowcomponents.push_back(Teuchos::rcp(new IntConditionComponent("inflow_vel_curve")));
  // inflow frequency of particles
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_freq")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("inflow_freq")));


  Teuchos::RCP<ConditionDefinition> particlecond =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INFLOW CONDITION",
                                         "ParticleInflow",
                                         "Particle Inflow Condition",
                                         DRT::Condition::ParticleInflow,
                                         false,
                                         DRT::Condition::Particle));


  for (unsigned i=0; i<particleinflowcomponents.size(); ++i)
  {
    particlecond->AddComponent(particleinflowcomponents[i]);
  }

  condlist.push_back(particlecond);

  /*--------------------------------------------------------------------*/
  // particle periodic boundary condition

  std::vector<Teuchos::RCP<ConditionComponent> > particlepbccomponents;
  // two vertices describing the bounding box for the pbc
  particlepbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  particlepbccomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("ONOFF",3)));
  particlepbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("boundaries")));
  particlepbccomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("boundaries",6)));

  Teuchos::RCP<ConditionDefinition> particlepbccond =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE PERIODIC BOUNDARY CONDITION",
                                         "ParticlePeriodic",
                                         "Particle Periodic Boundary Condition",
                                         DRT::Condition::ParticlePeriodic,
                                         false,
                                         DRT::Condition::Particle));


  for (unsigned i=0; i<particlepbccomponents.size(); ++i)
  {
    particlepbccond->AddComponent(particlepbccomponents[i]);
  }

  condlist.push_back(particlepbccond);

  /*--------------------------------------------------------------------*/
  // particle init radius condition

  std::vector<Teuchos::RCP<ConditionComponent> > particleinitradiuscomponents;

  particleinitradiuscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("SCALAR")));
  particleinitradiuscomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("SCALAR",1)));

  particleinitradiuscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  particleinitradiuscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("FUNCT",1)));

  Teuchos::RCP<ConditionDefinition> particleradiuscond =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INIT RADIUS CONDITIONS",
                                         "InitialParticleRadius",
                                         "Particle Initial Radius Condition",
                                         DRT::Condition::ParticleInitRadius,
                                         false,
                                         DRT::Condition::Particle));


  for (unsigned i=0; i<particleinitradiuscomponents.size(); ++i)
  {
    particleradiuscond->AddComponent(particleinitradiuscomponents[i]);
  }

  condlist.push_back(particleradiuscond);

  /*--------------------------------------------------------------------*/
  // particle wall condition

  std::vector<Teuchos::RCP<ConditionComponent> > particlewallcomponents;
  particlewallcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> surfpartwall =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PARTICLE WALL",
                                         "ParticleWall",
                                         "Wall for particle collisions",
                                         DRT::Condition::ParticleWall,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<particlewallcomponents.size(); ++i)
  {
    surfpartwall->AddComponent(particlewallcomponents[i]);
  }

  condlist.push_back(surfpartwall);

  /*--------------------------------------------------------------------*/
  // Surface current evaluation condition

  std::vector<Teuchos::RCP<ConditionComponent> > surfcurrcomponents;
  surfcurrcomponents.push_back(Teuchos::rcp(new IntConditionComponent("matching id")));

  Teuchos::RCP<ConditionDefinition> surfcurrcond =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE CURRENT EVALUATION CONDITION",
                                         "SurfaceCurrent",
                                         "Surface current",
                                         DRT::Condition::SurfaceCurrent,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<surfcurrcomponents.size(); ++i)
  {
    surfcurrcond->AddComponent(surfcurrcomponents[i]);
  }

  condlist.push_back(surfcurrcond);

  /*--------------------------------------------------------------------*/
  // level-set condition for contact points

  Teuchos::RCP<ConditionDefinition> linelscontact =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE LEVEL SET CONTACT CONDITION",
                                         "LsContact",
                                         "level-set condition for contact points",
                                         DRT::Condition::LsContact,
                                         false,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> pointlscontact =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT LEVEL SET CONTACT CONDITION",
                                         "LsContact",
                                         "level-set condition for contact points",
                                         DRT::Condition::LsContact,
                                         false,
                                         DRT::Condition::Point));

  condlist.push_back(linelscontact);
  condlist.push_back(pointlscontact);

  /*--------------------------------------------------------------------*/
  // absorbing boundary condition for acoustic problems
  // line
  Teuchos::RCP<ConditionDefinition> absorbing_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN ABSORBING LINE CONDITIONS",
                                         "Absorbing",
                                         "Absorbing line for acoustics",
                                         DRT::Condition::Absorb,
                                         true,
                                         DRT::Condition::Line));
  condlist.push_back(absorbing_line);

  // surface
  Teuchos::RCP<ConditionDefinition> absorbing_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN ABSORBING SURF CONDITIONS",
                                         "Absorbing",
                                         "Absorbing surface for acoustics",
                                         DRT::Condition::Absorb,
                                         true,
                                         DRT::Condition::Surface));
  condlist.push_back(absorbing_surface);

  /*--------------------------------------------------------------------*/
  // monitor condition for acoustic problems
  // line
  Teuchos::RCP<ConditionDefinition> pressmon_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN PRESSURE MONITOR LINE CONDITIONS",
                                         "PressureMonitor",
                                         "Pressure monitor line for acoustics",
                                         DRT::Condition::PressureMonitor,
                                         true,
                                         DRT::Condition::Line));
  condlist.push_back(pressmon_line);

  // surface
  Teuchos::RCP<ConditionDefinition> pressmon_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN PRESSURE MONITOR SURF CONDITIONS",
                                         "PressureMonitor",
                                         "Pressure monitor surface for acoustics",
                                         DRT::Condition::PressureMonitor,
                                         true,
                                         DRT::Condition::Surface));
  condlist.push_back(pressmon_surface);

  /*--------------------------------------------------------------------*/
  return vc;

}


