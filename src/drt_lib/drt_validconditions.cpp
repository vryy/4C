
#ifdef CCADISCRET

#include "drt_validconditions.H"
#include "drt_conditiondefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyConditionDefinitions(std::ostream& stream,
                                                std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  for (unsigned i=0; i<condlist.size(); ++i)
  {
    condlist[i]->Print(stream,NULL,true);
  }
}


/*----------------------------------------------------------------------*/
//! Print function to be called from C
/*----------------------------------------------------------------------*/
extern "C"
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
        Teuchos::tuple<std::string>("Live","Dead","PrescribedDomainLoad","constHydro_z","increaseHydro_z","orthopressure","LAS"),
        Teuchos::tuple<std::string>("neum_live","neum_dead","pres_domain_load","neum_consthydro_z","neum_increhydro_z","neum_orthopressure","neum_LAS"),
        true)));
  neumanncomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "surface","Mid",
        Teuchos::tuple<std::string>("Mid","Top","Bot"),
        Teuchos::tuple<std::string>("mid","top","bot"))));

  // optional
  neumanncomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",6,false,false,true)));

  Teuchos::RCP<ConditionDefinition> pointneumann =
    Teuchos::rcp(new ConditionDefinition("DESIGN POINT NEUMANN CONDITIONS",
                                         "PointNeumann",
                                         "Point Neumann",
                                         DRT::Condition::PointNeumann,
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

  for (unsigned i=0; i<neumanncomponents.size(); ++i)
  {
    pointneumann->AddComponent(neumanncomponents[i]);
    lineneumann->AddComponent(neumanncomponents[i]);
    surfneumann->AddComponent(neumanncomponents[i]);
    volneumann->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(pointneumann);
  condlist.push_back(lineneumann);
  condlist.push_back(surfneumann);
  condlist.push_back(volneumann);

  /*--------------------------------------------------------------------*/
  // Dirichlet

  std::vector<Teuchos::RCP<ConditionComponent> > dirichletcomponents;

  dirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",6)));
  dirichletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val",6)));
  dirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("curve",6,true,true)));

  // optional
  dirichletcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("funct",6,false,false,true)));

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

  for (unsigned i=0; i<dirichletcomponents.size(); ++i)
  {
    pointdirichlet->AddComponent(dirichletcomponents[i]);
    linedirichlet->AddComponent(dirichletcomponents[i]);
    surfdirichlet->AddComponent(dirichletcomponents[i]);
    voldirichlet->AddComponent(dirichletcomponents[i]);

    pointaledirichlet->AddComponent(dirichletcomponents[i]);
    linealedirichlet->AddComponent(dirichletcomponents[i]);
    surfaledirichlet->AddComponent(dirichletcomponents[i]);
  }

  condlist.push_back(pointdirichlet);
  condlist.push_back(linedirichlet);
  condlist.push_back(surfdirichlet);
  condlist.push_back(voldirichlet);

  condlist.push_back(pointaledirichlet);
  condlist.push_back(linealedirichlet);
  condlist.push_back(surfaledirichlet);

  /*--------------------------------------------------------------------*/
  // contact

  std::vector<Teuchos::RCP<ConditionComponent> > contactcomponents;

  contactcomponents.push_back(Teuchos::rcp(new IntConditionComponent("contact id")));
  contactcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Side","Master",
        Teuchos::tuple<std::string>("Master","Slave"),
        Teuchos::tuple<std::string>("Master","Slave"))));

  Teuchos::RCP<ConditionDefinition> linecontact =
    Teuchos::rcp(new ConditionDefinition("CONTACT CONDITIONS 2D",
                                         "Contact",
                                         "Line Contact",
                                         DRT::Condition::Contact,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfcontact =
    Teuchos::rcp(new ConditionDefinition("CONTACT CONDITIONS 3D",
                                         "Contact",
                                         "Surface Contact",
                                         DRT::Condition::Contact,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<contactcomponents.size(); ++i)
  {
    linecontact->AddComponent(contactcomponents[i]);
    surfcontact->AddComponent(contactcomponents[i]);
  }

  condlist.push_back(linecontact);
  condlist.push_back(surfcontact);

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
        Teuchos::tuple<std::string>("xy","yx","yz","zy","xz","zx"),
        Teuchos::tuple<std::string>("xy","xy","yz","yz","xz","xz"))));

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
  // FSI

  std::vector<Teuchos::RCP<ConditionComponent> > fsicomponents;

  fsicomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  fsicomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "field","structure",
        Teuchos::tuple<std::string>("structure","fluid","ale"),
        Teuchos::tuple<std::string>("structure","fluid","ale"))));

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
  // xfem

  std::vector<Teuchos::RCP<ConditionComponent> > xfemcomponents;

  xfemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

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

  for (unsigned i=0; i<xfemcomponents.size(); ++i)
  {
    linexfem->AddComponent(xfemcomponents[i]);
    surfxfem->AddComponent(xfemcomponents[i]);
  }

  condlist.push_back(linexfem);
  condlist.push_back(surfxfem);

  /*--------------------------------------------------------------------*/
  // surface tension

  Teuchos::RCP<ConditionDefinition> surftension =
    Teuchos::rcp(new ConditionDefinition("SURFACE TENSION CONDITIONS",
                                         "SurfaceStress",
                                         "Surface Stress (ideal water)",
                                         DRT::Condition::SurfaceTension,
                                         true,
                                         DRT::Condition::Surface));

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

  AddNamedReal(surfactant,"k1xCbulk");
  AddNamedReal(surfactant,"k2");
  AddNamedReal(surfactant,"m1");
  AddNamedReal(surfactant,"m2");
  AddNamedReal(surfactant,"gamma_0");
  AddNamedReal(surfactant,"gamma_min");
  AddNamedReal(surfactant,"gamma_min_eq");

  condlist.push_back(surfactant);

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
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE STRESS CALC CONDITIONS",
                                         "FluidStressCalc",
                                         "Line Fluid Stress Calculation",
                                         DRT::Condition::FluidStressCalc,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffluidstress =
    Teuchos::rcp(new ConditionDefinition("DESIGN FLUIDSTRESS COUPLING SURF CONDITIONS",
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
  liftdragcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("centerCoord",3)));

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

  condlist.push_back(volumeconstraint);

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

  condlist.push_back(areaconstraint);

  return vc;
}

#endif
