/*----------------------------------------------------------------------*/
/*!
\file drt_validconditions.cpp

\brief Setup of the list of valid conditions for input

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"


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

  for (unsigned i=0; i<neumanncomponents.size(); ++i)
  {
    pointneumann->AddComponent(neumanncomponents[i]);
    lineneumann->AddComponent(neumanncomponents[i]);
    surfneumann->AddComponent(neumanncomponents[i]);
    volneumann->AddComponent(neumanncomponents[i]);

    pointtransportneumann->AddComponent(neumanncomponents[i]);
    linetransportneumann->AddComponent(neumanncomponents[i]);
    surftransportneumann->AddComponent(neumanncomponents[i]);
    voltransportneumann->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(pointneumann);
  condlist.push_back(lineneumann);
  condlist.push_back(surfneumann);
  condlist.push_back(volneumann);

  condlist.push_back(pointtransportneumann);
  condlist.push_back(linetransportneumann);
  condlist.push_back(surftransportneumann);

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

  for (unsigned i=0; i<dirichletcomponents.size(); ++i)
  {
    pointdirichlet->AddComponent(dirichletcomponents[i]);
    linedirichlet->AddComponent(dirichletcomponents[i]);
    surfdirichlet->AddComponent(dirichletcomponents[i]);
    voldirichlet->AddComponent(dirichletcomponents[i]);

    pointaledirichlet->AddComponent(dirichletcomponents[i]);
    linealedirichlet->AddComponent(dirichletcomponents[i]);
    surfaledirichlet->AddComponent(dirichletcomponents[i]);

    pointtransportdirichlet->AddComponent(dirichletcomponents[i]);
    linetransportdirichlet->AddComponent(dirichletcomponents[i]);
    surftransportdirichlet->AddComponent(dirichletcomponents[i]);
  }

  condlist.push_back(pointdirichlet);
  condlist.push_back(linedirichlet);
  condlist.push_back(surfdirichlet);
  condlist.push_back(voldirichlet);

  condlist.push_back(pointaledirichlet);
  condlist.push_back(linealedirichlet);
  condlist.push_back(surfaledirichlet);

  condlist.push_back(pointtransportdirichlet);
  condlist.push_back(linetransportdirichlet);
  condlist.push_back(surftransportdirichlet);

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
  contactcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Initialization","Inactive",
        Teuchos::tuple<std::string>("Inactive","Active"),
        Teuchos::tuple<std::string>("Inactive","Active"))));

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
  // local coordinate systems

  std::vector<Teuchos::RCP<ConditionComponent> > locsyscomponents;

  locsyscomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("normal",3)));
  locsyscomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("tangent",3)));
 
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

  for (unsigned i=0; i<locsyscomponents.size(); ++i)
  {
    pointlocsys->AddComponent(locsyscomponents[i]);
    linelocsys->AddComponent(locsyscomponents[i]);
    surflocsys->AddComponent(locsyscomponents[i]);
    vollocsys->AddComponent(locsyscomponents[i]);
  }

  condlist.push_back(pointlocsys);
  condlist.push_back(linelocsys);
  condlist.push_back(surflocsys);
  condlist.push_back(vollocsys);

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
  // weak Dirichlet

  std::vector<Teuchos::RCP<ConditionComponent> > weakDirichletcomponents;

  // weak DBCs can be imposed adjoint consistent or adjoint inconsistent
  weakDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Choice of gamma parameter","adjoint-consistent",
        Teuchos::tuple<std::string>("adjoint-consistent","diffusive-optimal"),
        Teuchos::tuple<std::string>("adjoint-consistent","diffusive-optimal"))));

  // the penalty parameter could be computed dynamically (using Spaldings 
  // law of the wall) or using a fixed value
  weakDirichletcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "Definition of penalty parameter","constant",
        Teuchos::tuple<std::string>("constant","Spalding"),
        Teuchos::tuple<std::string>("constant","Spalding"))));


  // linearisation strategies --- the linearisation (i.e. the matrix 
  // contribution) of the convective term on the inflow could be
  // suppressed, since the flux is a kink function and including this one
  // might result in even worse convergence behaviour
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
  AddNamedReal(surfactant,"gamma_min_eq");

  condlist.push_back(surfactant);
  
  
  /*--------------------------------------------------------------------*/
  // Lennard Jones potential
  
  Teuchos::RCP<ConditionDefinition> lj_potential =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential",
                                         DRT::Condition::LJ_Potential,
                                         true,
                                         DRT::Condition::Surface));

  lj_potential->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential,"label");
  AddNamedReal(lj_potential,"depth");
  AddNamedReal(lj_potential,"rootDist");
  AddNamedReal(lj_potential,"cutOff");

  condlist.push_back(lj_potential);
  
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
  liftdragcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("centerCoord",3)));

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

  // axis of rotation for angular momentum calculation (only 3D)
  AddNamedVector(surfliftdrag,"axis",3);

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
  volumeconstraint->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));

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
  areaconstraint->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
  
  condlist.push_back(areaconstraint);

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
    Teuchos::tuple<std::string>("lung","artery","windkessel"),
    Teuchos::tuple<std::string>("lung","artery","windkessel"),
    true)));
  AddNamedReal(impedancebc,"termradius");
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("k1")));
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("k2")));
  impedancebc->AddComponent(Teuchos::rcp(new RealConditionComponent("k3")));

  condlist.push_back(impedancebc);

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
  // Multi point constraint in 3D, moving all constraint nodes synchronously 

   Teuchos::RCP<ConditionDefinition> nodemasterconst3D =
     Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE NORMALDIR MULTIPNT CONSTRAINT 3D",
                                          "MPC_NormalComponent_3D",
                                          "Node on Plane Constraint",
                                          DRT::Condition::MPC_NormalComponent_3D,
                                          false,
                                          DRT::Condition::Surface));

   nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new RealConditionComponent("activTime")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new IntConditionComponent("masterNode")));
   nodemasterconst3D->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("direction",3)));
   condlist.push_back(nodemasterconst3D);

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

  std::vector<Teuchos::RCP<ConditionComponent> > eleccomponents;

  eleccomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "electrode type","cathode",
        Teuchos::tuple<std::string>("cathode","anode"),
        Teuchos::tuple<std::string>("cathode","anode"))));
  eleccomponents.push_back(Teuchos::rcp(new IntConditionComponent("reactant id")));
  eleccomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "kinetic model","Butler-Volmer",
        Teuchos::tuple<std::string>("Butler-Volmer","Tafel","linear"),
        Teuchos::tuple<std::string>("Butler-Volmer","Tafel","linear"))));
  eleccomponents.push_back(Teuchos::rcp(new RealConditionComponent("pot0")));
  eleccomponents.push_back(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  eleccomponents.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
  eleccomponents.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
  eleccomponents.push_back(Teuchos::rcp(new RealConditionComponent("i0")));

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

  for (unsigned i=0; i<eleccomponents.size(); ++i)
  {
    lineelec->AddComponent(eleccomponents[i]);
    surfelec->AddComponent(eleccomponents[i]);
  }

  condlist.push_back(lineelec);
  condlist.push_back(surfelec);


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

  return vc;
}

#endif
