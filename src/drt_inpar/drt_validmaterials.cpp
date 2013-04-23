/*----------------------------------------------------------------------*/
/*!
\file drt_validmaterials.cpp

\brief Setup of the list of valid materials for input

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/
#include "drt_validmaterials.H"
#include "../drt_lib/drt_materialdefinition.H"
#include "inpar_material.H"
#include "../drt_lib/drt_colors.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyMaterialDefinitions(
  std::ostream& stream,
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >& matlist,
  bool color
  )
{
  std::string blue2light = "";
  std::string bluelight = "";
  std::string redlight = "";
  std::string yellowlight = "";
  std::string greenlight = "";
  std::string magentalight = "";
  std::string endcolor = "";

  if (color)
  {
    blue2light = BLUE2_LIGHT;
    bluelight = BLUE_LIGHT;
    redlight = RED_LIGHT;
    yellowlight = YELLOW_LIGHT;
    greenlight = GREEN_LIGHT;
    magentalight = MAGENTA_LIGHT;
    endcolor = END_COLOR;
  }

  const std::string sectionname = "MATERIALS";
  const unsigned l = sectionname.length();
  stream << redlight << "--";
  for (int i=0; i<std::max<int>(65-l,0); ++i) stream << '-';
  stream << greenlight << sectionname << endcolor << '\n';

  for (unsigned i=0; i<matlist.size(); ++i)
  {
    matlist[i]->Print(stream,NULL,color);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintMaterialDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> > > matlist = DRT::INPUT::ValidMaterials();
  DRT::INPUT::PrintEmptyMaterialDefinitions(std::cout,*matlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> > > DRT::INPUT::ValidMaterials()
{
  // a list containing all valid materials
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> > > vm
    = Teuchos::rcp(new std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >());

  // convenience
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >& matlist = *vm;


  /*----------------------------------------------------------------------*/
  // Newtonian fluid
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_fluid",
                                            "Newtonian fluid",
                                            INPAR::MAT::m_fluid));

    AddNamedReal(m,"DYNVISCOSITY","dynamic viscosity");
    AddNamedReal(m,"DENSITY","spatial mass density");
    AddNamedReal(m,"GAMMA","surface tension coeficient",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // fluid with non-linear viscosity according to Carreau-Yasuda
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_carreauyasuda",
                                            "fluid with non-linear viscosity according to Carreau-Yasuda",
                                            INPAR::MAT::m_carreauyasuda));

    AddNamedReal(m,"NU_0","zero-shear viscosity");
    AddNamedReal(m,"NU_INF","infinite-shear viscosity");
    AddNamedReal(m,"LAMBDA","characteristic time");
    AddNamedReal(m,"APARAM","constant parameter");
    AddNamedReal(m,"BPARAM","constant parameter");
    AddNamedReal(m,"DENSITY","density");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // fluid with nonlinear viscosity according to a modified power law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_modpowerlaw",
                                            "fluid with nonlinear viscosity according to a modified power law",
                                            INPAR::MAT::m_modpowerlaw));

    AddNamedReal(m,"MCONS","consistency");
    AddNamedReal(m,"DELTA","safety factor");
    AddNamedReal(m,"AEXP","exponent");
    AddNamedReal(m,"DENSITY","density");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // "yoghurt-type" fluid with nonlinear viscosity according to a power law
  // and extended by an Arrhenius-type term to account for temperature dependence
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_yoghurt",
                                            "yoghurt-type fluid with nonlinear viscosity",
                                            INPAR::MAT::m_yoghurt));

    AddNamedReal(m,"SHC","specific heat capacity at constant pressure (J/(kg*K))");
    AddNamedReal(m,"DENSITY","density");
    AddNamedReal(m,"THERMCOND","thermal conductivity (J/(m*K*s))");
    AddNamedReal(m,"STRAINRATEEXP","exponent of strain-rate term");
    AddNamedReal(m,"PREEXCON","pre-exponential constant (1/s)");
    AddNamedReal(m,"ACTENERGY","activation energy (J/kg)");
    AddNamedReal(m,"GASCON","specific gas constant R (J/(kg*K))");
    AddNamedReal(m,"DELTA","safety factor");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // fluid flow in a permeable material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_permeable",
                                            "permeability for flow in porous media",
                                            INPAR::MAT::m_permeable_fluid));

    AddNamedString(m,"TYPE","Problem type: Darcy or Darcy-Stokes","Darcy-Stokes");
    AddNamedReal(m,"DYNVISCOSITY","dynamic viscosity");
    AddNamedReal(m,"DENSITY","density");
    AddNamedReal(m,"PERMEABILITY","permeability of medium");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport material (with potential reaction coefficient)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra",
                                            "scalar transport material",
                                            INPAR::MAT::m_scatra));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"REACOEFF","reaction coefficient",true);
    AddNamedReal(m,"SCNUM","schmidt number",true);
    AddNamedReal(m,"FLDDENSITY","fluid density",true);

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // Myocard muscle material (with complicated reaction coefficient)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_myocard",
                                            "Myocard muscle material",
                                            INPAR::MAT::m_myocard));

    AddNamedReal(m,"MAIN_DIFFUSIVITY","conductivity in fiber direction");
    AddNamedReal(m,"OFF_DIFFUSIVITY","conductivity perpendicular to fiber direction");
    AddNamedReal(m,"PERTUBATION_DERIV","pertubation for calculation of reaction coefficient derivative");
    AddNamedReal(m,"U_O","base level potential");
    AddNamedReal(m,"U_U","maximum exited potential");
    AddNamedReal(m,"THETA_V","excitation threshold");
    AddNamedReal(m,"THETA_W","slow current threshold");
    AddNamedReal(m,"THETA_VM","v gate time constant threshold");
    AddNamedReal(m,"THETA_O","slow outward current time constant threshold");
    AddNamedReal(m,"TAU_V1M","v gate time constant");
    AddNamedReal(m,"TAU_V2M","v gate time constant");
    AddNamedReal(m,"TAU_VP","v gate time constant");
    AddNamedReal(m,"TAU_W1M","w gate time constant");
    AddNamedReal(m,"TAU_W2M","w gate time constant");
    AddNamedReal(m,"K_WM","w gate proportional factor");
    AddNamedReal(m,"U_WM","w gate time constant threshold");
    AddNamedReal(m,"TAU_WP","w gate time constant");
    AddNamedReal(m,"TAU_FI","fast inward current time constant");
    AddNamedReal(m,"TAU_O1","slow outward current time constant");
    AddNamedReal(m,"TAU_O2","slow outward current time constant");
    AddNamedReal(m,"TAU_SO1","slow outward current time constant");
    AddNamedReal(m,"TAU_SO2","slow outward current time constant");
    AddNamedReal(m,"K_SO","slow outward current proportional factor");
    AddNamedReal(m,"U_SO","slow outward current time constant threshold");
    AddNamedReal(m,"TAU_S1","s gate time constant");
    AddNamedReal(m,"TAU_S2","s gate time constant");
    AddNamedReal(m,"K_S","s gate proportional factor");
    AddNamedReal(m,"U_S","s gate time constant threshold");
    AddNamedReal(m,"TAU_SI","slow inward current time constant");
    AddNamedReal(m,"TAU_WINF","w gate infinity value time constant");
    AddNamedReal(m,"W_INFS","w gate infinity value");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material according to mixture-fraction approach
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_mixfrac",
                                            "material according to mixture-fraction approach",
                                            INPAR::MAT::m_mixfrac));

    AddNamedReal(m,"KINVISC","kinematic viscosity");
    AddNamedReal(m,"KINDIFF","kinematic diffusivity");
    AddNamedReal(m,"EOSFACA","equation-of-state factor a");
    AddNamedReal(m,"EOSFACB","equation-of-state factor b");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_sutherland",
                                            "material according to Sutherland law",
                                            INPAR::MAT::m_sutherland));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");
    AddNamedReal(m,"SHC","specific heat capacity at constant pressure (J/(kg*K))");
    AddNamedReal(m,"PRANUM","Prandtl number");
    AddNamedReal(m,"THERMPRESS","(initial) thermodynamic pressure (J/m^3)");
    AddNamedReal(m,"GASCON","specific gas constant R (J/(kg*K))");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with Arrhenius-type chemical
  // kinetics (species)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_arrhenius_spec",
                                            "Arrhenius-type chemical kinetics (species)",
                                            INPAR::MAT::m_arrhenius_spec));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");
    AddNamedReal(m,"SCHNUM","Schmidt number");
    AddNamedReal(m,"PREEXCON","pre-exponential constant (1/s)");
    AddNamedReal(m,"TEMPEXP","exponent of temperature dependence");
    AddNamedReal(m,"ACTEMP","activation temperature (K)");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with Arrhenius-type chemical
  // kinetics (temperature)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_arrhenius_temp",
                                            "Arrhenius-type chemical kinetics (temperature)",
                                            INPAR::MAT::m_arrhenius_temp));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");
    AddNamedReal(m,"SHC","specific heat capacity at constant pressure (J/(kg*K))");
    AddNamedReal(m,"PRANUM","Prandtl number");
    AddNamedReal(m,"REAHEAT","heat of reaction per unit mass (J/kg)");
    AddNamedReal(m,"PREEXCON","pre-exponential constant (1/s)");
    AddNamedReal(m,"TEMPEXP","exponent of temperature dependence");
    AddNamedReal(m,"ACTEMP","activation temperature (K)");
    AddNamedReal(m,"THERMPRESS","(initial) thermodynamic pressure (J/m^3)");
    AddNamedReal(m,"GASCON","specific gas constant R (J/(kg*K))");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with Arrhenius-type chemical
  // kinetics (progress variable)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_arrhenius_pv",
                                            "material with Arrhenius-type chemical kinetics (progress variable)",
                                            INPAR::MAT::m_arrhenius_pv));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");
    AddNamedReal(m,"PRANUM","Prandtl number");
    AddNamedReal(m,"PREEXCON","pre-exponential constant (1/s)");
    AddNamedReal(m,"TEMPEXP","exponent of temperature dependence");
    AddNamedReal(m,"ACTEMP","activation temperature (K)");
    AddNamedReal(m,"UNBSHC","specific heat capacity of unburnt phase (J/(kg*K))");
    AddNamedReal(m,"BURSHC","specific heat capacity of burnt phase (J/(kg*K))");
    AddNamedReal(m,"UNBTEMP","temperature of unburnt phase (K)");
    AddNamedReal(m,"BURTEMP","temperature of burnt phase (K)");
    AddNamedReal(m,"UNBDENS","density of unburnt phase (kg/m�)");
    AddNamedReal(m,"BURDENS","density of burnt phase (kg/m�)");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with simplified chemical
  // kinetics due to Ferziger and Echekki (1993) (original version and
  // modification by Poinsot and Veynante (2005)) (progress variable)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ferech_pv",
                                            "material with Ferziger-Echekki (1993) chemical kinetics (progress variable)",
                                            INPAR::MAT::m_ferech_pv));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");
    AddNamedReal(m,"PRANUM","Prandtl number");
    AddNamedReal(m,"REACRATECON","reaction-rate constant (1/s)");
    AddNamedReal(m,"PVCRIT","critical value of progress variable");
    AddNamedReal(m,"UNBSHC","specific heat capacity of unburnt phase (J/(kg*K))");
    AddNamedReal(m,"BURSHC","specific heat capacity of burnt phase (J/(kg*K))");
    AddNamedReal(m,"UNBTEMP","temperature of unburnt phase (K)");
    AddNamedReal(m,"BURTEMP","temperature of burnt phase (K)");
    AddNamedReal(m,"UNBDENS","density of unburnt phase (kg/m�)");
    AddNamedReal(m,"BURDENS","density of burnt phase (kg/m�)");
    AddNamedReal(m,"MOD","modification factor (0.0=original, 1.0=modified)");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrolyte solution (gjb 07/08)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ion",
                                            "material parameters for ion species in electrolyte solution",
                                            INPAR::MAT::m_ion));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"VALENCE","valence (= charge number)");
    AddNamedReal(m,"DENSIFICATION","densification coefficient",true);
    // via these two optional parameters we can bring the material parameters
    // of one eliminated ionic species into BACI if needed
    AddNamedReal(m,"ELIM_DIFFUSIVITY","kinematic diffusivity of elim. species",true);
    AddNamedReal(m,"ELIM_VALENCE","valence of elim. species",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrolyte solution (ehrl 07/12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_diffcond",
                                            "material parameters for ion species in electrolyte solution",
                                            INPAR::MAT::m_diffcond));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"VALENCE","valence (= charge number)");
    AddNamedReal(m,"TRANSFERENCE","transference number");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrolyte solution (ehrl 07/12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_newman",
                                            "material parameters for ion species in electrolyte solution",
                                            INPAR::MAT::m_newman));

    AddNamedReal(m,"VALENCE","valence (= charge number)");
    AddNamedInt(m,"CURVE_DIFF","curve number for kinematic diffusivity");
    AddNamedInt(m,"CURVE_TRANS","curve number for transference number");
    AddNamedReal(m,"A","constant for diffusion potential in current equation");
    AddNamedReal(m,"B","constant for diffusion potential in current equation");
    AddNamedReal(m,"C","constant for diffusion potential in current equation");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (gjb 07/08)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_matlist",
                                            "list/collection of materials, i.e. material IDs",
                                            INPAR::MAT::m_matlist));

    AddNamedBool(m,"LOCAL","individual materials allocated per element or only at global scope");
    //AddNamedInt(m,"LOCAL","individual materials allocated per element or only at global scope");
    AddNamedInt(m,"NUMMAT","number of materials in list");
    AddNamedIntVector(m,"MATIDS","the list material IDs","NUMMAT");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (ehrl 11/12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_elchmat",
                                            "specific list/collection of species and phases for elch applications",
                                            INPAR::MAT::m_elchmat));

    AddNamedBool(m,"CURRENT","current flow as a solution variable");
    AddNamedInt(m,"NUMSPEC","number of ionic species in electrolyte");
    AddNamedIntVector(m,"SPECIDS","the list material IDs","NUMSPEC");
    AddNamedInt(m,"NUMPHASE","number of phases in electrolyte");
    AddNamedIntVector(m,"PHASEIDS","the list phasel IDs","NUMPHASE");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (ehrl 11/12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_elchphase",
                                            "material parameters for ion species in electrolyte solution",
                                            INPAR::MAT::m_elchphase));

    AddNamedReal(m,"EPSILON","porousity of the phase");
    AddNamedReal(m,"TORTUOSITY","porousity of the phase");
    AddNamedReal(m,"CONDUCTIVITY","conductivity");
    AddNamedInt(m,"NR","conductivity depending on concentration: number of curve",0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // de St.Venant--Kirchhoff
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_StVenantKirchhoff",
                                            "de St.Venant--Kirchhoff material",
                                            INPAR::MAT::m_stvenant));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"THEXPANS","coefficient of linear thermal expansion",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // de St.Venant--Kirchhoff with temperature
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_ThrStVenantK",
                                            "Thermo St.Venant--Kirchhoff material",
                                            INPAR::MAT::m_thermostvenant));

    AddNamedInt(m,"YOUNGNUM","number of Young's modulus in list");
    AddNamedRealVector(m,"YOUNG","Young's modulus","YOUNGNUM");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"THEXPANS","coefficient of linear thermal expansion");
    AddNamedReal(m,"CAPA","capacity");
    AddNamedReal(m,"CONDUCT","conductivity");
    AddNamedReal(m,"INITTEMP","initial temperature");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Linear thermo-elastic St.Venant Kirchhoff / plastic von Mises
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_ThrPlasticLinElast",
                                            "Thermo-elastic St.Venant Kirchhoff / plastic von Mises material",
                                            INPAR::MAT::m_thermopllinelast));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"THEXPANS","coefficient of linear thermal expansion");
    AddNamedReal(m,"INITTEMP","initial temperature");
    AddNamedReal(m,"YIELD","yield stress");
    AddNamedReal(m,"ISOHARD","isotropic hardening modulus");
    AddNamedReal(m,"KINHARD","kinematic hardening modulus");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Thermo-hyperelasticity / finite strain von-Mises plasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_ThrPlasticHyperElast",
                                            "Thermo-hyperelastic / finite strain plastic von Mises material",
                                            INPAR::MAT::m_thermoplhyperelast));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"YIELD","yield stress");
    AddNamedReal(m,"ISOHARD","isotropic hardening modulus");
    AddNamedReal(m,"KINHARD","kinematic hardening modulus");

    AppendMaterialDefinition(matlist,m);
  }


  /*----------------------------------------------------------------------*/
  // Plastic Neo-Hooke / von Mises
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_PlasticNeoHooke",
                                            "elastic neo-Hooke / plastic von Mises material",
                                            INPAR::MAT::m_plneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"ISOHARD","isotropic hardening");
    AddNamedReal(m,"YIELD","yield stress");
    AddNamedReal(m,"INFYIELD","inf yield stress for nonlinear isotropic hardening");
    AddNamedReal(m,"EXP","exponent for nonlinear isotropic hardening");
    AddNamedReal(m,"KINHARD","kinematic hardening");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Hyperelasticity / finite strain von-Mises plasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_PlasticNlnLogNeoHooke",
                                            "hyperelastic / finite strain plastic von Mises material",
                                            INPAR::MAT::m_plnlnlogneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"YIELD","yield stress");
    AddNamedReal(m,"ISOHARD","isotropic hardening modulus");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Plastic linear elastic St.Venant Kirchhoff / von Mises
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_PlasticLinElast",
                                            "elastic St.Venant Kirchhoff / plastic von Mises material",
                                            INPAR::MAT::m_pllinelast));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"YIELD","yield stress");
    AddNamedReal(m,"ISOHARD","isotropic hardening modulus");
    AddNamedReal(m,"KINHARD","kinematic hardening modulus");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist,m);
  }


  /*----------------------------------------------------------------------*/
  // Robinson's visco-plastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Robinson",
                                            "Robinson's visco-plastic material",
                                            INPAR::MAT::m_vp_robinson));

    AddNamedString(m,"KIND","kind of Robinson material","arya_narloyz");
    AddNamedInt(m,"YOUNGNUM","number of Young's modulus in list");
    AddNamedRealVector(m,"YOUNG","Young's modulus","YOUNGNUM");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"THEXPANS","coefficient of linear thermal expansion");
    AddNamedReal(m,"INITTEMP","initial temperature");
    AddNamedReal(m,"HRDN_FACT","hardening factor 'A'");
    AddNamedReal(m,"HRDN_EXPO","hardening power 'n'");
    AddNamedInt(m,"SHRTHRSHLDNUM","number of shear stress threshold 'K^2'in list");
    AddNamedRealVector(m,"SHRTHRSHLD","Bingam-Prager shear stress threshold 'K^2'", "SHRTHRSHLDNUM");
    AddNamedReal(m,"RCVRY","recovery factor 'R_0'");
    AddNamedReal(m,"ACTV_ERGY","activation energy 'Q_0'");
    AddNamedReal(m,"ACTV_TMPR","activation temperature 'T_0'");
    AddNamedReal(m,"G0","'G_0'");
    AddNamedReal(m,"M_EXPO","'m'");
    AddNamedInt(m,"BETANUM","number of 'beta' in list");
    AddNamedRealVector(m,"BETA","beta", "BETANUM");
    AddNamedReal(m,"H_FACT","'H'");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Elasto-plastic material with damage, based on MAT_Struct_PlasticLinElast
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Damage",
                                            "elasto-plastic von Mises material with damage",
                                            INPAR::MAT::m_elpldamage));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"YIELD","yield stress");
    AddNamedReal(m,"KINHARD","kinematic hardening modulus");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // aneurysm wall material according to Raghavan and Vorp [2000]
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAANeoHooke",
                                            "aneurysm wall material according to Raghavan and Vorp [2000]",
                                            INPAR::MAT::m_aaaneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"BETA","2nd parameter");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }
  /*--------------------------------------------------------------------*/
    // aneurysm wall material according to Raghavan and Vorp [2000]
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAANeoHookeStopro",
                                              "aneurysm wall material according to Raghavan and Vorp [2000] with stochastic modelling of beta",
                                              INPAR::MAT::m_aaaneohooke_stopro));

      AddNamedReal(m,"YOUNG","Young's modulus");
      AddNamedReal(m,"BETA","2nd parameter");
      AddNamedReal(m,"NUE","Poisson's ratio");
      AddNamedReal(m,"DENS","mass density");
      // Stochastic properties are set via randomfield class

      AppendMaterialDefinition(matlist,m);
    }

  /*--------------------------------------------------------------------*/
  // AAA thrombus material according to GASSER et. al. [2008]
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAAGasser",
                                            "AAA thrombus material according to GASSER [2008]",
                                            INPAR::MAT::m_aaagasser));

    AddNamedReal(m,"DENS","mass density");
    AddNamedString(m,"VOL","Type of volumetric Strain Energy Density (OSM,SuBa,SiTa)","OSM");
    AddNamedReal(m,"NUE","Poisson's ratio (0.49)");
    AddNamedReal(m,"BETA","empiric constant for OSM (-2.0)");
    AddNamedReal(m,"CLUM","luminal stiffness parameter (2.62e3)");
    AddNamedReal(m,"CMED","medial stiffness parameter (2.62e3)");
    AddNamedReal(m,"CABLUM","abluminal stiffness parameter (2.62e3)");

    /*
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"KAPPA","dilatation modulus");
    AddNamedReal(m,"BETA","empiric constant");
    AddNamedReal(m,"CLUM","luminal stiffness parameter");
    AddNamedReal(m,"CMED","medial stiffness parameter");
    AddNamedReal(m,"CABLUM","abluminal stiffness parameter");
    */

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // aneurysm wall material according to Raghavan and Vorp [2000] with damage Simo
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Raghavan_Damage",
                                            "aneurysm wall material according to Raghavan and Vorp [2000] with damage",
                                            INPAR::MAT::m_aaaraghavanvorp_damage));

    AddNamedReal(m,"BULK","Bulk's modulus");
    AddNamedReal(m,"ALPHA","1nd parameter,alpha");
    AddNamedReal(m,"BETA","2nd parameter,beta");
    AddNamedReal(m,"EQSTRMIN","equivalent strain initial damage");
    AddNamedReal(m,"A","1st parameter, a");
    AddNamedReal(m,"B","2nd parameter, b");
    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }


  /*--------------------------------------------------------------------*/
  // aneurysm wall material SEF according  to Raghavan and Vorp [2000],
  // parameters according to mixed effects model
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAA_MixedEffects",
                                            "aneurysm wall material according to Mixed Effects Model",
                                            INPAR::MAT::m_aaa_mixedeffects));

    AddNamedReal(m,"AGE","age");
    AddNamedReal(m,"REFDIA","subrenal diameter");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Visco-elastic Neo-Hookean material law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_VISCONEOHOOKE",
                                            "visco-elastic neo-Hookean material law",
                                            INPAR::MAT::m_visconeohooke));
    AddNamedReal(m,"YOUNGS_SLOW","???");
    AddNamedReal(m,"POISSON","???");
    AddNamedReal(m,"DENS","???");
    AddNamedReal(m,"YOUNGS_FAST","???");
    AddNamedReal(m,"RELAX","???");
    AddNamedReal(m,"THETA","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Visco-elastic anisotropic fiber material law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_VISCOANISO",
                                            "visco-elastic anisotropic fibre material law",
                                            INPAR::MAT::m_viscoanisotropic));

    AddNamedReal(m,"KAPPA","dilatation modulus");
    AddNamedReal(m,"MUE","Shear Modulus");
    AddNamedReal(m,"DENS","Density");
    AddNamedReal(m,"K1","Parameter for linear fiber stiffness");
    AddNamedReal(m,"K2","Parameter for exponetial fiber stiffness");
    AddNamedReal(m,"GAMMA","angle between fibers");
    AddNamedReal(m,"BETA_ISO","ratio between elasticities in generalized Maxweel body");
    AddNamedReal(m,"BETA_ANISO","ratio between elasticities in generalized Maxweel body");
    AddNamedReal(m,"RELAX_ISO","isotropic relaxation time");
    AddNamedReal(m,"RELAX_ANISO","anisotropic relaxation time");
    AddNamedReal(m,"MINSTRETCH","minimal principal stretch fibers do respond to");
    AddNamedInt(m,"ELETHICKDIR","Element thickness direction applies also to fibers (only sosh)");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Structural micro-scale approach: material parameters are calculated from microscale simulation
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Multiscale",
                                            "Structural micro-scale approach: material parameters are calculated from microscale simulation",
                                            INPAR::MAT::m_struct_multiscale));

    AddNamedString(m,"MICROFILE","inputfile for microstructure","filename.dat");
    AddNamedInt(m,"MICRODIS_NUM","Number of microscale discretization");
    AddNamedReal(m,"INITVOL","Initial volume of RVE",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ElastHyper",
                                            "list/collection of hyperelastic materials, i.e. material IDs",
                                            INPAR::MAT::m_elasthyper));

    AddNamedInt(m,"NUMMAT","number of materials/potentials in list");
    AddNamedIntVector(m,"MATIDS","the list material/potential IDs","NUMMAT");
    AddNamedReal(m,"DENS","material mass density");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Generalized Maxwell Model compatible with elasthyper
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ViscoGenMax",
                                            "Generalized Maxwell model compatible with the collection of hyperelastic materials",
                                            INPAR::MAT::m_viscogenmax));

    AddNamedInt(m,"NUMMAT","number of materials/potentials in list");
    AddNamedIntVector(m,"MATIDS","the list material/potential IDs","NUMMAT");
    AddNamedReal(m,"DENS","material mass density");
    AddNamedReal(m,"GAMMA","fiber angle");
    AddNamedReal(m,"RELAX_ISOT_PRINC","relaxation time - isotropic not splitted formulation - Leave it to 0 if you don't use it");
    AddNamedReal(m,"BETA_ISOT_PRINC","viscous constant of the generalized maxwell model - isotropic not splitted formulation - Leave it to 0 if you don't use it");
    AddNamedReal(m,"RELAX_ISOT_MOD_VOL","relaxation time - isotropic splitted formulation - volumetric contribution - Leave it to 0 if you don't use it");
    AddNamedReal(m,"BETA_ISOT_MOD_VOL","viscous constant of the generalized maxwell model - isotropic splitted formulation - volumetric contribution - Leave it to 0 if you don't use it");
    AddNamedReal(m,"RELAX_ISOT_MOD_ISOC","relaxation time - isotropic splitted formulation - isochoric contribution - Leave it to 0 if you don't use it");
    AddNamedReal(m,"BETA_ISOT_MOD_ISOC","viscous constant of the generalized maxwell model - isotropic splitted formulation - isochoric contribution - Leave it to 0 if you don't use it");
    AddNamedReal(m,"RELAX_ANISOT_PRINC","relaxation time - anisotropic not splitted formulation - Leave it to 0 if you don't use it");
    AddNamedReal(m,"BETA_ANISOT_PRINC","viscous constant of the generalized maxwell model - anisotropic not splitted formulation - Leave it to 0 if you don't use it");

    // optional

    AddNamedInt(m,"INIT_MODE","initialization modus for fiber alignement", -1,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // logarithmic neo-Hooke material acc. to Bonet and Wood
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupLogNeoHooke",
                                            "logarithmic neo-Hooke material acc. to Bonet and Wood",
                                            INPAR::MAT::mes_couplogneohooke));

    AddNamedString(m,"MODE","parameter set: YN (Young's modulus and Poisson's ration) or Lame (mue and lambda)", "YN");
    AddNamedReal(m,"C1","E or mue");
    AddNamedReal(m,"C2","nue or lambda");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/

    //compressible neo-Hooke material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupNeoHooke",
                                            "compressible neo-Hooke material acc. to Holzapfel",
                                            INPAR::MAT::mes_coupneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus",true);
    AddNamedReal(m,"NUE","Poisson's ratio",true);

    AppendMaterialDefinition(matlist,m);
  }
    //Mooney Rivlin  material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupMooneyRivlin",
                                            "Mooney - Rivlin material acc. to Holzapfel",
                                            INPAR::MAT::mes_coupmooneyrivlin));

    AddNamedReal(m,"C1","material constant",true);
    AddNamedReal(m,"C2","material constant",true);
    AddNamedReal(m,"C3","material constant",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled Blatz and Ko material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupBlatzKo",
                                            "Blatz and Ko material acc. to Holzapfel",
                                            INPAR::MAT::mes_coupblatzko));

    AddNamedReal(m,"MUE","Shear modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"F","interpolation parameter");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of Neo-Hooke
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoNeoHooke",
                                            "isochoric part of  neo-Hooke material acc. to Holzapfel",
                                            INPAR::MAT::mes_isoneohooke));

    AddNamedReal(m,"MUE","Shear modulus");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric and volumetric contribution of HU dependent NeoHooke
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoVolHUDependentNeoHooke",
                                            "isochoric and volumetric part of HU dependent neo-Hooke material",
                                            INPAR::MAT::mes_isovolHUdependentneohooke));

    AddNamedReal(m,"ALPHA_MAX","");
    AddNamedReal(m,"CT_MIN","");
    AddNamedReal(m,"CT_MAX","");
    AddNamedReal(m,"NUE","");
    AddNamedReal(m,"BETA","");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric and volumetric contribution of AAAGasser
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoVolAAAGasser",
                                            "isochoric and volumetric part of AAAGasser material (thrombus)",
                                            INPAR::MAT::mes_isovolaaagasser));

    AddNamedReal(m,"CLUM","luminal stiffness parameter (2.62e3)");
    AddNamedReal(m,"CMED","medial stiffness parameter (2.62e3)");
    AddNamedReal(m,"CABLUM","abluminal stiffness parameter (2.62e3)");
    AddNamedReal(m,"NUE","");
    AddNamedReal(m,"BETA","");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of Yeoh
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoYeoh",
                                            "isochoric part of  Yeoh material acc. to Holzapfel",
                                            INPAR::MAT::mes_isoyeoh));

    AddNamedReal(m,"C1","Linear modulus");
    AddNamedReal(m,"C2","Quadratic modulus");
    AddNamedReal(m,"C3","Cubic modulus");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of Quad
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoQuad",
                                            "isochoric part of quadratic material",
                                            INPAR::MAT::mes_isoquad));

    AddNamedReal(m,"C","material parameter");

    AppendMaterialDefinition(matlist,m);
  }
  /*--------------------------------------------------------------------*/



  // isochoric contribution of Cub
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoCub",
                                            "isochoric part of cubic material",
                                            INPAR::MAT::mes_isocub));

    AddNamedReal(m,"C","material parameter");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/

  // isochoric contribution of iso1pow
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_Iso1Pow",
                                            "isochoric part of general power material",
                                            INPAR::MAT::mes_iso1pow));

    AddNamedReal(m,"C","material parameter");
    AddNamedInt(m,"D","exponent");
    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of iso2pow
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_Iso2Pow",
                                            "isochoric part of general power material",
                                            INPAR::MAT::mes_iso2pow));

    AddNamedReal(m,"C","material parameter");
    AddNamedInt(m,"D","exponent");
    AppendMaterialDefinition(matlist,m);
  }

    /*--------------------------------------------------------------------*/

  // contribution of coup1pow
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_Coup1Pow",
                                            "part of general power material",
                                            INPAR::MAT::mes_coup1pow));

    AddNamedReal(m,"C","material parameter");
    AddNamedInt(m,"D","exponent");
    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // contribution of iso2pow
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_Coup2Pow",
                                            "part of general power material",
                                            INPAR::MAT::mes_coup2pow));

    AddNamedReal(m,"C","material parameter");
    AddNamedInt(m,"D","exponent");
    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/



  // isochoric contribution of expo
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoExpoPow",
                                            "isochoric part of  exponential material acc. to Holzapfel",
                                            INPAR::MAT::mes_isoexpopow));

    AddNamedReal(m,"K1","material parameter");
    AddNamedReal(m,"K2","material parameter");
    AddNamedInt(m,"C","exponent");
    AppendMaterialDefinition(matlist,m);
  }



  /*--------------------------------------------------------------------*/
  // isochoric contribution of mooney rivlin
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoMooneyRivlin",
                                            "isochoric part of  Mooney-Rivlin material acc. to Holzapfel",
                                            INPAR::MAT::mes_isomooneyrivlin));

    AddNamedReal(m,"C1","Linear modulus for first invariant");
    AddNamedReal(m,"C2","Linear modulus for second invariant");
    AppendMaterialDefinition(matlist,m);
  }


  /*--------------------------------------------------------------------*/
  // volumetric contribution of Sussman Bathe
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_VolSussmanBathe",
                                            "volumetric part of  SussmanBathe material",
                                            INPAR::MAT::mes_volsussmanbathe));

    AddNamedReal(m,"KAPPA","dilatation modulus");

    AppendMaterialDefinition(matlist,m);
  }


  /*--------------------------------------------------------------------*/
  // volumetric penalty contribution
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_VolPenalty",
                                            "Penalty formulation for the volumetric part",
                                            INPAR::MAT::mes_volpenalty));

    AddNamedReal(m,"EPSILON","penalty parameter");
    AddNamedReal(m,"GAMMA","penalty parameter");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // volumetric contribution of Ogden
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_VolOgden",
                                            "Ogden formulation for the volumetric part",
                                            INPAR::MAT::mes_vologden));

    AddNamedReal(m,"KAPPA","dilatation modulus");
    AddNamedReal(m,"BETA","empiric constant");

    AppendMaterialDefinition(matlist,m);
  }


  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with one exponential fiber family
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpo",
                                            "anisotropic part with one exp. fiber",
                                            INPAR::MAT::mes_coupanisoexpo));

    AddNamedReal(m,"K1","linear constant");
    AddNamedReal(m,"K2","exponential constant");
    AddNamedReal(m,"GAMMA","angle");
    AddNamedReal(m,"K1COMP","linear constant");
    AddNamedReal(m,"K2COMP","exponential constant");
    AddNamedInt(m,"INIT","initialization modus for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with two exponential fiber families
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpoTwoCoup",
                                            "anisotropic part with two exp. fibers",
                                            INPAR::MAT::mes_coupanisoexpotwocoup));

    AddNamedReal(m,"A4","linear anisotropic constant for fiber 1");
    AddNamedReal(m,"B4","exponential anisotropic constant for fiber 1");
    AddNamedReal(m,"A6","linear anisotropic constant for fiber 2");
    AddNamedReal(m,"B6","exponential anisotropic constant for fiber 2");
    AddNamedReal(m,"A8","linear anisotropic constant for fiber 1 relating fiber 2");
    AddNamedReal(m,"B8","exponential anisotropic constant for fiber 1 relating fiber 2");
    AddNamedReal(m,"GAMMA","angle");
    AddNamedInt(m,"INIT","initialization modus for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with two exponential fiber families
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoNeoHooke",
                                            "anisotropic part with one neo Hookean fiber",
                                            INPAR::MAT::mes_coupanisoneohooke));

    AddNamedReal(m,"C","linear constant");
    AddNamedReal(m,"GAMMA","angle");
    AddNamedInt(m,"INIT","initialization modus for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist,m);
  }

   /*--------------------------------------------------------------------*/
  // coupled anisotropic material with the stress given by a the contraction law of Bestel-Clement-Sorine
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoNeoHooke_ActiveStress",
                                            "anisotropic part with one neo Hookean fiber with coefficient given by a simplification of the activation-contraction law of Bestel-Clement-Sorine-2001",
                                            INPAR::MAT::mes_coupanisoneohooke_activestress));

    AddNamedReal(m,"SIGMA","Contractility (maximal stress)");
    AddNamedReal(m,"TAUC0","Initial value for the active stress");
    AddNamedReal(m,"MAX_ACTIVATION","Maximal value for the rescaled activation");
    AddNamedReal(m,"MIN_ACTIVATION","Minimal value for the rescaled activation");
    AddNamedInt(m,"SOURCE_ACTIVATION","Where the activation comes from: 0=scatra , >0 Id for FUNCT");
    AddNamedReal(m,"GAMMA","azimuth angle", true);
    AddNamedReal(m,"THETA","polar angle", true);
    AddNamedInt(m,"INIT","initialization mode for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist,m);
  }


  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with variable stress coefficient
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoNeoHooke_VarProp",
                                            "anisotropic part with one neo Hookean fiber with variable coefficient",
                                            INPAR::MAT::mes_coupanisoneohooke_varprop));

    AddNamedReal(m,"C","linear constant");
    AddNamedInt(m,"SOURCE_ACTIVATION","Where the activation comes from: 0=scatra , >0 Id for FUNCT");
    AddNamedReal(m,"GAMMA","azimuth angle", true);
    AddNamedReal(m,"THETA","polar angle", true);
    AddNamedInt(m,"INIT","initialization mode for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist,m);
  }


  /*--------------------------------------------------------------------*/
  // isochoric anisotropic material with one exponential fiber family
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoAnisoExpo",
                                            "anisotropic part with one exp. fiber",
                                            INPAR::MAT::mes_isoanisoexpo));

    AddNamedReal(m,"K1","linear constant");
    AddNamedReal(m,"K2","exponential constant");
    AddNamedReal(m,"GAMMA","angle");
    AddNamedReal(m,"K1COMP","linear constant");
    AddNamedReal(m,"K2COMP","exponential constant");
    AddNamedInt(m,"INIT","initialization modus for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled Varga material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupVarga",
                                            "Varga material acc. to Holzapfel",
                                            INPAR::MAT::mes_coupvarga));

    AddNamedReal(m,"MUE","Shear modulus");
    AddNamedReal(m,"BETA","'Anti-modulus'");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric Varga material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoVarga",
                                            "Isochoric Varga material acc. to Holzapfel",
                                            INPAR::MAT::mes_isovarga));

    AddNamedReal(m,"MUE","Shear modulus");
    AddNamedReal(m,"BETA","'Anti-modulus'");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  // 1D Artery material with constant properties
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_CNST_ART",
                                            "artery with constant properties",
                                            INPAR::MAT::m_cnst_art));

    AddNamedReal(m,"VISCOSITY","viscosity of blood");
    AddNamedReal(m,"DENS","density of blood");
    AddNamedReal(m,"YOUNG","artery Youngs modulus of elasticity");
    AddNamedReal(m,"NUE","Poissons ratio of artery fiber");
    AddNamedReal(m,"DIAM","artery initial diameter");
    AddNamedReal(m,"TH","artery thickness");
    AddNamedReal(m,"PEXT1","artery fixed external pressure 1");
    AddNamedReal(m,"PEXT2","artery fixed external pressure 2");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  // Fourier's law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("THERM_FourierIso",
                                            "isotropic (linear) Fourier's law of heat conduction",
                                            INPAR::MAT::m_th_fourier_iso));

    AddNamedReal(m,"CAPA","capacity");
    AddNamedReal(m,"CONDUCT","conductivity");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // integration point based growth
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GROWTH",
                                            "integration point based growth",
                                            INPAR::MAT::m_growth));

    AddNamedReal(m,"DENS","Density");
    AddNamedInt(m,"IDMATELASTIC","number of elastic material in input file: MAT IDMATELASTIC ...");
    AddNamedReal(m,"STARTTIME","start growth after this time");
    AddNamedReal(m,"ENDTIME","end growth after this time");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");
    AddNamedReal(m,"KPLUS","growth law parameter kthetaplus");
    AddNamedReal(m,"MPLUS","growth law parameter mthetaplus");
    AddNamedReal(m,"KMINUS","growth law parameter kthetaminus");
    AddNamedReal(m,"MMINUS","growth law parameter mthetaminus");
    AddNamedReal(m,"HOMMANDEL","homeostatic value for mandelstress");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // growth and remodeling of arteries
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ConstraintMixture",
                                            "growth and remodeling of arteries",
                                            INPAR::MAT::m_constraintmixture));

    AddNamedReal(m,"DENS","Density");
    AddNamedReal(m,"MUE","Shear Modulus");
    AddNamedReal(m,"PHIE","mass fraction of elastin");
    AddNamedReal(m,"PREELA","prestretch of elastin");
    AddNamedReal(m,"K1","Parameter for linear collagen fiber stiffness");
    AddNamedReal(m,"K2","Parameter for exponential collagen fiber stiffness");
    AddNamedReal(m,"PRECOLL","prestretch of collagen fibers");
    AddNamedReal(m,"DAMAGE","damage stretch of collagen fibers");
    AddNamedReal(m,"K1M","Parameter for linear smooth muscle fiber stiffness");
    AddNamedReal(m,"K2M","Parameter for exponential smooth muscle fiber stiffness");
    AddNamedReal(m,"PHIM","mass fraction of smooth muscle");
    AddNamedReal(m,"PREMUS","prestretch of smooth muscle fibers");
    AddNamedReal(m,"SMAX","maximal active stress");
    AddNamedReal(m,"KAPPA","dilatation modulus");
    AddNamedReal(m,"LIFETIME","lifetime of collagen fibers");
    AddNamedReal(m,"HOMSTR","homeostatic target value of scalar stress measure");
    AddNamedReal(m,"GROWTHFAC","growth factor");
    AddNamedReal(m,"STARTTIME","at this time turnover of collagen starts");
    AddNamedString(m,"INTEGRATION","time integration scheme (Explicit, Implicit)","Explicit");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");
    AddNamedString(m,"GROWTHFORCE","driving force of growth (Single, All)","Single");
    AddNamedString(m,"INITSTRETCH","how to set stretches in the beginning (None, Homeo)","None");
    AddNamedInt(m,"CURVE","number of timecurve for increase of prestretch in time",0);
    AddNamedString(m,"DEGOPTION","which degradation function (Lin, Cos, Exp)","Lin");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // biofilm modeling (convection-diffusion-reaction equation)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_biofilm",
                                            "biofilm material",
                                            INPAR::MAT::m_biofilm));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"REARATE","substrate uptake rate coefficient");
    AddNamedReal(m,"SATCOEFF","substrate saturation coefficient");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // optimization modeling
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_opti",
                                            "optimization material",
                                            INPAR::MAT::m_opti_dens));

    AddNamedReal(m,"MINPORO","minimal porosity");
    AddNamedReal(m,"MAXPORO","maximal porosity");
    AddNamedReal(m,"SMEARFAC","smearing factor");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // hyperelastic material for poroelasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_StructPoro",
                                            "wrapper for structure porelastic material",
                                            INPAR::MAT::m_structporo));

    AddNamedInt(m,"MATID","ID of structure material");
    AddNamedReal(m,"INITPOROSITY","initial porosity of porous medium");
    AddNamedReal(m,"BULKMODULUS","bulk modulus of porous medium");
    AddNamedReal(m,"PENALTYPARAMETER","penalty paramter of porous medium");
  //  AddNamedBool(m,"REACTION","switch for reaction in porous medium");

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // hyperelastic material for poroelasticity with reaction
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_StructPoroReaction",
                                            "wrapper for structure porelastic material with reaction",
                                            INPAR::MAT::m_structpororeaction));

    AddNamedInt(m,"MATID","ID of structure material");
    AddNamedReal(m,"INITPOROSITY","initial porosity of porous medium");
    AddNamedReal(m,"BULKMODULUS","bulk modulus of porous medium");
    AddNamedReal(m,"PENALTYPARAMETER","penalty paramter of porous medium");

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // fluid flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoro",
                                            "flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo));

    AddNamedReal(m,"DYNVISCOSITY","dynamic viscosity");
    AddNamedReal(m,"DENSITY","density");
    AddNamedReal(m,"PERMEABILITY","permeability of medium");
    AddNamedString(m,"TYPE","Problem type: Darcy or Darcy-Brinkman","Darcy");
  //  AddNamedReal(m,"BULKMODULUS","bulk modulus of medium");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // elastic spring
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Spring",
                                            "elastic spring",
                                            INPAR::MAT::m_spring));

    AddNamedReal(m,"STIFFNESS","spring constant");
    AddNamedReal(m,"DENS","density");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // 0D Acinar material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS",
                                            "0D acinar material",
                                            INPAR::MAT::m_0d_maxwell_acinus));

    AddNamedReal(m,"Stiffness1","first stiffness");
    AddNamedReal(m,"Stiffness2","second stiffness");
    AddNamedReal(m,"Viscosity1","first viscosity");
    AddNamedReal(m,"Viscosity2","second viscosity");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*/
  // particle material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Particle",
                                            "particle material",
                                            INPAR::MAT::m_particlemat));

    AddNamedReal(m,"DENSITY","mass density");
    AddNamedReal(m,"INITRADIUS","initial radius of particle");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // deliver
  return vm;
}
