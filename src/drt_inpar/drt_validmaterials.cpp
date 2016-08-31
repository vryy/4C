/*----------------------------------------------------------------------*/
/*!
\file drt_validmaterials.cpp
\maintainer Martin Kronbichler

\brief Setup of the list of valid materials for input

<pre>
\level 1

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
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
    AddNamedReal(m,"GAMMA","surface tension coefficient",0.0,true);

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
  // fluid with non-linear viscosity according to Herschel-Bulkley
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_herschelbulkley",
                                            "fluid with non-linear viscosity according to Herschel-Bulkley",
                                            INPAR::MAT::m_herschelbulkley));

    AddNamedReal(m,"TAU_0","yield stress");
    AddNamedReal(m,"KFAC","constant factor");
    AddNamedReal(m,"NEXP","exponent");
    AddNamedReal(m,"MEXP","exponent");
    AddNamedReal(m,"LOLIMSHEARRATE","lower limit of shear rate");
    AddNamedReal(m,"UPLIMSHEARRATE","upper limit of shear rate");
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
  // cavitation fluid
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_cavitation",
                                            "Cavitation fluid",
                                            INPAR::MAT::m_cavitation));

    AddNamedReal(m,"DYNVISCOSITY","dynamic viscosity");
    AddNamedReal(m,"DENSITY","spatial mass density");
    AddNamedReal(m,"GAMMA","surface tension coefficient");
    AddNamedReal(m,"PVAPOR","vapor pressure");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // lubrication material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_lubrication",
                                            "lubrication material",
                                            INPAR::MAT::m_lubrication));

    AddNamedReal(m,"VISCOSITY","lubricant viscosity");
    AddNamedReal(m,"DENSITY","lubricant density");

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
    AddNamedReal(m,"REACOEFF","reaction coefficient",0.0,true);
    AddNamedReal(m,"SCNUM","schmidt number",0.0,true);
    AddNamedReal(m,"DENSIFICATION","densification coefficient",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }


  /*----------------------------------------------------------------------*/
  // scalar transport material (with potential reaction coefficient)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra_reaction_poro",
                                            "scalar transport material",
                                            INPAR::MAT::m_scatra_reaction_poroECM));

    AddNamedInt(m,"NUMSCAL","number of reactions for these elements");
    AddNamedIntVector(m,"STOICH","advanced reaction list","NUMSCAL");
    AddNamedReal(m,"REACCOEFF","reaction coefficient");
    AddNamedReal(m,"REACSCALE","scaling for reaction coefficient");
    AddNamedString(m,"COUPLING","type of coupling", "no_coupling",false);
    AddNamedRealVector(m,"ROLE","role in michaelis-menten like reactions","NUMSCAL",-1.0,false);
    AddNamedRealVector(m,"REACSTART","starting point of reaction","NUMSCAL",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // scalar transport reaction material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra_reaction",
                                            "advanced reaction material",
                                            INPAR::MAT::m_scatra_reaction));

    AddNamedInt(m,"NUMSCAL","number of reactions for these elements");
    AddNamedIntVector(m,"STOICH","advanced reaction list","NUMSCAL");
    AddNamedReal(m,"REACCOEFF","reaction coefficient");
    AddNamedString(m,"COUPLING","type of coupling", "no_coupling",false);
    AddNamedRealVector(m,"ROLE","role in michaelis-menten like reactions","NUMSCAL",-1.0,false);
    AddNamedRealVector(m,"REACSTART","starting point of reaction","NUMSCAL",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport reaction material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiporo",
                                            "advanced reaction material for multiphase porous flow",
                                            INPAR::MAT::m_scatra_multiporo));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedInt(m,"PHASEID","ID of fluid phase the scalar is associated with");
    AddNamedReal(m,"REACOEFF","reaction coefficient",0.0,true);
    AddNamedReal(m,"SCNUM","schmidt number",0.0,true);
    AddNamedReal(m,"DENSIFICATION","densification coefficient",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport chemotaxis material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra_chemotaxis",
                                            "chemotaxis material",
                                            INPAR::MAT::m_scatra_chemotaxis));

    AddNamedInt(m,"NUMSCAL","number of chemotactic pairs for these elements");
    AddNamedIntVector(m,"PAIR","chemotaxis pairing","NUMSCAL");
    AddNamedReal(m,"CHEMOCOEFF","chemotaxis coefficient");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // anisotropic scalar transport material (with potential reaction coefficient)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra_aniso",
                                            "anisotropic scalar transport material",
                                            INPAR::MAT::m_scatra_aniso));

    AddNamedReal(m,"DIFF1","kinematic diffusivity component 1");
    AddNamedReal(m,"DIFF2","kinematic diffusivity component 2");
    AddNamedReal(m,"DIFF3","kinematic diffusivity component 3");
    AddNamedReal(m,"REACOEFF","reaction coefficient",0.0,true);
//    AddNamedReal(m,"SCNUM","schmidt number",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport material for multi-scale approach
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiscale",
                                            "scalar transport material for multi-scale approach",
                                            INPAR::MAT::m_scatra_multiscale));

    AddNamedString(m,"MICROFILE","input file for micro scale","filename.dat");
    AddNamedInt(m,"MICRODIS_NUM","number of micro-scale discretization");
    AddNamedReal(m,"POROSITY","porosity");
    AddNamedReal(m,"TORTUOSITY","tortuosity");
    AddNamedReal(m,"A_s","specific micro-scale surface area");
    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"REACOEFF","reaction coefficient",0.0,true);
    AddNamedReal(m,"SCNUM","Schmidt number",0.0,true);
    AddNamedReal(m,"DENSIFICATION","densification coefficient",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Myocard muscle material (with complicated reaction coefficient)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_myocard",
                                            "Myocard muscle material",
                                            INPAR::MAT::m_myocard));

    AddNamedReal(m,"DIFF1","conductivity in fiber direction");
    AddNamedReal(m,"DIFF2","conductivity perpendicular to fiber direction");
    AddNamedReal(m,"DIFF3","conductivity perpendicular to fiber direction");
    AddNamedReal(m,"PERTUBATION_DERIV","pertubation for calculation of reaction coefficient derivative");
    AddNamedString(m,"MODEL","Model type: MV, FHN, TNNP, SAN or INADA", "MV");
    AddNamedString(m,"TISSUE","Tissue type: M, ENDO, EPI, AN, N or NH", "M");
    AddNamedReal(m, "TIME_SCALE","Scale factor for time units of Model");
    AddNamedInt(m, "GP","Number of Gauss Points for the evaluation of the material", 1, true);

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
    AddNamedReal(m,"GASCON","specific gas constant R (J/(kg*K))");

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
    AddNamedReal(m,"DENSIFICATION","densification coefficient",0.0,true);
    // via these two optional parameters we can bring the material parameters
    // of one eliminated ionic species into BACI if needed
    AddNamedReal(m,"ELIM_DIFFUSIVITY","kinematic diffusivity of elim. species",0.0,true);
    AddNamedReal(m,"ELIM_VALENCE","valence of elim. species",0.0,true);

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
    AddNamedInt(m,"DIFFCOEF","curve number for diffusion coefficient");
    AddNamedInt(m,"TRANSNR","curve number for transference number");
    AddNamedInt(m,"THERMFAC","curve number for thermodynamic factor");
    AddNamedInt(m,"COND","curve number for conductivity");
    // optional parameter for implemented concentration depending function
    AddNamedInt(m,"DIFF_PARA_NUM","number of parameters for diffusion coefficient",0,true);
    AddNamedRealVector(m,"DIFF_PARA","parameters for diffusion coefficient","DIFF_PARA_NUM",0.0,true);
    AddNamedInt(m,"TRANS_PARA_NUM","number of parameters for transference number",0,true);
    AddNamedRealVector(m,"TRANS_PARA","parameters for transference number","TRANS_PARA_NUM",0.0,true);
    AddNamedInt(m,"THERM_PARA_NUM","number of parameters for thermodynamic factor",0,true);
    AddNamedRealVector(m,"THERM_PARA","parameters for thermodynamic factor","THERM_PARA_NUM",0.0,true);
    AddNamedInt(m,"COND_PARA_NUM","number of parameters for conductivity",0,true);
    AddNamedRealVector(m,"COND_PARA","parameters for conductivity","COND_PARA_NUM",0.0,true);
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // electrode material (fang 02/15)
  {
    Teuchos::RCP<MaterialDefinition> matelectrode = Teuchos::rcp(new MaterialDefinition("MAT_electrode","electrode material",INPAR::MAT::m_electrode));

    // diffusivity and electronic conductivity
    AddNamedInt(matelectrode,"DIFFCOEF","curve number for diffusion coefficient");
    AddNamedInt(matelectrode,"COND","curve number for electronic conductivity");

    // optional parameters for concentration dependency of diffusivity and electronic conductivity
    AddNamedInt(matelectrode,"DIFF_PARA_NUM","number of parameters for diffusion coefficient",0,true);
    AddNamedRealVector(matelectrode,"DIFF_PARA","parameters for diffusion coefficient","DIFF_PARA_NUM",0.0,true);
    AddNamedInt(matelectrode,"COND_PARA_NUM","number of parameters for electronic conductivity",0,true);
    AddNamedRealVector(matelectrode,"COND_PARA","parameters for electronic conductivity","COND_PARA_NUM",0.0,true);

    // saturation value of intercalated Lithium concentration
    AddNamedReal(matelectrode,"C_MAX","saturation value of intercalated Lithium concentration");

    // model for half cell open circuit potential of electrode
    AddNamedString(matelectrode,"OCP_MODEL","model for half cell open circuit potential of electrode","none");

    // number of parameters underlying half cell open circuit potential model
    AddNamedInt(matelectrode,"OCP_PARA_NUM","number of parameters underlying half cell open circuit potential model");

    // parameters underlying half cell open circuit potential model
    AddNamedRealVector(matelectrode,"OCP_PARA","parameters underlying half cell open circuit potential model","OCP_PARA_NUM");

    // end of input line
    AddNamedSeparator(matelectrode,"END","indicating end of line");

    // add electrode material to global list of valid materials
    AppendMaterialDefinition(matlist,matelectrode);
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
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection with reactions (thon 09/14)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_matlist_reactions",
                                            "list/collection of materials, i.e. material IDs and list of reactions",
                                            INPAR::MAT::m_matlist_reactions));

    AddNamedBool(m,"LOCAL","individual materials allocated per element or only at global scope");
    AddNamedInt(m,"NUMMAT","number of materials in list");
    AddNamedIntVector(m,"MATIDS","the list material IDs","NUMMAT");
    AddNamedInt(m,"NUMREAC","number of reactions for these elements",0);
    AddNamedIntVector(m,"REACIDS","advanced reaction list","NUMREAC",0);
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection with chemotaxis (thon 06/15)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_matlist_chemotaxis",
                                            "list/collection of materials, i.e. material IDs and list of chemotactic pairs",
                                            INPAR::MAT::m_matlist_chemotaxis));

    AddNamedBool(m,"LOCAL","individual materials allocated per element or only at global scope");
    AddNamedInt(m,"NUMMAT","number of materials in list");
    AddNamedIntVector(m,"MATIDS","the list material IDs","NUMMAT");
    AddNamedInt(m,"NUMPAIR","number of pairs for these elements",0);
    AddNamedIntVector(m,"PAIRIDS","chemotaxis pairs list","NUMPAIR",0);
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection with reactions AND chemotaxis (thon 06/15)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_matlist_chemo_reac",
                                            "list/collection of materials, i.e. material IDs and list of reactive/chemotactic pairs",
                                            INPAR::MAT::m_matlist_chemoreac));

    AddNamedBool(m,"LOCAL","individual materials allocated per element or only at global scope");
    AddNamedInt(m,"NUMMAT","number of materials in list");
    AddNamedIntVector(m,"MATIDS","the list material IDs","NUMMAT");
    AddNamedInt(m,"NUMPAIR","number of pairs for these elements",0);
    AddNamedIntVector(m,"PAIRIDS","chemotaxis pairs list","NUMPAIR",0);
    AddNamedInt(m,"NUMREAC","number of reactions for these elements",0);
    AddNamedIntVector(m,"REACIDS","advanced reaction list","NUMREAC",0);
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (ehrl 11/12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_elchmat",
                                            "specific list/collection of species and phases for elch applications",
                                            INPAR::MAT::m_elchmat));

    AddNamedInt(m,"NUMDOF","number of dof's per node");
    AddNamedInt(m,"NUMSCAL","number of transported scalars per node");
    AddNamedInt(m,"NUMPHASE","number of phases in electrolyte");
    AddNamedIntVector(m,"PHASEIDS","the list phasel IDs","NUMPHASE");
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (ehrl 11/12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_elchphase",
                                            "material parameters for ion species in electrolyte solution",
                                            INPAR::MAT::m_elchphase));

    AddNamedReal(m,"EPSILON","phase porosity");
    AddNamedReal(m,"TORTUOSITY","inverse (!) of phase tortuosity");
    AddNamedInt(m,"NUMMAT","number of materials in electrolyte");
    AddNamedIntVector(m,"MATIDS","the list phasel IDs","NUMMAT");
    AddNamedSeparator(m,"END","indicating end of line");

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
    AddNamedReal(m,"THEXPANS","coefficient of linear thermal expansion",0.0,true);

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
    AddNamedInt(m,"SAMPLENUM","number of stress-strain pairs in list");
    AddNamedRealVector(m,"SIGMA_Y","yield stress", "SAMPLENUM");
    AddNamedRealVector(m,"EPSBAR_P","accumulated plastic strain corresponding to SIGMA_Y", "SAMPLENUM");
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
    AddNamedReal(m,"CTE","coefficient of thermal expansion",0.,true);
    AddNamedReal(m,"INITTEMP","initial, reference temperature",0.,true);
    AddNamedReal(m,"YIELD","initial yield stress");
    AddNamedReal(m,"ISOHARD","isotropic hardening modulus",0.,true);
    AddNamedReal(m,"SATHARDENING","saturation hardening",0.,true);
    AddNamedReal(m,"HARDEXPO","hardening exponent",0.,true);
    AddNamedReal(m,"YIELDSOFT","yield stress softening",0.,true);
    AddNamedReal(m,"HARDSOFT","hardening softening",0.,true);
    AddNamedReal(m,"TOL","tolerance for local Newton iteration",1.e-8,true);

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
    AddNamedReal(m,"SATHARDENING","saturation hardening");
    AddNamedReal(m,"HARDEXPO","hardening exponent");
    AddNamedReal(m,"VISC","VISCOSITY",0.,true);
    AddNamedReal(m,"RATE_DEPENDENCY","rate dependency",0.,true);

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
                                            "elasto-plastic von Mises material with ductile damage",
                                            INPAR::MAT::m_elpldamage));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedInt(m,"SAMPLENUM","number of stress-strain pairs in list");
    AddNamedRealVector(m,"SIGMA_Y","yield stress", "SAMPLENUM");
    AddNamedRealVector(m,"EPSBAR_P","accumulated plastic strain corresponding to SIGMA_Y", "SAMPLENUM");
    AddNamedReal(m,"DAMDEN","denominator of damage evoluation law");
    AddNamedReal(m,"DAMEXP","exponent of damage evoluation law");
    AddNamedReal(m,"DAMTHRESHOLD","damage threshold");
    AddNamedReal(m,"KINHARD","kinematic hardening modulus, stress-like variable");
    AddNamedReal(m,"KINHARD_REC","recovery factor, scalar-valued variable");
    AddNamedReal(m,"SATHARDENING","saturation hardening");
    AddNamedReal(m,"HARDEXPO","hardening exponent");
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
    AddNamedReal(m,"INITVOL","Initial volume of RVE",0.0,true);

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
    AddNamedInt(m,"POLYCONVEX","1.0 if polyconvexity of system is checked",0.,true);

    AppendMaterialDefinition(matlist,m);
  }

   /*----------------------------------------------------------------------*/
  // viscohyperelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ViscoElastHyper",
                                            "Viscohyperelastic material compatible with the collection of hyperelastic materials",
                                            INPAR::MAT::m_viscoelasthyper));

    AddNamedInt(m,"NUMMAT","number of materials/potentials in list");
    AddNamedIntVector(m,"MATIDS","the list material/potential IDs","NUMMAT");
    AddNamedReal(m,"DENS","material mass density");
    AddNamedInt(m,"POLYCONVEX","1.0 if polyconvexity of system is checked",0.,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials for finite strain plasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PlasticElastHyper",
                                            "list/collection of hyperelastic materials, i.e. material IDs",
                                            INPAR::MAT::m_plelasthyper));

    AddNamedInt(m,"NUMMAT","number of materials/potentials in list");
    AddNamedIntVector(m,"MATIDS","the list material/potential IDs","NUMMAT");
    AddNamedReal(m,"DENS","material mass density");
    AddNamedReal(m,"INITYIELD","initial yield stress");
    AddNamedInt(m,"POLYCONVEX","1.0 if polyconvexity of system is checked",0.,true);
    AddNamedReal(m,"ISOHARD","linear isotropic hardening modulus",0.,true);
    AddNamedReal(m,"EXPISOHARD","nonlinear isotropic hardening exponent",0.,true);
    AddNamedReal(m,"INFYIELD","saturation yield stress for nonlinear isotropic hardening",0.,true);
    AddNamedReal(m,"KINHARD","linear kinematic hardening modulus",0.,true);

    // visco-plasticity
    AddNamedReal(m,"VISC", "Visco-Plasticity parameter 'eta' in Perzyna model",0.,true);
    AddNamedReal(m,"RATE_DEPENDENCY", "Visco-Plasticity parameter 'eta' in Perzyna model",1.,true);
    AddNamedReal(m,"VISC_SOFT", "Visco-Plasticity temperature dependency (eta = eta_0 * (1-(T-T_0)*x)",0.,true);

    // optional pastic spin parameter
    AddNamedReal(m,"PL_SPIN_CHI","Plastic spin coupling parameter chi (often called eta)",0.0,true);

    // optional Hill yield parameters
    AddNamedReal(m,"rY_11","relative yield stress in fiber1-direction (Y_11/Y_0)",0.0,true);
    AddNamedReal(m,"rY_22","relative yield stress in fiber2-direction (Y_22/Y_0)",0.0,true);
    AddNamedReal(m,"rY_33","relative yield stress in fiber3-direction (Y_33/Y_0)",0.0,true);
    AddNamedReal(m,"rY_12","relative shear yield stress in 12-direction (Y_12/Y_0)",0.0,true);
    AddNamedReal(m,"rY_23","relative shear yield stress in 23-direction (Y_23/Y_0)",0.0,true);
    AddNamedReal(m,"rY_13","relative shear yield stress in 13-direction (Y_13/Y_0)",0.0,true);

    // optional TSI parameters
    AddNamedReal(m,"CTE","coefficient of thermal expansion",0.,true);
    AddNamedReal(m,"INITTEMP","initial, reference temperature",0.,true);
    AddNamedReal(m,"YIELDSOFT","yield stress softening",0.,true);
    AddNamedReal(m,"HARDSOFT","hardening softening",0.,true);
    AddNamedReal(m,"TAYLOR_QUINNEY","Taylor-Quinney factor for plastic heat conversion",1.,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials for finite strain plasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PlasticElastHyperVCU",
                                            "list/collection of hyperelastic materials, i.e. material IDs",
                                            INPAR::MAT::m_plelasthyperVCU));

    AddNamedInt(m,"NUMMAT","number of materials/potentials in list");
    AddNamedIntVector(m,"MATIDS","the list material/potential IDs","NUMMAT");
    AddNamedReal(m,"DENS","material mass density");
    AddNamedReal(m,"INITYIELD","initial yield stress");
    AddNamedReal(m,"ISOHARD","linear isotropic hardening modulus",0.,true);
    AddNamedReal(m,"EXPISOHARD","nonlinear isotropic hardening exponent",0.,true);
    AddNamedReal(m,"INFYIELD","saturation yield stress for nonlinear isotropic hardening",0.,true);
    AddNamedReal(m,"KINHARD","linear kinematic hardening modulus",0.,true);

    // visco-plasticity
    AddNamedReal(m,"VISC", "Visco-Plasticity parameter 'eta' in Perzyna model",0.,true);
    AddNamedReal(m,"RATE_DEPENDENCY", "Visco-Plasticity parameter 'eta' in Perzyna model",1.,true);
    AddNamedReal(m,"VISC_SOFT", "Visco-Plasticity temperature dependency (eta = eta_0 * (1-(T-T_0)*x)",0.,true);

    // optional pastic spin parameter
    AddNamedReal(m,"PL_SPIN_CHI","Plastic spin coupling parameter chi (often called eta)",0.0,true);

    // optional Hill yield parameters
    AddNamedReal(m,"rY_11","relative yield stress in fiber1-direction (Y_11/Y_0)",0.0,true);
    AddNamedReal(m,"rY_22","relative yield stress in fiber2-direction (Y_22/Y_0)",0.0,true);
    AddNamedReal(m,"rY_33","relative yield stress in fiber3-direction (Y_33/Y_0)",0.0,true);
    AddNamedReal(m,"rY_12","relative shear yield stress in 12-direction (Y_12/Y_0)",0.0,true);
    AddNamedReal(m,"rY_23","relative shear yield stress in 23-direction (Y_23/Y_0)",0.0,true);
    AddNamedReal(m,"rY_13","relative shear yield stress in 13-direction (Y_13/Y_0)",0.0,true);

    // optional TSI parameters
    AddNamedReal(m,"CTE","coefficient of thermal expansion",0.,true);
    AddNamedReal(m,"INITTEMP","initial, reference temperature",0.,true);
    AddNamedReal(m,"YIELDSOFT","yield stress softening",0.,true);
    AddNamedReal(m,"HARDSOFT","hardening softening",0.,true);
    AddNamedReal(m,"TAYLOR_QUINNEY","Taylor-Quinney factor for plastic heat conversion",1.,true);

    AddNamedInt(m,"POLYCONVEX","1.0 if polyconvexity of system is checked",0.,true);


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
  // Saint-Venant-Kirchhoff as elastic summand
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupSVK",
                                            "Saint-Venant-Kirchhoff as elastic summand",
                                            INPAR::MAT::mes_coupSVK));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // Simo-Pister type material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupSimoPister",
                                            "Simo-Pister type material",
                                            INPAR::MAT::mes_coupsimopister));

    AddNamedReal(m,"MUE","material constant");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // logarithmic mixed neo-Hooke material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupLogMixNeoHooke",
                                            "mixed logarithmic neo-Hooke material",
                                            INPAR::MAT::mes_couplogmixneohooke));

    AddNamedString(m,"MODE","parameter set: YN (Young's modulus and Poisson's ration) or Lame (mue and lambda)", "YN");
    AddNamedReal(m,"C1","E or mue");
    AddNamedReal(m,"C2","nue or lambda");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled exponential material for compressible material (according to Weikenmeier_2014)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupExpPol",
                                            "compressible, isochoric exponential material law for soft tissue",
                                            INPAR::MAT::mes_coupexppol));
    AddNamedReal(m,"A","material constant");
    AddNamedReal(m,"B","material constant linear I_1");
    AddNamedReal(m,"C","material constant linear J");

    AppendMaterialDefinition(matlist,m);
  }



  /*--------------------------------------------------------------------*/

    //compressible neo-Hooke material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupNeoHooke",
                                            "compressible neo-Hooke material acc. to Holzapfel",
                                            INPAR::MAT::mes_coupneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus",0.0,true);
    AddNamedReal(m,"NUE","Poisson's ratio",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }
    //Mooney Rivlin  material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupMooneyRivlin",
                                            "Mooney - Rivlin material acc. to Holzapfel",
                                            INPAR::MAT::mes_coupmooneyrivlin));

    AddNamedReal(m,"C1","material constant",0.0,true);
    AddNamedReal(m,"C2","material constant",0.0,true);
    AddNamedReal(m,"C3","material constant",0.0,true);

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
    // optional parameters for uncertainty quantification
    AddNamedReal(m,"MULUM","mu for luminal pdf, irrelevant for deterministic analysis",0.0,true);
    AddNamedReal(m,"MUMED","mu for medial pdf, irrelevant for deterministic analysis",0.0,true);
    AddNamedReal(m,"MUABLUM","mu for abluminal pdf, irrelevant for deterministic analysis",0.0,true);
    AddNamedReal(m,"SIGMALUM","std for luminal pdf, irrelevant for deterministic analysis",0.0,true);
    AddNamedReal(m,"SIGMAMED","std for medial pdf, irrelevant for deterministic analysis",0.0,true);
    AddNamedReal(m,"SIGMAABLUM","std for abluminal pdf, irrelevant for deterministic analysis",0.0,true);

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
  // test material to test elasthyper-toolbox
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_IsoTestMaterial",
                                            "test material to test elasthyper-toolbox",
                                            INPAR::MAT::mes_isotestmaterial));

    AddNamedReal(m,"C1","Modulus for first invariant");
    AddNamedReal(m,"C2","Modulus for second invariant");
    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // general fiber material for remodeling
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_RemodelFiber",
                                            "General fiber material for remodeling",
                                            INPAR::MAT::mes_remodelfiber));

    AddNamedInt(m,"NUMMAT","number of materials/potentials in list");
    AddNamedIntVector(m,"MATIDS","the list material/potential IDs","NUMMAT");
    AddNamedReal(m, "TDECAY", "decay time of Poisson (degradation) process");
    AddNamedReal(m, "GROWTHFAC", "time constant for collagen growth",0.0,true);
    AddNamedRealVector(m,"COLMASSFRAC","initial mass fraction of first collagen fiber family in constraint mixture","NUMMAT", 0.0, true);
    AddNamedReal(m, "DEPOSITIONSTRETCH", "deposition stretch");

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
  // volumetric power law contribution
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_VolPow",
                                            "Power law formulation for the volumetric part",
                                            INPAR::MAT::mes_volpow));

    AddNamedReal(m,"A","prefactor of power law");
    AddNamedReal(m,"EXPON","exponent of power law");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with one exponential fiber family
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpoActive",
                                            "anisotropic active fiber",
                                            INPAR::MAT::mes_coupanisoexpoactive));

    AddNamedReal(m,"K1","linear constant");
    AddNamedReal(m,"K2","exponential constant");
    AddNamedReal(m,"GAMMA","angle");
    AddNamedReal(m,"K1COMP","linear constant");
    AddNamedReal(m,"K2COMP","exponential constant");
    AddNamedInt(m,"INIT","initialization modus for fiber alignment", 1, true);
    AddNamedBool(m,"ADAPT_ANGLE","adapt angle during remodeling", false, true);
    AddNamedReal(m,"S","maximum contractile stress");
    AddNamedReal(m,"LAMBDAMAX","stretch at maximum active force generation");
    AddNamedReal(m,"LAMBDA0","stretch at zero active force generation");
    AddNamedReal(m,"DENS","total reference mass density of constrained mixture");

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
  // coupled anisotropic material with one pow-like fiber family
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoPow",
                                            "anisotropic part with one pow-like fiber",
                                            INPAR::MAT::mes_coupanisopow));

    AddNamedReal(m,"K","linear constant");
    AddNamedReal(m,"D1","exponential constant for fiber invariant");
    AddNamedReal(m,"D2","exponential constant for system");
    AddNamedReal(m,"ACTIVETHRES","Deformation threshold for activating fibers. Default:"
        " 1.0 (off at compression); If 0.0 (always active)",1.0,true);
    AddNamedInt(m,"FIBER","Number of the fiber family contained in the element",1, true);
    AddNamedReal(m,"GAMMA","angle",0.0,true);
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
    AddNamedBool(m,"FIB_COMP","fibers support compression: yes (true) or no (false)", true, true);
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
    AddNamedReal(m,"ACTIVATION_THRES","Threshold for activation (contraction starts when activation function is larger than this value, relaxes otherwise)");
    AddNamedReal(m,"GAMMA","azimuth angle",0.0,true);
    AddNamedReal(m,"THETA","polar angle",0.0,true);
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
    AddNamedReal(m,"GAMMA","azimuth angle",0.0,true);
    AddNamedReal(m,"THETA","polar angle",0.0,true);
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

  /*--------------------------------------------------------------------*/
  // isotropic viscous contribution of myocardial matrix (chapelle12)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("VISCO_CoupMyocard",
                                            "Isotropic viscous contribution of myocardial matrix",
                                            INPAR::MAT::mes_coupmyocard));

    AddNamedReal(m,"N","material parameter");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric rate dependent viscos material, modified from Pioletti,1997
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("VISCO_IsoRateDep",
                                            "Isochoric rate dependent viscous material",
                                            INPAR::MAT::mes_isoratedep));

    AddNamedReal(m,"N","material parameter");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // vicous contribution to visohyperelastic material according to SLS-Model
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("VISCO_GenMax",
                                            "Viscous contribution according to SLS-Model",
                                            INPAR::MAT::mes_genmax));

    AddNamedReal(m,"TAU","relaxation parameter");
    AddNamedReal(m,"BETA","emphasis of viscous to elastic part");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/


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

    AddNamedReal(m,"CAPA","volumetric heat capacity");
    AddNamedReal(m,"CONDUCT","thermal conductivity");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material for heat transport due to Fourier-type thermal conduction and the Soret effect (fang 06/15)
  {
    Teuchos::RCP<MaterialDefinition> matsoret = Teuchos::rcp(new MaterialDefinition("MAT_soret","material for heat transport due to Fourier-type thermal conduction and the Soret effect",INPAR::MAT::m_soret));

    // mandatory parameters
    AddNamedReal(matsoret,"CAPA","volumetric heat capacity");
    AddNamedReal(matsoret,"CONDUCT","thermal conductivity");
    AddNamedReal(matsoret,"SORET","Soret coefficient");

    // add Soret material to global list of valid materials
    AppendMaterialDefinition(matlist,matsoret);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // integration point based growth
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthVolumetric",
                                            "volumetric growth",
                                            INPAR::MAT::m_growth_volumetric));

    AddNamedInt(m,"GROWTHLAW","number of growth law in input file");
    AddNamedInt(m,"IDMATELASTIC","number of elastic material in input file: MAT IDMATELASTIC ...");
    AddNamedReal(m,"STARTTIME","start growth after this time");
    AddNamedReal(m,"ENDTIME","end growth after this time");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
   // growth and remodeling (homogenized constrained mixture model)
   {
     Teuchos::RCP<MaterialDefinition> m
       = Teuchos::rcp(new MaterialDefinition("MAT_GrowthRemodel_ElastHyper",
                                             "growth and remodeling",
                                             INPAR::MAT::m_growthremodel_elasthyper));

     AddNamedInt(m,"NUMMATRF","number of remodelfiber materials in list",0,true);
     AddNamedInt(m,"NUMMATEL","number of elastin matrix materials/potentials in list",0,true);
     AddNamedIntVector(m,"MATIDSRF","the list remodelfiber material IDs","NUMMATRF",-1,true);
     AddNamedIntVector(m,"MATIDSEL","the list elastin matrix material/potential IDs","NUMMATEL",-1,true);
     AddNamedInt(m,"MATIDPENALTY","penalty material ID",-1,true);
     AddNamedRealVector(m,"ELMASSFRAC","initial mass fraction of elastin matrix in constraint mixture","NUMMATEL",0.0, true);
     AddNamedReal(m,"DENS","material mass density");
     AddNamedReal(m,"PRESTRESSTIME","time of prestressing",false);
     AddNamedReal(m,"PRESTRETCHELASTINCIR","circumferential prestretch of elastin matrix",false);
     AddNamedReal(m,"PRESTRETCHELASTINAX","axial prestretch of elastin matrix",false);
     AddNamedInt(m,"DAMAGE","1: elastin damage after prestressing,0: no elastin damage",false);
     AddNamedInt(m,"GROWTHTYPE","flag to decide what type of collagen growth is used: 1: anisotropic growth; 0: isotropic growth",false);
     AddNamedInt(m,"LOCTIMEINT","flag to decide what type of local time integration scheme is used: 1: Backward Euler Method; 0: Forward Euler Method",false);
     AppendMaterialDefinition(matlist,m);
   }

  /*----------------------------------------------------------------------*/
   // integration point based and scalar dependent growth
   {
     Teuchos::RCP<MaterialDefinition> mm
       = Teuchos::rcp(new MaterialDefinition("MAT_GrowthVolumetricScd",
                                             "integration point based and scalar dependent growth",
                                             INPAR::MAT::m_growth_volumetric_scd));

     AddNamedInt(mm,"GROWTHLAW","growth function: linear(Default) or exponential");
     AddNamedInt(mm,"IDMATELASTIC","number of elastic material in input file: MAT IDMATELASTIC ...");
     AddNamedReal(mm,"STARTTIME","start growth after this time");
     AddNamedReal(mm,"ENDTIME","end growth after this time");
     AddNamedReal(mm,"REARATE","substrate uptake rate coefficient",-1.0,true);
     AddNamedReal(mm,"SATCOEFF","saturation coefficient for concentration dependent growth law",-1.0,true);

     AppendMaterialDefinition(matlist,mm);
   }

   /*----------------------------------------------------------------------*/
    // integration point based and scalar dependent interpolation between to materials
    {
      Teuchos::RCP<MaterialDefinition> mm
        = Teuchos::rcp(new MaterialDefinition("MAT_ScDepInterp",
                                              "integration point based and scalar dependent interpolation between to materials",
                                              INPAR::MAT::m_sc_dep_interp));

      AddNamedInt(mm,"IDZEROCONCMAT","material for zero concentrations");
      AddNamedInt(mm,"IDINFTYCONCMAT","material for infinite concentrations");
//      AddNamedReal(mm,"ALPHA","size of ",-1.0,true);

      AppendMaterialDefinition(matlist,mm);
    }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // linear growth law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthLinear",
                                            "linear growth law",
                                            INPAR::MAT::m_growth_linear));

    AddNamedReal(m,"THETAPLUS","maximal growth stretch");
    AddNamedReal(m,"KPLUS","growth law parameter kthetaplus");
    AddNamedReal(m,"MPLUS","growth law parameter mthetaplus");
    AddNamedReal(m,"THETAMINUS","minimal growth stretch");
    AddNamedReal(m,"KMINUS","growth law parameter kthetaminus");
    AddNamedReal(m,"MMINUS","growth law parameter mthetaminus");
    AddNamedReal(m,"HOMMANDEL","homeostatic value for mandelstress");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // exponential growth law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthExponential",
                                            "exponential growth law",
                                            INPAR::MAT::m_growth_exponential));

    AddNamedReal(m,"MANDEL","reference value for mandelstress");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // biofilm growth law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthBiofilm",
                                            "biofilm growth law",
                                            INPAR::MAT::m_growth_biofilm));

    AddNamedReal(m,"MANDEL","reference value for mandelstress");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // constant rate growth law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthConst",
                                            "constant growth law",
                                            INPAR::MAT::m_growth_const));

    AddNamedReal(m,"THETARATE","reference value for mandelstress");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // simple atherosclerosis growth law, scalar depended volumetric growth
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthAC",
                                            "scalar depended volumetric growth",
                                            INPAR::MAT::m_growth_ac));

    AddNamedInt(m,"SCALAR1","number of first growth inducing scalar");
    AddNamedReal(m,"ALPHA","volume per first scalar's mass density");
    AddNamedInt(m,"SCALAR2","number of second growth inducing scalar",1,true);
    AddNamedReal(m,"BETA","volume per second scalar's mass density",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // atherosclerosis growth law, scalar depended growth in radial direction
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthACRadial",
                                            "scalar depended growth in radial direction",
                                            INPAR::MAT::m_growth_ac_radial));

    AddNamedInt(m,"SCALAR1","number of first growth inducing scalar");
    AddNamedReal(m,"ALPHA","volume per first scalar's mass density");
    AddNamedInt(m,"SCALAR2","number of second growth inducing scalar",1,true);
    AddNamedReal(m,"BETA","volume per second scalar's mass density",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  // atherosclerosis growth law, scalar depended growth in radial direction
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_GrowthACRadialRefConc",
                                            "scalar depended growth in radial direction",
                                            INPAR::MAT::m_growth_ac_radial_refconc));

    AddNamedInt(m,"SCALAR1","number of first growth inducing scalar");
    AddNamedReal(m,"ALPHA","volume per first scalar's mass density");
    AddNamedInt(m,"SCALAR2","number of second growth inducing scalar",1,true);
    AddNamedReal(m,"BETA","volume per second scalar's mass density",0.0,true);

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
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"PHIE","mass fraction of elastin");
    AddNamedReal(m,"PREELA","prestretch of elastin");
    AddNamedReal(m,"K1","Parameter for linear collagen fiber stiffness");
    AddNamedReal(m,"K2","Parameter for exponential collagen fiber stiffness");
    AddNamedInt(m,"NUMHOM","Number of homeostatic parameters",1);
    AddNamedRealVector(m,"PRECOLL","prestretch of collagen fibers","NUMHOM");
    AddNamedReal(m,"DAMAGE","damage stretch of collagen fibers");
    AddNamedReal(m,"K1M","Parameter for linear smooth muscle fiber stiffness");
    AddNamedReal(m,"K2M","Parameter for exponential smooth muscle fiber stiffness");
    AddNamedReal(m,"PHIM","mass fraction of smooth muscle");
    AddNamedReal(m,"PREMUS","prestretch of smooth muscle fibers");
    AddNamedReal(m,"SMAX","maximal active stress");
    AddNamedReal(m,"KAPPA","dilatation modulus");
    AddNamedReal(m,"LIFETIME","lifetime of collagen fibers");
    AddNamedReal(m,"GROWTHFAC","growth factor for stress");
    AddNamedRealVector(m,"HOMSTR","homeostatic target value of scalar stress measure","NUMHOM");
    AddNamedReal(m,"SHEARGROWTHFAC","growth factor for shear");
    AddNamedReal(m,"HOMRAD","homeostatic target value of inner radius");
    AddNamedReal(m,"STARTTIME","at this time turnover of collagen starts");
    AddNamedString(m,"INTEGRATION","time integration scheme (Explicit, Implicit)","Explicit");
    AddNamedReal(m,"TOL","tolerance for local Newton iteration, only for implicit integration");
    AddNamedString(m,"GROWTHFORCE","driving force of growth (Single, All, ElaCol)","Single");
    AddNamedString(m,"ELASTINDEGRAD","how elastin is degraded (None, Rectangle, Time)","None");
    AddNamedString(m,"MASSPROD","how mass depends on driving force (Lin, CosCos)","Lin");
    AddNamedString(m,"INITSTRETCH","how to set stretches in the beginning (None, Homeo, UpdatePrestretch)","None");
    AddNamedInt(m,"CURVE","number of timecurve for increase of prestretch in time",0);
    AddNamedString(m,"DEGOPTION","which degradation function (Lin, Cos, Exp, ExpVar)","Lin");
    AddNamedReal(m,"MAXMASSPRODFAC","maximal factor of mass production");
    AddNamedReal(m,"ELASTINFAC","factor for elastin content",0.0,true);
    AddNamedBool(m,"STOREHISTORY","store all history variables, not recommended for forward simulations",false,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // biofilm modeling (convection-diffusion-reaction equation)
  // it contains a non-linear reaction term following the Monod Kinetic of the form rearate*phi/(satcoeff+phi)
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
                                            "wrapper for structure poroelastic material",
                                            INPAR::MAT::m_structporo));

    AddNamedInt(m,"MATID","ID of structure material");
    AddNamedInt(m,"POROLAWID","ID of porosity law");
    AddNamedReal(m,"INITPOROSITY","initial porosity of porous medium");

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // linear law for porosity in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawLinear",
                                            "linear constitutive law for porosity",
                                            INPAR::MAT::m_poro_law_linear));

    AddNamedReal(m,"BULKMODULUS","bulk modulus of porous medium");

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // constant law for porosity in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawConstant",
                                            "constant constitutive law for porosity",
                                            INPAR::MAT::m_poro_law_constant));

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // neo-hookean law for porosity in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawNeoHooke",
                                            "NeoHookean-like constitutive law for porosity",
                                            INPAR::MAT::m_poro_law_logNeoHooke_Penalty));

    AddNamedReal(m,"BULKMODULUS","bulk modulus of porous medium");
    AddNamedReal(m,"PENALTYPARAMETER","penalty paramter of porous medium");

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // incompressible skeleton law for porosity in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawIncompSkel",
                                            "porosity law for incompressible skeleton phase",
                                            INPAR::MAT::m_poro_law_incompr_skeleton));

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // incompressible skeleton law for porosity in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawLinBiot",
                                            "linear biot model for porosity law",
                                            INPAR::MAT::m_poro_law_linear_biot));

    AddNamedReal(m,"INVBIOTMODULUS","inverse Biot modulus of porous medium");
    AddNamedReal(m,"BIOTCEOFF","Biot coefficient of porous medium");

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
    AddNamedInt(m,"POROLAWID","ID of porosity law");
    AddNamedReal(m,"INITPOROSITY","initial porosity of porous medium");
    AddNamedInt(m,"DOFIDREACSCALAR","Id of DOF within scalar transport problem, which controls the reaction");

    AppendMaterialDefinition(matlist,m);
  }
  /*----------------------------------------------------------------------*/
  // hyperelastic material for poroelasticity with reaction
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_StructPoroReactionECM",
                                            "wrapper for structure porelastic material with reaction",
                                            INPAR::MAT::m_structpororeactionECM));

    AddNamedInt(m,"MATID","ID of structure material");
    AddNamedInt(m,"POROLAWID","ID of porosity law");
    AddNamedReal(m,"INITPOROSITY","initial porosity of porous medium");
    AddNamedReal(m,"DENSCOLLAGEN","density of collagen");
    AddNamedInt(m,"DOFIDREACSCALAR","Id of DOF within scalar transport problem, which controls the reaction");
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
    //optional parameter
    AddNamedString(m,"PERMEABILITYFUNCTION","Permeability function: Const(Default) or Kozeny_Carman","Const",true);
  //  AddNamedReal(m,"BULKMODULUS","bulk modulus of medium");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // multiphase flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroMultiPhase",
                                            "multi phase flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo_multiphase));

    AddNamedBool(m,"LOCAL","individual materials allocated per element or only at global scope");
    AddNamedReal(m,"PERMEABILITY","permeability of medium");
    AddNamedInt(m,"NUMMAT","number of materials in list");
    AddNamedIntVector(m,"MATIDS","the list material IDs","NUMMAT");
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // one phase for multiphase flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhase",
                                            "one phase for multiphase flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo_singlephase));

    AddNamedReal(m,"DYNVISCOSITY","dynamic viscosity");
    AddNamedReal(m,"DENSITY","density");
    AddNamedReal(m,"BULKMODULUS","bulk modulus of phase");
    AddNamedReal(m,"PERMEABILITY","relative permeability of phase");
    AddNamedInt(m,"DOFTYPEID","ID of dof definition");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofDiffPressure",
                                            "one degrree of freedom for multiphase flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo_phasedof_diffpressure));

    AddNamedInt(m,"PHASELAWID","ID of pressure-saturation law");
    AddNamedInt(m,"NUMDOF","number of DoFs",0);
    AddNamedIntVector(m,"PRESCOEFF","pressure IDs for differential pressure","NUMDOF",0);
    AddNamedSeparator(m,"END","indicating end of line");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofPressure",
                                            "one degrree of freedom for multiphase flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo_phasedof_pressure));

    AddNamedInt(m,"PHASELAWID","ID of pressure-saturation law");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofPressureSum",
                                            "one degrree of freedom for multiphase flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo_phasedof_pressuresum));

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofSaturation",
                                            "one degrree of freedom for multiphase flow in deformable porous media",
                                            INPAR::MAT::m_fluidporo_phasedof_saturation));

    AddNamedInt(m,"PHASELAWID","ID of pressure-saturation law");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // saturated law for pressure-saturation law in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PhaseLawLinear",
                                            "saturated fluid phase of porous medium",
                                            INPAR::MAT::m_fluidporo_phaselaw_linear));

    AddNamedReal(m,"RELTENSION","relative interface tensions");
    AddNamedReal(m,"SATURATION_0","saturation at zero differential pressure");
    AddNamedInt(m,"NUMDOF","number of DoFs",0,true);
    AddNamedIntVector(m,"PRESCOEFF","Coefficients for pressure dependence","NUMDOF",0,true);
    AddNamedSeparator(m,"END","indicating end of line",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // tangent law for pressure-saturation law in porous media problems
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PhaseLawTangent",
                                            "tangent fluid phase of porous medium",
                                            INPAR::MAT::m_fluidporo_phaselaw_tangent));

    AddNamedReal(m,"RELTENSION","relative interface tensions");
    AddNamedReal(m,"EXP","exponent in pressure-saturation law");
    AddNamedReal(m,"SATURATION_0","saturation at zero differential pressure");
    AddNamedInt(m,"NUMDOF","number of DoFs",0,true);
    AddNamedIntVector(m,"PRESCOEFF","Coefficients for pressure dependence","NUMDOF",0,true);
    AddNamedSeparator(m,"END","indicating end of line",true);

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
  // 0D Acinar material base
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


  /*--------------------------------------------------------------------*/
    // 0D NeoHookean Acinar material
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_NEOHOOKEAN",
                                              "0D acinar material neohookean",
                                              INPAR::MAT::m_0d_maxwell_acinus_neohookean));

      AddNamedReal(m,"Stiffness1","first stiffness");
      AddNamedReal(m,"Stiffness2","second stiffness");
      AddNamedReal(m,"Viscosity1","first viscosity");
      AddNamedReal(m,"Viscosity2","second viscosity");

      AppendMaterialDefinition(matlist,m);
    }

    /*--------------------------------------------------------------------*/
    // 0D Exponential Acinar material
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_EXPONENTIAL",
                                              "0D acinar material exponential",
                                              INPAR::MAT::m_0d_maxwell_acinus_exponential));

      AddNamedReal(m,"Stiffness1","first stiffness");
      AddNamedReal(m,"Stiffness2","second stiffness");
      AddNamedReal(m,"Viscosity1","first viscosity");
      AddNamedReal(m,"Viscosity2","second viscosity");

      AppendMaterialDefinition(matlist,m);
    }

    /*--------------------------------------------------------------------*/
    // 0D Exponential Acinar material
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_DOUBLEEXPONENTIAL",
                                              "0D acinar material doubleexponential",
                                              INPAR::MAT::m_0d_maxwell_acinus_doubleexponential));

      AddNamedReal(m,"Stiffness1","first stiffness");
      AddNamedReal(m,"Stiffness2","second stiffness");
      AddNamedReal(m,"Viscosity1","first viscosity");
      AddNamedReal(m,"Viscosity2","second viscosity");

      AppendMaterialDefinition(matlist,m);
    }

  /*--------------------------------------------------------------------*/
  // 0D Ogden Acinar material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_OGDEN",
                                            "0D acinar material ogden",
                                            INPAR::MAT::m_0d_maxwell_acinus_ogden));

    AddNamedReal(m,"Stiffness1","first stiffness");
    AddNamedReal(m,"Stiffness2","second stiffness");
    AddNamedReal(m,"Viscosity1","first viscosity");
    AddNamedReal(m,"Viscosity2","second viscosity");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // O2 hemoglobin saturation material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_0D_O2_HEMOGLOBIN_SATURATION",
                                            "0D O2 hemoglobin saturation material",
                                            INPAR::MAT::m_0d_o2_hemoglobin_saturation));

    AddNamedReal(m,"PerVolumeBlood","how much of blood satisfies this rule (usually 100ml)");
    AddNamedReal(m,"O2SaturationPerVolBlood","O2 saturation per volume blood (In healthy blood 21.36ml/100ml of blood)");
    AddNamedReal(m,"PressureHalf","PO2 of 50\% saturated O2 (In healthy blood 26mmHg)");
    AddNamedReal(m,"Power","Power of the Sigmoidal saturation curve (2.5)");
    AddNamedReal(m,"NumberOfO2PerVO2","Number of O2 moles per unit volume of O2");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // O2 air saturation material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_0D_O2_AIR_SATURATION",
                                            "0D O2 air saturation material",
                                            INPAR::MAT::m_0d_o2_air_saturation));

    AddNamedReal(m,"AtmosphericPressure","The atmospheric pressure");
    AddNamedReal(m,"NumberOfO2PerVO2","Number of O2 moles per unit volume of O2");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // particle material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Particle",
                                            "particle material",
                                            INPAR::MAT::m_particlemat));

    AddNamedReal(m,"DENSITY","initial mass density");
    AddNamedReal(m,"INITRADIUS","initial radius of particle");
    AddNamedReal(m,"NUE","poisson ratio",0.0,true);
    AddNamedReal(m,"YOUNG","youngs modulus",0.0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // multipurpose meshfree/particle material for additive manufacturing
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ExtParticle",
                                            "particle material",
                                            INPAR::MAT::m_extparticlemat));

    AddNamedReal(m,"DENSITY","initial mass density");
    AddNamedReal(m,"INITRADIUS","initial radius of particle");
    AddNamedReal(m,"NUE","poisson ratio",0.0,true);
    AddNamedReal(m,"YOUNG","youngs modulus",0.0,true);
    AddNamedReal(m,"INITTEMPERATURE","initial temperature");
    AddNamedReal(m,"CPS","specific heat - constant pressure - solid state");
    AddNamedReal(m,"CPL","specific heat - constant pressure - liquid state");
    AddNamedReal(m,"SL_LATENT_HEAT","specific latent heat - solid <-> liquid");
    AddNamedReal(m,"SL_TRANSITION_TEMPERATURE","transition temperature - solid <-> liquid");
    AddNamedReal(m,"S_THERMAL_EXPANSION","volumetric thermal expansion coefficient - solid state - usually 3 times the linear coefficient");
    AddNamedReal(m,"L_THERMAL_EXPANSION","volumetric thermal expansion coefficient - liquid state");
    AddNamedReal(m,"SL_THERMAL_EXPANSION","volumetric thermal expansion coefficient - solid <-> liquid - Beware! It is linked to latent heat instead of temperature!");

    AppendMaterialDefinition(matlist,m);
  }

   /*----------------------------------------------------------------------*/
   // nutrient diffusion modeling (diffusion-reaction equation). Contains a growth-dependent
   //reaction-term of the form 3*theta^2*time_derivative_theta*structure_density that is coupled with
   // scalar dependent growth law via 'theta'
   // and a non-linear reaction term following the Monod Kinetic of the form rearate*phi/(satcoeff+phi)
   {
     Teuchos::RCP<MaterialDefinition> m
       = Teuchos::rcp(new MaterialDefinition("MAT_Scatra_GrowthScd",
                                             "nutrientdiff material",
                                             INPAR::MAT::m_scatra_growth_scd));

     AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
     AddNamedReal(m,"STRDENSITY","density of structure material");
     AddNamedReal(m,"REARATE","substrate uptake rate coefficient",0.0);
     AddNamedReal(m,"SATCOEFF","substrate saturation coefficient",0.0);
     AddNamedString(m,"SOURCEMASS","source mass term","Standard", true);

     // AddNamedString(m,"KINETICS","Substrate consumption kinetics (SimpleMonod)","SimpleMonod");
     AppendMaterialDefinition(matlist,m);
   }

  /*----------------------------------------------------------------------*/
  // acoustic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Acoustic",
                                            "acoustic material",
                                            INPAR::MAT::m_acousticmat));

    AddNamedReal(m,"DENSITY","mass density");
    AddNamedReal(m,"C","speed of sound");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
    // acoustic solid material
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_AcousticSol",
                                              "acoustic solid material",
                                              INPAR::MAT::m_acousticsolmat));

      AddNamedReal(m,"DENSITY","mass density");
      AddNamedReal(m,"C","speed of sound");
      AddNamedReal(m,"VISC","viscosity mu");


      AppendMaterialDefinition(matlist,m);
    }

  /*----------------------------------------------------------------------*/

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    // active fiber formation for the modeling of living cells
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_ACTIVEFIBER",
                                              "active fiber formation for the modeling of living cells",
                                              INPAR::MAT::m_activefiber));

      AddNamedReal(m,"DENS","Density");
      AddNamedReal(m,"DECAY","decay constant of activation signal");
      AddNamedInt(m,"IDMATPASSIVE","number of passive material in input file: MAT IDMATPASSIVE ...");
      AddNamedReal(m,"KFOR","formation rate parameter kforwards");
      AddNamedReal(m,"KBACK","dissociation parameter kbackwards");
      AddNamedReal(m,"KVAR","fiber rate sensitivity");
      AddNamedReal(m,"SIGMAX","maximum tension exerted by stress fibres");
      AddNamedReal(m,"EPSNULL","reference strain rate of cross-bridge dynamics law");


      AppendMaterialDefinition(matlist,m);
    }

    /*----------------------------------------------------------------------*/
    // biochemo-mechano coupled active stress fiber formation for cells
    {
      Teuchos::RCP<MaterialDefinition> m
        = Teuchos::rcp(new MaterialDefinition("MAT_BIOCHEMOMECHANO",
                                              "biochemo-mechano coupled active stress fiber formation for cells",
                                               INPAR::MAT::m_biochemomechano));

          AddNamedReal(m,"DENS","Density");
          AddNamedInt(m,"IDMATPASSIVE","number of passive material in input file: MAT IDMATPASSIVE ...");
          AddNamedReal(m,"KFOR","formation rate parameter kforwards");
          AddNamedReal(m,"KBACK","dissociation rate parameter kbackwards");
          AddNamedReal(m,"KROCKETA","rate constant for Nfil to ROCK");
          AddNamedReal(m,"KACTIN","rate constant for Nfil to actin");
          AddNamedReal(m,"RATEMAX","maximum rate");
          AddNamedReal(m,"NMAX","maximum filament concentration");
          AddNamedReal(m,"KSTRESS","proportionality constant between Acto-Mysion-Activation and stress");
          AddNamedReal(m,"SOURCE","Constant for the Source for Stress-dependent Surface Scatra Condition");
          AddNamedString(m,"METHOD","Method for evaluating the material specific integral","2DGauss");
          AppendMaterialDefinition(matlist,m);
    }

  // deliver
  return vm;
}
