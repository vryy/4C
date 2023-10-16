/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid materials for input

\level 1

*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_validmaterials.H"

#include "baci_inpar_material.H"
#include "baci_lib_materialdefinition.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintEmptyMaterialDefinitions(
    std::ostream& stream, std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& matlist)
{
  const std::string sectionname = "MATERIALS";
  const unsigned l = sectionname.length();
  stream << "--" << std::string(std::max<int>(65 - l, 0), '-');
  stream << sectionname << '\n';

  for (auto& i : matlist)
  {
    i->Print(stream, nullptr);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintMaterialDatHeader()
{
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>> matlist =
      DRT::INPUT::ValidMaterials();
  DRT::INPUT::PrintEmptyMaterialDefinitions(std::cout, *matlist);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>> DRT::INPUT::ValidMaterials()
{
  using Teuchos::tuple;

  // a list containing all valid materials
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>> vm =
      Teuchos::rcp(new std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>());

  // convenience
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& matlist = *vm;


  /*----------------------------------------------------------------------*/
  // Newtonian fluid
  {
    auto m =
        Teuchos::rcp(new MaterialDefinition("MAT_fluid", "Newtonian fluid", INPAR::MAT::m_fluid));

    ::INPUT::AddNamedReal(m, "DYNVISCOSITY", "dynamic viscosity");
    ::INPUT::AddNamedReal(m, "DENSITY", "spatial mass density");
    ::INPUT::AddNamedReal(m, "GAMMA", "surface tension coefficient", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Weakly compressible fluid according to Murnaghan-Tait
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_fluid_murnaghantait",
        "Weakly compressible fluid according to Murnaghan-Tait",
        INPAR::MAT::m_fluid_murnaghantait));

    ::INPUT::AddNamedReal(m, "DYNVISCOSITY", "dynamic viscosity");
    ::INPUT::AddNamedReal(m, "REFDENSITY", "reference spatial mass density");
    ::INPUT::AddNamedReal(m, "REFPRESSURE", "reference pressure");
    ::INPUT::AddNamedReal(m, "REFBULKMODULUS", "reference bulk modulus");
    ::INPUT::AddNamedReal(m, "MATPARAMETER", "material parameter according to Murnaghan-Tait");
    ::INPUT::AddNamedReal(m, "GAMMA", "surface tension coefficient", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Linear law (pressure-dependent) for the density and the viscosity
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_fluid_linear_density_viscosity",
        "Linear law (pressure-dependent) for the density and the viscosity",
        INPAR::MAT::m_fluid_linear_density_viscosity));

    ::INPUT::AddNamedReal(m, "REFDENSITY", "reference density");
    ::INPUT::AddNamedReal(m, "REFVISCOSITY", "reference viscosity");
    ::INPUT::AddNamedReal(m, "REFPRESSURE", "reference pressure");
    ::INPUT::AddNamedReal(m, "COEFFDENSITY", "density-pressure coefficient");
    ::INPUT::AddNamedReal(m, "COEFFVISCOSITY", "viscosity-pressure coefficient");
    ::INPUT::AddNamedReal(m, "GAMMA", "surface tension coefficient", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Weakly compressible fluid
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_fluid_weakly_compressible",
        "Weakly compressible fluid", INPAR::MAT::m_fluid_weakly_compressible));

    ::INPUT::AddNamedReal(m, "VISCOSITY", "viscosity");
    ::INPUT::AddNamedReal(m, "REFDENSITY", "reference density");
    ::INPUT::AddNamedReal(m, "REFPRESSURE", "reference pressure");
    ::INPUT::AddNamedReal(m, "COMPRCOEFF", "compressibility coefficient");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // fluid with non-linear viscosity according to Carreau-Yasuda
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_carreauyasuda",
        "fluid with non-linear viscosity according to Carreau-Yasuda",
        INPAR::MAT::m_carreauyasuda));

    ::INPUT::AddNamedReal(m, "NU_0", "zero-shear viscosity");
    ::INPUT::AddNamedReal(m, "NU_INF", "infinite-shear viscosity");
    ::INPUT::AddNamedReal(m, "LAMBDA", "characteristic time");
    ::INPUT::AddNamedReal(m, "APARAM", "constant parameter");
    ::INPUT::AddNamedReal(m, "BPARAM", "constant parameter");
    ::INPUT::AddNamedReal(m, "DENSITY", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // fluid with nonlinear viscosity according to a modified power law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_modpowerlaw",
        "fluid with nonlinear viscosity according to a modified power law",
        INPAR::MAT::m_modpowerlaw));

    ::INPUT::AddNamedReal(m, "MCONS", "consistency");
    ::INPUT::AddNamedReal(m, "DELTA", "safety factor");
    ::INPUT::AddNamedReal(m, "AEXP", "exponent");
    ::INPUT::AddNamedReal(m, "DENSITY", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // fluid with non-linear viscosity according to Herschel-Bulkley
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_herschelbulkley",
        "fluid with non-linear viscosity according to Herschel-Bulkley",
        INPAR::MAT::m_herschelbulkley));

    ::INPUT::AddNamedReal(m, "TAU_0", "yield stress");
    ::INPUT::AddNamedReal(m, "KFAC", "constant factor");
    ::INPUT::AddNamedReal(m, "NEXP", "exponent");
    ::INPUT::AddNamedReal(m, "MEXP", "exponent");
    ::INPUT::AddNamedReal(m, "LOLIMSHEARRATE", "lower limit of shear rate");
    ::INPUT::AddNamedReal(m, "UPLIMSHEARRATE", "upper limit of shear rate");
    ::INPUT::AddNamedReal(m, "DENSITY", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // "yoghurt-type" fluid with nonlinear viscosity according to a power law
  // and extended by an Arrhenius-type term to account for temperature dependence
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_yoghurt", "yoghurt-type fluid with nonlinear viscosity", INPAR::MAT::m_yoghurt));

    ::INPUT::AddNamedReal(m, "SHC", "specific heat capacity at constant pressure (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "DENSITY", "density");
    ::INPUT::AddNamedReal(m, "THERMCOND", "thermal conductivity (J/(m*K*s))");
    ::INPUT::AddNamedReal(m, "STRAINRATEEXP", "exponent of strain-rate term");
    ::INPUT::AddNamedReal(m, "PREEXCON", "pre-exponential constant (1/s)");
    ::INPUT::AddNamedReal(m, "ACTENERGY", "activation energy (J/kg)");
    ::INPUT::AddNamedReal(m, "GASCON", "specific gas constant R (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "DELTA", "safety factor");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // fluid flow in a permeable material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_permeable", "permeability for flow in porous media", INPAR::MAT::m_permeable_fluid));

    ::INPUT::AddNamedString(
        m, "TYPE", "Problem type: Darcy, Darcy-Stokes (default)", "Darcy-Stokes");
    ::INPUT::AddNamedReal(m, "DYNVISCOSITY", "dynamic viscosity");
    ::INPUT::AddNamedReal(m, "DENSITY", "density");
    ::INPUT::AddNamedReal(m, "PERMEABILITY", "permeability of medium");

    AppendMaterialDefinition(matlist, m);
  }


  /*----------------------------------------------------------------------*/
  // lubrication material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_lubrication", "lubrication material", INPAR::MAT::m_lubrication));

    ::INPUT::AddNamedInt(m, "LUBRICATIONLAWID", "lubrication law id");
    ::INPUT::AddNamedReal(m, "DENSITY", "lubricant density");

    AppendMaterialDefinition(matlist, m);
  }


  /*----------------------------------------------------------------------*/
  // constant lubrication material law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_lubrication_law_constant",
        "constant lubrication material law", INPAR::MAT::m_lubrication_law_constant));

    ::INPUT::AddNamedReal(m, "VISCOSITY", "lubricant viscosity");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Barus viscosity lubrication material law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_lubrication_law_barus",
        "barus lubrication material law", INPAR::MAT::m_lubrication_law_barus));

    ::INPUT::AddNamedReal(m, "ABSViscosity", "absolute lubricant viscosity");
    ::INPUT::AddNamedReal(m, "PreVisCoeff", "pressure viscosity coefficient");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Roeland viscosity lubrication material law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_lubrication_law_roeland",
        "roeland lubrication material law", INPAR::MAT::m_lubrication_law_roeland));

    ::INPUT::AddNamedReal(m, "ABSViscosity", "absolute lubricant viscosity");
    ::INPUT::AddNamedReal(m, "PreVisCoeff", "pressure viscosity coefficient");
    ::INPUT::AddNamedReal(m, "RefVisc", "reference viscosity");
    ::INPUT::AddNamedReal(m, "RefPress", "reference Pressure");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport material (with potential reaction coefficient)
  {
    auto m = Teuchos::rcp(
        new MaterialDefinition("MAT_scatra", "scalar transport material", INPAR::MAT::m_scatra));

    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity");
    ::INPUT::AddNamedReal(m, "REACOEFF", "reaction coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "SCNUM", "schmidt number", 0.0, true);
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }


  /*----------------------------------------------------------------------*/
  // scalar transport material (with potential reaction coefficient)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_scatra_reaction_poro",
        "scalar transport material", INPAR::MAT::m_scatra_reaction_poroECM));

    ::INPUT::AddNamedInt(m, "NUMSCAL", "number of scalars for these elements");
    ::INPUT::AddNamedIntVector(m, "STOICH", "reaction stoichometrie list", "NUMSCAL");
    ::INPUT::AddNamedReal(m, "REACCOEFF", "reaction coefficient");
    ::INPUT::AddNamedReal(m, "REACSCALE", "scaling for reaction coefficient");
    // reacscale could now be done by constant distribution function
    ::INPUT::AddNamedInt(m, "DISTRFUNCT", "spatial distribution of reaction coefficient", 0, true);
    ::INPUT::AddNamedString(m, "COUPLING",
        "type of coupling: "
        "simple_multiplicative, power_multiplicative, constant, michaelis_menten, by_function, "
        "no_coupling (default)",
        "no_coupling", false);
    ::INPUT::AddNamedRealVector(m, "ROLE", "role in michaelis-menten like reactions", "NUMSCAL");
    ::INPUT::AddNamedRealVector(m, "REACSTART", "starting point of reaction", "NUMSCAL", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }
  /*----------------------------------------------------------------------*/
  // scalar transport reaction material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_scatra_reaction", "advanced reaction material", INPAR::MAT::m_scatra_reaction));

    ::INPUT::AddNamedInt(m, "NUMSCAL", "number of scalars for these elements");
    ::INPUT::AddNamedIntVector(m, "STOICH", "reaction stoichometrie list", "NUMSCAL");
    ::INPUT::AddNamedReal(m, "REACCOEFF", "reaction coefficient");
    ::INPUT::AddNamedInt(m, "DISTRFUNCT", "spatial distribution of reaction coefficient", 0, true);
    ::INPUT::AddNamedString(m, "COUPLING",
        "type of coupling: "
        "simple_multiplicative, power_multiplicative, constant, michaelis_menten, by_function, "
        "no_coupling (default)",
        "no_coupling", false);
    ::INPUT::AddNamedRealVector(m, "ROLE", "role in michaelis-menten like reactions", "NUMSCAL");
    ::INPUT::AddNamedRealVector(m, "REACSTART", "starting point of reaction", "NUMSCAL", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport reaction material (species in fluid)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiporo_fluid",
        "advanced reaction material for multiphase porous flow (species in fluid)",
        INPAR::MAT::m_scatra_multiporo_fluid));

    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity");
    ::INPUT::AddNamedInt(m, "PHASEID", "ID of fluid phase the scalar is associated with");
    ::INPUT::AddNamedReal(m, "REACOEFF", "reaction coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "SCNUM", "schmidt number", 0.0, true);
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "DELTA", "delta", 0.0, true);
    ::INPUT::AddNamedReal(m, "MIN_SAT",
        "minimum saturation under which also corresponding mass fraction is equal to zero", 1.0e-9,
        true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport reaction material (species in volume fraction)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiporo_volfrac",
        "advanced reaction material for multiphase porous flow (species in volfrac)",
        INPAR::MAT::m_scatra_multiporo_volfrac));

    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity");
    ::INPUT::AddNamedInt(m, "PHASEID", "ID of fluid phase the scalar is associated with");
    ::INPUT::AddNamedReal(m, "REACOEFF", "reaction coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "SCNUM", "schmidt number", 0.0, true);
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "DELTA", "delta", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport reaction material (species in solid)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiporo_solid",
        "advanced reaction material for multiphase "
        "porous flow (species in solid)",
        INPAR::MAT::m_scatra_multiporo_solid));

    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity");
    // no phaseID because only one solid phase
    ::INPUT::AddNamedReal(m, "REACOEFF", "reaction coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "SCNUM", "schmidt number", 0.0, true);
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "DELTA", "delta", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport reaction material (temperature)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiporo_temperature",
        "advanced reaction material for multiphase porous flow (temperature)",
        INPAR::MAT::m_scatra_multiporo_temperature));

    ::INPUT::AddNamedInt(m, "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE", "number of fluid dofs");
    ::INPUT::AddNamedRealVector(
        m, "CP_FLUID", "heat capacity fluid phases", "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE");
    ::INPUT::AddNamedInt(m, "NUMVOLFRAC", "number of volfrac dofs");
    ::INPUT::AddNamedRealVector(m, "CP_VOLFRAC", "heat capacity volfrac", "NUMVOLFRAC");
    ::INPUT::AddNamedReal(m, "CP_SOLID", "heat capacity solid");
    ::INPUT::AddNamedRealVector(m, "KAPPA_FLUID", "thermal diffusivity fluid phases",
        "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE");
    ::INPUT::AddNamedRealVector(m, "KAPPA_VOLFRAC", "thermal diffusivity volfrac", "NUMVOLFRAC");
    ::INPUT::AddNamedReal(m, "KAPPA_SOLID", "heat capacity solid");
    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity", 1.0, true);
    ::INPUT::AddNamedReal(m, "REACOEFF", "reaction coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "SCNUM", "schmidt number", 0.0, true);
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport chemotaxis material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_scatra_chemotaxis", "chemotaxis material", INPAR::MAT::m_scatra_chemotaxis));

    ::INPUT::AddNamedInt(m, "NUMSCAL", "number of chemotactic pairs for these elements");
    ::INPUT::AddNamedIntVector(m, "PAIR", "chemotaxis pairing", "NUMSCAL");
    ::INPUT::AddNamedReal(m, "CHEMOCOEFF", "chemotaxis coefficient");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // anisotropic scalar transport material (with potential reaction coefficient)
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_scatra_aniso", "anisotropic scalar transport material", INPAR::MAT::m_scatra_aniso));

    ::INPUT::AddNamedReal(m, "DIFF1", "kinematic diffusivity component 1");
    ::INPUT::AddNamedReal(m, "DIFF2", "kinematic diffusivity component 2");
    ::INPUT::AddNamedReal(m, "DIFF3", "kinematic diffusivity component 3");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // scalar transport material for multi-scale approach
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_scatra_multiscale",
        "scalar transport material for multi-scale approach", INPAR::MAT::m_scatra_multiscale));

    ::INPUT::AddNamedString(m, "MICROFILE", "input file for micro scale", "filename.dat");
    ::INPUT::AddNamedInt(m, "MICRODIS_NUM", "number of micro-scale discretization");
    ::INPUT::AddNamedReal(m, "POROSITY", "porosity");
    ::INPUT::AddNamedReal(m, "TORTUOSITY", "tortuosity");
    ::INPUT::AddNamedReal(m, "A_s", "specific micro-scale surface area");
    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity");
    ::INPUT::AddNamedReal(m, "REACOEFF", "reaction coefficient", 0.0, true);
    ::INPUT::AddNamedReal(m, "SCNUM", "Schmidt number", 0.0, true);
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Weickenmeier muscle material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Muscle_Weickenmeier",
        "Weickenmeier muscle material", INPAR::MAT::m_muscle_weickenmeier));

    ::INPUT::AddNamedReal(m, "ALPHA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "BETA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "GAMMA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "KAPPA", "material parameter for coupled volumetric contribution");
    ::INPUT::AddNamedReal(m, "OMEGA0", "weighting factor for isotropic tissue constituents");
    ::INPUT::AddNamedReal(
        m, "ACTMUNUM", "number of active motor units per undeformed muscle cross-sectional area");
    ::INPUT::AddNamedInt(m, "MUTYPESNUM", "number of motor unit types");
    ::INPUT::AddNamedRealVector(m, "INTERSTIM", "interstimulus interval", "MUTYPESNUM");
    ::INPUT::AddNamedRealVector(m, "FRACACTMU", "fraction of motor unit type", "MUTYPESNUM");
    ::INPUT::AddNamedRealVector(m, "FTWITCH", "twitch force of motor unit type", "MUTYPESNUM");
    ::INPUT::AddNamedRealVector(
        m, "TTWITCH", "twitch contraction time of motor unit type", "MUTYPESNUM");
    ::INPUT::AddNamedReal(m, "LAMBDAMIN", "minimal active fiber stretch");
    ::INPUT::AddNamedReal(
        m, "LAMBDAOPT", "optimal active fiber stretch related to active nominal stress maximum");
    ::INPUT::AddNamedReal(m, "DOTLAMBDAMIN", "minimal stretch rate");
    ::INPUT::AddNamedReal(m, "KE",
        "parameter controlling the curvature of the velocity dependent activation function in the "
        "eccentric case");
    ::INPUT::AddNamedReal(m, "KC",
        "parameter controlling the curvature of the velocity dependent activation function in the "
        "concentric case");
    ::INPUT::AddNamedReal(m, "DE",
        "parameter controlling the amplitude of the velocity dependent activation function in the "
        "eccentric case");
    ::INPUT::AddNamedReal(m, "DC",
        "parameter controlling the amplitude of the velocity dependent activation function in the "
        "concentric case");
    ::INPUT::AddNamedInt(m, "ACTTIMESNUM", "number of time boundaries to prescribe activation");
    ::INPUT::AddNamedRealVector(m, "ACTTIMES", "time boundaries between intervals", "ACTTIMESNUM");
    ::INPUT::AddNamedInt(m, "ACTINTERVALSNUM", "number of time intervals to prescribe activation");
    ::INPUT::AddNamedRealVector(m, "ACTVALUES",
        "scaling factor in intervals (1=full activation, 0=no activation)", "ACTINTERVALSNUM");
    ::INPUT::AddNamedReal(m, "DENS", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Combo muscle material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_Muscle_Combo", "Combo muscle material", INPAR::MAT::m_muscle_combo));

    ::INPUT::AddNamedReal(m, "ALPHA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "BETA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "GAMMA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "KAPPA", "material parameter for coupled volumetric contribution");
    ::INPUT::AddNamedReal(m, "OMEGA0", "weighting factor for isotropic tissue constituents");
    ::INPUT::AddNamedReal(m, "POPT", "tetanised optimal (maximal) active stress");
    ::INPUT::AddNamedReal(m, "LAMBDAMIN", "minimal active fiber stretch");
    ::INPUT::AddNamedReal(
        m, "LAMBDAOPT", "optimal active fiber stretch related to active nominal stress maximum");
    ::INPUT::AddNamedReal(m, "C", "constant scaling tanh-type activation function");
    ::INPUT::AddNamedReal(m, "ACTSTARTTIME", "starting time of muscle activation");
    ::INPUT::AddNamedReal(m, "DENS", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Active strain Giantesio muscle material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Muscle_Giantesio",
        "Giantesio active strain muscle material", INPAR::MAT::m_muscle_giantesio));

    ::INPUT::AddNamedReal(m, "ALPHA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "BETA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "GAMMA", "experimentally fitted material parameter");
    ::INPUT::AddNamedReal(m, "KAPPA", "material parameter for coupled volumetric contribution");
    ::INPUT::AddNamedReal(m, "OMEGA0", "weighting factor for isotropic tissue constituents");
    ::INPUT::AddNamedReal(
        m, "ACTMUNUM", "number of active motor units per undeformed muscle cross-sectional area");
    ::INPUT::AddNamedInt(m, "MUTYPESNUM", "number of motor unit types");
    ::INPUT::AddNamedRealVector(m, "INTERSTIM", "interstimulus interval", "MUTYPESNUM");
    ::INPUT::AddNamedRealVector(m, "FRACACTMU", "fraction of motor unit type", "MUTYPESNUM");
    ::INPUT::AddNamedRealVector(m, "FTWITCH", "twitch force of motor unit type", "MUTYPESNUM");
    ::INPUT::AddNamedRealVector(
        m, "TTWITCH", "twitch contraction time of motor unit type", "MUTYPESNUM");
    ::INPUT::AddNamedReal(m, "LAMBDAMIN", "minimal active fiber stretch");
    ::INPUT::AddNamedReal(
        m, "LAMBDAOPT", "optimal active fiber stretch related to active nominal stress maximum");
    ::INPUT::AddNamedReal(m, "DOTLAMBDAMIN", "minimal stretch rate");
    ::INPUT::AddNamedReal(m, "KE",
        "parameter controlling the curvature of the velocity dependent activation function in the "
        "eccentric case");
    ::INPUT::AddNamedReal(m, "KC",
        "parameter controlling the curvature of the velocity dependent activation function in the "
        "concentric case");
    ::INPUT::AddNamedReal(m, "DE",
        "parameter controlling the amplitude of the velocity dependent activation function in the "
        "eccentric case");
    ::INPUT::AddNamedReal(m, "DC",
        "parameter controlling the amplitude of the velocity dependent activation function in the "
        "concentric case");
    ::INPUT::AddNamedInt(m, "ACTTIMESNUM", "number of time boundaries to prescribe activation");
    ::INPUT::AddNamedRealVector(m, "ACTTIMES", "time boundaries between intervals", "ACTTIMESNUM");
    ::INPUT::AddNamedInt(m, "ACTINTERVALSNUM", "number of time intervals to prescribe activation");
    ::INPUT::AddNamedRealVector(m, "ACTVALUES",
        "scaling factor in intervals (1=full activation, 0=no activation)", "ACTINTERVALSNUM");
    ::INPUT::AddNamedReal(m, "DENS", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Myocard muscle material (with complicated reaction coefficient)
  {
    auto m = Teuchos::rcp(
        new MaterialDefinition("MAT_myocard", "Myocard muscle material", INPAR::MAT::m_myocard));

    ::INPUT::AddNamedReal(m, "DIFF1", "conductivity in fiber direction");
    ::INPUT::AddNamedReal(m, "DIFF2", "conductivity perpendicular to fiber direction");
    ::INPUT::AddNamedReal(m, "DIFF3", "conductivity perpendicular to fiber direction");
    ::INPUT::AddNamedReal(
        m, "PERTUBATION_DERIV", "pertubation for calculation of reaction coefficient derivative");
    ::INPUT::AddNamedString(m, "MODEL", "Model type: MV (default), FHN, TNNP, SAN or INADA", "MV");
    ::INPUT::AddNamedString(m, "TISSUE", "Tissue type: M (default), ENDO, EPI, AN, N or NH", "M");
    ::INPUT::AddNamedReal(m, "TIME_SCALE", "Scale factor for time units of Model");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material according to mixture-fraction approach
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_mixfrac", "material according to mixture-fraction approach", INPAR::MAT::m_mixfrac));

    ::INPUT::AddNamedReal(m, "KINVISC", "kinematic viscosity");
    ::INPUT::AddNamedReal(m, "KINDIFF", "kinematic diffusivity");
    ::INPUT::AddNamedReal(m, "EOSFACA", "equation-of-state factor a");
    ::INPUT::AddNamedReal(m, "EOSFACB", "equation-of-state factor b");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_sutherland", "material according to Sutherland law", INPAR::MAT::m_sutherland));

    ::INPUT::AddNamedReal(m, "REFVISC", "reference dynamic viscosity (kg/(m*s))");
    ::INPUT::AddNamedReal(m, "REFTEMP", "reference temperature (K)");
    ::INPUT::AddNamedReal(m, "SUTHTEMP", "Sutherland temperature (K)");
    ::INPUT::AddNamedReal(m, "SHC", "specific heat capacity at constant pressure (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "PRANUM", "Prandtl number");
    ::INPUT::AddNamedReal(m, "THERMPRESS", "(initial) thermodynamic pressure (J/m^3)");
    ::INPUT::AddNamedReal(m, "GASCON", "specific gas constant R (J/(kg*K))");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material for temperature-dependent water according to VDI Waermeatlas
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_tempdepwater",
        "material for temperature-dependent water", INPAR::MAT::m_tempdepwater));

    ::INPUT::AddNamedReal(m, "CRITDENS", "critical density (kg/m^3)");
    ::INPUT::AddNamedReal(m, "CRITTEMP", "critical temperature (K)");
    ::INPUT::AddNamedReal(m, "SHC", "specific heat capacity at constant pressure (J/(kg*K))");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with Arrhenius-type chemical
  // kinetics (species)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_arrhenius_spec",
        "Arrhenius-type chemical kinetics (species)", INPAR::MAT::m_arrhenius_spec));

    ::INPUT::AddNamedReal(m, "REFVISC", "reference dynamic viscosity (kg/(m*s))");
    ::INPUT::AddNamedReal(m, "REFTEMP", "reference temperature (K)");
    ::INPUT::AddNamedReal(m, "SUTHTEMP", "Sutherland temperature (K)");
    ::INPUT::AddNamedReal(m, "SCHNUM", "Schmidt number");
    ::INPUT::AddNamedReal(m, "PREEXCON", "pre-exponential constant (1/s)");
    ::INPUT::AddNamedReal(m, "TEMPEXP", "exponent of temperature dependence");
    ::INPUT::AddNamedReal(m, "ACTEMP", "activation temperature (K)");
    ::INPUT::AddNamedReal(m, "GASCON", "specific gas constant R (J/(kg*K))");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with Arrhenius-type chemical
  // kinetics (temperature)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_arrhenius_temp",
        "Arrhenius-type chemical kinetics (temperature)", INPAR::MAT::m_arrhenius_temp));

    ::INPUT::AddNamedReal(m, "REFVISC", "reference dynamic viscosity (kg/(m*s))");
    ::INPUT::AddNamedReal(m, "REFTEMP", "reference temperature (K)");
    ::INPUT::AddNamedReal(m, "SUTHTEMP", "Sutherland temperature (K)");
    ::INPUT::AddNamedReal(m, "SHC", "specific heat capacity at constant pressure (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "PRANUM", "Prandtl number");
    ::INPUT::AddNamedReal(m, "REAHEAT", "heat of reaction per unit mass (J/kg)");
    ::INPUT::AddNamedReal(m, "PREEXCON", "pre-exponential constant (1/s)");
    ::INPUT::AddNamedReal(m, "TEMPEXP", "exponent of temperature dependence");
    ::INPUT::AddNamedReal(m, "ACTEMP", "activation temperature (K)");
    ::INPUT::AddNamedReal(m, "GASCON", "specific gas constant R (J/(kg*K))");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with Arrhenius-type chemical
  // kinetics (progress variable)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_arrhenius_pv",
        "material with Arrhenius-type chemical kinetics (progress variable)",
        INPAR::MAT::m_arrhenius_pv));

    ::INPUT::AddNamedReal(m, "REFVISC", "reference dynamic viscosity (kg/(m*s))");
    ::INPUT::AddNamedReal(m, "REFTEMP", "reference temperature (K)");
    ::INPUT::AddNamedReal(m, "SUTHTEMP", "Sutherland temperature (K)");
    ::INPUT::AddNamedReal(m, "PRANUM", "Prandtl number");
    ::INPUT::AddNamedReal(m, "PREEXCON", "pre-exponential constant (1/s)");
    ::INPUT::AddNamedReal(m, "TEMPEXP", "exponent of temperature dependence");
    ::INPUT::AddNamedReal(m, "ACTEMP", "activation temperature (K)");
    ::INPUT::AddNamedReal(m, "UNBSHC", "specific heat capacity of unburnt phase (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "BURSHC", "specific heat capacity of burnt phase (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "UNBTEMP", "temperature of unburnt phase (K)");
    ::INPUT::AddNamedReal(m, "BURTEMP", "temperature of burnt phase (K)");
    ::INPUT::AddNamedReal(m, "UNBDENS", "density of unburnt phase (kg/m^3)");
    ::INPUT::AddNamedReal(m, "BURDENS", "density of burnt phase (kg/m^3)");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material according to Sutherland law with simplified chemical
  // kinetics due to Ferziger and Echekki (1993) (original version and
  // modification by Poinsot and Veynante (2005)) (progress variable)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ferech_pv",
        "material with Ferziger-Echekki (1993) chemical kinetics (progress variable)",
        INPAR::MAT::m_ferech_pv));

    ::INPUT::AddNamedReal(m, "REFVISC", "reference dynamic viscosity (kg/(m*s))");
    ::INPUT::AddNamedReal(m, "REFTEMP", "reference temperature (K)");
    ::INPUT::AddNamedReal(m, "SUTHTEMP", "Sutherland temperature (K)");
    ::INPUT::AddNamedReal(m, "PRANUM", "Prandtl number");
    ::INPUT::AddNamedReal(m, "REACRATECON", "reaction-rate constant (1/s)");
    ::INPUT::AddNamedReal(m, "PVCRIT", "critical value of progress variable");
    ::INPUT::AddNamedReal(m, "UNBSHC", "specific heat capacity of unburnt phase (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "BURSHC", "specific heat capacity of burnt phase (J/(kg*K))");
    ::INPUT::AddNamedReal(m, "UNBTEMP", "temperature of unburnt phase (K)");
    ::INPUT::AddNamedReal(m, "BURTEMP", "temperature of burnt phase (K)");
    ::INPUT::AddNamedReal(m, "UNBDENS", "density of unburnt phase (kg/m^3)");
    ::INPUT::AddNamedReal(m, "BURDENS", "density of burnt phase (kg/m^3)");
    ::INPUT::AddNamedReal(m, "MOD", "modification factor (0.0=original, 1.0=modified)");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrolyte solution (gjb 07/08)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ion",
        "material parameters for ion species in electrolyte solution", INPAR::MAT::m_ion));

    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "kinematic diffusivity");
    ::INPUT::AddNamedReal(m, "VALENCE", "valence (= charge number)");
    ::INPUT::AddNamedReal(m, "DENSIFICATION", "densification coefficient", 0.0, true);
    // via these two optional parameters we can bring the material parameters
    // of one eliminated ionic species into BACI if needed
    ::INPUT::AddNamedReal(
        m, "ELIM_DIFFUSIVITY", "kinematic diffusivity of elim. species", 0.0, true);
    ::INPUT::AddNamedReal(m, "ELIM_VALENCE", "valence of elim. species", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrolyte solution (ehrl 07/12)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_newman",
        "material parameters for ion species in electrolyte solution", INPAR::MAT::m_newman));

    ::INPUT::AddNamedReal(m, "VALENCE", "valence (= charge number)");
    ::INPUT::AddNamedInt(m, "DIFF_COEF_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of diffusion coefficient",
        0);
    ::INPUT::AddNamedInt(m, "DIFF_COEF_TEMP_SCALE_FUNCT",
        "FUNCT number describing temperature scaling of diffusion coefficient", 0);
    ::INPUT::AddNamedInt(m, "TRANSNR", "curve number for transference number");
    ::INPUT::AddNamedInt(m, "THERMFAC", "curve number for thermodynamic factor");
    ::INPUT::AddNamedInt(m, "COND_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of conductivity", 0);
    ::INPUT::AddNamedInt(m, "COND_TEMP_SCALE_FUNCT",
        "FUNCT number describing temperature scaling of conductivity", 0);
    // optional parameter for implemented concentration depending function
    ::INPUT::AddNamedInt(
        m, "DIFF_PARA_NUM", "number of parameters for diffusion coefficient", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "DIFF_PARA", "parameters for diffusion coefficient", "DIFF_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for scaling function describing temperature dependence of diffusion "
        "coefficient",
        0, true);
    ::INPUT::AddNamedRealVector(m, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA",
        "parameters for function describing temperature dependence of diffusion coefficient",
        "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        m, "TRANS_PARA_NUM", "number of parameters for transference number", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "TRANS_PARA", "parameters for transference number", "TRANS_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        m, "THERM_PARA_NUM", "number of parameters for thermodynamic factor", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "THERM_PARA", "parameters for thermodynamic factor", "THERM_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "COND_PARA_NUM", "number of parameters for conductivity", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "COND_PARA", "parameters for conductivity", "COND_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "COND_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for temperature scaling of conductivity", 0, true);
    ::INPUT::AddNamedRealVector(m, "COND_TEMP_SCALE_FUNCT_PARA",
        "parameters for temperature scaling of conductivity", "COND_TEMP_SCALE_FUNCT_PARA_NUM", 0.0,
        true);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrolyte solution for multi-scale approach (fang
  // 07/17)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_newman_multiscale",
        "material parameters for ion species in electrolyte solution for multi-scale approach",
        INPAR::MAT::m_newman_multiscale));

    ::INPUT::AddNamedReal(m, "VALENCE", "valence (= charge number)");
    ::INPUT::AddNamedInt(m, "DIFF_COEF_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of diffusion coefficient",
        0);
    ::INPUT::AddNamedInt(m, "DIFF_COEF_TEMP_SCALE_FUNCT",
        "FUNCT number describing temperature scaling of diffusion coefficient", 0);
    ::INPUT::AddNamedInt(m, "TRANSNR", "curve number for transference number");
    ::INPUT::AddNamedInt(m, "THERMFAC", "curve number for thermodynamic factor");
    ::INPUT::AddNamedInt(m, "COND_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of conductivity", 0);
    ::INPUT::AddNamedInt(m, "COND_TEMP_SCALE_FUNCT",
        "FUNCT number describing temperature scaling of conductivity", 0);
    ::INPUT::AddNamedReal(m, "SIGMA", "electronic conductivity");
    ::INPUT::AddNamedReal(m, "A_s", "specific micro-scale surface area");
    ::INPUT::AddNamedString(m, "MICROFILE", "input file for micro scale", "filename.dat");
    ::INPUT::AddNamedInt(m, "MICRODIS_NUM", "number of micro-scale discretization");
    // optional parameters for implemented concentration-depending functions
    ::INPUT::AddNamedInt(
        m, "DIFF_PARA_NUM", "number of parameters for diffusion coefficient", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "DIFF_PARA", "parameters for diffusion coefficient", "DIFF_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for scaling function describing temperature dependence of diffusion "
        "coefficient",
        0, true);
    ::INPUT::AddNamedRealVector(m, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA",
        "parameters for function describing temperature dependence of diffusion coefficient",
        "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        m, "TRANS_PARA_NUM", "number of parameters for transference number", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "TRANS_PARA", "parameters for transference number", "TRANS_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        m, "THERM_PARA_NUM", "number of parameters for thermodynamic factor", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "THERM_PARA", "parameters for thermodynamic factor", "THERM_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        m, "COND_PARA_NUM", "number of parameters for ionic conductivity", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "COND_PARA", "parameters for ionic conductivity", "COND_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "COND_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for temperature scaling of conductivity", 0, true);
    ::INPUT::AddNamedRealVector(m, "COND_TEMP_SCALE_FUNCT_PARA",
        "parameters for temperature scaling of conductivity", "COND_TEMP_SCALE_FUNCT_PARA_NUM", 0.0,
        true);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_scl", "material parameters for space charge layers", INPAR::MAT::m_scl));

    ::INPUT::AddNamedReal(m, "VALENCE", "valence/charge number");
    ::INPUT::AddNamedInt(m, "DIFF_COEF_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of diffusion coefficient",
        0);
    ::INPUT::AddNamedInt(m, "DIFF_COEF_TEMP_SCALE_FUNCT",
        "function number describing temperature scaling of diffusion coefficient", 0);
    ::INPUT::AddNamedInt(m, "TRANSNR", "curve number for transference number");
    ::INPUT::AddNamedInt(m, "COND_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of conductivity", 0);
    ::INPUT::AddNamedInt(m, "COND_TEMP_SCALE_FUNCT",
        "function number describing temperature scaling of conductivity", 0);
    ::INPUT::AddNamedInt(
        m, "DIFF_PARA_NUM", "number of parameters for diffusion coefficient", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "DIFF_PARA", "parameters for diffusion coefficient", "DIFF_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for scaling function describing temperature dependence of diffusion "
        "coefficient",
        0, true);
    ::INPUT::AddNamedRealVector(m, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA",
        "parameters for function describing temperature dependence of diffusion coefficient",
        "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        m, "TRANS_PARA_NUM", "number of parameters for transference number", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "TRANS_PARA", "parameters for transference number", "TRANS_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "COND_PARA_NUM", "number of parameters for conductivity", 0, true);
    ::INPUT::AddNamedRealVector(
        m, "COND_PARA", "parameters for conductivity", "COND_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(m, "COND_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for temperature scaling of conductivity", 0, true);
    ::INPUT::AddNamedRealVector(m, "COND_TEMP_SCALE_FUNCT_PARA",
        "parameters for temperature scaling of conductivity", "COND_TEMP_SCALE_FUNCT_PARA_NUM", 0.0,
        true);
    ::INPUT::AddNamedReal(m, "MAX_CONC", "maximum cation concentration", 1.0);
    ::INPUT::AddNamedInt(m, "EXTRAPOL_DIFF",
        "strategy for extrapolation of diffusion coefficient below 0 and above MAX_CONC (-1: "
        "disabled, 0: constant)",
        0);
    ::INPUT::AddNamedReal(m, "LIM_CONC", "limiting concentration for extrapolation", 1.0, true);
    ::INPUT::AddNamedReal(m, "BULK_CONC", "bulk ion concentration", 1.0);
    ::INPUT::AddNamedReal(m, "SUSCEPT", "susceptibility", 1.0);
    ::INPUT::AddNamedReal(
        m, "DELTA_NU", "difference of partial molar volumes (vacancy & cation)", 0.0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }


  /*----------------------------------------------------------------------*/
  // electrode material (fang 02/15)
  {
    auto matelectrode = Teuchos::rcp(
        new MaterialDefinition("MAT_electrode", "electrode material", INPAR::MAT::m_electrode));

    // diffusivity and electronic conductivity
    ::INPUT::AddNamedInt(matelectrode, "DIFF_COEF_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of diffusion coefficient",
        0);
    ::INPUT::AddNamedInt(matelectrode, "DIFF_COEF_TEMP_SCALE_FUNCT",
        "FUNCT number describing temperature scaling of diffusion coefficient", 0);
    ::INPUT::AddNamedInt(matelectrode, "COND_CONC_DEP_FUNCT",
        "function number of function describing concentration dependence of conductivity", 0);
    ::INPUT::AddNamedInt(matelectrode, "COND_TEMP_SCALE_FUNCT",
        "FUNCT number describing temperature scaling of conductivity", 0);

    // optional parameters for concentration dependency of diffusivity and electronic conductivity
    ::INPUT::AddNamedInt(
        matelectrode, "DIFF_PARA_NUM", "number of parameters for diffusion coefficient", 0, true);
    ::INPUT::AddNamedRealVector(matelectrode, "DIFF_PARA", "parameters for diffusion coefficient",
        "DIFF_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(matelectrode, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for scaling function describing temperature dependence of diffusion "
        "coefficient",
        0, true);
    ::INPUT::AddNamedRealVector(matelectrode, "DIFF_COEF_TEMP_SCALE_FUNCT_PARA",
        "parameters for function describing temperature dependence of diffusion coefficient",
        "DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(
        matelectrode, "COND_PARA_NUM", "number of parameters for electronic conductivity", 0, true);
    ::INPUT::AddNamedRealVector(matelectrode, "COND_PARA", "parameters for electronic conductivity",
        "COND_PARA_NUM", 0.0, true);
    ::INPUT::AddNamedInt(matelectrode, "COND_TEMP_SCALE_FUNCT_PARA_NUM",
        "number of parameters for temperature scaling of conductivity", 0, true);
    ::INPUT::AddNamedRealVector(matelectrode, "COND_TEMP_SCALE_FUNCT_PARA",
        "parameters for temperature scaling of conductivity", "COND_TEMP_SCALE_FUNCT_PARA_NUM", 0.0,
        true);
    // saturation value of intercalated Lithium concentration
    ::INPUT::AddNamedReal(
        matelectrode, "C_MAX", "saturation value of intercalated Lithium concentration");

    // lithiation value corresponding to saturation value of intercalated Lithium concentration
    ::INPUT::AddNamedReal(matelectrode, "CHI_MAX",
        "lithiation value corresponding to saturation value of intercalated Lithium concentration "
        "'C_MAX'");

    // model for half cell open circuit potential of electrode
    ::INPUT::AddNamedString(matelectrode, "OCP_MODEL",
        "model for half cell open circuit potential of electrode: "
        "Redlich-Kister, Taralov, Polynomial, csv",
        "none");

    // lower bound of range of validity as a fraction of C_MAX for ocp calculation model
    ::INPUT::AddNamedReal(matelectrode, "X_MIN",
        "lower bound of range of validity as a fraction of C_MAX for ocp calculation model", 2.0,
        false);

    // upper bound of range of validity as a fraction of C_MAX for ocp calculation model
    ::INPUT::AddNamedReal(matelectrode, "X_MAX",
        "upper bound of range of validity as a fraction of C_MAX for ocp calculation model", 2.0,
        false);

    // number of parameters underlying half cell open circuit potential model
    ::INPUT::AddNamedInt(matelectrode, "OCP_PARA_NUM",
        "number of parameters underlying half cell open circuit potential model", 0, true);

    // parameters underlying half cell open circuit potential model
    ::INPUT::AddNamedRealVector(matelectrode, "OCP_PARA",
        "parameters underlying half cell open circuit potential model", "OCP_PARA_NUM", 0., true);

    // *.csv file with data points for half cell open circuit potential
    ::INPUT::AddNamedString(matelectrode, "OCP_CSV",
        "\\*.csv file with data points for half cell open circuit potential", "", true);

    // end of input line
    ::INPUT::AddNamedSeparator(matelectrode, "END", "indicating end of line");

    // add electrode material to global list of valid materials
    AppendMaterialDefinition(matlist, matelectrode);
  }

  /*----------------------------------------------------------------------*/
  // material collection (gjb 07/08)
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_matlist", "list/collection of materials, i.e. material IDs", INPAR::MAT::m_matlist));

    ::INPUT::AddNamedBool(
        m, "LOCAL", "individual materials allocated per element or only at global scope");
    // ::INPUT::AddNamedInt(m,"LOCAL","individual materials allocated per element or only at global
    // scope");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material collection with reactions (thon 09/14)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_matlist_reactions",
        "list/collection of materials, i.e. material IDs and list of reactions",
        INPAR::MAT::m_matlist_reactions));

    ::INPUT::AddNamedBool(
        m, "LOCAL", "individual materials allocated per element or only at global scope");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");
    ::INPUT::AddNamedInt(m, "NUMREAC", "number of reactions for these elements", 0);
    ::INPUT::AddNamedIntVector(m, "REACIDS", "advanced reaction list", "NUMREAC", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material collection with chemotaxis (thon 06/15)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_matlist_chemotaxis",
        "list/collection of materials, i.e. material IDs and list of chemotactic pairs",
        INPAR::MAT::m_matlist_chemotaxis));

    ::INPUT::AddNamedBool(
        m, "LOCAL", "individual materials allocated per element or only at global scope");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");
    ::INPUT::AddNamedInt(m, "NUMPAIR", "number of pairs for these elements", 0);
    ::INPUT::AddNamedIntVector(m, "PAIRIDS", "chemotaxis pairs list", "NUMPAIR", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material collection with reactions AND chemotaxis (thon 06/15)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_matlist_chemo_reac",
        "list/collection of materials, i.e. material IDs and list of reactive/chemotactic pairs",
        INPAR::MAT::m_matlist_chemoreac));

    ::INPUT::AddNamedBool(
        m, "LOCAL", "individual materials allocated per element or only at global scope");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");
    ::INPUT::AddNamedInt(m, "NUMPAIR", "number of pairs for these elements", 0);
    ::INPUT::AddNamedIntVector(m, "PAIRIDS", "chemotaxis pairs list", "NUMPAIR", 0);
    ::INPUT::AddNamedInt(m, "NUMREAC", "number of reactions for these elements", 0);
    ::INPUT::AddNamedIntVector(m, "REACIDS", "advanced reaction list", "NUMREAC", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (ehrl 11/12)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_elchmat",
        "specific list/collection of species and phases for elch applications",
        INPAR::MAT::m_elchmat));

    ::INPUT::AddNamedBool(m, "LOCAL",
        "individual materials allocated per element or only at global scope", false, true);
    ::INPUT::AddNamedInt(m, "NUMDOF", "number of dof's per node");
    ::INPUT::AddNamedInt(m, "NUMSCAL", "number of transported scalars per node");
    ::INPUT::AddNamedInt(m, "NUMPHASE", "number of phases in electrolyte");
    ::INPUT::AddNamedIntVector(m, "PHASEIDS", "the list phasel IDs", "NUMPHASE");
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (ehrl 11/12)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_elchphase",
        "material parameters for ion species in electrolyte solution", INPAR::MAT::m_elchphase));

    ::INPUT::AddNamedBool(m, "LOCAL",
        "individual materials allocated per element or only at global scope", false, true);
    ::INPUT::AddNamedReal(m, "EPSILON", "phase porosity");
    ::INPUT::AddNamedReal(m, "TORTUOSITY", "inverse (!) of phase tortuosity");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in electrolyte");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list phasel IDs", "NUMMAT");
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // St.Venant--Kirchhoff
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_Struct_StVenantKirchhoff", "St.Venant--Kirchhoff material", INPAR::MAT::m_stvenant));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // St.Venant--Kirchhoff with temperature
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_ThrStVenantK",
        "Thermo St.Venant--Kirchhoff material", INPAR::MAT::m_thermostvenant));

    ::INPUT::AddNamedInt(m, "YOUNGNUM",
        "number of Young's modulus in list (if 1 Young is const, if >1 Young is temperature) "
        "dependent");
    ::INPUT::AddNamedRealVector(m, "YOUNG", "Young's modulus", "YOUNGNUM");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "THEXPANS", "constant coefficient of linear thermal expansion");
    ::INPUT::AddNamedReal(m, "CAPA", "capacity");
    ::INPUT::AddNamedReal(m, "CONDUCT", "conductivity");
    ::INPUT::AddNamedReal(m, "INITTEMP", "initial temperature");
    ::INPUT::AddNamedInt(m, "THERMOMAT", "mat id of thermal material part", -1, true);

    AppendMaterialDefinition(matlist, m);
  }
  /*----------------------------------------------------------------------*/
  // Plastic linear elastic St.Venant Kirchhoff / Drucker Prager plasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_DruckerPrager",
        "elastic St.Venant Kirchhoff / plastic drucker prager", INPAR::MAT::m_pldruckprag));

    AddNamedReal(m, "YOUNG", "Young's modulus");
    AddNamedReal(m, "NUE", "Poisson's ratio");
    AddNamedReal(m, "DENS", "Density");
    AddNamedReal(m, "ISOHARD", "linear isotropic hardening");
    AddNamedReal(m, "TOL", "Local Newton iteration tolerance");
    AddNamedReal(m, "C", "cohesion");
    AddNamedReal(m, "ETA", "Drucker Prager Constant Eta");
    AddNamedReal(m, "XI", "Drucker Prager Constant Xi");
    AddNamedReal(m, "ETABAR", "Drucker Prager Constant Etabar");
    AddNamedInt(m, "MAXITER", "Maximum Neuton Raphson Iterations", 50, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Linear thermo-elastic St.Venant Kirchhoff / plastic von Mises
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_ThrPlasticLinElast",
        "Thermo-elastic St.Venant Kirchhoff / plastic von Mises material",
        INPAR::MAT::m_thermopllinelast));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "THEXPANS", "coefficient of linear thermal expansion");
    ::INPUT::AddNamedReal(m, "INITTEMP", "initial temperature");
    ::INPUT::AddNamedReal(m, "YIELD", "yield stress");
    ::INPUT::AddNamedReal(m, "ISOHARD", "isotropic hardening modulus");
    ::INPUT::AddNamedReal(m, "KINHARD", "kinematic hardening modulus");
    ::INPUT::AddNamedInt(m, "SAMPLENUM", "number of stress-strain pairs in list");
    ::INPUT::AddNamedRealVector(m, "SIGMA_Y", "yield stress", "SAMPLENUM");
    ::INPUT::AddNamedRealVector(
        m, "EPSBAR_P", "accumulated plastic strain corresponding to SIGMA_Y", "SAMPLENUM");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Finite strain superelasticity of shape memory alloys
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_SuperElastSMA",
        "finite strain superelastic shape memory alloy", INPAR::MAT::m_superelast));

    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "EPSILON_L",
        "parameter representing the maximum deformation obtainable only by detwinning of the "
        "multiple-variant martensite");
    ::INPUT::AddNamedReal(m, "T_AS_s",
        "Temperature at which the phase transformation from austenite to martensite starts");
    ::INPUT::AddNamedReal(m, "T_AS_f",
        "Temperature at which the phase transformation from austenite to martensite finishes");
    ::INPUT::AddNamedReal(m, "T_SA_s",
        "Temperature at which the phase transformation from martensite to autenite starts");
    ::INPUT::AddNamedReal(m, "T_SA_f",
        "Temperature at which the phase transformation from martensite to autenite finishes");
    ::INPUT::AddNamedReal(m, "C_AS", "Coefficient of the linear temperature dependence of T_AS");
    ::INPUT::AddNamedReal(m, "C_SA", "Coefficient of the linear temperature dependence of T_SA");
    ::INPUT::AddNamedReal(m, "SIGMA_AS_s",
        "stress at which the phase transformation from austenite to martensite begins");
    ::INPUT::AddNamedReal(m, "SIGMA_AS_f",
        "stress at which the phase transformation from austenite to martensite finishes");
    ::INPUT::AddNamedReal(m, "SIGMA_SA_s",
        "stress at which the phase transformation from martensite to austenite begins");
    ::INPUT::AddNamedReal(m, "SIGMA_SA_f",
        "stress at which the phase transformation from martensite to austenite finishes");
    ::INPUT::AddNamedReal(m, "ALPHA", "pressure dependency in the drucker-prager-type loading");
    ::INPUT::AddNamedInt(m, "MODEL",
        "Model used for the evolution of martensitic fraction (1=exponential; 2=linear)");
    ::INPUT::AddNamedReal(m, "BETA_AS",
        "parameter, measuring the speed of the transformation from austenite to martensite", 0.,
        true);
    ::INPUT::AddNamedReal(m, "BETA_SA",
        "parameter, measuring the speed of the transformation from martensite to austenite", 0.,
        true);


    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Thermo-hyperelasticity / finite strain von-Mises plasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_ThrPlasticHyperElast",
        "Thermo-hyperelastic / finite strain plastic von Mises material "
        "with linear and exponential isotropic hardening",
        INPAR::MAT::m_thermoplhyperelast));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "CTE", "coefficient of thermal expansion", 0., true);
    ::INPUT::AddNamedReal(m, "INITTEMP", "initial, reference temperature", 0., true);
    ::INPUT::AddNamedReal(m, "YIELD", "initial yield stress");
    ::INPUT::AddNamedReal(m, "ISOHARD", "linear isotropic hardening modulus", 0., true);
    ::INPUT::AddNamedReal(m, "SATHARDENING", "saturation hardening", 0., true);
    ::INPUT::AddNamedReal(m, "HARDEXPO", "hardening exponent", 0., true);
    ::INPUT::AddNamedReal(m, "YIELDSOFT", "thermal yield stress softening", 0., true);
    ::INPUT::AddNamedReal(m, "HARDSOFT",
        "thermal hardening softening (acting on SATHARDENING and ISOHARD)", 0., true);
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration", 1.e-8, true);

    AppendMaterialDefinition(matlist, m);
  }


  /*----------------------------------------------------------------------*/
  // Hyperelasticity / finite strain von-Mises plasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_PlasticNlnLogNeoHooke",
        "hyperelastic / finite strain plastic von Mises material "
        "with linear and exponential isotropic hardening",
        INPAR::MAT::m_plnlnlogneohooke));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "YIELD", "yield stress");
    ::INPUT::AddNamedReal(m, "ISOHARD", "isotropic hardening modulus");
    ::INPUT::AddNamedReal(m, "SATHARDENING", "saturation hardening");
    ::INPUT::AddNamedReal(m, "HARDEXPO", "linear hardening exponent");
    ::INPUT::AddNamedReal(m, "VISC", "VISCOSITY", 0., true);
    ::INPUT::AddNamedReal(m, "RATE_DEPENDENCY", "rate dependency", 0., true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Plastic linear elastic St.Venant Kirchhoff / von Mises
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_PlasticLinElast",
        "elastic St.Venant Kirchhoff / plastic von Mises material "
        "with linear isotropic and kineamtic hardening",
        INPAR::MAT::m_pllinelast));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "YIELD", "yield stress");
    ::INPUT::AddNamedReal(m, "ISOHARD", "linear isotropic hardening modulus");
    ::INPUT::AddNamedReal(m, "KINHARD", "linear kinematic hardening modulus");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Elastic visco-plastic finite strain material law without yield surface
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Viscoplastic_No_Yield_Surface",
        "Elastic visco-plastic finite strain material law without yield surface",
        INPAR::MAT::m_vp_no_yield_surface));

    // elasticity parameters
    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "material mass density");
    // visco-plasticity parameters
    ::INPUT::AddNamedReal(m, "TEMPERATURE", "temperature in Kelvin");
    ::INPUT::AddNamedReal(
        m, "PRE_EXP_FAC", "pre-exponential factor of plastic shear strain rate 'A'");
    ::INPUT::AddNamedReal(m, "ACTIVATION_ENERGY", "activation energy 'Q'");
    ::INPUT::AddNamedReal(m, "GAS_CONSTANT", "gas constant 'R'");
    ::INPUT::AddNamedReal(m, "STRAIN_RATE_SENS", "strain-rate-sensitivity 'm'");
    ::INPUT::AddNamedReal(m, "INIT_FLOW_RES", "initial isotropic flow resistance 'S^0'");
    ::INPUT::AddNamedReal(m, "FLOW_RES_PRE_FAC", "flow resistance factor 'H_0'");
    ::INPUT::AddNamedReal(m, "FLOW_RES_EXP", "flow resistance exponential value 'a'");
    ::INPUT::AddNamedReal(m, "FLOW_RES_SAT_FAC", "flow resistance saturation factor 'S_*'");
    ::INPUT::AddNamedReal(m, "FLOW_RES_SAT_EXP", "flow resistance saturation exponent 'b'");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Robinson's visco-plastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_Struct_Robinson", "Robinson's visco-plastic material", INPAR::MAT::m_vp_robinson));

    ::INPUT::AddNamedString(m, "KIND",
        "kind of Robinson material: "
        "Butler, Arya, Arya_NarloyZ (default), Arya_CrMoSteel",
        "Arya_NarloyZ");
    ::INPUT::AddNamedInt(m, "YOUNGNUM", "number of Young's modulus in list");
    ::INPUT::AddNamedRealVector(m, "YOUNG", "Young's modulus", "YOUNGNUM");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "THEXPANS", "coefficient of linear thermal expansion");
    ::INPUT::AddNamedReal(m, "INITTEMP", "initial temperature");
    ::INPUT::AddNamedReal(m, "HRDN_FACT", "hardening factor 'A'");
    ::INPUT::AddNamedReal(m, "HRDN_EXPO", "hardening power 'n'");
    ::INPUT::AddNamedInt(m, "SHRTHRSHLDNUM", "number of shear stress threshold 'K^2'in list");
    ::INPUT::AddNamedRealVector(
        m, "SHRTHRSHLD", "Bingam-Prager shear stress threshold 'K^2'", "SHRTHRSHLDNUM");
    ::INPUT::AddNamedReal(m, "RCVRY", "recovery factor 'R_0'");
    ::INPUT::AddNamedReal(m, "ACTV_ERGY", "activation energy 'Q_0'");
    ::INPUT::AddNamedReal(m, "ACTV_TMPR", "activation temperature 'T_0'");
    ::INPUT::AddNamedReal(m, "G0", "'G_0'");
    ::INPUT::AddNamedReal(m, "M_EXPO", "'m'");
    ::INPUT::AddNamedInt(m, "BETANUM", "number of 'beta' in list");
    ::INPUT::AddNamedRealVector(m, "BETA", "beta", "BETANUM");
    ::INPUT::AddNamedReal(m, "H_FACT", "'H'");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Elasto-plastic material with damage, based on MAT_Struct_PlasticLinElast
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Damage",
        "elasto-plastic von Mises material with ductile damage", INPAR::MAT::m_elpldamage));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedInt(m, "SAMPLENUM", "number of stress-strain pairs in list");
    ::INPUT::AddNamedRealVector(m, "SIGMA_Y", "yield stress", "SAMPLENUM");
    ::INPUT::AddNamedRealVector(
        m, "EPSBAR_P", "accumulated plastic strain corresponding to SIGMA_Y", "SAMPLENUM");
    ::INPUT::AddNamedReal(m, "DAMDEN", "denominator of damage evoluation law");
    ::INPUT::AddNamedReal(m, "DAMEXP", "exponent of damage evoluation law");
    ::INPUT::AddNamedReal(m, "DAMTHRESHOLD", "damage threshold");
    ::INPUT::AddNamedReal(m, "KINHARD", "kinematic hardening modulus, stress-like variable");
    ::INPUT::AddNamedReal(m, "KINHARD_REC", "recovery factor, scalar-valued variable");
    ::INPUT::AddNamedReal(m, "SATHARDENING", "saturation hardening");
    ::INPUT::AddNamedReal(m, "HARDEXPO", "hardening exponent");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // aneurysm wall material according to Raghavan and Vorp [2000]
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAANeoHooke",
        "aneurysm wall material according to Raghavan and Vorp [2000]", INPAR::MAT::m_aaaneohooke));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "BETA", "2nd parameter");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // aneurysm wall material according to Raghavan and Vorp [2000]
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAANeoHookeStopro",
        "aneurysm wall material according to Raghavan and Vorp [2000] with stochastic "
        "modelling of beta",
        INPAR::MAT::m_aaaneohooke_stopro));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "BETA", "2nd parameter");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    // Stochastic properties are set via randomfield class

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // AAA thrombus material according to GASSER et. al. [2008]
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAAGasser",
        "AAA thrombus material according to GASSER [2008]", INPAR::MAT::m_aaagasser));

    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedString(
        m, "VOL", "Type of volumetric Strain Energy Density: OSM (default),SuBa,SiTa", "OSM");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio (0.49)");
    ::INPUT::AddNamedReal(m, "BETA", "empiric constant for OSM (-2.0)");
    ::INPUT::AddNamedReal(m, "CLUM", "luminal stiffness parameter (2.62e3)");
    ::INPUT::AddNamedReal(m, "CMED", "medial stiffness parameter (2.62e3)");
    ::INPUT::AddNamedReal(m, "CABLUM", "abluminal stiffness parameter (2.62e3)");

    /*
     ::INPUT::AddNamedReal(m,"DENS","mass density");
     ::INPUT::AddNamedReal(m,"KAPPA","dilatation modulus");
     ::INPUT::AddNamedReal(m,"BETA","empiric constant");
     ::INPUT::AddNamedReal(m,"CLUM","luminal stiffness parameter");
     ::INPUT::AddNamedReal(m,"CMED","medial stiffness parameter");
     ::INPUT::AddNamedReal(m,"CABLUM","abluminal stiffness parameter");
     */

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // aneurysm wall material according to Raghavan and Vorp [2000] with damage Simo
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Raghavan_Damage",
        "aneurysm wall material according to Raghavan and Vorp [2000] with damage",
        INPAR::MAT::m_aaaraghavanvorp_damage));

    ::INPUT::AddNamedReal(m, "BULK", "Bulk's modulus");
    ::INPUT::AddNamedReal(m, "ALPHA", "1nd parameter,alpha");
    ::INPUT::AddNamedReal(m, "BETA", "2nd parameter,beta");
    ::INPUT::AddNamedReal(m, "EQSTRMIN", "equivalent strain initial damage");
    ::INPUT::AddNamedReal(m, "A", "1st parameter, a");
    ::INPUT::AddNamedReal(m, "B", "2nd parameter, b");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // aneurysm wall material SEF according  to Raghavan and Vorp [2000],
  // parameters according to mixed effects model
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_AAA_MixedEffects",
        "aneurysm wall material according to Mixed Effects Model", INPAR::MAT::m_aaa_mixedeffects));

    ::INPUT::AddNamedReal(m, "AGE", "age");
    ::INPUT::AddNamedReal(m, "REFDIA", "subrenal diameter");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Visco-elastic Neo-Hookean material law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_VISCONEOHOOKE",
        "visco-elastic neo-Hookean material law", INPAR::MAT::m_visconeohooke));
    ::INPUT::AddNamedReal(m, "YOUNGS_SLOW", "???");
    ::INPUT::AddNamedReal(m, "POISSON", "???");
    ::INPUT::AddNamedReal(m, "DENS", "???");
    ::INPUT::AddNamedReal(m, "YOUNGS_FAST", "???");
    ::INPUT::AddNamedReal(m, "RELAX", "???");
    ::INPUT::AddNamedReal(m, "THETA", "???");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Visco-elastic anisotropic fiber material law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_VISCOANISO",
        "visco-elastic anisotropic fibre material law", INPAR::MAT::m_viscoanisotropic));

    ::INPUT::AddNamedReal(m, "KAPPA", "dilatation modulus");
    ::INPUT::AddNamedReal(m, "MUE", "Shear Modulus");
    ::INPUT::AddNamedReal(m, "DENS", "Density");
    ::INPUT::AddNamedReal(m, "K1", "Parameter for linear fiber stiffness");
    ::INPUT::AddNamedReal(m, "K2", "Parameter for exponetial fiber stiffness");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle between fibers");
    ::INPUT::AddNamedReal(m, "BETA_ISO", "ratio between elasticities in generalized Maxweel body");
    ::INPUT::AddNamedReal(
        m, "BETA_ANISO", "ratio between elasticities in generalized Maxweel body");
    ::INPUT::AddNamedReal(m, "RELAX_ISO", "isotropic relaxation time");
    ::INPUT::AddNamedReal(m, "RELAX_ANISO", "anisotropic relaxation time");
    ::INPUT::AddNamedReal(m, "MINSTRETCH", "minimal principal stretch fibers do respond to");
    ::INPUT::AddNamedInt(
        m, "ELETHICKDIR", "Element thickness direction applies also to fibers (only sosh)");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Structural micro-scale approach: material parameters are calculated from microscale simulation
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Multiscale",
        "Structural micro-scale approach: material parameters are calculated from microscale "
        "simulation",
        INPAR::MAT::m_struct_multiscale));

    ::INPUT::AddNamedString(m, "MICROFILE", "inputfile for microstructure", "filename.dat");
    ::INPUT::AddNamedInt(m, "MICRODIS_NUM", "Number of microscale discretization");
    ::INPUT::AddNamedReal(m, "INITVOL", "Initial volume of RVE", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ElastHyper",
        "list/collection of hyperelastic materials, i.e. material IDs", INPAR::MAT::m_elasthyper));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials/potentials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material/potential IDs", "NUMMAT");
    ::INPUT::AddNamedReal(m, "DENS", "material mass density");
    ::INPUT::AddNamedInt(m, "POLYCONVEX", "1.0 if polyconvexity of system is checked", 0., true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // viscohyperelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ViscoElastHyper",
        "Viscohyperelastic material compatible with the collection of hyperelastic materials",
        INPAR::MAT::m_viscoelasthyper));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials/potentials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material/potential IDs", "NUMMAT");
    ::INPUT::AddNamedReal(m, "DENS", "material mass density");
    ::INPUT::AddNamedInt(m, "POLYCONVEX", "1.0 if polyconvexity of system is checked", 0., true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials for finite strain plasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PlasticElastHyper",
        "list/collection of hyperelastic materials, i.e. material IDs",
        INPAR::MAT::m_plelasthyper));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials/potentials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material/potential IDs", "NUMMAT");
    ::INPUT::AddNamedReal(m, "DENS", "material mass density");
    ::INPUT::AddNamedReal(m, "INITYIELD", "initial yield stress");
    ::INPUT::AddNamedInt(m, "POLYCONVEX", "1.0 if polyconvexity of system is checked", 0., true);
    ::INPUT::AddNamedReal(m, "ISOHARD", "linear isotropic hardening modulus", 0., true);
    ::INPUT::AddNamedReal(m, "EXPISOHARD", "nonlinear isotropic hardening exponent", 0., true);
    ::INPUT::AddNamedReal(
        m, "INFYIELD", "saturation yield stress for nonlinear isotropic hardening", 0., true);
    ::INPUT::AddNamedReal(m, "KINHARD", "linear kinematic hardening modulus", 0., true);

    // visco-plasticity
    ::INPUT::AddNamedReal(m, "VISC", "Visco-Plasticity parameter 'eta' in Perzyna model", 0., true);
    ::INPUT::AddNamedReal(
        m, "RATE_DEPENDENCY", "Visco-Plasticity parameter 'eta' in Perzyna model", 1., true);
    ::INPUT::AddNamedReal(m, "VISC_SOFT",
        "Visco-Plasticity temperature dependency (eta = eta_0 * (1-(T-T_0)*x)", 0., true);

    // optional pastic spin parameter
    ::INPUT::AddNamedReal(
        m, "PL_SPIN_CHI", "Plastic spin coupling parameter chi (often called eta)", 0.0, true);

    // optional Hill yield parameters
    ::INPUT::AddNamedReal(
        m, "rY_11", "relative yield stress in fiber1-direction (Y_11/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_22", "relative yield stress in fiber2-direction (Y_22/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_33", "relative yield stress in fiber3-direction (Y_33/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_12", "relative shear yield stress in 12-direction (Y_12/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_23", "relative shear yield stress in 23-direction (Y_23/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_13", "relative shear yield stress in 13-direction (Y_13/Y_0)", 0.0, true);

    // optional TSI parameters
    ::INPUT::AddNamedReal(m, "CTE", "coefficient of thermal expansion", 0., true);
    ::INPUT::AddNamedReal(m, "INITTEMP", "initial, reference temperature", 0., true);
    ::INPUT::AddNamedReal(m, "YIELDSOFT", "yield stress softening", 0., true);
    ::INPUT::AddNamedReal(m, "HARDSOFT", "hardening softening", 0., true);
    ::INPUT::AddNamedReal(
        m, "TAYLOR_QUINNEY", "Taylor-Quinney factor for plastic heat conversion", 1., true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials for finite strain plasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PlasticElastHyperVCU",
        "list/collection of hyperelastic materials, i.e. material IDs",
        INPAR::MAT::m_plelasthyperVCU));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials/potentials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material/potential IDs", "NUMMAT");
    ::INPUT::AddNamedReal(m, "DENS", "material mass density");
    ::INPUT::AddNamedReal(m, "INITYIELD", "initial yield stress");
    ::INPUT::AddNamedReal(m, "ISOHARD", "linear isotropic hardening modulus", 0., true);
    ::INPUT::AddNamedReal(m, "EXPISOHARD", "nonlinear isotropic hardening exponent", 0., true);
    ::INPUT::AddNamedReal(
        m, "INFYIELD", "saturation yield stress for nonlinear isotropic hardening", 0., true);
    ::INPUT::AddNamedReal(m, "KINHARD", "linear kinematic hardening modulus", 0., true);

    // visco-plasticity
    ::INPUT::AddNamedReal(m, "VISC", "Visco-Plasticity parameter 'eta' in Perzyna model", 0., true);
    ::INPUT::AddNamedReal(
        m, "RATE_DEPENDENCY", "Visco-Plasticity parameter 'eta' in Perzyna model", 1., true);
    ::INPUT::AddNamedReal(m, "VISC_SOFT",
        "Visco-Plasticity temperature dependency (eta = eta_0 * (1-(T-T_0)*x)", 0., true);

    // optional pastic spin parameter
    ::INPUT::AddNamedReal(
        m, "PL_SPIN_CHI", "Plastic spin coupling parameter chi (often called eta)", 0.0, true);

    // optional Hill yield parameters
    ::INPUT::AddNamedReal(
        m, "rY_11", "relative yield stress in fiber1-direction (Y_11/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_22", "relative yield stress in fiber2-direction (Y_22/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_33", "relative yield stress in fiber3-direction (Y_33/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_12", "relative shear yield stress in 12-direction (Y_12/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_23", "relative shear yield stress in 23-direction (Y_23/Y_0)", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "rY_13", "relative shear yield stress in 13-direction (Y_13/Y_0)", 0.0, true);

    // optional TSI parameters
    ::INPUT::AddNamedReal(m, "CTE", "coefficient of thermal expansion", 0., true);
    ::INPUT::AddNamedReal(m, "INITTEMP", "initial, reference temperature", 0., true);
    ::INPUT::AddNamedReal(m, "YIELDSOFT", "yield stress softening", 0., true);
    ::INPUT::AddNamedReal(m, "HARDSOFT", "hardening softening", 0., true);
    ::INPUT::AddNamedReal(
        m, "TAYLOR_QUINNEY", "Taylor-Quinney factor for plastic heat conversion", 1., true);

    ::INPUT::AddNamedInt(m, "POLYCONVEX", "1.0 if polyconvexity of system is checked", 0., true);


    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // logarithmic neo-Hooke material acc. to Bonet and Wood
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupLogNeoHooke",
        "logarithmic neo-Hooke material acc. to Bonet and Wood", INPAR::MAT::mes_couplogneohooke));

    ::INPUT::AddNamedString(m, "MODE",
        "parameter set: YN (Young's modulus and Poisson's ration; default) or Lame (mue and "
        "lambda)",
        "YN");
    ::INPUT::AddNamedReal(m, "C1", "E or mue");
    ::INPUT::AddNamedReal(m, "C2", "nue or lambda");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // Saint-Venant-Kirchhoff as elastic summand
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_CoupSVK", "Saint-Venant-Kirchhoff as elastic summand", INPAR::MAT::mes_coupSVK));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // Simo-Pister type material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_CoupSimoPister", "Simo-Pister type material", INPAR::MAT::mes_coupsimopister));

    ::INPUT::AddNamedReal(m, "MUE", "material constant");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // logarithmic mixed neo-Hooke material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupLogMixNeoHooke",
        "mixed logarithmic neo-Hooke material", INPAR::MAT::mes_couplogmixneohooke));

    ::INPUT::AddNamedString(m, "MODE",
        "parameter set: YN (Young's modulus and Poisson's ration; default) or Lame (mue and "
        "lambda)",
        "YN");
    ::INPUT::AddNamedReal(m, "C1", "E or mue");
    ::INPUT::AddNamedReal(m, "C2", "nue or lambda");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled exponential material for compressible material (according to Weikenmeier_2014)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupExpPol",
        "compressible, isochoric exponential material law for soft tissue",
        INPAR::MAT::mes_coupexppol));
    ::INPUT::AddNamedReal(m, "A", "material constant");
    ::INPUT::AddNamedReal(m, "B", "material constant linear I_1");
    ::INPUT::AddNamedReal(m, "C", "material constant linear J");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // compressible neo-Hooke material acc. to Holzapfel
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupNeoHooke",
        "compressible neo-Hooke material acc. to Holzapfel", INPAR::MAT::mes_coupneohooke));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus", 0.0, true);
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }
  // Mooney Rivlin  material acc. to Holzapfel
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupMooneyRivlin",
        "Mooney - Rivlin material acc. to Holzapfel", INPAR::MAT::mes_coupmooneyrivlin));

    ::INPUT::AddNamedReal(m, "C1", "material constant", 0.0, true);
    ::INPUT::AddNamedReal(m, "C2", "material constant", 0.0, true);
    ::INPUT::AddNamedReal(m, "C3", "material constant", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled Blatz and Ko material acc. to Holzapfel
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupBlatzKo",
        "Blatz and Ko material acc. to Holzapfel", INPAR::MAT::mes_coupblatzko));

    ::INPUT::AddNamedReal(m, "MUE", "Shear modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "F", "interpolation parameter");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of Neo-Hooke
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoNeoHooke",
        "isochoric part of neo-Hooke material acc. to Holzapfel", INPAR::MAT::mes_isoneohooke));

    ::INPUT::AddNamedReal(m, "MUE", "Shear modulus");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of one-term Ogden material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoOgden",
        "isochoric part of the one-term Ogden material", INPAR::MAT::mes_isoogden));

    ::INPUT::AddNamedReal(m, "MUE", "Shear modulus");
    ::INPUT::AddNamedReal(m, "ALPHA", "Nonlinearity parameter");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric and volumetric contribution of AAAGasser
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoVolAAAGasser",
        "isochoric and volumetric part of AAAGasser material (thrombus)",
        INPAR::MAT::mes_isovolaaagasser));

    ::INPUT::AddNamedReal(m, "CLUM", "luminal stiffness parameter (2.62e3)");
    ::INPUT::AddNamedReal(m, "CMED", "medial stiffness parameter (2.62e3)");
    ::INPUT::AddNamedReal(m, "CABLUM", "abluminal stiffness parameter (2.62e3)");
    ::INPUT::AddNamedReal(m, "NUE", "");
    ::INPUT::AddNamedReal(m, "BETA", "");
    // optional parameters for uncertainty quantification
    ::INPUT::AddNamedReal(
        m, "MULUM", "mu for luminal pdf, irrelevant for deterministic analysis", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "MUMED", "mu for medial pdf, irrelevant for deterministic analysis", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "MUABLUM", "mu for abluminal pdf, irrelevant for deterministic analysis", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "SIGMALUM", "std for luminal pdf, irrelevant for deterministic analysis", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "SIGMAMED", "std for medial pdf, irrelevant for deterministic analysis", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "SIGMAABLUM", "std for abluminal pdf, irrelevant for deterministic analysis", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of Yeoh
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoYeoh",
        "isochoric part of  Yeoh material acc. to Holzapfel", INPAR::MAT::mes_isoyeoh));

    ::INPUT::AddNamedReal(m, "C1", "Linear modulus");
    ::INPUT::AddNamedReal(m, "C2", "Quadratic modulus");
    ::INPUT::AddNamedReal(m, "C3", "Cubic modulus");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of iso1pow
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_Iso1Pow", "isochoric part of general power material", INPAR::MAT::mes_iso1pow));

    ::INPUT::AddNamedReal(m, "C", "material parameter");
    ::INPUT::AddNamedInt(m, "D", "exponent");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of iso2pow
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_Iso2Pow", "isochoric part of general power material", INPAR::MAT::mes_iso2pow));

    ::INPUT::AddNamedReal(m, "C", "material parameter");
    ::INPUT::AddNamedInt(m, "D", "exponent");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // contribution of coup1pow
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_Coup1Pow", "part of general power material", INPAR::MAT::mes_coup1pow));

    ::INPUT::AddNamedReal(m, "C", "material parameter");
    ::INPUT::AddNamedInt(m, "D", "exponent");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // contribution of coup2pow
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_Coup2Pow", "part of general power material", INPAR::MAT::mes_coup2pow));

    ::INPUT::AddNamedReal(m, "C", "material parameter");
    ::INPUT::AddNamedInt(m, "D", "exponent");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // contribution of coup3pow
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_Coup3Pow", "part of general power material", INPAR::MAT::mes_coup3pow));

    ::INPUT::AddNamedReal(m, "C", "material parameter");
    ::INPUT::AddNamedInt(m, "D", "exponent");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // contribution of coup13apow
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_Coup13aPow",
        "hyperelastic potential summand for multiplicative coupled invariants I1 and I3",
        INPAR::MAT::mes_coup13apow));

    ::INPUT::AddNamedReal(m, "C", "material parameter");
    ::INPUT::AddNamedInt(m, "D", "exponent of all");
    ::INPUT::AddNamedReal(m, "A", "negative exponent of I3");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of expo
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoExpoPow",
        "isochoric part of  exponential material acc. to Holzapfel", INPAR::MAT::mes_isoexpopow));

    ::INPUT::AddNamedReal(m, "K1", "material parameter");
    ::INPUT::AddNamedReal(m, "K2", "material parameter");
    ::INPUT::AddNamedInt(m, "C", "exponent");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric contribution of mooney rivlin
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoMooneyRivlin",
        "isochoric part of  Mooney-Rivlin material acc. to Holzapfel",
        INPAR::MAT::mes_isomooneyrivlin));

    ::INPUT::AddNamedReal(m, "C1", "Linear modulus for first invariant");
    ::INPUT::AddNamedReal(m, "C2", "Linear modulus for second invariant");
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric anisotropic material with one exponential fiber family
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoMuscle_Blemker",
        "anisotropic Blemker muscle material", INPAR::MAT::mes_isomuscleblemker));

    ::INPUT::AddNamedReal(m, "G1", "muscle along fiber shear modulus");
    ::INPUT::AddNamedReal(m, "G2", "muscle cross fiber shear modulus");
    ::INPUT::AddNamedReal(m, "P1", "linear material parameter for passive along-fiber response");
    ::INPUT::AddNamedReal(
        m, "P2", "exponential material parameter for passive along-fiber response");
    ::INPUT::AddNamedReal(m, "SIGMAMAX", "maximal active isometric stress");
    ::INPUT::AddNamedReal(m, "LAMBDAOFL", "optimal fiber stretch");
    ::INPUT::AddNamedReal(
        m, "LAMBDASTAR", "stretch at which the normalized passive fiber force becomes linear");
    ::INPUT::AddNamedReal(m, "ALPHA", "tetanised activation level,");
    ::INPUT::AddNamedReal(m, "BETA", "constant scaling tanh-type activation function");
    ::INPUT::AddNamedReal(m, "ACTSTARTTIME", "starting time of muscle activation");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // test material to test elasthyper-toolbox
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoTestMaterial",
        "test material to test elasthyper-toolbox", INPAR::MAT::mes_isotestmaterial));

    ::INPUT::AddNamedReal(m, "C1", "Modulus for first invariant");
    ::INPUT::AddNamedReal(m, "C2", "Modulus for second invariant");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // general fiber material for remodeling
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_RemodelFiber",
        "General fiber material for remodeling", INPAR::MAT::mes_remodelfiber));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials/potentials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material/potential IDs", "NUMMAT");
    ::INPUT::AddNamedReal(m, "TDECAY", "decay time of Poisson (degradation) process");
    ::INPUT::AddNamedReal(m, "GROWTHFAC", "time constant for collagen growth", 0.0, true);
    ::INPUT::AddNamedRealVector(m, "COLMASSFRAC",
        "initial mass fraction of first collagen fiber family in constraint mixture", "NUMMAT", 0.0,
        true);
    ::INPUT::AddNamedReal(m, "DEPOSITIONSTRETCH", "deposition stretch");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // volumetric contribution of Sussman Bathe
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_VolSussmanBathe",
        "volumetric part of  SussmanBathe material", INPAR::MAT::mes_volsussmanbathe));

    ::INPUT::AddNamedReal(m, "KAPPA", "dilatation modulus");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // volumetric penalty contribution
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_VolPenalty",
        "Penalty formulation for the volumetric part", INPAR::MAT::mes_volpenalty));

    ::INPUT::AddNamedReal(m, "EPSILON", "penalty parameter");
    ::INPUT::AddNamedReal(m, "GAMMA", "penalty parameter");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // volumetric contribution of Ogden
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_VolOgden", "Ogden formulation for the volumetric part", INPAR::MAT::mes_vologden));

    ::INPUT::AddNamedReal(m, "KAPPA", "dilatation modulus");
    ::INPUT::AddNamedReal(m, "BETA", "empiric constant");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // volumetric power law contribution
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_VolPow", "Power law formulation for the volumetric part", INPAR::MAT::mes_volpow));

    ::INPUT::AddNamedReal(m, "A", "prefactor of power law");
    ::INPUT::AddNamedReal(m, "EXPON", "exponent of power law");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with one exponential fiber family
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpoActive",
        "anisotropic active fiber", INPAR::MAT::mes_coupanisoexpoactive));

    ::INPUT::AddNamedReal(m, "K1", "linear constant");
    ::INPUT::AddNamedReal(m, "K2", "exponential constant");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle");
    ::INPUT::AddNamedReal(m, "K1COMP", "linear constant");
    ::INPUT::AddNamedReal(m, "K2COMP", "exponential constant");
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);
    ::INPUT::AddNamedReal(m, "S", "maximum contractile stress");
    ::INPUT::AddNamedReal(m, "LAMBDAMAX", "stretch at maximum active force generation");
    ::INPUT::AddNamedReal(m, "LAMBDA0", "stretch at zero active force generation");
    ::INPUT::AddNamedReal(m, "DENS", "total reference mass density of constrained mixture");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with one exponential fiber family
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpo",
        "anisotropic part with one exp. fiber", INPAR::MAT::mes_coupanisoexpo));

    ::INPUT::AddNamedReal(m, "K1", "linear constant");
    ::INPUT::AddNamedReal(m, "K2", "exponential constant");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle");
    ::INPUT::AddNamedReal(m, "K1COMP", "linear constant");
    ::INPUT::AddNamedReal(m, "K2COMP", "exponential constant");
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);
    ::INPUT::AddNamedInt(
        m, "FIBER_ID", "Id of the fiber to be used (1 for first fiber, default)", 1, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with one exponential shear behavior between two fibers
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpoShear",
        "Exponential shear behavior between two fibers", INPAR::MAT::mes_coupanisoexposhear));

    ::INPUT::AddNamedReal(m, "K1", "linear constant");
    ::INPUT::AddNamedReal(m, "K2", "exponential constant");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle");
    ::INPUT::AddNamedReal(m, "K1COMP", "linear constant");
    ::INPUT::AddNamedReal(m, "K2COMP", "exponential constant");
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedIntVector(m, "FIBER_IDS",
        "Ids of the two fibers to be used (1 for the first fiber, 2 for the second, default)", 2);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with one pow-like fiber family
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoPow",
        "anisotropic part with one pow-like fiber", INPAR::MAT::mes_coupanisopow));

    ::INPUT::AddNamedReal(m, "K", "linear constant");
    ::INPUT::AddNamedReal(m, "D1", "exponential constant for fiber invariant");
    ::INPUT::AddNamedReal(m, "D2", "exponential constant for system");
    ::INPUT::AddNamedReal(m, "ACTIVETHRES",
        "Deformation threshold for activating fibers. Default:"
        " 1.0 (off at compression); If 0.0 (always active)",
        1.0, true);
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(
        m, "FIBER", "Number of the fiber family contained in the element", 1, true);
    ::INPUT::AddNamedReal(m, "GAMMA", "angle", 0.0, true);
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with two exponential fiber families
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpoTwoCoup",
        "anisotropic part with two exp. fibers", INPAR::MAT::mes_coupanisoexpotwocoup));

    ::INPUT::AddNamedReal(m, "A4", "linear anisotropic constant for fiber 1");
    ::INPUT::AddNamedReal(m, "B4", "exponential anisotropic constant for fiber 1");
    ::INPUT::AddNamedReal(m, "A6", "linear anisotropic constant for fiber 2");
    ::INPUT::AddNamedReal(m, "B6", "exponential anisotropic constant for fiber 2");
    ::INPUT::AddNamedReal(m, "A8", "linear anisotropic constant for fiber 1 relating fiber 2");
    ::INPUT::AddNamedReal(m, "B8", "exponential anisotropic constant for fiber 1 relating fiber 2");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle");
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(
        m, "FIB_COMP", "fibers support compression: yes (true) or no (false)", true, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with two exponential fiber families
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoNeoHooke",
        "anisotropic part with one neo Hookean fiber", INPAR::MAT::mes_coupanisoneohooke));

    ::INPUT::AddNamedReal(m, "C", "linear constant");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle");
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with the stress given by a simplified version of the contraction
  // law of Bestel-Clement-Sorine
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_AnisoActiveStress_Evolution",
        "anisotropic part with one fiber with coefficient given by a simplification of the "
        "activation-contraction law of Bestel-Clement-Sorine-2001",
        INPAR::MAT::mes_anisoactivestress_evolution));

    ::INPUT::AddNamedReal(m, "SIGMA", "Contractility (maximal stress)");
    ::INPUT::AddNamedReal(m, "TAUC0", "Initial value for the active stress");
    ::INPUT::AddNamedReal(m, "MAX_ACTIVATION", "Maximal value for the rescaled activation");
    ::INPUT::AddNamedReal(m, "MIN_ACTIVATION", "Minimal value for the rescaled activation");
    ::INPUT::AddNamedInt(
        m, "SOURCE_ACTIVATION", "Where the activation comes from: 0=scatra , >0 Id for FUNCT");
    ::INPUT::AddNamedReal(m, "ACTIVATION_THRES",
        "Threshold for activation (contraction starts when activation function is larger than this "
        "value, relaxes otherwise)");
    ::INPUT::AddNamedBool(m, "STRAIN_DEPENDENCY",
        "model strain dependency of contractility (Frank-Starling law): no (false) or yes (true)",
        false, true);
    ::INPUT::AddNamedReal(
        m, "LAMBDA_LOWER", "lower fiber stretch for Frank-Starling law", 1.0, true);
    ::INPUT::AddNamedReal(
        m, "LAMBDA_UPPER", "upper fiber stretch for Frank-Starling law", 1.0, true);
    ::INPUT::AddNamedReal(m, "GAMMA", "angle", 0.0, true);
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization mode for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled anisotropic material with variable stress coefficient
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoNeoHooke_VarProp",
        "anisotropic part with one neo Hookean fiber with variable coefficient",
        INPAR::MAT::mes_coupanisoneohooke_varprop));

    ::INPUT::AddNamedReal(m, "C", "linear constant");
    ::INPUT::AddNamedInt(
        m, "SOURCE_ACTIVATION", "Where the activation comes from: 0=scatra , >0 Id for FUNCT");
    ::INPUT::AddNamedReal(m, "GAMMA", "azimuth angle", 0.0, true);
    ::INPUT::AddNamedReal(m, "THETA", "polar angle", 0.0, true);
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization mode for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric anisotropic material with one exponential fiber family
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_IsoAnisoExpo",
        "anisotropic part with one exp. fiber", INPAR::MAT::mes_isoanisoexpo));

    ::INPUT::AddNamedReal(m, "K1", "linear constant");
    ::INPUT::AddNamedReal(m, "K2", "exponential constant");
    ::INPUT::AddNamedReal(m, "GAMMA", "angle");
    ::INPUT::AddNamedReal(m, "K1COMP", "linear constant");
    ::INPUT::AddNamedReal(m, "K2COMP", "exponential constant");
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);
    ::INPUT::AddNamedBool(m, "ADAPT_ANGLE", "adapt angle during remodeling", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // structural tensor
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_StructuralTensor",
        "Parameter for structural tensor strategy in anisotropic materials",
        INPAR::MAT::mes_structuraltensorstratgy));

    ::INPUT::AddNamedString(m, "STRATEGY",
        "Strategy for evaluation of structural tensor: "
        "Standard (default), ByDistributionFunction, DispersedTransverselyIsotropic",
        "Standard");

    // choose between:
    // "none"
    // "Bingham"
    // "vonMisesFisher"
    //  rauch 10/17
    ::INPUT::AddNamedString(m, "DISTR",
        "Type of distribution function around mean direction: "
        "none, Bingham, vonMisesFisher",
        "none", true);

    ::INPUT::AddNamedReal(m, "C1", "constant 1 for distribution function", 1.0, true);
    ::INPUT::AddNamedReal(m, "C2", "constant 2 for distribution function", 0.0, true);
    ::INPUT::AddNamedReal(m, "C3", "constant 3 for distribution function", 0.0, true);
    ::INPUT::AddNamedReal(m, "C4", "constant 4 for distribution function", 1e16, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // transversely isotropic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("ELAST_CoupTransverselyIsotropic",
        "transversely part of a simple othotropic, transversely "
        "isotropic hyperelastic constitutive equation",
        INPAR::MAT::mes_couptransverselyisotropic));

    ::INPUT::AddNamedReal(m, "ALPHA", "1-st constant");
    ::INPUT::AddNamedReal(m, "BETA", "2-nd constant");
    ::INPUT::AddNamedReal(m, "GAMMA", "3-rd constant");
    ::INPUT::AddNamedReal(m, "ANGLE", "fiber angle");
    ::INPUT::AddNamedInt(m, "STR_TENS_ID", "MAT ID for definition of Structural Tensor");
    ::INPUT::AddNamedInt(m, "FIBER", "exponential constant", 1, true);
    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for fiber alignment", 1, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // coupled Varga material acc. to Holzapfel
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_CoupVarga", "Varga material acc. to Holzapfel", INPAR::MAT::mes_coupvarga));

    ::INPUT::AddNamedReal(m, "MUE", "Shear modulus");
    ::INPUT::AddNamedReal(m, "BETA", "'Anti-modulus'");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric Varga material acc. to Holzapfel
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "ELAST_IsoVarga", "Isochoric Varga material acc. to Holzapfel", INPAR::MAT::mes_isovarga));

    ::INPUT::AddNamedReal(m, "MUE", "Shear modulus");
    ::INPUT::AddNamedReal(m, "BETA", "'Anti-modulus'");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isotropic viscous contribution of myocardial matrix (chapelle12)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("VISCO_CoupMyocard",
        "Isotropic viscous contribution of myocardial matrix", INPAR::MAT::mes_coupmyocard));

    ::INPUT::AddNamedReal(m, "N", "material parameter");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // isochoric rate dependent viscos material, modified from Pioletti,1997
  {
    auto m = Teuchos::rcp(new MaterialDefinition("VISCO_IsoRateDep",
        "Isochoric rate dependent viscous material", INPAR::MAT::mes_isoratedep));

    ::INPUT::AddNamedReal(m, "N", "material parameter");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // viscos contribution to visohyperelastic material according to SLS-Model
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "VISCO_GenMax", "Viscous contribution according to SLS-Model", INPAR::MAT::mes_genmax));

    ::INPUT::AddNamedReal(m, "TAU", "relaxation parameter");
    ::INPUT::AddNamedReal(m, "BETA", "emphasis of viscous to elastic part");
    ::INPUT::AddNamedString(m, "SOLVE",
        "Solution of evolution equation via: OST (default) or CONVOL (convolution integral)",
        "OST");


    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // viscos contribution to visohyperelastic material according to FSLS-Model
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "VISCO_Fract", "Viscous contribution according to FSLS-Model", INPAR::MAT::mes_fract));

    ::INPUT::AddNamedReal(m, "TAU", "relaxation parameter");
    ::INPUT::AddNamedReal(m, "ALPHA", "fractional order derivative");
    ::INPUT::AddNamedReal(m, "BETA", "emphasis of viscous to elastic part");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // viscous contribution of a branch of a generalized Maxwell model
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "VISCO_PART", "Viscous contribution of a viscoelastic Branch", INPAR::MAT::mes_viscopart));

    ::INPUT::AddNamedReal(m, "TAU", "dynamic viscosity divided by young's modulus of the branch");

    AppendMaterialDefinition(matlist, m);
  }
  /*--------------------------------------------------------------------*/
  // viscoelatic branches of a generalized Maxwell model
  {
    auto m = Teuchos::rcp(new MaterialDefinition("VISCO_GeneralizedGenMax",
        "Viscoelastic Branches of generalized Maxwell", INPAR::MAT::mes_generalizedgenmax));

    ::INPUT::AddNamedInt(m, "NUMBRANCH", "number of viscoelastic branches");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMBRANCH");
    ::INPUT::AddNamedString(m, "SOLVE",
        "Solution for evolution equation: OST (default) or CONVOL (convolution integral)",
        "CONVOL");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // description of a viscoelatic branch of a generalized Maxwell model
  {
    auto m = Teuchos::rcp(new MaterialDefinition("VISCO_BRANCH",
        "Viscoelastic Branch (viscous and elastic contribution)", INPAR::MAT::mes_viscobranch));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in the viscoelastic branch");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // 1D Artery material with constant properties
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_CNST_ART", "artery with constant properties", INPAR::MAT::m_cnst_art));

    ::INPUT::AddNamedReal(m, "VISCOSITY",
        "viscosity (for CONSTANT viscosity law taken as blood viscosity, for BLOOD viscosity law "
        "taken as the viscosity of blood plasma)");
    ::INPUT::AddNamedReal(m, "DENS", "density of blood");
    ::INPUT::AddNamedReal(m, "YOUNG", "artery Youngs modulus of elasticity");
    ::INPUT::AddNamedReal(m, "NUE", "Poissons ratio of artery fiber");
    ::INPUT::AddNamedReal(m, "TH", "artery thickness");
    ::INPUT::AddNamedReal(m, "PEXT1", "artery fixed external pressure 1");
    ::INPUT::AddNamedReal(m, "PEXT2", "artery fixed external pressure 2");
    ::INPUT::AddNamedString(
        m, "VISCOSITYLAW", "type of viscosity law, CONSTANT (default) or BLOOD", "CONSTANT", true);
    ::INPUT::AddNamedReal(m, "BLOOD_VISC_SCALE_DIAM_TO_MICRONS",
        "used to scale the diameter for blood viscosity law to microns if your problem is not "
        "given in microns, e.g., if you use mms, set this parameter to 1.0e3",
        1.0, true);
    ::INPUT::AddNamedString(m, "VARYING_DIAMETERLAW",
        "type of varying diameter law, CONSTANT (default) or BY_FUNCTION", "CONSTANT", true);
    ::INPUT::AddNamedInt(
        m, "VARYING_DIAMETER_FUNCTION", "function for varying diameter law", -1, true);
    ::INPUT::AddNamedReal(m, "COLLAPSE_THRESHOLD",
        "Collapse threshold for diameter (below this diameter element is assumed to be collapsed "
        "with zero diameter and is not evaluated)",
        -1.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // Fourier's law
  {
    auto m = Teuchos::rcp(new MaterialDefinition("THERM_FourierIso",
        "isotropic (linear) Fourier's law of heat conduction", INPAR::MAT::m_th_fourier_iso));

    ::INPUT::AddNamedReal(m, "CAPA", "volumetric heat capacity");
    ::INPUT::AddNamedReal(m, "CONDUCT", "thermal conductivity");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // material for heat transport due to Fourier-type thermal conduction and the Soret effect (fang
  // 06/15)
  {
    auto matsoret = Teuchos::rcp(new MaterialDefinition("MAT_soret",
        "material for heat transport due to Fourier-type thermal conduction and the Soret effect",
        INPAR::MAT::m_soret));

    // mandatory parameters
    ::INPUT::AddNamedReal(matsoret, "CAPA", "volumetric heat capacity");
    ::INPUT::AddNamedReal(matsoret, "CONDUCT", "thermal conductivity");
    ::INPUT::AddNamedReal(matsoret, "SORET", "Soret coefficient");

    // add Soret material to global list of valid materials
    AppendMaterialDefinition(matlist, matsoret);
  }

  /*----------------------------------------------------------------------*/
  // integration point based growth
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_GrowthVolumetric", "volumetric growth", INPAR::MAT::m_growth_volumetric));

    ::INPUT::AddNamedInt(m, "GROWTHLAW", "number of growth law in input file");
    ::INPUT::AddNamedInt(
        m, "IDMATELASTIC", "number of elastic material in input file: MAT IDMATELASTIC ...");
    ::INPUT::AddNamedReal(m, "STARTTIME", "start growth after this time");
    ::INPUT::AddNamedReal(m, "ENDTIME", "end growth after this time");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // collection of hyperelastic materials for membranes
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Membrane_ElastHyper",
        "list/collection of hyperelastic materials for membranes, i.e. material IDs",
        INPAR::MAT::m_membrane_elasthyper));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials/potentials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material/potential IDs", "NUMMAT");
    ::INPUT::AddNamedReal(m, "DENS", "material mass density");
    ::INPUT::AddNamedInt(m, "POLYCONVEX", "1.0 if polyconvexity of system is checked", 0., true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // active strain membrane material for gastric electromechanics
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_Membrane_ActiveStrain",
        "active strain membrane material", INPAR::MAT::m_membrane_activestrain));

    ::INPUT::AddNamedInt(m, "MATIDPASSIVE", "MATID for the passive material", false);
    ::INPUT::AddNamedInt(
        m, "SCALIDVOLTAGE", "ID of the scalar that represents the (SMC) voltage", false);
    ::INPUT::AddNamedReal(m, "DENS", "material mass density", false);
    ::INPUT::AddNamedReal(m, "BETA1", "Ca2+ dynamics", false);
    ::INPUT::AddNamedReal(m, "BETA2", "opening dynamics of the VDCC", false);
    ::INPUT::AddNamedReal(m, "VOLTHRESH", "voltage threshold for activation", false);
    ::INPUT::AddNamedReal(m, "ALPHA1", "intensity of contraction in fiber direction 1", false);
    ::INPUT::AddNamedReal(m, "ALPHA2", "intensity of contraction in fiber direction 2", false);
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // growth and remodeling (homogenized constrained mixture model)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthRemodel_ElastHyper",
        "growth and remodeling", INPAR::MAT::m_growthremodel_elasthyper));

    ::INPUT::AddNamedInt(m, "NUMMATRF", "number of remodelfiber materials in list", false);
    ::INPUT::AddNamedInt(
        m, "NUMMATEL3D", "number of 3d elastin matrix materials/potentials in list", 0, true);
    ::INPUT::AddNamedInt(
        m, "NUMMATEL2D", "number of 2d elastin matrix materials/potentials in list", false);
    ::INPUT::AddNamedIntVector(
        m, "MATIDSRF", "the list remodelfiber material IDs", "NUMMATRF", false);
    ::INPUT::AddNamedIntVector(m, "MATIDSEL3D", "the list 3d elastin matrix material/potential IDs",
        "NUMMATEL3D", -1, true);
    ::INPUT::AddNamedIntVector(
        m, "MATIDSEL2D", "the list 2d elastin matrix material/potential IDs", "NUMMATEL2D", false);
    ::INPUT::AddNamedInt(m, "MATIDELPENALTY", "penalty material ID", -1, true);
    ::INPUT::AddNamedReal(
        m, "ELMASSFRAC", "initial mass fraction of elastin matrix in constraint mixture", false);
    ::INPUT::AddNamedReal(m, "DENS", "material mass density", false);
    ::INPUT::AddNamedReal(
        m, "PRESTRETCHELASTINCIR", "circumferential prestretch of elastin matrix", false);
    ::INPUT::AddNamedReal(m, "PRESTRETCHELASTINAX", "axial prestretch of elastin matrix", false);
    ::INPUT::AddNamedReal(m, "THICKNESS",
        "reference wall thickness of the idealized cylindrical aneurysm [m]", -1, true);
    ::INPUT::AddNamedReal(m, "MEANPRESSURE", "mean blood pressure [Pa]", -1.0, true);
    ::INPUT::AddNamedReal(
        m, "RADIUS", "inner radius of the idealized cylindrical aneurysm [m]", -1.0, true);
    ::INPUT::AddNamedInt(
        m, "DAMAGE", "1: elastin damage after prestressing,0: no elastin damage", false);
    ::INPUT::AddNamedInt(m, "GROWTHTYPE",
        "flag to decide what type of collagen growth is used: 1: anisotropic growth; 0: isotropic "
        "growth",
        false);
    ::INPUT::AddNamedInt(m, "LOCTIMEINT",
        "flag to decide what type of local time integration scheme is used: 1: Backward Euler "
        "Method; 0: Forward Euler Method",
        false);
    ::INPUT::AddNamedInt(m, "MEMBRANE",
        "Flag whether Hex or Membrane elements are used ( Membrane: 1, Hex: Everything else )", -1,
        true);
    ::INPUT::AddNamedInt(m, "CYLINDER",
        "Flag that geometry is a cylinder. 1: aligned in x-direction; 2: y-direction; 3: "
        "z-direction",
        -1, true);
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // multiplicative split of deformation gradient in elastic and inelastic parts
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_MultiplicativeSplitDefgradElastHyper",
        "multiplicative split of deformation gradient",
        INPAR::MAT::m_multiplicative_split_defgrad_elasthyper));

    ::INPUT::AddNamedInt(m, "NUMMATEL", "number of elastic materials/potentials in list", 0, false);
    ::INPUT::AddNamedIntVector(
        m, "MATIDSEL", "the list of elastic material/potential IDs", "NUMMATEL", -1, false);
    ::INPUT::AddNamedInt(
        m, "NUMFACINEL", "number of factors of inelastic deformation gradient", false);
    ::INPUT::AddNamedIntVector(m, "INELDEFGRADFACIDS",
        "the list of inelastic deformation gradient factor IDs", "NUMFACINEL", false);
    ::INPUT::AddNamedReal(m, "DENS", "material mass density", false);
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // simple inelastic material law featuring no volume change
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradNoGrowth",
        "no volume change, i.e. the inelastic deformation gradient is the identity tensor",
        INPAR::MAT::mfi_no_growth));

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // simple isotropic, volumetric growth; growth is linearly dependent on scalar mapped to material
  // configuration, constant material density
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradLinScalarIso",
        "scalar dependent isotropic growth law; volume change linearly dependent on scalar (in "
        "material configuration)",
        INPAR::MAT::mfi_lin_scalar_iso));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of growth inducing scalar");
    ::INPUT::AddNamedReal(
        m, "SCALAR1_MolarGrowthFac", "isotropic molar growth factor due to scalar 1");
    ::INPUT::AddNamedReal(
        m, "SCALAR1_RefConc", "reference concentration of scalar 1 causing no strains");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // simple anisotropic, volumetric growth; growth direction prescribed in input-file;
  // growth is linearly dependent on scalar mapped to material configuration, constant material
  // density
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradLinScalarAniso",
        "scalar dependent anisotropic growth law; growth in direction as given in input-file; "
        "volume change linearly dependent on scalar (in material configuration)",
        INPAR::MAT::mfi_lin_scalar_aniso));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of growth inducing scalar");
    ::INPUT::AddNamedReal(
        m, "SCALAR1_MolarGrowthFac", "anisotropic molar growth factor due to scalar 1");
    ::INPUT::AddNamedReal(
        m, "SCALAR1_RefConc", "reference concentration of scalar 1 causing no strains");
    ::INPUT::AddNamedInt(m, "NUMSPACEDIM", "Number of space dimension (only 3 valid)");
    ::INPUT::AddNamedRealVector(
        m, "GrowthDirection", "vector that defines the growth direction", "NUMSPACEDIM");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // non-linear isotropic volumetric growth; growth is dependent on the degree of lithiation,
  // constant material density, nonlinear behavior prescribed by polynomial in input file
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradPolyIntercalFracIso",
        "scalar dependent isotropic growth law; volume change nonlinearly dependent on the "
        "intercalation fraction, that is calculated using the scalar concentration (in material "
        "configuration)",
        INPAR::MAT::mfi_poly_intercal_frac_iso));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of growth inducing scalar");
    ::INPUT::AddNamedReal(
        m, "SCALAR1_RefConc", "reference concentration of scalar 1 causing no strains");
    ::INPUT::AddNamedInt(m, "POLY_PARA_NUM", "number of polynomial coefficients");
    ::INPUT::AddNamedRealVector(m, "POLY_PARAMS", "coefficients of polynomial", "POLY_PARA_NUM");
    ::INPUT::AddNamedReal(m, "X_min", "lower bound of validity of polynomial");
    ::INPUT::AddNamedReal(m, "X_max", "upper bound of validity of polynomial");
    ::INPUT::AddNamedInt(m, "MATID", "material ID of the corresponding scatra material");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // non-linear anisotropic volumetric growth; growth direction prescribed in input-file;
  // growth is dependent on the degree of lithiation, constant material density, nonlinear behavior
  // prescribed by polynomial in input file
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradPolyIntercalFracAniso",
        "scalar dependent anisotropic growth law; growth in direction as given in input-file; "
        "volume change nonlinearly dependent on the intercalation fraction, that is calculated "
        "using the scalar concentration (in material configuration)",
        INPAR::MAT::mfi_poly_intercal_frac_aniso));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of growth inducing scalar");
    ::INPUT::AddNamedReal(
        m, "SCALAR1_RefConc", "reference concentration of scalar 1 causing no strains");
    ::INPUT::AddNamedInt(m, "NUMSPACEDIM", "Number of space dimension (only 3 valid)");
    ::INPUT::AddNamedRealVector(
        m, "GrowthDirection", "vector that defines the growth direction", "NUMSPACEDIM");
    ::INPUT::AddNamedInt(m, "POLY_PARA_NUM", "number of polynomial coefficients");
    ::INPUT::AddNamedRealVector(m, "POLY_PARAMS", "coefficients of polynomial", "POLY_PARA_NUM");
    ::INPUT::AddNamedReal(m, "X_min", "lower bound of validity of polynomial");
    ::INPUT::AddNamedReal(m, "X_max", "upper bound of validity of polynomial");
    ::INPUT::AddNamedInt(m, "MATID", "material ID of the corresponding scatra material");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradLinTempIso",
        "Temperature dependent growth law. Volume change linearly dependent on temperature",
        INPAR::MAT::mfi_lin_temp_iso));

    ::INPUT::AddNamedReal(m, "Temp_GrowthFac", "isotropic growth factor due to temperature");
    ::INPUT::AddNamedReal(m, "RefTemp", "reference temperature causing no strains");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_InelasticDefgradTimeFunct",
        "Time-dependent growth law. Determinant of volume change dependent on time function "
        "defined "
        "by 'FUNCT_NUM",
        INPAR::MAT::mfi_time_funct));

    ::INPUT::AddNamedInt(m, "FUNCT_NUM",
        "Time-dependent function of the determinant of the inelastic deformation gradient");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // integration point based and scalar dependent interpolation between to materials
  {
    auto mm = Teuchos::rcp(new MaterialDefinition("MAT_ScDepInterp",
        "integration point based and scalar dependent interpolation between to materials",
        INPAR::MAT::m_sc_dep_interp));

    ::INPUT::AddNamedInt(mm, "IDMATZEROSC", "material for lambda equal to zero");
    ::INPUT::AddNamedInt(mm, "IDMATUNITSC", "material for lambda equal to one");
    //      ::INPUT::AddNamedReal(mm,"ALPHA","size of ",-1.0,true);

    AppendMaterialDefinition(matlist, mm);
  }

  /*----------------------------------------------------------------------*/
  // anisotropic strain-dependent growth law (Goektepe et al., J Theor Biol 2010, Lee et al., BMMB
  // 2017)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthAnisoStrain",
        "growth law depending on elastic stretch in fiber direction, growth in fiber direction",
        INPAR::MAT::m_growth_aniso_strain));

    ::INPUT::AddNamedReal(m, "TAU", "growth time scale");
    ::INPUT::AddNamedReal(m, "TAU_REV", "reverse growth time scale");
    ::INPUT::AddNamedReal(m, "THETA_MIN", "lower limit for growth stretch");
    ::INPUT::AddNamedReal(m, "THETA_MAX", "upper limit for growth stretch");
    ::INPUT::AddNamedReal(m, "GAMMA", "growth non-linearity");
    ::INPUT::AddNamedReal(m, "GAMMA_REV", "reverse growth non-linearity");
    ::INPUT::AddNamedReal(m, "LAMBDA_CRIT", "critical fiber stretch");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // anisotropic strain-dependent growth law (Goektepe et al., J Theor Biol 2010, Lee et al., BMMB
  // 2017)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthAnisoStress",
        "growth law depending on elastic Mandel stress, growth perpendicular to fiber direction",
        INPAR::MAT::m_growth_aniso_stress));

    ::INPUT::AddNamedReal(m, "TAU", "growth time scale");
    ::INPUT::AddNamedReal(m, "TAU_REV", "reverse growth time scale");
    ::INPUT::AddNamedReal(m, "THETA_MIN", "lower limit for growth stretch");
    ::INPUT::AddNamedReal(m, "THETA_MAX", "upper limit for growth stretch");
    ::INPUT::AddNamedReal(m, "GAMMA", "growth non-linearity");
    ::INPUT::AddNamedReal(m, "GAMMA_REV", "reverse growth non-linearity");
    ::INPUT::AddNamedReal(m, "P_CRIT", "critical pressure");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // anisotropic strain-dependent growth law with constant prescribed trigger (for multiscale in
  // time)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthAnisoStrainConstTrig",
        "growth law depending on prescribed constant elastic stretch in fiber direction, "
        "growth in fiber direction",
        INPAR::MAT::m_growth_aniso_strain_const_trig));

    ::INPUT::AddNamedReal(m, "TAU", "growth time scale");
    ::INPUT::AddNamedReal(m, "TAU_REV", "reverse growth time scale");
    ::INPUT::AddNamedReal(m, "THETA_MIN", "lower limit for growth stretch");
    ::INPUT::AddNamedReal(m, "THETA_MAX", "upper limit for growth stretch");
    ::INPUT::AddNamedReal(m, "GAMMA", "growth non-linearity");
    ::INPUT::AddNamedReal(m, "GAMMA_REV", "reverse growth non-linearity");
    ::INPUT::AddNamedReal(m, "LAMBDA_CRIT", "critical fiber stretch");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // anisotropic strain-dependent growth law with constant prescribed trigger (for multiscale in
  // time)
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthAnisoStressConstTrig",
        "growth law depending on prescribed constant elastic Mandel stress, growth "
        "perpendicular to fiber direction",
        INPAR::MAT::m_growth_aniso_stress_const_trig));

    ::INPUT::AddNamedReal(m, "TAU", "growth time scale");
    ::INPUT::AddNamedReal(m, "TAU_REV", "reverse growth time scale");
    ::INPUT::AddNamedReal(m, "THETA_MIN", "lower limit for growth stretch");
    ::INPUT::AddNamedReal(m, "THETA_MAX", "upper limit for growth stretch");
    ::INPUT::AddNamedReal(m, "GAMMA", "growth non-linearity");
    ::INPUT::AddNamedReal(m, "GAMMA_REV", "reverse growth non-linearity");
    ::INPUT::AddNamedReal(m, "P_CRIT", "critical pressure");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // isotropic growth law (cf. Diss Tinkl 2015, LNM)
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_GrowthIsoStress", "stress-dependent growth law", INPAR::MAT::m_growth_iso_stress));

    ::INPUT::AddNamedReal(m, "THETAPLUS", "maximal growth stretch");
    ::INPUT::AddNamedReal(m, "KPLUS", "growth law parameter kthetaplus");
    ::INPUT::AddNamedReal(m, "MPLUS", "growth law parameter mthetaplus");
    ::INPUT::AddNamedReal(m, "THETAMINUS", "minimal growth stretch");
    ::INPUT::AddNamedReal(m, "KMINUS", "growth law parameter kthetaminus");
    ::INPUT::AddNamedReal(m, "MMINUS", "growth law parameter mthetaminus");
    ::INPUT::AddNamedReal(m, "HOMMANDEL", "homeostatic value for mandelstress");
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for local Newton iteration");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // simple atherosclerosis growth law, scalar-dependent volumetric growth
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_GrowthAC", "scalar depended volumetric growth", INPAR::MAT::m_growth_ac));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of first growth inducing scalar");
    ::INPUT::AddNamedReal(m, "ALPHA", "volume per first scalar's mass density");
    ::INPUT::AddNamedInt(m, "SCALAR2", "number of second growth inducing scalar", 1, true);
    ::INPUT::AddNamedReal(m, "BETA", "volume per second scalar's mass density", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // atherosclerosis growth law, scalar depended growth in radial direction
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthACRadial",
        "scalar depended growth in radial direction", INPAR::MAT::m_growth_ac_radial));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of first growth inducing scalar");
    ::INPUT::AddNamedReal(m, "ALPHA", "volume per first scalar's mass density");
    ::INPUT::AddNamedInt(m, "SCALAR2", "number of second growth inducing scalar", 1, true);
    ::INPUT::AddNamedReal(m, "BETA", "volume per second scalar's mass density", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // atherosclerosis growth law, scalar depended growth in radial direction
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_GrowthACRadialRefConc",
        "scalar depended growth in radial direction", INPAR::MAT::m_growth_ac_radial_refconc));

    ::INPUT::AddNamedInt(m, "SCALAR1", "number of first growth inducing scalar");
    ::INPUT::AddNamedReal(m, "ALPHA", "volume per first scalar's mass density");
    ::INPUT::AddNamedInt(m, "SCALAR2", "number of second growth inducing scalar", 1, true);
    ::INPUT::AddNamedReal(m, "BETA", "volume per second scalar's mass density", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // constant rate growth law
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_GrowthConst", "constant growth law", INPAR::MAT::m_growth_const));

    ::INPUT::AddNamedReal(m, "THETARATE", "reference value for mandelstress");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // growth and remodeling of arteries
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ConstraintMixture",
        "growth and remodeling of arteries", INPAR::MAT::m_constraintmixture));

    ::INPUT::AddNamedReal(m, "DENS", "Density");
    ::INPUT::AddNamedReal(m, "MUE", "Shear Modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "PHIE", "mass fraction of elastin");
    ::INPUT::AddNamedReal(m, "PREELA", "prestretch of elastin");
    ::INPUT::AddNamedReal(m, "K1", "Parameter for linear collagen fiber stiffness");
    ::INPUT::AddNamedReal(m, "K2", "Parameter for exponential collagen fiber stiffness");
    ::INPUT::AddNamedInt(m, "NUMHOM", "Number of homeostatic parameters", 1);
    ::INPUT::AddNamedRealVector(m, "PRECOLL", "prestretch of collagen fibers", "NUMHOM");
    ::INPUT::AddNamedReal(m, "DAMAGE", "damage stretch of collagen fibers");
    ::INPUT::AddNamedReal(m, "K1M", "Parameter for linear smooth muscle fiber stiffness");
    ::INPUT::AddNamedReal(m, "K2M", "Parameter for exponential smooth muscle fiber stiffness");
    ::INPUT::AddNamedReal(m, "PHIM", "mass fraction of smooth muscle");
    ::INPUT::AddNamedReal(m, "PREMUS", "prestretch of smooth muscle fibers");
    ::INPUT::AddNamedReal(m, "SMAX", "maximal active stress");
    ::INPUT::AddNamedReal(m, "KAPPA", "dilatation modulus");
    ::INPUT::AddNamedReal(m, "LIFETIME", "lifetime of collagen fibers");
    ::INPUT::AddNamedReal(m, "GROWTHFAC", "growth factor for stress");
    ::INPUT::AddNamedRealVector(
        m, "HOMSTR", "homeostatic target value of scalar stress measure", "NUMHOM");
    ::INPUT::AddNamedReal(m, "SHEARGROWTHFAC", "growth factor for shear");
    ::INPUT::AddNamedReal(m, "HOMRAD", "homeostatic target value of inner radius");
    ::INPUT::AddNamedReal(m, "STARTTIME", "at this time turnover of collagen starts");
    ::INPUT::AddNamedString(m, "INTEGRATION",
        "time integration scheme: "
        "Explicit (default), or Implicit",
        "Explicit");
    ::INPUT::AddNamedReal(
        m, "TOL", "tolerance for local Newton iteration, only for implicit integration");
    ::INPUT::AddNamedString(m, "GROWTHFORCE",
        "driving force of growth: "
        "Single (default), All, ElaCol",
        "Single");
    ::INPUT::AddNamedString(m, "ELASTINDEGRAD",
        "how elastin is degraded: "
        "None (default), Rectangle, Time",
        "None");
    ::INPUT::AddNamedString(m, "MASSPROD",
        "how mass depends on driving force: "
        "Lin (default), CosCos",
        "Lin");
    ::INPUT::AddNamedString(m, "INITSTRETCH",
        "how to set stretches in the beginning (None, Homeo, UpdatePrestretch)", "None");
    ::INPUT::AddNamedInt(m, "CURVE", "number of timecurve for increase of prestretch in time", 0);
    ::INPUT::AddNamedString(m, "DEGOPTION",
        "Type of degradation function: "
        "Lin (default), Cos, Exp, ExpVar",
        "Lin");
    ::INPUT::AddNamedReal(m, "MAXMASSPRODFAC", "maximal factor of mass production");
    ::INPUT::AddNamedReal(m, "ELASTINFAC", "factor for elastin content", 0.0, true);
    ::INPUT::AddNamedBool(m, "STOREHISTORY",
        "store all history variables, not recommended for forward simulations", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // hyperelastic material for poroelasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_StructPoro", "wrapper for structure poroelastic material", INPAR::MAT::m_structporo));

    ::INPUT::AddNamedInt(m, "MATID", "ID of structure material");
    ::INPUT::AddNamedInt(m, "POROLAWID", "ID of porosity law");
    ::INPUT::AddNamedReal(m, "INITPOROSITY", "initial porosity of porous medium");

    AppendMaterialDefinition(matlist, m);
  }
  /*----------------------------------------------------------------------*/
  // linear law for porosity in porous media problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawLinear",
        "linear constitutive law for porosity", INPAR::MAT::m_poro_law_linear));

    ::INPUT::AddNamedReal(m, "BULKMODULUS", "bulk modulus of porous medium");

    AppendMaterialDefinition(matlist, m);
  }
  /*----------------------------------------------------------------------*/
  // constant law for porosity in porous media problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawConstant",
        "constant constitutive law for porosity", INPAR::MAT::m_poro_law_constant));

    AppendMaterialDefinition(matlist, m);
  }
  /*----------------------------------------------------------------------*/
  // neo-hookean law for porosity in porous media problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawNeoHooke",
        "NeoHookean-like constitutive law for porosity",
        INPAR::MAT::m_poro_law_logNeoHooke_Penalty));

    ::INPUT::AddNamedReal(m, "BULKMODULUS", "bulk modulus of porous medium");
    ::INPUT::AddNamedReal(m, "PENALTYPARAMETER", "penalty paramter of porous medium");

    AppendMaterialDefinition(matlist, m);
  }
  /*----------------------------------------------------------------------*/
  // incompressible skeleton law for porosity in porous media problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawIncompSkel",
        "porosity law for incompressible skeleton phase", INPAR::MAT::m_poro_law_incompr_skeleton));

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // incompressible skeleton law for porosity in porous media problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawLinBiot",
        "linear biot model for porosity law", INPAR::MAT::m_poro_law_linear_biot));

    ::INPUT::AddNamedReal(m, "INVBIOTMODULUS", "inverse Biot modulus of porous medium");
    ::INPUT::AddNamedReal(m, "BIOTCEOFF", "Biot coefficient of porous medium");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // incompressible skeleton law for porosity depending on the density
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroLawDensityDependent",
        "porosity depending on the density", INPAR::MAT::m_poro_law_density_dependent));

    ::INPUT::AddNamedInt(m, "DENSITYLAWID", "material ID of density law");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // density law for constant density in porous multiphase medium
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroDensityLawConstant",
        "density law for constant density in porous multiphase medium",
        INPAR::MAT::m_poro_densitylaw_constant));

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // density law for constant density in porous multiphase medium
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PoroDensityLawExp",
        "density law for pressure dependent exponential function",
        INPAR::MAT::m_poro_densitylaw_exp));

    ::INPUT::AddNamedReal(m, "BULKMODULUS", "bulk modulus of porous medium");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // permeability law for constant permeability in porous multiphase medium
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroRelPermeabilityLawConstant",
        "permeability law for constant permeability in porous multiphase medium",
        INPAR::MAT::m_fluidporo_relpermeabilitylaw_constant));

    ::INPUT::AddNamedReal(m, "VALUE", "constant value of permeability");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // permeability law for permeability depending on saturation according to (saturation)^exp
  // in porous multiphase medium
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroRelPermeabilityLawExp",
        "permeability law depending on saturation in porous multiphase medium",
        INPAR::MAT::m_fluidporo_relpermeabilitylaw_exp));

    ::INPUT::AddNamedReal(m, "EXP", "exponent of the saturation of this phase");
    ::INPUT::AddNamedReal(m, "MIN_SAT", "minimum saturation which is used for calculation");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // viscosity law for constant viscosity in porous multiphase medium
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroViscosityLawConstant",
        "viscosity law for constant viscosity in porous multiphase medium",
        INPAR::MAT::m_fluidporo_viscositylaw_constant));

    ::INPUT::AddNamedReal(m, "VALUE", "constant value of viscosity");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // viscosity law for viscosity-dependency modelling cell adherence
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroViscosityLawCellAdherence",
        "visosity law depending on pressure gradient in porous multiphase medium",
        INPAR::MAT::m_fluidporo_viscositylaw_celladh));

    ::INPUT::AddNamedReal(m, "VISC_0", "Visc0 parameter for modelling cell adherence");
    ::INPUT::AddNamedReal(m, "XI", "xi parameter for modelling cell adherence");
    ::INPUT::AddNamedReal(m, "PSI", "psi parameter for modelling cell adherence");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // hyperelastic material for poroelasticity with reaction
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_StructPoroReaction",
        "wrapper for structure porelastic material with reaction",
        INPAR::MAT::m_structpororeaction));

    ::INPUT::AddNamedInt(m, "MATID", "ID of structure material");
    ::INPUT::AddNamedInt(m, "POROLAWID", "ID of porosity law");
    ::INPUT::AddNamedReal(m, "INITPOROSITY", "initial porosity of porous medium");
    ::INPUT::AddNamedInt(m, "DOFIDREACSCALAR",
        "Id of DOF within scalar transport problem, which controls the reaction");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // hyperelastic material for poroelasticity with reaction
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_StructPoroReactionECM",
        "wrapper for structure porelastic material with reaction",
        INPAR::MAT::m_structpororeactionECM));

    ::INPUT::AddNamedInt(m, "MATID", "ID of structure material");
    ::INPUT::AddNamedInt(m, "POROLAWID", "ID of porosity law");
    ::INPUT::AddNamedReal(m, "INITPOROSITY", "initial porosity of porous medium");
    ::INPUT::AddNamedReal(m, "DENSCOLLAGEN", "density of collagen");
    ::INPUT::AddNamedInt(m, "DOFIDREACSCALAR",
        "Id of DOF within scalar transport problem, which controls the reaction");
    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // fluid flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_FluidPoro", "flow in deformable porous media", INPAR::MAT::m_fluidporo));

    ::INPUT::AddNamedReal(m, "DYNVISCOSITY", "dynamic viscosity");
    ::INPUT::AddNamedReal(m, "DENSITY", "density");
    ::INPUT::AddNamedReal(m, "PERMEABILITY", "permeability of medium", 0.0, true);
    ::INPUT::AddNamedReal(
        m, "AXIALPERMEABILITY", "axial permeability for transverse isotropy", 0.0, true);
    ::INPUT::AddNamedReal(m, "ORTHOPERMEABILITY1", "first permeability for orthotropy", 0.0, true);
    ::INPUT::AddNamedReal(m, "ORTHOPERMEABILITY2", "second permeability for orthotropy", 0.0, true);
    ::INPUT::AddNamedReal(m, "ORTHOPERMEABILITY3", "third permeability for orthotropy", 0.0, true);
    ::INPUT::AddNamedString(m, "TYPE", "Problem type: Darcy (default) or Darcy-Brinkman", "Darcy");
    // optional parameter
    ::INPUT::AddNamedString(m, "PERMEABILITYFUNCTION",
        "Permeability function: Const(Default) or Kozeny_Carman", "Const", true);
    //  ::INPUT::AddNamedReal(m,"BULKMODULUS","bulk modulus of medium");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroMultiPhase",
        "multi phase flow in deformable porous media", INPAR::MAT::m_fluidporo_multiphase));

    ::INPUT::AddNamedBool(
        m, "LOCAL", "individual materials allocated per element or only at global scope");
    ::INPUT::AddNamedReal(m, "PERMEABILITY", "permeability of medium");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");
    ::INPUT::AddNamedInt(m, "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE", "number of fluid phases");
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // multiphase flow in a poroelastic material with reactions
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroMultiPhaseReactions",
        "multi phase flow in deformable porous media and list of reactions",
        INPAR::MAT::m_fluidporo_multiphase_reactions));

    ::INPUT::AddNamedBool(
        m, "LOCAL", "individual materials allocated per element or only at global scope");
    ::INPUT::AddNamedReal(m, "PERMEABILITY", "permeability of medium");
    ::INPUT::AddNamedInt(m, "NUMMAT", "number of materials in list");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "the list material IDs", "NUMMAT");
    ::INPUT::AddNamedInt(m, "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE", "number of fluid phases");
    ::INPUT::AddNamedInt(m, "NUMREAC", "number of reactions for these elements", 0);
    ::INPUT::AddNamedIntVector(m, "REACIDS", "advanced reaction list", "NUMREAC", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one reaction for multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSingleReaction",
        "advanced reaction material", INPAR::MAT::m_fluidporo_singlereaction));

    ::INPUT::AddNamedInt(m, "NUMSCAL", "number of scalars coupled with this problem");
    ::INPUT::AddNamedInt(m, "TOTALNUMDOF", "total number of multiphase-dofs");
    ::INPUT::AddNamedInt(m, "NUMVOLFRAC", "number of volfracs");
    ::INPUT::AddNamedIntVector(m, "SCALE", "advanced reaction list", "TOTALNUMDOF");
    ::INPUT::AddNamedString(m, "COUPLING",
        "type of coupling: "
        "scalar_by_function, no_coupling (default)",
        "no_coupling", false);
    ::INPUT::AddNamedInt(m, "FUNCTID", "function ID defining the reaction");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one phase for multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhase",
        "one phase for multiphase flow in deformable porous media",
        INPAR::MAT::m_fluidporo_singlephase));

    ::INPUT::AddNamedInt(m, "DENSITYLAWID", "ID of density law");
    ::INPUT::AddNamedReal(m, "DENSITY", "reference/initial density");
    ::INPUT::AddNamedInt(m, "RELPERMEABILITYLAWID", "ID of relative permeability law");
    ::INPUT::AddNamedInt(m, "VISCOSITYLAWID", "ID of viscosity law");
    ::INPUT::AddNamedInt(m, "DOFTYPEID", "ID of dof definition");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one volume fraction for multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSingleVolFrac",
        "one phase for multiphase flow in deformable porous media",
        INPAR::MAT::m_fluidporo_singlevolfrac));

    ::INPUT::AddNamedReal(m, "DENSITY", "reference/initial density");
    ::INPUT::AddNamedReal(m, "DIFFUSIVITY", "diffusivity of phase");
    ::INPUT::AddNamedBool(
        m, "AddScalarDependentFlux", "Is there additional scalar dependent flux (yes) or (no)");
    ::INPUT::AddNamedInt(m, "NUMSCAL", "Number of scalars", 0, true);
    ::INPUT::AddNamedRealVector(m, "SCALARDIFFS",
        "Diffusivities for additional scalar-dependent flux", "NUMSCAL", 0.0, true);
    ::INPUT::AddNamedRealVector(
        m, "OMEGA_HALF", "Constant for receptor kinetic law", "NUMSCAL", 1.0e13, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one volume fraction pressure for multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroVolFracPressure",
        "one volume fraction pressure for multiphase flow in deformable porous media",
        INPAR::MAT::m_fluidporo_volfracpressure));

    ::INPUT::AddNamedReal(m, "PERMEABILITY", "permeability of phase");
    ::INPUT::AddNamedInt(m, "VISCOSITYLAWID", "ID of viscosity law");
    ::INPUT::AddNamedReal(m, "MIN_VOLFRAC",
        "Minimum volume fraction under which we assume that VolfracPressure is zero", 1.0e-3, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofDiffPressure",
        "one degrree of freedom for multiphase flow in deformable porous media",
        INPAR::MAT::m_fluidporo_phasedof_diffpressure));

    ::INPUT::AddNamedInt(m, "PHASELAWID", "ID of pressure-saturation law");
    ::INPUT::AddNamedInt(m, "NUMDOF", "number of DoFs", 0);
    ::INPUT::AddNamedIntVector(
        m, "PRESCOEFF", "pressure IDs for differential pressure", "NUMDOF", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofPressure",
        "one degrree of freedom for multiphase flow in deformable porous media",
        INPAR::MAT::m_fluidporo_phasedof_pressure));

    ::INPUT::AddNamedInt(m, "PHASELAWID", "ID of pressure-saturation law");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // one degree of freedom for on single phase of a multiphase flow in a poroelastic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_FluidPoroSinglePhaseDofSaturation",
        "one degrree of freedom for multiphase flow in deformable porous media",
        INPAR::MAT::m_fluidporo_phasedof_saturation));

    ::INPUT::AddNamedInt(m, "PHASELAWID", "ID of pressure-saturation law");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // saturated law for pressure-saturation law in porous media problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PhaseLawLinear",
        "saturated fluid phase of porous medium", INPAR::MAT::m_fluidporo_phaselaw_linear));

    ::INPUT::AddNamedReal(m, "RELTENSION", "relative interface tensions");
    ::INPUT::AddNamedReal(m, "SATURATION_0", "saturation at zero differential pressure");
    ::INPUT::AddNamedInt(m, "NUMDOF", "number of DoFs", 0);
    ::INPUT::AddNamedIntVector(m, "PRESCOEFF", "Coefficients for pressure dependence", "NUMDOF", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // tangent law for pressure-saturation law in porous media multiphase problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PhaseLawTangent",
        "tangent fluid phase of porous medium", INPAR::MAT::m_fluidporo_phaselaw_tangent));

    ::INPUT::AddNamedReal(m, "RELTENSION", "relative interface tensions");
    ::INPUT::AddNamedReal(m, "EXP", "exponent in pressure-saturation law");
    ::INPUT::AddNamedReal(m, "SATURATION_0", "saturation at zero differential pressure");
    ::INPUT::AddNamedInt(m, "NUMDOF", "number of DoFs", 0);
    ::INPUT::AddNamedIntVector(m, "PRESCOEFF", "Coefficients for pressure dependence", "NUMDOF", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // constraint law for pressure-saturation law in porous media multiphase problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PhaseLawConstraint",
        "constraint fluid phase of porous medium", INPAR::MAT::m_fluidporo_phaselaw_constraint));

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // pressure-saturation law defined by functions in porous media multiphase problems
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_PhaseLawByFunction",
        "fluid phase of porous medium defined by functions",
        INPAR::MAT::m_fluidporo_phaselaw_byfunction));

    ::INPUT::AddNamedInt(m, "FUNCTPRES", "ID of function for differential pressure", 0);
    ::INPUT::AddNamedInt(m, "FUNCTSAT", "ID of function for saturation", 0);
    ::INPUT::AddNamedInt(m, "NUMDOF", "number of DoFs", 0);
    ::INPUT::AddNamedIntVector(m, "PRESCOEFF", "Coefficients for pressure dependence", "NUMDOF", 0);
    ::INPUT::AddNamedSeparator(m, "END", "indicating end of line");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // elastic spring
  {
    auto m = Teuchos::rcp(
        new MaterialDefinition("MAT_Struct_Spring", "elastic spring", INPAR::MAT::m_spring));

    ::INPUT::AddNamedReal(m, "STIFFNESS", "spring constant");
    ::INPUT::AddNamedReal(m, "DENS", "density");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // materials for beam elements (grill 02/17):

  /* The constitutive laws used in beam formulations are consistently
   * derived from a 3D solid continuum mechanics material law, e.g. a hyperelastic
   * stored energy function. The conceptual difference is that they are
   * formulated for stress and strain resultants, i.e. cross-section quantities.
   * Hence, the constitutive parameters that naturally occur in constitutive
   * relations of beam formulations are strongly related to the cross-section
   * specification (shape and dimensions) and can be identified as 'modal'
   * constitutive parameters (axial/shear/torsion/bending rigidity). See
   * Diss Meier, chapters 2.2.4 and 2.2.5 for formulae and details.
   *
   * This justifies the implementation and use of the following beam material
   * definitions. They combine cross-section specification and material definition
   * which can be done in two distinct ways:
   *
   * 1) by providing individual parameter values for cross-section specs
   *    (area, (polar) area moment of inertia, shear-correction factor, ...) and
   *    material (Young's modulus, Poisson's ratio).
   *
   * 2) by directly providing parameter values for modal constitutive parameters
   *    (axial/shear/torsion/bending rigidity).
   *    This is especially useful if experimentally determined values are used
   *    or artificial scaling of individual modes is desired in tests/debugging.
   *
   * The same logic applies to parameters required to model mass inertia.
   *
   * Reduced formulations such as Kirchhoff and isotropic/torsion-free Kirchhoff
   * beams of course require only a subset of parameters and hence use specific
   * material parameter definitions. Nevertheless, the material relations are
   * general enough such that only one class is used for the material relations of
   *  all types of beam formulations.
   */

  /*--------------------------------------------------------------------*/
  // material parameter definition for a Simo-Reissner type beam element
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition("MAT_BeamReissnerElastHyper",
        "material parameters for a Simo-Reissner type beam element based on "
        "hyperelastic stored energy function",
        INPAR::MAT::m_beam_reissner_elast_hyper));


    ::INPUT::AddNamedReal(matdef, "YOUNG", "Young's modulus");

    /* note: we define both of the two following (redundant) parameters to be optional.
     *       upon initialization of the material, we assure that one of them is
     *       properly defined. */
    ::INPUT::AddNamedReal(matdef, "SHEARMOD", "shear modulus", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "POISSONRATIO", "Poisson's ratio", -1.0, true);

    ::INPUT::AddNamedReal(matdef, "DENS", "mass density");

    ::INPUT::AddNamedReal(matdef, "CROSSAREA", "cross-section area");
    ::INPUT::AddNamedReal(matdef, "SHEARCORR", "shear correction factor");

    ::INPUT::AddNamedReal(matdef, "MOMINPOL", "polar/axial area moment of inertia");
    ::INPUT::AddNamedReal(matdef, "MOMIN2",
        "area moment of inertia w.r.t. first principal "
        "axis of inertia (i.e. second base vector)");
    ::INPUT::AddNamedReal(matdef, "MOMIN3",
        "area moment of inertia w.r.t. second principal "
        "axis of inertia (i.e. third base vector)");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);


    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }
  /*--------------------------------------------------------------------*/
  // material parameter definition for a Simo-Reissner type elasto-plastic beam element
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition("MAT_BeamReissnerElastPlastic",
        "material parameters for a Simo-Reissner type beam element based on "
        "hyperelastic stored energy function",
        INPAR::MAT::m_beam_reissner_elast_plastic));


    ::INPUT::AddNamedReal(matdef, "YOUNG", "Young's modulus");

    // optional parameters for plasticity
    ::INPUT::AddNamedReal(matdef, "YIELDN", "initial yield stress N", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "YIELDM", "initial yield stress M", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "ISOHARDN", "isotropic hardening modulus of forces", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "ISOHARDM", "isotropic hardening modulus of moments", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "TORSIONPLAST",
        "defines whether torsional moment contributes to plasticity", 0, true);

    /* note: we define both of the two following (redundant) parameters to be optional.
     *       upon initialization of the material, we assure that one of them is
     *       properly defined. */
    ::INPUT::AddNamedReal(matdef, "SHEARMOD", "shear modulus", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "POISSONRATIO", "Poisson's ratio", -1.0, true);

    ::INPUT::AddNamedReal(matdef, "DENS", "mass density");

    ::INPUT::AddNamedReal(matdef, "CROSSAREA", "cross-section area");
    ::INPUT::AddNamedReal(matdef, "SHEARCORR", "shear correction factor");

    ::INPUT::AddNamedReal(matdef, "MOMINPOL", "polar/axial area moment of inertia");
    ::INPUT::AddNamedReal(matdef, "MOMIN2",
        "area moment of inertia w.r.t. first principal "
        "axis of inertia (i.e. second base vector)");
    ::INPUT::AddNamedReal(matdef, "MOMIN3",
        "area moment of inertia w.r.t. second principal "
        "axis of inertia (i.e. third base vector)");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);


    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }
  /*--------------------------------------------------------------------*/
  // material parameter definition for a Simo-Reissner type beam element,
  // specified via 'modal' constitutive parameters (see comment above)
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition("MAT_BeamReissnerElastHyper_ByModes",
        "material parameters for a Simo-Reissner type beam element based on "
        "hyperelastic stored energy function, specified for individual "
        "deformation modes",
        INPAR::MAT::m_beam_reissner_elast_hyper_bymodes));


    ::INPUT::AddNamedReal(matdef, "EA", "axial rigidity");
    ::INPUT::AddNamedReal(matdef, "GA2", "shear rigidity w.r.t first principal axis of inertia");
    ::INPUT::AddNamedReal(matdef, "GA3", "shear rigidity w.r.t second principal axis of inertia");

    ::INPUT::AddNamedReal(matdef, "GI_T", "torsional rigidity");
    ::INPUT::AddNamedReal(matdef, "EI2",
        "flexural/bending rigidity w.r.t. first principal "
        "axis of inertia");
    ::INPUT::AddNamedReal(matdef, "EI3",
        "flexural/bending rigidity w.r.t. second principal "
        "axis of inertia");

    ::INPUT::AddNamedReal(
        matdef, "RhoA", "translational inertia: mass density * cross-section area");

    ::INPUT::AddNamedReal(matdef, "MASSMOMINPOL",
        "polar mass moment of inertia, i.e. w.r.t. "
        "rotation around beam axis");
    ::INPUT::AddNamedReal(matdef, "MASSMOMIN2",
        "mass moment of inertia w.r.t. first principal "
        "axis of inertia");
    ::INPUT::AddNamedReal(matdef, "MASSMOMIN3",
        "mass moment of inertia w.r.t. second principal "
        "axis of inertia");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);


    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }

  /*--------------------------------------------------------------------*/
  // material parameter definition for a Kirchhoff-Love type beam element
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition("MAT_BeamKirchhoffElastHyper",
        "material parameters for a Kirchhoff-Love type beam element based on "
        "hyperelastic stored energy function",
        INPAR::MAT::m_beam_kirchhoff_elast_hyper));


    ::INPUT::AddNamedReal(matdef, "YOUNG", "Young's modulus");

    /* note: we define both of the two following (redundant) parameters to be optional.
     *       upon initialization of the material, we assure that one of them is
     *       properly defined. */
    ::INPUT::AddNamedReal(matdef, "SHEARMOD", "shear modulus", -1.0, true);
    ::INPUT::AddNamedReal(matdef, "POISSONRATIO", "Poisson's ratio", -1.0, true);

    ::INPUT::AddNamedReal(matdef, "DENS", "mass density");

    ::INPUT::AddNamedReal(matdef, "CROSSAREA", "cross-section area");

    ::INPUT::AddNamedReal(matdef, "MOMINPOL", "polar/axial area moment of inertia");
    ::INPUT::AddNamedReal(matdef, "MOMIN2",
        "area moment of inertia w.r.t. first principal "
        "axis of inertia (i.e. second base vector)");
    ::INPUT::AddNamedReal(matdef, "MOMIN3",
        "area moment of inertia w.r.t. second principal "
        "axis of inertia (i.e. third base vector)");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);


    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }

  /*--------------------------------------------------------------------*/
  // material parameter definition for a Kirchhoff-Love type beam element,
  // specified via 'modal' constitutive parameters (see comment above)
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition("MAT_BeamKirchhoffElastHyper_ByModes",
        "material parameters for a Kirchhoff-Love type beam element based on "
        "hyperelastic stored energy function, specified for individual "
        "deformation modes",
        INPAR::MAT::m_beam_kirchhoff_elast_hyper_bymodes));


    ::INPUT::AddNamedReal(matdef, "EA", "axial rigidity");

    ::INPUT::AddNamedReal(matdef, "GI_T", "torsional rigidity");
    ::INPUT::AddNamedReal(matdef, "EI2",
        "flexural/bending rigidity w.r.t. first principal "
        "axis of inertia");
    ::INPUT::AddNamedReal(matdef, "EI3",
        "flexural/bending rigidity w.r.t. second principal "
        "axis of inertia");

    ::INPUT::AddNamedReal(
        matdef, "RhoA", "translational inertia: mass density * cross-section area");

    ::INPUT::AddNamedReal(matdef, "MASSMOMINPOL",
        "polar mass moment of inertia, i.e. w.r.t. "
        "rotation around beam axis");
    ::INPUT::AddNamedReal(matdef, "MASSMOMIN2",
        "mass moment of inertia w.r.t. first principal "
        "axis of inertia");
    ::INPUT::AddNamedReal(matdef, "MASSMOMIN3",
        "mass moment of inertia w.r.t. second principal "
        "axis of inertia");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);


    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }

  /*--------------------------------------------------------------------*/
  // material parameter definition for a torsion-free, isotropic
  // Kirchhoff-Love type beam element
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition("MAT_BeamKirchhoffTorsionFreeElastHyper",
        "material parameters for a torsion-free, isotropic Kirchhoff-Love "
        "type beam element based on hyperelastic stored energy function",
        INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper));


    ::INPUT::AddNamedReal(matdef, "YOUNG", "Young's modulus");

    ::INPUT::AddNamedReal(matdef, "DENS", "mass density");

    ::INPUT::AddNamedReal(matdef, "CROSSAREA", "cross-section area");

    ::INPUT::AddNamedReal(matdef, "MOMIN", "area moment of inertia");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);


    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }

  /*--------------------------------------------------------------------*/
  // material parameter definition for a torsion-free, isotropic
  // Kirchhoff-Love type beam element,
  // specified via 'modal' constitutive parameters (see comment above)
  {
    auto matdef =
        Teuchos::rcp(new MaterialDefinition("MAT_BeamKirchhoffTorsionFreeElastHyper_ByModes",
            "material parameters for a torsion-free, isotropic Kirchhoff-Love "
            "type beam element based on hyperelastic stored energy function, "
            "specified for individual deformation modes",
            INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes));


    ::INPUT::AddNamedReal(matdef, "EA", "axial rigidity");

    ::INPUT::AddNamedReal(matdef, "EI", "flexural/bending rigidity");


    ::INPUT::AddNamedReal(
        matdef, "RhoA", "translational inertia: mass density * cross-section area");
    ::INPUT::AddNamedBool(
        matdef, "FAD", "Does automatic differentiation have to be used", false, true);

    /* The following is optional because it is only required if we evaluate interactions
     * between beams such as contact, potential-based and whatever more to come.
     * For now, we always assume a circular cross-section if interactions are considered.
     *
     * This should be generalized to a type of cross-section shape (circular, rectangular,
     * elliptic, ...) and corresponding necessary dimensions (radius, sizes, ...) if needed. */
    ::INPUT::AddNamedReal(matdef, "INTERACTIONRADIUS",
        "radius of a circular cross-section which "
        "is EXCLUSIVELY used to evaluate interactions such as contact, potentials, ...",
        -1.0, true);

    AppendMaterialDefinition(matlist, matdef);
  }

  /*--------------------------------------------------------------------*/
  // material for a crosslinker in a biopolymer simulation
  {
    auto matdef = Teuchos::rcp(new MaterialDefinition(
        "MAT_Crosslinker", "material for a linkage between beams", INPAR::MAT::m_crosslinkermat));

    ::INPUT::AddNamedReal(matdef, "MATNUM", "number of beam elasthyper material");
    ::INPUT::AddNamedString(matdef, "JOINTTYPE",
        "type of joint: "
        "beam3rline2rigid (default), beam3rline2pin or truss",
        "beam3rline2rigid");
    ::INPUT::AddNamedReal(
        matdef, "LINKINGLENGTH", "distance between the two binding domains of a linker");
    ::INPUT::AddNamedReal(matdef, "LINKINGLENGTHTOL",
        "tolerance for linker length in the sense: length +- tolerance");
    ::INPUT::AddNamedReal(matdef, "LINKINGANGLE",
        "preferred binding angle enclosed by two filaments' axes in radians");
    ::INPUT::AddNamedReal(matdef, "LINKINGANGLETOL",
        "tolerance for preferred binding angle in radians in the sense of: angle +- tolerance");
    ::INPUT::AddNamedReal(matdef, "K_ON", "chemical association-rate");
    ::INPUT::AddNamedReal(matdef, "K_OFF", "chemical dissociation-rate");

    // optional parameter
    ::INPUT::AddNamedReal(
        matdef, "DELTABELLEQ", "deltaD in Bell's equation for force dependent off rate", 0.0, true);
    ::INPUT::AddNamedReal(matdef, "NOBONDDISTSPHERE",
        "distance to sphere elements in which no double bonded linker is allowed", 0.0, true);
    ::INPUT::AddNamedString(matdef, "TYPE",
        "type of crosslinker: "
        "arbitrary (default), actin, collagen, integrin",
        "arbitrary", true);

    AppendMaterialDefinition(matlist, matdef);
  }

  /*--------------------------------------------------------------------*/
  // 0D Acinar material base
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_0D_MAXWELL_ACINUS", "0D acinar material", INPAR::MAT::m_0d_maxwell_acinus));

    ::INPUT::AddNamedReal(m, "Stiffness1", "first stiffness");
    ::INPUT::AddNamedReal(m, "Stiffness2", "second stiffness");
    ::INPUT::AddNamedReal(m, "Viscosity1", "first viscosity");
    ::INPUT::AddNamedReal(m, "Viscosity2", "second viscosity");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // 0D NeoHookean Acinar material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_NEOHOOKEAN",
        "0D acinar material neohookean", INPAR::MAT::m_0d_maxwell_acinus_neohookean));

    ::INPUT::AddNamedReal(m, "Stiffness1", "first stiffness");
    ::INPUT::AddNamedReal(m, "Stiffness2", "second stiffness");
    ::INPUT::AddNamedReal(m, "Viscosity1", "first viscosity");
    ::INPUT::AddNamedReal(m, "Viscosity2", "second viscosity");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // 0D Exponential Acinar material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_EXPONENTIAL",
        "0D acinar material exponential", INPAR::MAT::m_0d_maxwell_acinus_exponential));

    ::INPUT::AddNamedReal(m, "Stiffness1", "first stiffness");
    ::INPUT::AddNamedReal(m, "Stiffness2", "second stiffness");
    ::INPUT::AddNamedReal(m, "Viscosity1", "first viscosity");
    ::INPUT::AddNamedReal(m, "Viscosity2", "second viscosity");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // 0D Exponential Acinar material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_DOUBLEEXPONENTIAL",
        "0D acinar material doubleexponential", INPAR::MAT::m_0d_maxwell_acinus_doubleexponential));

    ::INPUT::AddNamedReal(m, "Stiffness1", "first stiffness");
    ::INPUT::AddNamedReal(m, "Stiffness2", "second stiffness");
    ::INPUT::AddNamedReal(m, "Viscosity1", "first viscosity");
    ::INPUT::AddNamedReal(m, "Viscosity2", "second viscosity");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // 0D Ogden Acinar material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_0D_MAXWELL_ACINUS_OGDEN",
        "0D acinar material ogden", INPAR::MAT::m_0d_maxwell_acinus_ogden));

    ::INPUT::AddNamedReal(m, "Stiffness1", "first stiffness");
    ::INPUT::AddNamedReal(m, "Stiffness2", "second stiffness");
    ::INPUT::AddNamedReal(m, "Viscosity1", "first viscosity");
    ::INPUT::AddNamedReal(m, "Viscosity2", "second viscosity");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // O2 hemoglobin saturation material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_0D_O2_HEMOGLOBIN_SATURATION",
        "0D O2 hemoglobin saturation material", INPAR::MAT::m_0d_o2_hemoglobin_saturation));

    ::INPUT::AddNamedReal(
        m, "PerVolumeBlood", "how much of blood satisfies this rule (usually 100ml)");
    ::INPUT::AddNamedReal(m, "O2SaturationPerVolBlood",
        "O2 saturation per volume blood (In healthy blood 21.36ml/100ml of blood)");
    ::INPUT::AddNamedReal(m, "PressureHalf", "PO2 of 50\% saturated O2 (In healthy blood 26mmHg)");
    ::INPUT::AddNamedReal(m, "Power", "Power of the Sigmoidal saturation curve (2.5)");
    ::INPUT::AddNamedReal(m, "NumberOfO2PerVO2", "Number of O2 moles per unit volume of O2");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // O2 air saturation material
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_0D_O2_AIR_SATURATION",
        "0D O2 air saturation material", INPAR::MAT::m_0d_o2_air_saturation));

    ::INPUT::AddNamedReal(m, "AtmosphericPressure", "The atmospheric pressure");
    ::INPUT::AddNamedReal(m, "NumberOfO2PerVO2", "Number of O2 moles per unit volume of O2");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // particle material sph fluid
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ParticleSPHFluid",
        "particle material for SPH fluid", INPAR::MAT::m_particle_sph_fluid));

    ::INPUT::AddNamedReal(m, "INITRADIUS", "initial radius");
    ::INPUT::AddNamedReal(m, "INITDENSITY", "initial density");
    ::INPUT::AddNamedReal(m, "REFDENSFAC", "reference density factor in equation of state");
    ::INPUT::AddNamedReal(m, "EXPONENT", "exponent in equation of state");
    ::INPUT::AddNamedReal(
        m, "BACKGROUNDPRESSURE", "background pressure for transport velocity formulation");
    ::INPUT::AddNamedReal(m, "BULK_MODULUS", "bulk modulus");
    ::INPUT::AddNamedReal(m, "DYNAMIC_VISCOSITY", "dynamic shear viscosity");
    ::INPUT::AddNamedReal(m, "BULK_VISCOSITY", "bulk viscosity");
    ::INPUT::AddNamedReal(m, "ARTIFICIAL_VISCOSITY", "artificial viscosity");
    ::INPUT::AddNamedReal(m, "INITTEMPERATURE", "initial temperature", 0.0, true);
    ::INPUT::AddNamedReal(m, "THERMALCAPACITY", "thermal capacity", 0.0, true);
    ::INPUT::AddNamedReal(m, "THERMALCONDUCTIVITY", "thermal conductivity", 0.0, true);
    ::INPUT::AddNamedReal(m, "THERMALABSORPTIVITY", "thermal absorptivity", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // particle material sph boundary
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ParticleSPHBoundary",
        "particle material for SPH boundary", INPAR::MAT::m_particle_sph_boundary));

    ::INPUT::AddNamedReal(m, "INITRADIUS", "initial radius");
    ::INPUT::AddNamedReal(m, "INITDENSITY", "initial density");
    ::INPUT::AddNamedReal(m, "INITTEMPERATURE", "initial temperature", 0.0, true);
    ::INPUT::AddNamedReal(m, "THERMALCAPACITY", "thermal capacity", 0.0, true);
    ::INPUT::AddNamedReal(m, "THERMALCONDUCTIVITY", "thermal conductivity", 0.0, true);
    ::INPUT::AddNamedReal(m, "THERMALABSORPTIVITY", "thermal absorptivity", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // particle material dem
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_ParticleDEM", "particle material for DEM", INPAR::MAT::m_particle_dem));

    ::INPUT::AddNamedReal(m, "INITRADIUS", "initial radius of particle");
    ::INPUT::AddNamedReal(m, "INITDENSITY", "initial density of particle");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // particle wall material dem
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_ParticleWallDEM", "particle wall material for DEM", INPAR::MAT::m_particle_wall_dem));

    ::INPUT::AddNamedReal(
        m, "FRICT_COEFF_TANG", "friction coefficient for tangential contact", -1.0, true);
    ::INPUT::AddNamedReal(
        m, "FRICT_COEFF_ROLL", "friction coefficient for rolling contact", -1.0, true);
    ::INPUT::AddNamedReal(m, "ADHESION_SURFACE_ENERGY", "adhesion surface energy", -1.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // electromagnetic material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_Electromagnetic", "Electromagnetic material", INPAR::MAT::m_electromagneticmat));

    ::INPUT::AddNamedReal(m, "CONDUCTIVITY", "electrical conductivity");
    ::INPUT::AddNamedReal(m, "PERMITTIVITY", "Permittivity");
    ::INPUT::AddNamedReal(m, "PERMEABILITY", "Permeability");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // active fiber formation for the modeling of living cells
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_ACTIVEFIBER",
        "active fiber formation for the modeling of living cells", INPAR::MAT::m_activefiber));

    ::INPUT::AddNamedReal(m, "DENS", "Density");
    ::INPUT::AddNamedReal(m, "DECAY", "decay constant of activation signal");
    ::INPUT::AddNamedInt(
        m, "IDMATPASSIVE", "number of passive material in input file: MAT IDMATPASSIVE ...");
    ::INPUT::AddNamedReal(m, "KFOR", "formation rate parameter kforwards");
    ::INPUT::AddNamedReal(m, "KBACK", "dissociation parameter kbackwards");
    ::INPUT::AddNamedReal(m, "KVAR", "fiber rate sensitivity");
    ::INPUT::AddNamedReal(m, "SIGMAX", "maximum tension exerted by stress fibres");
    ::INPUT::AddNamedReal(m, "EPSNULL", "reference strain rate of cross-bridge dynamics law");


    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // General mixture models (used for prestretching and for homogenized constrained mixture models)
  {
    auto m = Teuchos::rcp(
        new MaterialDefinition("MAT_Mixture", "General mixture model", INPAR::MAT::m_mixture));

    ::INPUT::AddNamedInt(m, "NUMCONST", "number of mixture constituents");
    ::INPUT::AddNamedInt(m, "MATIDMIXTURERULE", "material id of the mixturerule");
    ::INPUT::AddNamedIntVector(
        m, "MATIDSCONST", "list material IDs of the mixture constituents", "NUMCONST");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for ElastHyper toolbox
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MIX_Constituent_ElastHyper", "ElastHyper toolbox", INPAR::MAT::mix_elasthyper));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of summands");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "list material IDs of the summands", "NUMMAT");
    ::INPUT::AddNamedInt(m, "PRESTRESS_STRATEGY",
        "Material id of the prestress strategy (optional, by default no prestretch)", 0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for ElastHyper toolbox with a damage process
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Constituent_ElastHyper_Damage",
        "ElastHyper toolbox with damage", INPAR::MAT::mix_elasthyper_damage));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of summands");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "list material IDs of the membrane summands", "NUMMAT");
    ::INPUT::AddNamedInt(m, "PRESTRESS_STRATEGY",
        "Material id of the prestress strategy (optional, by default no prestretch)", 0, true);
    ::INPUT::AddNamedInt(m, "DAMAGE_FUNCT",
        "Reference to the function that is a gain for the increase/decrease of the reference mass "
        "density.");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for ElastHyper toolbox with a damage process and a membrane constituent
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Constituent_ElastHyper_ElastinMembrane",
        "ElastHyper toolbox with damage and 2D membrane material",
        INPAR::MAT::mix_elasthyper_elastin_membrane));

    ::INPUT::AddNamedInt(m, "NUMMAT", "number of summands");
    ::INPUT::AddNamedIntVector(m, "MATIDS", "list material IDs of the membrane summands", "NUMMAT");
    ::INPUT::AddNamedInt(m, "MEMBRANENUMMAT", "number of summands");
    ::INPUT::AddNamedIntVector(
        m, "MEMBRANEMATIDS", "list material IDs of the membrane summands", "MEMBRANENUMMAT");
    ::INPUT::AddNamedInt(m, "PRESTRESS_STRATEGY",
        "Material id of the prestress strategy (optional, by default no prestretch)", 0, true);
    ::INPUT::AddNamedInt(m, "DAMAGE_FUNCT",
        "Reference to the function that is a gain for the increase/decrease of the reference mass "
        "density.");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for solid material
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MIX_Constituent_SolidMaterial", "Solid material", INPAR::MAT::mix_solid_material));

    ::INPUT::AddNamedInt(m, "MATID", "ID of the solid material");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Isotropic growth
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_GrowthStrategy_Isotropic", "isotropic growth",
        INPAR::MAT::mix_growth_strategy_isotropic));

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Anisotropic growth
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_GrowthStrategy_Anisotropic",
        "anisotropic growth", INPAR::MAT::mix_growth_strategy_anisotropic));


    ::INPUT::AddNamedInt(m, "INIT", "initialization modus for growth direction alignment", 1, true);
    ::INPUT::AddNamedInt(m, "FIBER_ID",
        "Id of the fiber to point the growth direction (1 for first fiber, default)", 1, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Extension of all constituents simultaneously -> Growth happens mainly in the direction with the
  // smallest stiffness
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_GrowthStrategy_Stiffness",
        "Extension of all constituents simultaneously", INPAR::MAT::mix_growth_strategy_stiffness));

    ::INPUT::AddNamedReal(
        m, "KAPPA", "Penalty parameter for the modified penalty term for incompressibility");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Constant predefined prestretch
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Prestress_Strategy_Constant",
        "Simple predefined prestress", INPAR::MAT::mix_prestress_strategy_constant));

    ::INPUT::AddNamedRealVector(m, "PRESTRETCH", "Definition of the prestretch as a 9x1 vector", 9);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Prestress strategy for a cylinder
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Prestress_Strategy_Cylinder",
        "Simple prestress strategy for a cylinder", INPAR::MAT::mix_prestress_strategy_cylinder));

    ::INPUT::AddNamedReal(m, "INNER_RADIUS", "Inner radius of the cylinder");
    ::INPUT::AddNamedReal(m, "WALL_THICKNESS", "Wall thickness of the cylinder");
    ::INPUT::AddNamedReal(m, "AXIAL_PRESTRETCH", "Prestretch in axial direction");
    ::INPUT::AddNamedReal(
        m, "CIRCUMFERENTIAL_PRESTRETCH", "Prestretch in circumferential direction");
    ::INPUT::AddNamedReal(m, "PRESSURE", "Pressure in the inner of the cylinder");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Iterative prestress strategy for any geometry
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Prestress_Strategy_Iterative",
        "Simple iterative prestress strategy for any geometry. Needed to be used within the "
        "mixture framework.",
        INPAR::MAT::mix_prestress_strategy_iterative));
    ::INPUT::AddNamedBool(m, "ACTIVE", "Flag whether prestretch tensor should be updated");
    ::INPUT::AddNamedBool(
        m, "ISOCHORIC", "Flag whether prestretch tensor is isochoric", false, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for a full constrained mixture fiber
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Constituent_FullConstrainedMixtureFiber",
        "A 1D constituent that grows with the full constrained mixture fiber theory",
        INPAR::MAT::mix_full_constrained_mixture_fiber));

    ::INPUT::AddNamedInt(m, "FIBER_ID", "Id of the fiber");
    ::INPUT::AddNamedInt(m, "FIBER_MATERIAL_ID", "Id of fiber material");
    ::INPUT::AddNamedBool(m, "GROWTH_ENABLED", "Switch for the growth", true, true);
    ::INPUT::AddNamedReal(m, "DECAY_TIME", "Decay time of deposited tissue");
    ::INPUT::AddNamedReal(m, "GROWTH_CONSTANT", "Growth constant of the tissue");
    ::INPUT::AddNamedReal(m, "DEPOSITION_STRETCH", "Stretch at which the fiber is deposited");
    ::INPUT::AddNamedInt(m, "INITIAL_DEPOSITION_STRETCH_TIMEFUNCT",
        "Id of the time function to scale the deposition stretch (Default: 0=None)", 0, true);
    ::INPUT::AddNamedInt(
        m, "INIT", "Initialization mode for fibers (1=element fibers, 3=nodal fibers)");
    ::INPUT::AddNamedBool(m, "ADAPTIVE_HISTORY",
        "Adaptively remove history snapshots based on a tolerance", false, true);
    ::INPUT::AddNamedReal(
        m, "ADAPTIVE_HISTORY_TOLERANCE", "Tolerance of the adaptive history", 1e-6, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for a remodel fiber
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Constituent_ExplicitRemodelFiber",
        "A 1D constituent that remodels", INPAR::MAT::mix_remodelfiber_expl));

    ::INPUT::AddNamedInt(m, "FIBER_ID", "Id of the fiber", 1, true);
    ::INPUT::AddNamedInt(m, "FIBER_MATERIAL_ID", "Id of fiber material");

    ::INPUT::AddNamedBool(m, "GROWTH_ENABLED", "Switch for the growth (default true)", true, true);
    ::INPUT::AddNamedReal(m, "DECAY_TIME", "Decay time of deposited tissue");
    ::INPUT::AddNamedReal(m, "GROWTH_CONSTANT", "Growth constant of the tissue");
    ::INPUT::AddNamedReal(m, "DEPOSITION_STRETCH", "Stretch at with the fiber is deposited");
    ::INPUT::AddNamedInt(m, "DEPOSITION_STRETCH_TIMEFUNCT",
        "Id of the time function to scale the deposition stretch (Default: 0=None)", 0, true);
    ::INPUT::AddNamedBool(
        m, "INELASTIC_GROWTH", "Mixture rule has inelastic growth (default false)", false, true);
    ::INPUT::AddNamedInt(
        m, "INIT", "Initialization mode for fibers (1=element fibers, 2=nodal fibers)");
    ::INPUT::AddNamedReal(
        m, "GAMMA", "Angle of fiber alignment in degree (default = 0.0 degrees)", 0.0, true);

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent for a remodel fiber
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_Constituent_ImplicitRemodelFiber",
        "A 1D constituent that remodels", INPAR::MAT::mix_remodelfiber_impl));

    ::INPUT::AddNamedInt(m, "FIBER_ID", "Id of the fiber");
    ::INPUT::AddNamedInt(m, "FIBER_MATERIAL_ID", "Id of fiber material");

    ::INPUT::AddNamedBool(m, "GROWTH_ENABLED", "Switch for the growth (default true)", true, true);
    ::INPUT::AddNamedReal(m, "DECAY_TIME", "Decay time of deposited tissue");
    ::INPUT::AddNamedReal(m, "GROWTH_CONSTANT", "Growth constant of the tissue");
    ::INPUT::AddNamedReal(m, "DEPOSITION_STRETCH", "Stretch at with the fiber is deposited");
    ::INPUT::AddNamedInt(m, "DEPOSITION_STRETCH_TIMEFUNCT",
        "Id of the time function to scale the deposition stretch (Default: 0=None)", 0, true);
    ::INPUT::AddNamedInt(
        m, "INIT", "Initialization mode for fibers (1=element fibers, 2=nodal fibers)");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent material for a remodel fiber with exponential strain energy function
  {
    auto m =
        Teuchos::rcp(new MaterialDefinition("MIX_Constituent_RemodelFiber_Material_Exponential",
            "An exponential strain energy function for the remodel fiber",
            INPAR::MAT::mix_remodelfiber_material_exponential));


    ::INPUT::AddNamedReal(m, "K1", "First parameter of exponential strain energy function");
    ::INPUT::AddNamedReal(m, "K2", "Second parameter of exponential strain energy function");
    ::INPUT::AddNamedBool(
        m, "COMPRESSION", "Bool, whether the fiber material also supports compressive forces.");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Mixture constituent material for a remodel fiber with exponential strain energy function and an
  // active contribution
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MIX_Constituent_RemodelFiber_Material_Exponential_Active",
        "An exponential strain energy function for the remodel fiber with an active contribution",
        INPAR::MAT::mix_remodelfiber_material_exponential_active));


    ::INPUT::AddNamedReal(m, "K1", "First parameter of exponential strain energy function");
    ::INPUT::AddNamedReal(m, "K2", "Second parameter of exponential strain energy function");
    ::INPUT::AddNamedBool(
        m, "COMPRESSION", "Bool, whether the fiber material also supports compressive forces.");
    ::INPUT::AddNamedReal(m, "SIGMA_MAX", "Maximum active Cauchy-stress");
    ::INPUT::AddNamedReal(m, "LAMBDAMAX", "Stretch at maximum active Cauchy-stress");
    ::INPUT::AddNamedReal(m, "LAMBDA0", "Stretch at zero active Cauchy-stress");
    ::INPUT::AddNamedReal(m, "LAMBDAACT", "Current stretch", 1.0, true);
    ::INPUT::AddNamedReal(m, "DENS", "Density of the whole mixture");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Base mixture rule for solid mixtures
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MIX_Rule_Simple", "Simple mixture rule", INPAR::MAT::mix_rule_simple));

    ::INPUT::AddNamedReal(m, "DENS", "");
    ::INPUT::AddNamedInt(m, "NUMCONST", "number of mixture constituents");
    ::INPUT::AddNamedRealVector(
        m, "MASSFRAC", "list mass fractions of the mixture constituents", "NUMCONST");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // Base mixture rule for solid mixtures
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MIX_GrowthRemodelMixtureRule",
        "Mixture rule for growth/remodel homogenized constrained mixture models",
        INPAR::MAT::mix_rule_growthremodel));

    ::INPUT::AddNamedInt(m, "GROWTH_STRATEGY", "Material id of the growth strategy");
    ::INPUT::AddNamedReal(m, "DENS", "");
    ::INPUT::AddNamedInt(m, "NUMCONST", "number of mixture constituents");
    ::INPUT::AddNamedRealVector(
        m, "MASSFRAC", "list mass fractions of the mixture constituents", "NUMCONST");

    AppendMaterialDefinition(matlist, m);
  }

  /*----------------------------------------------------------------------*/
  // crystal plasticity
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_crystal_plasticity", " Crystal plasticity ", INPAR::MAT::m_crystplast));
    ::INPUT::AddNamedReal(m, "TOL", "tolerance for internal Newton iteration");
    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "NUE", "Poisson's ratio");
    ::INPUT::AddNamedReal(m, "DENS", "Mass density");
    ::INPUT::AddNamedString(m, "LAT", "lattice type: FCC, BCC, HCP, D019 or L10", "FCC");
    ::INPUT::AddNamedReal(m, "CTOA", "c to a ratio of crystal unit cell");
    ::INPUT::AddNamedReal(m, "ABASE", "base length a of the crystal unit cell");
    ::INPUT::AddNamedInt(m, "NUMSLIPSYS", "number of slip systems");
    ::INPUT::AddNamedInt(m, "NUMSLIPSETS", "number of slip system sets");
    ::INPUT::AddNamedIntVector(m, "SLIPSETMEMBERS",
        "vector of NUMSLIPSYS indices ranging from 1 to NUMSLIPSETS that indicate to which set "
        "each slip system belongs",
        "NUMSLIPSYS");
    ::INPUT::AddNamedIntVector(m, "SLIPRATEEXP",
        "vector containing NUMSLIPSETS entries for the rate sensitivity exponent", "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "GAMMADOTSLIPREF",
        "vector containing NUMSLIPSETS entries for the reference slip shear rate", "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "DISDENSINIT",
        "vector containing NUMSLIPSETS entries for the initial dislocation density", "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "DISGENCOEFF",
        "vector containing NUMSLIPSETS entries for the dislocation generation coefficients",
        "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "DISDYNRECCOEFF",
        "vector containing NUMSLIPSETS entries for the coefficients for dynamic dislocation "
        "removal",
        "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "TAUY0",
        "vector containing NUMSLIPSETS entries for the lattice resistance to slip, e.g. the "
        "Peierls barrier",
        "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "MFPSLIP",
        "vector containing NUMSLIPSETS microstructural parameters that are relevant for Hall-Petch "
        "strengthening, e.g., grain size",
        "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "SLIPHPCOEFF",
        "vector containing NUMSLIPSETS entries for the Hall-Petch coefficients corresponding to "
        "the "
        "microstructural parameters given in MFPSLIP",
        "NUMSLIPSETS");
    ::INPUT::AddNamedRealVector(m, "SLIPBYTWIN",
        "(optional) vector containing NUMSLIPSETS entries for the work hardening coefficients by "
        "twinning on non-coplanar systems",
        "NUMSLIPSETS", 0., true);
    ::INPUT::AddNamedInt(m, "NUMTWINSYS", "(optional) number of twinning systems", 0, true);
    ::INPUT::AddNamedInt(
        m, "NUMTWINSETS", "(optional) number of sets of twinning systems", 0, true);
    ::INPUT::AddNamedIntVector(m, "TWINSETMEMBERS",
        "(optional) vector of NUMTWINSYS indices ranging from 1 to NUMTWINSETS that indicate to "
        "which set each slip system belongs",
        "NUMTWINSYS", 0, true);
    ::INPUT::AddNamedIntVector(m, "TWINRATEEXP",
        "(optional) vector containing NUMTWINSETS entries for the rate sensitivity exponent",
        "NUMTWINSETS", 0, true);
    ::INPUT::AddNamedRealVector(m, "GAMMADOTTWINREF",
        "(optional) vector containing NUMTWINSETS entries for the reference slip shear rate",
        "NUMTWINSETS", 0., true);
    ::INPUT::AddNamedRealVector(m, "TAUT0",
        "(optional) vector containing NUMTWINSETS entries for the lattice resistance to twinning, "
        "e.g. the Peierls "
        "barrier",
        "NUMTWINSETS", 0., true);
    ::INPUT::AddNamedRealVector(m, "MFPTWIN",
        "(optional) vector containing NUMTWINSETS microstructural parameters that are relevant for "
        "Hall-Petch "
        "strengthening of twins, e.g., grain size",
        "NUMTWINSETS", 0., true);
    ::INPUT::AddNamedRealVector(m, "TWINHPCOEFF",
        "(optional) vector containing NUMTWINSETS entries for the Hall-Petch coefficients "
        "corresponding to the "
        "microstructural parameters given in MFPTWIN",
        "NUMTWINSETS", 0., true);
    ::INPUT::AddNamedRealVector(m, "TWINBYSLIP",
        "(optional) vector containing NUMTWINSETS entries for the work hardening coefficients by "
        "slip",
        "NUMTWINSETS", 0., true);
    ::INPUT::AddNamedRealVector(m, "TWINBYTWIN",
        "(optional) vector containing NUMTWINSETS entries for the work hardening coefficients by "
        "twins on non-coplanar systems",
        "NUMTWINSETS", 0., true);
    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // linear elastic material in one direction
  {
    auto m = Teuchos::rcp(new MaterialDefinition(
        "MAT_LinElast1D", "linear elastic material in one direction", INPAR::MAT::m_linelast1D));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");

    AppendMaterialDefinition(matlist, m);
  }

  /*--------------------------------------------------------------------*/
  // linear elastic material with growth in one direction
  {
    auto m = Teuchos::rcp(new MaterialDefinition("MAT_LinElast1DGrowth",
        "linear elastic material with growth in one direction", INPAR::MAT::m_linelast1D_growth));

    ::INPUT::AddNamedReal(m, "YOUNG", "Young's modulus");
    ::INPUT::AddNamedReal(m, "DENS", "mass density");
    ::INPUT::AddNamedReal(m, "C0", "reference concentration");
    ::INPUT::AddNamedBool(m, "AOS_PROP_GROWTH",
        "growth proportional to amount of substance (AOS) if true or proportional to concentration "
        "if false");
    ::INPUT::AddNamedInt(m, "POLY_PARA_NUM", "number of polynomial coefficients");
    ::INPUT::AddNamedRealVector(m, "POLY_PARAMS", "coefficients of polynomial", "POLY_PARA_NUM");

    AppendMaterialDefinition(matlist, m);
  }

  // deliver
  return vm;
}
