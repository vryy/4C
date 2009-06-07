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
#ifdef CCADISCRET

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
extern "C"
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

    AddNamedReal(m,"VISCOSITY","kinematic or dynamic viscosity");
    AddNamedReal(m,"DENSITY","spatial mass density");
    AddNamedReal(m,"GAMMA","surface tension coeficient",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // fluid material according to Sutherland law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_sutherland_fluid",
                                            "fluid material according to Sutherland law ",
                                            INPAR::MAT::m_sutherland_fluid));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");

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
  // convection-diffusion material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_condif",
                                            "convection-diffusion material ",
                                            INPAR::MAT::m_condif));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"REACOEFF","reaction coefficient",true);
    AddNamedReal(m,"SHC","specific heat capacity",true);

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // convection-diffusion material according to Sutherland law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_sutherland_condif",
                                            "convection-diffusion material according to Sutherland law",
                                            INPAR::MAT::m_sutherland_condif));

    AddNamedReal(m,"REFVISC","reference dynamic viscosity (kg/(m*s))");
    AddNamedReal(m,"REFTEMP","reference temperature (K)");
    AddNamedReal(m,"SUTHTEMP","Sutherland temperature (K)");
    AddNamedReal(m,"SHC","specific heat capacity");
    AddNamedReal(m,"PRANUM","Prandtl number");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material parameters for ion species in electrlyte solution (gjb 07/08)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ion",
                                            "material parameters for ion species in electrolyte solution",
                                            INPAR::MAT::m_ion));

    AddNamedReal(m,"DIFFUSIVITY","kinematic diffusivity");
    AddNamedReal(m,"VALENCE","valence (= charge number)");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // material collection (gjb 07/08)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_matlist",
                                            "list/collection of materials, i.e. material IDs",
                                            INPAR::MAT::m_matlist));

    AddNamedInt(m,"NUMMAT","number of materials in list");
    AddNamedIntVector(m,"MATIDS","the list material IDs","NUMMAT");

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

  /*----------------------------------------------------------------------*/
  // Elastic orthotropic material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Orthotropic",
                                            "Elastic orthotropic material",
                                            INPAR::MAT::m_el_orth));

    AddNamedReal(m,"EMOD1","???");
    AddNamedReal(m,"EMOD2","???");
    AddNamedReal(m,"EMOD3","???");
    AddNamedReal(m,"GMOD12","???");
    AddNamedReal(m,"GMOD13","???");
    AddNamedReal(m,"GMOD23","???");
    AddNamedReal(m,"XNUE12","???");
    AddNamedReal(m,"XNUE13","???");
    AddNamedReal(m,"XNUE23","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Porous St.Venant--Kirchhoff material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_STVENPOR",
                                            "Porous St.Venant--Kirchhoff material",
                                            INPAR::MAT::m_stvenpor));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Possion's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedReal(m,"REFDENS","reference density");
    AddNamedReal(m,"EXPO","material parameter");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // neo-Hooke
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_NeoHooke",
                                            "neo-Hooke material",
                                            INPAR::MAT::m_neohooke));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");

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
  // logarithmic neo-Hooke material acc. to Bonet and Wood
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_LogNeoHooke",
                                            "logarithmic neo-Hooke material acc. to Bonet and Wood",
                                            INPAR::MAT::m_logneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");
    AddNamedReal(m,"DENS","mass density");
    AddNamedInt(m,"MODEL","sub model: 0=Bonet&Wood, 1=Volumetrically-isochorically decomposed",0,true);

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // Biological cell material model
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_BioCell",
                                            "Biological cell material model",
                                            INPAR::MAT::m_biocell));

    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // CHARMM API
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_CHARMM",
                                            "CHARmm API",
                                            INPAR::MAT::m_charmm));
    AddNamedInt(m,"ORIGIN","Evaluation at origin");
    AddNamedReal(m,"FCL","First characteristic length");
    AddNamedString(m,"FCD_TYPE","Type of the first characteristic direction","none");
    AddNamedRealVector(m,"FCD","First characteristic direction",3);
    AddNamedRealVector(m,"FCD_Space","First characteristic directional space",3);
    AddNamedReal(m,"SCL","Second characteristic length");
    AddNamedString(m,"SCD_TYPE","Type of the second charateristic direction","none");
    AddNamedRealVector(m,"SCD","Second characteristic direction",3);
    AddNamedRealVector(m,"SCD_Space","Second characteristic directional space",3);
    AddNamedInt(m,"FCD_Acceleration","Acceleration computation in FCD");
    AddNamedReal(m,"AtomicMass","Atomic mass [amu] of the moving part");
    AddNamedReal(m,"Facc_Scale","Scale factor from FE force to pN");
    AddNamedReal(m,"Time_AKMA","Scale factor from FE time to AKMA time");
    AddNamedReal(m,"Time_Scale","Linear scale factor for time span");
    AddNamedInt(m,"HARD","Use hard coded results");
    AddNamedReal(m,"c_Scale","Scale factor for c (Neo-Hookean)");
    AddNamedString(m,"PATH","Location of CHARMm problem case","none");
    AddNamedInt(m,"USE_OLD_RESULTS","Reuse previously computed results from CHARMm");
    AddNamedString(m,"SERPAR","Serial or parallel computations","ser");
    AddNamedString(m,"CHARMM","CHARMm binary location","none");
    AddNamedString(m,"INPUT","CHARMm input file","none");
    AddNamedReal(m,"NUE","Poisson ratio");
    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }
  /*--------------------------------------------------------------------*/
  // CHARMM API for proteins
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_PROTEIN",
                                            "CHARmm API for Proteins",
                                            INPAR::MAT::m_protein));

    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }
  /*--------------------------------------------------------------------*/
  // Itskov material law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ITSKOV",
                                            "Itskov material law",
                                            INPAR::MAT::m_itskov));

    AddNamedReal(m,"ALPHA","material parameter fibers");
    AddNamedReal(m,"BETA","material parameter fibers");
    AddNamedReal(m,"MU_FIBERS","mu fibers");
    AddNamedReal(m,"MU_GS","mu ground substance");
    AddNamedReal(m,"EPSILON","penalty parameter");
    AddNamedReal(m,"GAMMA","penalty parameter");
    AddNamedReal(m,"C","variable incompressibility");
    AddNamedReal(m,"DENS","mass density");

    AppendMaterialDefinition(matlist,m);
  }
  /*--------------------------------------------------------------------*/
  // MFOC
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_MFOC",
                                            "open cell foam material",
                                            INPAR::MAT::m_mfoc));

    AddNamedReal(m,"Es","Young's modulus (cell)");
    AddNamedReal(m,"pr","Possion's ratio");
    AddNamedReal(m,"dens","density foam");
    AddNamedReal(m,"denss","density (bulk)");
    AddNamedReal(m,"oce","exponent");
    AddNamedReal(m,"ocf","factor");
    AddNamedReal(m,"densmin","min. dens. foam (opti.)");
    AddNamedReal(m,"densmax","max. dens. foam (opti.)");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // MFCC
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_MFCC",
                                            "closed cell foam material",
                                            INPAR::MAT::m_mfcc));

    AddNamedReal(m,"Es","Young's modulus (cell)");
    AddNamedReal(m,"pr","Possion ratio");
    AddNamedReal(m,"dens","density foam");
    AddNamedReal(m,"denss","density (bulk)");
    AddNamedReal(m,"cce","exponent");
    AddNamedReal(m,"ccf","factor");
    AddNamedReal(m,"densmin","min. dens. foam (opti.)");
    AddNamedReal(m,"densmax","max. dens. foam (opti.)");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // NeoHMFCC
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_NeoHMFCC",
                                            "foam material",
                                            INPAR::MAT::m_nhmfcc));

    AddNamedReal(m,"Es","???");
    AddNamedReal(m,"pr","???");
    AddNamedReal(m,"dens","por. density");
    AddNamedReal(m,"denss","ref. density");
    AddNamedReal(m,"cce","???");
    AddNamedReal(m,"ccf","???");
    AddNamedReal(m,"densmin","???");
    AddNamedReal(m,"densmax","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Ogden
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Ogden",
                                            "???",
                                            INPAR::MAT::m_compogden));

    AddNamedReal(m,"NUE","???");
    AddNamedReal(m,"BETA","???");
    AddNamedReal(m,"ALFA1","???");
    AddNamedReal(m,"ALFA2","???");
    AddNamedReal(m,"ALFA3","???");
    AddNamedReal(m,"NU1","???");
    AddNamedReal(m,"NU2","???");
    AddNamedReal(m,"NU3","???");
    AddNamedReal(m,"DENS","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // visco-hyperelastic
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_Viscohyper",
                                            "???",
                                            INPAR::MAT::m_viscohyper));

    AddNamedReal(m,"NUE","???");
    AddNamedReal(m,"BETA","???");
    AddNamedReal(m,"ALFA1","???");
    AddNamedReal(m,"ALFA2","???");
    AddNamedReal(m,"ALFA3","???");
    AddNamedReal(m,"NU1","???");
    AddNamedReal(m,"NU2","???");
    AddNamedReal(m,"NU3","???");
    AddNamedReal(m,"DENS","???");
    AddNamedInt(m,"NMAXW","???");
    AddNamedReal(m,"TAU1","???");
    AddNamedReal(m,"TAU2","???");
    AddNamedReal(m,"TAU3","???");
    AddNamedReal(m,"TAU4","???");
    AddNamedReal(m,"BETA1","???");
    AddNamedReal(m,"BETA2","???");
    AddNamedReal(m,"BETA3","???");
    AddNamedReal(m,"BETA4","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // 3D von Mises plasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_3DMisesPlastic",
                                            "???",
                                            INPAR::MAT::m_pl_mises_3D));

    AddNamedReal(m,"YOUNG","???");
    AddNamedReal(m,"NUE","???");
    AddNamedReal(m,"ALFAT","???");
    AddNamedReal(m,"Sigy","???");
    AddNamedReal(m,"Hard","???");
    AddNamedReal(m,"GF","???");
    AddNamedReal(m,"BETAH","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // von Mises plasticity
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_MisesPlastic",
                                            "???",
                                            INPAR::MAT::m_pl_mises));

    AddNamedReal(m,"YOUNG","???");
    AddNamedReal(m,"NUE","???");
    AddNamedReal(m,"ALFAT","???");
    AddNamedReal(m,"Sigy","???");
    AddNamedReal(m,"Hard","???");
    AddNamedReal(m,"GF","???");
    AddNamedReal(m,"BETAH","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // Damage material
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Damage",
                                            "???",
                                            INPAR::MAT::m_damage));

    AddNamedReal(m,"YOUNG","???");
    AddNamedReal(m,"NUE","???");
    AddNamedInt(m,"Equival","???");
    AddNamedInt(m,"Damtyp","???");
    AddNamedReal(m,"Kappa_0","???");
    AddNamedReal(m,"Kappa_m","???");
    AddNamedReal(m,"Alpha","???");
    AddNamedReal(m,"Beta","???");
    AddNamedReal(m,"k_fac","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // Plastic foam
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_FoamPlastic",
                                            "???",
                                            INPAR::MAT::m_pl_foam));

    AddNamedReal(m,"YOUNG","???");
    AddNamedReal(m,"NUE","???");
    AddNamedReal(m,"ALFAT","???");
    AddNamedReal(m,"Sigy","???");
    AddNamedReal(m,"Hard","???");
    AddNamedReal(m,"GF","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // DP_Plastic
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_DP_Plastic",
                                            "???",
                                            INPAR::MAT::m_pl_dp));

    AddNamedReal(m,"YOUNG","???");
    AddNamedReal(m,"NUE","???");
    AddNamedReal(m,"ALFAT","???");
    AddNamedReal(m,"Sigy","???");
    AddNamedReal(m,"PHI","???");
    AddNamedReal(m,"Hard","???");
    AddNamedReal(m,"GF","???");
    AddNamedReal(m,"BETAH","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Lung Ogden
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_LungOgden",
                                            "lung Ogden",
                                            INPAR::MAT::m_lung_ogden));

    AddNamedReal(m,"C","???");
    AddNamedReal(m,"K1","???");
    AddNamedReal(m,"K2","???");
    AddNamedReal(m,"KAPPA","???");
    AddNamedReal(m,"BETA","???");
    AddNamedReal(m,"DENS","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Lung penalty
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_Struct_LungPenalty",
                                            "lung penalty",
                                            INPAR::MAT::m_lung_penalty));

    AddNamedReal(m,"C","???");
    AddNamedReal(m,"K1","???");
    AddNamedReal(m,"K2","???");
    AddNamedReal(m,"EPSILON","???");
    AddNamedReal(m,"GAMMA","???");
    AddNamedReal(m,"DENS","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Anisotropic Polyconvex Material Law based on Balzani et. al.
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ANISOTROPIC_BALZANI",
                                            "Anisotropic Polyconvex Material Law based on Balzani et. al.",
                                            INPAR::MAT::m_anisotropic_balzani));

    AddNamedReal(m,"C1","???");
    AddNamedReal(m,"EPS1","???");
    AddNamedReal(m,"EPS2","???");
    AddNamedReal(m,"ALPHA1","???");
    AddNamedReal(m,"ALPHA2","???");
    AddNamedReal(m,"DENS","???");
    AddNamedInt(m,"ALOC","???");
    AddNamedReal(m,"A1X","???");
    AddNamedReal(m,"A1Y","???");
    AddNamedReal(m,"A1Z","???");
    AddNamedReal(m,"ALPHA1_2","???");
    AddNamedReal(m,"ALPHA2_2","???");
    AddNamedReal(m,"A2X","???");
    AddNamedReal(m,"A2Y","???");
    AddNamedReal(m,"A2Z","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Mooney-Rivlin material law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_MOONEYRIVLIN",
                                            "Mooney-Rivlin material law",
                                            INPAR::MAT::m_mooneyrivlin));

    AddNamedReal(m,"C1","???");
    AddNamedReal(m,"C2","???");
    AddNamedReal(m,"KAPPA","???");
    AddNamedReal(m,"LAMBDA","???");
    AddNamedReal(m,"DENS","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Yeoh
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_YEOH",
                                            "hyperelastic material based on Yeoh",
                                            INPAR::MAT::m_yeoh));

    AddNamedReal(m,"C1","linear shear stiffness");
    AddNamedReal(m,"C2","quadratic shear stiffness");
    AddNamedReal(m,"C3","cubic shear stiffness");
    AddNamedReal(m,"KAPPA","volume dilatation modulus");
    AddNamedReal(m,"DENS","density");

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
  // Continuum chain network material law
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_CONTCHAINNETW",
                                            "continuum chain network material law",
                                            INPAR::MAT::m_contchainnetw));

    AddNamedReal(m,"LAMBDA","???");
    AddNamedReal(m,"MUE","???");
    AddNamedReal(m,"DENS","???");
    AddNamedReal(m,"NCHAIN","???");
    AddNamedReal(m,"ABSTEMP","???");
    AddNamedReal(m,"CONTL_L","???");
    AddNamedReal(m,"PERSL_A","???");
    AddNamedReal(m,"R0","???");
    AddNamedReal(m,"RELAX","???");
    AddNamedReal(m,"REMBEGT","???");
    AddNamedInt(m,"INITRAN","???");
    AddNamedInt(m,"UPDRATE","???");
    AddNamedReal(m,"DIFFTOL","???");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // Arterial wall material law (Holzapfel) with remodeling (Hariton)
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("MAT_ARTWALLREMOD",
                                            "Arterial wall material law (Holzapfel) with remodeling (Hariton)",
                                            INPAR::MAT::m_artwallremod));

    AddNamedReal(m,"MUE","???");
    AddNamedReal(m,"K1","???");
    AddNamedReal(m,"K2","???");
    AddNamedReal(m,"KAPPA","???");
    AddNamedReal(m,"DENS","???");
    AddNamedReal(m,"REMBEGT","???");
    AddNamedInt(m,"INIT","???");
    AddNamedReal(m,"GAMMA","???");
    AddNamedInt(m,"TENSION_ONLY","???");

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
    AddNamedReal(m,"GAMMA","fiber angle");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // logarithmic neo-Hooke material acc. to Bonet and Wood
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupLogNeoHooke",
                                            "logarithmic neo-Hooke material acc. to Bonet and Wood",
                                            INPAR::MAT::mes_couplogneohooke));

    AddNamedReal(m,"YOUNG","Young's modulus");
    AddNamedReal(m,"NUE","Poisson's ratio");

    AppendMaterialDefinition(matlist,m);
  }

  /*--------------------------------------------------------------------*/
  // coupled Blatz and Ko material acc. to Holzapfel
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupBlatzKo",
                                            "Blatz and Ko material acc. to Holtzapfel",
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
  // isochoric contribution of neohooke
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
  // coupled anisotropic material with two exponential fiber families
  {
    Teuchos::RCP<MaterialDefinition> m
      = Teuchos::rcp(new MaterialDefinition("ELAST_CoupAnisoExpoTwo",
                                            "anisotropic part with two exp. fibers",
                                            INPAR::MAT::mes_coupanisoexpotwo));

    AddNamedReal(m,"K1","linear constant fiber 1");
    AddNamedReal(m,"K2","exponential constant fiber 1");
    AddNamedReal(m,"K3","linear constant fiber 2");
    AddNamedReal(m,"K4","exponential constant fiber 2");

    AppendMaterialDefinition(matlist,m);
  }

  /*----------------------------------------------------------------------*/
  // deliver
  return vm;
}

#endif
