/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.cpp

\brief Setup of the list of valid input parameters

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <Teuchos_Array.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_any.hpp>

#include <AztecOO.h>


#include "drt_validparameters.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_inpar/inpar_solver.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_structure.H"


/*----------------------------------------------------------------------*/
//! Print function to be called from C
/*----------------------------------------------------------------------*/
extern "C"
void PrintValidParameters()
{
  Teuchos::RCP<const Teuchos::ParameterList> list = DRT::INPUT::ValidParameters();
  list->print(std::cout,
              Teuchos::ParameterList::PrintOptions()
              .showDoc(true)
              .showFlags(false)
              .indent(4)
              .showTypes(false));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintDatHeader(std::ostream& stream, const Teuchos::ParameterList& list, std::string parentname, bool color)
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

  // prevent invalid ordering of parameters caused by alphabetical output:
  // in the first run, print out all list elements that are not a sublist
  // in the second run, do the recursive call for all the sublists in the list
  for (int j=0; j<2; ++j)
  {
    for (Teuchos::ParameterList::ConstIterator i = list.begin();
    i!=list.end();
    ++i)
    {
      const Teuchos::ParameterEntry& entry = list.entry(i);
      if (entry.isList() && j==0) continue;
      if ((!entry.isList()) && j==1) continue;
      const std::string &name = list.name(i);
      Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();

      stream << blue2light << "//" << endcolor << '\n';

      std::string doc = entry.docString();
      if (doc!="")
      {
        Teuchos::StrUtils::printLines(stream,blue2light + "// ",doc);
        stream << endcolor;
      }

      if (entry.isList())
      {
        std::string secname = parentname;
        if (secname!="")
          secname += "/";
        secname += name;
        unsigned l = secname.length();
        stream << redlight << "--";
        for (int i=0; i<std::max<int>(65-l,0); ++i) stream << '-';
        stream << greenlight << secname << endcolor << '\n';
        PrintDatHeader(stream,list.sublist(name),secname,color);
      }
      else
      {
        if (validator!=Teuchos::null)
        {
          Teuchos::RCP<const Teuchos::Array<std::string> > values = validator->validStringValues();
          if (values!=Teuchos::null)
          {
            unsigned len = 0;
            for (unsigned i=0; i<values->size(); ++i)
            {
              len += (*values)[i].length()+1;
            }
            if (len<74)
            {
              stream << blue2light << "//     ";
              for (int i=0; i<static_cast<int>(values->size())-1; ++i)
              {
                stream << magentalight << (*values)[i] << blue2light << ",";
              }
              stream << magentalight << (*values)[values->size()-1] << endcolor << '\n';
            }
            else
            {
              for (unsigned i=0; i<values->size(); ++i)
              {
                stream << blue2light << "//     " << magentalight << (*values)[i] << endcolor << '\n';
              }
            }
          }
        }
        const Teuchos::any& v = entry.getAny(false);
        stream << bluelight << name << endcolor;
        unsigned l = name.length();
        for (int i=0; i<std::max<int>(31-l,0); ++i) stream << ' ';
        stream << ' ' << yellowlight << v << endcolor << '\n';
      }
    }
  }
}


/*----------------------------------------------------------------------*/
//! Print function to be called from C
/*----------------------------------------------------------------------*/
extern "C"
void PrintDefaultDatHeader()
{
  Teuchos::RCP<const Teuchos::ParameterList> list = DRT::INPUT::ValidParameters();
  DRT::INPUT::PrintDatHeader(std::cout,*list,"",true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::PrintDefaultParameters(std::ostream& stream, const Teuchos::ParameterList& list)
{
  bool hasDefault = false;
  for (Teuchos::ParameterList::ConstIterator i = list.begin();
       i!=list.end();
       ++i)
  {
    const Teuchos::ParameterEntry& entry = list.entry(i);
    if (entry.isDefault())
    {
      if (not hasDefault)
      {
        hasDefault = true;
        stream << "default parameters in list '" << list.name() << "':\n";
      }
      const Teuchos::any& v = entry.getAny(false);
      int l = list.name(i).length();
      stream << "    " << list.name(i);
      for (int i=0; i<std::max<int>(31-l,0); ++i) stream << ' ';
      stream << ' ' << v << '\n';
    }
  }
  if (hasDefault)
    stream << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::BoolParameter(std::string const& paramName,
                               std::string const& value,
                               std::string const& docString,
                               Teuchos::ParameterList* paramList)
{
  Teuchos::Array<std::string> yesnotuple = Teuchos::tuple<std::string>(
    "Yes",
    "No",
    "yes",
    "no",
    "YES",
    "NO");
  Teuchos::Array<int> yesnovalue = Teuchos::tuple<int>(
    true,
    false,
    true,
    false,
    true,
    false);
  Teuchos::setStringToIntegralParameter<int>(
    paramName,value,docString,
    yesnotuple,yesnovalue,
    paramList);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::IntParameter(std::string const &paramName,
                              int const value,
                              std::string const &docString,
                              Teuchos::ParameterList *paramList)
{
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
  validator.allowInt(true);
  Teuchos::setIntParameter(paramName,value,
                           docString,
                           paramList,validator);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::DoubleParameter(std::string const &paramName,
                                 double const &value,
                                 std::string const &docString,
                                 Teuchos::ParameterList *paramList)
{
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
  validator.allowDouble(true);
  validator.allowInt(true);
  Teuchos::setDoubleParameter(paramName,value,
                              docString,
                              paramList,validator);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> DRT::INPUT::ValidParameters()
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& discret = list->sublist("DISCRETISATION",false,"");

  IntParameter("NUMFLUIDDIS",1,"Number of meshes in fluid field",&discret);
  IntParameter("NUMSTRUCDIS",1,"Number of meshes in structural field",&discret);
  IntParameter("NUMALEDIS",1,"Number of meshes in ale field",&discret);
  IntParameter("NUMTHERMDIS",1,"Number of meshes in thermal field",&discret);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& size = list->sublist("PROBLEM SIZE",false,"");

  IntParameter("ELEMENTS",0,"Total number of elements",&size);
  IntParameter("NODES",0,"Total number of nodes",&size);
  IntParameter("NPATCHES",0,"number of nurbs patches",&size);
  IntParameter("DIM",3,"2d or 3d problem",&size);
  IntParameter("MATERIALS",0,"number of materials",&size);
  IntParameter("NUMDF",3,"maximum number of degrees of freedom",&size);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& type = list->sublist("PROBLEM TYP",false,"");

  setStringToIntegralParameter<PROBLEM_TYP>(
                               "PROBLEMTYP",
                               "Fluid_Structure_Interaction",
                               "",
                               tuple<std::string>(
                                 "Structure",
                                 "Fluid",
                                 "Fluid_XFEM",
                                 "Fluid_Ale",
                                 "Fluid_Freesurface",
                                 "Scalar_Transport",
                                 "Fluid_Structure_Interaction",
                                 "Fluid_Structure_Interaction_XFEM",
                                 "Ale",
                                 "Thermal_Structure_Interaction",
                                 "Structure_Multiscale",
                                 "Low_Mach_Number_Flow",
                                 "Electrochemistry",
                                 "Combustion"),
                               tuple<PROBLEM_TYP>(
                                 prb_structure,
                                 prb_fluid,
                                 prb_fluid_xfem,
                                 prb_fluid_ale,
                                 prb_freesurf,
                                 prb_scatra,
                                 prb_fsi,
                                 prb_fsi_xfem,
                                 prb_ale,
                                 prb_tsi,
                                 prb_struct_multi,
                                 prb_loma,
                                 prb_elch,
                                 prb_combust),
                                 &type);

  IntParameter("NUMFIELD",1,"",&type);
  setStringToIntegralParameter<TIME_TYP>("TIMETYP","Dynamic","",
                               tuple<std::string>("Static","Dynamic"),
                               tuple<TIME_TYP>(time_static,time_dynamic),
                               &type);
  //IntParameter("GRADERW",0,"",&type);
  IntParameter("MULTISC_STRUCT",0,"",&type);
  IntParameter("RESTART",0,"",&type);
  setStringToIntegralParameter<int>("ALGEBRA","Trilinos","outdated",
                               tuple<std::string>("Trilinos","ccarat"),
                               tuple<int>(1,0),
                               &type);
  setStringToIntegralParameter<int>("SHAPEFCT","Polynomial","Defines the function spaces for the spatial approximation",
                               tuple<std::string>("Polynomial","Nurbs"),
                               tuple<int>(1,0),
                               &type);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& io = list->sublist("IO",false,"");

  // are these needed?
  setStringToIntegralParameter<int>("OUTPUT_OUT","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("OUTPUT_GID","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("OUTPUT_BIN","No","",yesnotuple,yesnovalue,&io);

  setStringToIntegralParameter<int>("STRUCT_DISP","Yes","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<INPAR::STR::StressType>("STRUCT_STRESS","No","",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "Cauchy","cauchy",
                                                  "2PK", "2pk"),
                               tuple<INPAR::STR::StressType>(INPAR::STR::stress_none,INPAR::STR::stress_none,INPAR::STR::stress_none,
                                                             INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,
                                                             INPAR::STR::stress_cauchy,INPAR::STR::stress_cauchy,
                                                             INPAR::STR::stress_2pk,INPAR::STR::stress_2pk),
                               &io);
  setStringToIntegralParameter<INPAR::STR::StrainType>("STRUCT_STRAIN","No","",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "EA","ea",
                                                  "GL", "gl"),
                               tuple<INPAR::STR::StrainType>(INPAR::STR::strain_none,INPAR::STR::strain_none,INPAR::STR::strain_none,
                                                             INPAR::STR::strain_gl,INPAR::STR::strain_gl,INPAR::STR::strain_gl,
                                                             INPAR::STR::strain_ea,INPAR::STR::strain_ea,
                                                             INPAR::STR::strain_gl,INPAR::STR::strain_gl),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_SURFACTANT","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_SM_DISP","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_SOL","Yes","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_VIS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("ALE_DISP","No","",yesnotuple,yesnovalue,&io);

  setStringToIntegralParameter<int>("THERM_TEMPERATURE","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("THERM_HEATFLUX","No","",yesnotuple,yesnovalue,&io);

  IntParameter("FILESTEPS",1000,"",&io);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& design = list->sublist("DESIGN DESCRIPTION",false,"number of nodal clouds");

  IntParameter("NDPOINT",0,"number of points",&design);
  IntParameter("NDLINE",0,"number of line clouds",&design);
  IntParameter("NDSURF",0,"number of surface clouds",&design);
  IntParameter("NDVOL",0,"number of volume clouds",&design);

  /*----------------------------------------------------------------------*/
  // An empty list. The actual list is arbitrary and not validated.
  Teuchos::ParameterList& condition =
    list->sublist("CONDITION NAMES",false,
                "Names of conditions from exodus file.\n"
                "This section is not validated, any variable is allowed here.\n"
                "The names defined in this section can be used by all conditions instead of\n"
                "a design object number. This section assigns the respective numbers to\n"
                "the names.");

  condition.disableRecursiveValidation();

  //ParameterList& stat = list->sublist("STATIC",false,"");

  /*----------------------------------------------------------------------*/
  //Teuchos::ParameterList& eigen = list->sublist("EIGENVALUE ANALYSIS",false,"");


  /*--------------------------------------------------------------------*/
  /* parameters for NOX - non-linear solution */
  Teuchos::ParameterList& snox = list->sublist("STRUCT NOX",false,"");
  SetValidNoxParameters(snox);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& sdyn = list->sublist("STRUCTURAL DYNAMIC",false,"");

  setStringToIntegralParameter<INPAR::STR::DynamicType>("DYNAMICTYP","Gen_Alfa",
                               "type of time integration control",
                               tuple<std::string>(
                                 "Centr_Diff",
                                 "Gen_EMM",
                                 "Gen_Alfa",
                                 "Static",
                                 "Statics",
                                 "GenAlpha",
                                 "OneStepTheta",
                                 "GEMM",
                                 "AdamsBashforth2"),
                               tuple<INPAR::STR::DynamicType>(
                                 INPAR::STR::dyna_centr_diff,
                                 INPAR::STR::dyna_Gen_EMM,
                                 INPAR::STR::dyna_gen_alfa,
                                 INPAR::STR::dyna_gen_alfa_statics,
                                 INPAR::STR::dyna_statics,
                                 INPAR::STR::dyna_genalpha,
                                 INPAR::STR::dyna_onesteptheta,
                                 INPAR::STR::dyna_gemm,
                                 INPAR::STR::dyna_ab2),
                               &sdyn);
  // a temporary flag
  setStringToIntegralParameter<int>("ADAPTERDRIVE","No",
                                    "TEMPORARY FLAG: Switch on time integration driver based on ADAPTER::Structure rather than independent implementation",
                                    yesnotuple,yesnovalue,&sdyn);

  // Output type
  IntParameter("EIGEN",0,"EIGEN make eigenanalysis of the initial dynamic system",&sdyn);
  IntParameter("RESEVRYDISP",1,"save displacements and contact forces every RESEVRYDISP steps",&sdyn);
  IntParameter("RESEVRYSTRS",1,"save stresses every RESEVRYSTRS steps",&sdyn);
  IntParameter("RESEVRYERGY",0,"write system energies every requested step",&sdyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&sdyn);
  // Time loop control
  DoubleParameter("TIMESTEP",0.05,"time step size",&sdyn);
  IntParameter("NUMSTEP",200,"maximum number of steps",&sdyn);
  DoubleParameter("MAXTIME",5.0,"maximum time",&sdyn);
  // Generalised-alpha parameters
  DoubleParameter("BETA",0.25,"generalized alpha factors, also used by explicit time integration",&sdyn);
  DoubleParameter("DELTA",0.25,"generalized alpha factors",&sdyn);
  DoubleParameter("GAMMA",0.5,"generalized alpha factors, also used by explicit time integration",&sdyn);
  DoubleParameter("ALPHA_M",0.5,"generalized alpha factors",&sdyn);
  DoubleParameter("ALPHA_F",0.5,"generalized alpha factors",&sdyn);
  // Damping
  setStringToIntegralParameter<INPAR::STR::DampKind>("DAMPING","No",
                               "type of damping: (1) Rayleigh damping matrix and use it from M_DAMP x M + K_DAMP x K, (2) Material based and calculated in elements",
                               tuple<std::string>(
                                 "no",
                                 "No",
                                 "NO",
                                 "yes",
                                 "Yes",
                                 "YES",
                                 "Rayleigh",
                                 "Material"),
                               tuple<INPAR::STR::DampKind>(
                                 INPAR::STR::damp_none,
                                 INPAR::STR::damp_none,
                                 INPAR::STR::damp_none,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_material),
                               &sdyn);
  DoubleParameter("M_DAMP",0.5,"",&sdyn);
  DoubleParameter("K_DAMP",0.5,"",&sdyn);
  // Iteration
  setStringToIntegralParameter<int>("ITERATION","full","unused",
                               tuple<std::string>("full","Full","FULL"),
                               tuple<int>(1,1,1),
                               &sdyn);
  setStringToIntegralParameter<INPAR::STR::ConvCheck>("CONV_CHECK","AbsRes_Or_AbsDis","type of convergence check",
                               tuple<std::string>(
                                 "AbsRes_Or_AbsDis",
                                 "AbsRes_And_AbsDis",
                                 "RelRes_Or_AbsDis",
                                 "RelRes_And_AbsDis",
                                 "RelRes_Or_RelDis",
                                 "RelRes_And_RelDis"),
                               tuple<INPAR::STR::ConvCheck>(
                                 INPAR::STR::convcheck_absres_or_absdis,
                                 INPAR::STR::convcheck_absres_and_absdis,
                                 INPAR::STR::convcheck_relres_or_absdis,
                                 INPAR::STR::convcheck_relres_and_absdis,
                                 INPAR::STR::convcheck_relres_or_reldis,
                                 INPAR::STR::convcheck_relres_and_reldis),
                               &sdyn);

  DoubleParameter("TOLDISP",1.0E-10,
                  "tolerance in the displacement norm for the newton iteration",
                  &sdyn);
  DoubleParameter("TOLRES",1.0E-08,
                  "tolerance in the residual norm for the newton iteration",
                  &sdyn);
  DoubleParameter("TOLCONSTR",1.0E-08,
                  "tolerance in the constr error norm for the newton iteration",
                  &sdyn);
  IntParameter("MAXITER",50,
               "maximum number of iterations allowed for newton iteration before failure",
               &sdyn);

  setStringToIntegralParameter<INPAR::STR::NonlinSolTech>("NLNSOL","fullnewton","",
                               tuple<std::string>(
                                 "vague",
                                 "fullnewton",
                                 "lsnewton",
                                 "modnewton",
                                 "nlncg",
                                 "ptc",
                                 "newtonlinuzawa",
                                 "augmentedlagrange",
                                 "NoxNewtonLineSearch",
                                 "noxgeneral"),
                               tuple<INPAR::STR::NonlinSolTech>(
                                 INPAR::STR::soltech_vague,
                                 INPAR::STR::soltech_newtonfull,
                                 INPAR::STR::soltech_newtonls,
                                 INPAR::STR::soltech_newtonmod,
                                 INPAR::STR::soltech_nlncg,
                                 INPAR::STR::soltech_ptc,
                                 INPAR::STR::soltech_newtonuzawalin,
                                 INPAR::STR::soltech_newtonuzawanonlin,
                                 INPAR::STR::soltech_noxnewtonlinesearch,
                                 INPAR::STR::soltech_noxgeneral),
                               &sdyn);

  setStringToIntegralParameter<INPAR::STR::PredEnum>("PREDICT","ConstDis","",
                               tuple<std::string>(
                                 "Vague",
                                 "ConstDis",
                                 "ConstDisVelAcc",
                                 "TangDis"),
                               tuple<INPAR::STR::PredEnum>(
                                 INPAR::STR::pred_vague,
                                 INPAR::STR::pred_constdis,
                                 INPAR::STR::pred_constdisvelacc,
                                 INPAR::STR::pred_tangdis),
                               &sdyn);

  // time adaptivity (old style)
  IntParameter("TIMEADAPT",0,"",&sdyn);
  IntParameter("ITWANT",0,"",&sdyn);
  DoubleParameter("MAXDT",0.0,"",&sdyn);
  DoubleParameter("RESULTDT",0.0,"",&sdyn);
  // Uzawa iteration for constraint systems
  DoubleParameter("UZAWAPARAM",1.0,"Parameter for Uzawa algorithm dealing with lagrange multipliers",&sdyn);
  DoubleParameter("UZAWATOL",1.0E-8,"Tolerance for iterative solve with Uzawa algorithm",&sdyn);
  IntParameter("UZAWAMAXITER",50,"maximum number of iterations allowed for uzawa algorithm before failure going to next newton step",&sdyn);
  setStringToIntegralParameter<INPAR::STR::ConSolveAlgo>("UZAWAALGO","iterative","",
                                 tuple<std::string>(
                                   "iterative",
                                   "direct"),
                                 tuple<INPAR::STR::ConSolveAlgo>(
                                   INPAR::STR::consolve_iterative,
                                   INPAR::STR::consolve_direct),
                                 &sdyn);

  // convergence criteria adaptivity
  setStringToIntegralParameter<int>("ADAPTCONV","No",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&sdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&sdyn);

  /*--------------------------------------------------------------------*/
  /* parameters for time step size adaptivity in structural dynamics */
  Teuchos::ParameterList& tap = sdyn.sublist("TIMEADAPTIVITY",false,"");
  SetValidTimeAdaptivityParameters(tap);

  /*----------------------------------------------------------------------*/
  /* parameters for generalised-alpha structural integrator */
  Teuchos::ParameterList& genalpha = sdyn.sublist("GENALPHA",false,"");

  setStringToIntegralParameter<INPAR::STR::MidAverageEnum>("GENAVG","ImrLike",
                               "mid-average type of internal forces",
                               tuple<std::string>(
                                 "Vague",
                                 "ImrLike",
                                 "TrLike"),
                               tuple<INPAR::STR::MidAverageEnum>(
                                 INPAR::STR::midavg_vague,
                                 INPAR::STR::midavg_imrlike,
                                 INPAR::STR::midavg_trlike),
                               &genalpha);
  DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,1/2]",&genalpha);
  DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&genalpha);
  DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
  DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&genalpha);

  /*----------------------------------------------------------------------*/
  /* parameters for one-step-theta structural integrator */
  Teuchos::ParameterList& onesteptheta = sdyn.sublist("ONESTEPTHETA",false,"");

  DoubleParameter("THETA",0.5,"One-step-theta factor in (0,1]",&onesteptheta);


  /*----------------------------------------------------------------------*/
  /* parameters for generalised-energy-momentum structural integrator */
  Teuchos::ParameterList& gemm = sdyn.sublist("GEMM",false,"");

  DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&gemm);
  DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&gemm);
  DoubleParameter("XI",0.0,"generalisation factor in [0,1)",&gemm);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& iap = list->sublist("INVERSE ANALYSIS",false,"");

  // Inverse Analysis
  setStringToIntegralParameter<int>("INV_ANALYSIS","No",
                               "determines the material parameter for the lung material",
                               yesnotuple,yesnovalue,
                               &iap);
  // Measured displacement/load curve during the experiments  a1*(1-exp(-pow((a2*t), a3)))
  DoubleParameter("MEASURED_CURVE0",0.0,"measured displacment of the tension testing",&iap);
  DoubleParameter("MEASURED_CURVE1",0.0,"measured displacment of the tension testing",&iap);
  DoubleParameter("MEASURED_CURVE2",0.0,"measured displacment of the tension testing",&iap);

  // tolerance for inv_analysis
  DoubleParameter("INV_ANA_TOL",1.0,"tolerance for inverse analysis",&iap);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scontact = list->sublist("STRUCTURAL CONTACT",false,"");

  setStringToIntegralParameter<INPAR::CONTACT::ContactType>("CONTACT","None","Type of structural contact",
                               tuple<std::string>("None","none",
                                                  "Normal","normal",
                                                  "Frictional","frictional",
                                                  "MeshTying","Meshtying","meshtying"),
                               tuple<INPAR::CONTACT::ContactType>(INPAR::CONTACT::contact_none,INPAR::CONTACT::contact_none,
                                          INPAR::CONTACT::contact_normal,INPAR::CONTACT::contact_normal,
                                          INPAR::CONTACT::contact_frictional,INPAR::CONTACT::contact_frictional,
                                          INPAR::CONTACT::contact_meshtying,INPAR::CONTACT::contact_meshtying,
                                          INPAR::CONTACT::contact_meshtying),
                               &scontact);

  setStringToIntegralParameter<int>("BASISTRAFO","No","If chosen basis transformation is applied to displacements",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<INPAR::CONTACT::ContactFrictionType>("FRICTION","None","Type of friction law",
                                tuple<std::string>("None","none",
                                                   "Stick","stick",
                                                   "Tresca","tresca",
                                                   "Coulomb","coulomb"),
                                tuple<INPAR::CONTACT::ContactFrictionType>(INPAR::CONTACT::friction_none,INPAR::CONTACT::friction_none,
                                           INPAR::CONTACT::friction_stick,INPAR::CONTACT::friction_stick,
                                           INPAR::CONTACT::friction_tresca,INPAR::CONTACT::friction_tresca,
                                           INPAR::CONTACT::friction_coulomb,INPAR::CONTACT::friction_coulomb),
                                &scontact);

  DoubleParameter("FRBOUND",0.0,"Friction bound for Tresca friction",&scontact);
  DoubleParameter("FRCOEFF",0.0,"Friction coefficient for Coulomb friction",&scontact);

  setStringToIntegralParameter<int>("FULL_LINEARIZATION","No","If chosen full linearization of contact is applied",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON","No","If chosen semi-smooth Newton concept is applied",
                                 yesnotuple,yesnovalue,&scontact);

  DoubleParameter("SEMI_SMOOTH_CN",0.0,"Weighting factor cn for semi-smooth PDASS",&scontact);
  DoubleParameter("SEMI_SMOOTH_CT",0.0,"Weighting factor ct for semi-smooth PDASS",&scontact);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& bromotion = list->sublist("BROWNIAN MOTION",false,"");
  setStringToIntegralParameter<int>("BROWNIAN_MOTION","No",
                                "determines whether stochastical forces should be considered",
                                 yesnotuple,yesnovalue,
                                 &bromotion);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& statmech = list->sublist("STATISTICAL MECHANICS",false,"");

  //Reading kind of background fluid stream in the thermal bath
  setStringToIntegralParameter<INPAR::STATMECH::ThermalBathType>("THERMALBATH","None","Type of thermal bath applied to elements",
                               //listing possible strings in input file in category THERMALBATH
                               tuple<std::string>("None","none",
                                                  "Uniform","uniform",
                                                  "ShearFlow","shearflow","Shearflow"),
                               //translating input strings into BACI input parameters
                               tuple<INPAR::STATMECH::ThermalBathType>(INPAR::STATMECH::thermalbath_none,INPAR::STATMECH::thermalbath_none,
                                          INPAR::STATMECH::thermalbath_uniform,INPAR::STATMECH::thermalbath_uniform,
                                          INPAR::STATMECH::thermalbath_shearflow,INPAR::STATMECH::thermalbath_shearflow,INPAR::STATMECH::thermalbath_shearflow),
                               &statmech);
  //Reading which kind of special output should be written to files
  setStringToIntegralParameter<INPAR::STATMECH::StatOutput>("SPECIAL_OUTPUT","None","kind of special statistical output data written into files",
                                 //listing possible strings in input file in category SPECIAL_OUTPUT
                                 tuple<std::string>("None","none",
                                                    "EndToEnd_Log","endtoend_log","EndtoEnd_log",
                                                    "EndToEnd_Ergodicity","endtoend_ergodicity",
                                                    "Viscoelasticity","viscoelasticity","ViscoElasticity",
                                                    "Gmsh","gmsh"),
                                 //translating input strings into BACI input parameters
                                 tuple<INPAR::STATMECH::StatOutput>(INPAR::STATMECH::statout_none,INPAR::STATMECH::statout_none,
                                            INPAR::STATMECH::statout_endtoendlog,INPAR::STATMECH::statout_endtoendlog,INPAR::STATMECH::statout_endtoendlog,
                                            INPAR::STATMECH::statout_endtoendergodicity,INPAR::STATMECH::statout_endtoendergodicity,
                                            INPAR::STATMECH::statout_viscoelasticity,INPAR::STATMECH::statout_viscoelasticity,INPAR::STATMECH::statout_viscoelasticity,
                                            INPAR::STATMECH::statout_gmsh,INPAR::STATMECH::statout_gmsh),
                                 &statmech);
  //percentage of total simulation time after which writing of statistical output is started
  DoubleParameter("START_FACTOR",0.0,"Percentage of total simulation time after which writing of statistical output is started",&statmech);
  //Reading whether dynamics remodelling of cross linker distribution takes place
  setStringToIntegralParameter<int>("DYN_CROSSLINKERS","No","If chosen cross linker proteins are added and removed in each time step",
                               yesnotuple,yesnovalue,&statmech);
  //Reading double parameter for gradient of flow field
  DoubleParameter("GRADIENT",0.0,"Velocity gradient of shear flow",&statmech);
  //Reading double parameter for viscosity of background fluid
  DoubleParameter("ETA",0.0,"viscosity",&statmech);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&statmech);
  //Reading double parameter for crosslinker off-rate
  DoubleParameter("K_ON",0.0,"crosslinker on-rate",&statmech);
  //Reading double parameter for crosslinker off-rate
  DoubleParameter("K_OFF",0.0,"crosslinker off-rate",&statmech);
  //Reading double parameter for maximal cross linker protein length
  DoubleParameter("R_LINK",0.0,"Maximal distance between two nodes connected by a crosslinker",&statmech);
  //Reading double parameter for concentration of crosslinking protein
  DoubleParameter("C_CROSSLINKER",0.0,"Molar concentration of crosslinking protein",&statmech);
  //order of interpolation for stochastical fields
  IntParameter("STOCH_ORDER",0,"order of interpolation for stochastical fields",&statmech);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn = list->sublist("FLUID DYNAMIC",false,"");

  setStringToIntegralParameter<int>("LOWMACH","No",
                               "low-mach-number or incompressible flow",
                               tuple<std::string>(
                                 "No",
                                 "Yes"
                                 ),
                               tuple<int>(0,1),
                               &fdyn);

  setStringToIntegralParameter<int>("DYNAMICTYP","Nlin_Time_Int",
                               "Nonlinear Time Integraton Scheme",
                               tuple<std::string>(
                                 "Nlin_Time_Int",
                                 "Lin_Time_Int"
                                 ),
                               tuple<int>(
                                dyntyp_nln_time_int,
                                dyntyp_lin_time_int
                                ),
                               &fdyn);

  setStringToIntegralParameter<FLUID_TIMEINTTYPE>("TIMEINTEGR","One_Step_Theta",
                               "Time Integration Scheme",
                               tuple<std::string>(
                                 "Stationary",
                                 "Gen_Alfa",
                                 "Gen_Alpha",
                                 "Af_Gen_Alpha",
                                 "One_Step_Theta",
                                 "BDF2",
                                 "Inc_Acc_Gen_Alpha",
                                 "Theta_Adamsbashforth"
                                 ),
                               tuple<FLUID_TIMEINTTYPE>(
                                 timeint_stationary,
                                 timeint_gen_alpha,
                                 timeint_gen_alpha,
                                 timeint_afgenalpha,
                                 timeint_one_step_theta,
                                 timeint_bdf2,
                                 timeint_inc_acc_gen_alpha,
                                 timeint_theta_adamsbashforth
                                 ),
                               &fdyn);
  setStringToIntegralParameter<int>("STARTINGALGO","One_Step_Theta","",
                               tuple<std::string>(
                                 "One_Step_Theta"
                                 ),
                               tuple<int>(
                                 timeint_one_step_theta
                                 ),
                               &fdyn);
  setStringToIntegralParameter<int>("NONLINITER","fixed_point_like",
                               "Nonlinear iteration scheme",
                               tuple<std::string>(
                                 "fixed_point_like",
                                 "Newton",
                                 "minimal"
                                 ),
                               tuple<int>(1,2,3),
                               &fdyn);

  setStringToIntegralParameter<int>("PREDICTOR","steady_state_predictor",
                               "Predictor for first guess in nonlinear iteration",
                               tuple<std::string>(
                                 "steady_state_predictor",
                                 "zero_acceleration_predictor",
                                 "constant_acceleration_predictor",
                                 "constant_increment_predictor"
                                 ),
                               tuple<int>(1,2,3,4),
                               &fdyn);

  setStringToIntegralParameter<int>("CONVCHECK","L_2_norm",
                               "norm for convergence check",
                               tuple<std::string>(
                                 "No",
                                 "L_infinity_norm",
                                 "L_1_norm",
                                 "L_2_norm",
                                 "L_2_norm_without_residual_at_itemax"
                                 ),
                               tuple<std::string>(
                                 "do not check for convergence (ccarat)",
                                 "use max norm (ccarat)",
                                 "use abs. norm (ccarat)",
                                 "compute L2 errors of increments (reltive) and residuals (absolute)",
                                 "same as L_2_norm, only no residual norm is computed if itemax is reached (speedup for turbulence calculations, startup phase)"
                                 ),
                               tuple<int>(
                                 FLUID_DYNAMIC::fncc_no,
                                 FLUID_DYNAMIC::fncc_Linf,
                                 FLUID_DYNAMIC::fncc_L1,
                                 FLUID_DYNAMIC::fncc_L2,
                                 FLUID_DYNAMIC::fncc_L2_wo_res
                                 ),
                               &fdyn);
  setStringToIntegralParameter<int>("STEADYCHECK","L_2_norm",
                               "Norm of steady state check",
                               tuple<std::string>(
                                 "No",
                                 "L_infinity_norm",
                                 "L_1_norm",
                                 "L_2_norm"
                                 ),
                               tuple<int>(
                                 FLUID_DYNAMIC::fncc_no,
                                 FLUID_DYNAMIC::fncc_Linf,
                                 FLUID_DYNAMIC::fncc_L1,
                                 FLUID_DYNAMIC::fncc_L2
                                 ),
                               &fdyn);
  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Starting Field",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_from_file",
                                 "field_by_function",
                                 "disturbed_field_from_function",
                                 "SOLWAVE",
                                 "WAVEBREAKING",
                                 "BELTRAMI-FLOW",
                                 "KIM-MOIN-FLOW",
                                 "BREAKING-DAM"),
                               tuple<int>(0,1,2,3,6,7,8,9,10),
                               &fdyn);

  setStringToIntegralParameter<int>("LIFTDRAG","No",
                               "Calculate lift and drag forces along specified boundary",
                               tuple<std::string>(
                                 "No",
                                 "no",
                                 "Yes",
                                 "yes",
                                 "Nodeforce",
                                 "NODEFORCE",
                                 "nodeforce"
                                 ),
                               tuple<int>(
                                 FLUID_DYNAMIC::ld_none,
                                 FLUID_DYNAMIC::ld_none,
                                 FLUID_DYNAMIC::ld_nodeforce,
                                 FLUID_DYNAMIC::ld_nodeforce,
                                 FLUID_DYNAMIC::ld_nodeforce,
                                 FLUID_DYNAMIC::ld_nodeforce,
                                 FLUID_DYNAMIC::ld_nodeforce
                                 ),
                               &fdyn);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(0,1),
                               &fdyn);

  setStringToIntegralParameter<int>("FSSUGRVISC","No","fine-scale subgrid viscosity",
                               tuple<std::string>(
                                 "No",
                                 "artificial_all",
                                 "artificial_small",
                                 "Smagorinsky_all",
                                 "Smagorinsky_small"
                                 ),
                               tuple<int>(0,1,2,3,4),
                               &fdyn);

  setStringToIntegralParameter<int>("SIMPLER","no",
                               "Switch on SIMPLE family of solvers, needs additional FLUID PRESSURE SOLVER block!",
                               yesnotuple,yesnovalue,&fdyn);

  setStringToIntegralParameter<int>("ADAPTCONV","yes",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&fdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&fdyn);

  IntParameter("UPPSS",1,"Increment for visualisation (unused)",&fdyn);
  IntParameter("UPOUT",1,"Increment for writing solution to output file",&fdyn);
  IntParameter("UPRES",1,"Increment for writing solution",&fdyn);
  IntParameter("RESSTEP",0,"Restart Step",&fdyn);
  IntParameter("RESTARTEVRY",20,"Increment for writing restart",&fdyn);
  IntParameter("NUMSTEP",1,"Total number of Timesteps",&fdyn);
  IntParameter("STEADYSTEP",-1,"steady state check every step",&fdyn);
  IntParameter("NUMSTASTEPS",0,"Number of Steps for Starting Scheme",&fdyn);
  IntParameter("STARTFUNCNO",-1,"Function for Initial Starting Field",&fdyn);
  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&fdyn);
  IntParameter("GRIDVEL",1,"order of accuracy of mesh velocity determination",&fdyn);
  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&fdyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&fdyn);
  DoubleParameter("ALPHA_M",1.0,"Time integration factor",&fdyn);
  DoubleParameter("ALPHA_F",1.0,"Time integration factor",&fdyn);
  DoubleParameter("GAMMA",1.0,"Time integration factor",&fdyn);

  DoubleParameter("THETA",0.66,"Time integration factor",&fdyn);

  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&fdyn);
  DoubleParameter("STEADYTOL",1e-6,"Tolerance for steady state check",&fdyn);
  DoubleParameter("START_THETA",1.0,"Time integraton factor for starting scheme",&fdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_stab = fdyn.sublist("STABILIZATION",false,"");

  // this parameter seperates stabilized from unstabilized methods
  setStringToIntegralParameter<int>("STABTYPE",
                               "residual_based",
                               "Apply (un)stabilized fluid formulation",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "inconsistent",
                                 "residual_based"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> inf-sup stable elements required!",
                                 "Similar to residual based without second derivatives (i.e. only consistent for tau->0, but faster)",
                                 "Use a residual-based stabilization or, more generally, a stabilization \nbased on the concept of the residual-based variational multiscale method...\nExpecting additional input")  ,
                               tuple<int>(0,1,2),
                               &fdyn_stab);

  // the following parameters are necessary only if a residual based stabilized method is applied
  setStringToIntegralParameter<int>("TDS",
                               "quasistatic",
                               "Flag to allow time dependency of subscales for residual-based stabilization.",
                               tuple<std::string>(
                                 "quasistatic",
                                 "time_dependent"),
                               tuple<std::string>(
                                 "Use a quasi-static residual-based stabilization (standard case)",
                                 "Residual-based stabilization including time evolution equations for subscales"),
                               tuple<int>(0,1),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("TRANSIENT",
                               "no_transient",
                               "Specify how to treat the transient term.",
                               tuple<std::string>(
                                 "no_transient",
                                 "yes_transient",
                                 "transient_complete"),
                               tuple<std::string>(
                                 "Do not use transient term (currently only opportunity for quasistatic stabilization)",
                                 "Use transient term (recommended for time dependent subscales)",
                                 "Use transient term including a linearisation of 1/tau"),
                               tuple<int>(0,1,2),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("PSPG",
                               "yes_pspg",
                               "Flag to (de)activate PSPG.",
                               tuple<std::string>(
                                 "no_pspg",
                                 "yes_pspg"),
                               tuple<std::string>(
                                 "No PSPG -> inf-sup-stable elements mandatory",
                                 "Use PSPG -> allowing for equal-order interpolation"),
                               tuple<int>(0,1),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("SUPG",
                               "yes_supg",
                               "Flag to (de)activate SUPG.",
                               tuple<std::string>(
                                 "no_supg",
                                 "yes_supg"),
                               tuple<std::string>(
                                 "No SUPG",
                                 "Use SUPG."),
                               tuple<int>(0,1),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("VSTAB",
                               "no_vstab",
                               "Flag to (de)activate viscous term in residual-based stabilization.",
                               tuple<std::string>(
                                 "no_vstab",
                                 "vstab_gls",
                                 "vstab_gls_rhs",
                                 "vstab_usfem",
                                 "vstab_usfem_rhs"
                                 ),
                               tuple<std::string>(
                                 "No viscous term in stabilization",
                                 "Viscous stabilization of GLS type",
                                 "Viscous stabilization of GLS type, included only on the right hand side",
                                 "Viscous stabilization of USFEM type",
                                 "Viscous stabilization of USFEM type, included only on the right hand side"
                                 ),
                               tuple<int>(0,1,2,3,4),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("CSTAB",
                               "cstab_qs",
                               "Flag to (de)activate least-squares stabilization of continuity equation.",
                               tuple<std::string>(
                                 "no_cstab",
                                 "cstab_qs"),
                               tuple<std::string>(
                                 "No continuity stabilization",
                                 "Quasistatic continuity stabilization"),
                               tuple<int>(0,1),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("CROSS-STRESS",
                               "no_cross",
                               "Flag to (de)activate cross-stress term -> residual-based VMM.",
                               tuple<std::string>(
                                 "no_cross",
                                 "cross_complete",
                                 "cross_rhs"
                                 ),
                               tuple<std::string>(
                                 "No cross-stress term",
                                 "Include the cross-stress term with a linearization of the convective part",
                                 "Include cross-stress term, but only explicitly on right hand side"
                                 ),
                               tuple<int>(0,1,2),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("REYNOLDS-STRESS",
                               "no_reynolds",
                               "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
                               tuple<std::string>(
                                 "no_reynolds",
                                 "reynolds_rhs",
                                 "reynolds_complete"
                                 ),
                               tuple<std::string>(
                                 "No Reynolds-stress term",
                                 "Include Reynolds-stress term explicitly on right hand side",
                                 "Include Reynolds-stress term with linearisation"
                                 ),
                               tuple<int>(0,1,2),
                               &fdyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU",
                               "Barrenechea_Franca_Valentin_Wall",
                               "Definition of tau_M,C",
                               tuple<std::string>(
                                 "Barrenechea_Franca_Valentin_Wall",
                                 "Smoothed_FBVW",
                                 "Franca_Barrenechea_Valentin_Codina",
                                 "Bazilevs",
                                 "Codina"),
                               tuple<std::string>(
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea; tau_C: Wall",
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea (smoothed max opertaor using exp function); tau_C: Wall",
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea; tau_C: Codina"  ,
                                 "tau_M and tau_C (Bazilevs, based on G_ij and g_i)",
                                 "tau_M and tau_C: Codina")  ,
                                    tuple<int>(0,1,2,3,4),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("OUTFLOW_STAB",
                               "no_outstab",
                               "Flag to (de)activate outflow stabilization term",
                               tuple<std::string>(
                                 "no_outstab",
                                 "yes_outstab"),
                               tuple<std::string>(
                                 "No outflow stabilization term",
                                 "Add outflow stabilization term."),
                               tuple<int>(0,1),
                               &fdyn_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_turbu = fdyn.sublist("TURBULENCE MODEL",false,"");

  setStringToIntegralParameter<int>(
    "TURBULENCE_APPROACH",
    "DNS_OR_RESVMM_LES",
    "There are several options to deal with turbulent flows.",
    tuple<std::string>(
      "DNS_OR_RESVMM_LES",
      "CLASSICAL_LES",
      "RANS"),
    tuple<std::string>(
      "Try to solve flow as an underresolved DNS.\nMind that your stabilisation already acts as a kind of turbulence model!",
      "Perform a classical Large Eddy Simulation adding \naddititional turbulent viscosity. This may be based on various physical models.",
      "Solve Reynolds averaged Navier Stokes using an \nalgebraic, one- or two equation closure.\nNot implemented yet."),
    tuple<int>(0,1,2),
    &fdyn_turbu);

  setStringToIntegralParameter<int>(
    "PHYSICAL_MODEL",
    "no_model",
    "Classical LES approaches require an additional model for\nthe turbulent viscosity.",
    tuple<std::string>(
      "no_model",
      "Smagorinsky",
      "Smagorinsky_with_van_Driest_damping",
      "Dynamic_Smagorinsky"),
    tuple<std::string>(
      "If classical LES is our turbulence approach, this is a contradiction and should cause a dserror.",
      "Classical constant coefficient Smagorinsky model. Be careful if you \nhave a wall bounded flow domain!",
      "Use an exponential damping function for the turbulent viscosity \nclose to the wall. This is only implemented for a channel geometry of \nheight 2 in y direction. The viscous lengthscale l_tau is \nrequired as additional input.",
      "The solution is filtered and by comparison of the filtered \nvelocity field with the real solution, the Smagorinsky constant is \nestimated in each step --- mind that this procedure includes \nan averaging in the xz plane, hence this implementation will only work \nfor a channel flow."),
    tuple<int>(0,1,2,3),
    &fdyn_turbu);

  DoubleParameter("C_SMAGORINSKY",0.0,"Constant for the Smagorinsky model. Something between 0.1 to 0.24",&fdyn_turbu);

  setStringToIntegralParameter<int>(
    "CANONICAL_FLOW",
    "no",
    "Sampling is different for different canonical flows \n--- so specify what kind of flow you've got",
    tuple<std::string>(
      "no",
      "channel_flow_of_height_2",
      "lid_driven_cavity",
      "square_cylinder",
      "square_cylinder_nurbs",
      "loma_channel_flow_of_height_2",
      "loma_lid_driven_cavity"),
    tuple<std::string>(
      "The flow is not further specified, so spatial averaging \nand hence the standard sampling procedure is not possible",
      "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.",
      "For this flow, all statistical data are evaluated on the center lines of the xy-midplane, averaged only over time.",
      "For this flow, statistical data are evaluated on various lines of the xy-midplane, averaged only over time.",
      "For this flow, statistical data are evaluated on various lines of the xy-midplane, averaged over time and eventually in one hom.direction.",
      "For this low-Mach-number flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.",
      "For this low-Mach-number flow, all statistical data are evaluated on the center lines of the xy-midplane, averaged only over time."),
    tuple<int>(0,1,2,3,4,5,6),
    &fdyn_turbu);

  setStringToIntegralParameter<int>(
    "HOMDIR",
    "not_specified",
    "Specify the homogenous direction(s) of a flow",
    tuple<std::string>(
      "not_specified",
      "x"            ,
      "y"            ,
      "z"            ,
      "xy"           ,
      "xz"           ,
      "yz"           ),
    tuple<std::string>(
      "no homogeneous directions available, averaging is restricted to time averaging",
      "average along x-direction"                                                     ,
      "average along y-direction"                                                     ,
      "average along z-direction"                                                     ,
      "Wall normal direction is z, average in x and y direction"                      ,
      "Wall normal direction is y, average in x and z direction (standard case)"      ,
      "Wall normal direction is x, average in y and z direction"                      ),
    tuple<int>(0,1,2,3,4,5,6),
    &fdyn_turbu);

  DoubleParameter(
    "CHANNEL_L_TAU",
    0.0,
    "Used for normalisation of the wall normal distance in the Van \nDriest Damping function. May be taken from the output of \nthe apply_mesh_stretching.pl preprocessing script.",
    &fdyn_turbu);

  DoubleParameter(
    "CHAN_AMPL_INIT_DIST",
    0.1,
    "Max. amplitude of the random disturbance in percent of the initial value in mean flow direction.",
    &fdyn_turbu);

  IntParameter("SAMPLING_START",1,"Time step after when sampling shall be started",&fdyn_turbu);
  IntParameter("SAMPLING_STOP",1,"Time step when sampling shall be stopped",&fdyn_turbu);
  IntParameter("DUMPING_PERIOD",1,"Period of time steps after which statistical data shall be dumped",&fdyn_turbu);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC",false,"");

  DoubleParameter("TIMESTEP",0.1,"",&adyn);
  IntParameter("NUMSTEP",41,"",&adyn);
  DoubleParameter("MAXTIME",4.0,"",&adyn);
  setStringToIntegralParameter<int>("ALE_TYPE","classic_lin","ale mesh movement algorithm",
                               tuple<std::string>("classic_lin","incr_lin","springs"),
                               tuple<int>(ALE_DYNAMIC::classic_lin,
                                          ALE_DYNAMIC::incr_lin,
                                          ALE_DYNAMIC::springs),
                               &adyn);
  IntParameter("NUM_INITSTEP",0,"",&adyn);
  IntParameter("RESEVRYDISP",1,"",&adyn);

  setStringToIntegralParameter<int>("QUALITY","none","unused",
                               tuple<std::string>("none","NONE"),
                               tuple<int>(
                                 ALE_DYNAMIC::no_quality,
                                 ALE_DYNAMIC::no_quality),
                               &adyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn = list->sublist(
      "SCALAR TRANSPORT DYNAMIC",
      false,
      "control parameters for scalar transport problems\n");

  setStringToIntegralParameter<INPAR::SCATRA::TimeIntegrationScheme>("TIMEINTEGR","One_Step_Theta",
                               "Time Integration Scheme",
                               tuple<std::string>(
                                 "Stationary",
                                 "One_Step_Theta",
                                 "BDF2",
                                 "Gen_Alpha"
                                 ),
                               tuple<INPAR::SCATRA::TimeIntegrationScheme>(
                                   INPAR::SCATRA::timeint_stationary,
                                   INPAR::SCATRA::timeint_one_step_theta,
                                   INPAR::SCATRA::timeint_bdf2,
                                   INPAR::SCATRA::timeint_gen_alpha
                                 ),
                               &scatradyn);

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&scatradyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&scatradyn);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&scatradyn);
  IntParameter("ITEMAX",10,"Maximum number of nonlinear iterations",&scatradyn);
  DoubleParameter("THETA",0.5,"One-step-theta time integration factor",&scatradyn);
  DoubleParameter("ALPHA_M",0.5,"Generalized-alpha time integration factor",&scatradyn);
  DoubleParameter("ALPHA_F",0.5,"Generalized-alpha time integration factor",&scatradyn);
  DoubleParameter("GAMMA",0.5,"Generalized-alpha time integration factor",&scatradyn);
  //IntParameter("WRITESOLEVRY",1,"Increment for writing solution",&scatradyn);
  IntParameter("UPRES",1,"Increment for writing solution",&scatradyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&scatradyn);

  setStringToIntegralParameter<int>("VELOCITYFIELD","zero",
                               "type of velocity field used for scalar tranport problems",
                               tuple<std::string>(
                                 "zero",
                                 "function",
                                 "Navier_Stokes"
                                 ),
                               tuple<int>(0,1,2),
                               &scatradyn);

  IntParameter("VELFUNCNO",-1,"function number for scalar transport velocity field",&scatradyn);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Field for scalar transport problem",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function",
                                 "field_by_condition"),
                               tuple<int>(0,1,2),
                               &scatradyn);

  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&scatradyn);

  setStringToIntegralParameter<int>("CALCERROR","No",
                               "compute error compared to analytical solution",
                               tuple<std::string>(
                                 "No",
                                 "Kwok_Wu"
                                 ),
                               tuple<int>(0,1),
                               &scatradyn);

  setStringToIntegralParameter<int>("WRITEFLUX","No","output of diffusive/total flux vectors",
                               tuple<std::string>(
                                 "No",
                                 "totalflux_domain",
                                 "diffusiveflux_domain",
                                 "totalflux_boundary",
                                 "diffusiveflux_boundary"
                                 ),
                               tuple<int>(0,1,2,3,4),
                               &scatradyn);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(0,1),
                               &scatradyn);

  setStringToIntegralParameter<int>("FSSUGRVISC","No","fine-scale subgrid diffusivity",
                               tuple<std::string>(
                                 "No",
                                 "artificial_all",
                                 "artificial_small"
                                 ),
                               tuple<int>(0,1,2),
                               &scatradyn);

  setStringToIntegralParameter<int>("BLOCKPRECOND","no",
                               "Switch to block-preconditioned family of solvers, needs additional SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER block!",
                               yesnotuple,yesnovalue,&scatradyn);


  Teuchos::ParameterList& scatra_nonlin = scatradyn.sublist(
      "NONLINEAR",
      false,
      "control parameters for solving nonlinear SCATRA problems\n");

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&scatra_nonlin);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&scatra_nonlin);
  setStringToIntegralParameter<int>("EXPLPREDICT","yes",
                               "do an explicit predictor step before starting nonlinear iteration",
                               yesnotuple,yesnovalue,&scatra_nonlin);
  // convergence criteria adaptivity
  setStringToIntegralParameter<int>("ADAPTCONV","yes",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&scatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&scatra_nonlin);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn_stab = scatradyn.sublist("STABILIZATION",false,"");

  // this parameter seperates stabilized from unstabilized methods
  setStringToIntegralParameter<int>("STABTYPE",
                               "residual_based",
                               "Apply (un)stabilized scalar transport formulation",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "residual_based"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> only reasonable for low-Peclet-number flows",
                                 "Use a residual-based stabilization")  ,
                               tuple<int>(0,1),
                               &scatradyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU",
                               "Franca_Valentin",
                               "Definition of tau",
                               tuple<std::string>(
                                 "Franca_Valentin",
                                 "Bazilevs"),
                               tuple<std::string>(
                                 "tau according to Franca and Valentin (2000)",
                                 "tau according to Bazilevs et al. (2007) (based on G_ij and g_i)")  ,
                               tuple<int>(0,1),
                               &scatradyn_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& lomacontrol = list->sublist(
      "LOMA CONTROL",
      false,
      "control parameters for low-Mach-number flow problems\n");

  IntParameter("NUMSTEP",24,"Total number of time steps",&lomacontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&lomacontrol);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&lomacontrol);
  IntParameter("ITEMAX",10,"Maximum number of outer iterations",&lomacontrol);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&lomacontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&lomacontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&lomacontrol);
  DoubleParameter("THERMOPRESS",98100.0,"(initial) thermodynamic pressure",&lomacontrol);
  DoubleParameter("GASCONSTANT",287.0,"specific gas constant R (in J/(kg*K))",&lomacontrol);
  setStringToIntegralParameter<int>("CONSTHERMPRESS","Yes",
                               "treatment of thermodynamic pressure in time",
                               tuple<std::string>(
                                 "No_energy",
                                 "No_mass",
                                 "Yes"
                                 ),
                               tuple<int>(0,1,2),
                               &lomacontrol);
  setStringToIntegralParameter<int>("OUTMEAN","No",
                               "print out mean values of temperature/density",
                               tuple<std::string>(
                                 "No",
                                 "Yes"
                                 ),
                               tuple<int>(0,1),
                               &lomacontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& elchcontrol = list->sublist(
      "ELCH CONTROL",
      false,
      "control parameters for electrochemistry problems\n");

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&elchcontrol);
  IntParameter("NUMSTEP",24,"Total number of time steps",&elchcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&elchcontrol);
  IntParameter("ITEMAX",10,"Maximum number of nonlinear iterations",&elchcontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&elchcontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&elchcontrol);
  BoolParameter("OUTMEAN","No","Output of total and mean values",&elchcontrol);
  DoubleParameter("TEMPERATURE",298.0,"Constant temperature (Kelvin)",&elchcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrol = list->sublist("COMBUSTION CONTROL",false,"");

  DoubleParameter("MAXTIME",10.0,"Total simulation time",&combustcontrol);
  IntParameter("NUMSTEP",100,"Total number of Timesteps",&combustcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&combustcontrol);
  IntParameter("ITEMAX",10,"Total number of FG iterations",&combustcontrol);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields",&combustcontrol);
  IntParameter("RESTARTEVRY",20,"Increment for writing restart",&combustcontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&combustcontrol);
  setStringToIntegralParameter<FLUID_TIMEINTTYPE>("TIMEINTEGR","One_Step_Theta","Time Integration Scheme",
                               tuple<std::string>(
                                 "Stationary",
                                 "One_Step_Theta"),
                               tuple<FLUID_TIMEINTTYPE>(
                                 timeint_stationary,
                                 timeint_one_step_theta),
                               &combustcontrol);
  setStringToIntegralParameter<INPAR::COMBUST::ReInitialActionGfunc>("REINITGFUNCTION","Signed Distance Function","Type of reinitialization level set",
                               tuple<std::string>(
                                 "Function",
                                 "Signed Distance Function"),
                               tuple<INPAR::COMBUST::ReInitialActionGfunc>(
                                 INPAR::COMBUST::reinitialize_by_function,
                                 INPAR::COMBUST::compute_signeddistancefunction),
                               &combustcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolfluid = combustcontrol.sublist("COMBUSTION FLUID",false,"");

  DoubleParameter("LAMINAR_FLAMESPEED",1.0,"The laminar flamespeed incorporates all chemical kinetics into the problem for now",&combustcontrolfluid);
  DoubleParameter("MARKSTEIN_LENGTH",0.0,"The Markstein length takes flame curvature into account",&combustcontrolfluid);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolgfunc = combustcontrol.sublist("COMBUSTION GFUNCTION",false,"");

  IntParameter("REINITFUNCNO",-1,"function number for reinitialization of G-function field",&combustcontrolgfunc);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fsidyn = list->sublist(
    "FSI DYNAMIC",false,
    "Fluid Structure Interaction\n"
    "Partitioned FSI solver with various coupling methods"
    );

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_AITKEN_rel_param",
                               "Iteration Scheme over the fields",
                               tuple<std::string>(
                                 "basic_sequ_stagg",
                                 //"sequ_stagg_pred",
                                 //"sequ_stagg_shift",
                                 "iter_stagg_fixed_rel_param",
                                 "iter_stagg_AITKEN_rel_param",
                                 "iter_stagg_steep_desc",
                                 "iter_stagg_NLCG",
                                 "iter_stagg_MFNK_FD",
                                 "iter_stagg_MFNK_FSI",
                                 "iter_stagg_MPE",
                                 "iter_stagg_RRE",
                                 "iter_monolithic",
                                 "iter_monolithiclagrange",
                                 "iter_monolithicstructuresplit",
                                 "pseudo_structure"),
                               tuple<int>(
                                 fsi_basic_sequ_stagg,
                                 //fsi_sequ_stagg_pred,
                                 //fsi_sequ_stagg_shift,
                                 fsi_iter_stagg_fixed_rel_param,
                                 fsi_iter_stagg_AITKEN_rel_param,
                                 fsi_iter_stagg_steep_desc,
                                 fsi_iter_stagg_NLCG,
                                 fsi_iter_stagg_MFNK_FD,
                                 fsi_iter_stagg_MFNK_FSI,
                                 fsi_iter_stagg_MPE,
                                 fsi_iter_stagg_RRE,
                                 fsi_iter_monolithic,
                                 fsi_iter_monolithiclagrange,
                                 fsi_iter_monolithicstructuresplit,
                                 fsi_pseudo_structureale),
                                 &fsidyn);

  setStringToIntegralParameter<INPAR::FSI::PartitionedCouplingMethod>(
                               "PARTITIONED","DirichletNeumann",
                               "Coupling strategies for partitioned FSI solvers. Most of the time Dirichlet-Neumann is just right.",
                               tuple<std::string>(
                                 "DirichletNeumann",
                                 "RobinNeumann",
                                 "DirichletRobin",
                                 "RobinRobin"
                                 ),
                               tuple<INPAR::FSI::PartitionedCouplingMethod>(
                                 INPAR::FSI::DirichletNeumann,
                                 INPAR::FSI::RobinNeumann,
                                 INPAR::FSI::DirichletRobin,
                                 INPAR::FSI::RobinRobin
                                 ),
                               &fsidyn);

  DoubleParameter("ALPHA_F",-1.0,"Robin parameter fluid",&fsidyn);
  DoubleParameter("ALPHA_S",-1.0,"Robin parameter structure",&fsidyn);

  setStringToIntegralParameter<int>("DEBUGOUTPUT","No",
                               "Output of unconverged interface values during partitioned FSI iteration.\n"
                               "There will be a new control file for each time step.\n"
                               "This might be helpful to understand the coupling iteration.",
                               yesnotuple,yesnovalue,&fsidyn);

  setStringToIntegralParameter<int>("PREDICTOR","d(n)+dt*v(n)+0.5*dt^2*a(n)",
                               "Predictor for interface displacements",
                               tuple<std::string>(
                                 "d(n)",
                                 "d(n)+dt*(1.5*v(n)-0.5*v(n-1))",
                                 "d(n)+dt*v(n)",
                                 "d(n)+dt*v(n)+0.5*dt^2*a(n)"
                                 ),
                               tuple<int>(1,2,3,4),
                               &fsidyn);

  setStringToIntegralParameter<int>("CONVCRIT","||g(i)||:sqrt(neq)",
                               "Convergence criterium for iteration over fields (unused)",
                               tuple<std::string>(
                                 "||g(i)||:sqrt(neq)",
                                 "||g(i)||:||g(0)||"
                                 ),
                               tuple<int>(1,2),
                               &fsidyn);

  setStringToIntegralParameter<int>("COUPVARIABLE","Displacement",
                               "Coupling variable at the interface",
                               tuple<std::string>("Displacement","Force"),
                               tuple<int>(0,1),
                               &fsidyn);

  setStringToIntegralParameter<int>("ENERGYCHECK","No",
                               "Energy check for iteration over fields",
                               yesnotuple,yesnovalue,&fsidyn);

  setStringToIntegralParameter<int>("IALE","Pseudo_Structure",
                               "Treatment of ALE-field (outdated)",
                               tuple<std::string>(
                                 "Pseudo_Structure"
                                 ),
                               tuple<int>(1),
                               &fsidyn);

  setStringToIntegralParameter<int>("COUPMETHOD","conforming",
                               "Coupling Method Mortar (mtr) or conforming nodes at interface",
                               tuple<std::string>(
                                 "MTR",
                                 "Mtr",
                                 "mtr",
                                 "conforming"
                                 ),
                               tuple<int>(0,0,0,1),
                               &fsidyn);

  setStringToIntegralParameter<int>("COUPFORCE","nodeforce",
                               "Coupling force. Unused. We always couple with nodal forces.",
                               tuple<std::string>(
                                 "none",
                                 "stress",
                                 "nodeforce"
                                 ),
                               tuple<int>(
                                 FSI_DYNAMIC::cf_none,
                                 FSI_DYNAMIC::cf_stress,
                                 FSI_DYNAMIC::cf_nodeforce),
                               &fsidyn);

  setStringToIntegralParameter<int>("SECONDORDER","No",
                               "Second order coupling at the interface.",
                               yesnotuple,yesnovalue,&fsidyn);

  setStringToIntegralParameter<int>("SHAPEDERIVATIVES","No",
                               "Include linearization with respect to mesh movement in Navier Stokes equation.\n"
                               "Supported in monolithic FSI for now.",
                               yesnotuple,yesnovalue,&fsidyn);

  IntParameter("PRECONDREUSE",
               10,
               "Number of preconditioner reused in monolithic FSI",
               &fsidyn);

  setStringToIntegralParameter<INPAR::FSI::LinearBlockSolver>(
                               "LINEARBLOCKSOLVER","PreconditionedKrylov",
                               "Linear solver algorithm for monolithic block system in monolithic FSI.\n"
                               "Most of the time preconditioned Krylov is the right thing to choose. But there are\n"
                               "block Gauss-Seidel methods as well.",
                               tuple<std::string>(
                                 "PreconditionedKrylov",
                                 "PartitionedAitken",
                                 "PartitionedVectorExtrapolation",
                                 "PartitionedJacobianFreeNewtonKrylov",
                                 "BGSAitken",
                                 "BGSVectorExtrapolation",
                                 "BGSJacobianFreeNewtonKrylov"
                                 ),
                               tuple<INPAR::FSI::LinearBlockSolver>(
                                 INPAR::FSI::PreconditionedKrylov,
                                 INPAR::FSI::PartitionedAitken,
                                 INPAR::FSI::PartitionedVectorExtrapolation,
                                 INPAR::FSI::PartitionedJacobianFreeNewtonKrylov,
                                 INPAR::FSI::BGSAitken,
                                 INPAR::FSI::BGSVectorExtrapolation,
                                 INPAR::FSI::BGSJacobianFreeNewtonKrylov
                                 ),
                               &fsidyn);

  IntParameter("ITECHAPP",1,"unused",&fsidyn);
  IntParameter("ICHMAX",1,"unused",&fsidyn);
  IntParameter("ISDMAX",1,"not used up to now",&fsidyn);
  IntParameter("NUMSTEP",200,"Total number of Timesteps",&fsidyn);
  IntParameter("ITEMAX",100,"Maximum number of iterations over fields",&fsidyn);
  IntParameter("UPPSS",1,"Increment for visualisation",&fsidyn);
  IntParameter("UPRES",1,"Increment for writing solution",&fsidyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&fsidyn);

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&fsidyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&fsidyn);
  DoubleParameter("TOLENCHECK",1e-6,"Tolerance for energy check",&fsidyn);
  DoubleParameter("RELAX",1.0,"fixed relaxation parameter",&fsidyn);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields",&fsidyn);
  DoubleParameter("MAXOMEGA",0.0,"largest omega allowed for Aitken relaxation (0.0 means no constraint)",&fsidyn);

  DoubleParameter("BASETOL",1e-3,
                  "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
                  "This tolerance will be used for the linear solve of the FSI block system.\n"
                  "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
                  "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
                  "to the nonlinear convergance test.",
                  &fsidyn);

  DoubleParameter("ADAPTIVEDIST",0.0,
                  "Required distance for adaptive convergence check in Newton-type FSI.\n"
                  "This is the improvement we want to achieve in the linear extrapolation of the\n"
                  "adaptive convergence check. Set to zero to avoid the adaptive check altogether.",
                  &fsidyn);

  // monolithic preconditioner parameter

  DoubleParameter("STRUCTPCOMEGA",1.,
                  "Relaxation factor for Richardson iteration on structural block in MFSI block preconditioner",
                  &fsidyn);
  IntParameter("STRUCTPCITER",0,
               "Number of Richardson iterations on structural block in MFSI block preconditioner",
               &fsidyn);
  DoubleParameter("FLUIDPCOMEGA",1.,
                  "Relaxation factor for Richardson iteration on fluid block in MFSI block preconditioner",
                  &fsidyn);
  IntParameter("FLUIDPCITER",0,
               "Number of Richardson iterations on fluid block in MFSI block preconditioner",
               &fsidyn);

  setStringToIntegralParameter<int>("INFNORMSCALING","Yes","Scale Blocks in Mono-FSI with row infnorm?",
                                     yesnotuple,yesnovalue,&fsidyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL",false,"");

  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT","Yes","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("DLM_CONDENSATION","No","Do you want to condense the distributed Lagrange multiplier?",
                                 yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("CONDEST","No","Do you want to estimate the condition number? It is somewhat costly.",
                                   yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("EXP_INTERSECTION","No","Do you want to use the experimental intersection class?",
                                     yesnotuple,yesnovalue,&xfem_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fluidsolver = list->sublist("FLUID SOLVER",false,"");
  SetValidSolverParameters(fluidsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fluidpsolver = list->sublist("FLUID PRESSURE SOLVER",false,"pressure solver parameters for SIMPLE preconditioning");
  SetValidSolverParameters(fluidpsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& structsolver = list->sublist("STRUCT SOLVER",false,"");
  SetValidSolverParameters(structsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& alesolver = list->sublist("ALE SOLVER",false,"");
  SetValidSolverParameters(alesolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& thermalsolver = list->sublist("THERMAL SOLVER",false,"");
  SetValidSolverParameters(thermalsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatrasolver = list->sublist("SCALAR TRANSPORT SOLVER",false,"solver parameters for scalar transport problems");
  SetValidSolverParameters(scatrasolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatrapotsolver = list->sublist("SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER",false,"solver parameters for block-preconditioning");
  SetValidSolverParameters(scatrapotsolver);

  return list;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::SetValidSolverParameters(Teuchos::ParameterList& list)
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  setStringToIntegralParameter<INPAR::SOLVER::SolverType>(
    "SOLVER", "UMFPACK",
    "The solver to attack the system of linear equations arising of FE approach with.",
    tuple<std::string>("Amesos_KLU_sym",
                       "Amesos_KLU_nonsym",
                       "Superlu",
                       "vm3",
                       "Aztec_MSR",
                       "LAPACK_sym",
                       "LAPACK_nonsym",
                       "UMFPACK"),
    tuple<INPAR::SOLVER::SolverType>(INPAR::SOLVER::amesos_klu_sym,
                                     INPAR::SOLVER::amesos_klu_nonsym,
                                     INPAR::SOLVER::superlu,
                                     INPAR::SOLVER::vm3,
                                     INPAR::SOLVER::aztec_msr,
                                     INPAR::SOLVER::lapack_sym,
                                     INPAR::SOLVER::lapack_nonsym,
                                     INPAR::SOLVER::umfpack),
    &list
    );

  setStringToIntegralParameter<INPAR::SOLVER::AzSolverType>(
    "AZSOLVE", "GMRES",
    "Type of linear solver algorithm to use.",
    tuple<std::string>("CG",
                       "GMRES",
                       "CGS",
                       "TFQMR",
                       "BiCGSTAB",
                       "LU"),
    tuple<INPAR::SOLVER::AzSolverType>(INPAR::SOLVER::azsolv_CG,
                                       INPAR::SOLVER::azsolv_GMRES,
                                       INPAR::SOLVER::azsolv_CGS,
                                       INPAR::SOLVER::azsolv_TFQMR,
                                       INPAR::SOLVER::azsolv_BiCGSTAB,
                                       INPAR::SOLVER::azsolv_LU),
    &list
    );

  {
    // this one is longer than 15 and the tuple<> function does not support this,
    // so build the Tuple class directly (which can be any size)
    Teuchos::Tuple<std::string,17> name;
    Teuchos::Tuple<INPAR::SOLVER::AzPrecType,17>  number;

    name[0] = "none";                         number[0] = INPAR::SOLVER::azprec_none;
    name[1] = "ILU";                          number[1] = INPAR::SOLVER::azprec_ILU;
    name[2] = "ILUT";                         number[2] = INPAR::SOLVER::azprec_ILUT;
    name[3] = "Jacobi";                       number[3] = INPAR::SOLVER::azprec_Jacobi;
    name[4] = "SymmGaussSeidel";              number[4] = INPAR::SOLVER::azprec_SymmGaussSeidel;
    name[5] = "Least_Squares";                number[5] = INPAR::SOLVER::azprec_Least_Squares;
    name[6] = "Neumann";                      number[6] = INPAR::SOLVER::azprec_Neumann;
    name[7] = "ICC";                          number[7] = INPAR::SOLVER::azprec_ICC;
    name[8] = "LU";                           number[8] = INPAR::SOLVER::azprec_LU;
    name[9] = "RILU";                         number[9] = INPAR::SOLVER::azprec_RILU;
    name[10] = "BILU";                        number[10] = INPAR::SOLVER::azprec_BILU;
    name[11] = "ML";                          number[11] = INPAR::SOLVER::azprec_ML;
    name[12] = "MLFLUID";                     number[12] = INPAR::SOLVER::azprec_MLfluid;
    name[13] = "MLFLUID2";                    number[13] = INPAR::SOLVER::azprec_MLfluid2;
    name[14] = "MLAPI";                       number[14] = INPAR::SOLVER::azprec_MLAPI;
    name[15] = "GaussSeidel";                 number[15] = INPAR::SOLVER::azprec_GaussSeidel;
    name[16] = "DownwindGaussSeidel";         number[16] = INPAR::SOLVER::azprec_DownwindGaussSeidel;

    setStringToIntegralParameter<INPAR::SOLVER::AzPrecType>(
      "AZPREC", "ILU",
      "Type of internal preconditioner to use.\n"
      "Note! this preconditioner will only be used if the input operator\n"
      "supports the Epetra_RowMatrix interface and the client does not pass\n"
      "in an external preconditioner!",
      name,
      number,
      &list
      );
  }

  IntParameter(
    "AZOVERLAP", 0,
    "The amount of overlap used for the internal \"ilu\" and \"ilut\" preconditioners.",
    &list
    );
  IntParameter(
    "AZGFILL", 0,
    "The amount of fill allowed for the internal \"ilu\" preconditioner.",
    &list
    );
  DoubleParameter(
    "AZDROP", 0.0,
    "The tolerance below which an entry from the factors of an internal \"ilut\"\n"
    "preconditioner will be dropped.",
    &list
    );
  DoubleParameter(
    "AZFILL", 1.0,
    "The amount of fill allowed for an internal \"ilut\" preconditioner.",
    &list
    );
//   IntParameter(
//     Steps_name, 3,
//     "Number of steps taken for the \"Jacobi\" or the \"Symmetric Gauss-Seidel\"\n"
//     "internal preconditioners for each preconditioner application.",
//     &list
//     );
  IntParameter(
    "AZPOLY", 3,
    "The order for of the polynomials used for the \"Polynomial\" and\n"
    "\"Least-squares Polynomial\" internal preconditioners.",
    &list
    );
//   setStringToIntegralParameter(
//     RCMReordering_name, "Disabled",
//     "Determines if RCM reordering is used with the internal\n"
//     "\"ilu\" or \"ilut\" preconditioners.",
//     tuple<std::string>("Enabled","Disabled"),
//     tuple<int>(1,0),
//     &list
//     );
//   setStringToIntegralParameter(
//     Orthogonalization_name, "Classical",
//     "The type of orthogonalization to use with the \"GMRES\" solver.",
//     tuple<std::string>("Classical","Modified"),
//     tuple<int>(AZ_classic,AZ_modified),
//     &list
//     );
  IntParameter(
    "AZSUB", 300,
    "The maximum size of the Krylov subspace used with \"GMRES\" before\n"
    "a restart is performed.",
    &list
    );
  setStringToIntegralParameter<int>(
    "AZCONV", "AZ_r0", // Same as "rhs" when x=0
    "The convergence test to use for terminating the iterative solver.",
    tuple<std::string>(
      "AZ_r0",
      "AZ_rhs",
      "AZ_Anorm",
      "AZ_noscaled",
      "AZ_sol",
      "AZ_weighted",
      "AZ_expected_values",
      "AZTECOO_conv_test",
      "AZ_inf_noscaled"
      ),
    tuple<int>(
      AZ_r0,
      AZ_rhs,
      AZ_Anorm,
      AZ_noscaled,
      AZ_sol,
      AZ_weighted,
      AZ_expected_values,
      AZTECOO_conv_test,
      AZ_inf_noscaled
      ),
    &list
    );
//   DoubleParameter(
//     IllConditioningThreshold_name, 1e+11,
//     "The threshold tolerance above which a system is considered\n"
//     "ill conditioned.",
//     &list
//     );
  IntParameter(
    "AZOUTPUT", 0, // By default, no output from Aztec!
    "The number of iterations between each output of the solver's progress.",
    &list
    );

  IntParameter("AZREUSE", 0, "how often to recompute some preconditioners", &list);
  IntParameter("AZITER", 1000, "max iterations", &list);
  IntParameter("AZGRAPH", 0, "unused", &list);
  IntParameter("AZBDIAG", 0, "", &list);

  DoubleParameter("AZTOL", 1e-8, "tolerance in (un)scaled residual", &list);
  DoubleParameter("AZOMEGA", 0.0, "damping for GaussSeidel and jacobi type methods", &list);
  DoubleParameter("DWINDTAU",1.5,"threshold tau for downwinding", &list);

  setStringToIntegralParameter<int>(
    "AZSCAL","none","scaling of the system",
    tuple<std::string>("none","sym","infnorm"),
    tuple<int>(0,1,2),
    &list);

  // parameters of ML preconditioner

  IntParameter("ML_PRINT",0,
               "ML print-out level (0-10)",&list);
  IntParameter("ML_MAXCOARSESIZE",5000,
               "ML stop coarsening when coarse ndof smaller then this",&list);
  IntParameter("ML_MAXLEVEL",5,
               "ML max number of levels",&list);
  IntParameter("ML_AGG_SIZE",27,
               "objective size of an aggregate with METIS/VBMETIS, 2D: 9, 3D: 27",&list);

  DoubleParameter("ML_DAMPFINE",1.,"damping fine grid",&list);
  DoubleParameter("ML_DAMPMED",1.,"damping med grids",&list);
  DoubleParameter("ML_DAMPCOARSE",1.,"damping coarse grid",&list);
  DoubleParameter("ML_PROLONG_SMO",0.,"damping factor for prolongator smoother (usually 1.33 or 0.0)",&list);
  DoubleParameter("ML_PROLONG_THRES",0.,"threshold for prolongator smoother/aggregation",&list);

  setNumericStringParameter("ML_SMOTIMES","1 1 1 1 1 1",
                            "no. smoothing steps or polynomial order on each level (at least ML_MAXLEVEL numbers)",&list);

  setStringToIntegralParameter<int>(
    "ML_COARSEN","UC","",
    tuple<std::string>("UC","METIS","VBMETIS","MIS"),
    tuple<int>(0,1,2,3),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERFINE","ILU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERMED","ILU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERCOARSE","KLU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &list);

  // unused
  setStringToIntegralParameter<int>("PARTITION","Cut_Elements","unused",
                               tuple<std::string>("Cut_Elements"),
                               tuple<int>(0),
                               &list);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list)
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  setStringToIntegralParameter<INPAR::STR::TimAdaKind>(
    "KIND","None","Method for time step size adapivity",
    tuple<std::string>(
      "None",
      "ZienkiewiczXie",
      "AdamsBashforth2"),
    tuple<INPAR::STR::TimAdaKind>(
      INPAR::STR::timada_kind_none,
      INPAR::STR::timada_kind_zienxie,
      INPAR::STR::timada_kind_ab2),
    &list);

  DoubleParameter("OUTSYSPERIOD", 0.0, "Write system vectors (displacements, velocities, etc) every given period of time", &list);
  DoubleParameter("OUTSTRPERIOD", 0.0, "Write stress/strain every given period of time", &list);
  DoubleParameter("OUTENEPERIOD", 0.0, "Write energy every given period of time", &list);
  DoubleParameter("OUTRESTPERIOD", 0.0, "Write restart data every given period of time", &list);
  IntParameter("OUTSIZEEVERY", 0, "Write step size every given time step", &list);

  DoubleParameter("STEPSIZEMAX", 0.0, "Limit maximally permitted time step size (>0)", &list);
  DoubleParameter("STEPSIZEMIN", 0.0, "Limit minimally allowed time step size (>0)", &list);
  DoubleParameter("SIZERATIOMAX", 0.0, "Limit maximally permitted change of time step size compared to previous size, important for multi-step schemes (>0)", &list);
  DoubleParameter("SIZERATIOMIN", 0.0, "Limit minimally permitted change of time step size compared to previous size, important for multi-step schemes (>0)", &list);
  DoubleParameter("SIZERATIOSCALE", 0.9, "This is a safety factor to scale theretical optimal step size, should be lower than 1 and must be larger than 0", &list);

  setStringToIntegralParameter<INPAR::STR::VectorNorm>(
    "LOCERRNORM", "Vague", "Vector norm to treat error vector with",
    tuple<std::string>(
      "Vague",
      "L1",
      "L2",
      "Rms",
      "Inf"),
    tuple<INPAR::STR::VectorNorm>(
      INPAR::STR::norm_vague,
      INPAR::STR::norm_l1,
      INPAR::STR::norm_l2,
      INPAR::STR::norm_rms,
      INPAR::STR::norm_inf),
    &list);

  DoubleParameter("LOCERRTOL", 0.0, "Target local error tolerance (>0)", &list);
  IntParameter("ADAPTSTEPMAX", 0, "Limit maximally allowed step size reduction attempts (>0)", &list);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::SetValidNoxParameters(Teuchos::ParameterList& list)
{
  {
    Teuchos::Array<std::string> st = Teuchos::tuple<std::string>(
      "Line Search Based",
      "Trust Region Based",
      "Inexact Trust Region Based",
      "Tensor Based");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Nonlinear Solver","Line Search Based","",
      st,st,
      &list);  
  }

  // sub-list direction
  Teuchos::ParameterList& direction = list.sublist("Direction",false,"");

  {
    Teuchos::Array<std::string> st = Teuchos::tuple<std::string>(
      "Newton",
      "Steepest Descent",
      "NonlinearCG",
      "Broyden");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Method","Newton","",
      st,st,
      &direction);
  }

  // sub-sub-list "Newton"
  Teuchos::ParameterList& newton = direction.sublist("Newton",false,"");

  {
    Teuchos::Array<std::string> forcingtermmethod = Teuchos::tuple<std::string>(
      "Constant",
      "Type 1",
      "Type 2");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Forcing Term Method","Constant","",
      forcingtermmethod,forcingtermmethod,
      &newton);
    DoubleParameter("Forcing Term Initial Tolerance",0.1,"initial linear solver tolerance",&newton);
    DoubleParameter("Forcing Term Minimum Tolerance",1.0e-6,"",&newton);
    DoubleParameter("Forcing Term Maximum Tolerance",0.01,"",&newton);
    DoubleParameter("Forcing Term Alpha",1.5,"used only by \"Type 2\"",&newton);
    DoubleParameter("Forcing Term Gamma",0.9,"used only by \"Type 2\"",&newton);
    BoolParameter("Rescue Bad Newton Solver","Yes","If set to true, we will use the computed direction even if the linear solve does not achieve the tolerance specified by the forcing term",&newton);
  }

  // sub-sub-list "Steepest Descent"
  Teuchos::ParameterList& steepestdescent = direction.sublist("Steepest Descent",false,"");

  {
    Teuchos::Array<std::string> scalingtype = Teuchos::tuple<std::string>(
      "2-Norm",
      "Quadratic Model Min",
      "F 2-Norm",
      "None");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Scaling Type","None","",
      scalingtype,scalingtype,
      &steepestdescent);
  }

  // sub-list "Line Search"
  Teuchos::ParameterList& linesearch = list.sublist("Line Search",false,"");

  {
    Teuchos::Array<std::string> method = Teuchos::tuple<std::string>(
      "Full Step",
      "Backtrack" ,
      "Polynomial",
      "More'-Thuente",
      "User Defined");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Method","Full Step","",
      method,method,
      &linesearch);
  }

  // sub-sub-list "Full Step"
  Teuchos::ParameterList& fullstep = linesearch.sublist("Full Step",false,"");

  {
    DoubleParameter("Full Step",1.0,"length of a full step",&fullstep);
  }

  // sub-sub-list "Backtrack"
  Teuchos::ParameterList& backtrack = linesearch.sublist("Backtrack",false,"");

  {
    DoubleParameter("Default Step",1.0,"starting step length",&backtrack);
    DoubleParameter("Minimum Step",1.0e-12,"minimum acceptable step length",&backtrack);
    DoubleParameter("Recovery Step",1.0,"step to take when the line search fails (defaults to value for \"Default Step\")",&backtrack);
    IntParameter("Max Iters",50,"maximum number of iterations (i.e., RHS computations)",&backtrack);
    DoubleParameter("Reduction Factor",0.5,"A multiplier between zero and one that reduces the step size between line search iterations",&backtrack);
  }

  // sub-sub-list "Polynomial"
  Teuchos::ParameterList& polynomial = linesearch.sublist("Polynomial",false,"");

  {
    DoubleParameter("Default Step",1.0,"Starting step length",&polynomial);
    IntParameter("Max Iters",100,"Maximum number of line search iterations. The search fails if the number of iterations exceeds this value",&polynomial);
    DoubleParameter("Minimum Step",1.0e-12,"Minimum acceptable step length. The search fails if the computed $lambda_k$ is less than this value",&polynomial);
    Teuchos::Array<std::string> recoverysteptype = Teuchos::tuple<std::string>(
      "Constant",
      "Last Computed Step");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Recovery Step Type","Constant","Determines the step size to take when the line search fails",
      recoverysteptype,recoverysteptype,
      &polynomial);
    DoubleParameter("Recovery Step",1.0,"The value of the step to take when the line search fails. Only used if the \"Recovery Step Type\" is set to \"Constant\"",&polynomial);
    Teuchos::Array<std::string> interpolationtype = Teuchos::tuple<std::string>(
      "Quadratic",
      "Quadratic3",
      "Cubic");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Interpolation Type","Cubic","Type of interpolation that should be used",
      interpolationtype,interpolationtype,
      &polynomial);
    DoubleParameter("Min Bounds Factor",0.1,"Choice for $gamma_{min}$, i.e., the factor that limits the minimum size of the new step based on the previous step",&polynomial);
    DoubleParameter("Max Bounds Factor",0.5,"Choice for $gamma_{max}$, i.e., the factor that limits the maximum size of the new step based on the previous step",&polynomial);
    Teuchos::Array<std::string> sufficientdecreasecondition = Teuchos::tuple<std::string>(
      "Armijo-Goldstein",
      "Ared/Pred",
      "None");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Sufficient Decrease Condition","Armijo-Goldstein","Choice to use for the sufficient decrease condition",
      sufficientdecreasecondition,sufficientdecreasecondition,
      &polynomial);
    DoubleParameter("Alpha Factor",1.0e-4,"Parameter choice for sufficient decrease condition",&polynomial);
    BoolParameter("Force Interpolation","No","Set to true if at least one interpolation step should be used. The default is false which means that the line search will stop if the default step length satisfies the convergence criteria",&polynomial);
    BoolParameter("Use Counters","Yes","Set to true if we should use counters and then output the result to the paramter list as described in Output Parameters",&polynomial);
    IntParameter("Maximum Iteration for Increase",0,"Maximum index of the nonlinear iteration for which we allow a relative increase",&polynomial);
    DoubleParameter("Allowed Relative Increase",100,"",&polynomial);
  }

  // sub-sub-list "More'-Thuente"
  Teuchos::ParameterList& morethuente = linesearch.sublist("More'-Thuente",false,"");

  {
    DoubleParameter("Sufficient Decrease",1.0e-4,"The ftol in the sufficient decrease condition",&morethuente);
    DoubleParameter("Curvature Condition",0.9999,"The gtol in the curvature condition",&morethuente);
    DoubleParameter("Interval Width",1.0e-15,"The maximum width of the interval containing the minimum of the modified function",&morethuente);
    DoubleParameter("Maximum Step",1.0e6,"maximum allowable step length",&morethuente);
    DoubleParameter("Minimum Step",1.0e-12,"minimum allowable step length",&morethuente);
    IntParameter("Max Iters",20,"maximum number of right-hand-side and corresponding Jacobian evaluations",&morethuente);
    DoubleParameter("Default Step",1.0,"starting step length",&morethuente);
    Teuchos::Array<std::string> recoverysteptype = Teuchos::tuple<std::string>(
      "Constant",
      "Last Computed Step");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Recovery Step Type","Constant","Determines the step size to take when the line search fails",
      recoverysteptype,recoverysteptype,
      &morethuente);
    DoubleParameter("Recovery Step",1.0,"The value of the step to take when the line search fails. Only used if the \"Recovery Step Type\" is set to \"Constant\"",&morethuente);
    Teuchos::Array<std::string> sufficientdecreasecondition = Teuchos::tuple<std::string>(
      "Armijo-Goldstein",
      "Ared/Pred",
      "None");
    Teuchos::setStringToIntegralParameter<std::string>(
      "Sufficient Decrease Condition","Armijo-Goldstein","Choice to use for the sufficient decrease condition",
      sufficientdecreasecondition,sufficientdecreasecondition,
      &morethuente);    
    BoolParameter("Optimize Slope Calculation","No","Boolean value. If set to true the value of $s^T J^T F$ is estimated using a directional derivative in a call to NOX::LineSearch::Common::computeSlopeWithOutJac. If false the slope computation is computed with the NOX::LineSearch::Common::computeSlope method. Setting this to true eliminates having to compute the Jacobian at each inner iteration of the More'-Thuente line search",&morethuente);
  }

  // sub-list "Trust Region"
  Teuchos::ParameterList& trustregion = list.sublist("Trust Region",false,"");

  {
    DoubleParameter("Minimum Trust Region Radius",1.0e-6,"Minimum allowable trust region radius",&trustregion);
    DoubleParameter("Maximum Trust Region Radius",1.0e+9,"Maximum allowable trust region radius",&trustregion);
    DoubleParameter("Minimum Improvement Ratio",1.0e-4,"Minimum improvement ratio to accept the step",&trustregion);
    DoubleParameter("Contraction Trigger Ratio",0.1,"If the improvement ratio is less than this value, then the trust region is contracted by the amount specified by the \"Contraction Factor\". Must be larger than \"Minimum Improvement Ratio\"",&trustregion);
    DoubleParameter("Contraction Factor",0.25,"",&trustregion);
    DoubleParameter("Expansion Trigger Ratio",0.75,"If the improvement ratio is greater than this value, then the trust region is contracted by the amount specified by the \"Expansion Factor\"",&trustregion);
    DoubleParameter("Expansion Factor",4.0,"",&trustregion);
    DoubleParameter("Recovery Step",1.0,"",&trustregion);
  }

  // sub-list "Printing"
  Teuchos::ParameterList& printing = list.sublist("Printing",false,"");

  {
    BoolParameter("Error","No","",&printing);
    BoolParameter("Warning","Yes","",&printing);
    BoolParameter("Outer Iteration","Yes","",&printing);
    BoolParameter("Inner Iteration","Yes","",&printing);
    BoolParameter("Parameters","No","",&printing);
    BoolParameter("Details","No","",&printing);
    BoolParameter("Outer Iteration StatusTest","No","",&printing);
    BoolParameter("Linear Solver Details","No","",&printing);
    BoolParameter("Test Details","No","",&printing);
    /*  // for LOCA
    BoolParameter("Stepper Iteration","No","",&printing);
    BoolParameter("Stepper Details","No","",&printing);
    BoolParameter("Stepper Parameters","Yes","",&printing);
    */
    BoolParameter("Debug","No","",&printing);
  }
}

#endif
