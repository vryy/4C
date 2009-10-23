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
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_potential.H"
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_inpar/inpar_elch.H"

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
void DRT::INPUT::PrintDatHeader(std::ostream& stream,
                                const Teuchos::ParameterList& list,
                                std::string parentname,
                                bool color,
                                bool comment)
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
      if (name == PrintEqualSign()) continue;
      Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();

      if (comment)
      {
        stream << blue2light << "//" << endcolor << '\n';

        std::string doc = entry.docString();
        if (doc!="")
        {
          Teuchos::StrUtils::printLines(stream,blue2light + "// ",doc);
          stream << endcolor;
        }
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
        PrintDatHeader(stream,list.sublist(name),secname,color,comment);
      }
      else
      {
        if (comment)
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
        if (NeedToPrintEqualSign(list)) stream << " =";
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

  {
  Teuchos::Tuple<std::string,16> name;
  Teuchos::Tuple<PROBLEM_TYP,16> label;
  name[ 0] = "Structure";                                   label[ 0] = prb_structure;
  name[ 1] = "Fluid";                                       label[ 1] = prb_fluid;
  name[ 2] = "Fluid_XFEM";                                  label[ 2] = prb_fluid_xfem;
  name[ 3] = "Fluid_Ale";                                   label[ 3] = prb_fluid_ale;
  name[ 4] = "Fluid_Freesurface";                           label[ 4] = prb_freesurf;
  name[ 5] = "Scalar_Transport";                            label[ 5] = prb_scatra;
  name[ 6] = "Fluid_Structure_Interaction";                 label[ 6] = prb_fsi;
  name[ 7] = "Fluid_Structure_Interaction_XFEM";            label[ 7] = prb_fsi_xfem;
  name[ 8] = "Ale";                                         label[ 8] = prb_ale;
  name[ 9] = "Thermal_Structure_Interaction";               label[ 9] = prb_tsi;
  name[10] = "Thermo";                                      label[10] = prb_thermo;
  name[11] = "Structure_Multiscale";                        label[11] = prb_struct_multi;
  name[12] = "Low_Mach_Number_Flow";                        label[12] = prb_loma;
  name[13] = "Electrochemistry";                            label[13] = prb_elch;
  name[14] = "Combustion";                                  label[14] = prb_combust;
  name[15] = "ArterialNetwork";                             label[15] = prb_art_net;
  setStringToIntegralParameter<PROBLEM_TYP>(
                               "PROBLEMTYP",
                               "Fluid_Structure_Interaction",
                               "",
                               name,
                               label,
                               &type);
  }

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
  setStringToIntegralParameter<INPAR::THR::HeatFluxType>("THERM_HEATFLUX","None","",
                               tuple<std::string>("None",
                                                  "No",
                                                  "NO",
                                                  "no",
                                                  "Current",
                                                  "Initial"),
                               tuple<INPAR::THR::HeatFluxType>(INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_current,
                                                               INPAR::THR::heatflux_initial),
                               &io);
  setStringToIntegralParameter<INPAR::THR::TempGradType>("THERM_TEMPGRAD","None","",
                               tuple<std::string>("None",
                                                  "No",
                                                  "NO",
                                                  "no",
                                                  "Current",
                                                  "Initial"),
                               tuple<INPAR::THR::TempGradType>(INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_current,
                                                               INPAR::THR::tempgrad_initial),
                               &io);

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
                                 "AdamsBashforth2",
                                 "EulerMaruyama",
                                 "EulerImpStoch"),
                               tuple<INPAR::STR::DynamicType>(
                                 INPAR::STR::dyna_centr_diff,
                                 INPAR::STR::dyna_Gen_EMM,
                                 INPAR::STR::dyna_gen_alfa,
                                 INPAR::STR::dyna_gen_alfa_statics,
                                 INPAR::STR::dyna_statics,
                                 INPAR::STR::dyna_genalpha,
                                 INPAR::STR::dyna_onesteptheta,
                                 INPAR::STR::dyna_gemm,
                                 INPAR::STR::dyna_ab2,
                                 INPAR::STR::dyna_euma,
                                 INPAR::STR::dyna_euimsto),
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
                                 "Material",
                                 "BrownianMotion"),
                               tuple<INPAR::STR::DampKind>(
                                 INPAR::STR::damp_none,
                                 INPAR::STR::damp_none,
                                 INPAR::STR::damp_none,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_rayleigh,
                                 INPAR::STR::damp_material,
                                 INPAR::STR::damp_brownianmotion),
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
                                 "RelRes_And_RelDis",
                                 "MixRes_Or_MixDis",
                                 "MixRes_And_MixDis",
                                 "None"),
                               tuple<INPAR::STR::ConvCheck>(
                                 INPAR::STR::convcheck_absres_or_absdis,
                                 INPAR::STR::convcheck_absres_and_absdis,
                                 INPAR::STR::convcheck_relres_or_absdis,
                                 INPAR::STR::convcheck_relres_and_absdis,
                                 INPAR::STR::convcheck_relres_or_reldis,
                                 INPAR::STR::convcheck_relres_and_reldis,
                                 INPAR::STR::convcheck_mixres_or_mixdis,
                                 INPAR::STR::convcheck_mixres_and_mixdis,
                                 INPAR::STR::convcheck_vague),
                               &sdyn);

  DoubleParameter("TOLDISP",1.0E-10,
                  "tolerance in the displacement norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<INPAR::STR::ConvNorm>("NORM_DISP","Abs","type of norm for displacement convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<INPAR::STR::ConvNorm>(
                                 INPAR::STR::convnorm_abs,
                                 INPAR::STR::convnorm_rel,
                                 INPAR::STR::convnorm_mix),
                               &sdyn);

  DoubleParameter("TOLRES",1.0E-08,
                  "tolerance in the residual norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<INPAR::STR::ConvNorm>("NORM_RESF","Abs","type of norm for residual convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<INPAR::STR::ConvNorm>(
                                 INPAR::STR::convnorm_abs,
                                 INPAR::STR::convnorm_rel,
                                 INPAR::STR::convnorm_mix),
                               &sdyn);

  DoubleParameter("TOLPRE",1.0E-08,
                  "tolerance in pressure norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<INPAR::STR::ConvNorm>("NORM_PRES","Abs","type of norm for pressure convergence check",
                               tuple<std::string>(
                                 "Abs"),
                               tuple<INPAR::STR::ConvNorm>(
                                 INPAR::STR::convnorm_abs),
                               &sdyn);

  DoubleParameter("TOLINCO",1.0E-08,
                  "tolerance in the incompressible residual norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<INPAR::STR::ConvNorm>("NORM_INCO","Abs","type of norm for incompressible residual convergence check",
                               tuple<std::string>(
                                 "Abs"),
                               tuple<INPAR::STR::ConvNorm>(
                                 INPAR::STR::convnorm_abs),
                               &sdyn);

  setStringToIntegralParameter<INPAR::STR::BinaryOp>("NORMCOMBI_DISPPRES","And","binary operator to combine pressure and displacement values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<INPAR::STR::BinaryOp>(
                                 INPAR::STR::bop_and,
                                 INPAR::STR::bop_or),
                               &sdyn);

  setStringToIntegralParameter<INPAR::STR::BinaryOp>("NORMCOMBI_RESFINCO","And","binary operator to combine force and incompressible residual",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<INPAR::STR::BinaryOp>(
                                 INPAR::STR::bop_and,
                                 INPAR::STR::bop_or),
                               &sdyn);

  setStringToIntegralParameter<INPAR::STR::BinaryOp>("NORMCOMBI_RESFDISP","And","binary operator to combine displacement and residual force values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<INPAR::STR::BinaryOp>(
                                 INPAR::STR::bop_and,
                                 INPAR::STR::bop_or),
                               &sdyn);

  DoubleParameter("TOLCONSTR",1.0E-08,
                  "tolerance in the constr error norm for the newton iteration",
                  &sdyn);
  IntParameter("MAXITER",50,
               "maximum number of iterations allowed for Newton-Raphson iteration before failure",
               &sdyn);
  IntParameter("MINITER",0,
               "minimum number of iterations to be done within Newton-Raphson loop",
               &sdyn);
  setStringToIntegralParameter<INPAR::STR::VectorNorm>("ITERNORM","L2","type of norm to be applied to residuals",
                               tuple<std::string>(
                                 "L1",
                                 "L2",
                                 "Rms",
                                 "Inf"),
                               tuple<INPAR::STR::VectorNorm>(
                                 INPAR::STR::norm_l1,
                                 INPAR::STR::norm_l2,
                                 INPAR::STR::norm_rms,
                                 INPAR::STR::norm_inf),
                               &sdyn);
  setStringToIntegralParameter<int>("DIVERCONT","No",
                                    "Go on with time integration even if Newton-Raphson iteration failed",
                                    yesnotuple,yesnovalue,&sdyn);

  setStringToIntegralParameter<INPAR::STR::NonlinSolTech>("NLNSOL","fullnewton","Nonlinear solution technique",
                               tuple<std::string>(
                                 "vague",
                                 "fullnewton",
                                 "lsnewton",
                                 "oppnewton",
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
                                 INPAR::STR::soltech_newtonopp,
                                 INPAR::STR::soltech_newtonmod,
                                 INPAR::STR::soltech_nlncg,
                                 INPAR::STR::soltech_ptc,
                                 INPAR::STR::soltech_newtonuzawalin,
                                 INPAR::STR::soltech_newtonuzawanonlin,
                                 INPAR::STR::soltech_noxnewtonlinesearch,
                                 INPAR::STR::soltech_noxgeneral),
                               &sdyn);

  setStringToIntegralParameter<INPAR::STR::ControlType>("CONTROLTYPE","load","load, disp, arc1, arc2 control",
                               tuple<std::string>(
                                 "load",
                                 "Load",
                                 "disp",
                                 "Disp",
                                 "Displacement",
                                 "arc1",
                                 "Arc1",
                                 "arc2",
                                 "Arc2"),
                               tuple<INPAR::STR::ControlType>(
                                 INPAR::STR::control_load,
                                 INPAR::STR::control_load,
                                 INPAR::STR::control_disp,
                                 INPAR::STR::control_disp,
                                 INPAR::STR::control_disp,
                                 INPAR::STR::control_arc1,
                                 INPAR::STR::control_arc1,
                                 INPAR::STR::control_arc2,
                                 INPAR::STR::control_arc2),
                               &sdyn);

  setNumericStringParameter("CONTROLNODE","-1 -1 -1",
                            "for methods other than load control: [node(fortran numbering)] [dof(c-numbering)] [curve(fortran numbering)]",
                            &sdyn);

  setStringToIntegralParameter<int>("LOADLIN","no",
                                    "Use linearization of external follower load in Newton",
                                    yesnotuple,yesnovalue,&sdyn);

  setStringToIntegralParameter<INPAR::STR::PredEnum>("PREDICT","ConstDis","",
                               tuple<std::string>(
                                 "Vague",
                                 "ConstDis",
                                 "ConstDisVelAcc",
                                 "TangDis",
                                 "ConstDisPres",
                                 "ConstDisVelAccPres"),
                               tuple<INPAR::STR::PredEnum>(
                                 INPAR::STR::pred_vague,
                                 INPAR::STR::pred_constdis,
                                 INPAR::STR::pred_constdisvelacc,
                                 INPAR::STR::pred_tangdis,
                                 INPAR::STR::pred_constdispres,
                                 INPAR::STR::pred_constdisvelaccpres),
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
  DoubleParameter("MEASURED_CURVE0",0.0,"measured displacement of the tension testing",&iap);
  DoubleParameter("MEASURED_CURVE1",0.0,"measured displacement of the tension testing",&iap);
  DoubleParameter("MEASURED_CURVE2",0.0,"measured displacement of the tension testing",&iap);

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

  setStringToIntegralParameter<INPAR::CONTACT::SolvingStrategy>("STRATEGY","LagrangianMultipliers","Type of employed solving strategy",
        tuple<std::string>("LagrangianMultipliers","lagrange", "Lagrange",
            "PenaltyMethod","penalty", "Penalty",
            "AugmentedLagrange","augmented", "Augmented"),
            tuple<INPAR::CONTACT::SolvingStrategy>(INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult,
                INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty,
                INPAR::CONTACT::solution_auglag, INPAR::CONTACT::solution_auglag, INPAR::CONTACT::solution_auglag),
                &scontact);

  setStringToIntegralParameter<INPAR::CONTACT::ShapeFcn>("SHAPEFCN","Dual","Type of employed set of shape functions",
        tuple<std::string>("Dual", "dual",
            "Standard", "standard", "std"),
            tuple<INPAR::CONTACT::ShapeFcn>(INPAR::CONTACT::shape_dual, INPAR::CONTACT::shape_dual,
                INPAR::CONTACT::shape_standard, INPAR::CONTACT::shape_standard, INPAR::CONTACT::shape_standard),
                &scontact);

  DoubleParameter("PENALTYPARAM",0.0,"Penalty parameter for penalty / augmented solution strategy",&scontact);
  IntParameter("UZAWAMAXSTEPS",10,"Maximum no. of Uzawa steps for augmented / Uzawa solution strategy",&scontact);
  DoubleParameter("UZAWACONSTRTOL",1.0e-8,"Tolerance of constraint norm for augmented / Uzawa solution strategy",&scontact);
  DoubleParameter("FRBOUND",0.0,"Friction bound for Tresca friction",&scontact);
  DoubleParameter("FRCOEFF",0.0,"Friction coefficient for Coulomb friction",&scontact);

  setStringToIntegralParameter<int>("FULL_LINEARIZATION","Yes","If chosen full linearization of contact is applied",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON","Yes","If chosen semi-smooth Newton concept is applied",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("SEMI_SMOOTH_CN",1.0,"Weighting factor cn for semi-smooth PDASS",&scontact);
  DoubleParameter("SEMI_SMOOTH_CT",1.0,"Weighting factor ct for semi-smooth PDASS",&scontact);

  setStringToIntegralParameter<INPAR::CONTACT::ContactSearchAlgorithm>("SEARCH_ALGORITHM","Binarytree","Type of contact search",
                               tuple<std::string>("BruteForceNodeBased","bruteforcenodebased",
                                                  "BruteForceEleBased","bruteforceelebased",
                                                  "BinaryTree","Binarytree","binarytree"),
                               tuple<INPAR::CONTACT::ContactSearchAlgorithm>(INPAR::CONTACT::search_bfnode,INPAR::CONTACT::search_bfnode,
                                          INPAR::CONTACT::search_bfele,INPAR::CONTACT::search_bfele,
                                          INPAR::CONTACT::search_binarytree,INPAR::CONTACT::search_binarytree,
                                          INPAR::CONTACT::search_binarytree),
                               &scontact);

  DoubleParameter("SEARCH_PARAM",0.3,"Radius / Bounding volume inflation for contact search",&scontact);

  setStringToIntegralParameter<int>("COUPLING_AUXPLANE","Yes","If chosen auxiliary planes are used for 3D coupling",
                               yesnotuple,yesnovalue,&scontact);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& interaction_potential = list->sublist("INTERACTION POTENTIAL",false,"");

  // read if surfaces , volumes or both including fluid should be considered
  setStringToIntegralParameter<INPAR::POTENTIAL::PotentialType>("POTENTIAL_TYPE","Surface","Type of interaction potential",
                                tuple<std::string>("Surface",
                                                   "Volume",
                                                   "Surfacevolume",
                                                   "Surface_fsi",
                                                   "Volume_fsi",
                                                   "Surfacevolume_fsi"),
                                tuple<INPAR::POTENTIAL::PotentialType>(
                                   INPAR::POTENTIAL::potential_surface,
                                   INPAR::POTENTIAL::potential_volume,
                                   INPAR::POTENTIAL::potential_surfacevolume,
                                   INPAR::POTENTIAL::potential_surface_fsi,
                                   INPAR::POTENTIAL::potential_volume_fsi,
                                   INPAR::POTENTIAL::potential_surfacevolume_fsi),
                                &interaction_potential);

  // approximation method
  setStringToIntegralParameter<INPAR::POTENTIAL::ApproximationType>("APPROXIMATION_TYPE","None","Type of approximation",
                                tuple<std::string>("None",
                                                   "Surface_approx",
                                                   "Point_approx"),
                                tuple<INPAR::POTENTIAL::ApproximationType>(
                                           INPAR::POTENTIAL::approximation_none,
                                           INPAR::POTENTIAL::approximation_surface,
                                           INPAR::POTENTIAL::approximation_point),
                                &interaction_potential);

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
                                                    "anisotropic","Anisotropic",
                                                    "endtoend_const",
                                                    "Viscoelasticity","viscoelasticity","ViscoElasticity",
                                                    "Gmsh","gmsh"),
                                 //translating input strings into BACI input parameters
                                 tuple<INPAR::STATMECH::StatOutput>(INPAR::STATMECH::statout_none,INPAR::STATMECH::statout_none,
                                            INPAR::STATMECH::statout_endtoendlog,INPAR::STATMECH::statout_endtoendlog,INPAR::STATMECH::statout_endtoendlog,
                                            INPAR::STATMECH::statout_anisotropic,INPAR::STATMECH::statout_anisotropic,
                                            INPAR::STATMECH::statout_endtoendconst,INPAR::STATMECH::statout_endtoendconst,
                                            INPAR::STATMECH::statout_viscoelasticity,INPAR::STATMECH::statout_viscoelasticity,INPAR::STATMECH::statout_viscoelasticity,
                                            INPAR::STATMECH::statout_gmsh,INPAR::STATMECH::statout_gmsh),
                                 &statmech);
  //Reading which kind of friction model should be applied
  setStringToIntegralParameter<INPAR::STATMECH::FrictionModel>("FRICTION_MODEL","none","friction model for polymer dynamics",
                                 //listing possible strings in input file in category FRICTION_MODEL
                                 tuple<std::string>("none",
                                                    "isotropiclumped",
                                                    "isotropicconsistent",
                                                    "anisotropicconsistent"),
                                 //translating input strings into BACI input parameters
                                 tuple<INPAR::STATMECH::FrictionModel>(INPAR::STATMECH::frictionmodel_none,
                                                                    INPAR::STATMECH::frictionmodel_isotropiclumped,
                                                                    INPAR::STATMECH::frictionmodel_isotropicconsistent,
                                                                    INPAR::STATMECH::frictionmodel_anisotropicconsistent),
                                                                    &statmech);
  //percentage of total simulation time after which writing of statistical output is started
  DoubleParameter("START_FACTOR",0.0,"Percentage of total simulation time after which writing of statistical output is started",&statmech);
  //Reading whether dynamics remodelling of cross linker distribution takes place
  setStringToIntegralParameter<int>("DYN_CROSSLINKERS","No","If chosen cross linker proteins are added and removed in each time step",
                               yesnotuple,yesnovalue,&statmech);
  //Reading double parameters for shear flow field
  DoubleParameter("SHEARAMPLITUDE",0.0,"Shear amplitude of flow in z-direction; note: not amplitude of displacement, but of shear strain!",&statmech);
  DoubleParameter("SHEARFREQUENCY",0.0,"Shear frequency of flow in z-direction",&statmech);
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

  /*----------------------------------------------------------------------*/
   Teuchos::ParameterList& tdyn = list->sublist("THERMAL DYNAMIC",false,"");

   setStringToIntegralParameter<INPAR::THR::DynamicType>("DYNAMICTYP","OneStepTheta",
                                "type of time integration control",
                                tuple<std::string>(
                                  "Statics",
                                  "OneStepTheta",
                                  "GEMM",
                                  "GenAlpha"
                                  ),
                                tuple<INPAR::THR::DynamicType>(
                                   INPAR::THR::dyna_statics,
                                   INPAR::THR::dyna_onesteptheta,
                                   INPAR::THR::dyna_gemm,
                                   INPAR::THR::dyna_genalpha),
                                &tdyn);

   // Output type
   IntParameter("RESEVRYGLOB",1,"save temperature and other global quantities every RESEVRYGLOB steps",&tdyn);
   IntParameter("RESEVRYELEM",1,"save heat fluxes and other element quantities every RESEVRYELEM steps",&tdyn);
   IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&tdyn);
   // Time loop control
   DoubleParameter("TIMESTEP",0.05,"time step size",&tdyn);
   IntParameter("NUMSTEP",200,"maximum number of steps",&tdyn);
   DoubleParameter("MAXTIME",5.0,"maximum time",&tdyn);
   // Iterationparameters
   DoubleParameter("TOLTEMP",1.0E-10,
                   "tolerance in the temperature norm of the Newton iteration",
                   &tdyn);
   setStringToIntegralParameter<INPAR::THR::ConvNorm>("NORM_TEMP","Abs","type of norm for temperature convergence check",
                                tuple<std::string>(
                                  "Abs",
                                  "Rel",
                                  "Mix"),
                                tuple<INPAR::THR::ConvNorm>(
                                  INPAR::THR::convnorm_abs,
                                  INPAR::THR::convnorm_rel,
                                  INPAR::THR::convnorm_mix),
                                &tdyn);

   DoubleParameter("TOLRES",1.0E-08,
                   "tolerance in the residual norm for the Newton iteration",
                   &tdyn);
   setStringToIntegralParameter<INPAR::THR::ConvNorm>("NORM_RESF","Abs","type of norm for residual convergence check",
                                tuple<std::string>(
                                  "Abs",
                                  "Rel",
                                  "Mix"),
                                tuple<INPAR::THR::ConvNorm>(
                                  INPAR::THR::convnorm_abs,
                                  INPAR::THR::convnorm_rel,
                                  INPAR::THR::convnorm_mix),
                                &tdyn);

   setStringToIntegralParameter<INPAR::THR::BinaryOp>("NORMCOMBI_RESFTEMP","And","binary operator to combine temperature and residual force values",
                                tuple<std::string>(
                                  "And",
                                  "Or"),
                                tuple<INPAR::THR::BinaryOp>(
                                  INPAR::THR::bop_and,
                                  INPAR::THR::bop_or),
                                &tdyn);


   IntParameter("MAXITER",50,
                "maximum number of iterations allowed for Newton-Raphson iteration before failure",
                &tdyn);
   IntParameter("MINITER",0,
                "minimum number of iterations to be done within Newton-Raphson loop",
                &tdyn);
   setStringToIntegralParameter<INPAR::THR::VectorNorm>("ITERNORM","L2","type of norm to be applied to residuals",
                                tuple<std::string>(
                                  "L1",
                                  "L2",
                                  "Rms",
                                  "Inf"),
                                tuple<INPAR::THR::VectorNorm>(
                                  INPAR::THR::norm_l1,
                                  INPAR::THR::norm_l2,
                                  INPAR::THR::norm_rms,
                                  INPAR::THR::norm_inf),
                                &tdyn);
   setStringToIntegralParameter<int>("DIVERCONT","No",
                                     "Go on with time integration even if Newton-Raphson iteration failed",
                                     yesnotuple,yesnovalue,&tdyn);

   setStringToIntegralParameter<INPAR::THR::NonlinSolTech>("NLNSOL","fullnewton","Nonlinear solution technique",
                                tuple<std::string>(
                                  "vague",
                                  "fullnewton"),
                                tuple<INPAR::THR::NonlinSolTech>(
                                  INPAR::THR::soltech_vague,
                                  INPAR::THR::soltech_newtonfull),
                                &tdyn);

   setStringToIntegralParameter<INPAR::THR::PredEnum>("PREDICT","ConstTemp","Predictor of iterative solution techniques",
                                tuple<std::string>(
                                  "Vague",
                                  "ConstTemp",
                                  "ConstTempRate"
                                  "TangTemp"),
                                tuple<INPAR::THR::PredEnum>(
                                  INPAR::THR::pred_vague,
                                  INPAR::THR::pred_consttemp,
                                  INPAR::THR::pred_consttemprate,
                                  INPAR::THR::pred_tangtemp),
                                &tdyn);

   // convergence criteria solver adaptivity
   setStringToIntegralParameter<int>("ADAPTCONV","No",
                                "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                                yesnotuple,yesnovalue,&tdyn);
   DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&tdyn);

   /*----------------------------------------------------------------------*/
   /* parameters for generalised-alpha thermal integrator */
   Teuchos::ParameterList& tgenalpha = tdyn.sublist("GENALPHA",false,"");

   setStringToIntegralParameter<INPAR::THR::MidAverageEnum>("GENAVG","ImrLike",
                                "mid-average type of internal forces",
                                tuple<std::string>(
                                  "Vague",
                                  "ImrLike",
                                  "TrLike"),
                                tuple<INPAR::THR::MidAverageEnum>(
                                  INPAR::THR::midavg_vague,
                                  INPAR::THR::midavg_imrlike,
                                  INPAR::THR::midavg_trlike),
                                &tgenalpha);
   DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,1/2]",&tgenalpha);
   DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&tgenalpha);
   DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&tgenalpha);
   DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&tgenalpha);

   /*----------------------------------------------------------------------*/
   /* parameters for one-step-theta thermal integrator */
   Teuchos::ParameterList& tonesteptheta = tdyn.sublist("ONESTEPTHETA",false,"");

   DoubleParameter("THETA",0.5,"One-step-theta factor in (0,1]",&tonesteptheta);

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
                               "Nonlinear Time Integration Scheme",
                               tuple<std::string>(
                                 "Nlin_Time_Int",
                                 "Lin_Time_Int"
                                 ),
                               tuple<int>(
                                dyntyp_nln_time_int,
                                dyntyp_lin_time_int
                                ),
                               &fdyn);

  setStringToIntegralParameter<int>("FLUID_SOLVER", "Implicit",
							   "Solving strategy for fluid",
							   tuple<std::string>("Implicit","Pressure Correction","Pressure Correction SemiImplicit"),
							   tuple<int>(fluid_solver_implicit, fluid_solver_pressurecorrection, fluid_solver_pressurecorrection_semiimplicit),&fdyn);

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
                                 "compute L2 errors of increments (relative) and residuals (absolute)",
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
                                 "COUNTERVORT",
                                 "SOLWAVE",
                                 "WAVEBREAKING",
                                 "BELTRAMI-FLOW",
                                 "KIM-MOIN-FLOW",
                                 "BREAKING-DAM"),
                               tuple<int>(0,1,2,3,4,6,7,8,9,10),
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

  setStringToIntegralParameter<int>("NEUMANNINFLOW",
                               "no",
                               "Flag to (de)activate potential Neumann inflow term(s)",
                               tuple<std::string>(
                                 "no",
                                 "yes"),
                               tuple<std::string>(
                                 "No Neumann inflow term(s)",
                                 "Neumann inflow term(s) might occur"),
                               tuple<int>(0,1),
                               &fdyn);

  setStringToIntegralParameter<int>("FSSUGRVISC","No","fine-scale subgrid viscosity",
                               tuple<std::string>(
                                 "No",
                                 "Smagorinsky_all",
                                 "Smagorinsky_small"
                                 ),
                               tuple<int>(0,1,2),
                               &fdyn);

  setStringToIntegralParameter<int>("SIMPLER","no",
                               "Switch on SIMPLE family of solvers, needs additional FLUID PRESSURE SOLVER block!",
                               yesnotuple,yesnovalue,&fdyn);

  setStringToIntegralParameter<int>("ADAPTCONV","yes",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&fdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&fdyn);

  IntParameter("UPPSS",1,"Increment for visualisation (unused)",&fdyn);
  IntParameter("UPOUT",1,"Increment for writing solution to output file (unused)",&fdyn);
  IntParameter("UPRES",1,"Increment for writing solution",&fdyn);
  IntParameter("RESSTEP",0,"Restart Step (unused)",&fdyn);
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
  DoubleParameter("START_THETA",1.0,"Time integration factor for starting scheme",&fdyn);
 /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& andyn = list->sublist("ARTERIAL DYNAMIC",false,"");

  setStringToIntegralParameter<int>("DYNAMICTYP","ExpTaylorGalerkin",
                               "Explicit Taylor Galerkin Scheme",
                               tuple<std::string>(
                                 "ExpTaylorGalerkin"
                                 ),
                               tuple<int>(
                                typ_tay_gal
                                ),
                               &andyn);

  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&andyn);
  IntParameter("NUMSTEP",0,"Number of Time Steps",&andyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&andyn);
  IntParameter("UPRES",1,"Increment for writing solution",&andyn);
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
                                 "yes_cross",
                                 "cross_rhs",
                                 "cross_complete"
                                 ),
                               tuple<std::string>(
                                 "No cross-stress term",
                                 "Include the cross-stress term with a linearization of the convective part",
                                 "Include cross-stress term, but only explicitly on right hand side"
                                 ),
                               tuple<int>(0,1,2,3),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("REYNOLDS-STRESS",
                               "no_reynolds",
                               "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
                               tuple<std::string>(
                                 "no_reynolds",
                                 "yes_reynolds",
                                 "reynolds_rhs",
                                 "reynolds_complete"
                                 ),
                               tuple<std::string>(
                                 "No Reynolds-stress term",
                                 "Include Reynolds-stress term explicitly on right hand side",
                                 "Include Reynolds-stress term with linearisation"
                                 ),
                               tuple<int>(0,1,2,3),
                               &fdyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU",
                               "Barrenechea_Franca_Valentin_Wall",
                               "Definition of tau_M,C",
                               tuple<std::string>(
                                 "Barrenechea_Franca_Valentin_Wall",
                                 "BFVW_gradient_based_hk",
                                 "Smoothed_FBVW",
                                 "FBVW_without_dt",
                                 "Franca_Barrenechea_Valentin_Codina",
                                 "Bazilevs",
                                 "Codina"),
                               tuple<std::string>(
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea; tau_C: Wall",
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea; tau_C: Wall, gradien based element length",
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea (smoothed max operator using exp function); tau_C: Wall",
                                 "tau_M : Barrenechea, Valentin, Franca, Barrenechea; tau_C: Wall; no dt contribution",
                                 "tau_Mp: Barrenechea, Valentin; tau_M: Franca, Barrenechea; tau_C: Codina"  ,
                                 "tau_M and tau_C (Bazilevs, based on G_ij and g_i)",
                                 "tau_M and tau_C: Codina")  ,
                                    tuple<int>(0,1,2,3,4,5,6),
                               &fdyn_stab);

  setStringToIntegralParameter<INPAR::FLUID::TauType>("TAUTYPE","Franca_Barrenechea_Valentin_Wall",
                               "Type of definition of stabilization parameter",
                               tuple<std::string>(
                                 "Franca_Barrenechea_Valentin_Wall",
                                 "Bazilevs"),
                               tuple<INPAR::FLUID::TauType>(
                                 INPAR::FLUID::tautype_franca_barrenechea_valentin_wall,
                                 INPAR::FLUID::tautype_bazilevs),
                               &fdyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU",
                               "element_center",
                               "Location where tau is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate tau at element center",
                                 "evaluate tau at integration point")  ,
                                tuple<int>(0,1),
                               &fdyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT",
                               "element_center",
                               "Location where material law is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate material law at element center",
                                 "evaluate material law at integration point")  ,
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

  DoubleParameter("C_TURBPRANDTL",1.0,"(Constant) turbulent Prandtl number for the Smagorinsky model in scalar transport.",&fdyn_turbu);

  setStringToIntegralParameter<int>(
    "CANONICAL_FLOW",
    "no",
    "Sampling is different for different canonical flows \n--- so specify what kind of flow you've got",
    tuple<std::string>(
      "no",
      "channel_flow_of_height_2",
      "lid_driven_cavity",
      "backward_facing_step",
      "square_cylinder",
      "square_cylinder_nurbs",
      "rotating_circular_cylinder_nurbs",
      "loma_channel_flow_of_height_2",
      "loma_lid_driven_cavity",
      "loma_backward_facing_step"),
    tuple<std::string>(
      "The flow is not further specified, so spatial averaging \nand hence the standard sampling procedure is not possible",
      "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.",
      "For this flow, all statistical data are evaluated on the center lines of the xy-midplane, averaged only over time.",
      "For this flow, statistical data are evaluated on various lines, averaged over time and z.",
      "For this flow, statistical data are evaluated on various lines of the xy-midplane, averaged only over time.",
      "For this flow, statistical data are evaluated on various lines of the xy-midplane, averaged over time and eventually in one hom.direction.",
      "For this flow, statistical data is computed in concentric surfaces and averaged. in time and in one hom. direction",
      "For this low-Mach-number flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.",
      "For this low-Mach-number flow, all statistical data are evaluated on the center lines of the xy-midplane, averaged only over time.",
      "For this low-Mach-number flow, statistical data are evaluated on various lines, averaged over time and z."),
    tuple<int>(0,1,2,3,4,5,6,7,8,9),
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

  IntParameter("SAMPLING_START",10000000,"Time step after when sampling shall be started",&fdyn_turbu);
  IntParameter("SAMPLING_STOP",1,"Time step when sampling shall be stopped",&fdyn_turbu);
  IntParameter("DUMPING_PERIOD",1,"Period of time steps after which statistical data shall be dumped",&fdyn_turbu);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC",false,"");

  DoubleParameter("TIMESTEP",0.1,"",&adyn);
  IntParameter("NUMSTEP",41,"",&adyn);
  DoubleParameter("MAXTIME",4.0,"",&adyn);
  setStringToIntegralParameter<int>("ALE_TYPE","classic_lin","ale mesh movement algorithm",
                               tuple<std::string>("classic_lin","incr_lin","laplace","springs","springs_fixed_ref"),
                               tuple<int>(ALE_DYNAMIC::classic_lin,
                                          ALE_DYNAMIC::incr_lin,
                                          ALE_DYNAMIC::laplace,
                                          ALE_DYNAMIC::springs,
                                          ALE_DYNAMIC::springs_fixed_ref),
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

  setStringToIntegralParameter<INPAR::SCATRA::SolverType>("SOLVERTYPE","linear_full",
                               "type of scalar transport solver",
                               tuple<std::string>(
                                 "linear_full",
                                 "linear_incremental",
                                 "nonlinear"
                                 ),
                               tuple<INPAR::SCATRA::SolverType>(
                                   INPAR::SCATRA::solvertype_linear_full,
                                   INPAR::SCATRA::solvertype_linear_incremental,
                                   INPAR::SCATRA::solvertype_nonlinear),
                               &scatradyn);

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
  IntParameter("MATID",-1,"Material Id for automatic mesh creation",&scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::VelocityField>("VELOCITYFIELD","zero",
                               "type of velocity field used for scalar transport problems",
                               tuple<std::string>(
                                 "zero",
                                 "function",
                                 "Navier_Stokes"
                                 ),
                               tuple<INPAR::SCATRA::VelocityField>(
                                   INPAR::SCATRA::velocity_zero,
                                   INPAR::SCATRA::velocity_function,
                                   INPAR::SCATRA::velocity_Navier_Stokes),
                               &scatradyn);

  IntParameter("VELFUNCNO",-1,"function number for scalar transport velocity field",&scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::InitialField>("INITIALFIELD","zero_field",
                               "Initial Field for scalar transport problem",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function",
                                 "field_by_condition",
                                 "disturbed_field_by_function",
                                 "1D_DISCONTPV",
                                 "FVI_FERECHPRO",
                                 "RAYTAYMIXFRAC"),
                               tuple<INPAR::SCATRA::InitialField>(
                                   INPAR::SCATRA::initfield_zero_field,
                                   INPAR::SCATRA::initfield_field_by_function,
                                   INPAR::SCATRA::initfield_field_by_condition,
                                   INPAR::SCATRA::initfield_disturbed_field_by_function,
                                   INPAR::SCATRA::initfield_DISCONTPV_1D,
                                   INPAR::SCATRA::initfield_FVI_FERECHPRO,
                                   INPAR::SCATRA::initfield_RAYTAYMIXFRAC),
                               &scatradyn);

  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::CalcError>("CALCERROR","No",
                               "compute error compared to analytical solution",
                               tuple<std::string>(
                                 "No",
                                 "Kwok_Wu"
                                 ),
                               tuple<INPAR::SCATRA::CalcError>(
                                   INPAR::SCATRA::calcerror_no,
                                   INPAR::SCATRA::calcerror_Kwok_Wu),
                               &scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::FluxType>("WRITEFLUX","No","output of diffusive/total flux vectors",
                               tuple<std::string>(
                                 "No",
                                 "totalflux_domain",
                                 "diffusiveflux_domain",
                                 "totalflux_boundary",
                                 "diffusiveflux_boundary"
                                 ),
                               tuple<INPAR::SCATRA::FluxType>(
                                   INPAR::SCATRA::flux_no,
                                   INPAR::SCATRA::flux_total_domain,
                                   INPAR::SCATRA::flux_diffusive_domain,
                                   INPAR::SCATRA::flux_total_boundary,
                                   INPAR::SCATRA::flux_diffusive_boundary),
                               &scatradyn);

  BoolParameter("OUTMEAN","No","Output of mean values for scalars and density",&scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::ConvForm>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<INPAR::SCATRA::ConvForm>(
                                 INPAR::SCATRA::convform_convective,
                                 INPAR::SCATRA::convform_conservative),
                               &scatradyn);

  BoolParameter("NEUMANNINFLOW",
      "no","Flag to (de)activate potential Neumann inflow term(s)",&scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::FSSUGRDIFF>("FSSUGRDIFF",
                               "No",
                               "fine-scale subgrid diffusivity",
                               tuple<std::string>(
                                 "No",
                                 "artificial",
                                 "Smagorinsky_all",
                                 "Smagorinsky_small"
                                 ),
                               tuple<INPAR::SCATRA::FSSUGRDIFF>(
                                   INPAR::SCATRA::fssugrdiff_no,
                                   INPAR::SCATRA::fssugrdiff_artificial,
                                   INPAR::SCATRA::fssugrdiff_smagorinsky_all,
                                   INPAR::SCATRA::fssugrdiff_smagorinsky_small),
                               &scatradyn);

  BoolParameter("BLOCKPRECOND","NO",
      "Switch to block-preconditioned family of solvers, needs additional SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER block!",&scatradyn);

  setStringToIntegralParameter<INPAR::SCATRA::ScaTraType>("SCATRATYPE","Standard",
                               "Type of scalar transport problem",
                               tuple<std::string>(
                                 "Standard",
                                 "LevelSet"),
                               tuple<INPAR::SCATRA::ScaTraType>(
                                 INPAR::SCATRA::scatratype_standard,
                                 INPAR::SCATRA::scatratype_levelset),
                               &scatradyn);

  DoubleParameter("INITIALDENS",1.0,"Initial value for density",&scatradyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatra_nonlin = scatradyn.sublist(
      "NONLINEAR",
      false,
      "control parameters for solving nonlinear SCATRA problems\n");

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&scatra_nonlin);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&scatra_nonlin);
  BoolParameter("EXPLPREDICT","no","do an explicit predictor step before starting nonlinear iteration",&scatra_nonlin);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&scatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&scatra_nonlin);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn_levelset = scatradyn.sublist("LEVELSET",false,
      "control parameters for a level set function");

  setStringToIntegralParameter<INPAR::SCATRA::ReinitializationAction>("REINITIALIZATION","None",
                               "Type of reinitialization strategy for level set function",
                               tuple<std::string>(
                                 "None",
                                 "DirectDistance",
                                 "Sussman",
                                 "InterfaceProjection",
                                 "Function",
                                 "Signed_Distance_Function"),
                               tuple<INPAR::SCATRA::ReinitializationAction>(
                                 INPAR::SCATRA::reinitaction_none,
                                 INPAR::SCATRA::reinitaction_directdistance,
                                 INPAR::SCATRA::reinitaction_sussman,
                                 INPAR::SCATRA::reinitaction_interfaceprojection,
                                 INPAR::SCATRA::reinitaction_function,
                                 INPAR::SCATRA::reinitaction_signeddistancefunction),
                               &scatradyn_levelset);

  setStringToIntegralParameter<INPAR::SCATRA::MassCalculation>("MASSCALCULATION","No",
                               "Type of mass calculation",
                               tuple<std::string>(
                                 "No",
                                 "Squares",
                                 "Interpolated"),
                               tuple<INPAR::SCATRA::MassCalculation>(
                                 INPAR::SCATRA::masscalc_none,
                                 INPAR::SCATRA::masscalc_squares,
                                 INPAR::SCATRA::masscalc_interpolated),
                             &scatradyn_levelset);

/*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn_stab = scatradyn.sublist("STABILIZATION",false,"");

  // this parameter governs type of stabilization
  setStringToIntegralParameter<INPAR::SCATRA::StabType>("STABTYPE",
                                    "SUPG",
                                    "type of stabilization (if any)",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "SUPG",
                                 "GLS",
                                 "USFEM"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> only reasonable for low-Peclet-number flows",
                                 "Use SUPG",
                                 "Use GLS",
                                 "Use USFEM")  ,
                               tuple<INPAR::SCATRA::StabType>(
                                   INPAR::SCATRA::stabtype_no_stabilization,
                                   INPAR::SCATRA::stabtype_SUPG,
                                   INPAR::SCATRA::stabtype_GLS,
                                   INPAR::SCATRA::stabtype_USFEM),
                               &scatradyn_stab);

  // this parameter governs whether subgrid-scale velocity is included
  BoolParameter("SUGRVEL","no","potential incorporation of subgrid-scale velocity",&scatradyn_stab);

  // this parameter governs whether all-scale subgrid diffusivity is included
  BoolParameter("ASSUGRDIFF","no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) term",&scatradyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<INPAR::SCATRA::TauType>("DEFINITION_TAU",
                               "Franca_Valentin",
                               "Definition of tau",
                               tuple<std::string>(
                                 "Franca_Valentin",
                                 "Bazilevs",
                                 "Exact_1D",
                                 "Zero"),
                               tuple<std::string>(
                                 "tau according to Franca and Valentin (2000)",
                                 "tau according to Bazilevs et al. (2007) (based on G_ij and g_i)",
                                 "exact tau for stationary 1d problems and linear shape functions",
                                 "zero tau (no stabilizing effect)")  ,
                                tuple<INPAR::SCATRA::TauType>(
                                    INPAR::SCATRA::tau_franca_valentin,
                                    INPAR::SCATRA::tau_bazilevs,
                                    INPAR::SCATRA::tau_exact_1d,
                                    INPAR::SCATRA::tau_zero),
                               &scatradyn_stab);

  // this parameter selects the all-scale subgrid-diffusivity definition applied
  setStringToIntegralParameter<INPAR::SCATRA::AssgdType>("DEFINITION_ASSGD",
                               "artificial_linear",
                               "Definition of (all-scale) subgrid diffusivity",
                               tuple<std::string>(
                                 "artificial_linear",
                                 "Hughes_etal_86_nonlinear",
                                 "Tezduyar_Park_86_nonlinear",
                                 "doCarmo_Galeao_91_nonlinear",
                                 "Almeida_Silva_97_nonlinear"),
                               tuple<std::string>(
                                 "classical linear artificial subgrid-diffusivity",
                                 "nonlinear isotropic according to Hughes et al. (1986)",
                                 "nonlinear isotropic according to Tezduyar and Park (1986)",
                                 "nonlinear isotropic according to doCarmo and Galeao (1991)",
                                 "nonlinear isotropic according to Almeida and Silva (1997)")  ,
                                tuple<INPAR::SCATRA::AssgdType>(
                                    INPAR::SCATRA::assgd_artificial,
                                    INPAR::SCATRA::assgd_hughes,
                                    INPAR::SCATRA::assgd_tezduyar,
                                    INPAR::SCATRA::assgd_docarmo,
                                    INPAR::SCATRA::assgd_almeida),
                               &scatradyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<INPAR::SCATRA::EvalTau>("EVALUATION_TAU",
                               "element_center",
                               "Location where tau is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate tau at element center",
                                 "evaluate tau at integration point")  ,
                                tuple<INPAR::SCATRA::EvalTau>(
                                  INPAR::SCATRA::evaltau_element_center,
                                  INPAR::SCATRA::evaltau_integration_point),
                               &scatradyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<INPAR::SCATRA::EvalMat>("EVALUATION_MAT",
                               "element_center",
                               "Location where material law is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate material law at element center",
                                 "evaluate material law at integration point"),
                               tuple<INPAR::SCATRA::EvalMat>(
                                 INPAR::SCATRA::evalmat_element_center,
                                 INPAR::SCATRA::evalmat_integration_point),
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
  setStringToIntegralParameter<int>("CONSTHERMPRESS","Yes",
                               "treatment of thermodynamic pressure in time",
                               tuple<std::string>(
                                 "No_energy",
                                 "No_mass",
                                 "Yes"
                                 ),
                               tuple<int>(0,1,2),
                               &lomacontrol);
  setStringToIntegralParameter<int>(
    "CANONICAL_FLOW",
    "no",
    "Information on special flows",
    tuple<std::string>(
      "no",
      "loma_channel_flow_of_height_2",
      "loma_lid_driven_cavity",
      "loma_backward_facing_step"),
    tuple<std::string>(
      "The flow is not further specified.",
      "low-Mach-number in channel",
      "low-Mach-number flow in lid-driven cavity",
      "low-Mach-number flow over a backward-facing step"),
    tuple<int>(0,1,2,3),
    &lomacontrol);
  IntParameter("SAMPLING_START",1,"Time step after when sampling shall be started",&lomacontrol);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& elchcontrol = list->sublist(
      "ELCH CONTROL",
      false,
      "control parameters for electrochemistry problems\n");

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&elchcontrol);
  IntParameter("NUMSTEP",24,"Total number of time steps",&elchcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&elchcontrol);
  IntParameter("ITEMAX",10,"Maximum number of outer iterations",&elchcontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&elchcontrol);
  DoubleParameter("CONVTOL",1e-6,"Convergence check tolerance for outer loop",&elchcontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&elchcontrol);
  DoubleParameter("TEMPERATURE",298.0,"Constant temperature (Kelvin)",&elchcontrol);
  BoolParameter("MOVINGBOUNDARY","No","ELCH algorithm for deforming meshes",&elchcontrol);
  DoubleParameter("MOLARVOLUME",0.0,"Molar volume for electrode shape change computations",&elchcontrol);
  setStringToIntegralParameter<INPAR::ELCH::NatConv>("NATURAL_CONVECTION","No",
                               "Include natural convection effects",
                               tuple<string>(
                                 "No",
                                 "Natural_Convection_substance",
                                 "Natural_Convection_ion"
                                 ),
                                 tuple<INPAR::ELCH::NatConv>(
                                     INPAR::ELCH::natural_convection_no,
                                     INPAR::ELCH::natural_convection_substance,
                                     INPAR::ELCH::natural_convection_ion),
                               &elchcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrol = list->sublist("COMBUSTION CONTROL",false,
      "control parameters for a combustion problem");

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

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolfluid = combustcontrol.sublist("COMBUSTION FLUID",false,
      "control parameters for the fluid field of a combustion problem");

  setStringToIntegralParameter<INPAR::COMBUST::CombustionType>("COMBUSTTYPE","Premixed_Combustion",
                               "Type of combustion problem",
                               tuple<std::string>(
                                 "Premixed_Combustion",
                                 "Two_Phase_Flow"),
                               tuple<INPAR::COMBUST::CombustionType>(
                                 INPAR::COMBUST::combusttype_premixedcombustion,
                                 INPAR::COMBUST::combusttype_twophaseflow),
                               &combustcontrolfluid);
  DoubleParameter("LAMINAR_FLAMESPEED",1.0,"The laminar flamespeed incorporates all chemical kinetics into the problem for now",&combustcontrolfluid);
  DoubleParameter("MARKSTEIN_LENGTH",0.0,"The Markstein length takes flame curvature into account",&combustcontrolfluid);
  DoubleParameter("NITSCHE_VELOCITY",0.0,"Nitsche parameter to stabilize/penalize the velocity jump",&combustcontrolfluid);
  DoubleParameter("NITSCHE_PRESSURE",0.0,"Nitsche parameter to stabilize/penalize the pressure jump",&combustcontrolfluid);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolgfunc = combustcontrol.sublist("COMBUSTION GFUNCTION",false,
      "control parameters for the G-function (level set) field of a combustion problem");

  setStringToIntegralParameter<INPAR::COMBUST::ReInitialActionGfunc>("REINITIALIZATION","Signed_Distance_Function",
                               "Type of reinitialization strategy for level set function",
                               tuple<std::string>(
                                 "None",
                                 "Function",
                                 "Signed_Distance_Function"),
                               tuple<INPAR::COMBUST::ReInitialActionGfunc>(
                                 INPAR::COMBUST::reinitaction_none,
                                 INPAR::COMBUST::reinitaction_byfunction,
                                 INPAR::COMBUST::reinitaction_signeddistancefunction),
                               &combustcontrolgfunc);
  IntParameter("REINITFUNCNO",-1,"function number for reinitialization of level set (G-function) field",&combustcontrolgfunc);
  IntParameter("REINITINTERVAL",1,"reinitialization interval",&combustcontrolgfunc);
  setStringToIntegralParameter<int>("REFINEMENT","No","Turn refinement strategy for level set function on/off",
                                     yesnotuple,yesnovalue,&combustcontrolgfunc);
  IntParameter("REFINEMENTLEVEL",-1,"number of refinement level for refinement strategy",&combustcontrolgfunc);

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
                                 "iter_monolithicfluidsplit",
                                 "iter_monolithiclagrange",
                                 "iter_monolithicstructuresplit",
                                 "iter_monolithicxfem",
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
                                 fsi_iter_monolithicfluidsplit,
                                 fsi_iter_monolithiclagrange,
                                 fsi_iter_monolithicstructuresplit,
                                 fsi_iter_monolithicxfem,
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
                                 "FSIAMG",
                                 "PartitionedAitken",
                                 "PartitionedVectorExtrapolation",
                                 "PartitionedJacobianFreeNewtonKrylov",
                                 "BGSAitken",
                                 "BGSVectorExtrapolation",
                                 "BGSJacobianFreeNewtonKrylov"
                                 ),
                               tuple<INPAR::FSI::LinearBlockSolver>(
                                 INPAR::FSI::PreconditionedKrylov,
                                 INPAR::FSI::FSIAMG,
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
  IntParameter("UPPSS",1,"Increment for visualization (unused)",&fsidyn);
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
                  "to the nonlinear convergence test.",
                  &fsidyn);

  DoubleParameter("ADAPTIVEDIST",0.0,
                  "Required distance for adaptive convergence check in Newton-type FSI.\n"
                  "This is the improvement we want to achieve in the linear extrapolation of the\n"
                  "adaptive convergence check. Set to zero to avoid the adaptive check altogether.",
                  &fsidyn);

  // monolithic preconditioner parameter

  setNumericStringParameter("STRUCTPCOMEGA","1.0 1.0 1.0 1.0",
                  "Relaxation factor for Richardson iteration on structural block in MFSI block preconditioner",
                  &fsidyn);
  setNumericStringParameter("STRUCTPCITER","0 0 0 0",
               "Number of Richardson iterations on structural block in MFSI block preconditioner",
               &fsidyn);
  setNumericStringParameter("FLUIDPCOMEGA","1.0 1.0 1.0 1.0",
                  "Relaxation factor for Richardson iteration on fluid block in MFSI block preconditioner",
                  &fsidyn);
  setNumericStringParameter("FLUIDPCITER","0 0 0 0",
               "Number of Richardson iterations on fluid block in MFSI block preconditioner",
               &fsidyn);
  setNumericStringParameter("ALEPCOMEGA","1.0 1.0 1.0 1.0",
                  "Relaxation factor for Richardson iteration on ale block in MFSI block preconditioner",
                  &fsidyn);
  setNumericStringParameter("ALEPCITER","0 0 0 0",
               "Number of Richardson iterations on ale block in MFSI block preconditioner",
               &fsidyn);

  setNumericStringParameter("PCOMEGA","1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on whole MFSI block preconditioner",
                            &fsidyn);
  setNumericStringParameter("PCITER","1 1 3",
                            "Number of Richardson iterations on whole MFSI block preconditioner",
                            &fsidyn);

  //DoubleParameter("PCOMEGA",1.,
  //                "Relaxation factor for Richardson iteration on whole MFSI block preconditioner",
  //                &fsidyn);
  //IntParameter("PCITER",1,
  //             "Number of Richardson iterations on whole MFSI block preconditioner",
  //             &fsidyn);

  setStringToIntegralParameter<int>("INFNORMSCALING","Yes","Scale Blocks in Mono-FSI with row infnorm?",
                                     yesnotuple,yesnovalue,&fsidyn);
  setStringToIntegralParameter<int>("SYMMETRICPRECOND","No","Symmetric block GS preconditioner in monolithic FSI or ordinary GS",
                                     yesnotuple,yesnovalue,&fsidyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL",false,"");

  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT","Yes","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("DLM_CONDENSATION","Yes","Do you want to condense the distributed Lagrange multiplier?",
                                 yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("INCOMP_PROJECTION","No","Do you want to project the old velocity to an incompressible velocity field?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("FAST_INTEGRATION","No","Do you want to save gausspoints to speed up computation?",
                                   yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("CONDEST","No","Do you want to estimate the condition number? It is somewhat costly.",
                                   yesnotuple,yesnovalue,&xfem_general);
  DoubleParameter("volumeRatioLimit",1.0e-2,"don't enrich nodes of elements, when less than this fraction of the element is on one side of the interface",&xfem_general);
  DoubleParameter("boundaryRatioLimit",1.0e-4,"don't enrich element, when less than this area fraction is within this element",&xfem_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fluidsolver = list->sublist("FLUID SOLVER",false,"");
  SetValidSolverParameters(fluidsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fluidpsolver = list->sublist("FLUID PRESSURE SOLVER",false,"pressure solver parameters for SIMPLE preconditioning");
  SetValidSolverParameters(fluidpsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluidprojsolver = list->sublist("XFLUID PROJECTION SOLVER",false,"");
  SetValidSolverParameters(xfluidprojsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& structsolver = list->sublist("STRUCT SOLVER",false,"");
  SetValidSolverParameters(structsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& alesolver = list->sublist("ALE SOLVER",false,"");
  SetValidSolverParameters(alesolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& thermalsolver = list->sublist("THERMAL SOLVER",false,"linear solver for thermal problems");
  SetValidSolverParameters(thermalsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatrasolver = list->sublist("SCALAR TRANSPORT SOLVER",false,"solver parameters for scalar transport problems");
  SetValidSolverParameters(scatrasolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatrapotsolver = list->sublist("SCALAR TRANSPORT ELECTRIC POTENTIAL SOLVER",false,"solver parameters for block-preconditioning");
  SetValidSolverParameters(scatrapotsolver);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& artnetsolver = list->sublist("ARTERY NETWORK SOLVER",false,"");
  SetValidSolverParameters(artnetsolver);
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

  setNumericStringParameter("ML_SMOTIMES","1 1 1 1 1",
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
  SetPrintEqualSign(list,true);

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
  SetPrintEqualSign(direction,true);

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
  SetPrintEqualSign(newton,true);

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
  SetPrintEqualSign(steepestdescent,true);

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
  SetPrintEqualSign(linesearch,true);

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
  SetPrintEqualSign(fullstep,true);

  {
    DoubleParameter("Full Step",1.0,"length of a full step",&fullstep);
  }

  // sub-sub-list "Backtrack"
  Teuchos::ParameterList& backtrack = linesearch.sublist("Backtrack",false,"");
  SetPrintEqualSign(backtrack,true);

  {
    DoubleParameter("Default Step",1.0,"starting step length",&backtrack);
    DoubleParameter("Minimum Step",1.0e-12,"minimum acceptable step length",&backtrack);
    DoubleParameter("Recovery Step",1.0,"step to take when the line search fails (defaults to value for \"Default Step\")",&backtrack);
    IntParameter("Max Iters",50,"maximum number of iterations (i.e., RHS computations)",&backtrack);
    DoubleParameter("Reduction Factor",0.5,"A multiplier between zero and one that reduces the step size between line search iterations",&backtrack);
  }

  // sub-sub-list "Polynomial"
  Teuchos::ParameterList& polynomial = linesearch.sublist("Polynomial",false,"");
  SetPrintEqualSign(polynomial,true);

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
  SetPrintEqualSign(morethuente,true);

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
  SetPrintEqualSign(trustregion,true);

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
  SetPrintEqualSign(printing,true);

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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string DRT::INPUT::PrintEqualSign()
{
  return "*PrintEqualSign*";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::SetPrintEqualSign(Teuchos::ParameterList& list, const bool& pes)
{
  std::string printequalsign = PrintEqualSign();
  list.set<bool>(printequalsign,pes);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::INPUT::NeedToPrintEqualSign(const Teuchos::ParameterList& list)
{
  const std::string printequalsign = PrintEqualSign();
  bool pes = false;
  try
  {
    pes = list.get<bool>(printequalsign);
  }
  catch (Teuchos::Exceptions::InvalidParameter)
  {
    pes = false;
  }
  return pes;
}

#endif
