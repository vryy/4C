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

#include <Teuchos_Array.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_any.hpp>

#include "drt_validparameters.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem_enums.H"
#include "../drt_inpar/inpar_ale.H"
#include "../drt_inpar/inpar_artnet.H"
#include "../drt_inpar/inpar_solver.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_potential.H"
#include "../drt_inpar/inpar_problemtype.H"
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_inpar/inpar_turbulence.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_searchtree.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_inpar/inpar_poroscatra.H"
#include "../drt_inpar/inpar_ssi.H"

#include <AztecOO.h>

/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
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
//! Print help message
/*----------------------------------------------------------------------*/
void PrintHelpMessage()
{
#ifdef DEBUG
  char baci_build[] = "baci-debug";
#else
  char baci_build[] = "baci-release";
#endif

  std::cout << "NAME" << std::endl
            << "\t" << baci_build << " - simulate just about anything" << std::endl
            << std::endl
            << "SYNOPSIS" << std::endl
            << "\t" << baci_build << " [-h] [--help] [-p] [--parameters] [-d] [--datfile] [-ngroup=x] [-glayout=a,b,c,...] [-nptype=parallelism_type]" << std::endl
            << "\t\tdat_name output_name [restart=y] [restartfrom=restart_file_name] [ dat_name0 output_name0 [restart=y] [restartfrom=restart_file_name] ... ] [--interactive]" << std::endl
            << std::endl
            << "DESCRIPTION" << std::endl
            << "\tThe am besten simulation tool in the world." << std::endl
            << std::endl
            << "OPTIONS" << std::endl
            << "\t--help or -h" << std::endl
            << "\t\tPrint this message." << std::endl
            << std::endl
            << "\t--parameters or -p" << std::endl
            << "\t\tPrint a list of all available parameters for use in a dat_file." << std::endl
            << std::endl
            << "\t--datfile or -d" << std::endl
            << "\t\tPrint example dat_file with all available parameters." << std::endl
            << std::endl
            << "\t-ngroup=x" << std::endl
            << "\t\tSpecify the number of groups for nested parallelism. (default: 1)" << std::endl
            << std::endl
            << "\t-glayout=a,b,c,..." << std::endl
            << "\t\tSpecify the number of processors per group. Argument \"-ngroup\" is mandatory and must be preceding. (default: equal distribution)" << std::endl
            << std::endl
            << "\t-nptype=parallelism_type" << std::endl
            << "\t\tAvailable options: \"separateDatFiles\", \"everyGroupReadDatFile\" and \"copyDatFile\"; Must be set if \"-ngroup\" > 1." << std::endl
            << std::endl
            << "\tdat_name" << std::endl
            << "\t\tName of the input file (Usually *.dat)" << std::endl
            << std::endl
            << "\toutput_name" << std::endl
            << "\t\tPrefix of your output files." << std::endl
            << std::endl
            << "\trestart=y" << std::endl
            << "\t\tRestart the simulation from step y. It always refers to the previously defined dat_name and output_name. (default: 0 or from dat_name)" << std::endl
            << std::endl
            << "\trestartfrom=restart_file_name" << std::endl
            << "\t\tRestart the simulation from the files prefixed with restart_file_name. (default: output_name)" << std::endl
            << std::endl
            << "\t--interactive" << std::endl
            << "\t\tBaci waits at the beginning for keyboard input. Helpful for parallel debugging when attaching to a single job. Must be specified at the end in the command line." << std::endl
            << std::endl
            << "SEE ALSO" << std::endl
            << "\tguides/reports/global_report.pdf" << std::endl
            << std::endl
            << "BUGS" << std::endl
            << "\t100% bug free since 1964." << std::endl
            << std::endl
            << "TIPS" << std::endl
            << "\tCan be obtain from a friendly colleague." << std::endl
            << std::endl
            << "\tAlso, espresso may be donated to room MW1236." << std::endl;

  return;
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
            for (int i=0; i<(int)values->size(); ++i)
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
              for (int i=0; i<(int)values->size(); ++i)
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
//! Print function
/*----------------------------------------------------------------------*/
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

#if 0
namespace DRT
{
namespace INPUT
{

  template <class T>
  Teuchos::Array<T> tuple( const T & t1 )
  {
    Teuchos::Array<T> a;
    a.reserve(1);
    a.push_back(t1);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2 )
  {
    Teuchos::Array<T> a;
    a.reserve(2);
    a.push_back(t1);
    a.push_back(t2);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3 )
  {
    Teuchos::Array<T> a;
    a.reserve(3);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4 )
  {
    Teuchos::Array<T> a;
    a.reserve(4);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5 )
  {
    Teuchos::Array<T> a;
    a.reserve(5);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6 )
  {
    Teuchos::Array<T> a;
    a.reserve(6);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7 )
  {
    Teuchos::Array<T> a;
    a.reserve(7);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8 )
  {
    Teuchos::Array<T> a;
    a.reserve(8);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9 )
  {
    Teuchos::Array<T> a;
    a.reserve(9);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9, const T & t10 )
  {
    Teuchos::Array<T> a;
    a.reserve(10);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    a.push_back(t10);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9, const T & t10, const T & t11 )
  {
    Teuchos::Array<T> a;
    a.reserve(11);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    a.push_back(t10);
    a.push_back(t11);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9, const T & t10, const T & t11, const T & t12 )
  {
    Teuchos::Array<T> a;
    a.reserve(12);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    a.push_back(t10);
    a.push_back(t11);
    a.push_back(t12);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9, const T & t10, const T & t11, const T & t12, const T & t13 )
  {
    Teuchos::Array<T> a;
    a.reserve(13);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    a.push_back(t10);
    a.push_back(t11);
    a.push_back(t12);
    a.push_back(t13);
    return a;
  }

  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9, const T & t10, const T & t11, const T & t12, const T & t13, const T & t14 )
  {
    Teuchos::Array<T> a;
    a.reserve(14);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    a.push_back(t10);
    a.push_back(t11);
    a.push_back(t12);
    a.push_back(t13);
    a.push_back(t14);
    return a;
  }


  template <class T>
  Teuchos::Array<T> tuple( const T & t1, const T & t2, const T & t3, const T & t4, const T & t5, const T & t6, const T & t7, const T & t8, const T & t9, const T & t10, const T & t11, const T & t12, const T & t13, const T & t14, const T & t15 )
  {
    Teuchos::Array<T> a;
    a.reserve(15);
    a.push_back(t1);
    a.push_back(t2);
    a.push_back(t3);
    a.push_back(t4);
    a.push_back(t5);
    a.push_back(t6);
    a.push_back(t7);
    a.push_back(t8);
    a.push_back(t9);
    a.push_back(t10);
    a.push_back(t11);
    a.push_back(t12);
    a.push_back(t13);
    a.push_back(t14);
    a.push_back(t15);
    return a;
  }

}
}
#endif

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
  IntParameter("NUMARTNETDIS",1,"Number of meshes in arterial network field",&discret);
  IntParameter("NUMTHERMDIS",1,"Number of meshes in thermal field",&discret);
  IntParameter("NUMAIRWAYSDIS",1,"Number of meshes in reduced dimensional airways network field",&discret);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& size = list->sublist("PROBLEM SIZE",false,"");

  IntParameter("DIM",3,"2d or 3d problem",&size);

  // deactivate all the follwing (unused) parameters one day
  // they are nice as general info in the input file but should not
  // read into a parameter list. Misuse is possible
  IntParameter("ELEMENTS",0,"Total number of elements",&size);
  IntParameter("NODES",0,"Total number of nodes",&size);
  IntParameter("NPATCHES",0,"number of nurbs patches",&size);
  IntParameter("MATERIALS",0,"number of materials",&size);
  IntParameter("NUMDF",3,"maximum number of degrees of freedom",&size);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& type = list->sublist("PROBLEM TYP",false,"");

  {
    Teuchos::Array<std::string> name;
    Teuchos::Array<int> label;

    // fill the arrays
    {
      std::map<std::string,PROBLEM_TYP> map = DRT::StringToProblemTypeMap();
      std::map<std::string,PROBLEM_TYP>::const_iterator i;
      for (i = map.begin(); i != map.end();++i)
      {
        name. push_back(i->first);
        label.push_back(i->second);
      }
    }

    setStringToIntegralParameter<int>(
      "PROBLEMTYP",
      "Fluid_Structure_Interaction",
      "",
      name,
      label,
      &type);
  }

  IntParameter("NUMFIELD",1,"",&type); // unused. to be deleted
  IntParameter("RESTART",0,"",&type);
  setStringToIntegralParameter<int>("SHAPEFCT","Polynomial","Defines the function spaces for the spatial approximation",
                               tuple<std::string>("Polynomial","Nurbs"),
                               tuple<int>(1,0),
                               &type);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ps = list->sublist("PATIENT SPECIFIC",false,"");

  setStringToIntegralParameter<int>("PATSPEC","No",
                                    "Triggers application of patient specific tools in discretization construction",
                                    yesnotuple,yesnovalue,&ps);

  BoolParameter("REMODEL","No","Turn remodeling on/off",&ps);

  IntParameter("MAXHULUMEN",0,"max HU value within the blood lumen",&ps);
  setNumericStringParameter("CENTERLINEFILE","name.txt",
                            "filename of file containing centerline points",
                            &ps);

  setStringToIntegralParameter<int>("CALCSTRENGTH","No","Calculate strength on/off",yesnotuple,yesnovalue,&ps);
  DoubleParameter("AAA_SUBRENDIA",22.01,"subrenal diameter of the AAA",&ps);
  setStringToIntegralParameter<int>("FAMILYHIST","No","Does the patient have AAA family history",yesnotuple,yesnovalue,&ps);
  setStringToIntegralParameter<int>("MALE_PATIENT","Yes","Is the patient a male?",yesnotuple,yesnovalue,&ps);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& io = list->sublist("IO",false,"");

  // are these needed?
  setStringToIntegralParameter<int>("OUTPUT_OUT","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("OUTPUT_GID","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("OUTPUT_BIN","No","",yesnotuple,yesnovalue,&io);

  setStringToIntegralParameter<int>("STRUCT_DISP","Yes","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_SE","No","output of strain energy",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_STRESS","No","",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "Cauchy","cauchy",
                                                  "2PK", "2pk"),
                               tuple<int>(INPAR::STR::stress_none,INPAR::STR::stress_none,INPAR::STR::stress_none,
                                                             INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,
                                                             INPAR::STR::stress_cauchy,INPAR::STR::stress_cauchy,
                                                             INPAR::STR::stress_2pk,INPAR::STR::stress_2pk),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_STRAIN","No","",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "EA","ea",
                                                  "GL", "gl"),
                               tuple<int>(INPAR::STR::strain_none,INPAR::STR::strain_none,INPAR::STR::strain_none,
                                                             INPAR::STR::strain_gl,INPAR::STR::strain_gl,INPAR::STR::strain_gl,
                                                             INPAR::STR::strain_ea,INPAR::STR::strain_ea,
                                                             INPAR::STR::strain_gl,INPAR::STR::strain_gl),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_PLASTIC_STRAIN","No","",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "EA","ea",
                                                  "GL", "gl"),
                               tuple<int>(INPAR::STR::strain_none,INPAR::STR::strain_none,INPAR::STR::strain_none,
                                                             INPAR::STR::strain_gl,INPAR::STR::strain_gl,INPAR::STR::strain_gl,
                                                             INPAR::STR::strain_ea,INPAR::STR::strain_ea,
                                                             INPAR::STR::strain_gl,INPAR::STR::strain_gl),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_SURFACTANT","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_SM_DISP","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_SOL","Yes","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_WALL_SHEAR_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_VIS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("ALE_DISP","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("THERM_TEMPERATURE","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("THERM_HEATFLUX","None","",
                               tuple<std::string>("None",
                                                  "No",
                                                  "NO",
                                                  "no",
                                                  "Current",
                                                  "Initial"),
                               tuple<int>(INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_none,
                                                               INPAR::THR::heatflux_current,
                                                               INPAR::THR::heatflux_initial),
                               &io);
  setStringToIntegralParameter<int>("THERM_TEMPGRAD","None","",
                               tuple<std::string>("None",
                                                  "No",
                                                  "NO",
                                                  "no",
                                                  "Current",
                                                  "Initial"),
                               tuple<int>(INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_none,
                                                               INPAR::THR::tempgrad_current,
                                                               INPAR::THR::tempgrad_initial),
                               &io);

  IntParameter("FILESTEPS",1000,"",&io);
  IntParameter("STDOUTEVRY",1,"Print to screen every n step",&io);

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

  /*--------------------------------------------------------------------*/
  /* parameters for NOX - non-linear solution */
  Teuchos::ParameterList& snox = list->sublist("STRUCT NOX",false,"");
  SetValidNoxParameters(snox);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& sdyn = list->sublist("STRUCTURAL DYNAMIC",false,"");

  setStringToIntegralParameter<int>("DYNAMICTYP","GenAlpha",
                               "type of time integration control",
                               tuple<std::string>(
                                 "Statics",
                                 "GenAlpha",
                                 "OneStepTheta",
                                 "GEMM",
                                 "ExplicitEuler",
                                 "CentrDiff",
                                 "AdamsBashforth2",
                                 "EulerMaruyama",
                                 "EulerImpStoch",
                                 "StatMech"),
                               tuple<int>(
                                 INPAR::STR::dyna_statics,
                                 INPAR::STR::dyna_genalpha,
                                 INPAR::STR::dyna_onesteptheta,
                                 INPAR::STR::dyna_gemm,
                                 INPAR::STR::dyna_expleuler,
                                 INPAR::STR::dyna_centrdiff,
                                 INPAR::STR::dyna_ab2,
                                 INPAR::STR::dyna_euma,
                                 INPAR::STR::dyna_euimsto,
                                 INPAR::STR::dyna_statmech),
                               &sdyn);

  setStringToIntegralParameter<int>("PRESTRESS","none","prestressing takes values none mulf id",
                               tuple<std::string>("none","None","NONE",
                                                  "mulf","Mulf","MULF",
                                                  "id","Id","ID"),
                               tuple<int>(INPAR::STR::prestress_none,INPAR::STR::prestress_none,INPAR::STR::prestress_none,
                                                            INPAR::STR::prestress_mulf,INPAR::STR::prestress_mulf,INPAR::STR::prestress_mulf,
                                                            INPAR::STR::prestress_id,INPAR::STR::prestress_id,INPAR::STR::prestress_id),
                               &sdyn);

  DoubleParameter("PRESTRESSTIME",0.0,"time to switch from pre to post stressing",&sdyn);

  // a temporary flag
  setStringToIntegralParameter<int>("ADAPTERDRIVE","No",
                                    "TEMPORARY FLAG: Switch on time integration driver based on ADAPTER::Structure rather than independent implementation",
                                    yesnotuple,yesnovalue,&sdyn);

  // Output type
  IntParameter("EIGEN",0,"EIGEN make eigenanalysis of the initial dynamic system",&sdyn);
  IntParameter("RESULTSEVRY",1,"save displacements and contact forces every RESULTSEVRY steps",&sdyn);
  IntParameter("RESEVRYERGY",0,"write system energies every requested step",&sdyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&sdyn);
  // Time loop control
  DoubleParameter("TIMESTEP",0.05,"time step size",&sdyn);
  IntParameter("NUMSTEP",200,"maximum number of steps",&sdyn);
  DoubleParameter("MAXTIME",5.0,"maximum time",&sdyn);
  // Damping
  setStringToIntegralParameter<int>("DAMPING","No",
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
                               tuple<int>(
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

  DoubleParameter("TOLDISP",1.0E-10,
                  "tolerance in the displacement norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<int>("NORM_DISP","Abs","type of norm for displacement convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<int>(
                                 INPAR::STR::convnorm_abs,
                                 INPAR::STR::convnorm_rel,
                                 INPAR::STR::convnorm_mix),
                               &sdyn);

  DoubleParameter("TOLRES",1.0E-08,
                  "tolerance in the residual norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for residual convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<int>(
                                 INPAR::STR::convnorm_abs,
                                 INPAR::STR::convnorm_rel,
                                 INPAR::STR::convnorm_mix),
                               &sdyn);

  DoubleParameter("TOLPRE",1.0E-08,
                  "tolerance in pressure norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<int>("NORM_PRES","Abs","type of norm for pressure convergence check",
                               tuple<std::string>(
                                 "Abs"),
                               tuple<int>(
                                 INPAR::STR::convnorm_abs),
                               &sdyn);

  DoubleParameter("TOLINCO",1.0E-08,
                  "tolerance in the incompressible residual norm for the newton iteration",
                  &sdyn);
  setStringToIntegralParameter<int>("NORM_INCO","Abs","type of norm for incompressible residual convergence check",
                               tuple<std::string>(
                                 "Abs"),
                               tuple<int>(
                                 INPAR::STR::convnorm_abs),
                               &sdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_DISPPRES","And","binary operator to combine pressure and displacement values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 INPAR::STR::bop_and,
                                 INPAR::STR::bop_or),
                               &sdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINCO","And","binary operator to combine force and incompressible residual",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 INPAR::STR::bop_and,
                                 INPAR::STR::bop_or),
                               &sdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFDISP","And","binary operator to combine displacement and residual force values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 INPAR::STR::bop_and,
                                 INPAR::STR::bop_or),
                               &sdyn);

  setStringToIntegralParameter<int>("STC_SCALING","no",
      "Scaled director conditioning for thin shell structures",
      tuple<std::string>(
        "no",
        "No",
        "NO",
        "Symmetric",
        "Right"),
      tuple<int>(
        INPAR::STR::stc_none,
        INPAR::STR::stc_none,
        INPAR::STR::stc_none,
        INPAR::STR::stc_currsym,
        INPAR::STR::stc_curr),
      &sdyn);


  IntParameter("STC_LAYER",1,
               "number of STC layers for multilayer case",
               &sdyn);

  DoubleParameter("PTCDT",0.1,
                  "pseudo time step for pseudo transient continuation (PTC) stabilized Newton procedure",
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
  setStringToIntegralParameter<int>("ITERNORM","L2","type of norm to be applied to residuals",
                               tuple<std::string>(
                                 "L1",
                                 "L2",
                                 "Rms",
                                 "Inf"),
                               tuple<int>(
                                 INPAR::STR::norm_l1,
                                 INPAR::STR::norm_l2,
                                 INPAR::STR::norm_rms,
                                 INPAR::STR::norm_inf),
                               &sdyn);
  setStringToIntegralParameter<int>("DIVERCONT","No",
                                    "Go on with time integration even if Newton-Raphson iteration failed",
                                    yesnotuple,yesnovalue,&sdyn);

  setStringToIntegralParameter<int>("NLNSOL","fullnewton","Nonlinear solution technique",
                               tuple<std::string>(
                                 "vague",
                                 "fullnewton",
                                 "lsnewton",
                                 "oppnewton",
                                 "modnewton",
                                 "ptc",
                                 "newtonlinuzawa",
                                 "augmentedlagrange",
                                 "NoxNewtonLineSearch",
                                 "noxgeneral"),
                               tuple<int>(
                                 INPAR::STR::soltech_vague,
                                 INPAR::STR::soltech_newtonfull,
                                 INPAR::STR::soltech_newtonls,
                                 INPAR::STR::soltech_newtonopp,
                                 INPAR::STR::soltech_newtonmod,
                                 INPAR::STR::soltech_ptc,
                                 INPAR::STR::soltech_newtonuzawalin,
                                 INPAR::STR::soltech_newtonuzawanonlin,
                                 INPAR::STR::soltech_noxnewtonlinesearch,
                                 INPAR::STR::soltech_noxgeneral),
                               &sdyn);

  setStringToIntegralParameter<int>("CONTROLTYPE","load","load, disp, arc1, arc2 control",
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
                               tuple<int>(
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

  setStringToIntegralParameter<int>("PREDICT","ConstDis","Type of predictor",
                               tuple<std::string>(
                                 "Vague",
                                 "ConstDis",
                                 "ConstVel",
                                 "ConstAcc",
                                 "ConstDisVelAcc",
                                 "TangDis",
                                 "ConstDisPres",
                                 "ConstDisVelAccPres"),
                               tuple<int>(
                                 INPAR::STR::pred_vague,
                                 INPAR::STR::pred_constdis,
                                 INPAR::STR::pred_constvel,
                                 INPAR::STR::pred_constacc,
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
  setStringToIntegralParameter<int>("UZAWAALGO","direct","",
                                 tuple<std::string>(
                                   "uzawa",
                                   "simple",
                                   "direct"),
                                 tuple<int>(
                                   INPAR::STR::consolve_uzawa,
                                   INPAR::STR::consolve_simple,
                                   INPAR::STR::consolve_direct),
                                 &sdyn);

  // convergence criteria adaptivity
  setStringToIntegralParameter<int>("ADAPTCONV","No",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&sdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&sdyn);

  setStringToIntegralParameter<int>("LUMPMASS","No",
                               "Lump the mass matrix for explicit time integration",
                               yesnotuple,yesnovalue,&sdyn);

  setStringToIntegralParameter<int>("MODIFIEDEXPLEULER","Yes",
                               "Use the modified explicit Euler time integration scheme",
                               yesnotuple,yesnovalue,&sdyn);

  // linear solver id used for structural problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for structural problems",&sdyn);

  /*--------------------------------------------------------------------*/
  /* parameters for time step size adaptivity in structural dynamics */
  Teuchos::ParameterList& tap = sdyn.sublist("TIMEADAPTIVITY",false,"");
  SetValidTimeAdaptivityParameters(tap);

  /*----------------------------------------------------------------------*/
  /* parameters for generalised-alpha structural integrator */
  Teuchos::ParameterList& genalpha = sdyn.sublist("GENALPHA",false,"");

  setStringToIntegralParameter<int>("GENAVG","ImrLike",
                               "mid-average type of internal forces",
                               tuple<std::string>(
                                 "Vague",
                                 "ImrLike",
                                 "TrLike"),
                               tuple<int>(
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

  DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,0.5]",&gemm);
  DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&gemm);
  DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&gemm);
  DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&gemm);
  DoubleParameter("XI",0.0,"generalisation factor in [0,1)",&gemm);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& iap = list->sublist("INVERSE ANALYSIS",false,"");

  // Inverse Analysis
  setStringToIntegralParameter<int>("INV_ANALYSIS","none",
                               "types of inverse analysis and on/off switch",
                               tuple<std::string>(
                                 "none",
                                 "lung",
                                 "gen"),
                               tuple<int>(
                                 INPAR::STR::inv_none,
                                 INPAR::STR::inv_lung,
                                 INPAR::STR::inv_generalized),
                               &iap);


  DoubleParameter("MC_X_0",0.0,"measured displacment of the tension testing in x dir",&iap);
  DoubleParameter("MC_X_1",0.0,"measured displacment of the tension testing in x dir",&iap);
  DoubleParameter("MC_X_2",0.0,"measured displacment of the tension testing in x dir",&iap);
  DoubleParameter("MC_Y_0",0.0,"measured displacment of the tension testing in y dir",&iap);
  DoubleParameter("MC_Y_1",0.0,"measured displacment of the tension testing in y dir",&iap);
  DoubleParameter("MC_Y_2",0.0,"measured displacment of the tension testing in y dir",&iap);



  // tolerance for inv_analysis
  DoubleParameter("INV_ANA_TOL",1.0,"tolerance for inverse analysis",&iap);
  IntParameter("INV_ANA_MAX_RUN",100,"max iterations for inverse analysis",&iap);

  // perturbation parameters
  DoubleParameter("INV_ALPHA",1.0e-3,"perturbation parameters",&iap);
  DoubleParameter("INV_BETA",1.0e-3,"perturbation parameters",&iap);

  // initial regularization parameter
  DoubleParameter("INV_INITREG",1.0,"initial regularization parameter",&iap);


  setNumericStringParameter("MONITORFILE","none.monitor",
                            "filename of file containing measured displacements",
                            &iap);

  setNumericStringParameter("INV_LIST","-1",
                            "IDs of materials that have to be fitted",
                            &iap);

  setNumericStringParameter("INV_EH_LIST","-1",
                            "IDs of materials that have to be fitted",
                            &iap);

  setStringToIntegralParameter<int>("NEW_FILES","yes",
                                    "new result files for each run",
                                    yesnotuple,yesnovalue,&iap);
  /*----------------------------------------------------------------------*/

  /* parameters for multi-level monte carlo */
  Teuchos::ParameterList& mlmcp = list->sublist("MULTI LEVEL MONTE CARLO",false,"");

  setStringToIntegralParameter<int>("MLMC","no",
                                    "perform multi level monte carlo analysis",
                                    yesnotuple,yesnovalue,&mlmcp);
  IntParameter("NUMRUNS",200,"Number of Monte Carlo runs",&mlmcp);

  //IntParameter("NUMLEVELS",2,"Number of levels",&mlmcp);
  setStringToIntegralParameter<int>("DIFF_TO_LOWER_LEVEL","no","calculate difference to next lower level",yesnotuple,yesnovalue,&mlmcp);
  IntParameter("START_RUN",0,"Run to start calculating the difference to lower level", &mlmcp);
  //IntParameter("END_RUN",0,"Run to stop calculating the difference to lower level", &mlmcp);
  // NUMLEVEL additional inputfiles are read name must be standard_inputfilename+_level_i.dat
  setNumericStringParameter("DISCRETIZATION_FOR_PROLONGATION","filename.dat",
                            "filename of.dat file which contains discretization to which the results are prolongated",
                            &mlmcp);
  setNumericStringParameter("OUTPUT_FILE_OF_LOWER_LEVEL","level0",
                            "filename of controlfiles of next lower level",
                            &mlmcp);
  setStringToIntegralParameter<int>("PROLONGATERES","No",
                                    "Prolongate Displacements to finest Discretization",
                                    yesnotuple,yesnovalue,&mlmcp);
  //Parameter for Newton loop to find background element
  IntParameter("ITENODEINELE",20,"Number iteration in Newton loop to determine background element",&mlmcp);
  DoubleParameter("CONVTOL",10e-5,"Convergence tolerance for Newton loop",&mlmcp);
  IntParameter("INITRANDOMSEED",1000,"Random seed for first Monte Carlo run",&mlmcp);
  IntParameter("LEVELNUMBER",0,"Level number for Multi Level Monte Carlo", &mlmcp);
  IntParameter("WRITESTATS",1000,"Write statistics to file every WRITESTATS (only for polongated Dis",&mlmcp);
  setStringToIntegralParameter<int>("REDUCED_OUTPUT","NO",
                                    "Write reduced Coarse Level Output, i.e. no mesh stresses, just disp",
                                    yesnotuple,yesnovalue,&mlmcp);
  // Parameters to simulate random fields
  IntParameter("RANDOM_FIELD_DIMENSION",3,"Dimension of Random Field 2 or 3",&mlmcp);
  IntParameter("SIZE_PER_DIM",512,"Number of points per dimension",&mlmcp);
  DoubleParameter("CORRLENGTH",30,"Correlation length of Random Field",&mlmcp);
  IntParameter("NUM_COS_TERMS",64,"Number of terms in geometric row ",&mlmcp);
  setStringToIntegralParameter<int>("SPECTRAL_MATCHING","Yes",
                                    "Perform spectral matching with psd 1 yes ",
                                    yesnotuple,yesnovalue,&mlmcp);
  DoubleParameter("SIGMA",1.0,"sigma of random field",&mlmcp);
  DoubleParameter("MEAN",0.0,"Mean value of random field",&mlmcp);
  setStringToIntegralParameter<int>("CORRSTRUCT","gaussian",
                              "Correlation structure of random field",
                              tuple<std::string>("Gaussian","gaussian"),
                              tuple<int>(
                                  INPAR::MLMC::corr_gaussian,INPAR::MLMC::corr_gaussian),
                              &mlmcp);
  setStringToIntegralParameter<int>("MARGINALPDF","gaussian","Target marginal probability distribution function",
                                    tuple<std::string>("Gaussian","gaussian",
                                                        "Beta", "beta",
                                                        "lognormal", "Lognormal"),
                                    tuple<int>(
                                        INPAR::MLMC::pdf_gaussian,INPAR::MLMC::pdf_gaussian,
                                        INPAR::MLMC::pdf_beta,INPAR::MLMC::pdf_beta,
                                        INPAR::MLMC::pdf_lognormal,INPAR::MLMC::pdf_lognormal),
                                 &mlmcp);
  DoubleParameter("NONGAUSSPARAM1",0,"First parameter for non-gaussian pdf",&mlmcp);
  DoubleParameter("NONGAUSSPARAM2",0,"Second parameter for non-gaussian pdf",&mlmcp);
  DoubleParameter("KAPPA_U",6.283185307,"CUTOFF WAVE NUMBER FOR PSD",&mlmcp);
  setStringToIntegralParameter<int>("CALC_METHOD","FFT",
                                "Calculation method for the random field",
                                tuple<std::string>("FFT","COS"),
                                tuple<int>(
                                    INPAR::MLMC::calc_m_fft,INPAR::MLMC::calc_m_cos),
                                &mlmcp);
  setNumericStringParameter("OUTPUT_ELEMENT_IDS","-1",
                              "Set ID's of Output Elements, default is -1 which is none",
                              &mlmcp);

  /*----------------------------------------------------------------------*/
  /* parameters for meshtying and contact */
  Teuchos::ParameterList& scontact = list->sublist("MESHTYING AND CONTACT",false,"");

  // linear solver id used for contact/meshtying problems (both for fluid and structure meshtying as well as structural contact)
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for meshtying and contact",&scontact);

  // linear solver id used for contact/meshtying problems in saddlepoint formulation used for solving LAGRANGE multipliers,
  // used as SIMPLER preconditioner in FLUID DYNAMICS
  IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for meshtying and contact in saddlepoint formulation",&scontact);

  setStringToIntegralParameter<int>("APPLICATION","None","Type of contact or meshtying app",
       tuple<std::string>("None","none",
                          "MortarContact","mortarcontact",
                          "MortarMeshtying","mortarmeshtying",
                          "BeamContact","beamcontact"),
       tuple<int>(
                  INPAR::CONTACT::app_none,INPAR::CONTACT::app_none,
                  INPAR::CONTACT::app_mortarcontact,INPAR::CONTACT::app_mortarcontact,
                  INPAR::CONTACT::app_mortarmeshtying,INPAR::CONTACT::app_mortarmeshtying,
                  INPAR::CONTACT::app_beamcontact, INPAR::CONTACT::app_beamcontact),
       &scontact);

  setStringToIntegralParameter<int>("FRICTION","None","Type of friction law",
      tuple<std::string>("None","none",
                         "Stick","stick",
                         "Tresca","tresca",
                         "Coulomb","coulomb"),
      tuple<int>(
                 INPAR::CONTACT::friction_none,INPAR::CONTACT::friction_none,
                 INPAR::CONTACT::friction_stick,INPAR::CONTACT::friction_stick,
                 INPAR::CONTACT::friction_tresca,INPAR::CONTACT::friction_tresca,
                 INPAR::CONTACT::friction_coulomb,INPAR::CONTACT::friction_coulomb),
      &scontact);

  DoubleParameter("FRBOUND",0.0,"Friction bound for Tresca friction",&scontact);
  DoubleParameter("FRCOEFF",0.0,"Friction coefficient for Coulomb friction",&scontact);

  setStringToIntegralParameter<int>("WEAR","None","Type of wear law",
      tuple<std::string>("None","none",
                         "Archard","archard"),
      tuple<int>(
                 INPAR::CONTACT::wear_none,INPAR::CONTACT::wear_none,
                 INPAR::CONTACT::wear_archard,INPAR::CONTACT::wear_archard),
      &scontact);

  DoubleParameter("WEARCOEFF",0.0,"Wear coefficient",&scontact);

  setStringToIntegralParameter<int>("STRATEGY","LagrangianMultipliers","Type of employed solving strategy",
        tuple<std::string>("LagrangianMultipliers","lagrange", "Lagrange",
                           "PenaltyMethod","penalty", "Penalty",
                           "AugmentedLagrange","augmented", "Augmented"),
        tuple<int>(
                INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult,
                INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty,
                INPAR::CONTACT::solution_auglag, INPAR::CONTACT::solution_auglag, INPAR::CONTACT::solution_auglag),
        &scontact);

  setStringToIntegralParameter<int>("SHAPEFCN","Dual","Type of employed set of shape functions",
        tuple<std::string>("Dual", "dual",
                           "Standard", "standard", "std"),
        tuple<int>(
                INPAR::MORTAR::shape_dual, INPAR::MORTAR::shape_dual,
                INPAR::MORTAR::shape_standard, INPAR::MORTAR::shape_standard, INPAR::MORTAR::shape_standard),
        &scontact);

  setStringToIntegralParameter<int>("SYSTEM","Condensed","Type of linear system setup / solution",
        tuple<std::string>("Condensed","condensed", "cond",
                           "SaddlePointCoupled","saddlepointcoupled", "spcoupled",
                           "SaddlePointSimpler","saddlepointsimpler", "spsimpler"),
        tuple<int>(
                INPAR::CONTACT::system_condensed, INPAR::CONTACT::system_condensed, INPAR::CONTACT::system_condensed,
                INPAR::CONTACT::system_spcoupled, INPAR::CONTACT::system_spcoupled, INPAR::CONTACT::system_spcoupled,
                INPAR::CONTACT::system_spsimpler, INPAR::CONTACT::system_spsimpler, INPAR::CONTACT::system_spsimpler),
        &scontact);

  DoubleParameter("PENALTYPARAM",0.0,"Penalty parameter for penalty / augmented solution strategy",&scontact);
  DoubleParameter("PENALTYPARAMTAN",0.0,"Tangential penalty parameter for penalty / augmented solution strategy",&scontact);
  IntParameter("UZAWAMAXSTEPS",10,"Maximum no. of Uzawa steps for augmented / Uzawa solution strategy",&scontact);
  DoubleParameter("UZAWACONSTRTOL",1.0e-8,"Tolerance of constraint norm for augmented / Uzawa solution strategy",&scontact);

  setStringToIntegralParameter<int>("FULL_LINEARIZATION","Yes","If chosen full linearization of contact is applied",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON","Yes","If chosen semi-smooth Newton concept is applied",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("SEMI_SMOOTH_CN",1.0,"Weighting factor cn for semi-smooth PDASS",&scontact);
  DoubleParameter("SEMI_SMOOTH_CT",1.0,"Weighting factor ct for semi-smooth PDASS",&scontact);

  setStringToIntegralParameter<int>("SEARCH_ALGORITHM","Binarytree","Type of contact search",
       tuple<std::string>("BruteForce","bruteforce",
                          "BruteForceEleBased","bruteforceelebased",
                          "BinaryTree","Binarytree","binarytree"),
       tuple<int>(INPAR::MORTAR::search_bfele,INPAR::MORTAR::search_bfele,
                  INPAR::MORTAR::search_bfele,INPAR::MORTAR::search_bfele,
                  INPAR::MORTAR::search_binarytree,INPAR::MORTAR::search_binarytree,
                  INPAR::MORTAR::search_binarytree),
       &scontact);

  DoubleParameter("SEARCH_PARAM",0.3,"Radius / Bounding volume inflation for contact search",&scontact);

  setStringToIntegralParameter<int>("COUPLING_AUXPLANE","Yes","If chosen auxiliary planes are used for 3D coupling",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("LAGMULT_QUAD","undefined","Type of LM interpolation/weighting function",
       tuple<std::string>("undefined",
                          "quad_quad", "quadratic_quadratic",
                          "quad_pwlin", "quadratic_piecewiselinear",
                          "quad_lin", "quadratic_linear",
                          "pwlin_pwlin", "piecewiselinear_piecewiselinear",
                          "lin_lin","linear_linear"),
       tuple<int>(
                  INPAR::MORTAR::lagmult_undefined,
                  INPAR::MORTAR::lagmult_quad_quad, INPAR::MORTAR::lagmult_quad_quad,
                  INPAR::MORTAR::lagmult_quad_pwlin, INPAR::MORTAR::lagmult_quad_pwlin,
                  INPAR::MORTAR::lagmult_quad_lin, INPAR::MORTAR::lagmult_quad_lin,
                  INPAR::MORTAR::lagmult_pwlin_pwlin, INPAR::MORTAR::lagmult_pwlin_pwlin,
                  INPAR::MORTAR::lagmult_lin_lin, INPAR::MORTAR::lagmult_lin_lin),
       &scontact);

  setStringToIntegralParameter<int>("CROSSPOINTS","No","If chosen, multipliers are removed from crosspoints / edge nodes",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("EMOUTPUT","None","Type of energy and momentum output",
      tuple<std::string>("None","none", "No", "no",
                         "Screen", "screen",
                         "File", "file",
                         "Both", "both"),
      tuple<int>(
              INPAR::CONTACT::output_none, INPAR::CONTACT::output_none,
              INPAR::CONTACT::output_none, INPAR::CONTACT::output_none,
              INPAR::CONTACT::output_screen, INPAR::CONTACT::output_screen,
              INPAR::CONTACT::output_file, INPAR::CONTACT::output_file,
              INPAR::CONTACT::output_both, INPAR::CONTACT::output_both),
      &scontact);

  setStringToIntegralParameter<int>("ERROR_NORMS","None","Choice of analytical solution for error norm computation",
      tuple<std::string>("None","none", "No", "no",
                         "Zero", "zero",
                         "Bending", "bending",
                         "Sphere", "sphere",
                         "Thick", "thick"),
      tuple<int>(
              INPAR::MORTAR::errornorms_none, INPAR::MORTAR::errornorms_none,
              INPAR::MORTAR::errornorms_none, INPAR::MORTAR::errornorms_none,
              INPAR::MORTAR::errornorms_zero, INPAR::MORTAR::errornorms_zero,
              INPAR::MORTAR::errornorms_bending, INPAR::MORTAR::errornorms_bending,
              INPAR::MORTAR::errornorms_sphere, INPAR::MORTAR::errornorms_sphere,
              INPAR::MORTAR::errornorms_thicksphere, INPAR::MORTAR::errornorms_thicksphere),
      &scontact);

  setStringToIntegralParameter<int>("PARALLEL_REDIST","Static","Type of redistribution algorithm",
      tuple<std::string>("None","none", "No", "no",
                         "Static", "static",
                         "Dynamic", "dynamic"),
      tuple<int>(
              INPAR::MORTAR::parredist_none, INPAR::MORTAR::parredist_none,
              INPAR::MORTAR::parredist_none, INPAR::MORTAR::parredist_none,
              INPAR::MORTAR::parredist_static, INPAR::MORTAR::parredist_static,
              INPAR::MORTAR::parredist_dynamic, INPAR::MORTAR::parredist_dynamic),
      &scontact);

  DoubleParameter("MAX_BALANCE",2.0,"Maximum value of load balance measure before parallel redistribution",&scontact);
  IntParameter("MIN_ELEPROC",0,"Minimum no. of elements per processor for parallel redistribution",&scontact);

  DoubleParameter("HEATTRANSSLAVE",0.0,"Heat transfer parameter for slave side in thermal contact",&scontact);
  DoubleParameter("HEATTRANSMASTER",0.0,"Heat transfer parameter for master side in thermal contact",&scontact);

  setStringToIntegralParameter<int>("THERMOLAGMULT","Yes","Lagrange Multipliers are applied for thermo-contact",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("BEAMS_NEWGAP","No","choose between original or enhanced gapfunction",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("BEAMS_SMOOTHING","None","Application of smoothed tangent field",
       tuple<std::string>("None","none",
                          "Smoothed","smoothed",
                          "Partially","partially",
                          "Cpp", "cpp"),
       tuple<int>(
                  INPAR::CONTACT::bsm_none,INPAR::CONTACT::bsm_none,
                  INPAR::CONTACT::bsm_smoothed,INPAR::CONTACT::bsm_smoothed,
                  INPAR::CONTACT::bsm_partially,INPAR::CONTACT::bsm_partially,
                  INPAR::CONTACT::bsm_cpp,INPAR::CONTACT::bsm_cpp),
       &scontact);

  // enable octree search and determine type of bounding box (aabb = axis aligned, cobb = cylindrical oriented)
  setStringToIntegralParameter<int>("BEAMS_OCTREE","None","octree and bounding box type for octree search routine",
       tuple<std::string>("None","none","octree_axisaligned","octree_cylorient","octree_spherical"),
       tuple<int>(INPAR::CONTACT::boct_none,INPAR::CONTACT::boct_none,
                  INPAR::CONTACT::boct_aabb,INPAR::CONTACT::boct_cobb,INPAR::CONTACT::boct_spbb),
       &scontact);

  DoubleParameter("BEAMS_EXTFAC",1.05,"extrusion factor of the bounding box",&scontact);
  DoubleParameter("BEAMS_RADFAC",1.05,"radius extrusion factor of the bounding box",&scontact);
  IntParameter("BEAMS_TREEDEPTH",6,"max, tree depth of the octree",&scontact);
  IntParameter("BEAMS_BOXESINOCT",8,"max number of bounding boxes in any leaf octant",&scontact);

  setStringToIntegralParameter<int>("INITCONTACTBYGAP","No","Initialize init contact by weighted gap vector",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("INITCONTACTGAPVALUE",0.0,"Value for initialization of init contact set with gap vector",&scontact);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& interaction_potential = list->sublist("INTERACTION POTENTIAL",false,"");

  // read if surfaces , volumes or both including fluid should be considered
  setStringToIntegralParameter<int>("POTENTIAL_TYPE","Surface","Type of interaction potential",
                                tuple<std::string>("Surface",
                                                   "Volume",
                                                   "Surfacevolume",
                                                   "Surface_fsi",
                                                   "Volume_fsi",
                                                   "Surfacevolume_fsi"),
                                tuple<int>(
                                   INPAR::POTENTIAL::potential_surface,
                                   INPAR::POTENTIAL::potential_volume,
                                   INPAR::POTENTIAL::potential_surfacevolume,
                                   INPAR::POTENTIAL::potential_surface_fsi,
                                   INPAR::POTENTIAL::potential_volume_fsi,
                                   INPAR::POTENTIAL::potential_surfacevolume_fsi),
                                &interaction_potential);

  // approximation method
  setStringToIntegralParameter<int>("APPROXIMATION_TYPE","None","Type of approximation",
                                tuple<std::string>("None",
                                                   "Surface_approx",
                                                   "Point_approx"),
                                tuple<int>(
                                           INPAR::POTENTIAL::approximation_none,
                                           INPAR::POTENTIAL::approximation_surface,
                                           INPAR::POTENTIAL::approximation_point),
                                &interaction_potential);

  // switches on the analytical solution computation for two van der waals spheres or membranes
  setStringToIntegralParameter<int>("ANALYTICALSOLUTION","None", "Type of analytical solution"
                                 "computes analytical solutions for two Van der Waals spheres or membranes",
                                 tuple<std::string>("None",
                                                    "Sphere",
                                                    "Membrane"),
                                 tuple<int>(
                                            INPAR::POTENTIAL::solution_none,
                                            INPAR::POTENTIAL::solution_sphere,
                                            INPAR::POTENTIAL::solution_membrane),
                                            &interaction_potential);
  // use 2D integration for pseudo 3D
  setStringToIntegralParameter<int>("PSEUDO3D","no",
                                     "use 2D integration for pseudo 3D",
                                     yesnotuple,yesnovalue,&interaction_potential);

  // radius of can der Waals spheres for analytical testing
  DoubleParameter(  "VDW_RADIUS",0.0,
                    "radius of van der Waals spheres",
                    &interaction_potential);

  // thickness of sphericael Waals membranes for analytical testing
  DoubleParameter(  "THICKNESS",0.0,
                    "membrane thickness",
                    &interaction_potential);

  // number of atoms or molecules offset
  DoubleParameter(  "N_OFFSET",0.0,
                    "number of atoms or molecules offset",
                    &interaction_potential);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& statmech = list->sublist("STATISTICAL MECHANICS",false,"");

  //Reading kind of background fluid stream in the thermal bath
  setStringToIntegralParameter<int>("THERMALBATH","None","Type of thermal bath applied to elements",
                               //listing possible strings in input file in category THERMALBATH
                               tuple<std::string>("None","none",
                                                  "Uniform","uniform",
                                                  "ShearFlow","shearflow","Shearflow"),
                               //translating input strings into BACI input parameters
                               tuple<int>(INPAR::STATMECH::thermalbath_none,INPAR::STATMECH::thermalbath_none,
                                          INPAR::STATMECH::thermalbath_uniform,INPAR::STATMECH::thermalbath_uniform,
                                          INPAR::STATMECH::thermalbath_shearflow,INPAR::STATMECH::thermalbath_shearflow,INPAR::STATMECH::thermalbath_shearflow),
                               &statmech);
  //Reading which kind of special output should be written to files
  setStringToIntegralParameter<int>("SPECIAL_OUTPUT","None","kind of special statistical output data written into files",
                                 //listing possible strings in input file in category SPECIAL_OUTPUT
                                 tuple<std::string>("None","none",
                                                    "endtoend_log",
                                                    "anisotropic",
                                                    "orientationcorrelation",
                                                    "endtoend_const",
                                                    "viscoelasticity",
                                                    "densitydensitycorr",
                                                    "octree",
                                                    "avgdistloom",
                                                    "coverageloom",
                                                    "distandcoverloom",
                                                    "attractionloom"),
                                 //translating input strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::statout_none,INPAR::STATMECH::statout_none,
                                            INPAR::STATMECH::statout_endtoendlog,
                                            INPAR::STATMECH::statout_anisotropic,
                                            INPAR::STATMECH::statout_orientationcorrelation,
                                            INPAR::STATMECH::statout_endtoendconst,
                                            INPAR::STATMECH::statout_viscoelasticity,
                                            INPAR::STATMECH::statout_densitydensitycorr,
                                            INPAR::STATMECH::statout_octree,
                                            INPAR::STATMECH::statout_avgdistloom,
                                            INPAR::STATMECH::statout_coverageloom,
                                            INPAR::STATMECH::statout_distandcoverloom,
                                            INPAR::STATMECH::statout_attractionloom),
                                 &statmech);
  //Reading which kind of friction model should be applied
  setStringToIntegralParameter<int>("FRICTION_MODEL","none","friction model for polymer dynamics",
                                 //listing possible strings in input file in category FRICTION_MODEL
                                 tuple<std::string>("none",
                                                    "isotropiclumped",
                                                    "isotropicconsistent",
                                                    "anisotropicconsistent"),
                                 //translating input strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::frictionmodel_none,
                                            INPAR::STATMECH::frictionmodel_isotropiclumped,
                                            INPAR::STATMECH::frictionmodel_isotropicconsistent,
                                            INPAR::STATMECH::frictionmodel_anisotropicconsistent),
                                            &statmech);
  //Reading which kind of friction model should be applied
  setStringToIntegralParameter<int>("DBCTYPE","std","Dirichlet BC type applied",
                                 //listing possible strings in input file in category FRICTION_MODEL
                                 tuple<std::string>("none",
                                                    "std",
                                                    "shearfixed",
                                                    "sheartrans",
                                                    "pinnodes"),
                                 //translating input strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::dbctype_none,
                                            INPAR::STATMECH::dbctype_std,
                                            INPAR::STATMECH::dbctype_shearfixed,
                                            INPAR::STATMECH::dbctype_sheartrans,
                                            INPAR::STATMECH::dbctype_pinnodes),
                                            &statmech);
  //time after which writing of statistical output is started
  DoubleParameter("STARTTIMEOUT",0.0,"Time after which writing of statistical output is started",&statmech);
  //time after which certain action in simulation (e.g. DBCs in viscoelastic simulations) are started
  DoubleParameter("STARTTIMEACT",0.0,"Time after which certain action in simulation is started",&statmech);
  //time after which certain action in simulation (e.g. DBCs in viscoelastic simulations) are started
  DoubleParameter("EQUILIBTIME",0.0,"Time until which no crosslinkers are set",&statmech);
  //time after which certain action in simulation (e.g. DBCs in viscoelastic simulations) are started
  DoubleParameter("KTSWITCHTIME",0.0,"Time when KT value is changed to KTACT",&statmech);
  //alternative post-STARTTIME time step size
  DoubleParameter("DELTA_T_NEW",0.0,"A new time step size that comes into play once DBCs are have been activated",&statmech);
  //Reading whether dynamics remodelling of cross linker distribution takes place
  setStringToIntegralParameter<int>("DYN_CROSSLINKERS","No","If chosen cross linker proteins are added and removed in each time step",
                               yesnotuple,yesnovalue,&statmech);
  //Unlinking is dependent on internal forces of the linkers, too
  setStringToIntegralParameter<int>("FORCEDEPUNLINKING","No","Turns force-based unlinking of crosslinks on and off",
                               yesnotuple,yesnovalue,&statmech);
  //Toggles the use of internodal binding spots
  setStringToIntegralParameter<int>("INTERNODALBSPOTS","No","If yes, the four-noded beam element is applied which allows linker positions between FE nodes.",
                               yesnotuple,yesnovalue,&statmech);
  //Toggles helical binding spot structure of the actin filament
  setStringToIntegralParameter<int>("HELICALBINDINGSTRUCT","No","Turns double-helical binding spot geometry on and off",
                               yesnotuple,yesnovalue,&statmech);
  //Toggles helical binding spot structure of the actin filament
  setStringToIntegralParameter<int>("LOOMSETUP","No","Turns specific routines for loom network on and off",
                               yesnotuple,yesnovalue,&statmech);
  //Rise per monomer in the actin double helix according to Howard, p. 125
  DoubleParameter("RISEPERBSPOT",0.00277,"rise per monomer in the actin one-start helix",&statmech);
  //Rotation per monomer in the actin double helix according to Howard, p. 125
  DoubleParameter("ROTPERBSPOT",-2.8999,"rotation per monomer in the actin double-helix",&statmech);
  //angular offset of the binding spot orientation (constant for each filament)
  DoubleParameter("BSPOTOFFSET",0.0,"angular offset of the binding spot orientation (constant for each filament)",&statmech);
  //angle between binding spot orientation and the surface of the cone-shaped binding spot reactive volume
  DoubleParameter("PHIBSPOT",0.524,"angle between binding spot orientation and the surface of the cone-shaped binding spot reactive volume",&statmech);
  //Reading double parameter for shear flow field
  DoubleParameter("SHEARAMPLITUDE",0.0,"Shear amplitude of flow in z-direction; note: not amplitude of displacement, but of shear strain!",&statmech);
  //Reading double parameter for viscosity of background fluid
  DoubleParameter("ETA",0.0,"viscosity",&statmech);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&statmech);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KTACT",0.0,"thermal energy for t>=STARTTIMEACT",&statmech);
  //Reading double parameter for crosslinker on-rate at the beginning
  DoubleParameter("K_ON_start",0.0,"crosslinker on-rate at the end",&statmech);
  //Reading double parameter for crosslinker on-rate at the end
  DoubleParameter("K_ON_end",0.0,"crosslinker on-rate at the end",&statmech);
  //Reading double parameter for crosslinker off-rate at the beginning
  DoubleParameter("K_OFF_start",0.0,"crosslinker off-rate at the beginning",&statmech);
  //Reading double parameter for crosslinker off-rate at the end
  DoubleParameter("K_OFF_end",0.0,"crosslinker off-rate at the end",&statmech);
  //Reading double parameter for crosslinker off-rate at the end
  DoubleParameter("K_ON_SELF",0.0,"crosslinker on-rate for crosslinkers with both bonds on same filament",&statmech);
  //number of overall crosslink molecules in the boundary volume
  IntParameter("N_crosslink",0,"number of crosslinkers for switching on- and off-rates; if molecule diffusion model is used: number of crosslink molecules",&statmech);
  //number of overall crosslink molecules in the boundary volume
  IntParameter("INITOCCUPIEDBSPOTS",0,"binding spots occupied by (singly-bound) crosslinkers before the first time step",&statmech);
  //number by which the number of crosslinkers is reduced.
  IntParameter("REDUCECROSSLINKSBY",0,"number of crosslinker elements by which the overall number of crosslinker is reduced.",&statmech);
  //Reading double parameter for crosslinker protein mean length
  DoubleParameter("R_LINK",0.0,"Mean distance between two nodes connected by a crosslinker",&statmech);
  //Absolute value of difference between maximal/minimal and mean cross linker length
  DoubleParameter("DeltaR_LINK",0.0,"Absolute value of difference between maximal/minimal and mean cross linker length",&statmech);
  // upper bound of the interval within which uniformly distributed random numbers are generated
  DoubleParameter("MaxRandValue",0.0,"Upper bound of the interval within which uniformly distributed random numbers are generated (usually equal to PeriodLength)",&statmech);
  // Three values representing the size of the periodic box in each spatial direction
  setNumericStringParameter("PERIODLENGTH","0.0 0.0 0.0",
                            "Values representing the size of the periodic box in each spatial direction",&statmech);
  //angle between filament axes at crosslinked points with zero potential energy
  DoubleParameter("PHI0",0.0,"equilibrium angle between crosslinker axis and filament at each binding site",&statmech);
  //only angles in the range PHI0 +/- PHIODEV are admitted at all for the angle PHI between filament axes at crosslinked points; the default value for this parameter is 2*pi so that by default any value is admitted
  DoubleParameter("PHI0DEV",6.28,"only angles in the range PHI0 +/- PHIODEV",&statmech);
  //stiffness of orientation potential of crosslinkers
  DoubleParameter("CORIENT",0.0,"stiffness of orientation potential of crosslinkers",&statmech);
 //Young's modulus of crosslinkers
  DoubleParameter("ELINK",0.0,"Moment of inertia of area of crosslinkers",&statmech);
  //Moment of inertia of area of crosslinkers
  DoubleParameter("ILINK",0.0,"Moment of inertia of area of crosslinkers",&statmech);
  //Polar moment of inertia of area of crosslinkers
  DoubleParameter("IPLINK",0.0,"Polar moment of inertia of area of crosslinkers",&statmech);
  //Cross section of crosslinkers
  DoubleParameter("ALINK",0.0,"Cross section of crosslinkers",&statmech);
  //Parameter for PTC according to Cyron,Wall (2011):Numerical method for the simulation of the Brownian dynamics of rod-like microstructures with three dimensional nonlinear beam elements
  DoubleParameter("CTRANSPTC0",0.0,"PTC factor for translational DOF in first iteration step",&statmech);
  //Parameter for PTC according to Cyron,Wall (2011):Numerical method for the simulation of the Brownian dynamics of rod-like microstructures with three dimensional nonlinear beam elements
  DoubleParameter("CROTPTC0",0.145,"PTC factor for rotational DOF in first iteration step",&statmech);
  //Parameter for PTC according to Cyron,Wall (2011):Numerical method for the simulation of the Brownian dynamics of rod-like microstructures with three dimensional nonlinear beam elements
  DoubleParameter("ALPHAPTC",6.0,"exponent of power law for reduction of PTC factor",&statmech);
  //Makes filaments and crosslinkers be plotted by that factor thicker than they are acutally
  DoubleParameter("PlotFactorThick",0.0,"Makes filaments and crosslinkers be plotted by that factor thicker than they are acutally",&statmech);
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("CHECKORIENT","No","If chosen crosslinkers are set only after check of orientation of linked filaments",
                               yesnotuple,yesnovalue,&statmech);
  //Gmsh Output switch
  setStringToIntegralParameter<int>("GMSHOUTPUT","No","If chosen gmsh output is generated.",
                               yesnotuple,yesnovalue,&statmech);
  // toggling Gmsh Output for structure detection
  setStringToIntegralParameter<int>("GMSHNETSTRUCT","No","If chosen, special gmsh visualization for network structure types is generated.",
                               yesnotuple,yesnovalue,&statmech);
  setStringToIntegralParameter<int>("FIXEDDIRICHNODES","Yes","If chosen, the set of Dirichlet Nodes is fixed and is not updated/changed anymore.",
                               yesnotuple,yesnovalue,&statmech);
  //Number of time steps between two special outputs written
  IntParameter("OUTPUTINTERVALS",1,"Number of time steps between two special outputs written",&statmech);
  //Number of time steps between two gmsh outputs written
  IntParameter("GMSHOUTINTERVALS",100,"Number of time steps between two gmsh outputs written",&statmech);
  //Reading direction of oscillatory motion that DBC nodes are subjected to (we need this when using periodic BCs)
  IntParameter("OSCILLDIR",-1,"Global spatial direction of oscillatory motion by Dirichlet BCs",&statmech);
  //Reading time curve number for oscillatory motion
  IntParameter("CURVENUMBER",0,"Specifies Time Curve number of oscillatory motion",&statmech);
  //Reading number of elements that are taken into account when applying Dirichlet Conditions (useful to avoid redundant evaluation)
  // when Crosslink elements are added or the bead-spring-model is used
  IntParameter("NUM_EVAL_ELEMENTS",-1,"number of elements that are taken into account when applying Dirichlet Conditions",&statmech);
  // number of partitions along the edge length of the volume determining the resolution of the search grid
  IntParameter("SEARCHRES",1,"leads to the indexing of SEARCHRES^3 cubic volume partitions",&statmech);
  //Reading direction of oscillatory motion that DBC nodes are subjected to (we need this when using periodic BCs)
  IntParameter("INITIALSEED",0,"Integer value which guarantuees reproducable random number, default 0",&statmech);
  // number of histogram bins for post-analysis
  IntParameter("HISTOGRAMBINS",1,"number of bins for histograms showing the density-density-correlation-function",&statmech);
  // number of raster point along for ddcorr output boundary box shift -> n³ points in volume
  IntParameter("NUMRASTERPOINTS",3,"number of bins for histograms showing the density-density-correlation-function",&statmech);
  //Reading whether DBCs shall be applied to broken elements
  setStringToIntegralParameter<int>("PERIODICDBC","No","If chosen, Point DBCs are applied to the nodes of discontinuous elements",
                               yesnotuple,yesnovalue,&statmech);
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("FIXEDSEED","No","If chosen fixed seed for random numbers in each time step is applied",
                               yesnotuple,yesnovalue,&statmech);
  //Reading whether beam contact is switched on or not
  setStringToIntegralParameter<int>("BEAMCONTACT","No","If chosen beam contact is calculated",yesnotuple,yesnovalue,&statmech);
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& tdyn = list->sublist("THERMAL DYNAMIC",false,"");

  setStringToIntegralParameter<int>("DYNAMICTYP","OneStepTheta",
                               "type of time integration control",
                               tuple<std::string>(
                                 "Statics",
                                 "OneStepTheta",
                                 "GEMM",
                                 "GenAlpha",
                                 "ExplicitEuler"
                                 ),
                               tuple<int>(
                                  INPAR::THR::dyna_statics,
                                  INPAR::THR::dyna_onesteptheta,
                                  INPAR::THR::dyna_gemm,
                                  INPAR::THR::dyna_genalpha,
                                  INPAR::THR::dyna_expleuler),
                               &tdyn);

  // Output type
  IntParameter("RESEVRYGLOB",1,"save temperature and other global quantities every RESEVRYGLOB steps",&tdyn);
  IntParameter("RESEVRYERGY",0,"write system energies every requested step",&tdyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&tdyn);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Field for thermal problem",
                               tuple<std::string>(
                                "zero_field",
                                 "field_by_function",
                                 "field_by_condition"
                                 ),
                               tuple<int>(
                                   INPAR::THR::initfield_zero_field,
                                   INPAR::THR::initfield_field_by_function,
                                   INPAR::THR::initfield_field_by_condition
                                   ),
                               &tdyn);
  IntParameter("INITFUNCNO",-1,"function number for thermal initial field",&tdyn);

  // Time loop control
  DoubleParameter("TIMESTEP",0.05,"time step size",&tdyn);
  IntParameter("NUMSTEP",200,"maximum number of steps",&tdyn);
  DoubleParameter("MAXTIME",5.0,"maximum time",&tdyn);

  // Iterationparameters
  DoubleParameter("TOLTEMP",1.0E-10,
                  "tolerance in the temperature norm of the Newton iteration",
                  &tdyn);
  setStringToIntegralParameter<int>("NORM_TEMP","Abs","type of norm for temperature convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<int>(
                                 INPAR::THR::convnorm_abs,
                                 INPAR::THR::convnorm_rel,
                                 INPAR::THR::convnorm_mix),
                               &tdyn);

  DoubleParameter("TOLRES",1.0E-08,
                  "tolerance in the residual norm for the Newton iteration",
                  &tdyn);
  setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for residual convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<int>(
                                 INPAR::THR::convnorm_abs,
                                 INPAR::THR::convnorm_rel,
                                 INPAR::THR::convnorm_mix),
                               &tdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFTEMP","And","binary operator to combine temperature and residual force values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 INPAR::THR::bop_and,
                                 INPAR::THR::bop_or),
                               &tdyn);

  IntParameter("MAXITER",50,
               "maximum number of iterations allowed for Newton-Raphson iteration before failure",
               &tdyn);
  IntParameter("MINITER",0,
               "minimum number of iterations to be done within Newton-Raphson loop",
               &tdyn);
  setStringToIntegralParameter<int>("ITERNORM","L2","type of norm to be applied to residuals",
                               tuple<std::string>(
                                 "L1",
                                 "L2",
                                 "Rms",
                                 "Inf"),
                               tuple<int>(
                                 INPAR::THR::norm_l1,
                                 INPAR::THR::norm_l2,
                                 INPAR::THR::norm_rms,
                                 INPAR::THR::norm_inf),
                               &tdyn);
  setStringToIntegralParameter<int>("DIVERCONT","No",
                                    "Go on with time integration even if Newton-Raphson iteration failed",
                                     yesnotuple,yesnovalue,&tdyn);

  setStringToIntegralParameter<int>("NLNSOL","fullnewton","Nonlinear solution technique",
                               tuple<std::string>(
                                 "vague",
                                 "fullnewton"),
                                 tuple<int>(
                                     INPAR::THR::soltech_vague,
                                     INPAR::THR::soltech_newtonfull),
                                     &tdyn);

  setStringToIntegralParameter<int>("PREDICT","ConstTemp","Predictor of iterative solution techniques",
                               tuple<std::string>(
                                 "Vague",
                                 "ConstTemp",
                                 "ConstTempRate",
                                 "TangTemp"),
                               tuple<int>(
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

  setStringToIntegralParameter<int>("LUMPCAPA","No",
                               "Lump the capacity matrix for explicit time integration",
                               yesnotuple,yesnovalue,&tdyn);

  // number of linear solver used for thermal problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for thermal problems",&tdyn);

  /*----------------------------------------------------------------------*/
  /* parameters for generalised-alpha thermal integrator */
  Teuchos::ParameterList& tgenalpha = tdyn.sublist("GENALPHA",false,"");

  setStringToIntegralParameter<int>("GENAVG","ImrLike",
                              "mid-average type of internal forces",
                              tuple<std::string>(
                                "Vague",
                                "ImrLike",
                                "TrLike"),
                              tuple<int>(
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
  Teuchos::ParameterList& tsidyn = list->sublist(
   "TSI DYNAMIC",false,
   "Thermo Structure Interaction\n"
   "Dynamic section for TSI solver with various coupling methods"
   );

  // Coupling strategy for (partitioned and monolithic) TSI solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","tsi_monolithic",
                              "Coupling strategies for TSI solvers",
                              tuple<std::string>(
                                "tsi_oneway",
                                "tsi_sequstagg",
                                "tsi_iterstagg",
                                "tsi_iterstagg_aitken",
                                "tsi_iterstagg_aitkenirons",
                                "tsi_monolithic"
                                ),
                              tuple<int>(
                                INPAR::TSI::OneWay,
                                INPAR::TSI::SequStagg,
                                INPAR::TSI::IterStagg,
                                INPAR::TSI::IterStaggAitken,
                                INPAR::TSI::IterStaggAitkenIrons,
                                INPAR::TSI::Monolithic
                                ),
                              &tsidyn);

  // decide in partitioned TSI which one-way coupling should be used
  setStringToIntegralParameter<int>("COUPVARIABLE","Displacement",
                              "Coupling variable",
                              tuple<std::string>(
                                "Displacement",
                                "Temperature"
                                ),
                              tuple<int>(0,1),
                              &tsidyn);

  // Output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&tsidyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&tsidyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&tsidyn);
  DoubleParameter("TIMESTEP",0.05,"time step size dt",&tsidyn);
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check",&tsidyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&tsidyn);
  IntParameter("ITEMIN",1,"minimal number of iterations over fields",&tsidyn);
  IntParameter("UPRES",1,"increment for writing solution",&tsidyn);

  // Iterationparameters
  setStringToIntegralParameter<int>("NORM_INC","Abs","type of norm for primary variables convergence check",
                               tuple<std::string>(
                                 "Abs"
                                 ),
                               tuple<int>(
                                 INPAR::TSI::convnorm_abs
                                 ),
                               &tsidyn);

  setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for residual convergence check",
                                 tuple<std::string>(
                                   "Abs",
                                   "Rel",
                                   "RelIter0",
                                   "Mix"
                                   ),
                                 tuple<int>(
                                   INPAR::TSI::convnorm_abs,
                                   INPAR::TSI::convnorm_rel,
                                   INPAR::TSI::convnorm_reliter0,
                                   INPAR::TSI::convnorm_mix
                                   ),
                                 &tsidyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                               tuple<std::string>("And"),
                               tuple<int>(
                                 INPAR::TSI::bop_and
                                 ),
                               &tsidyn);

  // number of linear solver used for monolithic TSI
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for monolithic TSI problems",&tsidyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& poroelastdyn = list->sublist(
   "POROELASTICITY DYNAMIC",false,
   "Poroelasticity"
   );

  // Coupling strategy for (monolithic) porous media solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","poroelast_monolithic",
                              "Coupling strategies for poroelasticity solvers",
                              tuple<std::string>(
                                 "poroelast_monolithic",
                                 "poroelast_monolithicstructuresplit",
                                 "poroelast_monolithicfluidsplit"
                                ),
                              tuple<int>(
                                INPAR::POROELAST::Monolithic,
                                INPAR::POROELAST::Monolithic_structuresplit,
                                INPAR::POROELAST::Monolithic_fluidsplit
                                ),
                              &poroelastdyn);

  // Output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&poroelastdyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&poroelastdyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&poroelastdyn);
  DoubleParameter("TIMESTEP",0.05,"time step size dt",&poroelastdyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&poroelastdyn);
  IntParameter("ITEMIN",1,"minimal number of iterations over fields",&poroelastdyn);
  IntParameter("UPRES",1,"increment for writing solution",&poroelastdyn);

  // Iterationparameters
  DoubleParameter("RESTOL",1e-8,"tolerance in the residual norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("INCTOL",1e-8,"tolerance in the increment norm for the Newton iteration",&poroelastdyn);

  setStringToIntegralParameter<int>("NORM_INC","Abs","type of norm for primary variables convergence check",
                               tuple<std::string>(
                                 "Abs"
                                 ),
                               tuple<int>(
                                 INPAR::POROELAST::convnorm_abs
                                 ),
                               &poroelastdyn);

  setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for residual convergence check",
                                 tuple<std::string>(
                                   "Abs"
                                   ),
                                 tuple<int>(
                                   INPAR::POROELAST::convnorm_abs
                                   ),
                                 &poroelastdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                               tuple<std::string>(
                                     "And",
                                     "Or"),
                                     tuple<int>(
                                       INPAR::POROELAST::bop_and,
                                       INPAR::POROELAST::bop_or),
                               &poroelastdyn);

  setStringToIntegralParameter<int>("SECONDORDER","No",
                               "Second order coupling at the interface.",
                               yesnotuple,yesnovalue,&poroelastdyn);

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for poroelasticity problems",&poroelastdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& poroscatradyn = list->sublist(
   "POROSCATRA CONTROL",false,
   "Control paramters for scatra porous media coupling"
   );

  // Output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&poroscatradyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&poroscatradyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&poroscatradyn);
  DoubleParameter("TIMESTEP",0.05,"time step size dt",&poroscatradyn);
  IntParameter("UPRES",1,"increment for writing solution",&poroscatradyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidyn = list->sublist(
   "SSI CONTROL",false,
   "Control paramters for scatra structure interaction"
   );

  // Output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&ssidyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&ssidyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&ssidyn);
  DoubleParameter("TIMESTEP",0.05,"time step size dt",&ssidyn);
  IntParameter("UPRES",1,"increment for writing solution",&ssidyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&ssidyn);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&ssidyn);

  // Coupling strategy for SSI solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","one_way",
                              "Coupling strategies for SSI solvers",
                              tuple<std::string>(
                                "one_way",
                                "two_way"
                                ),
                              tuple<int>(
                                INPAR::SSI::Part_OneWay,
                                INPAR::SSI::Part_TwoWay
                                ),
                              &ssidyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& flucthydro = list->sublist("FLUCTUATING HYDRODYNAMICS",false,"");
  DoubleParameter("TEMPERATURE",300,"Temperature in K",&flucthydro);
  DoubleParameter("BOLTZMANNCONST",1.380650424e-23,"Boltzmann constant",&flucthydro);
  setStringToIntegralParameter<int>("SEEDCONTROL","No",
                                      "control seeding with given unsigned integer",
                                      yesnotuple,yesnovalue,&flucthydro);
  IntParameter("SEEDVARIABLE",0,"seed variable",&flucthydro);
  IntParameter("SAMPLEPERIOD",1,"sample period",&flucthydro);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn = list->sublist("FLUID DYNAMIC",false,"");

  // physical type of fluid flow (incompressible, varying density, loma, Boussinesq approximation)
  setStringToIntegralParameter<int>("PHYSICAL_TYPE","Incompressible",
                               "Physical Type",
                               tuple<string>(
                                 "Incompressible",
                                 "Varying_density",
                                 "Loma",
                                 "Boussinesq",
                                 "Poro",
                                 "Topology_optimization"
                                 ),
                               tuple<int>(
                                     INPAR::FLUID::incompressible,
                                     INPAR::FLUID::varying_density,
                                     INPAR::FLUID::loma,
                                     INPAR::FLUID::boussinesq,
                                     INPAR::FLUID::poro,
                                     INPAR::FLUID::topopt),
                               &fdyn);

  // number of linear solver used for fluid problem
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for fluid dynamics",&fdyn);

  // number of linear solver used for fluid problem (former fluid pressure solver for SIMPLER preconditioning with fluid)
  IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for fluid dynamics (fluid pressure solver within SIMPLER preconditioner)",&fdyn);

  setStringToIntegralParameter<int>("TIMEINTEGR","One_Step_Theta",
      "Time Integration Scheme",
      tuple<std::string>(
          "Stationary",
          "Np_Gen_Alpha",
          "Gen_Alpha",
          "Af_Gen_Alpha",
          "One_Step_Theta",
          "BDF2"),
      tuple<int>(
          INPAR::FLUID::timeint_stationary,
          INPAR::FLUID::timeint_npgenalpha,        // fluid3's implementation (solving for velocity increment)
          INPAR::FLUID::timeint_gen_alpha,         // Peter's implementation (solving for acceleration increment)
          INPAR::FLUID::timeint_afgenalpha,
          INPAR::FLUID::timeint_one_step_theta,
          INPAR::FLUID::timeint_bdf2),
          &fdyn);

  setStringToIntegralParameter<int>("STARTINGALGO","One_Step_Theta","",
                               tuple<std::string>(
                                 "One_Step_Theta"
                                 ),
                               tuple<int>(
                                 INPAR::FLUID::timeint_one_step_theta
                                 ),
                               &fdyn);
  setStringToIntegralParameter<int>("NONLINITER","fixed_point_like",
                               "Nonlinear iteration scheme",
                               tuple<std::string>(
                                 "fixed_point_like",
                                 "Newton",
                                 "minimal"
                                 ),
                               tuple<int>(
                                     INPAR::FLUID::fixed_point_like,
                                     INPAR::FLUID::Newton,
                                     INPAR::FLUID::minimal),
                               &fdyn);

  setStringToIntegralParameter<int>("PREDICTOR","default",
                                    "Predictor for first guess in nonlinear iteration",
                                    tuple<std::string>(
                                      "disabled",
                                      "default",
                                      "steady_state_predictor",
                                      "zero_acceleration_predictor",
                                      "constant_acceleration_predictor",
                                      "constant_increment_predictor",
                                      "explicit_second_order_midpoint"
                                      ),
                                    tuple<int>(1,2,3,4,5,6,7),
                                    &fdyn);

  setStringToIntegralParameter<int>("CONVCHECK","L_2_norm",
                               "norm for convergence check",
                               tuple<std::string>(
                                 //"L_infinity_norm",
                                 //"L_1_norm",
                                 "L_2_norm",
                                 "L_2_norm_without_residual_at_itemax"
                                 ),
                               tuple<std::string>(
                                 //"use max norm (ccarat)",
                                 //"use abs. norm (ccarat)",
                                 "compute L2 errors of increments (relative) and residuals (absolute)",
                                 "same as L_2_norm, only no residual norm is computed if itemax is reached (speedup for turbulence calculations, startup phase)"
                                 ),
                               tuple<int>(
                                 //INPAR::FLUID::fncc_Linf,
                                 //INPAR::FLUID::fncc_L1,
                                 INPAR::FLUID::fncc_L2,
                                 INPAR::FLUID::fncc_L2_wo_res
                                 ),
                               &fdyn);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial field for fluid problem",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function",
                                 "disturbed_field_from_function",
                                 "FLAME_VORTEX_INTERACTION",
                                 "BELTRAMI-FLOW",
                                 "KIM-MOIN-FLOW",
                                 "BOCHEV-TEST"),
                               tuple<int>(
                                     INPAR::FLUID::initfield_zero_field,
                                     INPAR::FLUID::initfield_field_by_function,
                                     INPAR::FLUID::initfield_disturbed_field_from_function,
                                     INPAR::FLUID::initfield_flame_vortex_interaction,
                                     INPAR::FLUID::initfield_beltrami_flow,
                                     INPAR::FLUID::initfield_kim_moin_flow,
                                     INPAR::FLUID::initfield_bochev_test),
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
                                 INPAR::FLUID::liftdrag_none,
                                 INPAR::FLUID::liftdrag_none,
                                 INPAR::FLUID::liftdrag_nodeforce,
                                 INPAR::FLUID::liftdrag_nodeforce,
                                 INPAR::FLUID::liftdrag_nodeforce,
                                 INPAR::FLUID::liftdrag_nodeforce,
                                 INPAR::FLUID::liftdrag_nodeforce
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

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
                                  tuple<std::string>(
                                    "no",
                                    "Condensed_Smat",
                                    "Condensed_Bmat",
                                    "Condensed_Bmat_merged",
                                    "SaddlePointSystem_coupled",
                                    "SaddlePointSystem_pc",
                                    "Coupling_ionTransport_Laplace"),
                                  tuple<int>(
                                      INPAR::FLUID::no_meshtying,
                                      INPAR::FLUID::condensed_smat,
                                      INPAR::FLUID::condensed_bmat,
                                      INPAR::FLUID::condensed_bmat_merged,
                                      INPAR::FLUID::sps_coupled,
                                      INPAR::FLUID::sps_pc,
                                      INPAR::FLUID::coupling_iontransport_laplace),
                                  &fdyn);

  setStringToIntegralParameter<int>("CALCERROR",
                               "no",
                               "Flag to (de)activate error calculation",
                               tuple<std::string>(
                                 "no",
                                 "beltrami_flow",
                                 "channel2D",
                                 "gravitation",
                                 "shear_flow",
                                 "jeffery_hamel_flow"),
                               tuple<int>(
                                   INPAR::FLUID::no_error_calculation,
                                   INPAR::FLUID::beltrami_flow,
                                   INPAR::FLUID::channel2D,
                                   INPAR::FLUID::gravitation,
                                   INPAR::FLUID::shear_flow,
                                   INPAR::FLUID::jeffery_hamel_flow),
                               &fdyn);

  setStringToIntegralParameter<int>("SIMPLER","no",
                               "Switch on SIMPLE family of solvers, needs additional FLUID PRESSURE SOLVER block!",
                               yesnotuple,yesnovalue,&fdyn);

/*  setStringToIntegralParameter<int>("SPLITFLUID","no",
                               "If yes, the fluid matrix is splitted into a block sparse matrix for velocity and pressure degrees of freedom (similar to SIMPLER flag)",
                               yesnotuple,yesnovalue,&fdyn);*/

  setStringToIntegralParameter<int>("ADAPTCONV","yes",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&fdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&fdyn);

  setStringToIntegralParameter<int>("INFNORMSCALING","no",
                               "Scale blocks of matrix with row infnorm?",
                               yesnotuple,yesnovalue,&fdyn);

  BoolParameter("GMSH_OUTPUT","No","write output to gmsh files",&fdyn);
  IntParameter("UPRES",1,"Increment for writing solution",&fdyn);
  IntParameter("RESTARTEVRY",20,"Increment for writing restart",&fdyn);
  IntParameter("NUMSTEP",1,"Total number of Timesteps",&fdyn);
  IntParameter("STEADYSTEP",-1,"steady state check every step",&fdyn);
  IntParameter("NUMSTASTEPS",0,"Number of Steps for Starting Scheme",&fdyn);
  IntParameter("STARTFUNCNO",-1,"Function for Initial Starting Field",&fdyn);
  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&fdyn);
  IntParameter("INITSTATITEMAX",5,"max number of nonlinear iterations for initial stationary solution",&fdyn);
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

  setStringToIntegralParameter<int>("STRONG_REDD_3D_COUPLING_TYPE",
                               "no",
                               "Flag to (de)activate potential Strong 3D redD coupling",
                               tuple<std::string>(
                                 "no",
                                 "yes"),
                               tuple<std::string>(
                                 "Weak coupling",
                                 "Strong coupling"),
                               tuple<int>(0,1),
                               &fdyn);

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

  // number of linear solver used for arterial dynamics
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for arterial dynamics",&andyn);

 /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& redawdyn = list->sublist("REDUCED DIMENSIONAL AIRWAYS DYNAMIC",false,"");

  setStringToIntegralParameter<int>("DYNAMICTYP","CrankNicolson",
                               "CrankNicolson Scheme",
                               tuple<std::string>(
                                 "CrankNicolson"
                                 ),
                               tuple<int>(
                                typ_crank_nicolson
                                ),
                               &redawdyn);


  setStringToIntegralParameter<int>("SOLVERTYPE","Linear",
                               "Solver type",
                               tuple<std::string>(
                                 "Linear",
                                 "Nonlinear"
                                 ),
                               tuple<int>(
                                 linear,
                                 nonlinear
                                ),
                               &redawdyn);

  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&redawdyn);
  IntParameter("NUMSTEP",0,"Number of Time Steps",&redawdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&redawdyn);
  IntParameter("UPRES",1,"Increment for writing solution",&redawdyn);

  IntParameter("MAXITERATIONS",1,"maximum iteration steps",&redawdyn);
  DoubleParameter("TOLERANCE",1.0E-6,"tolerance",&redawdyn);

  // number of linear solver used for reduced dimensional airways dynamic
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for reduced dim arterial dynamics",&redawdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_stab = fdyn.sublist("STABILIZATION",false,"");

  // this parameter seperates stabilized from unstabilized methods
  setStringToIntegralParameter<int>("STABTYPE",
                               "residual_based",
                               "Apply (un)stabilized fluid formulation",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "inconsistent",
                                 "residual_based",
                                 "edge_based"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> inf-sup stable elements required!",
                                 "Similar to residual based without second derivatives (i.e. only consistent for tau->0, but faster)",
                                 "Use a residual-based stabilization or, more generally, a stabilization \nbased on the concept of the residual-based variational multiscale method...\nExpecting additional input",
                                 "Use an edge-based stabilization, especially for XFEM")  ,
                               tuple<int>(0,1,2,3),
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
                                 tuple<int>(
                                   INPAR::FLUID::subscales_quasistatic,
                                   INPAR::FLUID::subscales_time_dependent),
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
                               tuple<int>(
                                   INPAR::FLUID::inertia_stab_drop,
                                   INPAR::FLUID::inertia_stab_keep,
                                   INPAR::FLUID::inertia_stab_keep_complete),
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
                               tuple<int>(
                                 INPAR::FLUID::pstab_assume_inf_sup_stable,   //No PSPG -> inf-sup-stable elements mandatory
                                 INPAR::FLUID::pstab_use_pspg),               // Use PSPG -> allowing for equal-order interpolation
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
                               tuple<int>(
                                 INPAR::FLUID::convective_stab_none,  // no SUPG stabilization
                                 INPAR::FLUID::convective_stab_supg), // use SUPG stabilization
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
                                 tuple<int>(
                                     INPAR::FLUID::viscous_stab_none,
                                     INPAR::FLUID::viscous_stab_gls,
                                     INPAR::FLUID::viscous_stab_gls_only_rhs,
                                     INPAR::FLUID::viscous_stab_usfem,
                                     INPAR::FLUID::viscous_stab_usfem_only_rhs
                                   ),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("RSTAB",
                               "no_rstab",
                               "Flag to (de)activate reactive term in residual-based stabilization.",
                               tuple<std::string>(
                                 "no_rstab",
                                 "rstab_gls",
                                 "rstab_usfem"
                                 ),
                               tuple<std::string>(
                                 "no reactive term in stabilization",
                                 "reactive stabilization of GLS type",
                                 "reactive stabilization of USFEM type"
                                 ),
                                 tuple<int>(
                                     INPAR::FLUID::reactive_stab_none,
                                     INPAR::FLUID::reactive_stab_gls,
                                     INPAR::FLUID::reactive_stab_usfem
                                   ),
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
                               tuple<int>(
                                   INPAR::FLUID::continuity_stab_none,
                                   INPAR::FLUID::continuity_stab_yes
                                 ),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("CROSS-STRESS",
                               "no_cross",
                               "Flag to (de)activate cross-stress term -> residual-based VMM.",
                               tuple<std::string>(
                                 "no_cross",
                                 "yes_cross",
                                 "cross_rhs"
                                 //"cross_complete"
                                 ),
                               tuple<std::string>(
                                 "No cross-stress term",
                                 "Include the cross-stress term with a linearization of the convective part",
                                 "Include cross-stress term, but only explicitly on right hand side"
                                 //""
                                 ),
                               tuple<int>(
                                   INPAR::FLUID::cross_stress_stab_none,
                                   INPAR::FLUID::cross_stress_stab,
                                   INPAR::FLUID::cross_stress_stab_only_rhs
                                ),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("REYNOLDS-STRESS",
                               "no_reynolds",
                               "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
                               tuple<std::string>(
                                 "no_reynolds",
                                 "yes_reynolds",
                                 "reynolds_rhs"
                                 //"reynolds_complete"
                                 ),
                               tuple<std::string>(
                                 "No Reynolds-stress term",
                                 "Include Reynolds-stress term with linearisation",
                                 "Include Reynolds-stress term explicitly on right hand side"
                                 //""
                                 ),
                               tuple<int>(
                                   INPAR::FLUID::reynolds_stress_stab_none,
                                   INPAR::FLUID::reynolds_stress_stab,
                                   INPAR::FLUID::reynolds_stress_stab_only_rhs
                               ),
                               &fdyn_stab);

  // this parameter selects the tau definition applied
  // the options "tau_not_defined" are only available in Peter's Genalpha code
  // (there, which_tau is set via the string, not via INPAR::FLUID::TauType)
  setStringToIntegralParameter<int>("DEFINITION_TAU",
                               "Franca_Barrenechea_Valentin_Frey_Wall",
                               "Definition of tau_M and Tau_C",
                               tuple<std::string>(
                                 "Taylor_Hughes_Zarins",
                                 "Taylor_Hughes_Zarins_wo_dt",
                                 "Taylor_Hughes_Zarins_Whiting_Jansen",
                                 "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt",
                                 "Taylor_Hughes_Zarins_scaled",
                                 "Taylor_Hughes_Zarins_scaled_wo_dt",
                                 "Franca_Barrenechea_Valentin_Frey_Wall",
                                 "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt",
                                 "Shakib_Hughes_Codina",
                                 "Shakib_Hughes_Codina_wo_dt",
                                 "Codina",
                                 "Codina_wo_dt",
                                 "Franca_Madureira_Valentin_Badia_Codina",
                                 "Franca_Madureira_Valentin_Badia_Codina_wo_dt",
                                 //"BFVW_gradient_based_hk",
                                 "Smoothed_FBVW"),
                               tuple<int>(
                                   INPAR::FLUID::tau_taylor_hughes_zarins,
                                   INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt,
                                   INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen,
                                   INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt,
                                   INPAR::FLUID::tau_taylor_hughes_zarins_scaled,
                                   INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt,
                                   INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall,
                                   INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt,
                                   INPAR::FLUID::tau_shakib_hughes_codina,
                                   INPAR::FLUID::tau_shakib_hughes_codina_wo_dt,
                                   INPAR::FLUID::tau_codina,
                                   INPAR::FLUID::tau_codina_wo_dt,
                                   INPAR::FLUID::tau_franca_madureira_valentin_badia_codina,
                                   INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt,
                                   //INPAR::FLUID::tau_not_defined,
                                   INPAR::FLUID::tau_not_defined),
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

  // these parameters active additional terms in loma continuity equation
  // which might be identified as SUPG-/cross- and Reynolds-stress term
  setStringToIntegralParameter<int>("LOMA_CONTI_SUPG",
                               "no_supg",
                               "Flag to (de)activate SUPG-term in loma continuity equation.",
                               tuple<std::string>(
                                 "no_supg",
                                 "yes_supg"),
                               tuple<std::string>(
                                 "No SUPG",
                                 "Use SUPG."),
                               tuple<int>(
                                 INPAR::FLUID::convective_stab_none,  // no SUPG stabilization
                                 INPAR::FLUID::convective_stab_supg), // use SUPG stabilization
                               &fdyn_stab);

  setStringToIntegralParameter<int>("LOMA_CONTI_CROSS_STRESS",
                               "no_cross",
                               "Flag to (de)activate cross-stress term loma continuity equation-> residual-based VMM.",
                               tuple<std::string>(
                                 "no_cross",
                                 "yes_cross",
                                 "cross_rhs"
                                 //"cross_complete"
                                 ),
                               tuple<std::string>(
                                 "No cross-stress term",
                                 "Include the cross-stress term with a linearization of the convective part",
                                 "Include cross-stress term, but only explicitly on right hand side"
                                 //""
                                 ),
                               tuple<int>(
                                   INPAR::FLUID::cross_stress_stab_none,
                                   INPAR::FLUID::cross_stress_stab,
                                   INPAR::FLUID::cross_stress_stab_only_rhs
                                ),
                               &fdyn_stab);

  setStringToIntegralParameter<int>("LOMA_CONTI_REYNOLDS_STRESS",
                               "no_reynolds",
                               "Flag to (de)activate Reynolds-stress term loma continuity equation-> residual-based VMM.",
                               tuple<std::string>(
                                 "no_reynolds",
                                 "yes_reynolds",
                                 "reynolds_rhs"
                                 //"reynolds_complete"
                                 ),
                               tuple<std::string>(
                                 "No Reynolds-stress term",
                                 "Include Reynolds-stress term with linearisation",
                                 "Include Reynolds-stress term explicitly on right hand side"
                                 //""
                                 ),
                               tuple<int>(
                                   INPAR::FLUID::reynolds_stress_stab_none,
                                   INPAR::FLUID::reynolds_stress_stab,
                                   INPAR::FLUID::reynolds_stress_stab_only_rhs
                               ),
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
      "Dynamic_Smagorinsky",
      "Scale_Similarity",
      "Scale_Similarity_basic",
      "Multifractal_Subgrid_Scales"),
    tuple<std::string>(
      "If classical LES is our turbulence approach, this is a contradiction and should cause a dserror.",
      "Classical constant coefficient Smagorinsky model. Be careful if you \nhave a wall bounded flow domain!",
      "Use an exponential damping function for the turbulent viscosity \nclose to the wall. This is only implemented for a channel geometry of \nheight 2 in y direction. The viscous lengthscale l_tau is \nrequired as additional input.",
      "The solution is filtered and by comparison of the filtered \nvelocity field with the real solution, the Smagorinsky constant is \nestimated in each step --- mind that this procedure includes \nan averaging in the xz plane, hence this implementation will only work \nfor a channel flow.",
      "Scale Similarity Model coherent with the variational multiscale formulation",
      "Scale Similarity Model according to liu, meneveau, katz",
      "Multifractal Subgrid-Scale Modeling based on the work of burton"),
    tuple<int>(0,1,2,3,4,5,6),
    &fdyn_turbu);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    // Otherwise BACI DEBUG version will crash during runtime!
    Teuchos::Tuple<std::string,15> name;
    Teuchos::Tuple<int,15> label;
    name[ 0] = "no";                                      label[ 0] = 0;
    name[ 1] = "time_averaging";                          label[ 1] = 1;
    name[ 2] = "channel_flow_of_height_2";                label[ 2] = 2;
    name[ 3] = "lid_driven_cavity";                       label[ 3] = 3;
    name[ 4] = "backward_facing_step";                    label[ 4] = 4;
    name[ 5] = "square_cylinder";                         label[ 5] = 5;
    name[ 6] = "square_cylinder_nurbs";                   label[ 6] = 6;
    name[ 7] = "rotating_circular_cylinder_nurbs";        label[ 7] = 7;
    name[ 8] = "rotating_circular_cylinder_nurbs_scatra"; label[ 8] = 8;
    name[ 9] = "loma_channel_flow_of_height_2";           label[ 9] = 9;
    name[10] = "loma_lid_driven_cavity";                  label[10] = 10;
    name[11] = "loma_backward_facing_step";               label[11] = 11;
    name[12] = "combust_oracles";                         label[12] = 12;
    name[13] = "bubbly_channel_flow";                     label[13] = 13;
    name[14] = "scatra_channel_flow_of_height_2";         label[14] = 14;

    Teuchos::Tuple<std::string,15> description;
    description[0]="The flow is not further specified, so spatial averaging \nand hence the standard sampling procedure is not possible";
    description[1]="The flow is not further specified, but time averaging of velocity and pressure field is performed";
    description[2]="For this flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.";
    description[3]="For this flow, all statistical data are evaluated on the center lines of the xy-midplane, averaged only over time.";
    description[4]="For this flow, statistical data are evaluated on various lines, averaged over time and z.";
    description[5]="For this flow, statistical data are evaluated on various lines of the xy-midplane, averaged only over time.";
    description[6]="For this flow, statistical data are evaluated on various lines of the xy-midplane, averaged over time and eventually in one hom.direction.";
    description[7]="For this flow, statistical data is computed in concentric surfaces and averaged. in time and in one hom. direction";
    description[8]="For this flow with mass transport, statistical data is computed in concentric surfaces and averaged. in time and in one hom. direction";
    description[9]="For this low-Mach-number flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.";
    description[10]="For this low-Mach-number flow, all statistical data are evaluated on the center lines of the xy-midplane, averaged only over time.";
    description[11]="For this low-Mach-number flow, statistical data are evaluated on various lines, averaged over time and z.";
    description[12]="ORACLES test rig for turbulent premixed combustion.";
    description[13]="Turbulent two-phase flow: bubbly channel flow, statistical data are averaged in homogeneous planse and over time.";
    description[14]="For this flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.";

    setStringToIntegralParameter<int>(
        "CANONICAL_FLOW",
        "no",
        "Sampling is different for different canonical flows \n--- so specify what kind of flow you've got",
        name,
        description,
        label,
        &fdyn_turbu);
  }

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
    "CHAN_AMPL_INIT_DIST",
    0.1,
    "Max. amplitude of the random disturbance in percent of the initial value in mean flow direction.",
    &fdyn_turbu);

  IntParameter("SAMPLING_START",10000000,"Time step after when sampling shall be started",&fdyn_turbu);
  IntParameter("SAMPLING_STOP",1,"Time step when sampling shall be stopped",&fdyn_turbu);
  IntParameter("DUMPING_PERIOD",1,"Period of time steps after which statistical data shall be dumped",&fdyn_turbu);
  BoolParameter("SUBGRID_DISSIPATION","No","Flag to (de)activate estimation of subgrid-scale dissipation.",&fdyn_turbu);

  setStringToIntegralParameter<int>("FSSUGRVISC","No","fine-scale subgrid viscosity",
                               tuple<std::string>(
                                 "No",
                                 "Smagorinsky_all",
                                 "Smagorinsky_small"
                                 ),
                               tuple<int>(
                                   INPAR::FLUID::no_fssgv,
                                   INPAR::FLUID::smagorinsky_all,
                                   INPAR::FLUID::smagorinsky_small
                                   ),
                               &fdyn_turbu);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for smagorinsky model
  Teuchos::ParameterList& fdyn_turbsgv = fdyn.sublist("SUBGRID VISCOSITY",false,"");

  DoubleParameter("C_SMAGORINSKY",0.0,"Constant for the Smagorinsky model. Something between 0.1 to 0.24",&fdyn_turbsgv);
  BoolParameter("C_SMAGORINSKY_AVERAGED","No","Flag to (de)activate averaged Smagorinksy constant",&fdyn_turbsgv);

  DoubleParameter(
    "CHANNEL_L_TAU",
    0.0,
    "Used for normalisation of the wall normal distance in the Van \nDriest Damping function. May be taken from the output of \nthe apply_mesh_stretching.pl preprocessing script.",
    &fdyn_turbsgv);

  DoubleParameter("C_TURBPRANDTL",1.0,"(Constant) turbulent Prandtl number for the Smagorinsky model in scalar transport.",&fdyn_turbsgv);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for multifractal subgrid-scales
  Teuchos::ParameterList& fdyn_turbmfs = fdyn.sublist("MULTIFRACTAL SUBGRID SCALES",false,"");

  DoubleParameter(
    "CSGS",
    0.0,
    "Modelparameter of multifractal subgrid-scales.",
    &fdyn_turbmfs);

  setStringToIntegralParameter<int>(
    "SCALE_SEPARATION",
    "no_scale_sep",
    "Specify the filter type for scale separation in LES",
    tuple<std::string>(
      "no_scale_sep",
      "box_filter",
      "algebraic_multigrid_operator",
      "geometric_multigrid_operator"),
    tuple<std::string>(
      "no scale separation",
      "classical box filter",
      "scale separation by algebraic multigrid operator",
      "scale separation by geometric multigrid operator"),
    tuple<int>(0,1,2,3),
    &fdyn_turbmfs);

  BoolParameter("CALC_N","No","Flag to (de)activate calculation of N from the Reynolds number.",&fdyn_turbmfs);

  DoubleParameter(
    "N",
    1.0,
    "Set grid to viscous scale ratio..",
    &fdyn_turbmfs);

  setStringToIntegralParameter<int>(
    "REF_LENGTH",
    "cube_edge",
    "Specify the reference length for Re-dependent N.",
    tuple<std::string>(
      "cube_edge",
      "sphere_diameter",
      "streamlength",
      "gradient_based",
      "metric_tensor"),
    tuple<std::string>(
      "edge length of volume equivalent cube",
      "diameter of volume equivalent sphere",
      "streamlength taken from stabilization",
      "gradient based length taken from stabilization",
      "metric tensor taken from stabilization"),
    tuple<int>(0,1,2,3,4),
    &fdyn_turbmfs);

  setStringToIntegralParameter<int>(
    "REF_VELOCITY",
    "strainrate",
    "Specify the reference velocity for Re-dependent N.",
    tuple<std::string>(
      "strainrate",
      "resolved",
      "fine_scale"),
    tuple<std::string>(
      "norm of strain rate",
      "resolved velocity",
      "fine-scale velocity"),
    tuple<int>(0,1,2),
    &fdyn_turbmfs);

  DoubleParameter(
    "C_NU",
    1.0,
    "Proportionality constant between Re and ratio viscous scale to element length.",
    &fdyn_turbmfs);

  BoolParameter("NEAR_WALL_LIMIT","No","Flag to (de)activate near-wall limit.",&fdyn_turbmfs);

  setStringToIntegralParameter<int>("EVALUATION_B",
                               "element_center",
                               "Location where B is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate B at element center",
                                 "evaluate B at integration point")  ,
                                tuple<int>(0,1),
                               &fdyn_turbmfs);

  DoubleParameter(
    "BETA",
    0.0,
    "Cross- and Reynolds-stress terms only on right-hand-side.",
    &fdyn_turbmfs);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(0,1),
                               &fdyn_turbmfs);

  DoubleParameter("C_SCALE_SIMILARITY",1.0,"Constant for the scale similarity model. Something between 0.45 +- 0.15 or 1.0.", &fdyn_turbmfs);

  DoubleParameter(
    "CSGS_PHI",
    0.0,
    "Modelparameter of multifractal subgrid-scales for scalar transport.",
    &fdyn_turbmfs);

  DoubleParameter(
    "C_DIFF",
    1.0,
    "Proportionality constant between Re*Pr and ratio dissipative scale to element length. Usually equal cnu.",
    &fdyn_turbmfs);

  // activate cross- and Reynolds-stress terms in loma continuity equation
  BoolParameter("LOMA_CONTI","No","Flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation.",&fdyn_turbmfs);

  /*----------------------------------------------------------------------*/
   Teuchos::ParameterList& fdyn_turbinf = fdyn.sublist("TURBULENT INFLOW",false,"");

   BoolParameter("TURBULENTINFLOW","No","Flag to (de)activate potential separate turbulent inflow section",&fdyn_turbinf);

   setStringToIntegralParameter<int>("INITIALINFLOWFIELD","zero_field",
                                "Initial field for inflow section",
                                tuple<std::string>(
                                  "zero_field",
                                  "field_by_function",
                                  "disturbed_field_from_function"),
                                tuple<int>(
                                      INPAR::FLUID::initfield_zero_field,
                                      INPAR::FLUID::initfield_field_by_function,
                                      INPAR::FLUID::initfield_disturbed_field_from_function),
                                &fdyn_turbinf);

   IntParameter("INFLOWFUNC",-1,"Function number for initial flow field in inflow section",&fdyn_turbinf);

   DoubleParameter(
     "INFLOW_INIT_DIST",
     0.1,
     "Max. amplitude of the random disturbance in percent of the initial value in mean flow direction.",
     &fdyn_turbinf);

   IntParameter("NUMINFLOWSTEP",1,"Total number of time steps for development of turbulent flow",&fdyn_turbinf);

   setStringToIntegralParameter<int>(
       "CANONICAL_INFLOW",
       "no",
       "Sampling is different for different canonical flows \n--- so specify what kind of flow you've got",
       tuple<std::string>(
       "no",
       "time_averaging",
       "channel_flow_of_height_2",
       "loma_channel_flow_of_height_2",
       "scatra_channel_flow_of_height_2"),
       tuple<std::string>(
       "The flow is not further specified, so spatial averaging \nand hence the standard sampling procedure is not possible",
       "The flow is not further specified, but time averaging of velocity and pressure field is performed",
       "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.",
       "For this low-Mach-number flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow.",
       "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it is essentially a statistically one dimensional flow."),
       tuple<int>(0,1,2,3,4),
       &fdyn_turbinf);

   DoubleParameter(
     "INFLOW_CHA_SIDE",
     0.0,
     "Most right side of inflow channel. Necessary to define sampling domain.",
     &fdyn_turbinf);

   setStringToIntegralParameter<int>(
     "INFLOW_HOMDIR",
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
     &fdyn_turbinf);

   IntParameter("INFLOW_SAMPLING_START",10000000,"Time step after when sampling shall be started",&fdyn_turbinf);
   IntParameter("INFLOW_SAMPLING_STOP",1,"Time step when sampling shall be stopped",&fdyn_turbinf);
   IntParameter("INFLOW_DUMPING_PERIOD",1,"Period of time steps after which statistical data shall be dumped",&fdyn_turbinf);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC",false,"");

  DoubleParameter("TIMESTEP",0.1,"",&adyn);
  IntParameter("NUMSTEP",41,"",&adyn);
  DoubleParameter("MAXTIME",4.0,"",&adyn);
  setStringToIntegralParameter<int>("ALE_TYPE","classic_lin","ale mesh movement algorithm",
                               tuple<std::string>("classic_lin","incr_lin","laplace","springs","springs_fixed_ref"),
                               tuple<int>(INPAR::ALE::classic_lin,
                                   INPAR::ALE::incr_lin,
                                   INPAR::ALE::laplace,
                                   INPAR::ALE::springs,
                                   INPAR::ALE::springs_fixed_ref),
                               &adyn);
  IntParameter("NUM_INITSTEP",0,"",&adyn);
  IntParameter("RESULTSEVRY",1,"",&adyn);

  // linear solver id used for scalar ale problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for ale problems...",&adyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn = list->sublist(
      "SCALAR TRANSPORT DYNAMIC",
      false,
      "control parameters for scalar transport problems\n");

  setStringToIntegralParameter<int>("SOLVERTYPE","linear_full",
                               "type of scalar transport solver",
                               tuple<std::string>(
                                 "linear_full",
                                 "linear_incremental",
                                 "nonlinear"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::solvertype_linear_full,
                                   INPAR::SCATRA::solvertype_linear_incremental,
                                   INPAR::SCATRA::solvertype_nonlinear),
                               &scatradyn);

  setStringToIntegralParameter<int>("TIMEINTEGR","One_Step_Theta",
                               "Time Integration Scheme",
                               tuple<std::string>(
                                 "Stationary",
                                 "One_Step_Theta",
                                 "BDF2",
                                 "Gen_Alpha",
                                 "Taylor_Galerkin_2",
                                 "Taylor_Galerkin_2_LW",
                                 "Taylor_Galerkin_3",
                                 "Taylor_Galerkin_4_Leapfrog",
                                 "Taylor_Galerkin_4_onestep"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::timeint_stationary,
                                   INPAR::SCATRA::timeint_one_step_theta,
                                   INPAR::SCATRA::timeint_bdf2,
                                   INPAR::SCATRA::timeint_gen_alpha,
                                   INPAR::SCATRA::timeint_tg2,
                                   INPAR::SCATRA::timeint_tg2_LW,
                                   INPAR::SCATRA::timeint_tg3,                     //schott 05/11
                                   INPAR::SCATRA::timeint_tg4_leapfrog,
                                   INPAR::SCATRA::timeint_tg4_onestep
                                 ),
                               &scatradyn);

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&scatradyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&scatradyn);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&scatradyn);
  DoubleParameter("THETA",0.5,"One-step-theta time integration factor",&scatradyn);
  DoubleParameter("ALPHA_M",0.5,"Generalized-alpha time integration factor",&scatradyn);
  DoubleParameter("ALPHA_F",0.5,"Generalized-alpha time integration factor",&scatradyn);
  DoubleParameter("GAMMA",0.5,"Generalized-alpha time integration factor",&scatradyn);
  //IntParameter("WRITESOLEVRY",1,"Increment for writing solution",&scatradyn);
  IntParameter("UPRES",1,"Increment for writing solution",&scatradyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&scatradyn);
  IntParameter("MATID",-1,"Material ID for automatic mesh generation",&scatradyn);

  setStringToIntegralParameter<int>("VELOCITYFIELD","zero",
                               "type of velocity field used for scalar transport problems",
                               tuple<std::string>(
                                 "zero",
                                 "function",
                                 "function_and_curve",
                                 "Navier_Stokes"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::velocity_zero,
                                   INPAR::SCATRA::velocity_function,
                                   INPAR::SCATRA::velocity_function_and_curve,
                                   INPAR::SCATRA::velocity_Navier_Stokes),
                               &scatradyn);

  IntParameter("VELFUNCNO",-1,"function number for scalar transport velocity field",&scatradyn);

  IntParameter("VELCURVENO",-1,"curve number for time-dependent scalar transport velocity field",&scatradyn);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Field for scalar transport problem",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function",
                                 "field_by_condition",
                                 "disturbed_field_by_function",
                                 "1D_DISCONTPV",
                                 "FLAME_VORTEX_INTERACTION",
                                 "RAYTAYMIXFRAC",
                                 "L_shaped_domain",
                                 "facing_flame_fronts",
                                 "oracles_flame"),
                               tuple<int>(
                                   INPAR::SCATRA::initfield_zero_field,
                                   INPAR::SCATRA::initfield_field_by_function,
                                   INPAR::SCATRA::initfield_field_by_condition,
                                   INPAR::SCATRA::initfield_disturbed_field_by_function,
                                   INPAR::SCATRA::initfield_discontprogvar_1D,
                                   INPAR::SCATRA::initfield_flame_vortex_interaction,
                                   INPAR::SCATRA::initfield_raytaymixfrac,
                                   INPAR::SCATRA::initfield_Lshapeddomain,
                                   INPAR::SCATRA::initfield_facing_flame_fronts,
                                   INPAR::SCATRA::initfield_oracles_flame),
                               &scatradyn);

  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&scatradyn);

  setStringToIntegralParameter<int>("CALCERROR","No",
                               "compute error compared to analytical solution",
                               tuple<std::string>(
                                 "No",
                                 "Kwok_Wu",
                                 "ConcentricCylinders",
                                 "Electroneutrality"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::calcerror_no,
                                   INPAR::SCATRA::calcerror_Kwok_Wu,
                                   INPAR::SCATRA::calcerror_cylinder,
                                   INPAR::SCATRA::calcerror_electroneutrality
                                   ),
                               &scatradyn);

  setStringToIntegralParameter<int>("WRITEFLUX","No","output of diffusive/total flux vectors",
                               tuple<std::string>(
                                 "No",
                                 "totalflux_domain",
                                 "diffusiveflux_domain",
                                 "totalflux_boundary",
                                 "diffusiveflux_boundary",
                                 "convectiveflux_boundary"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::flux_no,
                                   INPAR::SCATRA::flux_total_domain,
                                   INPAR::SCATRA::flux_diffusive_domain,
                                   INPAR::SCATRA::flux_total_boundary,
                                   INPAR::SCATRA::flux_diffusive_boundary,
                                   INPAR::SCATRA::flux_convective_boundary),
                               &scatradyn);

  setNumericStringParameter("WRITEFLUX_IDS","-1",
      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
      &scatradyn);

  BoolParameter("OUTMEAN","No","Output of mean values for scalars and density",&scatradyn);
  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&scatradyn);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(
                                 INPAR::SCATRA::convform_convective,
                                 INPAR::SCATRA::convform_conservative),
                               &scatradyn);

  BoolParameter("NEUMANNINFLOW",
      "no","Flag to (de)activate potential Neumann inflow term(s)",&scatradyn);

  BoolParameter("CONV_HEAT_TRANS",
      "no","Flag to (de)activate potential convective heat transfer boundary conditions",&scatradyn);

  BoolParameter("SKIPINITDER",
      "no","Flag to skip computation of initial time derivative",&scatradyn);

  BoolParameter("INITPOTCALC","no",
      "Automatically calculate initial field for electric potential",&scatradyn);

  setStringToIntegralParameter<int>("FSSUGRDIFF",
                               "No",
                               "fine-scale subgrid diffusivity",
                               tuple<std::string>(
                                 "No",
                                 "artificial",
                                 "Smagorinsky_all",
                                 "Smagorinsky_small"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::fssugrdiff_no,
                                   INPAR::SCATRA::fssugrdiff_artificial,
                                   INPAR::SCATRA::fssugrdiff_smagorinsky_all,
                                   INPAR::SCATRA::fssugrdiff_smagorinsky_small),
                               &scatradyn);

  BoolParameter("BLOCKPRECOND","NO",
      "Switch to block-preconditioned family of solvers, needs additional SIMPLER SOLVER block!",&scatradyn);

  setStringToIntegralParameter<int>("SCATRATYPE","Undefined",
                               "Type of scalar transport problem",
                               tuple<std::string>(
                                 "Undefined",
                                 "ConvectionDiffusion",
                                 "LowMachNumberFlow",
                                 "Elch_ENC",
                                 "Elch_ENC_PDE",
                                 "Elch_ENC_PDE_ELIM",
                                 "Elch_Poisson",
                                 "Elch_Laplace",
                                 "LevelSet",
                                 "TurbulentPassiveScalar",
                                 "Poroscatra"),
                               tuple<int>(
                                 INPAR::SCATRA::scatratype_undefined,
                                 INPAR::SCATRA::scatratype_condif,
                                 INPAR::SCATRA::scatratype_loma,
                                 INPAR::SCATRA::scatratype_elch_enc,
                                 INPAR::SCATRA::scatratype_elch_enc_pde,
                                 INPAR::SCATRA::scatratype_elch_enc_pde_elim,
                                 INPAR::SCATRA::scatratype_elch_poisson,
                                 INPAR::SCATRA::scatratype_elch_laplace,
                                 INPAR::SCATRA::scatratype_levelset,
                                 INPAR::SCATRA::scatratype_turbpassivesca,
                                 INPAR::SCATRA::scatratype_poro),
                                 &scatradyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
                                  tuple<std::string>(
                                      "no",
                                      "Condensed_Smat",
                                      "Condensed_Bmat",
                                      "Condensed_Bmat_merged",
                                      "SaddlePointSystem_coupled",
                                      "SaddlePointSystem_pc",
                                      "Coupling_ionTransport_Laplace"), //use the condensed_bmat_merged strategy
                                    tuple<int>(
                                        INPAR::FLUID::no_meshtying,
                                        INPAR::FLUID::condensed_smat,
                                        INPAR::FLUID::condensed_bmat,
                                        INPAR::FLUID::condensed_bmat_merged,
                                        INPAR::FLUID::sps_coupled,
                                        INPAR::FLUID::sps_pc,
                                        INPAR::FLUID::coupling_iontransport_laplace),   //use the condensed_bmat_merged strategy
                                    &scatradyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for scalar transport/elch...",&scatradyn);
  IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with SIMPLER)...",&scatradyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatra_nonlin = scatradyn.sublist(
      "NONLINEAR",
      false,
      "control parameters for solving nonlinear SCATRA problems\n");

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&scatra_nonlin);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&scatra_nonlin);
  BoolParameter("EXPLPREDICT","no","do an explicit predictor step before starting nonlinear iteration",&scatra_nonlin);
  DoubleParameter("ABSTOLRES",1e-14,"Absolute tolerance for deciding if residual of nonlinear problem is already zero",&scatra_nonlin);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&scatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&scatra_nonlin);


/*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn_stab = scatradyn.sublist("STABILIZATION",false,"");

  // this parameter governs type of stabilization
  setStringToIntegralParameter<int>("STABTYPE",
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
                               tuple<int>(
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
  setStringToIntegralParameter<int>("DEFINITION_TAU",
                               "Franca_Valentin",
                               "Definition of tau",
                               tuple<std::string>(
                                 "Taylor_Hughes_Zarins",
                                 "Taylor_Hughes_Zarins_wo_dt",
                                 "Franca_Valentin",
                                 "Franca_Valentin_wo_dt",
                                 "Shakib_Hughes_Codina",
                                 "Shakib_Hughes_Codina_wo_dt",
                                 "Codina",
                                 "Codina_wo_dt",
                                 "Franca_Madureira_Valentin",
                                 "Franca_Madureira_Valentin_wo_dt",
                                 "Exact_1D",
                                 "Zero"),
                                tuple<int>(
                                    INPAR::SCATRA::tau_taylor_hughes_zarins,
                                    INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt,
                                    INPAR::SCATRA::tau_franca_valentin,
                                    INPAR::SCATRA::tau_franca_valentin_wo_dt,
                                    INPAR::SCATRA::tau_shakib_hughes_codina,
                                    INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt,
                                    INPAR::SCATRA::tau_codina,
                                    INPAR::SCATRA::tau_codina_wo_dt,
                                    INPAR::SCATRA::tau_franca_madureira_valentin,
                                    INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt,
                                    INPAR::SCATRA::tau_exact_1d,
                                    INPAR::SCATRA::tau_zero),
                               &scatradyn_stab);

  // this parameter selects the all-scale subgrid-diffusivity definition applied
  setStringToIntegralParameter<int>("DEFINITION_ASSGD",
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
                                tuple<int>(
                                    INPAR::SCATRA::assgd_artificial,
                                    INPAR::SCATRA::assgd_hughes,
                                    INPAR::SCATRA::assgd_tezduyar,
                                    INPAR::SCATRA::assgd_docarmo,
                                    INPAR::SCATRA::assgd_almeida),
                               &scatradyn_stab);

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
                                tuple<int>(
                                  INPAR::SCATRA::evaltau_element_center,
                                  INPAR::SCATRA::evaltau_integration_point),
                               &scatradyn_stab);

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
                                 "evaluate material law at integration point"),
                               tuple<int>(
                                 INPAR::SCATRA::evalmat_element_center,
                                 INPAR::SCATRA::evalmat_integration_point),
                               &scatradyn_stab);

  // this parameter selects methods for improving consistency of stabilization terms
  setStringToIntegralParameter<int>("CONSISTENCY",
                               "no",
                               "improvement of consistency for stabilization",
                               tuple<std::string>(
                                 "no",
                                 "L2_projection_lumped"),
                               tuple<std::string>(
                                 "inconsistent",
                                 "L2 projection with lumped mass matrix")  ,
                                tuple<int>(
                                  INPAR::SCATRA::consistency_no,
                                  INPAR::SCATRA::consistency_l2_projection_lumped),
                               &scatradyn_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fs3idyn = list->sublist(
      "FS3I CONTROL",
      false,
      "control parameters for FS3I problems\n");

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&fs3idyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&fs3idyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&fs3idyn);
  IntParameter("ITEMAX",10,"Maximum number of outer iterations",&fs3idyn);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&fs3idyn);
  IntParameter("UPRES",1,"Increment for writing solution",&fs3idyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&fs3idyn);
  setStringToIntegralParameter<int>("SCATRA_SOLVERTYPE","nonlinear",
                               "type of scalar transport solver",
                               tuple<std::string>(
                                 "linear",
                                 "nonlinear"
                                 ),
                               tuple<int>(
                                 INPAR::SCATRA::solvertype_linear_incremental,
                                 INPAR::SCATRA::solvertype_nonlinear),
                               &fs3idyn);
  BoolParameter("INF_PERM","yes","Flag for infinite permeability",&fs3idyn);
  setStringToIntegralParameter<int>("CONSTHERMPRESS","Yes",
                               "treatment of thermodynamic pressure in time",
                               tuple<std::string>(
                                 "No_energy",
                                 "No_mass",
                                 "Yes"
                                 ),
                               tuple<int>(0,1,2),
                               &fs3idyn);

  // number of linear solver used for fs3i problems
  IntParameter("COUPLED_LINEAR_SOLVER",-1,"number of linear solver used for fs3i problem",&fs3idyn);
  IntParameter("LINEAR_SOLVER1",-1,"number of linear solver used for fluid problem",&fs3idyn);
  IntParameter("LINEAR_SOLVER2",-1,"number of linear solver used for structural problem",&fs3idyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& lomacontrol = list->sublist(
      "LOMA CONTROL",
      false,
      "control parameters for low-Mach-number flow problems\n");

  BoolParameter("MONOLITHIC","no","monolithic solver",&lomacontrol);
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
  BoolParameter("SGS_MATERIAL_UPDATE","no","update material by adding subgrid-scale scalar field",&lomacontrol);

  // number of linear solver used for LOMA solver
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for LOMA problem",&lomacontrol);

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
  // parameter for possible types of ELCH algorithms for deforming meshes
  setStringToIntegralParameter<int>("MOVINGBOUNDARY",
                              "No",
                              "ELCH algorithm for deforming meshes",
                               tuple<std::string>(
                                 "No",
                                 "pseudo-transient",
                                 "fully-transient"),
                               tuple<std::string>(
                                 "no moving boundary algorithm",
                                 "pseudo-transient moving boundary algorithm",
                                 "full moving boundary algorithm including fluid solve")  ,
                                tuple<int>(
                                  INPAR::ELCH::elch_mov_bndry_no,
                                  INPAR::ELCH::elch_mov_bndry_pseudo_transient,
                                  INPAR::ELCH::elch_mov_bndry_fully_transient),
                                 &elchcontrol);
  DoubleParameter("MOLARVOLUME",0.0,"Molar volume for electrode shape change computations",&elchcontrol);
  DoubleParameter("MOVBOUNDARYTHETA",0.0,"One-step-theta factor for electrode shape change computations",&elchcontrol);
  BoolParameter("NATURAL_CONVECTION","No","Include natural convection effects",&elchcontrol);
  BoolParameter("GALVANOSTATIC","No","flag for galvanostatic mode",&elchcontrol);
  IntParameter("GSTATCONDID_CATHODE",0,"condition id of electrode kinetics for cathode",&elchcontrol);
  IntParameter("GSTATCONDID_ANODE",1,"condition id of electrode kinetics for anode",&elchcontrol);
  DoubleParameter("GSTATCONVTOL",1.e-5,"Convergence check tolerance for galvanostatic mode",&elchcontrol);
  DoubleParameter("GSTATCURTOL",1.e-15,"Current Tolerance",&elchcontrol);
  IntParameter("GSTATCURVENO",-1,"function number defining the imposed current curve",&elchcontrol);
  IntParameter("GSTATITEMAX",10,"maximum number of iterations for galvanostatic mode",&elchcontrol);
  DoubleParameter("GSTAT_LENGTH_CURRENTPATH",0.0,"average length of the current path",&elchcontrol);
  IntParameter("MAGNETICFIELD_FUNCNO",-1,"function number defining an externally imposed magnetic field",&elchcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& biofilmcontrol = list->sublist(
      "BIOFILM CONTROL",
      false,
      "control parameters for biofilm problems\n");

  BoolParameter("SURFACEGROWTH","No","Scatra algorithm for surface growth",&biofilmcontrol);
  DoubleParameter("GROWNVOLUME",0.0,"Volume for surface growth",&biofilmcontrol);
  DoubleParameter("BIOTIMESTEP",0.05,"Time step size for surface grown",&biofilmcontrol);
  IntParameter("BIONUMSTEP",0,"Maximum number of steps for surface grown",&biofilmcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptcontrol = list->sublist("TOPOLOGY OPTIMIZATION CONTROL",false,
      "control parameters for topology optimization problems");

  DoubleParameter("MAXTIME",10.0,"Total simulation time",&topoptcontrol);
  IntParameter("NUMSTEP",100,"Total number of timesteps",&topoptcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&topoptcontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&topoptcontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&topoptcontrol);

  setStringToIntegralParameter<int>("RESTART_ACTION","Fluid_Time_Step","Startint field of Restart",
      tuple<std::string>(
          "Fluid_Time_Step",
          "Adjoint_Time_Step",
          "Evaluated_Gradient",
          "Finished_Optimization_Step"),
          tuple<int>(
              INPAR::TOPOPT::fluid,
              INPAR::TOPOPT::adjoint,
              INPAR::TOPOPT::gradient,
              INPAR::TOPOPT::opti_step),
              &topoptcontrol);

  setStringToIntegralParameter<int>("CONV_CHECK_TYPE","Residuum","Convergence check due to given fields",
      tuple<std::string>(
          "Increment",
          "Residuum",
          "Increment_and_Residuum"),
          tuple<int>(
              INPAR::TOPOPT::inc,
              INPAR::TOPOPT::res,
              INPAR::TOPOPT::inc_and_res),
              &topoptcontrol);

  DoubleParameter("RESTOL",1e-5,"Convergence tolerance of the objective function",&topoptcontrol);
  DoubleParameter("INCTOL",1e-5,"Convergence tolerance of the optimized variable",&topoptcontrol);

  BoolParameter("OBJECTIVE_DISSIPATION","No","dissipation part of the objective function",&topoptcontrol);
  BoolParameter("OBJECTIVE_PRESSURE_DROP","No","pressure drop part of the objective function",&topoptcontrol);

  DoubleParameter("DISSIPATION_FAC",0.0,"factor for the dissipation part of the objective",&topoptcontrol);
  DoubleParameter("PRESSURE_DROP_FAC",0.0,"factor for the mean pressure drop part of the objective",&topoptcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptoptimizer = topoptcontrol.sublist("TOPOLOGY OPTIMIZER",false,
      "control parameters for the optimizer of a topology optimization problem");

  IntParameter("MAX_ITER",-1,"Maximal number of optimization steps",&topoptoptimizer);
  IntParameter("MATID",-1,"Material ID for automatic mesh generation",&topoptoptimizer);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Field for density = topology optimization's optimization variable",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function"),
                               tuple<int>(
                                   INPAR::TOPOPT::initdensfield_zero_field,
                                   INPAR::TOPOPT::initdensfield_field_by_function),
                               &topoptoptimizer);
  IntParameter("INITFUNCNO",-1,"function number for initial density field in topology optimization",&topoptoptimizer);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptadjointfluiddyn = topoptcontrol.sublist("TOPOLOGY ADJOINT FLUID",false,
      "control parameters for the adjoint fluid of a topology optimization problem");

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field","Initial field for adjoint problem",
      tuple<std::string>(
          "zero_field",
          "field_by_function"),
          tuple<int>(
              INPAR::TOPOPT::initadjointfield_zero_field,
              INPAR::TOPOPT::initadjointfield_field_by_function),
              &topoptadjointfluiddyn);

  setStringToIntegralParameter<int>("TESTCASE","test_no","test case for adjoint problem",
      tuple<std::string>(
          "test_no",
          "test_stat_const_vel_lin_pres",
          "test_stat_lin_vel_quad_pres",
          "test_stat_quad_vel_lin_pres",
          "test_stat_all_terms_all_constants",
          "test_instat_varying_theta",
          "test_instat_all_terms_all_constants",
          "test_instat_primal_and_dual"),
          tuple<int>(
              INPAR::TOPOPT::adjointtest_no,
              INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres,
              INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres,
              INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres,
              INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants,
              INPAR::TOPOPT::adjointtest_instat_varying_theta,
              INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants,
              INPAR::TOPOPT::adjointtest_instat_primal_and_dual),
              &topoptadjointfluiddyn);

  IntParameter("INITFUNCNO",-1,"Function for initial field",&topoptadjointfluiddyn);

  DoubleParameter("THETA_PRES",1.0,"One-Step-Theta-factor for pressure terms",&topoptadjointfluiddyn);
  DoubleParameter("THETA_DIV",1.0,"One-Step-Theta-factor for divergence terms",&topoptadjointfluiddyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrol = list->sublist("COMBUSTION CONTROL",false,
      "control parameters for a combustion problem");

  DoubleParameter("MAXTIME",10.0,"Total simulation time",&combustcontrol);
  IntParameter("NUMSTEP",100,"Total number of timesteps",&combustcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&combustcontrol);
  IntParameter("ITEMAX",1,"Total number of FG iterations",&combustcontrol);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields",&combustcontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&combustcontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&combustcontrol);
  DoubleParameter("PARALLEL_REDIST_RATIO_FAC",0.0,"Factor by which the max to min element evaluation time has to increase to trigger a redistribution",&combustcontrol);
  setStringToIntegralParameter<int>("TIMEINT","One_Step_Theta","Time Integration Scheme",
      tuple<std::string>(
          "Stationary",
          "One_Step_Theta",
          "Gen_Alpha"),
          tuple<int>(
              INPAR::FLUID::timeint_stationary,
              INPAR::FLUID::timeint_one_step_theta,
              INPAR::FLUID::timeint_afgenalpha),
              &combustcontrol);

  BoolParameter("RESTART_FROM_FLUID","No","Restart from a standard fluid problem (no scalar transport field). No XFEM dofs allowed!",&combustcontrol);
  BoolParameter("RESTART_SCATRA_INPUT","No","Use ScaTra field from .dat-file instead",&combustcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolfluid = combustcontrol.sublist("COMBUSTION FLUID",false,
      "control parameters for the fluid field of a combustion problem");

  setStringToIntegralParameter<int>("COMBUSTTYPE","Premixed_Combustion",
      "Type of combustion problem",
      tuple<std::string>(
          "Premixed_Combustion",
          "Two_Phase_Flow",
          "Two_Phase_Flow_Surf",
          "Two_Phase_Flow_Jumps"),
          tuple<int>(
              INPAR::COMBUST::combusttype_premixedcombustion,
              INPAR::COMBUST::combusttype_twophaseflow,
              INPAR::COMBUST::combusttype_twophaseflow_surf,
              INPAR::COMBUST::combusttype_twophaseflowjump),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMINTEGRATION","Cut",
      "Type of integration strategy for intersected elements",
      tuple<std::string>(
          "Cut",
          "Tetrahedra",
          "Hexahedra"),
          tuple<int>(
              INPAR::COMBUST::xfemintegration_cut,
              INPAR::COMBUST::xfemintegration_tetrahedra,
              INPAR::COMBUST::xfemintegration_hexahedra),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMTIMEINT","DoNothing",
      "Type of time integration strategy for standard degrees of freedom",
      tuple<std::string>(
          "DoNothing",
          "SemiLagrange",
          "Extrapolation",
          "MixedSemiLagrangeExtrapolation",
          "MixedGhostvalSemiLagrange",
          "MixedGhostvalExtrapolation",
          "MixedGhostvalSemiLagrangeExtrapolation"),
          tuple<int>(
              INPAR::COMBUST::xfemtimeint_donothing,
              INPAR::COMBUST::xfemtimeint_semilagrange,
              INPAR::COMBUST::xfemtimeint_extrapolation,
              INPAR::COMBUST::xfemtimeint_mixedSLExtrapol,
              INPAR::COMBUST::xfemtimeint_mixedghostSL,
              INPAR::COMBUST::xfemtimeint_mixedghostExtrapol,
              INPAR::COMBUST::xfemtimeint_mixedghostSLExtrapol),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMTIMEINT_ENR","DoNothing",
      "Type of time integration strategy for enrichment degrees of freedom",
      tuple<std::string>(
          "DoNothing",
          "QuasiStatic",
          "Projection",
          "ProjectionScalar",
          "Extrapolation",
          "ExtrapolationScalar"),
          tuple<int>(
              INPAR::COMBUST::xfemtimeintenr_donothing,
              INPAR::COMBUST::xfemtimeintenr_quasistatic,
              INPAR::COMBUST::xfemtimeintenr_project,
              INPAR::COMBUST::xfemtimeintenr_project_scalar,
              INPAR::COMBUST::xfemtimeintenr_extrapolate,
              INPAR::COMBUST::xfemtimeintenr_extrapolate_scalar),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("XFEMTIMEINT_ENR_COMP","Standard",
      "Type of time integration strategy for enrichment computation",
      tuple<std::string>(
          "Standard",
          "Full",
          "Minimal"),
          tuple<int>(
              INPAR::COMBUST::xfemtimeintenr_standard,
              INPAR::COMBUST::xfemtimeintenr_full,
              INPAR::COMBUST::xfemtimeintenr_minimal),
              &combustcontrolfluid);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field","Initial field for fluid problem",
      tuple<std::string>(
          "zero_field",
          "field_by_function",
          "disturbed_function_by_function",
          "flame_vortex_interaction",
          "beltrami_flow"),
          tuple<int>(
              INPAR::COMBUST::initfield_zero_field,
              INPAR::COMBUST::initfield_field_by_function,
              INPAR::COMBUST::initfield_disturbed_field_by_function,
              INPAR::COMBUST::initfield_flame_vortex_interaction,
              INPAR::COMBUST::initfield_beltrami_flow),
              &combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_ERROR","nitsche_error_none","To which analyt. solution do we compare?",
      tuple<std::string>(
          "nitsche_error_none",
          "nitsche_error_static_bubble_nxnx1",
          "nitsche_error_static_bubble_nxnxn",
          "nitsche_error_shear",
          "nitsche_error_couette_20x20x1",
          "nitsche_error_straight_bodyforce",
          "nitsche_error_ellipsoid_bubble_2D",
          "nitsche_error_ellipsoid_bubble_3D",
          "nitsche_error_beltrami"),
          tuple<int>(
              INPAR::COMBUST::nitsche_error_none,
              INPAR::COMBUST::nitsche_error_static_bubble_nxnx1,
              INPAR::COMBUST::nitsche_error_static_bubble_nxnxn,
              INPAR::COMBUST::nitsche_error_shear,
              INPAR::COMBUST::nitsche_error_couette_20x20x1,
              INPAR::COMBUST::nitsche_error_straight_bodyforce,
              INPAR::COMBUST::nitsche_error_ellipsoid_bubble_2D,
              INPAR::COMBUST::nitsche_error_ellipsoid_bubble_3D,
              INPAR::COMBUST::nitsche_error_beltrami),
              &combustcontrolfluid);
  setStringToIntegralParameter<int>("SURFTENSAPPROX","surface_tension_approx_none","Type of surface tension approximation",
      tuple<std::string>(
          "surface_tension_approx_none",
          "surface_tension_approx_fixed_curvature",
          "surface_tension_approx_divgrad",
          "surface_tension_approx_divgrad_normal",
          "surface_tension_approx_nodal_curvature",
          "surface_tension_approx_laplacebeltrami",
          "surface_tension_approx_laplacebeltrami_smoothed"),
          tuple<int>(
              INPAR::COMBUST::surface_tension_approx_none,
              INPAR::COMBUST::surface_tension_approx_fixed_curvature,
              INPAR::COMBUST::surface_tension_approx_divgrad,
              INPAR::COMBUST::surface_tension_approx_divgrad_normal,
              INPAR::COMBUST::surface_tension_approx_nodal_curvature,
              INPAR::COMBUST::surface_tension_approx_laplacebeltrami,
              INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed),
              &combustcontrolfluid);
  DoubleParameter("VARIABLESURFTENS",0.0,"Variable surface tension coefficient",&combustcontrolfluid);
  setStringToIntegralParameter<int>("SMOOTHGRADPHI","smooth_grad_phi_none","Type of smoothing for grad(phi)",
      tuple<std::string>(
          "smooth_grad_phi_none",
          "smooth_grad_phi_meanvalue",
          "smooth_grad_phi_leastsquares_3D",
          "smooth_grad_phi_leastsquares_2Dx",
          "smooth_grad_phi_leastsquares_2Dy",
          "smooth_grad_phi_leastsquares_2Dz"),
          tuple<int>(
              INPAR::COMBUST::smooth_grad_phi_none,
              INPAR::COMBUST::smooth_grad_phi_meanvalue,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_3D,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz),
              &combustcontrolfluid);
  // set parameters VELOCITY_JUMP_TYPE and FLUX_JUMP_TYPE in case of CombustType Premixed_Combustion
  // Two_Phase_Flow_Jumps is equal to Premixed_Combustion & vel_jump_none & flux_jump_surface_tension
  setStringToIntegralParameter<int>("VELOCITY_JUMP_TYPE","vel_jump_none","Type of velocity jump",
      tuple<std::string>(
          "vel_jump_none",
          "vel_jump_const",
          "vel_jump_premixed_combustion"),
          tuple<int>(
              INPAR::COMBUST::vel_jump_none,
              INPAR::COMBUST::vel_jump_const,
              INPAR::COMBUST::vel_jump_premixed_combustion),
              &combustcontrolfluid);
  setStringToIntegralParameter<int>("FLUX_JUMP_TYPE","flux_jump_none","Type of flux jump",
      tuple<std::string>(
          "flux_jump_none",
          "flux_jump_const",
          "flux_jump_premixed_combustion",
          "flux_jump_surface_tension"),
          tuple<int>(
              INPAR::COMBUST::flux_jump_none,
              INPAR::COMBUST::flux_jump_const,
              INPAR::COMBUST::flux_jump_premixed_combustion,
              INPAR::COMBUST::flux_jump_surface_tension),
              &combustcontrolfluid);
  IntParameter("INITFUNCNO",-1,"Function for initial field",&combustcontrolfluid);
  // unused
  //DoubleParameter("PHI_MODIFY_TOL",1.0E-10,"We modify GfuncValues near zero",&combustcontrolfluid);
  IntParameter("ITE_MAX_FRS",1,"The maximal number of iterations between fluid and recomputation of reference solution",&combustcontrolfluid);
  DoubleParameter("LAMINAR_FLAMESPEED",1.0,"The laminar flamespeed incorporates all chemical kinetics into the problem for now",&combustcontrolfluid);
  DoubleParameter("MOL_DIFFUSIVITY",0.0,"Molecular diffusivity",&combustcontrolfluid);
  DoubleParameter("MARKSTEIN_LENGTH",0.0,"The Markstein length takes flame curvature into account",&combustcontrolfluid);
  DoubleParameter("NITSCHE_VELOCITY",100.0,"Nitsche parameter to stabilize/penalize the velocity jump",&combustcontrolfluid);
  DoubleParameter("NITSCHE_PRESSURE",0.0,"Nitsche parameter to stabilize/penalize the pressure jump",&combustcontrolfluid);
  setStringToIntegralParameter<int>("CONNECTED_INTERFACE","No","Turn refinement strategy for level set function on/off",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("SMOOTHED_BOUNDARY_INTEGRATION","No","Turn on/off type of boundary integration",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("INITSTATSOL","No","Compute stationary solution as initial solution",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  // for selection of enriched fields (velocity, pressure, velocity+pressure)
  setStringToIntegralParameter<int>("SELECTED_ENRICHMENT","both",
       "select fields which get enriched dofs",
       tuple<std::string>(
           "both",
           "velocity",
           "pressure",
           "none"),
           tuple<int>(
               INPAR::COMBUST::selectedenrichment_both,
               INPAR::COMBUST::selectedenrichment_velocity,
               INPAR::COMBUST::selectedenrichment_pressure,
               INPAR::COMBUST::selectedenrichment_none),
               &combustcontrolfluid);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolgfunc = combustcontrol.sublist("COMBUSTION GFUNCTION",false,
      "control parameters for the G-function (level set) field of a combustion problem");

  setStringToIntegralParameter<int>("REINITIALIZATION","Signed_Distance_Function",
                               "Type of reinitialization strategy for level set function",
                               tuple<std::string>(
                                 "None",
                                 "Function",
                                 "Signed_Distance_Function",
                                 "Fast_Signed_Distance_Function",
                                 "Sussman"),
                               tuple<int>(
                                 INPAR::COMBUST::reinitaction_none,
                                 INPAR::COMBUST::reinitaction_byfunction,
                                 INPAR::COMBUST::reinitaction_signeddistancefunction,
                                 INPAR::COMBUST::reinitaction_fastsigneddistancefunction,
                                 INPAR::COMBUST::reinitaction_sussman),
                               &combustcontrolgfunc);

  BoolParameter("REINIT_INITIAL","Yes","Has to level set field to be reinitialized before first time step?",&combustcontrolgfunc);

  BoolParameter("REINIT_OUTPUT","No","gmsh output for reinitialized pseudo timesteps?",&combustcontrolgfunc);
  IntParameter("REINITFUNCNO",-1,"function number for reinitialization of level set (G-function) field",&combustcontrolgfunc);
  IntParameter("REINITINTERVAL",1,"reinitialization interval",&combustcontrolgfunc);
  BoolParameter("REINITBAND","No","reinitialization only within a band around the interface, or entire domain?",&combustcontrolgfunc);
  DoubleParameter("REINITBANDWIDTH",1.0,"G-function value defining band width for reinitialization",&combustcontrolgfunc);
  BoolParameter("REINITVOLCORRECTION","No","volume correction after reinitialization",&combustcontrolgfunc);

  setStringToIntegralParameter<int>("TRANSPORT_VEL","Compute_flame_vel",
                               "Type of reinitialization strategy for level set function",
                               tuple<std::string>(
                                 "Compute_flame_vel",
                                 "Convective_vel",
                                 "Vel_by_funct"),
                               tuple<int>(
                                 INPAR::COMBUST::transport_vel_flamevel,
                                 INPAR::COMBUST::transport_vel_convvel,
                                 INPAR::COMBUST::transport_vel_byfunct),
                               &combustcontrolgfunc);

  IntParameter("TRANSPORT_VEL_FUNC",-1,"number of function to overwrite the transport velocity",&combustcontrolgfunc);

  setStringToIntegralParameter<int>("REFINEMENT","No","Turn refinement strategy for level set function on/off",
                                     yesnotuple,yesnovalue,&combustcontrolgfunc);
  IntParameter("REFINEMENTLEVEL",-1,"number of refinement level for refinement strategy",&combustcontrolgfunc);

  BoolParameter("EXTRACT_INTERFACE_VEL","No","replace computed velocity at nodes of given distance of interface by approximated interface velocity",&combustcontrolgfunc);
  IntParameter("NUM_CONVEL_LAYERS",0,"number of layers around the interface which keep their computed convective velocity",&combustcontrolgfunc);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolpdereinit = combustcontrol.sublist("COMBUSTION PDE REINITIALIZATION",false,"");

  setStringToIntegralParameter<int>("REINIT_METHOD", "PDE_Based_Characteristic_Galerkin",
                                 "Type of reinitialization method",
                                 tuple<std::string>(
                                  "PDE_Based_Characteristic_Galerkin",
                                  "PDE_Based_Linear_Convection",
                                  "None"),
                                tuple<int>(
                                  INPAR::SCATRA::reinitstrategy_pdebased_characteristic_galerkin,
                                  INPAR::SCATRA::reinitstrategy_pdebased_linear_convection,
                                  INPAR::SCATRA::reinitstrategy_none),
                                  &combustcontrolpdereinit);

  setStringToIntegralParameter<int>("REINIT_TIMEINTEGR","Taylor_Galerkin_2",
                               "Time Integration Scheme for PDE-based reinitialization",
                               tuple<std::string>(
                                 "Taylor_Galerkin_2",
                                 "One_Step_Theta"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::timeint_tg2,
                                   INPAR::SCATRA::timeint_one_step_theta
                                 ),
                               &combustcontrolpdereinit);

  setStringToIntegralParameter<int>("STATIONARY_CHECK", "Integrated_L1_Norm",
                                 "Type of check for stationary solution",
                                 tuple<std::string>(
                                  "Integrated_L1_Norm",
                                  "Num_Steps"),
                                tuple<int>(
                                  INPAR::SCATRA::reinit_stationarycheck_L1normintegrated,
                                  INPAR::SCATRA::reinit_stationarycheck_numsteps),
                                  &combustcontrolpdereinit);

  setStringToIntegralParameter<int>("SMOOTHED_SIGN_TYPE", "LinEtAl2005",
                                 "Type of check for stationary solution",
                                 tuple<std::string>(
                                  "NonSmoothed",
                                  "Nagrath2005",
                                  "LinEtAl2005",
                                  "LinEtAl_Normalized"),
                                tuple<int>(
                                  INPAR::SCATRA::signtype_nonsmoothed,
                                  INPAR::SCATRA::signtype_Nagrath2005,
                                  INPAR::SCATRA::signtype_LinEtAl2005,
                                  INPAR::SCATRA::signtype_LinEtAl_normalized),
                                  &combustcontrolpdereinit);

  setStringToIntegralParameter<int>("PENALTY_METHOD", "None",
                                 "Type of interface penalizing",
                                 tuple<std::string>(
                                  "None",
                                  "Intersection_Points",
                                  "Akkerman"),
                                tuple<int>(
                                  INPAR::SCATRA::penalty_method_none,
                                  INPAR::SCATRA::penalty_method_intersection_points,
                                  INPAR::SCATRA::penalty_method_akkerman),
                                  &combustcontrolpdereinit);

  IntParameter("NUMPSEUDOSTEPS", 5, "number of pseudo time steps for pde based reinitialization", &combustcontrolpdereinit);
  DoubleParameter("PSEUDOTIMESTEP_FACTOR", 0.1, "factor for pseudo time step size for pde based reinitialization (factor*meshsize)", &combustcontrolpdereinit);
  DoubleParameter("PENALTY_INTERFACE", 100.0, "factor for penalizing interface moving during ", &combustcontrolpdereinit);
  DoubleParameter("EPSILON_BANDWIDTH", 2.0, "factor for interface thickness (multiplied by meshsize)", &combustcontrolpdereinit);
  DoubleParameter("SHOCK_CAPTURING_DIFFUSIVITY", 1e-003, "diffusivity used for shock capturing operator", &combustcontrolpdereinit);
  BoolParameter("SHOCK_CAPTURING","no","Switch on shock capturing",&combustcontrolpdereinit);
  /*----------------------------------------------------------------------*/


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fsidyn = list->sublist(
    "FSI DYNAMIC",false,
    "Fluid Structure Interaction\n"
    "Partitioned FSI solver with various coupling methods"
    );

  Teuchos::Tuple<std::string,20> name;
  Teuchos::Tuple<int,20> label;

  name[ 0] = "basic_sequ_stagg";                   label[ 0] = fsi_basic_sequ_stagg;
  name[ 1] = "iter_stagg_fixed_rel_param";         label[ 1] = fsi_iter_stagg_fixed_rel_param;
  name[ 2] = "iter_stagg_AITKEN_rel_param";        label[ 2] = fsi_iter_stagg_AITKEN_rel_param;
  name[ 3] = "iter_stagg_steep_desc";              label[ 3] = fsi_iter_stagg_steep_desc;
  name[ 4] = "iter_stagg_NLCG";                    label[ 4] = fsi_iter_stagg_NLCG;
  name[ 5] = "iter_stagg_MFNK_FD";                 label[ 5] = fsi_iter_stagg_MFNK_FD;
  name[ 6] = "iter_stagg_MFNK_FSI";                label[ 6] = fsi_iter_stagg_MFNK_FSI;
  name[ 7] = "iter_stagg_MPE";                     label[ 7] = fsi_iter_stagg_MPE;
  name[ 8] = "iter_stagg_RRE";                     label[ 8] = fsi_iter_stagg_RRE;
  name[ 9] = "iter_monolithicfluidsplit";          label[ 9] = fsi_iter_monolithicfluidsplit;
  name[10] = "iter_monolithicstructuresplit";      label[10] = fsi_iter_monolithicstructuresplit;
  name[11] = "iter_lung_monolithicstructuresplit"; label[11] = fsi_iter_lung_monolithicstructuresplit;
  name[12] = "iter_lung_monolithicfluidsplit";     label[12] = fsi_iter_lung_monolithicfluidsplit;
  name[13] = "iter_monolithicxfem";                label[13] = fsi_iter_monolithicxfem;
  name[14] = "pseudo_structure";                   label[14] = fsi_pseudo_structureale;
  name[15] = "iter_constr_monolithicfluidsplit";     label[15] = fsi_iter_constr_monolithicfluidsplit;
  name[16] = "iter_constr_monolithicstructuresplit";     label[16] = fsi_iter_constr_monolithicstructuresplit;
  name[17] = "iter_mortar_monolithicstructuresplit";     label[17] = fsi_iter_mortar_monolithicstructuresplit;
  name[18] = "iter_mortar_monolithicfluidsplit";     label[18] = fsi_iter_mortar_monolithicfluidsplit;
  name[19] = "iter_fluidfluid_monolithicstructuresplit";     label[19] = fsi_iter_fluidfluid_monolithicstructuresplit;

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_AITKEN_rel_param",
                                    "Iteration Scheme over the fields",
                                    name,
                                    label,
                                    &fsidyn);

  setStringToIntegralParameter<int>(
                               "PARTITIONED","DirichletNeumann",
                               "Coupling strategies for partitioned FSI solvers. Most of the time Dirichlet-Neumann is just right.",
                               tuple<std::string>(
                                 "DirichletNeumann",
                                 "DirichletNeumannSlideALE"
                                 ),
                               tuple<int>(
                                 INPAR::FSI::DirichletNeumann,
                                 INPAR::FSI::DirichletNeumannSlideale
                                 ),
                               &fsidyn);

  setStringToIntegralParameter<int>("DEBUGOUTPUT","No",
                                    "Output of unconverged interface values during partitioned FSI iteration.\n"
                                    "There will be a new control file for each time step.\n"
                                    "This might be helpful to understand the coupling iteration.",
                                    tuple<std::string>(
                                      "No",
                                      "Yes",
                                      "Interface",
                                      "Preconditioner",
                                      "All"
                                      ),
                                    tuple<int>(
                                      0,
                                      1,
                                      1,
                                      2,
                                      256
                                      ),
                                    &fsidyn);

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

  setStringToIntegralParameter<int>(
                               "LINEARBLOCKSOLVER","PreconditionedKrylov",
                               "Linear solver algorithm for monolithic block system in monolithic FSI.\n"
                               "Most of the time preconditioned Krylov is the right thing to choose. But there are\n"
                               "block Gauss-Seidel methods as well.",
                               tuple<std::string>(
                                 "PreconditionedKrylov",
                                 "FSIAMG"
                                 ),
                               tuple<int>(
                                 INPAR::FSI::PreconditionedKrylov,
                                 INPAR::FSI::FSIAMG
                                 ),
                               &fsidyn);

  setStringToIntegralParameter<int>("FSIAMGANALYZE","No",
                               "run analysis on fsiamg multigrid scheme\n"
                               "Supported in monolithic FSI for now.",
                               yesnotuple,yesnovalue,&fsidyn);

  IntParameter("NUMSTEP",200,"Total number of Timesteps",&fsidyn);
  IntParameter("ITEMAX",100,"Maximum number of iterations over fields",&fsidyn);
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

  // Iterationparameters for convergence check of newton loop
  setStringToIntegralParameter<int>("NORM_INC","Rel","type of norm for primary variables convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"
                                 ),
                               tuple<int>(
                                 INPAR::FSI::convnorm_abs,
                                 INPAR::FSI::convnorm_rel,
                                 INPAR::FSI::convnorm_mix
                                 ),
                               &fsidyn);

  setStringToIntegralParameter<int>("NORM_RESF","Rel","type of norm for residual convergence check",
                                 tuple<std::string>(
                                   "Abs",
                                   "Rel",
                                   "Mix"
                                   ),
                                 tuple<int>(
                                   INPAR::FSI::convnorm_abs,
                                   INPAR::FSI::convnorm_rel,
                                   INPAR::FSI::convnorm_mix
                                   ),
                                 &fsidyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                               tuple<std::string>("And"),
                               tuple<int>(
                                 INPAR::FSI::bop_and
                                 ),
                               &fsidyn);

  // monolithic preconditioner parameter

  setNumericStringParameter("STRUCTPCOMEGA","1.0 1.0 1.0 1.0",
                  "Relaxation factor for Richardson iteration on structural block in MFSI block preconditioner",
                  &fsidyn);
  setNumericStringParameter("STRUCTPCITER","1 1 1 1",
               "Number of Richardson iterations on structural block in MFSI block preconditioner",
               &fsidyn);
  setNumericStringParameter("FLUIDPCOMEGA","1.0 1.0 1.0 1.0",
                  "Relaxation factor for Richardson iteration on fluid block in MFSI block preconditioner",
                  &fsidyn);
  setNumericStringParameter("FLUIDPCITER","1 1 1 1",
               "Number of Richardson iterations on fluid block in MFSI block preconditioner",
               &fsidyn);
  setNumericStringParameter("ALEPCOMEGA","1.0 1.0 1.0 1.0",
                  "Relaxation factor for Richardson iteration on ale block in MFSI block preconditioner",
                  &fsidyn);
  setNumericStringParameter("ALEPCITER","1 1 1 1",
               "Number of Richardson iterations on ale block in MFSI block preconditioner",
               &fsidyn);

  setNumericStringParameter("PCOMEGA","1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on whole MFSI block preconditioner",
                            &fsidyn);
  setNumericStringParameter("PCITER","1 1 1",
                            "Number of Richardson iterations on whole MFSI block preconditioner",
                            &fsidyn);

  setNumericStringParameter("BLOCKSMOOTHER","BGS BGS BGS",
                            "Type of block smoother, can be BGS or Schur",
                            &fsidyn);

  setNumericStringParameter("SCHUROMEGA","0.001 0.01 0.1",
                            "Damping factor for Schur complement construction",
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

  setStringToIntegralParameter<int>("SLIDEALEPROJ","None",
                                 "Projection method to use for sliding FSI.",
                                 tuple<std::string>("None","Curr","Ref","RotZ","RotZSphere"),
                                 tuple<int>(
                                     INPAR::FSI::ALEprojection_none,
                                     INPAR::FSI::ALEprojection_curr,
                                     INPAR::FSI::ALEprojection_ref,
                                     INPAR::FSI::ALEprojection_rot_z,
                                     INPAR::FSI::ALEprojection_rot_zsphere),
                                 &fsidyn);

  setStringToIntegralParameter<int> ("DIVPROJECTION", "no", "Project velocity into divergence-free subspace",
                                     yesnotuple,yesnovalue,&fsidyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& constrfsi = fsidyn.sublist("CONSTRAINT",false,"");

  setStringToIntegralParameter<int> ("PRECONDITIONER","Simple","preconditioner to use",
      tuple<std::string>("Simple","Simplec"),tuple<int>(INPAR::FSI::Simple,INPAR::FSI::Simplec),&constrfsi);
  IntParameter("SIMPLEITER",2,"Number of iterations for simple pc",&constrfsi);
  DoubleParameter("ALPHA",0.8,"alpha parameter for simple pc",&constrfsi);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& search_tree = list->sublist("SEARCH TREE",false,"");

  setStringToIntegralParameter<int>("TREE_TYPE","notree","set tree type",
                                   tuple<std::string>("notree","octree3d","quadtree3d","quadtree2d"),
                                   tuple<int>(
                                     INPAR::GEO::Notree,
                                     INPAR::GEO::Octree3D,
                                     INPAR::GEO::Quadtree3D,
                                     INPAR::GEO::Quadtree2D),
                                     &search_tree);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL",false,"");


  // OUTPUT options
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT","Yes","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT_SCREEN","No","Do you want to be informed, if Gmsh output is written?",
                                 yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_SOL_OUT","Yes","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_DISCRET_OUT","Yes","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_CUT_OUT","Yes","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);


  IntParameter("MAX_NUM_DOFSETS",3,"Maximum number of volumecells in the XFEM element",&xfem_general);

  // Integration options
  setStringToIntegralParameter<int>("VOLUME_GAUSS_POINTS_BY","Tessellation","how to find Gauss Points for the cut volumes",
                               tuple<std::string>("Tessellation","MomentFitting","DirectDivergence"),
                               tuple<int>(0,1,2),
                               &xfem_general);

  setStringToIntegralParameter<int>("BOUNDARY_GAUSS_POINTS_BY","Tessellation","how to find Gauss Points for the boundary cells",
                                 tuple<std::string>("Tessellation","MomentFitting","DirectDivergence"),
                                 tuple<int>(0,1,2),
                                 &xfem_general);

  setStringToIntegralParameter<int>("2DFLOW",
                               "no",
                               "Flag needed for pseudo 2D-simulations for fluid-fluid-Coupling",
                               tuple<std::string>(
                                 "no",
                                 "yes"),
                               tuple<std::string>(
                                 "No 2D-Simulation",
                                 "2D-Simulation"),
                               tuple<int>(0,1),
                               &xfem_general);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_dyn = list->sublist("XFLUID DYNAMIC",false,"");

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_general = xfluid_dyn.sublist("GENERAL",false,"");

  setStringToIntegralParameter<int>("XFLUID_BOUNDARY","xfluid_moving_boundary","moving boundary or stationary boundary or xfsi boundary",
                               tuple<std::string>("xfluid_moving_boundary","xfluid_stationary_boundary", "xfsi_moving_boundary"),
                               tuple<int>(
                                   INPAR::XFEM::XFluidMovingBoundary,       // moving boundary
                                   INPAR::XFEM::XFluidStationaryBoundary,    // stationary boundary
                                   INPAR::XFEM::XFSIMovingBoundary
                                   ),
                               &xfluid_general);


  setStringToIntegralParameter<int>("INTERFACE_VEL_INITIAL","interface_vel_init_zero","how to compute or define the initial interface velocity",
                               tuple<std::string>("interface_vel_init_by_funct", "interface_vel_init_zero"),
                               tuple<int>(
                                   INPAR::XFEM::interface_vel_init_by_funct,   // define interface velocity by function
                                   INPAR::XFEM::interface_vel_init_zero        // zero interface velocity function
                                   ),
                               &xfluid_general);


  setStringToIntegralParameter<int>("INTERFACE_VEL","interface_vel_by_disp","how to compute or define the interface velocity",
                               tuple<std::string>("interface_vel_by_disp", "interface_vel_by_funct", "interface_vel_zero"),
                               tuple<int>(
                                   INPAR::XFEM::interface_vel_by_disp,    // define interface velocity by displacement of solid
                                   INPAR::XFEM::interface_vel_by_funct,   // define interface velocity by function
                                   INPAR::XFEM::interface_vel_zero        // zero interface velocity function
                                   ),
                               &xfluid_general);

  setStringToIntegralParameter<int>("INTERFACE_DISP","interface_disp_zero","how to define the interface displacement",
                               tuple<std::string>("interface_disp_by_fsi", "interface_disp_by_funct", "interface_disp_zero"),
                               tuple<int>(
                                   INPAR::XFEM::interface_disp_by_fsi,     // define interface displacement by structure solution of fsi algo
                                   INPAR::XFEM::interface_disp_by_funct,   // define interface displacement by function
                                   INPAR::XFEM::interface_disp_zero        // zero interface displacement function
                                   ),
                               &xfluid_general);

  IntParameter("DISP_FUNCT_NO",-1,"funct number for interface displacement",&xfluid_general);
  IntParameter("DISP_CURVE_NO",-1,"curve number for interface displacement",&xfluid_general);

  IntParameter("VEL_FUNCT_NO",-1,"funct number for WDBC or Neumann Condition at embedded boundary/interface",&xfluid_general);
  IntParameter("VEL_INIT_FUNCT_NO",-1,"funct number for initial interface velocity",&xfluid_general);

  // How many monolithic steps we keep the fluidfluid-boundary fixed
  IntParameter("RELAXING_ALE_EVERY",1,"Relaxing Ale after how many monolithic steps",&xfluid_general);

  BoolParameter("RELAXING_ALE","yes","switch on/off for relaxing Ale in monolithic fluid-fluid-fsi",&xfluid_general);

  // xfluidfluid-fsi-monolithic approach
  setStringToIntegralParameter<int>("MONOLITHIC_XFFSI_APPROACH","xffsi_fixedALE_partitioned","The monolithic apporach for xfluidfluid-fsi",
                                    tuple<std::string>("xffsi_full_newton", "xffsi_fixedALE_interpolation", "xffsi_fixedALE_partitioned"),
                                    tuple<int>(
                                      INPAR::XFEM::XFFSI_Full_Newton,    //xffsi with no fixed xfem-coupling
                                      INPAR::XFEM::XFFSI_FixedALE_Interpolation,  // xffsi with fixed xfem-coupling in every newtonstep
                                                                                  // and interpolations for embedded-dis afterwards
                                      INPAR::XFEM::XFFSI_FixedALE_Partitioned      // xffsi with fixed xfem-coupling in every newtonstep
                                                                                  // and solving fluid-field again
                                      ),
                                    &xfluid_general);

  // xfluidfluid time integration approach
  setStringToIntegralParameter<int>("XFLUIDFLUID_TIMEINT","Xff_TimeInt_FullProj","The xfluidfluid-timeintegration approach",
                                    tuple<std::string>("Xff_TimeInt_FullProj", "Xff_TimeInt_ProjIfMoved","Xff_TimeInt_KeepGhostValues"),
                                    tuple<int>(
                                      INPAR::XFEM::Xff_TimeInt_FullProj   ,    //always project nodes from embedded to background nodes
                                      INPAR::XFEM::Xff_TimeInt_ProjIfMoved,     //project nodes just if the status of background nodes changed
                                      INPAR::XFEM::Xff_TimeInt_KeepGhostValues  //always keep the ghost values of the background discretization
                                      ),
                                    &xfluid_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_stab = xfluid_dyn.sublist("STABILIZATION",false,"");

  // Boundary-Coupling options
  setStringToIntegralParameter<int>("EMBEDDED_BOUNDARY","BoundaryTypeSigma","method how to enforce embedded boundary/coupling conditions at the interface",
                               tuple<std::string>("BoundaryTypeSigma", "BoundaryTypeNitsche", "BoundaryTypeNeumann"),
                               tuple<int>(
                                   INPAR::XFEM::BoundaryTypeSigma,       // stress/hybrid formulation
                                   INPAR::XFEM::BoundaryTypeNitsche,     // Nitsche's formulation
                                   INPAR::XFEM::BoundaryTypeNeumann      // interior Neumann condition
                                   ),
                               &xfluid_stab);


  setStringToIntegralParameter<int>("COUPLING_STRATEGY","xfluid_sided_mortaring","on which side is the main part of enforcing",
                               tuple<std::string>("xfluid_sided_mortaring", "embedded_sided_mortaring", "two_sided_mortaring"),
                               tuple<int>(
                                   INPAR::XFEM::Xfluid_Sided_Mortaring,    // coupling on cut mesh at interface
                                   INPAR::XFEM::Embedded_Sided_Mortaring,  // coupling on embedded mesh at interface
                                   INPAR::XFEM::Two_Sided_Mortaring        // coupling on cut mesh and embedded mesh at interface
                                   ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("MSH_L2_PROJ","part_ele_proj","perform the L2 projection between stress fields on whole element or on fluid part?",
                               tuple<std::string>("full_ele_proj", "part_ele_proj"),
                               tuple<int>(
                                   INPAR::XFEM::MSH_L2_Proj_full,   // L2 stress projection on whole fluid element
                                   INPAR::XFEM::MSH_L2_Proj_part    // L2 stress projection on partial fluid element volume
                                   ),
                               &xfluid_stab);

  // viscous and convective Nitsche/MSH stabilization parameter
  DoubleParameter("VISC_STAB_FAC",      1.0, "define stabilization parameter for viscous part of interface stabilization (Nitsche, MSH)",&xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_SCALING","visc_div_by_hk","scaling factor for viscous interface stabilization (Nitsche, MSH)",
                               tuple<std::string>("visc_div_by_hk", "inv_hk", "const"),
                               tuple<int>(
                                   INPAR::XFEM::ViscStabScaling_visc_div_by_hk,   // scaling with mu/hk
                                   INPAR::XFEM::ViscStabScaling_inv_hk,           // scaling with 1/hk
                                   INPAR::XFEM::ViscStabScaling_const             // scaling with 1.0=const
                                   ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_HK","vol_equivalent","how to define the characteristic element length in cut elements",
                                 tuple<std::string>("vol_equivalent", "vol_div_by_surf", "longest_ele_length"),
                                 tuple<int>(
                                     INPAR::XFEM::ViscStab_hk_vol_equivalent,     // volume equivalent element diameter
                                     INPAR::XFEM::ViscStab_hk_vol_div_by_surf,    // scaling with hk ~ vol/surf
                                     INPAR::XFEM::ViscStab_hk_longest_ele_length  // longest element length
                                     ),
                                 &xfluid_stab);

  DoubleParameter("CONV_STAB_FAC",       1.0, "define stabilization parameter for convective part of interface stabilization (inflow, inflow and outflow)",&xfluid_stab);

  setStringToIntegralParameter<int>("CONV_STAB_SCALING","abs_normal_vel","scaling factor for viscous interface stabilization (Nitsche, MSH)",
                               tuple<std::string>("inflow", "abs_normal_vel", "const", "none"),
                               tuple<int>(
                                   INPAR::XFEM::ConvStabScaling_inflow,           // scaling with max(0,-u*n)
                                   INPAR::XFEM::ConvStabScaling_abs_normal_vel,   // scaling with |u*n|
                                   INPAR::XFEM::ConvStabScaling_const,            // scaling with 1.0=const
                                   INPAR::XFEM::ConvStabScaling_none              // no convective stabilization
                                   ),
                               &xfluid_stab);

  BoolParameter("GHOST_PENALTY_STAB","no","switch on/off ghost penalty interface stabilization",&xfluid_stab);

  DoubleParameter("GHOST_PENALTY_FAC",       0.1, "define stabilization parameter ghost penalty interface stabilization",&xfluid_stab);

  setStringToIntegralParameter<int>("EOS_GP_PATTERN","u-v-w-p-diagonal-block","which matrix pattern shall be assembled for 'Edgebased' fluid stabilization and 'GhostPenalty' stabilization?",
                               tuple<std::string>("u-v-w-p-diagonal-block", "u-p-block", "full"),
                               tuple<int>(
                                   INPAR::XFEM::EOS_GP_Pattern_uvwp,    // u-v-w-p-diagonal-block matrix pattern
                                   INPAR::XFEM::EOS_GP_Pattern_up,      // u-p-block matrix pattern
                                   INPAR::XFEM::EOS_GP_Pattern_full     // full matrix pattern
                                   ),
                               &xfluid_stab);


  setStringToIntegralParameter<int>("DLM_CONDENSATION","Yes","Do you want to condense the discontinuous stress field?",
                               yesnotuple,yesnovalue,&xfluid_stab);

  /*----------------------------------------------------------------------*/
  // set valid parameters for solver blocks

  // Note: the maximum number of solver blocks is hardwired here. If you change this,
  // don't forget to edit the corresponding parts in globalproblems.cpp, too.
  for (int i = 1; i<10; i++) {
    std::stringstream ss;
    ss << "SOLVER " << i;
    std::stringstream ss_description;
    ss_description << "solver parameters for solver block " << i;
    Teuchos::ParameterList& solverlist = list->sublist(ss.str(),false,ss_description.str());
    SetValidSolverParameters(solverlist);
  }

  /*----------------------------------------------------------------------*/
  // UMFPACK solver section
  // some people just need a solver quickly. We provide a special paramter set
  // for UMFPACK that users can just use temporarely without introducing a
  // separate solver block.
  Teuchos::ParameterList& solver_u = list->sublist("UMFPACK SOLVER",false,"solver parameters for UMFPACK");
  SetValidSolverParameters(solver_u);

  return list;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::SetValidSolverParameters(Teuchos::ParameterList& list)
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Tuple<std::string,13> solver_name;
  Teuchos::Tuple<int,13>  solver_number;

  solver_name[0] = "Amesos_KLU_sym";               solver_number[0] = INPAR::SOLVER::amesos_klu_sym;
  solver_name[1] = "Amesos_KLU_nonsym";            solver_number[1] = INPAR::SOLVER::amesos_klu_nonsym;
  solver_name[2] = "Superlu";                      solver_number[2] = INPAR::SOLVER::superlu;
  solver_name[3] = "vm3";                          solver_number[3] = INPAR::SOLVER::vm3;
  solver_name[4] = "Aztec_MSR";                    solver_number[4] = INPAR::SOLVER::aztec_msr;
  solver_name[5] = "LAPACK_sym";                   solver_number[5] = INPAR::SOLVER::lapack_sym;
  solver_name[6] = "LAPACK_nonsym";                solver_number[6] = INPAR::SOLVER::lapack_nonsym;
  solver_name[7] = "UMFPACK";                      solver_number[7] = INPAR::SOLVER::umfpack;
  solver_name[8] = "Belos";                        solver_number[8] = INPAR::SOLVER::belos;
  solver_name[9] = "Stratimikos_Amesos";           solver_number[9] = INPAR::SOLVER::stratimikos_amesos;
  solver_name[10]= "Stratimikos_Aztec";            solver_number[10]= INPAR::SOLVER::stratimikos_aztec;
  solver_name[11]= "Stratimikos_Belos";            solver_number[11]= INPAR::SOLVER::stratimikos_belos;
  solver_name[12]= "undefined";                    solver_number[12]= INPAR::SOLVER::undefined;

  setStringToIntegralParameter<int>(
    "SOLVER", "undefined",
    "The solver to attack the system of linear equations arising of FE approach with.",
    solver_name,
    solver_number,
    &list
    );

  setStringToIntegralParameter<int>(
    "AZSOLVE", "GMRES",
    "Type of linear solver algorithm to use.",
    tuple<std::string>("CG",
                       "GMRES",
                       "GMRESR",
                       "CGS",
                       "TFQMR",
                       "BiCGSTAB",
                       "LU",
                       "FGMRES"),
    tuple<int>(INPAR::SOLVER::azsolv_CG,
               INPAR::SOLVER::azsolv_GMRES,
               INPAR::SOLVER::azsolv_GMRESR,
               INPAR::SOLVER::azsolv_CGS,
               INPAR::SOLVER::azsolv_TFQMR,
               INPAR::SOLVER::azsolv_BiCGSTAB,
               INPAR::SOLVER::azsolv_LU,
               INPAR::SOLVER::belos_FGMRES),
    &list
    );

  {
    // this one is longer than 15 and the tuple<> function does not support this,
    // so build the Tuple class directly (which can be any size)
    Teuchos::Tuple<std::string,26> name;
    Teuchos::Tuple<int,26>  number;

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
    name[17] = "AMG(Braess-Sarazin)";         number[17] = INPAR::SOLVER::azprec_AMGBS;
    name[18] = "AMG";                         number[18] = INPAR::SOLVER::azprec_AMG;
    name[19] = "BGS2x2";                      number[19] = INPAR::SOLVER::azprec_BGS2x2;
    name[20] = "BGSnxn";                      number[20] = INPAR::SOLVER::azprec_BGSnxn;
    name[21] = "TekoSIMPLE";                  number[21] = INPAR::SOLVER::azprec_TekoSIMPLE;
    name[22] = "CheapSIMPLE";                 number[22] = INPAR::SOLVER::azprec_CheapSIMPLE;
    name[23] = "MueLu_sym";                   number[23] = INPAR::SOLVER::azprec_MueLuAMG_sym;
    name[24] = "MueLu_nonsym";                number[24] = INPAR::SOLVER::azprec_MueLuAMG_nonsym;
    name[25] = "MueLu_contact";               number[25] = INPAR::SOLVER::azprec_MueLuAMG_contact;

    setStringToIntegralParameter<int>(
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
    "IFPACKOVERLAP", 0,
    "The amount of overlap used for the ifpack \"ilu\" and \"ilut\" preconditioners.",
    &list);

  IntParameter(
    "IFPACKGFILL", 0,
    "The amount of fill allowed for the internal \"ilu\" preconditioner.",
    &list);

  DoubleParameter(
    "IFPACKFILL", 1.0,
    "The amount of fill allowed for an internal \"ilut\" preconditioner.",
    &list);

  setStringToIntegralParameter<int>(
    "IFPACKCOMBINE","Add","Combine mode for Ifpack Additive Schwarz",
    tuple<std::string>("Add","Insert","Zero"),
    tuple<int>(0,1,2),
    &list);

  DoubleParameter(
    "AZDROP", 0.0,
    "The tolerance below which an entry from the factors of an internal \"ilut\"\n"
    "preconditioner will be dropped.",
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
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS","Umfpack"),
    tuple<int>(0,1,2,3,4,5,6,7,8,9),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERMED","ILU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS","Umfpack"),
    tuple<int>(0,1,2,3,4,5,6,7,8,9),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERCOARSE","Umfpack","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS","Umfpack"),
    tuple<int>(0,1,2,3,4,5,6,7,8,9),
    &list);

  // parameters for AMG(BS)
  setNumericStringParameter("AMGBS_BS_DAMPING","1.3 1.3 1.3",
                            "Relaxation factor for Braess-Sarazin smoother within AMGBS method",
                            &list);

  setNumericStringParameter("AMGBS_BS_PCSWEEPS","2 2 2",
                            "number of jacobi/sgs sweeps for smoothing/solving pressure correction equation within Braess-Sarazin. only necessary for jacobi/gauss seidel",
                            &list);

  setNumericStringParameter("AMGBS_BS_PCDAMPING","1.0 1.0 1.0",
                              "jacobi damping factors for smoothing/solving pressure correction equation within Braess-Sarazin. only necessary for jacobi/gauss seidel",
                              &list);


  setStringToIntegralParameter<int>(
    "AMGBS_PSMOOTHER_VEL","PA-AMG","Prolongation/Restriction smoothing strategy (velocity part in AMGBS preconditioner)",
    tuple<std::string>("PA-AMG","SA-AMG","PG-AMG","PG2-AMG"),
    tuple<int>(INPAR::SOLVER::PA_AMG,INPAR::SOLVER::SA_AMG,INPAR::SOLVER::PG_AMG,INPAR::SOLVER::PG2_AMG),
    &list);
  setStringToIntegralParameter<int>(
    "AMGBS_PSMOOTHER_PRE","PA-AMG","Prolongation/Restriction smoothing strategy (pressure part in AMGBS preconditioner)",
    tuple<std::string>("PA-AMG","SA-AMG","PG-AMG","PG2-AMG"),
    tuple<int>(INPAR::SOLVER::PA_AMG,INPAR::SOLVER::SA_AMG,INPAR::SOLVER::PG_AMG,INPAR::SOLVER::PG2_AMG),
    &list);

  setStringToIntegralParameter<int>(
    "AMGBS_BS_PCCOARSE","Umfpack","approximation algorithm for solving pressure correction equation (coarsest level)",
    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel stand-alone"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &list);

  setStringToIntegralParameter<int>(
    "AMGBS_BS_PCMEDIUM","Umfpack","approximation algorithm for solving pressure correction equation (medium level)",
    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel stand-alone"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &list);

  setStringToIntegralParameter<int>(
    "AMGBS_BS_PCFINE","Umfpack","approximation algorithm for solving pressure correction equation (finest level)",
    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel stand-alone"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &list);

  // switch order of blocks in BGS2x2 preconditioner
  setStringToIntegralParameter<int>(
    "BGS2X2_FLIPORDER","block0_block1_order","BGS2x2 flip order parameter",
    tuple<std::string>("block0_block1_order","block1_block0_order"),
    tuple<int>(0,1),
    &list);

  // verbosity flag (for Belos)
  IntParameter("VERBOSITY",0,"verbosity level (0=no output,... 10=extreme), for Belos only",&list);

  // the only one stratimikos specific parameter
  setNumericStringParameter("STRATIMIKOS_XMLFILE","",
                              "xml file for stratimikos parameters",
                              &list);

  // user-given name of solver block (just for beauty)
  setNumericStringParameter("NAME","No_name",
                              "User specified name for solver block",
                              &list);

  // damping parameter for SIMPLE
  DoubleParameter("SIMPLE_DAMPING",1.,"damping parameter for SIMPLE preconditioner",&list);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list)
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  setStringToIntegralParameter<int>(
    "KIND","None","Method for time step size adapivity",
    tuple<std::string>(
      "None",
      "ZienkiewiczXie",
      "AdamsBashforth2"),
    tuple<int>(
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

  setStringToIntegralParameter<int>(
    "LOCERRNORM", "Vague", "Vector norm to treat error vector with",
    tuple<std::string>(
      "Vague",
      "L1",
      "L2",
      "Rms",
      "Inf"),
    tuple<int>(
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
    Teuchos::setStringToIntegralParameter<int>(
      "Nonlinear Solver","Line Search Based","",
      st,Teuchos::tuple<int>( 0, 1, 2, 3 ),
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
    Teuchos::setStringToIntegralParameter<int>(
      "Method","Newton","",
      st,Teuchos::tuple<int>( 0, 1, 2, 3 ),
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
    Teuchos::setStringToIntegralParameter<int>(
      "Forcing Term Method","Constant","",
      forcingtermmethod,Teuchos::tuple<int>( 0, 1, 2 ),
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
    Teuchos::setStringToIntegralParameter<int>(
      "Scaling Type","None","",
      scalingtype,Teuchos::tuple<int>( 0, 1, 2, 3 ),
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
    Teuchos::setStringToIntegralParameter<int>(
      "Method","Full Step","",
      method,Teuchos::tuple<int>( 0, 1, 2, 3, 4 ),
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
    Teuchos::setStringToIntegralParameter<int>(
      "Recovery Step Type","Constant","Determines the step size to take when the line search fails",
      recoverysteptype,Teuchos::tuple<int>( 0, 1 ),
      &polynomial);
    DoubleParameter("Recovery Step",1.0,"The value of the step to take when the line search fails. Only used if the \"Recovery Step Type\" is set to \"Constant\"",&polynomial);
    Teuchos::Array<std::string> interpolationtype = Teuchos::tuple<std::string>(
      "Quadratic",
      "Quadratic3",
      "Cubic");
    Teuchos::setStringToIntegralParameter<int>(
      "Interpolation Type","Cubic","Type of interpolation that should be used",
      interpolationtype,Teuchos::tuple<int>( 0, 1, 2 ),
      &polynomial);
    DoubleParameter("Min Bounds Factor",0.1,"Choice for $gamma_{min}$, i.e., the factor that limits the minimum size of the new step based on the previous step",&polynomial);
    DoubleParameter("Max Bounds Factor",0.5,"Choice for $gamma_{max}$, i.e., the factor that limits the maximum size of the new step based on the previous step",&polynomial);
    Teuchos::Array<std::string> sufficientdecreasecondition = Teuchos::tuple<std::string>(
      "Armijo-Goldstein",
      "Ared/Pred",
      "None");
    Teuchos::setStringToIntegralParameter<int>(
      "Sufficient Decrease Condition","Armijo-Goldstein","Choice to use for the sufficient decrease condition",
      sufficientdecreasecondition,Teuchos::tuple<int>( 0, 1, 2 ),
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
    Teuchos::setStringToIntegralParameter<int>(
      "Recovery Step Type","Constant","Determines the step size to take when the line search fails",
      recoverysteptype,Teuchos::tuple<int>( 0, 1 ),
      &morethuente);
    DoubleParameter("Recovery Step",1.0,"The value of the step to take when the line search fails. Only used if the \"Recovery Step Type\" is set to \"Constant\"",&morethuente);
    Teuchos::Array<std::string> sufficientdecreasecondition = Teuchos::tuple<std::string>(
      "Armijo-Goldstein",
      "Ared/Pred",
      "None");
    Teuchos::setStringToIntegralParameter<int>(
      "Sufficient Decrease Condition","Armijo-Goldstein","Choice to use for the sufficient decrease condition",
      sufficientdecreasecondition,Teuchos::tuple<int>( 0, 1, 2 ),
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

