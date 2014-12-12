/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.cpp

\brief Setup of the list of valid input parameters

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_Array.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_any.hpp>

#include "drt_validparameters.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem_enums.H"
#include "inpar.H"
#include "inpar_meshfree.H"
#include "inpar_ale.H"
#include "inpar_artnet.H"
#include "inpar_solver.H"
#include "inpar_fluid.H"
#include "inpar_cut.H"
#include "inpar_combust.H"
#include "inpar_mortar.H"
#include "inpar_contact.H"
#include "inpar_statmech.H"
#include "inpar_fsi.H"
#include "inpar_topopt.H"
#include "inpar_scatra.H"
#include "inpar_structure.H"
#include "inpar_potential.H"
#include "inpar_problemtype.H"
#include "inpar_thermo.H"
#include "inpar_tsi.H"
#include "inpar_turbulence.H"
#include "inpar_elch.H"
#include "inpar_invanalysis.H"
#include "inpar_statinvanalysis.H"
#include "inpar_searchtree.H"
#include "inpar_xfem.H"
#include "inpar_mlmc.H"
#include "inpar_poroelast.H"
#include "inpar_poroscatra.H"
#include "inpar_immersed.H"
#include "inpar_fpsi.H"
#include "inpar_ssi.H"
#include "inpar_cavitation.H"
#include "inpar_crack.H"
#include "inpar_levelset.H"
#include "inpar_wear.H"
#include "inpar_beamcontact.H"
#include "inpar_beampotential.H"
#include "inpar_acou.H"
#include "inpar_volmortar.H"

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

  std::cout << "NAME\n"
            << "\t" << baci_build << " - simulate just about anything\n"
            << "\n"
            << "SYNOPSIS\n"
            << "\t" << baci_build << " [-h] [--help] [-p] [--parameters] [-d] [--datfile] [-ngroup=x] [-glayout=a,b,c,...] [-nptype=parallelism_type]\n"
            << "\t\tdat_name output_name [restart=y] [restartfrom=restart_file_name] [ dat_name0 output_name0 [restart=y] [restartfrom=restart_file_name] ... ] [--interactive]\n"
            << "\n"
            << "DESCRIPTION\n"
            << "\tThe am besten simulation tool in the world.\n"
            << "\n"
            << "OPTIONS\n"
            << "\t--help or -h\n"
            << "\t\tPrint this message.\n"
            << "\n"
            << "\t--parameters or -p\n"
            << "\t\tPrint a list of all available parameters for use in a dat_file.\n"
            << "\n"
            << "\t--datfile or -d\n"
            << "\t\tPrint example dat_file with all available parameters.\n"
            << "\n"
            << "\t-ngroup=x\n"
            << "\t\tSpecify the number of groups for nested parallelism. (default: 1)\n"
            << "\n"
            << "\t-glayout=a,b,c,...\n"
            << "\t\tSpecify the number of processors per group. Argument \"-ngroup\" is mandatory and must be preceding. (default: equal distribution)\n"
            << "\n"
            << "\t-nptype=parallelism_type\n"
            << "\t\tAvailable options: \"separateDatFiles\", \"everyGroupReadDatFile\" and \"copyDatFile\"; Must be set if \"-ngroup\" > 1.\n"
            << "\t\t\"diffgroupx\" can be used to compare results from separate but parallel baci runs; x must be 0 and 1 for the respective run"
            << "\n"
            << "\tdat_name\n"
            << "\t\tName of the input file (Usually *.dat)\n"
            << "\n"
            << "\toutput_name\n"
            << "\t\tPrefix of your output files.\n"
            << "\n"
            << "\trestart=y\n"
            << "\t\tRestart the simulation from step y. It always refers to the previously defined dat_name and output_name. (default: 0 or from dat_name)\n"
            << "\n"
            << "\trestartfrom=restart_file_name\n"
            << "\t\tRestart the simulation from the files prefixed with restart_file_name. (default: output_name)\n"
            << "\n"
            << "\t--interactive\n"
            << "\t\tBaci waits at the beginning for keyboard input. Helpful for parallel debugging when attaching to a single job. Must be specified at the end in the command line.\n"
            << "\n"
            << "SEE ALSO\n"
            << "\tguides/reports/global_report.pdf\n"
            << "\n"
            << "BUGS\n"
            << "\t100% bug free since 1964.\n"
            << "\n"
            << "TIPS\n"
            << "\tCan be obtain from a friendly colleague.\n"
            << "\n"
            << "\tAlso, espresso may be donated to room MW1236.\n";

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
void DRT::INPUT::PrintDefaultParameters(IO::Pstream& stream, const Teuchos::ParameterList& list)
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
void DRT::INPUT::StringParameter(std::string const &paramName,
                                 std::string const &value,
                                 std::string const &docString,
                                 Teuchos::ParameterList *paramList)
{

  // The method Teuchos::setNumericStringParameter() cannot be used for arbitrary
  // std::string parameters, since the validate() method of the underlying
  // AnyNumberParameterEntryValidator always tries to convert a given std::string to DOUBLE(s)!
  // This may cause error messages in valgrind.
  // Thus, for arbitrary std::strings, such as needed for specifying a file or solver name, for instance,
  // this method which uses a StringValidator has to be used!

  Teuchos::RCP<Teuchos::StringValidator> validator
    = Teuchos::rcp(new Teuchos::StringValidator());

  paramList->set(paramName, value, docString, validator);
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

  // define some tuples that are often used to account for different writing of certain key words
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

  IntParameter("RESTART",0,"",&type);
  DoubleParameter("RESTARTTIME", -1.0, "Used defined restart time", &type);
  setStringToIntegralParameter<int>("SHAPEFCT","Polynomial","Defines the function spaces for the spatial approximation",
                                    tuple<std::string>("Polynomial","Nurbs","Meshfree","HDG"),
                                    tuple<int>(1,0,2,3),
                                    &type);
  IntParameter("RANDSEED",-1,"Set the random seed. If < 0 use current time.",&type);

#if 0 // currently not in use
//  BoolParameter("BANDWITHOPT","No","Do bandwith optimization of dof numbering",&type);
  setStringToIntegralParameter<int>("BANDWIDTHOPT","No",
                                    "Do bandwith optimization of dof numbering",
                                    yesnotuple,yesnovalue,&type);
#endif

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& meshfree = list->sublist("MESHFREE",false,"");
  setStringToIntegralParameter<int>("TYPE","MaxEnt","Type of meshfree discretisation.",
                                    tuple<std::string>("MaxEnt","Particle","GeoDecoupled"),
                                    tuple<int>(INPAR::MESHFREE::maxent,
                                               INPAR::MESHFREE::particle,
                                               INPAR::MESHFREE::geo_decoupled),
                                    &meshfree);
  setStringToIntegralParameter<int>("NODEPOINTASSIGNMENT","procwise","Type of assignment of nodes to cells/knots.",
                                    tuple<std::string>("procwise","blockwise"),
                                    tuple<int>(INPAR::MESHFREE::procwise,
                                               INPAR::MESHFREE::blockwise),
                                    &meshfree);
  DoubleParameter("NEWTON_TOL",1e-6,"Tolerance at which Newton is considered to be converged.",&meshfree);
  DoubleParameter("NEWTON_MAXITER",10,"Maximum number of Newton steps.",&meshfree);
  IntParameter("DBC_SOLVER",-1,"Solver number for solving non-constant Dirichlet BC if necessary.",&meshfree);
  BoolParameter("PARTITION_OF_UNITY","Yes","Enforcement of partition of unity constraint",&meshfree);
  DoubleParameter("NEGATIVITY",0.0,"Decides if and to which degree negativity is allowed",&meshfree);
  DoubleParameter("VARIANCE",1,"Variance of the basis solution function prior.",&meshfree);
  DoubleParameter("RANGE_TOL",1e-6,"Threshhold at which basis solution function prior is considered nmuerically zero.",&meshfree);
  DoubleParameter("CUTOFF_RADIUS",-1.0,"Cutoff radius for influence of meshfree points on each other.",&meshfree);
  setNumericStringParameter("BIN_PER_DIR","-1 -1 -1",
                            "Number of bins per direction (x, y, z) in particle simulations.",
                            &meshfree);
  setNumericStringParameter("BOUNDINGBOX","-1e12 -1e12 -1e12 1e12 1e12 1e12",
                            "Bounding box for binning strategy in particle simulations.",
                            &meshfree);
  setStringToIntegralParameter<int>("PRIOR","Gauss","Defines the prior type of the basis solution function.",
                                    tuple<std::string>("Gauss"),
                                    tuple<int>(INPAR::MESHFREE::p_gauss),
                                    &meshfree);
  setStringToIntegralParameter<int>("S_COMPLIANCE","linear","Defines the compliance type enforced for max-ent scheme of the basis solution functions.",
                                    tuple<std::string>("linear","stream","freespace"),
                                    tuple<int>(
                                      INPAR::MESHFREE::c_linear,
                                      INPAR::MESHFREE::c_stream,
                                      INPAR::MESHFREE::c_freesp),
                                    &meshfree);

  setStringToIntegralParameter<int>("W_COMPLIANCE","linear","Defines the compliance type enforced for max-ent scheme of the basis solution functions.",
                                    tuple<std::string>("linear","stream","freespace"),
                                    tuple<int>(
                                      INPAR::MESHFREE::c_linear,
                                      INPAR::MESHFREE::c_stream,
                                      INPAR::MESHFREE::c_freesp),
                                    &meshfree);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ps = list->sublist("PATIENT SPECIFIC",false,"");

  setStringToIntegralParameter<int>("PATSPEC","No",
                                    "Triggers application of patient specific tools in discretization construction",
                                    yesnotuple,yesnovalue,&ps);

  BoolParameter("REMODEL","No","Turn remodeling on/off",&ps);
  BoolParameter("CALCINNERRADIUS","No","Compute inner radius for pure structural wall shear stress",&ps);
  BoolParameter("LINEARCENTERLINE","No","Is the centerline linear? Only important when CALCINNERRADIUS is turned on",&ps);
  setNumericStringParameter("CENTERLINEDIRECTION","-1","direction of linear centerline",&ps);
  setNumericStringParameter("CENTERLINEPOINT","-1","point on linear centerline",&ps);

  IntParameter("MAXHULUMEN",0,"max HU value within the blood lumen",&ps);
  StringParameter("CENTERLINEFILE","name.txt",
                  "filename of file containing centerline points",
                  &ps);

  setStringToIntegralParameter<int>("CALCSTRENGTH","No","Calculate strength on/off",yesnotuple,yesnovalue,&ps);
  DoubleParameter("AAA_SUBRENDIA",22.01,"subrenal diameter of the AAA",&ps);
  setStringToIntegralParameter<int>("FAMILYHIST","No","Does the patient have AAA family history",yesnotuple,yesnovalue,&ps);
  setStringToIntegralParameter<int>("MALE_PATIENT","Yes","Is the patient a male?",yesnotuple,yesnovalue,&ps);
  // historically the maximum ilt thickness was computed based on distance to orthopressure/fsi surface on luminal side of
  // the ilt. From the maximum an approximate wall thickness of 1.0 mm is subtrated (hardcoded in patspec).
  // This obviously can cause problems when the wall thickness is not constant e.g. during UQ analysis.
  // Therefore, a new more accurate method was added which needs the luminal and the outer surface of the ILT
  // as AAA surface condition. To use this mehtod set the flag below to yes.
  // The old way is kept here only to allow evaluation of the AAA database.
  setStringToIntegralParameter<int>("CALC_ACCURATE_MAX_ILT_THICK","no","Method with which the Max ILT thickness is calculated"
                                    ,yesnotuple,yesnovalue,&ps);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& io = list->sublist("IO",false,"");

  setStringToIntegralParameter<int>("OUTPUT_GMSH","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("OUTPUT_BIN","yes","Do you want to have binary output?",yesnotuple,yesnovalue,&io);

  // Output every iteration (for debugging purposes)
  setStringToIntegralParameter<int>("OUTPUT_EVERY_ITER","no","Do you desire structural displ. output every Newton iteration",
                                    yesnotuple,yesnovalue,&io);
  IntParameter("OEI_FILE_COUNTER",0,
                "Add an output name affix by introducing a additional number",
                &io);

  // Structural output
  setStringToIntegralParameter<int>("STRUCT_DISP","Yes","Output of displacements",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_VEL_ACC","No","Output of velocity and acceleration",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_SE","No","Output of strain energy",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("STRUCT_STRESS","No","Output of stress",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "Cauchy","cauchy",
                                                  "2PK", "2pk"),
                               tuple<int>(INPAR::STR::stress_none,INPAR::STR::stress_none,INPAR::STR::stress_none,
                                          INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,
                                          INPAR::STR::stress_cauchy,INPAR::STR::stress_cauchy,
                                          INPAR::STR::stress_2pk,INPAR::STR::stress_2pk),
                               &io);
  // in case of a coupled problem (e.g. TSI) the additional stresses are
  // (TSI: thermal stresses) are printed here
  setStringToIntegralParameter<int>("STRUCT_COUPLING_STRESS","No","",
                               tuple<std::string>(
                                 "No","no","NO",
                                 "Yes","yes","YES",
                                 "Cauchy","cauchy",
                                 "2PK", "2pk"),
                               tuple<int>(
                                 INPAR::STR::stress_none,INPAR::STR::stress_none,INPAR::STR::stress_none,
                                 INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,INPAR::STR::stress_2pk,
                                 INPAR::STR::stress_cauchy,INPAR::STR::stress_cauchy,
                                 INPAR::STR::stress_2pk,INPAR::STR::stress_2pk),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_STRAIN","No","Output of strains",
                               tuple<std::string>("No","no","NO",
                                                  "Yes","yes","YES",
                                                  "EA","ea",
                                                  "GL", "gl",
                                                  "LOG","log"),
                               tuple<int>(INPAR::STR::strain_none,INPAR::STR::strain_none,INPAR::STR::strain_none,
                                          INPAR::STR::strain_gl,INPAR::STR::strain_gl,INPAR::STR::strain_gl,
                                          INPAR::STR::strain_ea,INPAR::STR::strain_ea,
                                          INPAR::STR::strain_gl,INPAR::STR::strain_gl,
                                          INPAR::STR::strain_log,INPAR::STR::strain_log),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_PLASTIC_STRAIN","No","",
                               tuple<std::string>(
                                 "No","no","NO",
                                 "Yes","yes","YES",
                                 "EA","ea",
                                 "GL", "gl"),
                               tuple<int>(
                                 INPAR::STR::strain_none,INPAR::STR::strain_none,INPAR::STR::strain_none,
                                 INPAR::STR::strain_gl,INPAR::STR::strain_gl,INPAR::STR::strain_gl,
                                 INPAR::STR::strain_ea,INPAR::STR::strain_ea,
                                 INPAR::STR::strain_gl,INPAR::STR::strain_gl),
                               &io);
  setStringToIntegralParameter<int>("STRUCT_SURFACTANT","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_SOL","Yes","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_WALL_SHEAR_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_ELEDATA_EVRY_STEP","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("FLUID_VIS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("THERM_TEMPERATURE","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter<int>("THERM_HEATFLUX","None","",
                               tuple<std::string>(
                                 "None",
                                 "No",
                                 "NO",
                                 "no",
                                 "Current",
                                 "Initial"),
                               tuple<int>(
                                 INPAR::THR::heatflux_none,
                                 INPAR::THR::heatflux_none,
                                 INPAR::THR::heatflux_none,
                                 INPAR::THR::heatflux_none,
                                 INPAR::THR::heatflux_current,
                                 INPAR::THR::heatflux_initial),
                               &io);
  setStringToIntegralParameter<int>("THERM_TEMPGRAD","None","",
                               tuple<std::string>(
                                 "None",
                                 "No",
                                 "NO",
                                 "no",
                                 "Current",
                                 "Initial"),
                               tuple<int>(
                                 INPAR::THR::tempgrad_none,
                                 INPAR::THR::tempgrad_none,
                                 INPAR::THR::tempgrad_none,
                                 INPAR::THR::tempgrad_none,
                                 INPAR::THR::tempgrad_current,
                                 INPAR::THR::tempgrad_initial),
                               &io);

  IntParameter("FILESTEPS",1000,"Amount of timesteps written to a single result file",&io);
  IntParameter("STDOUTEVRY",1,"Print to screen every n step",&io);

  BoolParameter("WRITE_TO_SCREEN",  "Yes","Write screen output",                       &io);
  BoolParameter("WRITE_TO_FILE",    "No", "Write the output into a file",              &io);
  BoolParameter("PREFIX_GROUP_ID",  "No", "Put a <GroupID>: in front of every line",   &io);
  IntParameter("LIMIT_OUTP_TO_PROC", -1,  "Only the specified procs will write output",&io);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& design = list->sublist("DESIGN DESCRIPTION",false,"number of nodal clouds");

  IntParameter("NDPOINT",0,"number of points",&design);
  IntParameter("NDLINE",0,"number of line clouds",&design);
  IntParameter("NDSURF",0,"number of surface clouds",&design);
  IntParameter("NDVOL",0,"number of volume clouds",&design);
  IntParameter("NDPARTICLE",0,"number of particle clouds",&design);

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

  // Output type
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
  DoubleParameter("M_DAMP",-1.0,"",&sdyn);
  DoubleParameter("K_DAMP",-1.0,"",&sdyn);

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
  DoubleParameter("TOLWINDKESSEL",1.0E-08,
                  "tolerance in the Windkessel error norm for the newton iteration",
                  &sdyn);
  DoubleParameter("WINDKESSEL_TIMINT_THETA",0.5,
                  "theta for one-step-theta time-integration scheme of Windkessel",
                  &sdyn);
  setStringToIntegralParameter<int>("RESTART_WITH_WINDKESSEL","No","Must be chosen if a non-Windkessel simulation is to be restarted with Windkessel",
                                 yesnotuple,yesnovalue,&sdyn);
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

  setStringToIntegralParameter<int>("DIVERCONT","stop","What to do with time integration when Newton-Raphson iteration failed",
                                tuple<std::string>(
                                  "stop",
                                  "continue",
                                  "repeat_step",
                                  "halve_step",
                                  "repeat_simulation"),
                                tuple<int>(
                                  INPAR::STR::divcont_stop,
                                  INPAR::STR::divcont_continue,
                                  INPAR::STR::divcont_repeat_step,
                                  INPAR::STR::divcont_halve_step,
                                  INPAR::STR::divcont_repeat_simulation),
                                &sdyn);

  setStringToIntegralParameter<int>("NLNSOL","fullnewton","Nonlinear solution technique",
                               tuple<std::string>(
                                 "vague",
                                 "fullnewton",
                                 "modnewton",
                                 "lsnewton",
                                 "ptc",
                                 "newtonlinuzawa",
                                 "augmentedlagrange",
                                 "NoxNewtonLineSearch",
                                 "noxgeneral",
                                 "NLNSOL"),
                               tuple<int>(
                                 INPAR::STR::soltech_vague,
                                 INPAR::STR::soltech_newtonfull,
                                 INPAR::STR::soltech_newtonmod,
                                 INPAR::STR::soltech_newtonls,
                                 INPAR::STR::soltech_ptc,
                                 INPAR::STR::soltech_newtonuzawalin,
                                 INPAR::STR::soltech_newtonuzawanonlin,
                                 INPAR::STR::soltech_noxnewtonlinesearch,
                                 INPAR::STR::soltech_noxgeneral,
                                 INPAR::STR::soltech_nlnsol),
                               &sdyn);

  IntParameter("LSMAXITER",30,
               "maximum number of line search steps",
               &sdyn);
  DoubleParameter("ALPHA_LS",0.5,
                  "step reduction factor alpha in (Newton) line search scheme",
                  &sdyn);
  DoubleParameter("SIGMA_LS",1.e-4,
                  "sufficient descent factor in (Newton) line search scheme",
                  &sdyn);

  setStringToIntegralParameter<int>("MATERIALTANGENT","analytical","way of evaluating the constitutive matrix",
                               tuple<std::string>(
                                 "analytical",
                                 "finitedifferences"),
                               tuple<int>(
                                 0,1),
                               &sdyn);

  // Currently not used, but structure will be kept if someone wants to reimplement
  // AN 2013_05
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
  // Currently not used, but structure will be kept if someone wants to reimplement
  // AN 2013_05
  setNumericStringParameter("CONTROLNODE","-1 -1 -1",
                            "for methods other than load control: [node(fortran numbering)] [dof(c-numbering)] [curve(fortran numbering)]",
                            &sdyn);

  setStringToIntegralParameter<int>("LOADLIN","yes",
                                    "Use linearization of external follower load in Newton",
                                    yesnotuple,yesnovalue,&sdyn);

  setStringToIntegralParameter<int>("MASSLIN","No","Application of nonlinear inertia terms",
  tuple<std::string>("No","no",
                     "Standard", "standard",
                     "Rotations", "rotations"),

  tuple<int>(INPAR::STR::ml_none,INPAR::STR::ml_none,
             INPAR::STR::ml_standard,INPAR::STR::ml_standard,
             INPAR::STR::ml_rotations,INPAR::STR::ml_rotations),
             &sdyn);


// Since predicor "none" would be misleading, the usage of no predictor is called vague.
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

  // flag decides if young's modulus is temperature dependent, so far only available
  // for temperature-dependent St.Venant Kirchhoff material
  setStringToIntegralParameter<int>("YOUNG_IS_TEMP_DEPENDENT","No",
                               "Use temperature-dependent Young's modulus",
                               yesnotuple,yesnovalue,&sdyn);

  // where the geometry comes from
  setStringToIntegralParameter<int>(
    "GEOMETRY","full",
    "How the geometry is specified",
    tuple<std::string>(
      "full",
      "box",
      "file"),
    tuple<int>(
      INPAR::geometry_full,
      INPAR::geometry_box,
      INPAR::geometry_file),
    &sdyn);

  /*--------------------------------------------------------------------*/
  /* parameters for time step size adaptivity in structural dynamics */
  Teuchos::ParameterList& tap = sdyn.sublist("TIMEADAPTIVITY",false,"");
  SetValidTimeAdaptivityParameters(tap);

  /*----------------------------------------------------------------------*/
  /* parameters for generalised-alpha structural integrator */
  Teuchos::ParameterList& genalpha = sdyn.sublist("GENALPHA",false,"");

  setStringToIntegralParameter<int>("GENAVG","TrLike",
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
  Teuchos::ParameterList& statinvp = list->sublist("STAT INVERSE ANALYSIS",false,"");

  // Statistical Inverse Analysis switch
  setStringToIntegralParameter<int>("STAT_INV_ANALYSIS","none",
                                    "types of statistical inverse analysis and on/off switch",
                                    tuple<std::string>(
                                      "none",
                                      "GradientDescent",
                                      "MonteCarlo",
                                      "LBFGS"),
                                    tuple<int>(
                                      INPAR::INVANA::stat_inv_none,
                                      INPAR::INVANA::stat_inv_graddesc,
                                      INPAR::INVANA::stat_inv_mc,
                                      INPAR::INVANA::stat_inv_lbfgs),
                                    &statinvp);

  // initial scaling for the LBFGS algorithm
  BoolParameter("LBFGSINITSCAL","yes","want initial scaling for the LBFGS?", &statinvp);

  // step to restart from
  IntParameter("FPRESTART",0,"forward problem restart",&statinvp);

  // write restart info every so often
  IntParameter("RESTARTEVRY",1,"write restart information every x-th step",&statinvp);

  // decide which parametrization of material parameters to use
  setStringToIntegralParameter<int>("PARAMETRIZATION","none",
                                      "how to parametrize the parameter field",
                                    tuple<std::string>(
                                      "none",
                                      "kernelsmoothing",
                                      "elementwise",
                                      "uniform"),
                                    tuple<int>(
                                      INPAR::INVANA::stat_inv_mp_none,
                                      INPAR::INVANA::stat_inv_mp_smoothkernel,
                                      INPAR::INVANA::stat_inv_mp_elementwise,
                                      INPAR::INVANA::stat_inv_mp_uniform),
                                    &statinvp);

  // want some regularization
  setStringToIntegralParameter<int>("REGULARIZATION","none",
                                    "want regularization? ('tikhonov', 'totalvariation,'none')",
                                    tuple<std::string>(
                                      "none",
                                      "tikhonov",
                                      "totalvariation"),
                                    tuple<int>(
                                      INPAR::INVANA::stat_inv_reg_none,
                                      INPAR::INVANA::stat_inv_reg_tikhonov,
                                      INPAR::INVANA::stat_inv_reg_totalvariation),
                                    &statinvp);

  // want some regularization
  setStringToIntegralParameter<int>("OBJECTIVEFUNCT","none",
                                    "choose type of objective function ('displacements', 'surfcurr')",
                                    tuple<std::string>(
                                      "none",
                                      "displacements",
                                      "surfcurr"),
                                    tuple<int>(
                                      INPAR::INVANA::stat_inv_obj_none,
                                      INPAR::INVANA::stat_inv_obj_disp,
                                      INPAR::INVANA::stat_inv_obj_surfcurr),
                                    &statinvp);

  // scaling of objective function
  BoolParameter("OBJECTIVEFUNCTSCAL","No","want scaling of objective function?", &statinvp);

  // monitorfile to provide measurements
  StringParameter("MONITORFILE","none.monitor",
                  "filename of file containing measured displacements",
                  &statinvp);

  // target discretization for surface currents
  StringParameter("TARGETDISCRETIZATION", "none.dat",
                  "datfile containing target discretization",
                  &statinvp);

  // list of parameters for the respective material
  StringParameter("PARAMLIST","none",
                  "list of std::string of parameters to be optimized, order as in INV_LIST e.g. 1 YOUNG BETA",
                  &statinvp);

  // number of optimization steps
  IntParameter("MAXITER",100,"max iterations for inverse analysis",&statinvp);

  // number of optimization steps before using
  // parameter continuation in the forward problem
  IntParameter("ITERTOPC",10,"iterations before parameter continuation in the forward problem",&statinvp);

  // stepsize for deterministic gradient based schemes
  DoubleParameter("STEPSIZE",1.0,"stepsize for the gradient descent scheme",&statinvp);

  // convergence criterion tolerance
  DoubleParameter("CONVTOL",1.0e-06,"stop optimizaiton iterations for convergence criterion below this value",&statinvp);

  // weight of the regularization
  DoubleParameter("REG_WEIGHT",0.1,"weight of the regularization",&statinvp);

  // regulaization of the totalvariation functional
  DoubleParameter("TVD_EPS",0.001,"differentiation epsilon for total variation",&statinvp);

  // number of optimization steps
  IntParameter("SIZESTORAGE",20,"number of vectors to keep in storage; defaults to 20 (lbfgs usage only)",&statinvp);

  // meta parametrization of material parameters
  BoolParameter("METAPARAMS","Yes","want metaparametrization of material parameters to kept them in range?", &statinvp);

  // scale of the kernel functions used in surface currents
  DoubleParameter("KERNELSCALE", -1.0, "scale of the kernel function", &statinvp);

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

  // strategy to update regularization parameter
  setStringToIntegralParameter<int>("UPDATE_REG","RES","Update strategy for regularization parameter ",
                                    tuple<std::string>("RES","res",
                                                        "GRAD","grad"),
                                    tuple<int>(
                                        INPAR::STR::reg_update_res,INPAR::STR::reg_update_res,
                                        INPAR::STR::reg_update_grad,INPAR::STR::reg_update_grad
                                  ),
                                 &iap);


  StringParameter("MONITORFILE","none.monitor",
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
  setStringToIntegralParameter<int>("PARAM_BOUNDS","no",
                                      "Reset parameters if optstep predicts negative values",
                                      yesnotuple,yesnovalue,&iap);

  BoolParameter("PATCHES","No","Do you want to use smoothed patches?",&iap);
  StringParameter("DEFINEPATCHES","MaterialNumber",
      "define how the patches are defined: MaterialNumber or Uniform",
      &iap);
  IntParameter("NUMPATCHES",0,"number of patches",&iap);
  setNumericStringParameter("INV_LIST_PATCHES","-1",
                      "IDs of materials that are included in the patches",
                      &iap);
  IntParameter("SMOOTHINGSTEPSPATCHES",1,"number of smoothing steps that are performed",&iap);
  setNumericStringParameter("STARTVALUESFORPATCHES","1.0",
                  "startvalues for the patches, only needed for Uniform Patches",
                  &iap);

  /*----------------------------------------------------------------------*/

  /* parameters for multi-level monte carlo */
  Teuchos::ParameterList& mlmcp = list->sublist("MULTI LEVEL MONTE CARLO",
      false, "");

  //setStringToIntegralParameter<int>("MLMC","no",
  // "perform multi level monte carlo analysis",
  // yesnotuple,yesnovalue,&mlmcp);
  IntParameter("NUMRUNS", 200, "Number of Monte Carlo runs", &mlmcp);

  IntParameter("START_RUN", 0,
      "Run to start calculating the difference to lower level", &mlmcp);

  IntParameter("INITRANDOMSEED", 1000, "Random seed for first Monte Carlo run",
      &mlmcp);


  setStringToIntegralParameter<int>("REDUCED_OUTPUT", "NO",
      "Write reduced Coarse Level Output, i.e. no mesh stresses, just disp",
      yesnotuple, yesnovalue, &mlmcp);

  IntParameter("CONTNUMMAXTRIALS", 8,
      "Half stepsize CONTNUMMAXTRIALS times before giving up", &mlmcp);

  setStringToIntegralParameter<int>("FWDPROBLEM", "structure",
      "WHICH KIND OF FORWARD DO WE WANT",
      tuple<std::string>("structure", "red_airways"),
      tuple<int>(INPAR::MLMC::structure, INPAR::MLMC::red_airways), &mlmcp);

  setStringToIntegralParameter<int>(
      "UQSTRATEGY",
      "MC_PLAIN",
      "WHICH UQ STRATEGY WILL BE USED",
      tuple<std::string>("MC_PLAIN", "mc_plain", "MC_PARAMETERCONTINUATION",
          "mc_parametercontinuation", "MC_SCALEDTHICKNESS",
          "mc_scaledthickness"),
      tuple<int>(INPAR::MLMC::mc_plain, INPAR::MLMC::mc_plain,
          INPAR::MLMC::mc_paracont, INPAR::MLMC::mc_paracont,
          INPAR::MLMC::mc_scaledthick, INPAR::MLMC::mc_scaledthick), &mlmcp);

  IntParameter("NUMCONTSTEPS",2,"Number of continuation steps",&mlmcp);

   // list of materials with stochastic constitutive parameters and
  // the stochastic parameters for the respective material
  StringParameter(
      "PARAMLIST",
      "none",
      "list of std::string of parameters that are to be modelled as random fields, 1 YOUNG BETA",
      &mlmcp);
  setNumericStringParameter("OUTPUT_ELEMENT_IDS", "-1",
      "Set ID's of Output Elements, default is -1 which is none", &mlmcp);

  setStringToIntegralParameter<int>("WRITE_RV_TO_FILE", "NO",
      "Write random variables used for random field generation to file",
      yesnotuple, yesnovalue, &mlmcp);

  // For variable geometry/wall thickness
  setStringToIntegralParameter<int>("RANDOMGEOMETRY", "No",
      "Do consider random geometry/wall thickness", yesnotuple, yesnovalue,
      &mlmcp);

  // For variable geometry/wall thickness
  setStringToIntegralParameter<int>("MEAN_GEO_FROM_ALE_POINTDBC", "No",
      "Read in mean from point dbc conds", yesnotuple, yesnovalue, &mlmcp);

  IntParameter("NUMALESTEPS", 1,
      "How many ALE steps do we want to use to compute uncertain geometry",
      &mlmcp);

  DoubleParameter("INITIALTHICKNESS", 10., "wall thickness in input file",
      &mlmcp);

  DoubleParameter("Z_POS_AAA_START_RF", -1078., "Location of bifurcation of AAA including an offset ",
      &mlmcp);

  DoubleParameter("TRANSITION_WIDTH", 15., "Transition domain to blend in random wall thickness",
      &mlmcp);
  // For variable geometry/wall thickness
  setStringToIntegralParameter<int>("START_RF_ABOVE_BIFURCATION", "No",
      "Start random geometry above bifurcation", yesnotuple, yesnovalue, &mlmcp);

  // Legacy input parameters for multilevel mc
  IntParameter("WRITESTATS", 1000,
      "Write statistics to file every WRITESTATS (only for polongated Dis)",
      &mlmcp);

  // keep this input parameter incase we want to read in another discretization (thats all it does at the moment)
  setStringToIntegralParameter<int>("PROLONGATERES", "No",
      "Prolongate Displacements to finest Discretization", yesnotuple,
      yesnovalue, &mlmcp);

  // additional inputfile which should be read (currently not used, but we want to keep the feature for now)
  StringParameter(
      "DISCRETIZATION_FOR_PROLONGATION",
      "filename.dat",
      "filename of.dat file which contains discretization to which the results are prolongated",
      &mlmcp);


  /*----------------------------------------------------------------------*/
  // set valid parameters for random fields

  // Note: the maximum number of random fields is hardwired here. If you change this,
  // don't forget to edit the corresponding parts in globalproblems.cpp, too.
  for (int i = 1; i<4; i++) {
    std::stringstream ss;
    ss << "RANDOM FIELD " << i;
    std::stringstream ss_description;
    ss_description << "random field parameters for uncertainty quantification " << i;
    Teuchos::ParameterList& randomfieldlist = list->sublist(ss.str(),false,ss_description.str());
    SetValidRandomFieldParameters(randomfieldlist);
  }



  /*----------------------------------------------------------------------*/
  /* parameters for mortar coupling */
  Teuchos::ParameterList& mortar = list->sublist("MORTAR COUPLING",false,"");

  setStringToIntegralParameter<int>("LM_SHAPEFCN","Dual","Type of employed set of shape functions",
        tuple<std::string>("Dual", "dual",
                           "Standard", "standard", "std",
                           "PetrovGalerkin", "petrovgalerkin", "pg"),
        tuple<int>(
                INPAR::MORTAR::shape_dual, INPAR::MORTAR::shape_dual,
                INPAR::MORTAR::shape_standard, INPAR::MORTAR::shape_standard, INPAR::MORTAR::shape_standard,
                INPAR::MORTAR::shape_petrovgalerkin,INPAR::MORTAR::shape_petrovgalerkin,INPAR::MORTAR::shape_petrovgalerkin),
        &mortar);

  setStringToIntegralParameter<int>("SEARCH_ALGORITHM","Binarytree","Type of contact search",
       tuple<std::string>("BruteForce","bruteforce",
                          "BruteForceEleBased","bruteforceelebased",
                          "BinaryTree","Binarytree","binarytree"),
       tuple<int>(INPAR::MORTAR::search_bfele,INPAR::MORTAR::search_bfele,
                  INPAR::MORTAR::search_bfele,INPAR::MORTAR::search_bfele,
                  INPAR::MORTAR::search_binarytree,INPAR::MORTAR::search_binarytree,
                  INPAR::MORTAR::search_binarytree),
       &mortar);

  DoubleParameter("SEARCH_PARAM",0.3,"Radius / Bounding volume inflation for contact search",&mortar);

  setStringToIntegralParameter<int>("SEARCH_USE_AUX_POS","Yes","If chosen auxiliary position is used for computing dops",
                               yesnotuple,yesnovalue,&mortar);

  setStringToIntegralParameter<int>("LM_QUAD","undefined","Type of LM interpolation for quadratic FE",
       tuple<std::string>("undefined",
                          "quad", "quadratic",
                          "pwlin", "piecewiselinear",
                          "lin","linear"),
       tuple<int>(
                  INPAR::MORTAR::lagmult_undefined,
                  INPAR::MORTAR::lagmult_quad, INPAR::MORTAR::lagmult_quad,
                  INPAR::MORTAR::lagmult_pwlin, INPAR::MORTAR::lagmult_pwlin,
                  INPAR::MORTAR::lagmult_lin, INPAR::MORTAR::lagmult_lin),
       &mortar);

  setStringToIntegralParameter<int>("CROSSPOINTS","No","If chosen, multipliers are removed from crosspoints / edge nodes",
                               yesnotuple,yesnovalue,&mortar);

  setStringToIntegralParameter<int>("LM_DUAL_CONSISTENT","Yes","If chosen consistent dual shape functions are calculated",
                               yesnotuple,yesnovalue,&mortar);

  setStringToIntegralParameter<int>("LM_NODAL_SCALE","No","If chosen a nodal scaling factor is calculated for each LM",
                               yesnotuple,yesnovalue,&mortar);

  setStringToIntegralParameter<int>("HERMITE_SMOOTHING","No","If chosen hermite interface smoothing is activated for line2 elements",
                               yesnotuple,yesnovalue,&mortar);

  setStringToIntegralParameter<int>("MESH_RELOCATION","Initial","Type of mesh relocation",
      tuple<std::string>("Initial","initial",
                         "Every_Timestep", "every_timestep",
                         "No", "no"),
      tuple<int>(
              INPAR::MORTAR::relocation_initial,  INPAR::MORTAR::relocation_initial,
              INPAR::MORTAR::relocation_timestep, INPAR::MORTAR::relocation_timestep,
              INPAR::MORTAR::relocation_none,     INPAR::MORTAR::relocation_none),
      &mortar);

  setStringToIntegralParameter<int>("REDUNDANT_STORAGE","Master","Type of redundancy in interface storage",
      tuple<std::string>("All","all",
                         "Master", "master",
                         "None", "none"),
      tuple<int>(
              INPAR::MORTAR::redundant_all,    INPAR::MORTAR::redundant_all,
              INPAR::MORTAR::redundant_master, INPAR::MORTAR::redundant_master,
              INPAR::MORTAR::redundant_none,   INPAR::MORTAR::redundant_none),
      &mortar);

  setStringToIntegralParameter<int>("PARALLEL_STRATEGY","redundant_ghosting","Type of parallel interface evaluation",
      tuple<std::string>("rg","redundant_ghosting", "ghosting",
                         "rre", "roundrobinevaluate","RoundRobinEvaluate",
                         "rrg", "roundrobinghost", "RoundRobinGhost",
                         "bs", "binningstrategy", "binning"),
      tuple<int>(
              INPAR::MORTAR::ghosting_redundant, INPAR::MORTAR::ghosting_redundant,INPAR::MORTAR::ghosting_redundant,
              INPAR::MORTAR::roundrobinevaluate, INPAR::MORTAR::roundrobinevaluate,INPAR::MORTAR::roundrobinevaluate,
              INPAR::MORTAR::roundrobinghost, INPAR::MORTAR::roundrobinghost,INPAR::MORTAR::roundrobinghost,
              INPAR::MORTAR::binningstrategy, INPAR::MORTAR::binningstrategy, INPAR::MORTAR::binningstrategy),
      &mortar);

  setStringToIntegralParameter<int>("PARALLEL_REDIST","Static","Type of redistribution algorithm",
      tuple<std::string>("None","none", "No", "no",
                         "Static", "static",
                         "Dynamic", "dynamic"),
      tuple<int>(
              INPAR::MORTAR::parredist_none, INPAR::MORTAR::parredist_none,
              INPAR::MORTAR::parredist_none, INPAR::MORTAR::parredist_none,
              INPAR::MORTAR::parredist_static, INPAR::MORTAR::parredist_static,
              INPAR::MORTAR::parredist_dynamic, INPAR::MORTAR::parredist_dynamic),
      &mortar);

  setStringToIntegralParameter<int>("ALGORITHM","Mortar","Type of meshtying/contact algorithm",
      tuple<std::string>("mortar","Mortar",
                         "nts","NTS"),
      tuple<int>(
                 INPAR::MORTAR::algorithm_mortar,INPAR::MORTAR::algorithm_mortar,
                 INPAR::MORTAR::algorithm_nts,INPAR::MORTAR::algorithm_nts),
      &mortar);

  DoubleParameter("MAX_BALANCE",2.0,"Maximum value of load balance measure before parallel redistribution",&mortar);
  IntParameter("MIN_ELEPROC",0,"Minimum no. of elements per processor for parallel redistribution",&mortar);

  setStringToIntegralParameter<int>("INTTYPE","Segments","Type of numerical integration scheme",
    tuple<std::string>("Segments","segments",
             "Elements","elements",
             "Elements_BS", "elements_BS"),
    tuple<int>(
        INPAR::MORTAR::inttype_segments, INPAR::MORTAR::inttype_segments,
        INPAR::MORTAR::inttype_elements, INPAR::MORTAR::inttype_elements,
        INPAR::MORTAR::inttype_elements_BS, INPAR::MORTAR::inttype_elements_BS),
    &mortar);

  IntParameter("NUMGP_PER_DIM",0,"Number of employed integration points per dimension",&mortar);

  /*----------------------------------------------------------------------*/
  /* parameters for structural meshtying and contact */
  Teuchos::ParameterList& scontact = list->sublist("CONTACT DYNAMIC",false,"");

  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for meshtying and contact",&scontact);

  setStringToIntegralParameter<int>("RESTART_WITH_CONTACT","No","Must be chosen if a non-contact simulation is to be restarted with contact",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("ADHESION","None","Type of adhesion law",
      tuple<std::string>("None","none",
                         "bounded","b"),
      tuple<int>(
                 INPAR::CONTACT::adhesion_none,INPAR::CONTACT::adhesion_none,
                 INPAR::CONTACT::adhesion_bound,INPAR::CONTACT::adhesion_bound),
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

  setStringToIntegralParameter<int>("FRLESS_FIRST","No",
      "If chosen the first time step of a newly in contact slave node is regarded as frictionless",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("GP_SLIP_INCR","No",
      "If chosen the slip increment is computed gp-wise which results to a non-objective quantity, but this would be consistent to wear and tsi calculations.",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("STRATEGY","LagrangianMultipliers","Type of employed solving strategy",
        tuple<std::string>("LagrangianMultipliers","lagrange", "Lagrange",
                           "PenaltyMethod","penalty", "Penalty",
                           "UzawaAugementedLagrange","uzawa","Uzawa",
                           "AugmentedLagrange","augmented", "Augmented"),
        tuple<int>(
                INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult,
                INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty,
                INPAR::CONTACT::solution_uzawa, INPAR::CONTACT::solution_uzawa, INPAR::CONTACT::solution_uzawa,
                INPAR::CONTACT::solution_augmented, INPAR::CONTACT::solution_augmented, INPAR::CONTACT::solution_augmented),
        &scontact);

  setStringToIntegralParameter<int>("SYSTEM","Condensed","Type of linear system setup / solution",
        tuple<std::string>("Condensed","condensed", "cond",
                           "SaddlePoint","Saddlepoint","saddlepoint", "sp"),
        tuple<int>(
                INPAR::CONTACT::system_condensed, INPAR::CONTACT::system_condensed, INPAR::CONTACT::system_condensed,
                INPAR::CONTACT::system_saddlepoint, INPAR::CONTACT::system_saddlepoint,
                INPAR::CONTACT::system_saddlepoint, INPAR::CONTACT::system_saddlepoint),
        &scontact);

  DoubleParameter("PENALTYPARAM",0.0,"Penalty parameter for penalty / Uzawa augmented solution strategy",&scontact);
  DoubleParameter("PENALTYPARAMTAN",0.0,"Tangential penalty parameter for penalty / Uzawa augmented solution strategy",&scontact);
  IntParameter("UZAWAMAXSTEPS",10,"Maximum no. of Uzawa steps for Uzawa solution strategy",&scontact);
  DoubleParameter("UZAWACONSTRTOL",1.0e-8,"Tolerance of constraint norm for Uzawa solution strategy",&scontact);

  setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON","Yes","If chosen semi-smooth Newton concept is applied",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("SEMI_SMOOTH_CN",1.0,"Weighting factor cn for semi-smooth PDASS",&scontact);
  DoubleParameter("SEMI_SMOOTH_CT",1.0,"Weighting factor ct for semi-smooth PDASS",&scontact);

  setStringToIntegralParameter<int>(
      "CONTACTFORCE_ENDTIME",
      "No",
      "If chosen, the contact force is not evaluated at the generalized midpoint, but at the end of the time step",
      yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("VELOCITY_UPDATE","No","If chosen, velocity update method is applied",
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
              INPAR::CONTACT::errornorms_none, INPAR::CONTACT::errornorms_none,
              INPAR::CONTACT::errornorms_none, INPAR::CONTACT::errornorms_none,
              INPAR::CONTACT::errornorms_zero, INPAR::CONTACT::errornorms_zero,
              INPAR::CONTACT::errornorms_bending, INPAR::CONTACT::errornorms_bending,
              INPAR::CONTACT::errornorms_sphere, INPAR::CONTACT::errornorms_sphere,
              INPAR::CONTACT::errornorms_thicksphere, INPAR::CONTACT::errornorms_thicksphere),
      &scontact);

  setStringToIntegralParameter<int>("INITCONTACTBYGAP","No","Initialize init contact by weighted gap vector",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("INITCONTACTGAPVALUE",0.0,"Value for initialization of init contact set with gap vector",&scontact);

  // solver convergence test parameters for contact/meshtying in saddlepoint formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFCONTCONSTR","And",
    "binary operator to combine contact constraints and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::STR::bop_and,
      INPAR::STR::bop_or),
    &scontact
    );

  setStringToIntegralParameter<int>("NORMCOMBI_DISPLAGR","And",
      "binary operator to combine displacement increments and Lagrange multiplier increment values",
      tuple<std::string>(
        "And",
        "Or"),
      tuple<int>(
        INPAR::STR::bop_and,
        INPAR::STR::bop_or),
      &scontact
      );

  DoubleParameter("TOLCONTCONSTR",1.0E-6,
                  "tolerance in the contact constraint norm for the newton iteration (saddlepoint formulation only)",
                  &scontact);
  DoubleParameter("TOLLAGR",1.0E-6,
                  "tolerance in the LM norm for the newton iteration (saddlepoint formulation only)",
                  &scontact);

  setStringToIntegralParameter<int>("MESH_ADAPTIVE_CN","no",
                                     "use a scaling of cn with the local mesh size",
                                     yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("MESH_ADAPTIVE_CT","no",
                                     "use a scaling of ct with the local mesh size",
                                     yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("CONSTRAINT_DIRECTIONS","ntt",
      "formulation of constraints in normal/tangential or xyz-direction",
      tuple<std::string>(
        "ntt",
        "xyz"),
      tuple<int>(
        INPAR::CONTACT::constr_ntt,
        INPAR::CONTACT::constr_xyz),
      &scontact
      );

  /*----------------------------------------------------------------------*/
  /* parameters for volmortar */
  Teuchos::ParameterList& volmortar = list->sublist("VOLMORTAR COUPLING",false,"");

  setStringToIntegralParameter<int>("INTTYPE","Elements","Type of numerical integration scheme",
    tuple<std::string>("Elements","elements",
                       "Segments","segments"),
    tuple<int>(
        INPAR::VOLMORTAR::inttype_elements, INPAR::VOLMORTAR::inttype_elements,
        INPAR::VOLMORTAR::inttype_segments, INPAR::VOLMORTAR::inttype_segments),
    &volmortar);

  setStringToIntegralParameter<int>("COUPLINGTYPE","Volmortar","Type of coupling",
    tuple<std::string>("Volmortar","volmortar",
                       "consistentinterpolation","consint"),
    tuple<int>(
        INPAR::VOLMORTAR::couplingtype_volmortar, INPAR::VOLMORTAR::couplingtype_volmortar,
        INPAR::VOLMORTAR::couplingtype_coninter, INPAR::VOLMORTAR::couplingtype_coninter),
    &volmortar);

  setStringToIntegralParameter<int>("CUTTYPE","dd","Type of cut procedure/ integration point calculation",
    tuple<std::string>("dd","directdivergence","DirectDivergence",
                       "tessellation","t","Tessellation"),
    tuple<int>(
        INPAR::VOLMORTAR::cuttype_directdivergence, INPAR::VOLMORTAR::cuttype_directdivergence, INPAR::VOLMORTAR::cuttype_directdivergence,
        INPAR::VOLMORTAR::cuttype_tessellation, INPAR::VOLMORTAR::cuttype_tessellation, INPAR::VOLMORTAR::cuttype_tessellation),
    &volmortar);

  setStringToIntegralParameter<int>("DUALQUAD","nomod","Type of dual shape function for weighting function for quadr. problems",
    tuple<std::string>("nm","nomod",
                       "lm","lin_mod",
                       "qm","quad_mod"),
    tuple<int>(
        INPAR::VOLMORTAR::dualquad_no_mod,   INPAR::VOLMORTAR::dualquad_no_mod,
        INPAR::VOLMORTAR::dualquad_lin_mod,  INPAR::VOLMORTAR::dualquad_lin_mod,
        INPAR::VOLMORTAR::dualquad_quad_mod, INPAR::VOLMORTAR::dualquad_quad_mod),
    &volmortar);

  setStringToIntegralParameter<int>("MESH_INIT","No","If chosen, mesh initialization procedure is performed",
                               yesnotuple,yesnovalue,&volmortar);

  setStringToIntegralParameter<int>("KEEP_EXTENDEDGHOSTING","Yes","If chosen, extended ghosting is kept for simulation",
                               yesnotuple,yesnovalue,&volmortar);

  /*----------------------------------------------------------------------*/
  /* parameters for wear */
  Teuchos::ParameterList& wear = list->sublist("WEAR",false,"");

  setStringToIntegralParameter<int>("WEARLAW","None","Type of wear law",
      tuple<std::string>("None","none",
                         "Archard","archard"),
      tuple<int>(
                 INPAR::WEAR::wear_none, INPAR::WEAR::wear_none,
                 INPAR::WEAR::wear_archard, INPAR::WEAR::wear_archard),
      &wear);

  setStringToIntegralParameter<int>("WEARCOEFF_CONF","material","configuration in which wcoeff is defined",
      tuple<std::string>("material","mat",
                         "spatial","sp"),
      tuple<int>(
                 INPAR::WEAR::wear_conf_mat, INPAR::WEAR::wear_conf_mat,
                 INPAR::WEAR::wear_conf_sp, INPAR::WEAR::wear_conf_sp),
      &wear);

  setStringToIntegralParameter<int>("WEAR_SHAPEFCN","std","Type of employed set of shape functions for wear",
        tuple<std::string>("Dual", "dual",
                           "Standard", "standard", "std"),
        tuple<int>(
                INPAR::WEAR::wear_shape_dual, INPAR::WEAR::wear_shape_dual,
                INPAR::WEAR::wear_shape_standard, INPAR::WEAR::wear_shape_standard, INPAR::WEAR::wear_shape_standard),
        &wear);

  DoubleParameter("WEARCOEFF",0.0,"Wear coefficient for slave surface",&wear);
  DoubleParameter("WEARCOEFF_MASTER",0.0,"Wear coefficient for master surface",&wear);
  DoubleParameter("WEAR_TIMERATIO",1.0,"Time step ratio between wear and spatial time scale",&wear);
  DoubleParameter("SSSLIP",1.0,"Fixed slip for steady state wear",&wear);

  setStringToIntegralParameter<int>("SSWEAR","No","flag for steady state wear",
                               yesnotuple,yesnovalue,&wear);

  setStringToIntegralParameter<int>("VOLMASS_OUTPUT","No","flag for output of mass/volume in ref,mat and cur. conf.",
                               yesnotuple,yesnovalue,&wear);

  setStringToIntegralParameter<int>("BOTH_SIDED_WEAR","No","Definition of wear side",
        tuple<std::string>("No","no", "none" ,
                           "Mapping","mapping", "map",
                           "bothdiscr","bd", "sm"),
        tuple<int>(
                INPAR::WEAR::wear_slave, INPAR::WEAR::wear_slave, INPAR::WEAR::wear_slave,
                INPAR::WEAR::wear_both_map, INPAR::WEAR::wear_both_map,INPAR::WEAR::wear_both_map,
                INPAR::WEAR::wear_both_discr, INPAR::WEAR::wear_both_discr,INPAR::WEAR::wear_both_discr),
        &wear);

  setStringToIntegralParameter<int>("WEARTYPE","internal_state_expl","Definition of wear algorithm",
        tuple<std::string>("intstate_impl","is_impl", "internal_state_impl" ,
                           "intstate_expl","is_expl", "internal_state_expl",
                           "primvar","pv", "primary_variable"),
        tuple<int>(
                INPAR::WEAR::wear_intstate_impl, INPAR::WEAR::wear_intstate_impl, INPAR::WEAR::wear_intstate_impl,
                INPAR::WEAR::wear_intstate_expl, INPAR::WEAR::wear_intstate_expl, INPAR::WEAR::wear_intstate_expl,
                INPAR::WEAR::wear_primvar, INPAR::WEAR::wear_primvar, INPAR::WEAR::wear_primvar),
        &wear);

  setStringToIntegralParameter<int>("WEAR_COUPALGO","stagg","Definition of wear (ALE) coupling algorithm",
        tuple<std::string>("stagg","s",
                           "iterstagg","is",
                           "monolithic","mono"),
        tuple<int>(
                INPAR::WEAR::wear_stagg, INPAR::WEAR::wear_stagg,
                INPAR::WEAR::wear_iterstagg, INPAR::WEAR::wear_iterstagg,
                INPAR::WEAR::wear_monolithic, INPAR::WEAR::wear_monolithic),
        &wear);

  setStringToIntegralParameter<int>("WEAR_TIMESCALE","equal","Definition wear time scale compares to std. time scale",
        tuple<std::string>("equal","e",
                           "different","d"),
        tuple<int>(
                INPAR::WEAR::wear_time_equal, INPAR::WEAR::wear_time_equal,
                INPAR::WEAR::wear_time_different, INPAR::WEAR::wear_time_different),
        &wear);

  /*----------------------------------------------------------------------*/
  /* parameters for tsi contact */
  Teuchos::ParameterList& tsic = list->sublist("TSI CONTACT",false,"");

  DoubleParameter("HEATTRANSSLAVE",0.0,"Heat transfer parameter for slave side in thermal contact",&tsic);
  DoubleParameter("HEATTRANSMASTER",0.0,"Heat transfer parameter for master side in thermal contact",&tsic);

  setStringToIntegralParameter<int>("THERMOLAGMULT","Yes","Lagrange Multipliers are applied for thermo-contact",
                               yesnotuple,yesnovalue,&tsic);

  /*----------------------------------------------------------------------*/
  /* parameters for beam contact */
  Teuchos::ParameterList& beamcontact = list->sublist("BEAM CONTACT",false,"");

  setStringToIntegralParameter<int>("BEAMS_STRATEGY","None","Type of employed solving strategy",
        tuple<std::string>("None","none",
                           "Penalty", "penalty",
                           "Uzawa","uzawa"),
        tuple<int>(
                INPAR::BEAMCONTACT::bstr_none, INPAR::BEAMCONTACT::bstr_none,
                INPAR::BEAMCONTACT::bstr_penalty, INPAR::BEAMCONTACT::bstr_penalty,
                INPAR::BEAMCONTACT::bstr_uzawa, INPAR::BEAMCONTACT::bstr_uzawa),
        &beamcontact);

  setStringToIntegralParameter<int>("BEAMS_NEWGAP","No","choose between original or enhanced gapfunction",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_SEGCON","No","choose between beam contact with and without subsegment generation",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_DEBUG","No","This flag can be used for testing purposes. When it is switched on, some sanity checks are not performed!",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_INACTIVESTIFF","No","Always apply contact stiffness in first Newton step for pairs which have active in last time step",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_BTSOL","No","decide, if also the contact between beams and solids is possible",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_BTSPH","No","decide, if also the contact between beams and spheres is possible",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_SMOOTHING","None","Application of smoothed tangent field",
       tuple<std::string>("None","none",
                          "Cpp", "cpp"),
       tuple<int>(
                  INPAR::BEAMCONTACT::bsm_none,INPAR::BEAMCONTACT::bsm_none,
                  INPAR::BEAMCONTACT::bsm_cpp,INPAR::BEAMCONTACT::bsm_cpp),
       &beamcontact);

  setStringToIntegralParameter<int>("BEAMS_DAMPING","No","Application of a contact damping force",
       tuple<std::string>("No","no",
                          "Yes", "yes"),
       tuple<int>(
                  INPAR::BEAMCONTACT::bd_no,INPAR::BEAMCONTACT::bd_no,
                  INPAR::BEAMCONTACT::bd_yes,INPAR::BEAMCONTACT::bd_yes),
       &beamcontact);

  DoubleParameter("BEAMS_BTBPENALTYPARAM",0.0,"Penalty parameter for beam-to-beam point contact",&beamcontact);
  DoubleParameter("BEAMS_BTBLINEPENALTYPARAM",-1.0,"Penalty parameter per unit length for beam-to-beam line contact",&beamcontact);
  DoubleParameter("BEAMS_BTSPENALTYPARAM",0.0,"Penalty parameter for beam-to-solid contact",&beamcontact);
  DoubleParameter("BEAMS_BTSPH_PENALTYPARAM",0.0,"Penalty parameter for beam-to-rigidsphere penalty / Uzawa augmented solution strategy",&beamcontact);
  IntParameter("BEAMS_BTBUZAWAMAXSTEPS",10,"Maximum no. of Uzawa steps for Uzawa solution strategy",&beamcontact);
  DoubleParameter("BEAMS_BTBUZAWACONSTRTOL",1.0e-8,"Tolerance of constraint norm for Uzawa solution strategy",&beamcontact);
  DoubleParameter("BEAMS_DAMPINGPARAM",-1000.0,"Damping parameter for contact damping force",&beamcontact);
  DoubleParameter("BEAMS_DAMPREGPARAM1",-1000.0,"First (at gap1, with gap1>gap2) regularization parameter for contact damping force",&beamcontact);
  DoubleParameter("BEAMS_DAMPREGPARAM2",-1000.0,"Second (at gap2, with gap1>gap2) regularization parameter for contact damping force",&beamcontact);
  DoubleParameter("BEAMS_MAXDISISCALEFAC",-1.0,"Scale factor in order to limit maximal iterative displacement increment (resiudal displacement)",&beamcontact);
  DoubleParameter("BEAMS_MAXDELTADISSCALEFAC",1.0,"Scale factor in order to limit maximal displacement per time step",&beamcontact);

  DoubleParameter("BEAMS_PERPSHIFTANGLE1",-1.0,"Lower shift angle (in degrees) for penalty scaling of large-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_PERPSHIFTANGLE2",-1.0,"Upper shift angle (in degrees) for penalty scaling of large-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_PARSHIFTANGLE1",-1.0,"Lower shift angle (in degrees) for penalty scaling of small-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_PARSHIFTANGLE2",-1.0,"Upper shift angle (in degrees) for penalty scaling of small-angle-contact",&beamcontact);
  IntParameter("BEAMS_NUMINTEGRATIONINTERVAL",1,"Number of integration intervals per element",&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_PENALTYLAW","LinPen","Applied Penalty Law",
       tuple<std::string>("LinPen",
                          "QuadPen",
                          "LinNegQuadPen",
                          "LinPosQuadPen",
                          "LinPosCubPen",
                          "LinPosDoubleQuadPen",
                          "LinPosExpPen"),
       tuple<int>(
                  INPAR::BEAMCONTACT::pl_lp,
                  INPAR::BEAMCONTACT::pl_qp,
                  INPAR::BEAMCONTACT::pl_lnqp,
                  INPAR::BEAMCONTACT::pl_lpqp,
                  INPAR::BEAMCONTACT::pl_lpcp,
                  INPAR::BEAMCONTACT::pl_lpdqp,
                  INPAR::BEAMCONTACT::pl_lpep),
       &beamcontact);

  DoubleParameter("BEAMS_PENREGPARAM_G0",-1.0,"First penalty regularization parameter G0 >=0: For gap<G0 contact is active!",&beamcontact);
  DoubleParameter("BEAMS_PENREGPARAM_F0",-1.0,"Second penalty regularization parameter F0 >=0: F0 represents the force at the transition point between regularized and linear force law!",&beamcontact);
  DoubleParameter("BEAMS_PENREGPARAM_C0",-1.0,"Third penalty regularization parameter C0 >=0: C0 has different physical meanings for the different penalty laws!",&beamcontact);
  DoubleParameter("BEAMS_BASICSTIFFGAP",-1.0,"For gaps > -BEAMS_BASICSTIFFGAP, only the basic part of the contact linearization is applied!",&beamcontact);

  // enable octree search and determine type of bounding box (aabb = axis aligned, cobb = cylindrical oriented)
  setStringToIntegralParameter<int>("BEAMS_OCTREE","None","octree and bounding box type for octree search routine",
       tuple<std::string>("None","none","octree_axisaligned","octree_cylorient","octree_spherical"),
       tuple<int>(INPAR::BEAMCONTACT::boct_none,INPAR::BEAMCONTACT::boct_none,
                  INPAR::BEAMCONTACT::boct_aabb,INPAR::BEAMCONTACT::boct_cobb,INPAR::BEAMCONTACT::boct_spbb),
       &beamcontact);

  setStringToIntegralParameter<int>("BEAMS_ADDITEXT","Yes","Switch between No==multiplicative extrusion factor and Yes==additive extrusion factor",
                               yesnotuple,yesnovalue,&beamcontact);
  setNumericStringParameter("BEAMS_EXTVAL","-1.0", "extrusion value(s) of the bounding box, Depending on BEAMS_ADDITIVEEXTFAC is either additive or multiplicative. Give one or two values.",&beamcontact);
  IntParameter("BEAMS_TREEDEPTH",6,"max, tree depth of the octree",&beamcontact);
  IntParameter("BEAMS_BOXESINOCT",8,"max number of bounding boxes in any leaf octant",&beamcontact);


  /*----------------------------------------------------------------------*/
  /* parameters for potential-based beam interaction */
  Teuchos::ParameterList& beampotential = list->sublist("BEAM POTENTIAL",false,"");

  setNumericStringParameter("POT_LAW_EXPONENT","", "negative(!) exponent(s) m_i of potential law Phi(r) = sum_i (k_i * r^(-m_i)).",&beampotential);
  setNumericStringParameter("POT_LAW_PREFACTOR","", "prefactor(s) k_i of potential law Phi(r) = sum_i (k_i * r^(-m_i)).",&beampotential);
  DoubleParameter("CUTOFFRADIUS",-1.0,"cutoff radius for search of potential-based interaction pairs",&beampotential);

  setStringToIntegralParameter<int>("BEAMPOTENTIAL_TYPE","Surface","Type of potential interaction: surface (default) or volume potential",
       tuple<std::string>("Surface","surface",
                          "Volume", "volume"),
       tuple<int>(
                  INPAR::BEAMPOTENTIAL::beampot_surf,INPAR::BEAMPOTENTIAL::beampot_surf,
                  INPAR::BEAMPOTENTIAL::beampot_vol,INPAR::BEAMPOTENTIAL::beampot_vol),
       &beampotential);

  setStringToIntegralParameter<int>("BEAMPOT_BTSOL","No","decide, whether potential-based interaction between beams and solids is considered",
                               yesnotuple,yesnovalue,&beampotential);

  setStringToIntegralParameter<int>("BEAMPOT_BTSPH","No","decide, whether potential-based interaction between beams and spheres is considered",
                               yesnotuple,yesnovalue,&beampotential);

  // enable octree search and determine type of bounding box (aabb = axis aligned, spbb = spherical)
  setStringToIntegralParameter<int>("BEAMPOT_OCTREE","None","octree and bounding box type for octree search routine",
       tuple<std::string>("None","none","octree_axisaligned","octree_cylorient","octree_spherical"),
       tuple<int>(INPAR::BEAMCONTACT::boct_none,INPAR::BEAMCONTACT::boct_none,
                  INPAR::BEAMCONTACT::boct_aabb,INPAR::BEAMCONTACT::boct_cobb,INPAR::BEAMCONTACT::boct_spbb),
       &beampotential);

  IntParameter("BEAMPOT_TREEDEPTH",6,"max, tree depth of the octree",&beampotential);
  IntParameter("BEAMPOT_BOXESINOCT",8,"max number of bounding boxes in any leaf octant",&beampotential);

  /*----------------------------------------------------------------------*/
  /* parameters for semi-smooth Newton plasticity algorithm */
  Teuchos::ParameterList& iplast = list->sublist("SEMI-SMOOTH PLASTICITY",false,"");

  DoubleParameter("SEMI_SMOOTH_CPL",1.0,"Weighting factor cpl for semi-smooth PDASS",&iplast);
  DoubleParameter("STABILIZATION_S",1.0,"Stabilization factor s for semi-smooth PDASS",&iplast);

  // solver convergence test parameters for semi-smooth plasticity formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFPLASTCONSTR","And",
    "binary operator to combine plasticity constraints and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::STR::bop_and,
      INPAR::STR::bop_or),
    &iplast
    );

  setStringToIntegralParameter<int>("NORMCOMBI_DISPPLASTINCR","And",
      "binary operator to combine displacement increments and plastic flow (Delta Lp) increment values",
      tuple<std::string>(
        "And",
        "Or"),
      tuple<int>(
        INPAR::STR::bop_and,
        INPAR::STR::bop_or),
      &iplast
      );

  DoubleParameter("TOLPLASTCONSTR",1.0E-8,
                  "tolerance in the plastic constraint norm for the newton iteration",
                  &iplast);
  DoubleParameter("TOLDELTALP",1.0E-8,
                  "tolerance in the plastic flow (Delta Lp) norm for the Newton iteration",
                  &iplast);

  setStringToIntegralParameter<int>("NORMCOMBI_EASRES","And",
    "binary operator to combine EAS-residual and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::STR::bop_and,
      INPAR::STR::bop_or),
    &iplast
    );

  setStringToIntegralParameter<int>("NORMCOMBI_EASINCR","And",
      "binary operator to combine displacement increments and EAS increment values",
      tuple<std::string>(
        "And",
        "Or"),
      tuple<int>(
        INPAR::STR::bop_and,
        INPAR::STR::bop_or),
      &iplast
      );

  DoubleParameter("TOLEASRES",1.0E-8,
                  "tolerance in the EAS residual norm for the newton iteration",
                  &iplast);
  DoubleParameter("TOLEASINCR",1.0E-8,
                  "tolerance in the EAS increment norm for the Newton iteration",
                  &iplast);

  setStringToIntegralParameter<int>("DISSIPATION_MODE","pl_multiplier",
      "method to calculate the plastic dissipation",
      tuple<std::string>(
        "pl_multiplier",
        "pl_flow"),
      tuple<int>(
        INPAR::TSI::pl_multiplier,
        INPAR::TSI::pl_flow),
      &iplast
      );

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
                               //listing possible std::strings in input file in category THERMALBATH
                               tuple<std::string>("None","none",
                                                  "Uniform","uniform",
                                                  "ShearFlow","shearflow","Shearflow"),
                               //translating input std::strings into BACI input parameters
                               tuple<int>(INPAR::STATMECH::thermalbath_none,INPAR::STATMECH::thermalbath_none,
                                          INPAR::STATMECH::thermalbath_uniform,INPAR::STATMECH::thermalbath_uniform,
                                          INPAR::STATMECH::thermalbath_shearflow,INPAR::STATMECH::thermalbath_shearflow,INPAR::STATMECH::thermalbath_shearflow),
                               &statmech);
  //Reading which kind of special output should be written to files
  setStringToIntegralParameter<int>("SPECIAL_OUTPUT","None","kind of special statistical output data written into files",
                                 //listing possible std::strings in input file in category SPECIAL_OUTPUT
                                 tuple<std::string>("None","none",
                                                    "endtoend_log",
                                                    "anisotropic",
                                                    "orientationcorrelation",
                                                    "endtoend_const",
                                                    "viscoelasticity",
                                                    "networkcreep",
                                                    "networkrelax",
                                                    "networkdispfield",
                                                    "structanaly",
                                                    "octree",
                                                    "loom",
                                                    "loomelnrg",
                                                    "motassay"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::statout_none,INPAR::STATMECH::statout_none,
                                            INPAR::STATMECH::statout_endtoendlog,
                                            INPAR::STATMECH::statout_anisotropic,
                                            INPAR::STATMECH::statout_orientationcorrelation,
                                            INPAR::STATMECH::statout_endtoendconst,
                                            INPAR::STATMECH::statout_viscoelasticity,
                                            INPAR::STATMECH::statout_networkcreep,
                                            INPAR::STATMECH::statout_networkrelax,
                                            INPAR::STATMECH::statout_networkdispfield,
                                            INPAR::STATMECH::statout_structanaly,
                                            INPAR::STATMECH::statout_octree,
                                            INPAR::STATMECH::statout_loom,
                                            INPAR::STATMECH::statout_loomelnrg,
                                            INPAR::STATMECH::statout_motassay),
                                 &statmech);
  //Reading which kind of friction model should be applied
  setStringToIntegralParameter<int>("FRICTION_MODEL","none","friction model for polymer dynamics",
                                 //listing possible std::strings in input file in category FRICTION_MODEL
                                 tuple<std::string>("none",
                                                    "isotropiclumped",
                                                    "isotropicconsistent",
                                                    "anisotropicconsistent"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::frictionmodel_none,
                                            INPAR::STATMECH::frictionmodel_isotropiclumped,
                                            INPAR::STATMECH::frictionmodel_isotropicconsistent,
                                            INPAR::STATMECH::frictionmodel_anisotropicconsistent),
                                            &statmech);
  //Reading which kind of Dirichlet boundary condition should be applied
  setStringToIntegralParameter<int>("DBCTYPE","std","Dirichlet BC type applied",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("none",
                                                    "std",
                                                    "shearfixed",
                                                    "shearfixeddel",
                                                    "sheartrans",
                                                    "pinnodes" ,
                                                    "affineshear",
                                                    "affinesheardel",
                                                    "movablesupport1d"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::dbctype_none,
                                            INPAR::STATMECH::dbctype_std,
                                            INPAR::STATMECH::dbctype_shearfixed,
                                            INPAR::STATMECH::dbctype_shearfixeddel,
                                            INPAR::STATMECH::dbctype_sheartrans,
                                            INPAR::STATMECH::dbctype_pinnodes,
                                            INPAR::STATMECH::dbctype_affineshear,
                                            INPAR::STATMECH::dbctype_affinesheardel,
                                            INPAR::STATMECH::dbctype_movablesupport1d),
                                            &statmech);
  //Reading which kind of Dirichlet boundary condition should be applied
  setStringToIntegralParameter<int>("NBCTYPE","std","Neumann BC type applied",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("std",
                                                    "constcreep",
                                                    "randompointforce"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::nbctype_std,
                                            INPAR::STATMECH::nbctype_constcreep,
                                            INPAR::STATMECH::nbctype_randompointforce),
                                            &statmech);
  //Reading which kind of biopolymer network will be simulated
  setStringToIntegralParameter<int>("NETWORKTYPE","std","Network type simulated",
                                 //listing possible std::strings in input file in category FILAMENTMODEL
                                 tuple<std::string>("std",
                                                    "loom",
                                                    "motassay"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::networktype_std,
                                            INPAR::STATMECH::networktype_casimir,
                                            INPAR::STATMECH::networktype_motassay),
                                            &statmech);
  //Reading which kind of linker model should be applied
  setStringToIntegralParameter<int>("LINKERMODEL","none","Linker model applied in Statmech simulations",
                                 //listing possible std::strings in input file in category LINKERMODEL
                                 tuple<std::string>("none",
                                                    "std",
                                                    "stdintpol",
                                                    "bellseq",
                                                    "bellseqintpol",
                                                    "active",
                                                    "activeintpol",
                                                    "myosinthick"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::linkermodel_none,
                                            INPAR::STATMECH::linkermodel_std,
                                            INPAR::STATMECH::linkermodel_stdintpol,
                                            INPAR::STATMECH::linkermodel_bellseq,
                                            INPAR::STATMECH::linkermodel_bellseqintpol,
                                            INPAR::STATMECH::linkermodel_active,
                                            INPAR::STATMECH::linkermodel_activeintpol,
                                            INPAR::STATMECH::linkermodel_myosinthick),
                                            &statmech);
  setStringToIntegralParameter<int>("PLANELINKERMOTION","No",
                                 "Plane Brownian Motion of linkers",
                                 yesnotuple,yesnovalue,&statmech);
  setStringToIntegralParameter<int>("CROSSBRIDGEMODEL","No",
                                 "swinging cross bridge model for active linkers",
                                 yesnotuple,yesnovalue,&statmech);
  //Reading which kind of filament model should be applied
  setStringToIntegralParameter<int>("FILAMENTMODEL","std","Filament model applied in Statmech simulations",
                                 //listing possible std::strings in input file in category FILAMENTMODEL
                                 tuple<std::string>("std",
                                                    "helical"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::filamentmodel_std,
                                            INPAR::STATMECH::filamentmodel_helical),
                                            &statmech);
  //Reading which kind of search routine is used
  setStringToIntegralParameter<int>("BINDINGSITESEARCH","volpart","binding site search method",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("volpart",
                                                    "binning",
                                                    "octree"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(INPAR::STATMECH::bsstype_volpart,
                                            INPAR::STATMECH::bsstype_binning,
                                            INPAR::STATMECH::bsstype_octree),
                                            &statmech);

  //time after which writing of statistical output is started
  DoubleParameter("STARTTIMEOUT",0.0,"Time after which writing of statistical output is started",&statmech);
  // Time values at which certain actions are carried out
  setNumericStringParameter("ACTIONTIME","-1.0","Points in time (corresponding to ACTIONDT values), where certain actions are carried out. Order: [t_equilib; t_ktswitch; ...; t_act]",&statmech);
  // time step sizes corresponding to ACTIONTIME
  setNumericStringParameter("ACTIONDT","-1.0","Time step sizes corresponding to ACTIONTIME values.",&statmech);
  // index controlling the start of BC application (see ACTIONTIME)
  IntParameter("BCTIMEINDEX",-1,"Integer refers to the n-th entry of ACTIONTIME. States beginning of BC application. Starts counting at '1' !",&statmech);
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("FILAMENTPOLARITY","No","toggles filament polarity",yesnotuple,yesnovalue,&statmech);
  //Rise per monomer in the actin double helix according to Howard, p. 125
  setNumericStringParameter("RISEPERBSPOT","0.00277","rise per monomer in the actin one-start helix",&statmech);
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
  // chemical rate for contractile conformation change
  DoubleParameter("K_ACT_SHORT_start",0.0,"rate",&statmech);
  // chemical rate for contractile conformation change
  DoubleParameter("K_ACT_SHORT_end",0.0,"rate",&statmech);
  // chemical rate for extensional conformation change
  DoubleParameter("K_ACT_LONG_start",0.0,"rate",&statmech);
  // chemical rate for extensional conformation change
  DoubleParameter("K_ACT_LONG_end",0.0,"rate",&statmech);
  // active linker stroke distance
  DoubleParameter("STROKEDISTANCE",0.005,"active linker stroke distance",&statmech);
  // scaling factor for linker length changes
  DoubleParameter("LINKERSCALEFACTOR",0.0,"Scaling factor for active linker length changes. No effect on BEAM3CL elements!",&statmech);
  //displacement in the reaction coordinate used in Bell's euqations
  DoubleParameter("DELTABELLSEQ",0.0,"displacement in the reaction coordinate used in Bell's eqation (<0.0 -> catch bond, >0 -> std bond, 0 == no Bell's equation ",&statmech);
  // active linker fraction
  DoubleParameter("ACTIVELINKERFRACTION",0.0,"Fraction of linkers that show active behavior", &statmech);
  // cycle time
  DoubleParameter("ACTIVELINKERCYCLE",0.04,"duration of a work cycle of an active linker",&statmech);
  // time fraction during which no bonding is possible due to the linker being in its recovery state
  DoubleParameter("ACTIVERECOVERYFRACTION",0.95,"fraction of ACTIVELINKERCYCLE during which the linker recovers (i.e. is unbound)",&statmech);
  //number of overall crosslink molecules in the boundary volume
  IntParameter("N_crosslink",0,"number of crosslinkers for switching on- and off-rates; if molecule diffusion model is used: number of crosslink molecules",&statmech);
  //number of overall crosslink molecules in the boundary volume
  IntParameter("INITOCCUPIEDBSPOTS",0,"binding spots occupied by (singly-bound) crosslinkers before the first time step",&statmech);
  //number of filaments used as substrate for motility assay setups
  IntParameter("NUMSUBSTRATEFIL",0,"Number of filaments used as substrate filaments",&statmech);
  //number by which the number of crosslinkers is reduced.
  IntParameter("REDUCECROSSLINKSBY",0,"number of crosslinker elements by which the overall number of crosslinker is reduced.",&statmech);
  IntParameter("BSPOTINTERVAL",1,"determines every n-th binding spot available for crosslinking",&statmech);
  //Reading double parameter for crosslinker protein mean length
  DoubleParameter("R_LINK",0.0,"Mean distance between two nodes connected by a crosslinker",&statmech);
  //Absolute value of difference between maximal/minimal and mean cross linker length
  DoubleParameter("DeltaR_LINK",0.0,"Absolute value of difference between maximal/minimal and mean cross linker length",&statmech);
  // Three values representing the size of the periodic box in each spatial direction
  setNumericStringParameter("PERIODLENGTH","0.0 0.0 0.0", "Values representing the size of the periodic box in each spatial direction",&statmech);
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
  //Torsional stiffness of crosslinkers
  DoubleParameter("KTOR1_LINK",0.0,"Torsional stiffness at binding spots between Truss linkers and filaments",&statmech);
  //Torsional stiffness of crosslinkers
  DoubleParameter("KTOR2_LINK",0.0,"Torsional stiffness between tangents of filaments",&statmech);
  //Parameter for PTC according to Cyron,Wall (2011):Numerical method for the simulation of the Brownian dynamics of rod-like microstructures with three dimensional nonlinear beam elements
  DoubleParameter("CTRANSPTC0",0.0,"PTC factor for translational DOF in first iteration step",&statmech);
  //Parameter for PTC according to Cyron,Wall (2011):Numerical method for the simulation of the Brownian dynamics of rod-like microstructures with three dimensional nonlinear beam elements
  DoubleParameter("CROTPTC0",0.145,"PTC factor for rotational DOF in first iteration step",&statmech);
  //Parameter for PTC (rigid spherical particle)
  DoubleParameter("CSPHEREPTC0",0.0,"PTC factor for DOFs of rigid spherical particle in first iteration step",&statmech);
  //Parameter for PTC according to Cyron,Wall (2011):Numerical method for the simulation of the Brownian dynamics of rod-like microstructures with three dimensional nonlinear beam elements
  DoubleParameter("ALPHAPTC",6.0,"exponent of power law for reduction of PTC factor",&statmech);
  //Number of iterations after which PTC is turned off
  IntParameter("MAXITERPTC",5,"Number of iterations after which PTC is turned off! 0 entails standard Newton-Raphson.",&statmech);
  // fraction of initial residual below which PTC is turned off
  DoubleParameter("RESLOWPTC",0.001,"fraction of initial residual below which PTC is turned off!",&statmech);
  //Makes filaments and crosslinkers be plotted by that factor thicker than they are acutally
  DoubleParameter("PlotFactorThick",0.0,"Makes filaments and crosslinkers be plotted by that factor thicker than they are acutally",&statmech);
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("CHECKORIENT","No","If chosen crosslinkers are set only after check of orientation of linked filaments", yesnotuple, yesnovalue, &statmech);
  //Gmsh Output switch
  setStringToIntegralParameter<int>("GMSHOUTPUT","No","If chosen gmsh output is generated.", yesnotuple,yesnovalue,&statmech);
  // toggling Gmsh Output for structure detection
  setStringToIntegralParameter<int>("GMSHNETSTRUCT","No","If chosen, special gmsh visualization for network structure types is generated.", yesnotuple,yesnovalue,&statmech);
  //Number of time steps between two special outputs written
  IntParameter("OUTPUTINTERVALS",1,"Number of time steps between two special outputs written",&statmech);
  //Number of time steps between two gmsh outputs written
  IntParameter("GMSHOUTINTERVALS",100,"Number of time steps between two gmsh outputs written",&statmech);
  //Reading direction of oscillatory motion that DBC nodes are subjected to (we need this when using periodic BCs)
  IntParameter("DBCDISPDIR",0,"Global spatial direction of oscillatory motion by Dirichlet BCs",&statmech);
  //Reading time curve number for Dirichlet boundary conditions
  IntParameter("CURVENUMBER",0,"Specifies Time Curve number of imposed Dirichlet BCs",&statmech);
  //Reading time curve number for Neumann boundary conditions
  IntParameter("NBCCURVENUMBER",0,"Specifies Time Curve number of Neumann BCs",&statmech);
  // amplitude of Neumann boundary force
  DoubleParameter("NBCFORCEAMP",0.0,"constant creep force in NBCs",&statmech);
  // amplitude of Neumann boundary force
  IntParameter("NUMNBCNODES",0,"Number of nodes to which Neumann point forces are applied sequentially.",&statmech);
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
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("FIXEDSEED","No","If chosen fixed seed for random numbers in each time step is applied", yesnotuple,yesnovalue,&statmech);
  //time interval in which random numbers are constant
  DoubleParameter("RANDNUMTIMEINT",-1.0,"Within this time interval the random numbers remain constant. -1.0 means no prescribed time interval.'",&statmech);
  //cutoff for random forces, which determines the maximal value
  DoubleParameter("MAXRANDFORCE",-1.0,"Any random force beyond MAXRANDFORCE*(standard dev.) will be omitted and redrawn. -1.0 means no bounds.'",&statmech);
  // vector of constant background velocity
  setNumericStringParameter("BACKGROUNDVELOCITY","0.0 0.0 0.0", "Vector of constant background velocity (x,y,z)",&statmech);
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& tdyn = list->sublist("THERMAL DYNAMIC",false,"");

  setStringToIntegralParameter<int>("DYNAMICTYP","OneStepTheta",
    "type of time integration control",
    tuple<std::string>(
      "Statics",
      "OneStepTheta",
      "GenAlpha",
      "ExplicitEuler"),
    tuple<int>(
      INPAR::THR::dyna_statics,
      INPAR::THR::dyna_onesteptheta,
      INPAR::THR::dyna_genalpha,
      INPAR::THR::dyna_expleuler),
    &tdyn
    );

  // output type
  IntParameter("RESEVRYGLOB",1,"save temperature and other global quantities every RESEVRYGLOB steps",&tdyn);
  IntParameter("RESEVRYERGY",0,"write system energies every requested step",&tdyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&tdyn);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
    "Initial Field for thermal problem",
    tuple<std::string>(
      "zero_field",
      "field_by_function",
      "field_by_condition"),
    tuple<int>(
      INPAR::THR::initfield_zero_field,
      INPAR::THR::initfield_field_by_function,
      INPAR::THR::initfield_field_by_condition),
    &tdyn
    );

  IntParameter("INITFUNCNO",-1,"function number for thermal initial field",&tdyn);

  // Time loop control
  DoubleParameter("TIMESTEP",0.05,"time step size",&tdyn);
  IntParameter("NUMSTEP",200,"maximum number of steps",&tdyn);
  DoubleParameter("MAXTIME",5.0,"maximum time",&tdyn);

  // Iterationparameters
  DoubleParameter("TOLTEMP",1.0E-10,
    "tolerance in the temperature norm of the Newton iteration",
    &tdyn
    );

  setStringToIntegralParameter<int>("NORM_TEMP","Abs",
    "type of norm for temperature convergence check",
    tuple<std::string>(
      "Abs",
      "Rel",
      "Mix"),
    tuple<int>(
      INPAR::THR::convnorm_abs,
      INPAR::THR::convnorm_rel,
      INPAR::THR::convnorm_mix),
    &tdyn
    );

  DoubleParameter("TOLRES",1.0E-08,
    "tolerance in the residual norm for the Newton iteration",
    &tdyn
    );

  setStringToIntegralParameter<int>("NORM_RESF","Abs",
    "type of norm for residual convergence check",
    tuple<std::string>(
      "Abs",
      "Rel",
      "Mix"),
    tuple<int>(
      INPAR::THR::convnorm_abs,
      INPAR::THR::convnorm_rel,
      INPAR::THR::convnorm_mix),
    &tdyn
    );

  setStringToIntegralParameter<int>("NORMCOMBI_RESFTEMP","And",
    "binary operator to combine temperature and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::THR::bop_and,
      INPAR::THR::bop_or),
    &tdyn
    );

  IntParameter("MAXITER",50,
    "maximum number of iterations allowed for Newton-Raphson iteration before failure",
    &tdyn
    );

  IntParameter("MINITER",0,
    "minimum number of iterations to be done within Newton-Raphson loop",
    &tdyn
    );

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
    &tdyn
    );

  setStringToIntegralParameter<int>("DIVERCONT","stop","What to do with time integration when Newton-Raphson iteration failed",
    tuple<std::string>(
    "stop",
    "continue",
    "halve_step",
    "repeat_step",
    "repeat_simulation"),
    tuple<int>(
    INPAR::THR::divcont_stop,
    INPAR::THR::divcont_continue,
    INPAR::THR::divcont_halve_step,
    INPAR::THR::divcont_repeat_step,
    INPAR::THR::divcont_repeat_simulation),
  &tdyn);

  setStringToIntegralParameter<int>("NLNSOL","fullnewton","Nonlinear solution technique",
    tuple<std::string>(
      "vague",
      "fullnewton"),
    tuple<int>(
      INPAR::THR::soltech_vague,
      INPAR::THR::soltech_newtonfull),
    &tdyn
    );

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
    &tdyn
    );

  // convergence criteria solver adaptivity
  setStringToIntegralParameter<int>("ADAPTCONV","No",
    "Switch on adaptive control of linear solver tolerance for nonlinear solution",
    yesnotuple,yesnovalue,&tdyn
    );
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&tdyn);

  setStringToIntegralParameter<int>("LUMPCAPA","No",
    "Lump the capacity matrix for explicit time integration",
    yesnotuple,yesnovalue,&tdyn
    );

  // number of linear solver used for thermal problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for thermal problems",&tdyn);

  // where the geometry comes from
  setStringToIntegralParameter<int>(
    "GEOMETRY","full",
    "How the geometry is specified",
    tuple<std::string>(
      "full",
      "box",
      "file"),
    tuple<int>(
      INPAR::geometry_full,
      INPAR::geometry_box,
      INPAR::geometry_file),
    &tdyn);

  /*----------------------------------------------------------------------*/
  /* parameters for generalised-alpha thermal integrator */
  Teuchos::ParameterList& tgenalpha = tdyn.sublist("GENALPHA",false,"");

  setStringToIntegralParameter<int>("GENAVG","TrLike",
    "mid-average type of internal forces",
    tuple<std::string>(
      "Vague",
      "ImrLike",
      "TrLike"),
    tuple<int>(
      INPAR::THR::midavg_vague,
      INPAR::THR::midavg_imrlike,
      INPAR::THR::midavg_trlike),
    &tgenalpha
    );

  // default values correspond to midpoint-rule
  DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&tgenalpha);
  DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0.5,1)",&tgenalpha);
  DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0.5,1)",&tgenalpha);

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

  // coupling strategy for (partitioned and monolithic) TSI solvers
  setStringToIntegralParameter<int>("COUPALGO","tsi_monolithic",
    "Coupling strategies for TSI solvers",
    tuple<std::string>(
      "tsi_oneway",
      "tsi_sequstagg",
      "tsi_iterstagg",
      "tsi_iterstagg_aitken",
      "tsi_iterstagg_aitkenirons",
      "tsi_iterstagg_fixedrelax",
      "tsi_monolithic"),
    tuple<int>(
      INPAR::TSI::OneWay,
      INPAR::TSI::SequStagg,
      INPAR::TSI::IterStagg,
      INPAR::TSI::IterStaggAitken,
      INPAR::TSI::IterStaggAitkenIrons,
      INPAR::TSI::IterStaggFixedRel,
      INPAR::TSI::Monolithic),
    &tsidyn
    );

  BoolParameter("MATCHINGGRID","Yes","is matching grid",&tsidyn);

  // coupling strategy for BACI-INCA coupling
  setStringToIntegralParameter<int>("TFSI_COUPALGO","tfsi",
    "Coupling strategies for BACI-INCA coupling (TFSI)",
    tuple<std::string>(
      "tfsi",
      "fsi",
      "conj_heat_transfer",
      "no_inca_fsi"),
    tuple<int>(
        INPAR::TSI::TFSI,
        INPAR::TSI::FSI,
        INPAR::TSI::ConjHeatTransfer,
        INPAR::TSI::NoIncaFSI),
    &tsidyn
    );

  // scaling factor for AeroTFSI problems when length unit other than SI [m] is used
  setStringToIntegralParameter<double>("TFSI_length_unit","m",
    "Used unit for extension in the structural model in AeroTFSI",
    tuple<std::string>(
      "m",
      "cm",
      "mm"),
    tuple<double>(
      1.0,
      100.0,
      1000.0),
    &tsidyn
    );

  // scaling factor for AeroTFSI problems when time unit other than SI [s] is used
  setStringToIntegralParameter<double>("TFSI_time_unit","s",
    "Used unit for time in the structural model in AeroTFSI",
    tuple<std::string>(
      "s",
      "ms",
      "mikros"),
    tuple<double>(
      1.0,
      1.0e3,
      1.0e6),
    &tsidyn
    );

  // output type
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&tsidyn);

  // time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&tsidyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&tsidyn);
  DoubleParameter("TIMESTEP",0.05,"time step size dt",&tsidyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&tsidyn);
  IntParameter("ITEMIN",1,"minimal number of iterations over fields",&tsidyn);
  IntParameter("UPRES",1,"increment for writing solution",&tsidyn);

  setStringToIntegralParameter<int>("NORM_INC","Abs",
    "type of norm for convergence check of primary variables in TSI",
    tuple<std::string>(
      "Abs",
      "Rel",
      "Mix"),
    tuple<int>(
      INPAR::TSI::convnorm_abs,
      INPAR::TSI::convnorm_rel,
      INPAR::TSI::convnorm_mix),
    &tsidyn
    );

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic TSI */
  Teuchos::ParameterList& tsidynmono = tsidyn.sublist("MONOLITHIC",false,
    "Monolithic Thermo Structure Interaction\n"
    "Dynamic section for monolithic TSI"
    );

  // convergence tolerance of tsi residual
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of TSI",&tsidynmono);
  // Iterationparameters
  DoubleParameter("TOLINC",1.0e-6,"tolerance for convergence check of TSI-increment in monolithic TSI",&tsidynmono);

  setStringToIntegralParameter<int>("NORM_RESF","Abs",
    "type of norm for residual convergence check",
    tuple<std::string>(
      "Abs",
      "Rel",
      "Mix"),
    tuple<int>(
      INPAR::TSI::convnorm_abs,
      INPAR::TSI::convnorm_rel,
      INPAR::TSI::convnorm_mix),
    &tsidynmono
    );

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","Coupl_And_Singl",
    "binary operator to combine primary variables and residual force values",
    tuple<std::string>(
      "And",
      "Or",
      "Coupl_Or_Singl",
      "Coupl_And_Singl",
      "And_Singl",
      "Or_Singl"),
    tuple<int>(
      INPAR::TSI::bop_and,
      INPAR::TSI::bop_or,
      INPAR::TSI::bop_coupl_or_singl,
      INPAR::TSI::bop_coupl_and_singl,
      INPAR::TSI::bop_and_singl,
      INPAR::TSI::bop_or_singl),
    &tsidynmono
    );

  setStringToIntegralParameter<int>("ITERNORM","Rms",
    "type of norm to be applied to residuals",
    tuple<std::string>(
      "L1",
      "L1_Scaled",
      "L2",
      "Rms",
      "Inf"),
    tuple<int>(
      INPAR::TSI::norm_l1,
      INPAR::TSI::norm_l1_scaled,
      INPAR::TSI::norm_l2,
      INPAR::TSI::norm_rms,
      INPAR::TSI::norm_inf),
    &tsidynmono
    );

  setStringToIntegralParameter<int>("NLNSOL","fullnewton",
    "Nonlinear solution technique",
    tuple<std::string>(
      "fullnewton",
      "ptc"),
    tuple<int>(
      INPAR::TSI::soltech_newtonfull,
      INPAR::TSI::soltech_ptc),
    &tsidynmono
    );

  DoubleParameter("PTCDT",0.1,
    "pseudo time step for pseudo-transient continuation (PTC) stabilised Newton procedure",
    &tsidynmono
    );

  // number of linear solver used for monolithic TSI
  IntParameter("LINEAR_SOLVER",-1,
    "number of linear solver used for monolithic TSI problems",
    &tsidynmono
    );

  // convergence criteria adaptivity of monolithic TSI solver
  setStringToIntegralParameter<int>("ADAPTCONV","No",
    "Switch on adaptive control of linear solver tolerance for nonlinear solution",
    yesnotuple,
    yesnovalue,
    &tsidynmono
    );
  DoubleParameter("ADAPTCONV_BETTER",0.1,
    "The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",
    &tsidynmono
    );

  setStringToIntegralParameter<int>("INFNORMSCALING","yes",
    "Scale blocks of matrix with row infnorm?",
    yesnotuple,
    yesnovalue,
    &tsidynmono
    );

  // merge TSI block matrix to enable use of direct solver in monolithic TSI
  // default: "No", i.e. use block matrix
  BoolParameter("MERGE_TSI_BLOCK_MATRIX","No","Merge TSI block matrix",&tsidynmono);

  // in case of monolithic TSI nodal values (displacements, temperatures and
  // reaction forces) at fix points of the body can be calculated
  // default: "No", i.e. nothing is calculated
  BoolParameter("CALC_NECKING_TSI_VALUES","No",
    "Calculate nodal values for evaluation and validation of necking",
    &tsidynmono
    );

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned TSI */
  Teuchos::ParameterList& tsidynpart = tsidyn.sublist(
      "PARTITIONED",false,
      "Partitioned Thermo Structure Interaction\n"
      "Dynamic section for partitioned TSI"
       );

  // decide in partitioned TSI which one-way coupling or predictor should be used
  setStringToIntegralParameter<int>("COUPVARIABLE","Displacement",
    "Coupling variable",
    tuple<std::string>(
      "Displacement",
      "Temperature"),
    tuple<int>(0,1),
    &tsidynpart
    );

  // Solver parameter for relaxation of iterative staggered partitioned TSI
  DoubleParameter("MAXOMEGA",0.0,"largest omega allowed for Aitken relaxation (0.0 means no constraint)",&tsidynpart);
  DoubleParameter("FIXEDOMEGA",1.0,"fixed relaxation parameter",&tsidynpart);

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of outer iteraiton within partitioned TSI",&tsidynpart);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& poroelastdyn = list->sublist(
   "POROELASTICITY DYNAMIC",false,
   "Poroelasticity"
   );

  // Coupling strategy for (monolithic) porous media solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","poro_monolithic",
                              "Coupling strategies for poroelasticity solvers",
                              tuple<std::string>(
                                 "poro_partitioned",
                                 "poro_monolithic",
                                 "poro_monolithicstructuresplit",
                                 "poro_monolithicfluidsplit",
                                 "poro_monolithicnopenetrationsplit"
                                ),
                              tuple<int>(
                                INPAR::POROELAST::Partitioned,
                                INPAR::POROELAST::Monolithic,
                                INPAR::POROELAST::Monolithic_structuresplit,
                                INPAR::POROELAST::Monolithic_fluidsplit,
                                INPAR::POROELAST::Monolithic_nopenetrationsplit
                                ),
                              &poroelastdyn);

  // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq approximation)
  setStringToIntegralParameter<int>("PHYSICAL_TYPE","Poro",
                                    "Physical Type of Porofluid",
                                    tuple<std::string>(
                                      "Poro",
                                      "Poro_P1",
                                      "Poro_P2"
                                      ),
                                    tuple<int>(
                                      INPAR::FLUID::poro,
                                      INPAR::FLUID::poro_p1,
                                      INPAR::FLUID::poro_p2),
                                    &poroelastdyn);

  // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq approximation)
  setStringToIntegralParameter<int>("TIME_DISTYPE_CONTI","pressure",
                                    "type of time discretization for continuity equation",
                                    tuple<std::string>(
                                      "pressure",
                                      "pres",
                                      "porosity"
                                      ),
                                    tuple<int>(
                                      INPAR::POROELAST::pressure,
                                      INPAR::POROELAST::pressure,
                                      INPAR::POROELAST::porosity),
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
  DoubleParameter("TOLRES_GLOBAL",1e-8,"tolerance in the residual norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLINC_GLOBAL",1e-8,"tolerance in the increment norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLRES_DISP",1e-8,"tolerance in the residual norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLINC_DISP",1e-8,"tolerance in the increment norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLRES_PORO",1e-8,"tolerance in the residual norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLINC_PORO",1e-8,"tolerance in the increment norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLRES_VEL",1e-8,"tolerance in the residual norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLINC_VEL",1e-8,"tolerance in the increment norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLRES_PRES",1e-8,"tolerance in the residual norm for the Newton iteration",&poroelastdyn);
  DoubleParameter("TOLINC_PRES",1e-8,"tolerance in the increment norm for the Newton iteration",&poroelastdyn);

  setStringToIntegralParameter<int>("NORM_INC","AbsSingleFields","type of norm for primary variables convergence check",
                               tuple<std::string>(
                                   "AbsGlobal",
                                   "AbsSingleFields"
                                 ),
                               tuple<int>(
                                   INPAR::POROELAST::convnorm_abs_global,
                                   INPAR::POROELAST::convnorm_abs_singlefields
                                 ),
                               &poroelastdyn);

  setStringToIntegralParameter<int>("NORM_RESF","AbsSingleFields","type of norm for residual convergence check",
                                 tuple<std::string>(
                                     "AbsGlobal",
                                     "AbsSingleFields"
                                   ),
                                 tuple<int>(
                                   INPAR::POROELAST::convnorm_abs_global,
                                   INPAR::POROELAST::convnorm_abs_singlefields
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

  setStringToIntegralParameter<int>("VECTORNORM_RESF","L2",
                                "type of norm to be applied to residuals",
                                tuple<std::string>(
                                  "L1",
                                  "L1_Scaled",
                                  "L2",
                                  "Rms",
                                  "Inf"),
                                tuple<int>(
                                  INPAR::POROELAST::norm_l1,
                                  INPAR::POROELAST::norm_l1_scaled,
                                  INPAR::POROELAST::norm_l2,
                                  INPAR::POROELAST::norm_rms,
                                  INPAR::POROELAST::norm_inf),
                                &poroelastdyn
                                );

  setStringToIntegralParameter<int>("VECTORNORM_INC","L2",
                              "type of norm to be applied to residuals",
                              tuple<std::string>(
                                "L1",
                                "L1_Scaled",
                                "L2",
                                "Rms",
                                "Inf"),
                              tuple<int>(
                                INPAR::POROELAST::norm_l1,
                                INPAR::POROELAST::norm_l1_scaled,
                                INPAR::POROELAST::norm_l2,
                                INPAR::POROELAST::norm_rms,
                                INPAR::POROELAST::norm_inf),
                              &poroelastdyn
                              );

  setStringToIntegralParameter<int>("SECONDORDER","Yes",
                               "Second order coupling at the interface.",
                               yesnotuple,yesnovalue,&poroelastdyn);

  setStringToIntegralParameter<int>("CONTIPARTINT","No",
                               "Partial integration of porosity gradient in continuity equation",
                               yesnotuple,yesnovalue,&poroelastdyn);

  setStringToIntegralParameter<int>("CONTACTNOPEN","No",
                               "No-Penetration Condition on active contact surface in case of poro contact problem!",
                               yesnotuple,yesnovalue,&poroelastdyn);

  BoolParameter("MATCHINGGRID","Yes","is matching grid",&poroelastdyn);

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
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&poroscatradyn);
  IntParameter("ITEMIN",1,"minimal number of iterations over fields",&poroscatradyn);

  // Iterationparameters
  DoubleParameter("TOLRES_GLOBAL",1e-8,"tolerance in the residual norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLINC_GLOBAL",1e-8,"tolerance in the increment norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLRES_DISP",1e-8,"tolerance in the residual norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLINC_DISP",1e-8,"tolerance in the increment norm for the Newton iteration",&poroscatradyn);
//  DoubleParameter("TOLRES_PORO",1e-8,"tolerance in the residual norm for the Newton iteration",&poroscatradyn);
//  DoubleParameter("TOLINC_PORO",1e-8,"tolerance in the increment norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLRES_VEL",1e-8,"tolerance in the residual norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLINC_VEL",1e-8,"tolerance in the increment norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLRES_PRES",1e-8,"tolerance in the residual norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLINC_PRES",1e-8,"tolerance in the increment norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLRES_SCALAR",1e-8,"tolerance in the residual norm for the Newton iteration",&poroscatradyn);
  DoubleParameter("TOLINC_SCALAR",1e-8,"tolerance in the increment norm for the Newton iteration",&poroscatradyn);

  setStringToIntegralParameter<int>("NORM_INC","AbsSingleFields","type of norm for primary variables convergence check",
                               tuple<std::string>(
                                   "AbsGlobal",
                                   "AbsSingleFields"
                                 ),
                               tuple<int>(
                                   INPAR::POROELAST::convnorm_abs_global,
                                   INPAR::POROELAST::convnorm_abs_singlefields
                                 ),
                               &poroscatradyn);

  setStringToIntegralParameter<int>("NORM_RESF","AbsSingleFields","type of norm for residual convergence check",
                                 tuple<std::string>(
                                     "AbsGlobal",
                                     "AbsSingleFields"
                                   ),
                                 tuple<int>(
                                   INPAR::POROELAST::convnorm_abs_global,
                                   INPAR::POROELAST::convnorm_abs_singlefields
                                   ),
                                 &poroscatradyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                               tuple<std::string>(
                                     "And",
                                     "Or"),
                                     tuple<int>(
                                       INPAR::POROELAST::bop_and,
                                       INPAR::POROELAST::bop_or),
                               &poroscatradyn);

  setStringToIntegralParameter<int>("VECTORNORM_RESF","L2",
                                "type of norm to be applied to residuals",
                                tuple<std::string>(
                                  "L1",
                                  "L1_Scaled",
                                  "L2",
                                  "Rms",
                                  "Inf"),
                                tuple<int>(
                                  INPAR::POROELAST::norm_l1,
                                  INPAR::POROELAST::norm_l1_scaled,
                                  INPAR::POROELAST::norm_l2,
                                  INPAR::POROELAST::norm_rms,
                                  INPAR::POROELAST::norm_inf),
                                &poroscatradyn
                                );

  setStringToIntegralParameter<int>("VECTORNORM_INC","L2",
                              "type of norm to be applied to residuals",
                              tuple<std::string>(
                                "L1",
                                "L1_Scaled",
                                "L2",
                                "Rms",
                                "Inf"),
                              tuple<int>(
                                INPAR::POROELAST::norm_l1,
                                INPAR::POROELAST::norm_l1_scaled,
                                INPAR::POROELAST::norm_l2,
                                INPAR::POROELAST::norm_rms,
                                INPAR::POROELAST::norm_inf),
                              &poroscatradyn
                              );

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for monolithic poroscatra problems",&poroscatradyn);

  // Coupling strategy for poroscatra solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","solid_to_scatra",
                              "Coupling strategies for poroscatra solvers",
                              tuple<std::string>(
                                "monolithic",
                                "scatra_to_solid",
                                "solid_to_scatra",
                                "two_way"
                                ),
                              tuple<int>(
                                INPAR::PORO_SCATRA::Monolithic,
                                INPAR::PORO_SCATRA::Part_ScatraToPoro,
                                INPAR::PORO_SCATRA::Part_PoroToScatra,
                                INPAR::PORO_SCATRA::Part_TwoWay
                                ),
                              &poroscatradyn);

  /*----------------------------------------------------------------------*/
  /* parameters for SSI */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidyn = list->sublist(
   "SSI CONTROL",false,
   "Control paramters for scatra structure interaction"
   );

  // Output type
  DoubleParameter("RESTARTEVRYTIME",0,"write restart possibility every RESTARTEVRY steps",&ssidyn);
  IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&ssidyn);
  // Time loop control
  IntParameter("NUMSTEP",200,"maximum number of Timesteps",&ssidyn);
  DoubleParameter("MAXTIME",1000.0,"total simulation time",&ssidyn);
  DoubleParameter("TIMESTEP",-1,"time step size dt",&ssidyn);
  BoolParameter("DIFFTIMESTEPSIZE","No","use different step size for scatra and solid",&ssidyn);
  DoubleParameter("UPRESTIME",0,"increment for writing solution",&ssidyn);
  IntParameter("UPRES",1,"increment for writing solution",&ssidyn);
  IntParameter("ITEMAX",10,"maximum number of iterations over fields",&ssidyn);
  BoolParameter("SCATRA_FROM_RESTART_FILE","No","read scatra result from restart files (use option 'restartfromfile' during execution of baci)",&ssidyn);
  StringParameter("SCATRA_FILENAME","nil","Control-file name for reading scatra results in SSI",&ssidyn);


  // Coupling strategy for SSI solvers
  setStringToIntegralParameter<int>(
                              "COUPALGO","ssi_IterStagg",
                              "Coupling strategies for SSI solvers",
                              tuple<std::string>(
                                "ssi_OneWay_ScatraToSolid",
                                "ssi_OneWay_SolidToScatra",
//                                "ssi_SequStagg_ScatraToSolid",
//                                "ssi_SequStagg_SolidToScatra",
                                "ssi_IterStagg",
                                "ssi_IterStaggFixedRel_ScatraToSolid",
                                "ssi_IterStaggFixedRel_SolidToScatra",
                                "ssi_IterStaggAitken_ScatraToSolid",
                                "ssi_IterStaggAitken_SolidToScatra"
                                ),
                              tuple<int>(
                                INPAR::SSI::ssi_OneWay_ScatraToSolid,
                                INPAR::SSI::ssi_OneWay_SolidToScatra,
//                                INPAR::SSI::ssi_SequStagg_ScatraToSolid,
//                                INPAR::SSI::ssi_SequStagg_SolidToScatra,
                                INPAR::SSI::ssi_IterStagg,
                                INPAR::SSI::ssi_IterStaggFixedRel_ScatraToSolid,
                                INPAR::SSI::ssi_IterStaggFixedRel_SolidToScatra,
                                INPAR::SSI::ssi_IterStaggAitken_ScatraToSolid,
                                INPAR::SSI::ssi_IterStaggAitken_SolidToScatra
                                ),
                              &ssidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned SSI */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidynpart = ssidyn.sublist(
      "PARTITIONED",false,
      "Partitioned Structure Scalar Interaction\n"
      "Control section for partitioned SSI"
       );

  // Solver parameter for relaxation of iterative staggered partitioned SSI
  DoubleParameter("MAXOMEGA",10.0,"largest omega allowed for Aitken relaxation",&ssidynpart);
  DoubleParameter("MINOMEGA",0.1,"smallest omega allowed for Aitken relaxation",&ssidynpart);
  DoubleParameter("STARTOMEGA",1.0,"fixed relaxation parameter",&ssidynpart);

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL",1e-6,"tolerance for convergence check of outer iteraiton within partitioned TSI",&ssidynpart);

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
                                    tuple<std::string>(
                                      "Incompressible",
                                      "Artificial_compressibility",
                                      "Varying_density",
                                      "Loma",
                                      "Boussinesq",
                                      "Topology_optimization",
                                      "Stokes",
                                      "Oseen"
                                      ),
                                    tuple<int>(
                                      INPAR::FLUID::incompressible,
                                      INPAR::FLUID::artcomp,
                                      INPAR::FLUID::varying_density,
                                      INPAR::FLUID::loma,
                                      INPAR::FLUID::boussinesq,
                                      INPAR::FLUID::topopt,
                                      INPAR::FLUID::stokes,
                                      INPAR::FLUID::oseen),
                                    &fdyn);

  // number of linear solver used for fluid problem
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for fluid dynamics",&fdyn);

  // number of linear solver used for fluid problem (former fluid pressure solver for SIMPLER preconditioning with fluid)
  IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for fluid dynamics (ONLY NECESSARY FOR BlockGaussSeidel solver block within fluid mehstying case any more!!!!)",&fdyn);

  // Set ML-solver number for smooting of residual-based calculated wallshearstress via plain aggregation.
  IntParameter("WSS_ML_AGR_SOLVER",-1,"Set ML-solver number for smoothing of residual-based calculated wallshearstress via plain aggregation.",&fdyn);

  setStringToIntegralParameter<int>(
    "TIMEINTEGR","One_Step_Theta",
    "Time Integration Scheme",
    tuple<std::string>(
      "Stationary",
      "Np_Gen_Alpha",
      "Af_Gen_Alpha",
      "One_Step_Theta",
      "BDF2"),
    tuple<int>(
      INPAR::FLUID::timeint_stationary,
      INPAR::FLUID::timeint_npgenalpha,
      INPAR::FLUID::timeint_afgenalpha,
      INPAR::FLUID::timeint_one_step_theta,
      INPAR::FLUID::timeint_bdf2),
    &fdyn);

  setStringToIntegralParameter<int>(
    "OST_CONT_PRESS","Cont_normal_Press_normal",
    "One step theta option for time discretization of continuity eq. and pressure",
    tuple<std::string>(
      "Cont_normal_Press_normal",
      "Cont_impl_Press_normal",
      "Cont_impl_Press_impl"),
    tuple<int>(
      INPAR::FLUID::Cont_normal_Press_normal,
      INPAR::FLUID::Cont_impl_Press_normal,
      INPAR::FLUID::Cont_impl_Press_impl),
    &fdyn);

  setStringToIntegralParameter<int>(
    "GEOMETRY","full",
    "How the geometry is specified",
    tuple<std::string>(
      "full",
      "box",
      "file"),
    tuple<int>(
      INPAR::geometry_full,
      INPAR::geometry_box,
      INPAR::geometry_file),
    &fdyn);

  setStringToIntegralParameter<int>(
    "NONLINITER","fixed_point_like",
    "Nonlinear iteration scheme",
    tuple<std::string>(
      "fixed_point_like",
      "Newton"
      ),
    tuple<int>(
      INPAR::FLUID::fixed_point_like,
      INPAR::FLUID::Newton),
    &fdyn);

  setStringToIntegralParameter<int>("PREDICTOR","steady_state",
                                    "Predictor for first guess in nonlinear iteration",
                                    tuple<std::string>(
                                      "steady_state",
                                      "zero_acceleration",
                                      "constant_acceleration",
                                      "constant_increment",
                                      "explicit_second_order_midpoint",
                                      "TangVel"
                                      ),
                                    tuple<int>(1,2,3,4,5,6),
                                    &fdyn);

  setStringToIntegralParameter<int>("CONVCHECK","L_2_norm",
                               "norm for convergence check",
                               tuple<std::string>(
                                 //"L_infinity_norm",
                                 //"L_1_norm",
                                 "L_2_norm"
                                 //"L_2_norm_without_residual_at_itemax"
                                 ),
                               tuple<std::string>(
                                 //"use max norm (ccarat)",
                                 //"use abs. norm (ccarat)",
                                 "compute L2 errors of increments (relative) and residuals (absolute)"
                                 //"same as L_2_norm, only no residual norm is computed if itemax is reached (speedup for turbulence calculations, startup phase)"
                                 ),
                               tuple<int>(
                                 //INPAR::FLUID::fncc_Linf,
                                 //INPAR::FLUID::fncc_L1,
                                 INPAR::FLUID::fncc_L2
                                 //INPAR::FLUID::fncc_L2_wo_res
                                 ),
                               &fdyn);

  BoolParameter("INCONSISTENT_RESIDUAL","No","do not evaluate residual after solution has converged (->faster)",&fdyn);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string,11> name;
    Teuchos::Tuple<int,11> label;
    name[ 0] = "zero_field";                             label[ 0] = INPAR::FLUID::initfield_zero_field;
    name[ 1] = "field_by_function";                      label[ 1] = INPAR::FLUID::initfield_field_by_function;
    name[ 2] = "disturbed_field_from_function";          label[ 2] = INPAR::FLUID::initfield_disturbed_field_from_function;
    name[ 3] = "FLAME_VORTEX_INTERACTION";               label[ 3] = INPAR::FLUID::initfield_flame_vortex_interaction;
    name[ 4] = "BELTRAMI-FLOW";                          label[ 4] = INPAR::FLUID::initfield_beltrami_flow;
    name[ 5] = "KIM-MOIN-FLOW";                          label[ 5] = INPAR::FLUID::initfield_kim_moin_flow;
    name[ 6] = "BOCHEV-TEST";                            label[ 6] = INPAR::FLUID::initfield_bochev_test;
    name[ 7] = "hit_comte_bellot_corrsin_initial_field"; label[ 7] = INPAR::FLUID::initfield_hit_comte_bellot_corrsin;
    name[ 8] = "forced_hit_simple_algebraic_spectrum";   label[ 8] = INPAR::FLUID::initfield_forced_hit_simple_algebraic_spectrum;
    name[ 9] = "forced_hit_numeric_spectrum";            label[ 9] = INPAR::FLUID::initfield_forced_hit_numeric_spectrum;
    name[10] = "forced_hit_passive";                     label[10] = INPAR::FLUID::initfield_passive_hit_const_input;

    setStringToIntegralParameter<int>(
        "INITIALFIELD",
        "zero_field",
        "Initial field for fluid problem",
        name,
        label,
        &fdyn);
  }

  IntParameter("OSEENFIELDFUNCNO",-1,"function number of Oseen advective field",&fdyn);

  BoolParameter("LIFTDRAG","No","Calculate lift and drag forces along specified boundary",&fdyn);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(0,1),
                               &fdyn);

  setStringToIntegralParameter<int>("NONLINEARBC",
                               "no",
                               "Flag to activate check for potential nonlinear boundary conditions",
                               tuple<std::string>(
                                 "no",
                                 "yes"),
                               tuple<std::string>(
                                 "no nonlinear boundary conditions",
                                 "nonlinear boundary conditions might be present"),
                               tuple<int>(0,1),
                               &fdyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
                                  tuple<std::string>(
                                    "no",
                                    "Condensed_Smat",
                                    "Condensed_Bmat",
                                    "Condensed_Bmat_merged"),
                                  tuple<int>(
                                      INPAR::FLUID::no_meshtying,
                                      INPAR::FLUID::condensed_smat,
                                      INPAR::FLUID::condensed_bmat,
                                      INPAR::FLUID::condensed_bmat_merged),
                                  &fdyn);

  setStringToIntegralParameter<int>("GRIDVEL", "BE", "scheme for determination of gridvelocity from displacements",
                                  tuple<std::string>(
                                    "BE",
                                    "BDF2",
                                    "OST"),
                                  tuple<int>(
                                      INPAR::FLUID::BE,
                                      INPAR::FLUID::BDF2,
                                      INPAR::FLUID::OST),
                                  &fdyn);

  BoolParameter("ALLDOFCOUPLED","Yes","all dof (incl. pressure) are coupled",&fdyn);

  {
   Teuchos::Tuple<std::string,17> name;
   Teuchos::Tuple<int,17> label;

   name[ 0] = "no";                             label[ 0] = INPAR::FLUID::no_error_calculation;
   name[ 1] = "beltrami_flow";                  label[ 1] = INPAR::FLUID::beltrami_flow;
   name[ 2] = "channel2D";                      label[ 2] = INPAR::FLUID::channel2D;
   name[ 3] = "gravitation";                    label[ 3] = INPAR::FLUID::gravitation;
   name[ 4] = "shear_flow";                     label[ 4] = INPAR::FLUID::shear_flow;
   name[ 5] = "jeffery_hamel_flow";             label[ 5] = INPAR::FLUID::jeffery_hamel_flow;
   name[ 6] = "byfunct1";                       label[ 6] = INPAR::FLUID::byfunct1;
   name[ 7] = "beltrami_stat_stokes";           label[ 7] = INPAR::FLUID::beltrami_stat_stokes;
   name[ 8] = "beltrami_stat_navier_stokes";    label[ 8] = INPAR::FLUID::beltrami_stat_navier_stokes;
   name[ 9] = "beltrami_instat_stokes";         label[ 9] = INPAR::FLUID::beltrami_instat_stokes;
   name[10] = "beltrami_instat_navier_stokes";  label[10] = INPAR::FLUID::beltrami_instat_navier_stokes;
   name[11] = "kimmoin_stat_stokes";            label[11] = INPAR::FLUID::kimmoin_stat_stokes;
   name[12] = "kimmoin_stat_navier_stokes";     label[12] = INPAR::FLUID::kimmoin_stat_navier_stokes;
   name[13] = "kimmoin_instat_stokes";          label[13] = INPAR::FLUID::kimmoin_instat_stokes;
   name[14] = "kimmoin_instat_navier_stokes";   label[14] = INPAR::FLUID::kimmoin_instat_navier_stokes;
   name[15] = "fsi_fluid_pusher";               label[15] = INPAR::FLUID::fsi_fluid_pusher;
   name[16] = "topopt_channel";                 label[16] = INPAR::FLUID::topoptchannel;

   setStringToIntegralParameter<int>("CALCERROR","no",
                                     "Flag to (de)activate error calculations",
                                     name,
                                     label,
                                     &fdyn);
  }

  setStringToIntegralParameter<int>("SIMPLER","no",
                               "Switch on SIMPLE family of solvers, only works with block preconditioners like CheapSIMPLE!",
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
  BoolParameter("COMPUTE_DIVU","No","Compute divergence of velocity field at the element center",&fdyn);
  BoolParameter("COMPUTE_EKIN","No","Compute kinetic energy at the end of each time step and write it to file.",&fdyn);
  BoolParameter("NEW_OST","No","Solve the Navier-Stokes equation with the new One Step Theta algorithm",&fdyn);  //TODO: To be removed.
  IntParameter("UPRES",1,"Increment for writing solution",&fdyn);
  IntParameter("RESTARTEVRY",20,"Increment for writing restart",&fdyn);
  IntParameter("NUMSTEP",1,"Total number of Timesteps",&fdyn);
  IntParameter("STEADYSTEP",-1,"steady state check every step",&fdyn);
  IntParameter("NUMSTASTEPS",0,"Number of Steps for Starting Scheme",&fdyn);
  IntParameter("STARTFUNCNO",-1,"Function for Initial Starting Field",&fdyn);
  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&fdyn);
  IntParameter("INITSTATITEMAX",5,"max number of nonlinear iterations for initial stationary solution",&fdyn);
  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&fdyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&fdyn);
  DoubleParameter("ALPHA_M",1.0,"Time integration factor",&fdyn);
  DoubleParameter("ALPHA_F",1.0,"Time integration factor",&fdyn);
  DoubleParameter("GAMMA",1.0,"Time integration factor",&fdyn);
  DoubleParameter("THETA",0.66,"Time integration factor",&fdyn);

  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&fdyn);
  DoubleParameter("STEADYTOL",1e-6,"Tolerance for steady state check",&fdyn);
  DoubleParameter("START_THETA",1.0,"Time integration factor for starting scheme",&fdyn);

  DoubleParameter("CFL_NUMBER",-1.0,"CFL number for adaptive time step",&fdyn);
  IntParameter("FREEZE_ADAPTIVE_DT_AT",1000000,"keep time step constant after this step, otherwise turbulence statistics sampling is not consistent",&fdyn);

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

  IntParameter("VELGRAD_PROJ_SOLVER",-1,"Number of linear solver used for L2 projection",&fdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_stab = fdyn.sublist("RESIDUAL-BASED STABILIZATION",false,"");

  // this parameter defines various stabilized methods
  setStringToIntegralParameter<int>("STABTYPE",
                               "residual_based",
                               "Apply (un)stabilized fluid formulation",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "residual_based",
                                 "edge_based",
                                 "pressure_projection"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> inf-sup stable elements required!",
                                 "Use a residual-based stabilization or, more generally, a stabilization \nbased on the concept of the residual-based variational multiscale method...\nExpecting additional input",
                                 "Use an edge-based stabilization, especially for XFEM",
                                 "Element/cell based polynomial pressure projection, see Dohrmann/Bochev 2004, IJNMF")  ,
                               tuple<int>(
                                   INPAR::FLUID::stabtype_nostab,
                                   INPAR::FLUID::stabtype_residualbased,
                                   INPAR::FLUID::stabtype_edgebased,
                                   INPAR::FLUID::stabtype_pressureprojection),
                               &fdyn_stab);

  BoolParameter("INCONSISTENT","No","residual based without second derivatives (i.e. only consistent for tau->0, but faster)",&fdyn_stab);

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

  BoolParameter("PSPG","Yes","Flag to (de)activate PSPG stabilization.",&fdyn_stab);
  BoolParameter("SUPG","Yes","Flag to (de)activate SUPG stabilization.",&fdyn_stab);
  BoolParameter("GRAD_DIV","Yes","Flag to (de)activate grad-div term.",&fdyn_stab);

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

  {
    // this parameter selects the tau definition applied
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string,16> name;
    Teuchos::Tuple<int,16> label;
    name[ 0] = "Taylor_Hughes_Zarins";                        label[ 0] = INPAR::FLUID::tau_taylor_hughes_zarins;
    name[ 1] = "Taylor_Hughes_Zarins_wo_dt";                  label[ 1] = INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt;
    name[ 2] = "Taylor_Hughes_Zarins_Whiting_Jansen";         label[ 2] = INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen;
    name[ 3] = "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt";   label[ 3] = INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt;
    name[ 4] = "Taylor_Hughes_Zarins_scaled";                 label[ 4] = INPAR::FLUID::tau_taylor_hughes_zarins_scaled;
    name[ 5] = "Taylor_Hughes_Zarins_scaled_wo_dt";           label[ 5] = INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt;
    name[ 6] = "Franca_Barrenechea_Valentin_Frey_Wall";       label[ 6] = INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall;
    name[ 7] = "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt"; label[ 7] = INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt;
    name[ 8] = "Shakib_Hughes_Codina";                        label[ 8] = INPAR::FLUID::tau_shakib_hughes_codina;
    name[ 9] = "Shakib_Hughes_Codina_wo_dt";                  label[ 9] = INPAR::FLUID::tau_shakib_hughes_codina_wo_dt;
    name[10] = "Codina";                                      label[10] = INPAR::FLUID::tau_codina;
    name[11] = "Codina_wo_dt";                                label[11] = INPAR::FLUID::tau_codina_wo_dt;
    name[12] = "Codina_convscaled";                           label[12] = INPAR::FLUID::tau_codina_convscaled;
    name[13] = "Franca_Madureira_Valentin_Badia_Codina";      label[13] = INPAR::FLUID::tau_franca_madureira_valentin_badia_codina;
    name[14] = "Franca_Madureira_Valentin_Badia_Codina_wo_dt";label[14] = INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt;
    name[15] = "Hughes_Franca_Balestra_wo_dt";                label[15] = INPAR::FLUID::tau_hughes_franca_balestra_wo_dt;

    setStringToIntegralParameter<int>(
        "DEFINITION_TAU",
        "Franca_Barrenechea_Valentin_Frey_Wall",
        "Definition of tau_M and Tau_C",
        name,
        label,
        &fdyn_stab);
  }

  // this parameter selects the characteristic element length for tau_Mu for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_U",
                               "streamlength",
                               "Characteristic element length for tau_Mu",
                               tuple<std::string>(
                                 "streamlength",
                                 "volume_equivalent_diameter",
                                 "root_of_volume"),
                               tuple<int>(
                                   INPAR::FLUID::streamlength_u,
                                   INPAR::FLUID::volume_equivalent_diameter_u,
                                   INPAR::FLUID::root_of_volume_u),
                               &fdyn_stab);

  // this parameter selects the characteristic element length for tau_Mp and tau_C for
  // all stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_PC",
                               "volume_equivalent_diameter",
                               "Characteristic element length for tau_Mp/tau_C",
                               tuple<std::string>(
                                 "streamlength",
                                 "volume_equivalent_diameter",
                                 "root_of_volume"),
                               tuple<int>(
                                   INPAR::FLUID::streamlength_pc,
                                   INPAR::FLUID::volume_equivalent_diameter_pc,
                                   INPAR::FLUID::root_of_volume_pc),
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
  BoolParameter("LOMA_CONTI_SUPG","No","Flag to (de)activate SUPG stabilization in loma continuity equation.",&fdyn_stab);

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
  Teuchos::ParameterList& fdyn_edge_based_stab = fdyn.sublist("EDGE-BASED STABILIZATION",false,"");

  //! Flag to (de)activate edge-based (EOS) pressure stabilization
  setStringToIntegralParameter<int>("EOS_PRES",
                               "none",
                               "Flag to (de)activate pressure edge-based stabilization.",
                               tuple<std::string>(
                                 "none",
                                 "std_eos",
                                 "xfem_gp"),
                               tuple<std::string>(
                                 "do not use pressure edge-based stabilization",
                                 "use pressure edge-based stabilization as standard edge-based stabilization on the entire domain",
                                 "use pressure edge-based stabilization as xfem ghost-penalty stabilization just around cut elements"),
                               tuple<int>(
                                 INPAR::FLUID::EOS_PRES_none,       // no pressure edge-based stabilization
                                 INPAR::FLUID::EOS_PRES_std_eos,    // pressure edge-based stabilization on the entire domain
                                 INPAR::FLUID::EOS_PRES_xfem_gp     // pressure edge-based stabilization as ghost penalty around cut elements
                               ),
                               &fdyn_edge_based_stab);

  //! Flag to (de)activate edge-based (EOS) convective streamline stabilization
  setStringToIntegralParameter<int>("EOS_CONV_STREAM",
                               "none",
                               "Flag to (de)activate convective streamline edge-based stabilization.",
                               tuple<std::string>(
                                 "none",
                                 "std_eos",
                                 "xfem_gp"),
                               tuple<std::string>(
                                 "do not use convective streamline edge-based stabilization",
                                 "use convective streamline edge-based stabilization as standard edge-based stabilization on the entire domain",
                                 "use convective streamline edge-based stabilization as xfem ghost-penalty stabilization just around cut elements"),
                               tuple<int>(
                                 INPAR::FLUID::EOS_CONV_STREAM_none,       // no convective streamline edge-based stabilization
                                 INPAR::FLUID::EOS_CONV_STREAM_std_eos,    // convective streamline edge-based stabilization on the entire domain
                                 INPAR::FLUID::EOS_CONV_STREAM_xfem_gp     // pressure edge-based stabilization as ghost penalty around cut elements
                               ),
                               &fdyn_edge_based_stab);

  //! Flag to (de)activate edge-based (EOS) convective crosswind stabilization
  setStringToIntegralParameter<int>("EOS_CONV_CROSS",
                               "none",
                               "Flag to (de)activate convective crosswind edge-based stabilization.",
                               tuple<std::string>(
                                 "none",
                                 "std_eos",
                                 "xfem_gp"),
                               tuple<std::string>(
                                 "do not use convective crosswind edge-based stabilization",
                                 "use convective crosswind edge-based stabilization as standard edge-based stabilization on the entire domain",
                                 "use convective crosswind edge-based stabilization as xfem ghost-penalty stabilization just around cut elements"),
                               tuple<int>(
                                 INPAR::FLUID::EOS_CONV_CROSS_none,       // no convective crosswind edge-based stabilization
                                 INPAR::FLUID::EOS_CONV_CROSS_std_eos,    // convective crosswind edge-based stabilization on the entire domain
                                 INPAR::FLUID::EOS_CONV_CROSS_xfem_gp     // convective crosswind edge-based stabilization as ghost penalty around cut elements
                               ),
                               &fdyn_edge_based_stab);

  //! Flag to (de)activate edge-based (EOS) divergence stabilization
  setStringToIntegralParameter<int>("EOS_DIV",
                               "none",
                               "Flag to (de)activate divergence edge-based stabilization.",
                               tuple<std::string>(
                                 "none",
                                 "vel_jump_std_eos",
                                 "vel_jump_xfem_gp",
                                 "div_jump_std_eos",
                                 "div_jump_xfem_gp"),
                               tuple<std::string>(
                                 "do not use divergence edge-based stabilization",
                                 "divergence edge-based stabilization based on velocity jump on the entire domain",
                                 "divergence edge-based stabilization based on divergence jump just around cut elements",
                                 "divergence edge-based stabilization based on velocity jump on the entire domain",
                                 "divergence edge-based stabilization based on divergence jump just around cut elements"),
                               tuple<int>(
                                 INPAR::FLUID::EOS_DIV_none,                       // no convective edge-based stabilization
                                 INPAR::FLUID::EOS_DIV_vel_jump_std_eos,           // streamline convective edge-based stabilization
                                 INPAR::FLUID::EOS_DIV_vel_jump_xfem_gp,           // streamline convective edge-based stabilization
                                 INPAR::FLUID::EOS_DIV_div_jump_std_eos,           // crosswind convective edge-based stabilization
                                 INPAR::FLUID::EOS_DIV_div_jump_xfem_gp            // crosswind convective edge-based stabilization
                               ),
                               &fdyn_edge_based_stab);

  //! special least-squares condition for pseudo 2D examples where pressure level is determined via Krylov-projection
  BoolParameter("PRES_KRYLOV_2Dz","No","residual based without second derivatives (i.e. only consistent for tau->0, but faster)",&fdyn_edge_based_stab);

  //! this parameter selects the definition of Edge-based stabilization parameter
  setStringToIntegralParameter<int>("EOS_DEFINITION_TAU",
                                    "Burman_Hansbo_DAngelo_Zunino",
                                    "Definition of stabilization parameter for edge-based stabilization",
                                    tuple<std::string>(
                                    "Burman_Fernandez_Hansbo",
                                    "Burman_Fernandez_Hansbo_wo_dt",
                                    "Braack_Burman_John_Lube",
                                    "Braack_Burman_John_Lube_wo_divjump",
                                    "Franca_Barrenechea_Valentin_Wall",
                                    "Burman_Fernandez",
                                    "Burman_Hansbo_DAngelo_Zunino",
                                    "Burman_Hansbo_DAngelo_Zunino_wo_dt",
                                    "Burman",
                                    "Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling",
                                    "tau_not_defined"
                                    ),
                                    tuple<std::string>(
                                      "definition of burman_fernandez_hansbo",
                                      "definition of burman_fernandez_hansbo for stationary problems",
                                      "definition of braack_burman_john_lube",
                                      "definition of braack_burman_john_lube without explicit inclusion of divergence jump",
                                      "definition of tau_franca_barrenechea_valentin_wall",
                                      "definition of EOS_tau_burman_fernandez",
                                      "definition of EOS_tau_burman_hansbo_dangelo_zunino",
                                      "definition of EOS_tau_burman_hansbo_dangelo_zunino for stationary problems",
                                      "definition of EOS_tau_burman",
                                      "definition of EOS_tau related to residual-based stabilization",
                                      "no chosen definition"
                                      ),
                                    tuple<int>(
                                      INPAR::FLUID::EOS_tau_burman_fernandez_hansbo,
                                      INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt,
                                      INPAR::FLUID::EOS_tau_braack_burman_john_lube,
                                      INPAR::FLUID::EOS_tau_braack_burman_john_lube_wo_divjump,
                                      INPAR::FLUID::EOS_tau_franca_barrenechea_valentin_wall,
                                      INPAR::FLUID::EOS_tau_burman_fernandez,
                                      INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino,
                                      INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt,
                                      INPAR::FLUID::EOS_tau_burman,
                                      INPAR::FLUID::EOS_tau_Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling,
                                      INPAR::FLUID::EOS_tau_not_defined) ,
                                    &fdyn_edge_based_stab);

  //! this parameter selects how the element length of Edge-based stabilization is defined
  setStringToIntegralParameter<int>("EOS_H_DEFINITION",
                                    "EOS_he_max_diameter_to_opp_surf",
                                    "Definition of element length for edge-based stabilization",
                                    tuple<std::string>(
                                    "EOS_he_max_diameter_to_opp_surf",
                                    "EOS_he_max_dist_to_opp_surf",
                                    "EOS_he_surf_with_max_diameter",
                                    "EOS_hk_max_diameter",
                                    "EOS_he_surf_diameter",
                                    "EOS_he_vol_eq_diameter"),
                                    tuple<std::string>(
                                      "take the maximal (nsd-1)D diameter of faces that connect the internal face to its opposite faces",
                                      "take the maximal 1D distance along 1D edge to opposite surface for both parent elements",
                                      "take the maximal (nsd-1)D face diameter of all faces for both parent elements",
                                      "maximal nD diameter of the neighboring elements",
                                      "maximal (n-1)D diameter of the internal face/edge",
                                      "take the maximal volume eqivalent diameter of adjecent elements"
                                      ),
                                    tuple<int>(
                                    INPAR::FLUID::EOS_he_max_diameter_to_opp_surf,
                                    INPAR::FLUID::EOS_he_max_dist_to_opp_surf,
                                    INPAR::FLUID::EOS_he_surf_with_max_diameter,
                                    INPAR::FLUID::EOS_hk_max_diameter,
                                    INPAR::FLUID::EOS_he_surf_diameter,
                                    INPAR::FLUID::EOS_he_vol_eq_diameter) ,
                                    &fdyn_edge_based_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_porostab = fdyn.sublist("POROUS-FLOW STABILIZATION",false,"");

  // this parameter defines various stabilized methods
  setStringToIntegralParameter<int>("STABTYPE",
                               "residual_based",
                               "Apply (un)stabilized fluid formulation",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "residual_based",
                                 "edge_based"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> inf-sup stable elements required!",
                                 "Use a residual-based stabilization or, more generally, a stabilization \nbased on the concept of the residual-based variational multiscale method...\nExpecting additional input",
                                 "Use an edge-based stabilization, especially for XFEM")  ,
                               tuple<int>(
                                   INPAR::FLUID::stabtype_nostab,
                                   INPAR::FLUID::stabtype_residualbased,
                                   INPAR::FLUID::stabtype_edgebased),
                               &fdyn_porostab);

  BoolParameter("INCONSISTENT","No","residual based without second derivatives (i.e. only consistent for tau->0, but faster)",&fdyn_porostab);

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
                               &fdyn_porostab);

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
                               &fdyn_porostab);

  BoolParameter("PSPG","Yes","Flag to (de)activate PSPG stabilization.",&fdyn_porostab);
  BoolParameter("SUPG","Yes","Flag to (de)activate SUPG stabilization.",&fdyn_porostab);
  BoolParameter("GRAD_DIV","Yes","Flag to (de)activate grad-div term.",&fdyn_porostab);

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
                               &fdyn_porostab);

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
                               &fdyn_porostab);

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
                               &fdyn_porostab);

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
                               &fdyn_porostab);

  // this parameter selects the tau definition applied
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
                                 "Franca_Madureira_Valentin_Badia_Codina_wo_dt"),
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
                                   INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt),
                               &fdyn_porostab);

  // this parameter selects the characteristic element length for tau_Mu for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_U",
                               "streamlength",
                               "Characteristic element length for tau_Mu",
                               tuple<std::string>(
                                 "streamlength",
                                 "volume_equivalent_diameter",
                                 "root_of_volume"),
                               tuple<int>(
                                   INPAR::FLUID::streamlength_u,
                                   INPAR::FLUID::volume_equivalent_diameter_u,
                                   INPAR::FLUID::root_of_volume_u),
                               &fdyn_porostab);

  // this parameter selects the characteristic element length for tau_Mp and tau_C for
  // all stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_PC",
                               "volume_equivalent_diameter",
                               "Characteristic element length for tau_Mp/tau_C",
                               tuple<std::string>(
                                 "streamlength",
                                 "volume_equivalent_diameter",
                                 "root_of_volume"),
                               tuple<int>(
                                   INPAR::FLUID::streamlength_pc,
                                   INPAR::FLUID::volume_equivalent_diameter_pc,
                                   INPAR::FLUID::root_of_volume_pc),
                               &fdyn_porostab);

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
                               &fdyn_porostab);

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
                               &fdyn_porostab);

  // these parameters active additional terms in loma continuity equation
  // which might be identified as SUPG-/cross- and Reynolds-stress term
  BoolParameter("LOMA_CONTI_SUPG","No","Flag to (de)activate SUPG stabilization in loma continuity equation.",&fdyn_porostab);

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
                               &fdyn_porostab);

  setStringToIntegralParameter<int>("LOMA_CONTI_REYNOLDS_STRESS",
                               "no_reynolds",
                               "Flag to (de)activate Reynolds-stress term loma continuity equation-> residual-based VMM.",
                               tuple<std::string>(
                                 "no_reynolds",
                                 "yes_reynolds",
                                 "reynolds_rhs"
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
                               &fdyn_porostab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_turbu = fdyn.sublist("TURBULENCE MODEL",false,"");

  //----------------------------------------------------------------------
  // modeling strategies
  //----------------------------------------------------------------------

  setStringToIntegralParameter<int>(
    "TURBULENCE_APPROACH",
    "DNS_OR_RESVMM_LES",
    "There are several options to deal with turbulent flows.",
    tuple<std::string>(
      "DNS_OR_RESVMM_LES",
      "CLASSICAL_LES"),
    tuple<std::string>(
      "Try to solve flow as an underresolved DNS.\nMind that your stabilisation already acts as a kind of turbulence model!",
      "Perform a classical Large Eddy Simulation adding \naddititional turbulent viscosity. This may be based on various physical models."),
    tuple<int>(0,1),
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
      "Multifractal_Subgrid_Scales",
      "Vreman",
      "Dynamic_Vreman"),
    tuple<std::string>(
      "If classical LES is our turbulence approach, this is a contradiction and should cause a dserror.",
      "Classical constant coefficient Smagorinsky model. Be careful if you \nhave a wall bounded flow domain!",
      "Use an exponential damping function for the turbulent viscosity \nclose to the wall. This is only implemented for a channel geometry of \nheight 2 in y direction. The viscous lengthscale l_tau is \nrequired as additional input.",
      "The solution is filtered and by comparison of the filtered \nvelocity field with the real solution, the Smagorinsky constant is \nestimated in each step --- mind that this procedure includes \nan averaging in the xz plane, hence this implementation will only work \nfor a channel flow.",
      "Scale Similarity Model coherent with the variational multiscale formulation",
      "Scale Similarity Model according to liu, meneveau, katz",
      "Multifractal Subgrid-Scale Modeling based on the work of burton",
      "Vremans constant model",
      "Dynamic Vreman model according to You and Moin (2007)"),
    tuple<int>(0,1,2,3,4,5,6,7,8),
    &fdyn_turbu);

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

  //----------------------------------------------------------------------
  // turbulence specific output and statistics
  //----------------------------------------------------------------------

  IntParameter("SAMPLING_START",10000000,"Time step after when sampling shall be started",&fdyn_turbu);
  IntParameter("SAMPLING_STOP",1,"Time step when sampling shall be stopped",&fdyn_turbu);
  IntParameter("DUMPING_PERIOD",1,"Period of time steps after which statistical data shall be dumped",&fdyn_turbu);

  BoolParameter("SUBGRID_DISSIPATION","No","Flag to (de)activate estimation of subgrid-scale dissipation (only for seclected flows).",&fdyn_turbu);

  BoolParameter("OUTMEAN","No","Flag to (de)activate averaged paraview output",&fdyn_turbu);

  BoolParameter("TURBMODEL_LS","Yes","Flag to (de)activate turbulence model in level-set equation",&fdyn_turbu);

  //----------------------------------------------------------------------
  // turbulent flow problem and general characteristics
  //----------------------------------------------------------------------

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    // Otherwise BACI DEBUG version will crash during runtime!
    Teuchos::Tuple<std::string,20> name;
    Teuchos::Tuple<int,20> label;
    name[ 0] = "no";                                             label[ 0] = 0;
    name[ 1] = "time_averaging";                                 label[ 1] = 1;
    name[ 2] = "channel_flow_of_height_2";                       label[ 2] = 2;
    name[ 3] = "lid_driven_cavity";                              label[ 3] = 3;
    name[ 4] = "backward_facing_step";                           label[ 4] = 4;
    name[ 5] = "square_cylinder";                                label[ 5] = 5;
    name[ 6] = "square_cylinder_nurbs";                          label[ 6] = 6;
    name[ 7] = "rotating_circular_cylinder_nurbs";               label[ 7] = 7;
    name[ 8] = "rotating_circular_cylinder_nurbs_scatra";        label[ 8] = 8;
    name[ 9] = "loma_channel_flow_of_height_2";                  label[ 9] = 9;
    name[10] = "loma_lid_driven_cavity";                         label[10] = 10;
    name[11] = "loma_backward_facing_step";                      label[11] = 11;
    name[12] = "combust_oracles";                                label[12] = 12;
    name[13] = "bubbly_channel_flow";                            label[13] = 13;
    name[14] = "scatra_channel_flow_of_height_2";                label[14] = 14;
    name[15] = "decaying_homogeneous_isotropic_turbulence";      label[15] = 15;
    name[16] = "forced_homogeneous_isotropic_turbulence";        label[16] = 16;
    name[17] = "scatra_forced_homogeneous_isotropic_turbulence"; label[17] = 17;
    name[18] = "periodic_hill";                                  label[18] = 18;
    name[19] = "blood_fda_flow";                                 label[19] = 19;

    Teuchos::Tuple<std::string,20> description;
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
    description[15]="For this flow, all statistical data could be averaged in \nthe in all homogeneous directions  --- it is essentially a statistically zero dimensional flow.";
    description[16]="For this flow, all statistical data could be averaged in \nthe in all homogeneous directions  --- it is essentially a statistically zero dimensional flow.";
    description[17]="For this flow, all statistical data could be averaged in \nthe in all homogeneous directions  --- it is essentially a statistically zero dimensional flow.";
    description[18]="For this flow, statistical data is evaluated on various lines, averaged over time and z.";
    description[19]="For this flow, statistical data is evaluated on various planes.";

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
      "yz"           ,
      "xyz"          ),
    tuple<std::string>(
      "no homogeneous directions available, averaging is restricted to time averaging",
      "average along x-direction"                                                     ,
      "average along y-direction"                                                     ,
      "average along z-direction"                                                     ,
      "Wall normal direction is z, average in x and y direction"                      ,
      "Wall normal direction is y, average in x and z direction (standard case)"      ,
      "Wall normal direction is x, average in y and z direction"                      ,
      "averageing in all directions"                                                 ),
    tuple<int>(0,1,2,3,4,5,6,7),
    &fdyn_turbu);

  //---------------------------------------
  // further problem-specific parameters

  // CHANNEL FLOW
  //--------------

  DoubleParameter(
    "CHAN_AMPL_INIT_DIST",
    0.1,
    "Max. amplitude of the random disturbance in percent of the initial value in mean flow direction.",
    &fdyn_turbu);

  setStringToIntegralParameter<int>("FORCING_TYPE","linear_compensation_from_intermediate_spectrum","forcing strategy",
                               tuple<std::string>(
                                 "linear_compensation_from_intermediate_spectrum",
                                 "fixed_power_input"
                                 ),
                               tuple<int>(
                                   INPAR::FLUID::linear_compensation_from_intermediate_spectrum,
                                   INPAR::FLUID::fixed_power_input
                                   ),
                               &fdyn_turbu);

  IntParameter("CHA_NUMSUBDIVISIONS",5,"Number of homogenious sampling planes in element",&fdyn_turbu);

  // HIT
  //--------------

  IntParameter(
    "FORCING_TIME_STEPS",
    0,
    "Number of time steps during which forcing is applied. Decaying homogeneous isotropic turbulence only.",
    &fdyn_turbu);

  DoubleParameter(
    "THRESHOLD_WAVENUMBER",
    0.0,
    "Forcing is only applied to wave numbers lower or equal than the given threshold wave number.",
    &fdyn_turbu);

  DoubleParameter(
    "POWER_INPUT", 0.0, "power of forcing", &fdyn_turbu);

  setStringToIntegralParameter<int>(
    "SCALAR_FORCING",
    "no",
    "Define forcing for scalar field.",
    tuple<std::string>(
      "no",
      "isotropic",
      "mean_scalar_gradient"),
    tuple<std::string>(
      "Do not force the scalar field",
      "Force scalar field isotropically such as the fluid field.",
      "Force scalar field by imposed mean-scalar gradient."),
    tuple<int>(0,1,2),
    &fdyn_turbu);

  DoubleParameter(
    "MEAN_SCALAR_GRADIENT",
    0.0,
    "Value of imposed mean-scalar gradient to force scalar field.",
    &fdyn_turbu);

  // filtering with xfem
  //--------------

  BoolParameter("EXCLUDE_XFEM","No","Flag to (de)activate XFEM dofs in calculation of fine-scale velocity.",&fdyn_turbu);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for Smagorinsky model
  Teuchos::ParameterList& fdyn_turbsgv = fdyn.sublist("SUBGRID VISCOSITY",false,"");

  DoubleParameter("C_SMAGORINSKY",0.0,"Constant for the Smagorinsky model. Something between 0.1 to 0.24. Vreman constant if the constant vreman model is applied (something between 0.07 and 0.01).",&fdyn_turbsgv);
  DoubleParameter("C_YOSHIZAWA",-1.0,"Constant for the compressible Smagorinsky model: isotropic part of subgrid-stress tensor. About 0.09 or 0.0066. Ci will not be squared!",&fdyn_turbsgv);
  BoolParameter("C_SMAGORINSKY_AVERAGED","No","Flag to (de)activate averaged Smagorinksy constant",&fdyn_turbsgv);
  BoolParameter("C_INCLUDE_CI","No","Flag to (de)inclusion of Yoshizawa model",&fdyn_turbsgv);
  //remark: following Moin et al. 1991, the extension of the dynamic Smagorinsky model to compressibel flow
  //        also contains a model for the isotropic part of the subgrid-stress tensor according to Yoshizawa 1989
  //        although used in literature for turbulent variable-density flow at low Mach number, this additional term
  //        seems to destabilize the simulation when the flow is only weakly compressible
  //        therefore C_INCLUDE_CI allows to exclude this term
  //        if C_SMAGORINSKY_AVERAGED == true
  //           if C_INCLUDE_CI==true and C_YOSHIZAWA>=0.0 then the given value C_YOSHIZAWA is used
  //           if C_INCLUDE_CI==true and C_YOSHIZAWA<0.0 then C_YOSHIZAWA is determined dynamically
  //        else all values are taken from input

  DoubleParameter(
    "CHANNEL_L_TAU",
    0.0,
    "Used for normalisation of the wall normal distance in the Van \nDriest Damping function. May be taken from the output of \nthe apply_mesh_stretching.pl preprocessing script.",
    &fdyn_turbsgv);

  DoubleParameter("C_TURBPRANDTL",1.0,"(Constant) turbulent Prandtl number for the Smagorinsky model in scalar transport.",&fdyn_turbsgv);

  setStringToIntegralParameter<int>("FILTER_WIDTH","CubeRootVol","The Vreman model requires a filter width.",
                               tuple<std::string>(
                                 "CubeRootVol",
                                 "Direction_dependent",
                                 "Minimum_length"),
                                 tuple<int>(
                                     INPAR::FLUID::cuberootvol,
                                     INPAR::FLUID::dir_dep,
                                     INPAR::FLUID::min_len),
                               &fdyn_turbsgv);


  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for Smagorinsky model
  Teuchos::ParameterList& fdyn_wallmodel = fdyn.sublist("WALL MODEL",false,"");

  BoolParameter("X_WALL","No","Flag to switch on the xwall model",&fdyn_wallmodel);

  setStringToIntegralParameter<int>(
    "Tauw_Type",
    "constant",
    "Methods for calculating/updating the wall shear stress necessary for Spalding's law.",
    tuple<std::string>(
      "constant",
      "mean_between_steps",
      "mean_iter",
      "between_steps",
      "fix_point_iter_with_step_control",
      "fully_linearized"),
    tuple<std::string>(
      "Use the constant wall shear stress given in the input file for the whole simulation.",
      "Calculate wall shear stress in between time steps and use the mean value.",
      "Calculate wall shear stress at every iteration and use the mean value.",
      "Calculate wall shear stress in between time steps.",
      "Calculate wall shear stress in between every non-linear iteration and use fix point iteration to converge.\n Since this results in an unstable behavior, the increment of tauw is reduced if necessary.",
      "Fully linearized wall shear stress."),
    tuple<int>(0,1,2,3,4,5),
    &fdyn_wallmodel);

  setStringToIntegralParameter<int>(
    "Tauw_Calc_Type",
    "residual",
    "Methods for calculating the wall shear stress necessary for Spalding's law.",
    tuple<std::string>(
      "residual",
      "spalding",
      "gradient",
      "gradient_to_residual"),
    tuple<std::string>(
      "Residual (force) devided by area.",
      "Shear stress by inverse calculation via Spalding's law.",
      "Gradient via shape functions and nodal values.",
      "First gradient, then residual."),
    tuple<int>(0,1,2,3),
    &fdyn_wallmodel);

  IntParameter(
    "Switch_Step",
    -1,
    "Switch from gradient to residual based tauw.",
  &fdyn_wallmodel);

  setStringToIntegralParameter<int>(
    "Projection",
    "No",
    "Flag to switch projection of the enriched dofs after updating tauw, alternatively with or without continuity constraint.",
    tuple<std::string>(
      "No",
      "onlyl2projection",
      "l2projectionwithcontinuityconstraint"),
    tuple<std::string>(
      "Switch off projection.",
      "Only l2 projection.",
      "L2 projection with continuity constraint."),
    tuple<int>(0,1,2),
    &fdyn_wallmodel);

  DoubleParameter("C_Tauw",1.0,"Constant wall shear stress for Spalding's law, if applicable",&fdyn_wallmodel);

  DoubleParameter("Min_Tauw",2.0e-9,"Minimum wall shear stress preventing system to become singular",&fdyn_wallmodel);

  DoubleParameter("Inc_Tauw",1.0,"Increment of Tauw of full step, between 0.0 and 1.0",&fdyn_wallmodel);

  DoubleParameter("Penalty_Param",1000.0,"Penalty parameter for divergence free projection",&fdyn_wallmodel);

  setStringToIntegralParameter<int>(
    "Blending_Type",
    "none",
    "Methods for blending the enrichment space.",
    tuple<std::string>(
      "none",
      "ramp_function",
      "tauw_transformation"),
    tuple<std::string>(
      "No ramp function, does not converge!",
      "Enrichment is multiplied with linear ramp function resulting in zero enrichment at the interface",
      "Wall shear stress is modified at blending nodes such that y+ is constant at the interface"),
    tuple<int>(0,1,2),
    &fdyn_wallmodel);

  IntParameter("GP_Wall_Normal",3,"Gauss points in wall normal direction",&fdyn_wallmodel);

  IntParameter("GP_Wall_Parallel",3,"Gauss points in wall parallel direction",&fdyn_wallmodel);

  BoolParameter("Enr_MFS_Fine_Scale","No","Flag to add the enrichment part to the fine scale velocity. if no: enrichment part of coarse scales",&fdyn_wallmodel);

  IntParameter(
    "PROJECTION_SOLVER",
    -1,
    "Set solver number for l2-projection.",
  &fdyn_wallmodel);

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

  IntParameter(
    "ML_SOLVER",
    -1,
    "Set solver number for scale separation via level set transfer operators from plain aggregation.",
  &fdyn_turbmfs);

  BoolParameter("CALC_N","No","Flag to (de)activate calculation of N from the Reynolds number.",&fdyn_turbmfs);

  DoubleParameter(
    "N",
    1.0,
    "Set grid to viscous scale ratio.",
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

  BoolParameter("ADAPT_CSGS_PHI","No","Flag to (de)activate adaption of CsgsD to CsgsB.",&fdyn_turbmfs);

  BoolParameter("NEAR_WALL_LIMIT_CSGS_PHI","No","Flag to (de)activate near-wall limit for scalar field.",&fdyn_turbmfs);

  BoolParameter("CONSISTENT_FLUID_RESIDUAL","No","Flag to (de)activate the consistency term for residual-based stabilization.",&fdyn_turbmfs);

  DoubleParameter(
    "C_DIFF",
    1.0,
    "Proportionality constant between Re*Pr and ratio dissipative scale to element length. Usually equal cnu.",
    &fdyn_turbmfs);

  BoolParameter("SET_FINE_SCALE_VEL","No","Flag to set fine-scale velocity for parallel nightly tests.",&fdyn_turbmfs);

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
   Teuchos::ParameterList& twophasedyn = list->sublist("TWO PHASE FLOW",false,"");
   DoubleParameter("INTERFACE_THICKNESS",0.0,"Thickness of interface for multiphase flow",&twophasedyn);
   IntParameter("NUMSTEP",10,"Number of Time Steps",&twophasedyn);
   DoubleParameter("TIMESTEP",0.01,"Time increment dt",&twophasedyn);
   DoubleParameter("MAXTIME",0.0,"Total simulation time",&twophasedyn);
   DoubleParameter("CONVTOL",1E-6,"Tolerance for convergence check",&twophasedyn);
   IntParameter("UPRES",1,"Increment for writing solution",&twophasedyn);
   IntParameter("RESTARTEVRY",1,"Increment for writing restart",&twophasedyn);
   IntParameter("ITEMAX",1,"Maximum number of iterations in levelset-fluid loop",&twophasedyn);
   BoolParameter("WRITE_CENTER_OF_MASS","No","Write center of mass to file",&twophasedyn);
   BoolParameter("RESTART_SCATRA_INPUT","No","Use ScaTra field from .dat-file instead",&twophasedyn);
   BoolParameter("ENHANCED_GAUSSRULE","No","Set higher order gaussrule within the interface layer.",&twophasedyn);

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
    setStringToIntegralParameter<int>("SOLVESCATRA",
                               "no",
                               "Flag to (de)activate solving scalar transport in blood",
                               tuple<std::string>(
                                 "no",
                                 "yes"),
                               tuple<std::string>(
                                 "do not solve scatra",
                                 "solve scatra"),
                               tuple<int>(0,1),
                               &andyn);

    // number of linear solver used for arterial dynamics
    IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for arterial dynamics",&andyn);

    /*----------------------------------------------------------------------*/
     Teuchos::ParameterList& redtisdyn = list->sublist("COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC",false,"");
     DoubleParameter("CONVTOL_P",1E-6,"Coupled red_airway and tissue iteration convergence for pressure",&redtisdyn);
     DoubleParameter("CONVTOL_Q",1E-6,"Coupled red_airway and tissue iteration convergence for flux",&redtisdyn);
     IntParameter("MAXITER",5,"Maximum coupling iterations",&redtisdyn);
     setStringToIntegralParameter<int>("RELAXTYPE","norelaxation","Dynamic Relaxation Type",
                                   tuple<std::string>(
                                       "norelaxation",
                                       "fixedrelaxation",
                                       "Aitken",
                                       "SD"),
                                   tuple<int>(
                                       INPAR::ARTNET::norelaxation,
                                       INPAR::ARTNET::fixedrelaxation,
                                       INPAR::ARTNET::Aitken,
                                       INPAR::ARTNET::SD),
                                   &redtisdyn);
     DoubleParameter("TIMESTEP",0.01,"Time increment dt",&redtisdyn);
     IntParameter("NUMSTEP",1,"Number of Time Steps",&redtisdyn);
     DoubleParameter("MAXTIME",4.0,"",&redtisdyn);
     DoubleParameter("NORMAL",1.0,"",&redtisdyn);

    /*----------------------------------------------------------------------*/
     Teuchos::ParameterList& redawdyn = list->sublist("REDUCED DIMENSIONAL AIRWAYS DYNAMIC",false,"");

     setStringToIntegralParameter<int>("DYNAMICTYP","OneStepTheta",
                                  "OneStepTheta Scheme",
                                  tuple<std::string>(
                                    "OneStepTheta"
                                    ),
                                  tuple<int>(
                                   one_step_theta
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
     DoubleParameter("THETA",1.0,"One-step-theta time integration factor",&redawdyn);

     IntParameter("MAXITERATIONS",1,"maximum iteration steps",&redawdyn);
     DoubleParameter("TOLERANCE",1.0E-6,"tolerance",&redawdyn);

     // number of linear solver used for reduced dimensional airways dynamic
     IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for reduced dim arterial dynamics",&redawdyn);

     // Solve scatra flag
     setStringToIntegralParameter<int>("SOLVESCATRA",
                                       "no",
                                       "Flag to (de)activate solving scalar transport in blood",
                                       tuple<std::string>(
                                         "no",
                                         "yes"),
                                       tuple<std::string>(
                                         "do not solve scatra",
                                         "solve scatra"),
                                       tuple<int>(0,1),
                                       &redawdyn);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC",false,"");

  DoubleParameter("TIMESTEP",0.1,"",&adyn);
  IntParameter("NUMSTEP",41,"",&adyn);
  DoubleParameter("MAXTIME",4.0,"",&adyn);
  setStringToIntegralParameter<int>("ALE_TYPE","solid","ale mesh movement algorithm",
                               tuple<std::string>("solid","laplace","springs"),
                               tuple<int>(INPAR::ALE::solid,
                                   INPAR::ALE::laplace,
                                   INPAR::ALE::springs),
                               &adyn);
  IntParameter("MAXITER",1,"Maximum number of newton iterations.",&adyn);
  DoubleParameter("TOLRES",1.0e-06,"Absolute tolerance for length scaled L2 residual norm ",&adyn);
  DoubleParameter("TOLDISP",1.0e-06,"Absolute tolerance for length scaled L2 increment norm ",&adyn);

  IntParameter("NUM_INITSTEP",0,"",&adyn);
  IntParameter("RESTARTEVRY",1,"write restart data every RESTARTEVRY steps",&adyn);
  IntParameter("RESULTSEVRY",0,"write results every RESULTSTEVRY steps",&adyn);
  setStringToIntegralParameter<int>("DIVERCONT", "continue",
                                    "What to do if nonlinear solver does not converge?",
                                    tuple<std::string>(
                                        "stop",
                                        "continue"),
                                    tuple<int>(
                                        INPAR::ALE::divcont_stop,
                                        INPAR::ALE::divcont_continue),
                                    &adyn);

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
                                 "Gen_Alpha"
                                 ),
                               tuple<int>(
                                   INPAR::SCATRA::timeint_stationary,
                                   INPAR::SCATRA::timeint_one_step_theta,
                                   INPAR::SCATRA::timeint_bdf2,
                                   INPAR::SCATRA::timeint_gen_alpha
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

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string,12> name;
    Teuchos::Tuple<int,12> label;
    name[ 0] = "zero_field";                   label[ 0] = INPAR::SCATRA::initfield_zero_field;
    name[ 1] = "field_by_function";            label[ 1] = INPAR::SCATRA::initfield_field_by_function;
    name[ 2] = "field_by_condition";           label[ 2] = INPAR::SCATRA::initfield_field_by_condition;
    name[ 3] = "disturbed_field_by_function";  label[ 3] = INPAR::SCATRA::initfield_disturbed_field_by_function;
    name[ 4] = "1D_DISCONTPV";                 label[ 4] = INPAR::SCATRA::initfield_discontprogvar_1D;
    name[ 5] = "FLAME_VORTEX_INTERACTION";     label[ 5] = INPAR::SCATRA::initfield_flame_vortex_interaction;
    name[ 6] = "RAYTAYMIXFRAC";                label[ 6] = INPAR::SCATRA::initfield_raytaymixfrac;
    name[ 7] = "L_shaped_domain";              label[ 7] = INPAR::SCATRA::initfield_Lshapeddomain;
    name[ 8] = "facing_flame_fronts";          label[ 8] = INPAR::SCATRA::initfield_facing_flame_fronts;
    name[ 9] = "oracles_flame";                label[ 9] = INPAR::SCATRA::initfield_oracles_flame;
    name[10] = "high_forced_hit";              label[10] = INPAR::SCATRA::initialfield_forced_hit_high_Sc;
    name[11] = "low_forced_hit";               label[11] = INPAR::SCATRA::initialfield_forced_hit_low_Sc;

    setStringToIntegralParameter<int>(
        "INITIALFIELD",
        "zero_field",
        "Initial Field for scalar transport problem",
        name,
        label,
        &scatradyn);
  }

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

  // Parameters for reaction-diffusion systems (for example cardiac electrophysiology)
  IntParameter("WRITEMAXINTSTATE",0,"number of maximal internal state variables to be postprocessed",&scatradyn);
  IntParameter("WRITEMAXIONICCURRENTS",0,"number of maximal ionic currents to be postprocessed",&scatradyn);
  DoubleParameter("ACTTHRES",1.0,"threshold for the potential for computing and postprocessing activation time ",&scatradyn);

  setNumericStringParameter("WRITEFLUX_IDS","-1",
      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
      &scatradyn);

  BoolParameter("OUTMEAN","No","Output of mean values for scalars and density",&scatradyn);
  BoolParameter("OUTINTEGRREAC","No","Output of integral reaction values",&scatradyn);
  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&scatradyn);

  BoolParameter("MATLAB_STATE_OUTPUT","No","Do you want to write the state solution to Matlab file?",&scatradyn);

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

  setStringToIntegralParameter<int>("INITPOTCALC","no",
                                    "Automatically calculate initial field for electric potential",
                                    tuple<std::string>(
                                        "no",
                                        "yes"),
                                    tuple<int>(
                                        INPAR::SCATRA::initpotcalc_no,
                                        INPAR::SCATRA::initpotcalc_yes),
                                    &scatradyn);

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
      "Switch to block-preconditioned family of solvers, only works with block preconditioners like CheapSIMPLE!",&scatradyn);

  setStringToIntegralParameter<int>("SCATRATYPE","Undefined",
                               "Type of scalar transport problem",
                               tuple<std::string>(
                                 "Undefined",
                                 "ConvectionDiffusion",
                                 "LowMachNumberFlow",
                                 "Elch",
                                 "LevelSet",
                                 "Poroscatra",
                                 "Advanced_Reaction",
                                 "Poro_Scatra_Reaction",
                                 "AnisotropicDiffusion",
                                 "Cardiac_Monodomain"),
                               tuple<int>(
                                 INPAR::SCATRA::scatratype_undefined,
                                 INPAR::SCATRA::scatratype_condif,
                                 INPAR::SCATRA::scatratype_loma,
                                 INPAR::SCATRA::scatratype_elch,
                                 INPAR::SCATRA::scatratype_levelset,
                                 INPAR::SCATRA::scatratype_poro,
                                 INPAR::SCATRA::scatratype_advreac,
                                 INPAR::SCATRA::scatratype_pororeac,
                                 INPAR::SCATRA::scatratype_anisotrop,
                                 INPAR::SCATRA::scatratype_cardiac_monodomain),
                                 &scatradyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
                                  tuple<std::string>(
                                      "no",
                                      "Condensed_Smat",
                                      "Condensed_Bmat",
                                      "Condensed_Bmat_merged"), //use the condensed_bmat_merged strategy
                                    tuple<int>(
                                        INPAR::FLUID::no_meshtying,
                                        INPAR::FLUID::condensed_smat,
                                        INPAR::FLUID::condensed_bmat,
                                        INPAR::FLUID::condensed_bmat_merged),   //use the condensed_bmat_merged strategy
                                    &scatradyn);

  BoolParameter("ONLYPOTENTIAL","no",
        "Coupling of general ion transport equation with Laplace equation",&scatradyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for scalar transport/elch...",&scatradyn);
  //IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with SIMPLER)...",&scatradyn);

  // parameters for natural convection effects
  BoolParameter("NATURAL_CONVECTION","No","Include natural convection effects",&scatradyn);
  IntParameter("NATCONVITEMAX",10,"Maximum number of outer iterations for natural convection",&scatradyn);
  DoubleParameter("NATCONVCONVTOL",1e-6,"Convergence check tolerance for outer loop for natural convection",&scatradyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none, local, or global",
                                    tuple<std::string>(
                                      "none",
                                      "local",    // perform finite difference check on element level
                                      "global"),  // perform finite difference check on time integrator level
                                    tuple<int>(
                                        INPAR::SCATRA::fdcheck_none,
                                        INPAR::SCATRA::fdcheck_local,
                                        INPAR::SCATRA::fdcheck_global),
                                    &scatradyn);
  DoubleParameter("FDCHECKEPS",1.e-6,"dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, whereas smaller values don't)",&scatradyn);
  DoubleParameter("FDCHECKTOL",1.e-6,"relative tolerance for finite difference check",&scatradyn);

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

  // this parameter selects the characteristic element length for tau for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH",
                               "streamlength",
                               "Characteristic element length for tau",
                               tuple<std::string>(
                                 "streamlength",
                                 "volume_equivalent_diameter",
                                 "root_of_volume"),
                               tuple<int>(
                                   INPAR::SCATRA::streamlength,
                                   INPAR::SCATRA::volume_equivalent_diameter,
                                   INPAR::SCATRA::root_of_volume),
                               &scatradyn_stab);

  // this parameter selects the all-scale subgrid-diffusivity definition applied
  setStringToIntegralParameter<int>("DEFINITION_ASSGD",
                               "artificial_linear",
                               "Definition of (all-scale) subgrid diffusivity",
                               tuple<std::string>(
                                 "artificial_linear",
                                 "Hughes_etal_86_nonlinear",
                                 "Tezduyar_Park_86_nonlinear",
                                 "Tezduyar_Park_86_nonlinear_wo_phizero",
                                 "doCarmo_Galeao_91_nonlinear",
                                 "Almeida_Silva_97_nonlinear"),
                               tuple<std::string>(
                                 "classical linear artificial subgrid-diffusivity",
                                 "nonlinear isotropic according to Hughes et al. (1986)",
                                 "nonlinear isotropic according to Tezduyar and Park (1986)",
                                 "nonlinear isotropic according to Tezduyar and Park (1986) without user parameter phi_zero",
                                 "nonlinear isotropic according to doCarmo and Galeao (1991)",
                                 "nonlinear isotropic according to Almeida and Silva (1997)")  ,
                                tuple<int>(
                                    INPAR::SCATRA::assgd_artificial,
                                    INPAR::SCATRA::assgd_hughes,
                                    INPAR::SCATRA::assgd_tezduyar,
                                    INPAR::SCATRA::assgd_tezduyar_wo_phizero,
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
  IntParameter("ITEMAX_BEFORE_SAMPLING",1,"Maximum number of outer iterations before samling (for turbulent flows only)",&lomacontrol);
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

  setStringToIntegralParameter<int>("ELCHTYPE","Undefined",
                                 "Type of elch formulation",
                                 tuple<std::string>(
                                   "Undefined",
                                   "Nernst_Planck",
                                   "DiffCond"),
                                 tuple<int>(
                                   INPAR::ELCH::elchtype_undefined,
                                   INPAR::ELCH::elchtype_nernst_planck,
                                   INPAR::ELCH::elchtype_diffcond),
                                   &elchcontrol);
  IntParameter("MOVBOUNDARYITEMAX",10,"Maximum number of outer iterations in electrode shape change computations",&elchcontrol);
  DoubleParameter("MOVBOUNDARYCONVTOL",1e-6,"Convergence check tolerance for outer loop in electrode shape change computations",&elchcontrol);
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
  DoubleParameter("MOVBOUNDARYTHETA",0.0,"One-step-theta factor in electrode shape change computations",&elchcontrol);
  BoolParameter("GALVANOSTATIC","No","flag for galvanostatic mode",&elchcontrol);
  setStringToIntegralParameter<int>("GSTAT_APPROX_ELECT_RESIST",
                                 "relation_pot_cur",
                                 "relation of potential and current flow",
                                  tuple<std::string>(
                                    "relation_pot_cur",
                                    "effective_length_with_initial_cond",
                                    "effective_length_with_integrated_cond"),
                                   tuple<int>(
                                     INPAR::ELCH::approxelctresist_relpotcur,
                                     INPAR::ELCH::approxelctresist_effleninitcond,
                                     INPAR::ELCH::approxelctresist_efflenintegcond),
                                    &elchcontrol);
  IntParameter("GSTATCONDID_CATHODE",0,"condition id of electrode kinetics for cathode",&elchcontrol);
  IntParameter("GSTATCONDID_ANODE",1,"condition id of electrode kinetics for anode",&elchcontrol);
  DoubleParameter("GSTATCONVTOL",1.e-5,"Convergence check tolerance for galvanostatic mode",&elchcontrol);
  DoubleParameter("GSTATCURTOL",1.e-15,"Current Tolerance",&elchcontrol);
  IntParameter("GSTATCURVENO",-1,"function number defining the imposed current curve",&elchcontrol);
  IntParameter("GSTATITEMAX",10,"maximum number of iterations for galvanostatic mode",&elchcontrol);
  DoubleParameter("GSTAT_LENGTH_CURRENTPATH",0.0,"average length of the current path",&elchcontrol);
//  IntParameter("MAGNETICFIELD_FUNCNO",-1,"function number defining an externally imposed magnetic field",&elchcontrol);
  setStringToIntegralParameter<int>("EQUPOT",
                               "Undefined",
                               "type of closing equation for electric potential",
                                tuple<std::string>(
                                  "Undefined",
                                  "ENC",
                                  "ENC_PDE",
                                  "ENC_PDE_ELIM",
                                  "Poisson",
                                  "Laplace",
                                  "divi"),
                                 tuple<int>(
                                   INPAR::ELCH::equpot_undefined,
                                   INPAR::ELCH::equpot_enc,
                                   INPAR::ELCH::equpot_enc_pde,
                                   INPAR::ELCH::equpot_enc_pde_elim,
                                   INPAR::ELCH::equpot_poisson,
                                   INPAR::ELCH::equpot_laplace,
                                   INPAR::ELCH::equpot_divi),
                                  &elchcontrol);

  /*----------------------------------------------------------------------*/
  // attention: this list is a sublist of elchcontrol

    Teuchos::ParameterList& elchdiffcondcontrol = elchcontrol.sublist(
        "DIFFCOND",
        false,
        "control parameters for electrochemical diffusion conduction problems\n");

    BoolParameter("DIFFCOND_FORMULATION","No",
          "Activation of diffusion-conduction formulation",&elchdiffcondcontrol);

    BoolParameter("CURRENT_SOLUTION_VAR","No","Current as a solution variable",&elchdiffcondcontrol);

    BoolParameter("MAT_DIFFCOND_DIFFBASED","Yes",
        "Coupling terms of chemical diffusion for current equation are based on t and kappa",&elchdiffcondcontrol);

    /// dilute solution theory (diffusion potential in current equation):
    ///    A          B
    ///   |--|  |----------|
    ///   z_1 + (z_2 - z_1) t_1
    /// ------------------------ (RT/F kappa (1+f+-) 1/c_k grad c_k)
    ///      z_1 z_2
    ///     |________|
    ///         C
    //
    // default: concentrated solution theory according to Newman
    DoubleParameter("MAT_NEWMAN_CONST_A",2.0,"Constant A for the Newman model(term for the concentration overpotential)",&elchdiffcondcontrol);
    DoubleParameter("MAT_NEWMAN_CONST_B",-2.0,"Constant B for the Newman model(term for the concentration overpotential)",&elchdiffcondcontrol);
    DoubleParameter("MAT_NEWMAN_CONST_C",-1.0,"Constant C for the Newman model(term for the concentration overpotential)",&elchdiffcondcontrol);

    /*----------------------------------------------------------------------*/
    Teuchos::ParameterList& levelsetcontrol = list->sublist(
        "LEVEL-SET CONTROL",
        false,
        "control parameters for level-set problems\n");

    IntParameter("NUMSTEP",24,"Total number of time steps",&levelsetcontrol);
    DoubleParameter("TIMESTEP",0.1,"Time increment dt",&levelsetcontrol);
    DoubleParameter("MAXTIME",1000.0,"Total simulation time",&levelsetcontrol);
    IntParameter("UPRES",1,"Increment for writing solution",&levelsetcontrol);
    IntParameter("RESTARTEVRY",1,"Increment for writing restart",&levelsetcontrol);

    setStringToIntegralParameter<int>("CALCERROR","No",
                                 "compute error compared to analytical solution",
                                 tuple<std::string>(
                                   "No",
                                   "ZalesaksDisk"
                                   ),
                                 tuple<int>(
                                     INPAR::SCATRA::calcerror_no_ls,
                                     INPAR::SCATRA::calcerror_initial_field
                                     ),
                                     &levelsetcontrol);

    BoolParameter("EXTRACT_INTERFACE_VEL","No","replace computed velocity at nodes of given distance of interface by approximated interface velocity",&levelsetcontrol);
    IntParameter("NUM_CONVEL_LAYERS",-1,"number of layers around the interface which keep their computed convective velocity",&levelsetcontrol);


    /*----------------------------------------------------------------------*/
    Teuchos::ParameterList& ls_reinit = levelsetcontrol.sublist("REINITIALIZATION",false,"");

    setStringToIntegralParameter<int>("REINITIALIZATION","None",
                                 "Type of reinitialization strategy for level set function",
                                 tuple<std::string>(
                                   "None",
                                   "Signed_Distance_Function",
                                   "Sussman"),
                                 tuple<int>(
                                   INPAR::SCATRA::reinitaction_none,
                                   INPAR::SCATRA::reinitaction_signeddistancefunction,
                                   INPAR::SCATRA::reinitaction_sussman),
                                 &ls_reinit);

    BoolParameter("REINIT_INITIAL","No","Has level set field to be reinitialized before first time step?",&ls_reinit);
    IntParameter("REINITINTERVAL",1,"reinitialization interval",&ls_reinit);

    // parameters for signed distance reinitialization
    BoolParameter("REINITBAND","No","reinitialization only within a band around the interface, or entire domain?",&ls_reinit);
    DoubleParameter("REINITBANDWIDTH",1.0,"level-set value defining band width for reinitialization",&ls_reinit);

    // parameters for reinitialization equation
    IntParameter("NUMSTEPSREINIT",1,"(maximal) number of pseudo-time steps",&ls_reinit);
    // this parameter selects the tau definition applied
    setStringToIntegralParameter<int>("LINEARIZATIONREINIT",
                                 "fixed_point",
                                 "linearization of reinitialization equation",
                                 tuple<std::string>(
                                   "newton",
                                   "fixed_point"),
                                  tuple<int>(
                                      INPAR::SCATRA::newton,
                                      INPAR::SCATRA::fixed_point),
                                 &ls_reinit);
    DoubleParameter("TIMESTEPREINIT",1.0,"pseudo-time step length (usually a * characteristic element length of discretization with a>0)",&ls_reinit);
    DoubleParameter("THETAREINIT",1.0,"theta for time discretization of reinitialization equation",&ls_reinit);
    setStringToIntegralParameter<int>("STABTYPEREINIT",
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
                                 &ls_reinit);
    // this parameter selects the tau definition applied
    setStringToIntegralParameter<int>("DEFINITION_TAU_REINIT",
                                 "Taylor_Hughes_Zarins",
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
                                      INPAR::SCATRA::tau_exact_1d,
                                      INPAR::SCATRA::tau_zero),
                                 &ls_reinit);
    // this parameter governs whether all-scale subgrid diffusivity is included
    setStringToIntegralParameter<int>("ARTDIFFREINIT",
                                 "no",
                                 "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) term",
                                 tuple<std::string>(
                                   "no",
                                   "isotropic",
                                   "crosswind"),
                                 tuple<std::string>(
                                   "no artificial diffusion",
                                   "homogeneous artificial diffusion",
                                   "artificial diffusion in crosswind directions only")  ,
                                  tuple<int>(
                                      INPAR::SCATRA::artdiff_none,
                                      INPAR::SCATRA::artdiff_isotropic,
                                      INPAR::SCATRA::artdiff_crosswind),
                                 &ls_reinit);
    // this parameter selects the all-scale subgrid-diffusivity definition applied
    setStringToIntegralParameter<int>("DEFINITION_ARTDIFFREINIT",
                                 "artificial_linear",
                                 "Definition of (all-scale) subgrid diffusivity",
                                 tuple<std::string>(
                                   "artificial_linear",
                                   "artificial_linear_reinit",
                                   "Hughes_etal_86_nonlinear",
                                   "Tezduyar_Park_86_nonlinear",
                                   "Tezduyar_Park_86_nonlinear_wo_phizero",
                                   "doCarmo_Galeao_91_nonlinear",
                                   "Almeida_Silva_97_nonlinear",
                                   "YZbeta_nonlinear",
                                   "Codina_nonlinear"),
                                 tuple<std::string>(
                                   "simple linear artificial diffusion",
                                   "simple linear artificial diffusion const*h",
                                   "nonlinear isotropic according to Hughes et al. (1986)",
                                   "nonlinear isotropic according to Tezduyar and Park (1986)",
                                   "nonlinear isotropic according to Tezduyar and Park (1986) without user parameter phi_zero",
                                   "nonlinear isotropic according to doCarmo and Galeao (1991)",
                                   "nonlinear isotropic according to Almeida and Silva (1997)",
                                   "nonlinear YZ beta model",
                                   "nonlinear isotropic according to Codina")  ,
                                  tuple<int>(
                                      INPAR::SCATRA::assgd_artificial,
                                      INPAR::SCATRA::assgd_lin_reinit,
                                      INPAR::SCATRA::assgd_hughes,
                                      INPAR::SCATRA::assgd_tezduyar,
                                      INPAR::SCATRA::assgd_tezduyar_wo_phizero,
                                      INPAR::SCATRA::assgd_docarmo,
                                      INPAR::SCATRA::assgd_almeida,
                                      INPAR::SCATRA::assgd_yzbeta,
                                      INPAR::SCATRA::assgd_codina),
                                 &ls_reinit);

    setStringToIntegralParameter<int>("SMOOTHED_SIGN_TYPE",
                                      "SussmanSmerekaOsher1994",
                                      "sign function for reinitialization equation",
                                      tuple<std::string>(
                                      "NonSmoothed",
                                      "SussmanFatemi1999", // smeared-out Heaviside function
                                      "SussmanSmerekaOsher1994",
                                      "PengEtAl1999"),
                                      tuple<int>(
                                      INPAR::SCATRA::signtype_nonsmoothed,
                                      INPAR::SCATRA::signtype_SussmanFatemi1999,
                                      INPAR::SCATRA::signtype_SussmanSmerekaOsher1994,
                                      INPAR::SCATRA::signtype_PengEtAl1999),
                                      &ls_reinit);
    setStringToIntegralParameter<int>("CHARELELENGTHREINIT",
                                          "root_of_volume",
                                          "characteristic element length for sign function",
                                          tuple<std::string>(
                                          "root_of_volume",
                                          "streamlength"),
                                          tuple<int>(
                                          INPAR::SCATRA::root_of_volume_reinit,
                                          INPAR::SCATRA::streamlength_reinit),
                                          &ls_reinit);
    DoubleParameter("INTERFACE_THICKNESS", 1.0, "factor for interface thickness (multiplied by element length)", &ls_reinit);
    setStringToIntegralParameter<int>("VELREINIT",
                                          "integration_point_based",
                                          "evaluate velocity at integration point or compute node-based velocity",
                                          tuple<std::string>(
                                          "integration_point_based",
                                          "node_based"),
                                          tuple<int>(
                                          INPAR::SCATRA::vel_reinit_integration_point_based,
                                          INPAR::SCATRA::vel_reinit_node_based),
                                          &ls_reinit);
    setStringToIntegralParameter<int>("LINEARIZATIONREINIT",
                                              "newton",
                                              "linearization scheme for nonlinear convective term of reinitalization equation",
                                              tuple<std::string>(
                                              "newton",
                                              "fixed_point"),
                                              tuple<int>(
                                              INPAR::SCATRA::newton,
                                              INPAR::SCATRA::fixed_point),
                                              &ls_reinit);
    BoolParameter("CORRECTOR_STEP", "yes", "correction of interface position via volume constraint according to Sussman & Fatemi", &ls_reinit);
    DoubleParameter("CONVTOL_REINIT", -1.0, "tolerance for convergence check according to Sussman et al. 1994 (turned off negative)", &ls_reinit);

    BoolParameter("REINITVOLCORRECTION","No","volume correction after reinitialization",&ls_reinit);

    /*----------------------------------------------------------------------*/
    Teuchos::ParameterList& ls_particle = levelsetcontrol.sublist("PARTICLE",false,"");

    BoolParameter("INCLUDE_PARTICLE","No","Activate a hybrid particle-level-set method",&ls_particle);

    setStringToIntegralParameter<int>("DIMENSION","3D",
                                 "number of space dimensions for handling of quasi-2D problems with 3D particles",
                                 tuple<std::string>(
                                   "3D",
                                   "2Dx",
                                   "2Dy",
                                   "2Dz"),
                                 tuple<int>(
                                   INPAR::PARTICLE::particle_3D,
                                   INPAR::PARTICLE::particle_2Dx,
                                   INPAR::PARTICLE::particle_2Dy,
                                   INPAR::PARTICLE::particle_2Dz),
                                 &ls_particle);

    IntParameter("NUMPARTICLE",64,"number of particles in bins around interface (usually 3D: 64, 2D: 32)",&ls_particle);
    DoubleParameter("PARTICLEBANDWIDTH",3.0,"multiple of maximal element length defining band around interface filled with particles, i.e., alpha*max(dx,dy,dz): here we give alpha, max(dx,dy,dz) is defined by the cut_off radius for the bins!",&ls_particle);
    DoubleParameter("MIN_RADIUS",0.1,"minimal radius of particles, usually a multiple of min(dx,dy,dz)",&ls_particle);
    DoubleParameter("MAX_RADIUS",0.5,"maximal radius of particles, usually a multiple of min(dx,dy,dz)",&ls_particle);
    IntParameter("RESEEDING",10000,"reseeding interval",&ls_particle);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& biofilmcontrol = list->sublist(
      "BIOFILM CONTROL",
      false,
      "control parameters for biofilm problems\n");

  BoolParameter("BIOFILMGROWTH","No","Scatra algorithm for biofilm growth",&biofilmcontrol);
  BoolParameter("AVGROWTH","No","The calculation of growth parameters is based on averaged values",&biofilmcontrol);
  DoubleParameter("FLUXCOEF",0.0,"Coefficient for growth due to scalar flux",&biofilmcontrol);
  DoubleParameter("NORMFORCEPOSCOEF",0.0,"Coefficient for erosion due to traction normal surface forces",&biofilmcontrol);
  DoubleParameter("NORMFORCENEGCOEF",0.0,"Coefficient for erosion due to compression normal surface forces",&biofilmcontrol);
  DoubleParameter("TANGONEFORCECOEF",0.0,"Coefficient for erosion due to the first tangential surface force",&biofilmcontrol);
  DoubleParameter("TANGTWOFORCECOEF",0.0,"Coefficient for erosion due to the second tangential surface force",&biofilmcontrol);
  DoubleParameter("BIOTIMESTEP",0.05,"Time step size for biofilm growth",&biofilmcontrol);
  IntParameter("BIONUMSTEP",0,"Maximum number of steps for biofilm growth",&biofilmcontrol);
  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&biofilmcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptcontrol = list->sublist("TOPOLOGY OPTIMIZATION CONTROL",false,
      "control parameters for topology optimization problems");

  DoubleParameter("MAXTIME",10.0,"Total simulation time",&topoptcontrol);
  IntParameter("NUMSTEP",100,"Total number of timesteps",&topoptcontrol);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&topoptcontrol);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&topoptcontrol);
  IntParameter("UPRES",1,"Increment for writing solution",&topoptcontrol);

  setStringToIntegralParameter<int>("DENS_TYPE","node_based","type of optimization = density = porosity field",
      tuple<std::string>(
          "node_based",
          "element_based"),
          tuple<int>(
              INPAR::TOPOPT::dens_node_based,
              INPAR::TOPOPT::dens_ele_based),
              &topoptcontrol);

  setStringToIntegralParameter<int>("GRADIENT_TYPE","adjoints","basic type of adjoint equations",
      tuple<std::string>(
          "adjoints",
          "FD1",
          "FD2"),
          tuple<int>(
              INPAR::TOPOPT::gradientByAdjoints,
              INPAR::TOPOPT::gradientByFD1,
              INPAR::TOPOPT::gradientByFD2),
              &topoptcontrol);

  setStringToIntegralParameter<int>("RESTART_ACTION","Finished_Optimization_Step","Startint field of Restart",
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

  setStringToIntegralParameter<int>("OBJECTIVE_DISSIPATION","No","Cdissipation part of the objective function",
      tuple<std::string>(
          "No",
          "Yes",
          "Physical"),
          tuple<int>(
              INPAR::TOPOPT::obj_diss_no,
              INPAR::TOPOPT::obj_diss_yes,
              INPAR::TOPOPT::obj_diss_physical),
              &topoptcontrol);
  DoubleParameter("DISSIPATION_FAC",-1.0,"factor for the dissipation part of the objective",&topoptcontrol);

  BoolParameter("OBJECTIVE_PRESSURE_DROP","No","pressure drop part of the objective function",&topoptcontrol);
  DoubleParameter("PRESSURE_DROP_FAC",-1.0,"factor for the mean pressure drop part of the objective",&topoptcontrol);

  BoolParameter("OUTPUT_EVERY_ITER","No","write output of every iteration",&topoptcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptoptimizer = topoptcontrol.sublist("TOPOLOGY OPTIMIZER",false,
      "control parameters for the optimizer of a topology optimization problem");

  DoubleParameter("THETA",0.5,"theta for temporal integration of objective function",&topoptoptimizer);
  IntParameter("MAX_ITER",100,"Maximal number of optimization steps",&topoptoptimizer);
  IntParameter("MAX_GRAD_ITER",100,"Maximal number of optimization steps containing the gradient",&topoptoptimizer);
  IntParameter("MAX_INNER_ITER",20,"Maximal number of inner optimization steps",&topoptoptimizer);
  IntParameter("MAX_SUB_ITER",200,"Maximal iteration number within subproblem",&topoptoptimizer);
  IntParameter("MAX_INNER_SUB_ITER",50,"Maximal iteration number within inner subproblem routine",&topoptoptimizer);
  IntParameter("MATID",-1,"Material ID for automatic mesh generation",&topoptoptimizer);

  setStringToIntegralParameter<int>("INITIALFIELD","zero_field",
                               "Initial Field for density = topology optimization's optimization variable",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_by_function",
                                 "channel0",
                                 "channel05",
                                 "channel1",
                                 "channelstep0",
                                 "channelstep05",
                                 "channelstep1"),
                               tuple<int>(
                                   INPAR::TOPOPT::initdensfield_zero_field,
                                   INPAR::TOPOPT::initdensfield_field_by_function,
                                   INPAR::TOPOPT::initdensfield_channelflow0,
                                   INPAR::TOPOPT::initdensfield_channelflow05,
                                   INPAR::TOPOPT::initdensfield_channelflow1,
                                   INPAR::TOPOPT::initdensfield_channelstepflow0,
                                   INPAR::TOPOPT::initdensfield_channelstepflow05,
                                   INPAR::TOPOPT::initdensfield_channelstepflow1),
                               &topoptoptimizer);

  setStringToIntegralParameter<int>("TESTCASE","test_no","test case for optimizer",
      tuple<std::string>(
          "test_no",
          "test_snake_one_constr",
          "test_snake_multiple_constr",
          "test_workflow_without_fluiddata",
          "test_channel",
          "test_channel_with_step",
          "test_cornerflow",
          "test_lin_poro",
          "test_quad_poro",
          "test_cubic_poro"),
          tuple<int>(
              INPAR::TOPOPT::optitest_no,
              INPAR::TOPOPT::optitest_snake_one_constr,
              INPAR::TOPOPT::optitest_snake_multiple_constr,
              INPAR::TOPOPT::optitest_workflow_without_fluiddata,
              INPAR::TOPOPT::optitest_channel,
              INPAR::TOPOPT::optitest_channel_with_step,
              INPAR::TOPOPT::optitest_cornerflow,
              INPAR::TOPOPT::optitest_lin_poro,
              INPAR::TOPOPT::optitest_quad_poro,
              INPAR::TOPOPT::optitest_cub_poro),
              &topoptoptimizer);

  IntParameter("INITFUNCNO",-1,"function number for initial density field in topology optimization",&topoptoptimizer);
  DoubleParameter("VOLUME_BOUNDARY",0.7,"maximal percentage of fluid volume in background domain",&topoptoptimizer);
  DoubleParameter("TOL_KKT",1.0e-5,"tolerance of optimization problem (for KKT-conditions)",&topoptoptimizer);
  DoubleParameter("TOL_SUB",1.0e-9,"tolerance of subproblem",&topoptoptimizer);
  DoubleParameter("X_DIFF_MIN",1.0e-5,"minimal difference of upper and lower boundary of optimization variable",&topoptoptimizer);
  DoubleParameter("RHO_INIT",1.0e-2,"initial rho value",&topoptoptimizer);
  DoubleParameter("RHOMIN",1.0e-6,"minimal parameter value",&topoptoptimizer);
  DoubleParameter("FACMIN",1.0e-10,"minimal parameter value",&topoptoptimizer);
  IntParameter("UPRES",1,"Increment for writing solution",&topoptoptimizer);
  BoolParameter("GMSH_OUTPUT","No","Write Gmsh files",&topoptoptimizer);
  DoubleParameter("c_init",1000.0,"initial value for solver parameter",&topoptoptimizer);
  DoubleParameter("tol_sub_fac",1.001,"convergence factor for subproblem",&topoptoptimizer);
  DoubleParameter("tol_reducefac",0.1,"reduction factor for subproblem tolerance",&topoptoptimizer);
  DoubleParameter("resfac_sub",0.9,"residuum reduction factor for subproblem",&topoptoptimizer);
  DoubleParameter("fac_stepsize",-1.01,"factor for adjusting step size in every optimization step",&topoptoptimizer);
  DoubleParameter("RHO_FAC1",1.1,"factor for updating rho",&topoptoptimizer);
  DoubleParameter("RHO_FAC2",10.0,"factor for updating rho",&topoptoptimizer);
  DoubleParameter("asymptotes_fac1",10.0,"factor for updating asymptotes",&topoptoptimizer);
  DoubleParameter("asymptotes_fac2",0.01,"factor for updating asymptotes",&topoptoptimizer);
  DoubleParameter("fac_x_boundaries",0.1,"unsensible factor for computation of boundaries for optimization variable",&topoptoptimizer);
  DoubleParameter("fac_sub_reg",0.001,"regularisation factor in subproblem",&topoptoptimizer);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptadjointfluiddyn = topoptcontrol.sublist("TOPOLOGY ADJOINT FLUID",false,
      "control parameters for the adjoint fluid of a topology optimization problem");

  setStringToIntegralParameter<int>("ADJOINT_TYPE","discrete_adjoint","basic type of adjoint equations",
      tuple<std::string>(
          "discrete_adjoint",
          "cont_adjoint"),
          tuple<int>(
              INPAR::TOPOPT::discrete_adjoint,
              INPAR::TOPOPT::cont_adjoint),
              &topoptadjointfluiddyn);

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
          "test_primal",
          "test_stat_const_vel_lin_pres",
          "test_stat_lin_vel_quad_pres",
          "test_stat_quad_vel_lin_pres",
          "test_stat_all_terms_all_constants",
          "test_instat_varying_theta",
          "test_instat_all_terms_all_constants",
          "test_instat_primal_and_dual"),
          tuple<int>(
              INPAR::TOPOPT::adjointtest_no,
              INPAR::TOPOPT::adjointtest_primal,
              INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres,
              INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres,
              INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres,
              INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants,
              INPAR::TOPOPT::adjointtest_instat_varying_theta,
              INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants,
              INPAR::TOPOPT::adjointtest_instat_primal_and_dual),
              &topoptadjointfluiddyn);

  IntParameter("INITFUNCNO",-1,"Function for initial field",&topoptadjointfluiddyn);

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

  BoolParameter("GENERATE_FLOW_FIELD","No","Do not solve g-function, since we merely want to set up a fluid field",&combustcontrol);

  BoolParameter("RESTART_FROM_FLUID","No","Restart from a standard fluid problem (no scalar transport field). No XFEM dofs allowed!",&combustcontrol);
  BoolParameter("RESTART_SCATRA_INPUT","No","Use ScaTra field from .dat-file instead",&combustcontrol);

  BoolParameter("WRITE_CENTER_OF_MASS","No","write center of mass to file",&combustcontrol);

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

  BoolParameter("XFEMSTABILIZATION","Yes","Switch on/off face integrals based on edge-based stabilization, i.e., ghost penalty terms",&combustcontrolfluid);

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
          "ExtrapolationOld",
          "ExtrapolationNew",
          "MixedSemiLagrangeExtrapolation",
          "MixedSemiLagrangeExtrapolationNew",
          "MixedGhostvalSemiLagrange",
          "MixedGhostvalExtrapolation",
          "MixedGhostvalSemiLagrangeExtrapolation"),
          tuple<int>(
              INPAR::COMBUST::xfemtimeint_donothing,
              INPAR::COMBUST::xfemtimeint_semilagrange,
              INPAR::COMBUST::xfemtimeint_extrapolationold,
              INPAR::COMBUST::xfemtimeint_extrapolationnew,
              INPAR::COMBUST::xfemtimeint_mixedSLExtrapol,
              INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew,
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
          "disturbed_field_by_function",
          "flame_vortex_interaction",
          "darrieus_landau_instability",
          "beltrami_flow"),
          tuple<int>(
              INPAR::COMBUST::initfield_zero_field,
              INPAR::COMBUST::initfield_field_by_function,
              INPAR::COMBUST::initfield_disturbed_field_by_function,
              INPAR::COMBUST::initfield_flame_vortex_interaction,
              INPAR::COMBUST::initfield_darrieus_landau_instability,
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
          "smooth_grad_phi_leastsquares_2Dz",
          "smooth_grad_phi_l2_projection"),
          tuple<int>(
              INPAR::COMBUST::smooth_grad_phi_none,
              INPAR::COMBUST::smooth_grad_phi_meanvalue,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_3D,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy,
              INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz,
              INPAR::COMBUST::smooth_grad_phi_l2_projection),
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

  IntParameter("ITE_MAX_FRS",1,"The maximal number of iterations between fluid and recomputation of reference solution",&combustcontrolfluid);
  DoubleParameter("LAMINAR_FLAMESPEED",1.0,"The laminar flamespeed incorporates all chemical kinetics into the problem for now",&combustcontrolfluid);
  DoubleParameter("MOL_DIFFUSIVITY",0.0,"Molecular diffusivity",&combustcontrolfluid);
  DoubleParameter("MARKSTEIN_LENGTH",0.0,"The Markstein length takes flame curvature into account",&combustcontrolfluid);
  DoubleParameter("NITSCHE_VELOCITY",100.0,"Nitsche parameter to stabilize/penalize the velocity jump",&combustcontrolfluid);
  DoubleParameter("NITSCHE_PRESSURE",0.0,"Nitsche parameter to stabilize/penalize the pressure jump",&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_CONVFLUX","Yes","(De)activate Nitsche convective flux term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_CONVSTAB","Yes","(De)activate Nitsche convective stabilization term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_CONVPENALTY","No","(De)activate Nitsche convective penalty term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_MASS","No","(De)activate Nitsche mass conservation term",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("NITSCHE_WEIGHT","intersection_visc_based_harmonic","Definition of weighting",
                                    tuple<std::string>(
                                    "visc_based_harmonic",
                                    "intersection_visc_based_harmonic"),
                                    tuple<int>(
                                    INPAR::COMBUST::weight_visc_based_harmonic,
                                    INPAR::COMBUST::weight_intersection_visc_based_harmonic),
                                    &combustcontrolfluid);
  setStringToIntegralParameter<int>("CONNECTED_INTERFACE","No","laplace-beltrami surface-tension evaluation only: consider boundary integrals if interface is not closed",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("SMOOTHED_BOUNDARY_INTEGRATION","No","Turn on/off usage of smoothed normal vectors at interface",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);
  setStringToIntegralParameter<int>("INITSTATSOL","No","Compute stationary solution as initial solution",
                                     yesnotuple,yesnovalue,&combustcontrolfluid);

  BoolParameter("L2_PROJECTION_SECOND_DERIVATIVES","No","L2 Projection Second Derivatives of Level Set",&combustcontrolfluid);
  setStringToIntegralParameter<int>("NODAL_CURVATURE","l2_projected","Type of calculation of nodal curvature value",
      tuple<std::string>(
          "l2_projected",
          "averaged"),
          tuple<int>(
              INPAR::COMBUST::l2_projected,
              INPAR::COMBUST::averaged),
              &combustcontrolfluid);
  DoubleParameter("SMOOTHING_PARAMETER",0.0,"Diffusion Coefficient for Smoothing",&combustcontrolfluid);

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

  BoolParameter("REPELLANT_FORCE","No","Activate repellant force for turbulent bubbly channel flow",&combustcontrolfluid);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& combustcontrolgfunc = combustcontrol.sublist("COMBUSTION GFUNCTION",false,
      "control parameters for the G-function (level set) field of a combustion problem");

  setStringToIntegralParameter<int>("REFINEMENT","No","Turn refinement strategy for level set function on/off",
                                     yesnotuple,yesnovalue,&combustcontrolgfunc);
  IntParameter("REFINEMENTLEVEL",-1,"number of refinement level for refinement strategy",&combustcontrolgfunc);
  /*----------------------------------------------------------------------*/


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fsidyn = list->sublist(
    "FSI DYNAMIC",false,
    "Fluid Structure Interaction\n"
    "FSI solver with various coupling methods"
    );

  Teuchos::Tuple<std::string,23> name;
  Teuchos::Tuple<int,23> label;

  name[ 0] = "basic_sequ_stagg";                              label[ 0] = fsi_basic_sequ_stagg;
  name[ 1] = "iter_stagg_fixed_rel_param";                    label[ 1] = fsi_iter_stagg_fixed_rel_param;
  name[ 2] = "iter_stagg_AITKEN_rel_param";                   label[ 2] = fsi_iter_stagg_AITKEN_rel_param;
  name[ 3] = "iter_stagg_steep_desc";                         label[ 3] = fsi_iter_stagg_steep_desc;
  name[ 4] = "iter_stagg_NLCG";                               label[ 4] = fsi_iter_stagg_NLCG;
  name[ 5] = "iter_stagg_MFNK_FD";                            label[ 5] = fsi_iter_stagg_MFNK_FD;
  name[ 6] = "iter_stagg_MFNK_FSI";                           label[ 6] = fsi_iter_stagg_MFNK_FSI;
  name[ 7] = "iter_stagg_MPE";                                label[ 7] = fsi_iter_stagg_MPE;
  name[ 8] = "iter_stagg_RRE";                                label[ 8] = fsi_iter_stagg_RRE;
  name[ 9] = "iter_monolithicfluidsplit";                     label[ 9] = fsi_iter_monolithicfluidsplit;
  name[10] = "iter_monolithicstructuresplit";                 label[10] = fsi_iter_monolithicstructuresplit;
  name[11] = "iter_lung_monolithicstructuresplit";            label[11] = fsi_iter_lung_monolithicstructuresplit;
  name[12] = "iter_lung_monolithicfluidsplit";                label[12] = fsi_iter_lung_monolithicfluidsplit;
  name[13] = "iter_xfem_monolithic";                          label[13] = fsi_iter_xfem_monolithic;
  name[14] = "pseudo_structure";                              label[14] = fsi_pseudo_structureale;
  name[15] = "iter_constr_monolithicfluidsplit";              label[15] = fsi_iter_constr_monolithicfluidsplit;
  name[16] = "iter_constr_monolithicstructuresplit";          label[16] = fsi_iter_constr_monolithicstructuresplit;
  name[17] = "iter_mortar_monolithicstructuresplit";          label[17] = fsi_iter_mortar_monolithicstructuresplit;
  name[18] = "iter_mortar_monolithicfluidsplit";              label[18] = fsi_iter_mortar_monolithicfluidsplit;
  name[19] = "iter_fluidfluid_monolithicstructuresplit";      label[19] = fsi_iter_fluidfluid_monolithicstructuresplit;
  name[20] = "iter_fluidfluid_monolithicfluidsplit";          label[20] = fsi_iter_fluidfluid_monolithicfluidsplit;
  name[21] = "iter_fluidfluid_monolithicstructuresplit_nox";  label[21] = fsi_iter_fluidfluid_monolithicstructuresplit_nox;
  name[22] = "iter_fluidfluid_monolithicfluidsplit_nox";      label[22] = fsi_iter_fluidfluid_monolithicfluidsplit_nox;

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_AITKEN_rel_param",
                                    "Iteration Scheme over the fields",
                                    name,
                                    label,
                                    &fsidyn);

  setStringToIntegralParameter<int>("DEBUGOUTPUT","No",
                                    "Output of unconverged interface values during FSI iteration.\n"
                                    "There will be a new control file for each time step.\n"
                                    "This might be helpful to understand the coupling iteration.",
                                    tuple<std::string>(
                                      "No",
                                      "Yes",
                                      "no",
                                      "yes",
                                      "NO",
                                      "YES",
                                      "Interface",
                                      "Preconditioner",
                                      "All"
                                      ),
                                    tuple<int>(
                                      0,
                                      1,
                                      0,
                                      1,
                                      0,
                                      1,
                                      1,
                                      2,
                                      256
                                      ),
                                    &fsidyn);

  setStringToIntegralParameter<int>("SECONDORDER","No",
                               "Second order displacement-velocity conversion at the interface.",
                               yesnotuple,yesnovalue,&fsidyn);

  IntParameter("NUMSTEP",200,"Total number of Timesteps",&fsidyn);
  IntParameter("UPRES",1,"Increment for writing solution",&fsidyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&fsidyn);

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&fsidyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&fsidyn);

  setStringToIntegralParameter<int>("SLIDEALEPROJ","None",
                                 "Projection method to use for sliding FSI.",
                                 tuple<std::string>(
                                     "None",
                                     "Curr",
                                     "Ref",
                                     "RotZ",
                                     "RotZSphere"),
                                 tuple<int>(
                                     INPAR::FSI::ALEprojection_none,
                                     INPAR::FSI::ALEprojection_curr,
                                     INPAR::FSI::ALEprojection_ref,
                                     INPAR::FSI::ALEprojection_rot_z,
                                     INPAR::FSI::ALEprojection_rot_zsphere),
                                 &fsidyn);

  setStringToIntegralParameter<int> ("DIVPROJECTION", "no", "Project velocity into divergence-free subspace for partitioned fsi",
                                     yesnotuple,yesnovalue,&fsidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for time step size adaptivity in fsi dynamics */
  Teuchos::ParameterList& fsiadapt = fsidyn.sublist("TIMEADAPTIVITY",false,"");

  DoubleParameter("DTMAX", 0.1, "Limit maximally permitted time step size (>0)", &fsiadapt);
  DoubleParameter("DTMIN", 1.0e-4, "Limit minimally allowed time step size (>0)", &fsiadapt);

  DoubleParameter("LOCERRTOLFLUID", 1.0e-3, "Tolerance for the norm of local velocity error", &fsiadapt);

  IntParameter("ADAPTSTEPMAX", 5, "Maximum number of repetitions of one time step for adapting/reducing the time step size (>0)", &fsiadapt);
  DoubleParameter("SIZERATIOMAX", 2.0, "Limit maximally permitted change of time step size compared to previous size (>0).", &fsiadapt);
  DoubleParameter("SIZERATIOMIN", 0.5, "Limit minimally permitted change of time step size compared to previous size (>0).", &fsiadapt);
  DoubleParameter("SAFETYFACTOR", 0.9, "This is a safety factor to scale theoretical optimal step size, should be lower than 1 and must be larger than 0", &fsiadapt);

  IntParameter("NUMINCREASESTEPS", 0, "Number of consecutive steps that want to increase time step size before\n"
                                      "actually increasing it. Set 0 to deactivate this feature.", &fsiadapt);

  setNumericStringParameter("AVERAGINGDT", "0.3 0.7",
                            "Averaging of time step sizes in case of increasing time step size.\n"
                            "Parameters are ordered from most recent weight to the most historic one.\n"
                            "Number of parameters determines the number of previous time steps that are involved\n"
                            "in the averaging procedure.",
                            &fsiadapt);

  setStringToIntegralParameter<int>("TIMEADAPTON", "No",
                                    "Activate or deactivate time step size adaptivity",
                                    yesnotuple,yesnovalue, &fsiadapt);

  setStringToIntegralParameter<int>("AUXINTEGRATORFLUID", "AB2",
                                    "Method for error estimation in the fluid field",
                                    tuple<std::string>(
                                        "None",
                                        "ExplicitEuler",
                                        "AB2"),
                                    tuple<int>(
                                        INPAR::FSI::timada_fld_none,
                                        INPAR::FSI::timada_fld_expleuler,
                                        INPAR::FSI::timada_fld_adamsbashforth2),
                                    &fsiadapt);

  setStringToIntegralParameter<int>("DIVERCONT", "stop",
                                    "What to do if nonlinear solver does not converge?",
                                    tuple<std::string>(
                                        "stop",
                                        "continue",
                                        "halve_step",
                                        "revert_dt"),
                                    tuple<int>(
                                        INPAR::FSI::divcont_stop,
                                        INPAR::FSI::divcont_continue,
                                        INPAR::FSI::divcont_halve_step,
                                        INPAR::FSI::divcont_revert_dt),
                                    &fsiadapt);

  /*----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic FSI solvers */
  Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER",false,"");

  setStringToIntegralParameter<int>("SHAPEDERIVATIVES","No",
                                    "Include linearization with respect to mesh movement in Navier Stokes equation.",
                                    yesnotuple,yesnovalue,&fsimono);

  setStringToIntegralParameter<int>("ENERGYFILE","No",
                                    "Write artificial interface energie due to temporal discretization to file",
                                    yesnotuple,yesnovalue,&fsimono);

  IntParameter("ITEMAX",100,"Maximum allowed number of nonlinear iterations",&fsimono);

  setStringToIntegralParameter<int>("INFNORMSCALING","Yes","Scale Blocks with row infnorm?",
                                    yesnotuple,yesnovalue,&fsimono);

  setStringToIntegralParameter<int>("SYMMETRICPRECOND","No","Symmetric block GS preconditioner or ordinary GS",
                                    yesnotuple,yesnovalue,&fsimono);

  IntParameter("PRECONDREUSE", 0, "Number of preconditioner reused in monolithic FSI", &fsimono);

  setStringToIntegralParameter<int>(
                               "LINEARBLOCKSOLVER","PreconditionedKrylov",
                               "Linear solver algorithm for monolithic block system in monolithic FSI.\n"
                               "Most of the time preconditioned Krylov is the right thing to choose. But there are\n"
                               "block Gauss-Seidel methods as well.",
                               tuple<std::string>(
                                 "PreconditionedKrylov",
                                 "FSIAMG",
                                 "AMGnxn"
                                 ),
                               tuple<int>(
                                 INPAR::FSI::PreconditionedKrylov,
                                 INPAR::FSI::FSIAMG,
                                 INPAR::FSI::AMGnxn
                                 ),
                               &fsimono);

  setStringToIntegralParameter<int>("FSIAMGANALYZE","No",
                                   "run analysis on fsiamg multigrid scheme",
                                   yesnotuple,yesnovalue,&fsimono);

  // monolithic preconditioner parameter
  setNumericStringParameter("STRUCTPCOMEGA","1.0 1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on structural block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("STRUCTPCITER","1 1 1 1",
                            "Number of Richardson iterations on structural block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("FLUIDPCOMEGA","1.0 1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on fluid block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("FLUIDPCITER","1 1 1 1",
                            "Number of Richardson iterations on fluid block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("ALEPCOMEGA","1.0 1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on ale block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("ALEPCITER","1 1 1 1",
                            "Number of Richardson iterations on ale block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);

  setNumericStringParameter("PCOMEGA","1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on whole MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("PCITER","1 1 1",
                            "Number of Richardson iterations on whole MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);

  StringParameter("BLOCKSMOOTHER","BGS BGS BGS",
                  "Type of block smoother, can be BGS or Schur",
                  &fsimono);

  setNumericStringParameter("SCHUROMEGA","0.001 0.01 0.1",
                            "Damping factor for Schur complement construction",
                            &fsimono);

  DoubleParameter("ADAPTIVEDIST",0.0,
                  "Required distance for adaptive convergence check in Newton-type FSI.\n"
                  "This is the improvement we want to achieve in the linear extrapolation of the\n"
                  "adaptive convergence check. Set to zero to avoid the adaptive check altogether.",
                  &fsimono);

  DoubleParameter("BASETOL",1e-3,
                  "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
                  "This tolerance will be used for the linear solve of the FSI block system.\n"
                  "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
                  "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
                  "to the nonlinear convergence test using a absolute residual norm.",
                  &fsimono);

  DoubleParameter("CONVTOL",1e-6,"Nonlinear tolerance for lung/constraint/fluid-fluid FSI",&fsimono); // ToDo remove

  // Iteration parameters for convergence check of newton loop
  // for implementations without NOX
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
                               &fsimono);

  // for implementations without NOX
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
                                 &fsimono);

  // for implementations without NOX
  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                                 tuple<std::string>(
                                   "And"),
                                 tuple<int>(
                                 INPAR::FSI::bop_and),
                               &fsimono);

  // tolerances for convergence check of nonlinear solver in monolithic FSI
  // structure displacements
  DoubleParameter("TOL_DIS_RES_L2",1e-6,"Absolute tolerance for structure displacement residual in L2-norm",&fsimono);
  DoubleParameter("TOL_DIS_RES_INF",1e-6,"Absolute tolerance for structure displacement residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_DIS_INC_L2",1e-6,"Absolute tolerance for structure displacement increment in L2-norm",&fsimono);
  DoubleParameter("TOL_DIS_INC_INF",1e-6,"Absolute tolerance for structure displacement increment in Inf-norm",&fsimono);
  // interface tolerances
  DoubleParameter("TOL_FSI_RES_L2",1e-6,"Absolute tolerance for interface residual in L2-norm",&fsimono);
  DoubleParameter("TOL_FSI_RES_INF",1e-6,"Absolute tolerance for interface residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_FSI_INC_L2",1e-6,"Absolute tolerance for interface increment in L2-norm",&fsimono);
  DoubleParameter("TOL_FSI_INC_INF",1e-6,"Absolute tolerance for interface increment in Inf-norm",&fsimono);
  // fluid pressure
  DoubleParameter("TOL_PRE_RES_L2",1e-6,"Absolute tolerance for fluid pressure residual in L2-norm",&fsimono);
  DoubleParameter("TOL_PRE_RES_INF",1e-6,"Absolute tolerance for fluid pressure residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_PRE_INC_L2",1e-6,"Absolute tolerance for fluid pressure increment in L2-norm",&fsimono);
  DoubleParameter("TOL_PRE_INC_INF",1e-6,"Absolute tolerance for fluid pressure increment in Inf-norm",&fsimono);
  // fluid velocities
  DoubleParameter("TOL_VEL_RES_L2",1e-6,"Absolute tolerance for fluid velocity residual in L2-norm",&fsimono);
  DoubleParameter("TOL_VEL_RES_INF",1e-6,"Absolute tolerance for fluid velocity residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_VEL_INC_L2",1e-6,"Absolute tolerance for fluid velocity increment in L2-norm",&fsimono);
  DoubleParameter("TOL_VEL_INC_INF",1e-6,"Absolute tolerance for fluid velocity increment in Inf-norm",&fsimono);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic FSI solvers */
  Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER",false,"");

  setStringToIntegralParameter<int>(
                                 "PARTITIONED","DirichletNeumann",
                                 "Coupling strategies for partitioned FSI solvers.",
                                 tuple<std::string>(
                                   "DirichletNeumann",
                                   "DirichletNeumannSlideALE"
                                   ),
                                 tuple<int>(
                                   INPAR::FSI::DirichletNeumann,
                                   INPAR::FSI::DirichletNeumannSlideale
                                   ),
                                 &fsipart);

  setStringToIntegralParameter<int>("PREDICTOR","d(n)",
                                 "Predictor for interface displacements",
                                 tuple<std::string>(
                                   "d(n)",
                                   "d(n)+dt*(1.5*v(n)-0.5*v(n-1))",
                                   "d(n)+dt*v(n)",
                                   "d(n)+dt*v(n)+0.5*dt^2*a(n)"
                                   ),
                                 tuple<int>(1,2,3,4),
                                 &fsipart);

    setStringToIntegralParameter<int>("COUPVARIABLE","Displacement",
                                 "Coupling variable at the interface",
                                 tuple<std::string>("Displacement","Force"),
                                 tuple<int>(0,1),
                                 &fsipart);

    setStringToIntegralParameter<int>("COUPMETHOD","conforming",
                                 "Coupling Method Mortar (mtr) or conforming nodes at interface",
                                 tuple<std::string>(
                                   "MTR",
                                   "Mtr",
                                   "mtr",
                                   "conforming",
                                   "immersed"
                                   ),
                                 tuple<int>(0,0,0,1,2),
                                 &fsipart);

    DoubleParameter("BASETOL",1e-3,
                    "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
                    "This tolerance will be used for the linear solve of the FSI block system.\n"
                    "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
                    "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
                    "to the nonlinear convergence test using a absolute residual norm.",
                    &fsipart);

    DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields in case of partitioned scheme",&fsipart);
    DoubleParameter("RELAX",1.0,"fixed relaxation parameter for partitioned FSI solvers",&fsipart);
    DoubleParameter("MAXOMEGA",0.0,"largest omega allowed for Aitken relaxation (0.0 means no constraint)",&fsipart);
    IntParameter("ITEMAX",100,"Maximum number of iterations over fields",&fsipart);

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
    Teuchos::ParameterList& immersedmethod = list->sublist(
      "IMMERSED METHOD",false,
      "General parameters for any immersed problem"
      );

    setStringToIntegralParameter<int>(
                                 "COUPALGO","partitioned",
                                 "Coupling strategies for immersed method.",
                                 tuple<std::string>(
                                   "partitioned",
                                   "monolithic"),
                                   tuple<int>(
                                   INPAR::IMMERSED::partitioned,
                                   INPAR::IMMERSED::monolithic),
                                   &immersedmethod);

    setStringToIntegralParameter<int>(
                                 "SCHEME","dirichletneumann",
                                 "Coupling schemes for partitioned immersed method.",
                                 tuple<std::string>(
                                   "neumannneumann",
                                   "dirichletneumann"),
                                   tuple<int>(
                                   INPAR::IMMERSED::neumannneumann,
                                   INPAR::IMMERSED::dirichletneumann),
                                   &immersedmethod);


    setStringToIntegralParameter<int>(
                                 "PROJECTION","shapefunctions",
                                 "Projection of nodal values between the coupled fields.",
                                 tuple<std::string>(
                                   "shapefunctions",
                                   "mortar"),
                                   tuple<int>(
                                   INPAR::IMMERSED::shapefunctions,
                                   INPAR::IMMERSED::mortar),
                                   &immersedmethod);

    setStringToIntegralParameter<int>(
                                 "APPLY_FORCE_RELAX","globally",
                                 "Relax whole force vector or not.",
                                 tuple<std::string>(
                                   "globally",
                                   "selectively"),
                                   tuple<int>(
                                   INPAR::IMMERSED::globally,
                                   INPAR::IMMERSED::selectively),
                                   &immersedmethod);

    setStringToIntegralParameter<int>(
                                 "APPLY_VEL_RELAX","globally",
                                 "Relax whole velocity vector or not.",
                                 tuple<std::string>(
                                   "globally",
                                   "selectively"),
                                   tuple<int>(
                                   INPAR::IMMERSED::globally,
                                   INPAR::IMMERSED::selectively),
                                   &immersedmethod);

    setStringToIntegralParameter<int>(
                                 "DIVERCONT","stop",
                                 "What to do after maxiter is reached.",
                                 tuple<std::string>(
                                   "stop",
                                   "continue"),
                                   tuple<int>(
                                   INPAR::IMMERSED::nlnsolver_stop,
                                   INPAR::IMMERSED::nlnsolver_continue),
                                   &immersedmethod);

    setStringToIntegralParameter<int>(
                                 "DETECT_VEL_SIGNIFICANCE","no",
                                 "set projected velocity to veln if not of significant magnitude",
                                 tuple<std::string>(
                                   "yes",
                                   "no"),
                                   tuple<int>(
                                   1,
                                   0),
                                   &immersedmethod);

    DoubleParameter("FORCE_RELAX",1.0,"Force Relaxaton Parameter"    ,&immersedmethod);
    DoubleParameter("VEL_RELAX"  ,1.0,"Velocity Relaxation Parameter",&immersedmethod);
    DoubleParameter("FLD_SRCHRADIUS_FAC",1.0,"fac times fluid ele. diag. length",&immersedmethod);
    DoubleParameter("STRCT_SRCHRADIUS_FAC",0.5,"fac times structure bounding box diagonal",&immersedmethod);

  /*----------------------------------------------------------------------*/
    Teuchos::ParameterList& fpsidyn = list->sublist(
      "FPSI DYNAMIC",false,
      "Fluid Porous Structure Interaction\n"
      "FPSI solver with various coupling methods"
      );

    Teuchos::Tuple<std::string,1> fpsiname;
    Teuchos::Tuple<int,1> fpsilabel;

    name[0] = "fpsi_monolithic_plain"; label[0] = fpsi_monolithic_plain;
    setStringToIntegralParameter<int>("COUPALGO","fpsi_monolithic_plain",
                                      "Iteration Scheme over the fields",
                                      name,
                                      label,
                                      &fpsidyn);

    setStringToIntegralParameter<int>("SHAPEDERIVATIVES","No",
                                 "Include linearization with respect to mesh movement in Navier Stokes equation.\n"
                                 "Supported in monolithic FPSI for now.",
                                 yesnotuple,yesnovalue,&fpsidyn);

    setStringToIntegralParameter<int>("USESHAPEDERIVATIVES","No",
                                 "Add linearization with respect to mesh movement in Navier Stokes equation to stiffness matrix.\n"
                                 "Supported in monolithic FPSI for now.",
                                 yesnotuple,yesnovalue,&fpsidyn);

    setStringToIntegralParameter<int>(
                                 "PARTITIONED","RobinNeumann",
                                 "Coupling strategies for partitioned FPSI solvers.",
                                 tuple<std::string>(
                                   "RobinNeumann",
                                   "DirichletNeumann",
                                   "monolithic",
                                   "nocoupling"),
                                   tuple<int>(
                                   INPAR::FPSI::RobinNeumann,
                                   INPAR::FPSI::DirichletNeumann,
                                   INPAR::FPSI::monolithic,
                                   INPAR::FPSI::nocoupling),
                                   &fpsidyn);

    setStringToIntegralParameter<int>("SECONDORDER","No",
                                 "Second order coupling at the interface.",
                                 yesnotuple,yesnovalue,&fpsidyn);

    // Iterationparameters
    setNumericStringParameter("RESTOL","1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
                              "tolerances for single fields in the residual norm for the Newton iteration \n"
                              "for NORM_RESF != *_split only the first value is used for all fields \n"
                              "order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, fluidpressure, ale",
                              &fpsidyn);

    setNumericStringParameter("INCTOL","1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
                              "tolerance in the increment norm for the Newton iteration \n"
                              "for NORM_INC != *_split only the first value is used for all fields \n"
                              "order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, fluidpressure, ale",
                              &fpsidyn);

    setStringToIntegralParameter<int>("NORM_INC","Abs","type of norm for primary variables convergence check \n"
                                 "Abs: absolute values, Abs_sys_split: absolute values with correction of systemsize for every field seperate, Rel_sys: relative values with correction of systemsize",
                                 tuple<std::string>(
                                   "Abs","Abs_sys_split","Rel_sys"
                                   ),
                                 tuple<int>(
                                     INPAR::FPSI::absoluteconvergencenorm,INPAR::FPSI::absoluteconvergencenorm_sys_split,INPAR::FPSI::relativconvergencenorm_sys
                                   ),
                                 &fpsidyn);

    setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for primary variables convergence check \n"
                                   "Abs: absolute values, Abs_sys_split: absolute values with correction of systemsize for every field seperate, Rel_sys: relative values with correction of systemsize",
                                   tuple<std::string>(
                                     "Abs","Abs_sys_split","Rel_sys"
                                     ),
                                   tuple<int>(
                                     INPAR::FPSI::absoluteconvergencenorm,INPAR::FPSI::absoluteconvergencenorm_sys_split,INPAR::FPSI::relativconvergencenorm_sys
                                     ),
                                   &fpsidyn);

    setStringToIntegralParameter<int>("NORMCOMBI_RESFINC","And","binary operator to combine primary variables and residual force values",
                                 tuple<std::string>(
                                       "And",
                                       "Or"),
                                       tuple<int>(
                                         INPAR::FPSI::bop_and,
                                         INPAR::FPSI::bop_or),
                                 &fpsidyn);

    setStringToIntegralParameter<int>("LineSearch","No","adapt increment in case of non-monotonic residual convergence or residual oscillations",
                                 tuple<std::string>(
                                       "Yes",
                                       "No"),
                                       tuple<int>(
                                         1,
                                         0),
                                 &fpsidyn);

    setStringToIntegralParameter<int>("FDCheck","No","perform FPSIFDCheck() finite difference check",
                                 tuple<std::string>(
                                       "Yes",
                                       "No"),
                                       tuple<int>(
                                         1,
                                         0),
                                 &fpsidyn);

    // number of linear solver used for poroelasticity
    IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for FPSI problems",&fpsidyn);

    IntParameter("ITEMAX",10,"maximum number of iterations over fields",&fpsidyn);
    IntParameter("ITEMIN",1,"minimal number of iterations over fields",&fpsidyn);
    IntParameter("NUMSTEP",200,"Total number of Timesteps",&fpsidyn);
    IntParameter("ITEMAX",100,"Maximum number of iterations over fields",&fpsidyn);
    IntParameter("UPRES",1,"Increment for writing solution",&fpsidyn);
    IntParameter("RESTARTEVRY",1,"Increment for writing restart",&fpsidyn);

    IntParameter("FDCheck_row",0,"print row value during FDCheck",&fpsidyn);
    IntParameter("FDCheck_column",0,"print column value during FDCheck",&fpsidyn);

    DoubleParameter("TIMESTEP",0.1,"Time increment dt",&fpsidyn);
    DoubleParameter("MAXTIME",1000.0,"Total simulation time",&fpsidyn);
    DoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields",&fpsidyn);
    DoubleParameter("ALPHABJ",1.0,"Beavers-Joseph-Coefficient for Slip-Boundary-Condition at Fluid-Porous-Interface (0.1-4)",&fpsidyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL",false,"");


  // OUTPUT options
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT_SCREEN","No","Do you want to be informed, if Gmsh output is written?",
                                 yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_SOL_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_EOS_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_DISCRET_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_CUT_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);


  IntParameter("MAX_NUM_DOFSETS",3,"Maximum number of volumecells in the XFEM element",&xfem_general);

  // Integration options
  setStringToIntegralParameter<int>("VOLUME_GAUSS_POINTS_BY","Tessellation","how to find Gauss Points for the cut volumes",
                               tuple<std::string>("Tessellation","MomentFitting","DirectDivergence"),
                               tuple<int>(
                                   INPAR::CUT::VCellGaussPts_Tessellation,
                                   INPAR::CUT::VCellGaussPts_MomentFitting,
                                   INPAR::CUT::VCellGaussPts_DirectDivergence
                                   ),
                               &xfem_general);

  setStringToIntegralParameter<int>("BOUNDARY_GAUSS_POINTS_BY","Tessellation","how to find Gauss Points for the boundary cells",
                                 tuple<std::string>("Tessellation","MomentFitting","DirectDivergence"),
                                 tuple<int>(
                                     INPAR::CUT::BCellGaussPts_Tessellation,
                                     INPAR::CUT::BCellGaussPts_MomentFitting,
                                     INPAR::CUT::BCellGaussPts_DirectDivergence
                                     ),
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
                               tuple<std::string>("interface_vel_by_disp", "interface_vel_by_funct", "interface_vel_by_curve", "interface_vel_zero"),
                               tuple<int>(
                                   INPAR::XFEM::interface_vel_by_disp,    // define interface velocity by displacement of solid
                                   INPAR::XFEM::interface_vel_by_funct,   // define interface velocity by function
                                   INPAR::XFEM::interface_vel_by_curve,   // define interface velocity by curve
                                   INPAR::XFEM::interface_vel_zero        // zero interface velocity function
                                   ),
                               &xfluid_general);

  setStringToIntegralParameter<int>("INTERFACE_DISP","interface_disp_zero","how to define the interface displacement",
                               tuple<std::string>("interface_disp_by_fsi", "interface_disp_by_funct", "interface_disp_by_curve", "interface_disp_zero", "interface_disp_by_implementation"),
                               tuple<int>(
                                   INPAR::XFEM::interface_disp_by_fsi,     // define interface displacement by structure solution of fsi algo
                                   INPAR::XFEM::interface_disp_by_funct,   // define interface displacement by function
                                   INPAR::XFEM::interface_disp_by_curve,   // define interface displacement by curve
                                   INPAR::XFEM::interface_disp_zero,       // zero interface displacement function
                                   INPAR::XFEM::interface_disp_by_implementation // interface displacement by implementation
                                   ),
                               &xfluid_general);

  IntParameter("DISP_FUNCT_NO",-1,"funct number for interface displacement",&xfluid_general);

  setNumericStringParameter("DISP_CURVE_NO","-1 -1 -1",
                              "Curve numbers for interface displacement in (x, y, z) direction",
                              &xfluid_general);

  IntParameter("VEL_FUNCT_NO",-1,"funct number for WDBC or Neumann Condition at embedded boundary/interface",&xfluid_general);
  setNumericStringParameter("VEL_CURVE_NO","-1 -1 -1",
                            "Curve numbers for interface velocity in (x, y, z) direction",
                            &xfluid_general);
  IntParameter("VEL_INIT_FUNCT_NO",-1,"funct number for initial interface velocity",&xfluid_general);

  // How many monolithic steps we keep the fluidfluid-boundary fixed
  IntParameter("RELAXING_ALE_EVERY",1,"Relaxing Ale after how many monolithic steps",&xfluid_general);

  BoolParameter("RELAXING_ALE","yes","switch on/off for relaxing Ale in monolithic fluid-fluid-fsi",&xfluid_general);

  DoubleParameter("XFLUIDFLUID_SEARCHRADIUS",1.0,"Radius of the search tree",&xfluid_general);

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
                                    tuple<std::string>("Xff_TimeInt_FullProj", "Xff_TimeInt_ProjIfMoved","Xff_TimeInt_KeepGhostValues","Xff_TimeInt_IncompProj"),
                                    tuple<int>(
                                      INPAR::XFEM::Xff_TimeInt_FullProj   ,      //always project nodes from embedded to background nodes
                                      INPAR::XFEM::Xff_TimeInt_ProjIfMoved,      //project nodes just if the status of background nodes changed
                                      INPAR::XFEM::Xff_TimeInt_KeepGhostValues,  //always keep the ghost values of the background discretization
                                      INPAR::XFEM::Xff_TimeInt_IncompProj        //after projecting nodes do a incompressibility projection
                                      ),
                                    &xfluid_general);

  setStringToIntegralParameter<int>("XFLUID_TIMEINT","STD=COPY/SL_and_GHOST=COPY/GP","The xfluid time integration approach",
                               tuple<std::string>("STD=COPY_and_GHOST=COPY/GP", "STD=COPY/SL_and_GHOST=COPY/GP", "STD=SL(boundary-zone)_and_GHOST=GP"),
                               tuple<int>(
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP,       // STD= only copy, GHOST= copy or ghost penalty reconstruction
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP, // STD= copy or SL, GHOST= copy or ghost penalty reconstruction
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP         // STD= only SL on whole boundary zone, GHOST= ghost penalty reconstruction
                                   ),
                               &xfluid_general);

  BoolParameter("ALE_XFluid","no","XFluid is Ale Fluid?",&xfluid_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_stab = xfluid_dyn.sublist("STABILIZATION",false,"");

  // Boundary-Coupling options
  setStringToIntegralParameter<int>("COUPLING_METHOD","Nitsche","method how to enforce embedded boundary/coupling conditions at the interface",
                               tuple<std::string>("Hybrid_LM_Cauchy_stress", "Hybrid_LM_viscous_stress", "Nitsche"),
                               tuple<int>(
                                   INPAR::XFEM::Hybrid_LM_Cauchy_stress,  // Cauchy stress-based mixed/hybrid formulation
                                   INPAR::XFEM::Hybrid_LM_viscous_stress, // viscous stress-based mixed/hybrid formulation
                                   INPAR::XFEM::Nitsche                   // Nitsche's formulation
                                   ),
                               &xfluid_stab);


  setStringToIntegralParameter<int>("COUPLING_STRATEGY","xfluid_sided_weak_DBC","enforce weak DBC or coupling conditions with average w.r.t which side?",
                               tuple<std::string>("xfluid_sided_weak_DBC", "xfluid_sided_coupling", "embedded_sided_coupling", "two_sided_coupling"),
                               tuple<int>(
                                   INPAR::XFEM::Xfluid_Sided_weak_DBC,     // one sided weak DBC at the interface w.r.t background mesh
                                   INPAR::XFEM::Xfluid_Sided_Coupling,     // coupling with average w.r.t background mesh
                                   INPAR::XFEM::Embedded_Sided_Coupling,   // coupling with average w.r.t embedded mesh
                                   INPAR::XFEM::Two_Sided_Coupling         // coupling with mean average between background and embedded mesh
                                   ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("HYBRID_LM_L2_PROJ","part_ele_proj","perform the L2 projection between stress fields on whole element or on fluid part?",
                               tuple<std::string>("full_ele_proj", "part_ele_proj"),
                               tuple<int>(
                                   INPAR::XFEM::Hybrid_LM_L2_Proj_full,   // L2 stress projection on whole fluid element
                                   INPAR::XFEM::Hybrid_LM_L2_Proj_part    // L2 stress projection on partial fluid element volume
                                   ),
                               &xfluid_stab);

  BoolParameter("VISC_ADJOINT_SYMMETRY","yes","viscous and adjoint viscous interface terms with matching sign?",&xfluid_stab);

  // viscous and convective Nitsche/MSH stabilization parameter
  DoubleParameter("NIT_STAB_FAC", 35.0, " ( stabilization parameter for Nitsche's penalty term",&xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_TRACE_ESTIMATE","CT_div_by_hk","how to estimate the scaling from the trace inequality in Nitsche's method",
                               tuple<std::string>("CT_div_by_hk", "eigenvalue"),
                               tuple<int>(
                                   INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk,   // estimate the trace inequality by a trace-constant CT and hk: CT/hk
                                   INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue      // estimate the trace inequality by solving a local element-wise eigenvalue problem
                                   ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_HK","ele_vol_div_by_max_ele_surf","how to define the characteristic element length in cut elements",
                                 tuple<std::string>(
                                     "vol_equivalent",
                                     "cut_vol_div_by_cut_surf",
                                     "ele_vol_div_by_cut_surf",
                                     "ele_vol_div_by_ele_surf",
                                     "ele_vol_div_by_max_ele_surf"
                                     ),
                                 tuple<int>(
                                     INPAR::XFEM::ViscStab_hk_vol_equivalent,             /// volume equivalent element diameter
                                     INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf,    /// physical partial/cut volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, not useful for sliver and/or dotted cut situations)
                                     INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf,    /// full element volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, however, avoids problems with sliver cuts, not useful for dotted cuts)
                                     INPAR::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf,    /// full element volume divided by surface measure ( used for uncut situations, standard weak Dirichlet boundary/coupling conditions)
                                     INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf /// default: full element volume divided by maximal element surface measure ( used to estimate the trace inequality for stretched elements in combination with ghost-penalties)
                                 ),
                                 &xfluid_stab);


  setStringToIntegralParameter<int>("CONV_STAB_SCALING","none","scaling factor for viscous interface stabilization (Nitsche, MSH)",
                                    tuple<std::string>("inflow", "abs_inflow", "none"),
                                    tuple<int>(
                                      INPAR::XFEM::ConvStabScaling_inflow,                // scaling with max(0,-u*n)
                                      INPAR::XFEM::ConvStabScaling_abs_inflow,            // scaling with |u*n|
                                      INPAR::XFEM::ConvStabScaling_none                   // no convective stabilization
                                      ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("XFF_CONV_STAB_SCALING","none","scaling factor for convective interface stabilization of fluid-fluid Coupling",
                                    tuple<std::string>("inflow", "averaged", "none"),
                                    tuple<int>(
                                      INPAR::XFEM::XFF_ConvStabScaling_upwinding,          // one-sided inflow stabilization
                                      INPAR::XFEM::XFF_ConvStabScaling_only_averaged,      // averaged inflow stabilization
                                      INPAR::XFEM::XFF_ConvStabScaling_none                // no convective stabilization
                                      ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_COMBO","max","choose the maximum from viscous and convective contributions or just sum both up",
                                    tuple<std::string>("max", "sum"),
                                    tuple<int>(
                                      INPAR::XFEM::MassConservationCombination_max,        /// use the maximum contribution
                                      INPAR::XFEM::MassConservationCombination_sum         /// sum viscous and convective contributions
                                      ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_SCALING","only_visc","apply additional scaling of penalty term to enforce mass conservation for convection-dominated flow",
                                    tuple<std::string>("full", "only_visc"),
                                    tuple<int>(
                                      INPAR::XFEM::MassConservationScaling_full,           /// apply mass-conserving convective scaling additionally
                                      INPAR::XFEM::MassConservationScaling_only_visc       /// use only the viscous scaling
                                      ),
                               &xfluid_stab);

  BoolParameter("GHOST_PENALTY_STAB","no","switch on/off ghost penalty interface stabilization",&xfluid_stab);

  BoolParameter("GHOST_PENALTY_TRANSIENT_STAB","no","switch on/off ghost penalty transient interface stabilization",&xfluid_stab);

  BoolParameter("GHOST_PENALTY_2nd_STAB","no","switch on/off ghost penalty interface stabilization for 2nd order derivatives",&xfluid_stab);

  DoubleParameter("GHOST_PENALTY_FAC",       0.1, "define stabilization parameter ghost penalty interface stabilization",&xfluid_stab);

  DoubleParameter("GHOST_PENALTY_TRANSIENT_FAC",       0.001, "define stabilization parameter ghost penalty transient interface stabilization",&xfluid_stab);

  BoolParameter("XFF_EOS_PRES_EMB_LAYER","no","switch on/off edge-based pressure stabilization on interface-contributing elements of the embedded fluid",&xfluid_stab);

  BoolParameter("IS_PSEUDO_2D","no","modify viscous interface stabilization due to the vanishing polynomial in third dimension when using strong Dirichlet conditions to block polynomials in one spatial dimension",&xfluid_stab);


  /*----------------------------------------------------------------------*/
   Teuchos::ParameterList& particledyn = list->sublist(
       "PARTICLE DYNAMIC",
       false,
       "control parameters for particle problems\n");

   setStringToIntegralParameter<int>("DYNAMICTYP","CentrDiff",
                                "type of time integration control",
                                tuple<std::string>(
                                  "ExplicitEuler",
                                  "CentrDiff",
                                  "RK2",
                                  "RK4"),
                                tuple<int>(
                                  INPAR::PARTICLE::dyna_expleuler,
                                  INPAR::PARTICLE::dyna_centrdiff,
                                  INPAR::PARTICLE::dyna_rk2,
                                  INPAR::PARTICLE::dyna_rk4
                                ),
                                &particledyn);

   // Output type
   IntParameter("RESULTSEVRY",1,"save displacements and contact forces every RESULTSEVRY steps",&particledyn);
   IntParameter("RESEVRYERGY",0,"write system energies every requested step",&particledyn);
   IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&particledyn);
   // Time loop control
   DoubleParameter("TIMESTEP",0.05,"time step size",&particledyn);
   IntParameter("NUMSTEP",200,"maximum number of steps",&particledyn);
   DoubleParameter("MAXTIME",5.0,"maximum time",&particledyn);

   setStringToIntegralParameter<int>(
                               "CONTACT_STRATEGY","None",
                               "Contact strategies for particle problems",
                               tuple<std::string>(
                                 "None",
                                 "NormalContact_DEM",
                                 "NormalContact_MD",
                                 "NormalAndTangentialContact_DEM"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::None,
                                 INPAR::PARTICLE::Normal_DEM,
                                 INPAR::PARTICLE::Normal_MD,
                                 INPAR::PARTICLE::NormalAndTang_DEM
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "NORMAL_CONTACT_LAW","LinearSpringDamp",
                               "contact law for normal contact of particles",
                               tuple<std::string>(
                                 "LinearSpring",
                                 "Hertz",
                                 "LinearSpringDamp",
                                 "LeeHerrmann",
                                 "KuwabaraKono",
                                 "Tsuji"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::LinSpring,
                                 INPAR::PARTICLE::Hertz,
                                 INPAR::PARTICLE::LinSpringDamp,
                                 INPAR::PARTICLE::LeeHerrmann,
                                 INPAR::PARTICLE::KuwabaraKono,
                                 INPAR::PARTICLE::Tsuji
                                 ),
                               &particledyn);
   DoubleParameter("MIN_RADIUS",-1.0,"smallest particle radius",&particledyn);
   DoubleParameter("MAX_RADIUS",-1.0,"largest particle radius",&particledyn);
   DoubleParameter("REL_PENETRATION",-1.0,"relative particle penetration",&particledyn);
   DoubleParameter("MAX_VELOCITY",-1.0,"highest particle velocity",&particledyn);
   DoubleParameter("COEFF_RESTITUTION",-1.0,"coefficient of restitution",&particledyn);
   DoubleParameter("COEFF_RESTITUTION_WALL",-1.0,"coefficient of restitution (wall)",&particledyn);
   DoubleParameter("FRICT_COEFF_WALL",-1.0,"friction coefficient for contact particle-wall",&particledyn);
   DoubleParameter("FRICT_COEFF",-1.0,"dynamic friction coefficient for contact particle-particle",&particledyn);
   DoubleParameter("NORMAL_STIFF",-1.0,"stiffness for normal contact force",&particledyn);
   DoubleParameter("TANG_STIFF",-1.0,"stiffness for tangential contact force",&particledyn);
   DoubleParameter("NORMAL_DAMP",-1.0,"damping coefficient for normal contact force",&particledyn);
   DoubleParameter("TANG_DAMP",-1.0,"damping coefficient for tangential contact force",&particledyn);
   DoubleParameter("NORMAL_DAMP_WALL",-1.0,"damping coefficient for normal contact force (wall)",&particledyn);
   DoubleParameter("TANG_DAMP_WALL",-1.0,"damping coefficient for tangential contact force (wall)",&particledyn);
   BoolParameter("TENSION_CUTOFF","no","switch on/off tension cutoff",&particledyn);
   BoolParameter("MOVING_WALLS","no","switch on/off moving walls",&particledyn);
   DoubleParameter("RANDOM_AMPLITUDE",0.0,"random value for initial position",&particledyn);

   setNumericStringParameter("GRAVITY_ACCELERATION","0.0 0.0 0.0",
                             "Acceleration due to gravity in particle/cavitation simulations.",
                             &particledyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cavitationdyn = list->sublist(
      "CAVITATION DYNAMIC",
      false,
      "control parameters for cavitation problems\n");

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&cavitationdyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&cavitationdyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&cavitationdyn);
  IntParameter("UPRES",1,"Increment for writing solution",&cavitationdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&cavitationdyn);

  // Coupling strategy
  setStringToIntegralParameter<int>(
                              "COUPALGO","cavitation_twowaymomentum",
                              "Coupling strategies for cavitation problems",
                              tuple<std::string>(
                                "cavitation_oneway",
                                "cavitation_twowaymomentum",
                                "cavitation_twowayfull"
                                ),
                              tuple<int>(
                                INPAR::CAVITATION::OneWay,
                                INPAR::CAVITATION::TwoWayMomentum,
                                INPAR::CAVITATION::TwoWayFull
                                ),
                              &cavitationdyn);
  setStringToIntegralParameter<int>(
                              "VOID_FRACTION_CALC","analytical_quadraticpoly",
                              "Void fraction calculation strategy",
                              tuple<std::string>(
                                "analytical_constpoly",
                                "analytical_quadraticpoly",
                                "gaussian_integration"
                                ),
                              tuple<int>(
                                INPAR::CAVITATION::analytical_constpoly,
                                INPAR::CAVITATION::analytical_quadraticpoly,
                                INPAR::CAVITATION::gaussian_integration
                                ),
                              &cavitationdyn);
  IntParameter("NUM_GP_VOID_FRACTION",4,"Number of gauss points in each direction for void fraction calculation",&cavitationdyn);

  BoolParameter("APPROX_ELECOORDS_INIT","no","switch on/off approximate initial guess for computing element coordinates",&cavitationdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& crackdyn = list->sublist("CRACK",false,"");

  // type of crack propagation model -- either using linear elastic fracture mechanics concepts, or cohesive crack models
  setStringToIntegralParameter<int>("CRACK_MODEL","none",
                                      "type of crack propagation modeling",
                                      tuple<std::string>(
                                        "none",
                                        "lefm",
                                        "cohesive"),
                                      tuple<int>(
                                        INPAR::CRACK::crack_none,
                                        INPAR::CRACK::crack_lefm,
                                        INPAR::CRACK::crack_cohesive),
                                      &crackdyn);

  DoubleParameter("CRITICAL_K1",50000.0,"Critical stress intensity factor in normal mode",&crackdyn);

  DoubleParameter("CRITICAL_K2",50000.0,"Critical stress intensity factor in shear mode",&crackdyn);

  StringParameter("THICKNESS_ASSUMPTION","plane_strain",
                  "Whether is this plane strain or plane stress problem?",
                   &crackdyn);

  // type of cohesive crack propagation modeling
  setStringToIntegralParameter<int>("COHESIVE_CRACK_MODEL","none",
                                  "type of cohesive crack propagation modeling",
                                  tuple<std::string>(
                                    "none",
                                    "dczm",
                                    "ddzm"),
                                  tuple<int>(
                                    INPAR::CRACK::cohesive_none,
                                    INPAR::CRACK::cohesive_dczm,
                                    INPAR::CRACK::cohesive_ddzm),
                                  &crackdyn);

  setStringToIntegralParameter<int>("TRACTION_SEPARATION_LAW","exponential",
                                    "type of traction-separation law used for cohesive elements",
                                    tuple<std::string>(
                                      "linear",
                                      "trapezoidal",
                                      "exponential",
                                      "sinusoidal",
                                      "ppr"),
                                    tuple<int>(
                                      INPAR::CRACK::linear,
                                      INPAR::CRACK::trapezoidal,
                                      INPAR::CRACK::exponential,
                                      INPAR::CRACK::sinusoidal,
                                      INPAR::CRACK::ppr),
                                    &crackdyn);

  // are we modeling cracks with known propagation direction?
  setStringToIntegralParameter<int>("IS_CRACK_PREDEFINED","No","Have you already predefined the crack path?",
                               yesnotuple,yesnovalue,&crackdyn);

  DoubleParameter("NORMAL_COHESIVE_STRENGTH",50000.0,"Cohesive strength for normal separation",&crackdyn);
  DoubleParameter("SHEAR_COHESIVE_STRENGTH",500000000.0,"Cohesive strength for shear separation",&crackdyn);
  DoubleParameter("G_I",1.0,"Model I fracture energy (normal separation)",&crackdyn);
  DoubleParameter("G_II",1.0,"Model II fracture energy (shear separation)",&crackdyn);

  DoubleParameter("ALFA_PPR",3.0,"Constant alpha in PPR model",&crackdyn);
  DoubleParameter("BETA_PPR",3.0,"Constant beta in PPR model",&crackdyn);

  DoubleParameter("SLOPE_NORMAL",0.02,"Initial slope indicator in normal direction for PPR model",&crackdyn);
  DoubleParameter("SLOPE_SHEAR",0.02,"Initial slope indicator in normal direction for PPR model",&crackdyn);

  setStringToIntegralParameter<int>("GMSH_OUT","No","Do you want to write Gmsh output of displacement each timestep?",
                                 yesnotuple,yesnovalue,&crackdyn);

  IntParameter("START_NEW_NODE_ID",0,"Id of first node that will be introduced into discretization while propagating crack. This "
      "should be set greater than total no of nodes in the initial discretization",&crackdyn);

  IntParameter("START_NEW_ELE_ID",0,"Id of first wedge element that will be introduced into discretization while propagating crack. This "
        "should be set greater than total no of elements in the initial discretization",&crackdyn);

  // type of crack propagation model -- either using linear elastic fracture mechanics concepts, or cohesive crack models
  setStringToIntegralParameter<int>("CRACK_PROPAGATION_CRITERION","displacement_correlation",
                                      "Crack propagation criterion used for LEFM",
                                      tuple<std::string>(
                                        "displacement_correlation",
                                        "InteractionIntegral"),
                                      tuple<int>(
                                        INPAR::CRACK::displacementCorrelation,
                                        INPAR::CRACK::InteractionIntegral),
                                      &crackdyn);

  IntParameter("NO_LAYERS_J_INT",4,"No of element layers used for J-integral calculation",&crackdyn);

  DoubleParameter("CRITICAL_J",1.0,"Critical energy release rate",&crackdyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fsi_crackdyn = list->sublist("FSI CRACK",false,"");

  setStringToIntegralParameter<int>("CHECK_CONDITION","No","Do you want to check crack mouth opening condition?",
                                     yesnotuple,yesnovalue,&fsi_crackdyn);

  DoubleParameter("CRACK_OPENING_DIST",0.01,"Critial crack mouth opening distance",&fsi_crackdyn);



  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& acousticdyn = list->sublist("ACOUSTIC DYNAMIC",false,"control parameters for acoustic or photoacoustic problems\n");

  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&acousticdyn);
  IntParameter("NUMSTEP",100,"Total number of time steps",&acousticdyn);
  DoubleParameter("MAXTIME",1.0,"Total simulation time",&acousticdyn);
  setStringToIntegralParameter<int>("CALCERROR","No",
                               "compute error compared to analytical solution",
                               tuple<std::string>(
                                 "No",
                                 "1D",
                                 "2D"
                                 ),
                               tuple<int>(
                                   INPAR::ACOU::calcerror_no,
                                   INPAR::ACOU::calcerror_1d,
                                   INPAR::ACOU::calcerror_2d
                                   ),
                               &acousticdyn);


  IntParameter("UPRES",1,"Increment for writing solution",&acousticdyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&acousticdyn);
  IntParameter("LINEAR_SOLVER",-1,"Number of linear solver used for acoustical problem",&acousticdyn);
  IntParameter("STARTFUNCNO",-1,"Function for Initial Starting Field",&acousticdyn);
  IntParameter("SOURCETERMFUNCNO",-1,"Function for source term in volume",&acousticdyn);

  // distinguish viscous and lossless flows
  setStringToIntegralParameter<int>("PHYSICAL_TYPE","lossless",
                    "fluid properties",
                    tuple<std::string>(
                    "lossless",
                    "viscous"),
                    tuple<int>(
                    INPAR::ACOU::acou_lossless,
                    INPAR::ACOU::acou_viscous),
                    &acousticdyn);


  // for nonlinear time integration
  DoubleParameter("CONVTOL",1.0e-9,"Convergence tolerance for Newton loop",&acousticdyn);
  IntParameter("ITEMAX",10,"Maximum number of iterations for Newton loop",&acousticdyn);

  // photoacoustics
  DoubleParameter("PULSEDURATION",15e-9,"Laser pulse duration",&acousticdyn);
  BoolParameter("PHOTOACOU","No","Coupling with Scatra for Diffusive Light Transport",&acousticdyn);
  BoolParameter("MESHCONFORM","No","Conformity of scatra and acoustical mesh",&acousticdyn);

  // local postprocessing and p-adaptivity
  BoolParameter("ERRORMAPS","No","Output of error maps obtained by local postprocessing",&acousticdyn);
  BoolParameter("P_ADAPTIVITY","No","p-adaptivity in time integration",&acousticdyn);
  DoubleParameter("P_ADAPT_TOL",1.0e-15,"Error tolerance for p-adaptivity",&acousticdyn);

  // time integration
  setStringToIntegralParameter<int>("TIMEINT","impl",
                    "Type of time integration scheme",
                    tuple<std::string>(
                    "impl",
                    "trap",
                    "dirk23",
                    "dirk33",
                    "dirk34",
                    "dirk54",
                    "bdf2",
                    "bdf3",
                    "bdf4"),
                    tuple<int>(
                    INPAR::ACOU::acou_impleuler,
                    INPAR::ACOU::acou_trapezoidal,
                    INPAR::ACOU::acou_dirk23,
                    INPAR::ACOU::acou_dirk33,
                    INPAR::ACOU::acou_dirk34,
                    INPAR::ACOU::acou_dirk54,
                    INPAR::ACOU::acou_bdf2,
                    INPAR::ACOU::acou_bdf3,
                    INPAR::ACOU::acou_bdf4),
                    &acousticdyn);

  setStringToIntegralParameter<int>("INV_ANALYSIS","none",
                 "Types of inverse analysis and on/off switch",
                 tuple<std::string>(
                   "none",
                   "pat"), // here, backprojection could be added
                 tuple<int>(
                   INPAR::ACOU::inv_none,
                   INPAR::ACOU::inv_pat),
                 &acousticdyn);


  Teuchos::ParameterList& acou_inv = acousticdyn.sublist("PA IMAGE RECONSTRUCTION",false,"");

  setStringToIntegralParameter<int>("OPTIMIZATION","LBFGS",
                                    "types of optimization algorithm",
                                    tuple<std::string>(
                                      "GradientDescent",
                                      "LBFGS"),
                                    tuple<int>(
                                      INPAR::ACOU::inv_gd,
                                      INPAR::ACOU::inv_lbfgs),
                                    &acou_inv);

  StringParameter("MONITORFILE","none.monitor","Filename of file containing measured pressure values",&acou_inv);
  BoolParameter("FDCHECK","No","Finite difference check",&acou_inv);
  DoubleParameter("INV_TOL",1e-16,"Tolerance for objective function of inverse pat analysis",&acou_inv);
  BoolParameter("INV_TOL_GRAD_YN","No","Flag to indicate check of the norm of the gradient",&acou_inv);
  DoubleParameter("INV_TOL_GRAD",0.0,"Tolerance for norm of gradient of inverse pat analysis",&acou_inv);
  BoolParameter("ELE_SCALING","No","Should gradient be scaled by element size?",&acou_inv);
  IntParameter("INV_MAX_RUN",10,"Maximal run number for inverse pat analysis",&acou_inv);
  IntParameter("INV_LS_MAX_RUN",10,"Maximal run number for line search in inverse pat analysis",&acou_inv);
  DoubleParameter("LS_DECREASECOND",0.0,"coefficient for calculation of sufficient decrease condition",&acou_inv);
  DoubleParameter("LS_STEPLENGTHRED",0.5,"step length is multiplied by this value if line search not yet sufficient",&acou_inv);
  DoubleParameter("ALPHA_MUA",0.0,"Regularization parameter for absorption coefficient",&acou_inv);
  DoubleParameter("BETA_MUA",0.0,"Regularization parameter for gradient of absorption coefficient",&acou_inv);
  IntParameter("SIZESTORAGE",10,"number of vectors to keep in storage; defaults to 10 (lbfgs usage only)",&acou_inv);

  // decide which parametrization of material parameters to use
  setStringToIntegralParameter<int>("PARAMETRIZATION","none",
                                      "how to parametrize the parameter field",
                                    tuple<std::string>(
                                      "none",
                                      "elementwise",
                                      "uniform"),
                                    tuple<int>(
                                      INPAR::ACOU::inv_mat_none,
                                      INPAR::ACOU::inv_mat_elementwise,
                                      INPAR::ACOU::inv_mat_uniform),
                                    &acou_inv);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& nlnsol = list->sublist("NONLINEAR SOLVER", false, "Configuration of nonlinear solver package");

  StringParameter("XML_FILE", "none", "Filename of XML file with configuration of nonlinear solver", &nlnsol);

  /*----------------------------------------------------------------------*/
  // set valid parameters for solver blocks

  // Note: the maximum number of solver blocks is hardwired here. If you change this,
  // don't forget to edit the corresponding parts in globalproblems.cpp, too.
  for (int i = 1; i<10; i++)
  {
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

  Teuchos::Tuple<std::string,12> solver_name;
  Teuchos::Tuple<int,12>  solver_number;

  solver_name[0] = "Amesos_KLU_sym";               solver_number[0] = INPAR::SOLVER::amesos_klu_sym;
  solver_name[1] = "Amesos_KLU_nonsym";            solver_number[1] = INPAR::SOLVER::amesos_klu_nonsym;
  solver_name[2] = "Superlu";                      solver_number[2] = INPAR::SOLVER::superlu;
  solver_name[3] = "Aztec_MSR";                    solver_number[3] = INPAR::SOLVER::aztec_msr;
  solver_name[4] = "LAPACK_sym";                   solver_number[4] = INPAR::SOLVER::lapack_sym;
  solver_name[5] = "LAPACK_nonsym";                solver_number[5] = INPAR::SOLVER::lapack_nonsym;
  solver_name[6] = "UMFPACK";                      solver_number[6] = INPAR::SOLVER::umfpack;
  solver_name[7] = "Belos";                        solver_number[7] = INPAR::SOLVER::belos;
  solver_name[8] = "Stratimikos_Amesos";           solver_number[8] = INPAR::SOLVER::stratimikos_amesos;
  solver_name[9] = "Stratimikos_Aztec";            solver_number[9]= INPAR::SOLVER::stratimikos_aztec;
  solver_name[10]= "Stratimikos_Belos";            solver_number[10]= INPAR::SOLVER::stratimikos_belos;
  solver_name[11]= "undefined";                    solver_number[11]= INPAR::SOLVER::undefined;


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
    Teuchos::Tuple<std::string,28> name;
    Teuchos::Tuple<int,28>  number;

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
    name[10] = "ML";                          number[10] = INPAR::SOLVER::azprec_ML;
    name[11] = "MLFLUID";                     number[11] = INPAR::SOLVER::azprec_MLfluid;
    name[12] = "MLFLUID2";                    number[12] = INPAR::SOLVER::azprec_MLfluid2;
    name[13] = "MLAPI";                       number[13] = INPAR::SOLVER::azprec_MLAPI;
    name[14] = "GaussSeidel";                 number[14] = INPAR::SOLVER::azprec_GaussSeidel;
    name[15] = "DownwindGaussSeidel";         number[15] = INPAR::SOLVER::azprec_DownwindGaussSeidel;
    name[16] = "BGS2x2";                      number[16] = INPAR::SOLVER::azprec_BGS2x2;
    name[17] = "BGSnxn";                      number[17] = INPAR::SOLVER::azprec_BGSnxn;
    name[18] = "TekoSIMPLE";                  number[18] = INPAR::SOLVER::azprec_TekoSIMPLE;
    name[19] = "CheapSIMPLE";                 number[19] = INPAR::SOLVER::azprec_CheapSIMPLE;
    name[20] = "MueLu_sym";                   number[20] = INPAR::SOLVER::azprec_MueLuAMG_sym;
    name[21] = "MueLu_nonsym";                number[21] = INPAR::SOLVER::azprec_MueLuAMG_nonsym;
    name[22] = "MueLu_contact";               number[22] = INPAR::SOLVER::azprec_MueLuAMG_contact;
    name[23] = "MueLu_contact2";              number[23] = INPAR::SOLVER::azprec_MueLuAMG_contact2;
    name[24] = "MueLu_contact3";              number[24] = INPAR::SOLVER::azprec_MueLuAMG_contact3;
    name[25] = "MueLu_contactSP";             number[25] = INPAR::SOLVER::azprec_MueLuAMG_contactSP;
    name[26] = "MueLu_contactPenalty";        number[26] = INPAR::SOLVER::azprec_MueLuAMG_contactPen;
    name[27] = "AMGnxn";                      number[27] = INPAR::SOLVER::azprec_AMGnxn;

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
    "AZSUB", 50,
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
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS","Umfpack","BS","SIMPLE","SIMPLEC","IBD","Uzawa"),
    tuple<int>(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERMED","ILU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS","Umfpack","BS","SIMPLE","SIMPLEC","IBD","Uzawa"),
    tuple<int>(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),
    &list);

  setStringToIntegralParameter<int>(
    "ML_SMOOTHERCOARSE","Umfpack","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu","GS","DGS","Umfpack","BS","SIMPLE","SIMPLEC","IBD","Uzawa"),
    tuple<int>(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),
    &list);

  // TODO remove this
  setStringToIntegralParameter<int>(
    "MueLu_INITSMOOTHER","SGS","",
    tuple<std::string>("SGS","Jacobi","Chebychev","ILU","GS"),
    tuple<int>(0,1,2,4,7),
    &list);

  IntParameter("SUB_SOLVER1", -1  ,"sub solver/smoother block number (SIMPLE/C: used for prediction of primary variable on all levels, BS: used for fine and intermedium BraessSarazin (BS) level smoother)",&list);
  IntParameter("SUB_SOLVER2", -1  ,"sub solver/smoother block number (SIMPLE/C: used for SchurComplement eq. on all levels, BS: used for coarse BraessSarazin (BS) level smoother)",&list);

  // TODO remove this
  IntParameter("MueLu_INITSMOO_SWEEPS", 1  ,"number of sweeps for adaptive SA smoother (initialization phase). For Chebyshev it is used as polynomial degree",&list);
  DoubleParameter("MueLu_INITSMOO_DAMPING",1.,"damping parameter for adaptive SA smoother (initialization phase). For Chebyshev it is used as alpha parameter",&list);

  IntParameter("MueLu_MIN_AGG_SIZE",6,
               "Minimal objective size of an aggregate (to influence the coarsening rate)",&list);

  setStringToIntegralParameter<int>(
    "MueLu_REBALANCE","No","activate rebalancing using Zoltan/Isorropia",
    tuple<std::string>("NO", "No","no","YES","Yes","yes"),
    tuple<int>(0,1,2,3,4,5),
    &list);
  DoubleParameter("MueLu_REBALANCE_NONZEROIMBALANCE",1.2,"maximum allowed nonzero imbalance factor",&list);
  IntParameter("MueLu_REBALANCE_MINROWS", 1000  ,"minimum numbers of rows per processor before rebalancing is necessary",&list);

  {
    // MueLu_Reuse
    Teuchos::Tuple<std::string,3> name;
    Teuchos::Tuple<int,3>  number;

    name[0] = "nothing";                      number[0] = INPAR::SOLVER::Reuse_nothing;
    name[1] = "Ptent";                        number[1] = INPAR::SOLVER::Reuse_Ptent;
    name[2] = "full";                         number[2] = INPAR::SOLVER::Reuse_full;

    setStringToIntegralParameter<int>(
      "MueLu_REUSE", "nothing",
      "Reuse strategy in MueLu contact preconditioner 2.\n"
      "The options are 'none', 'Ptent' and 'full'.\n"
      "'full' means: reuse the full multigrid hierarchy.\n"
      "'Ptent': reuse aggregates and nonsmoothed transfer operator.\n"
      "The MueLu_Reuse parameter only makes sense together with AZREUSE.\n",
      name,
      number,
      &list
      );
  }

//  // parameters for AMG(BS)
//  setNumericStringParameter("AMGBS_BS_DAMPING","1.3 1.3 1.3",
//                            "Relaxation factor for Braess-Sarazin smoother within AMGBS method",
//                            &list);
//
//  setNumericStringParameter("AMGBS_BS_PCSWEEPS","2 2 2",
//                            "number of jacobi/sgs sweeps for smoothing/solving pressure correction equation within Braess-Sarazin. only necessary for jacobi/gauss seidel",
//                            &list);
//
//  setNumericStringParameter("AMGBS_BS_PCDAMPING","1.0 1.0 1.0",
//                              "jacobi damping factors for smoothing/solving pressure correction equation within Braess-Sarazin. only necessary for jacobi/gauss seidel",
//                              &list);
//
//
//  setStringToIntegralParameter<int>(
//    "AMGBS_PSMOOTHER_VEL","PA-AMG","Prolongation/Restriction smoothing strategy (velocity part in AMGBS preconditioner)",
//    tuple<std::string>("PA-AMG","SA-AMG","PG-AMG","PG2-AMG"),
//    tuple<int>(INPAR::SOLVER::PA_AMG,INPAR::SOLVER::SA_AMG,INPAR::SOLVER::PG_AMG,INPAR::SOLVER::PG2_AMG),
//    &list);
//  setStringToIntegralParameter<int>(
//    "AMGBS_PSMOOTHER_PRE","PA-AMG","Prolongation/Restriction smoothing strategy (pressure part in AMGBS preconditioner)",
//    tuple<std::string>("PA-AMG","SA-AMG","PG-AMG","PG2-AMG"),
//    tuple<int>(INPAR::SOLVER::PA_AMG,INPAR::SOLVER::SA_AMG,INPAR::SOLVER::PG_AMG,INPAR::SOLVER::PG2_AMG),
//    &list);
//
//  setStringToIntegralParameter<int>(
//    "AMGBS_BS_PCCOARSE","Umfpack","approximation algorithm for solving pressure correction equation (coarsest level)",
//    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel stand-alone"),
//    tuple<int>(0,1,2,3,4,5,6,7,8),
//    &list);
//
//  setStringToIntegralParameter<int>(
//    "AMGBS_BS_PCMEDIUM","Umfpack","approximation algorithm for solving pressure correction equation (medium level)",
//    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel stand-alone"),
//    tuple<int>(0,1,2,3,4,5,6,7,8),
//    &list);
//
//  setStringToIntegralParameter<int>(
//    "AMGBS_BS_PCFINE","Umfpack","approximation algorithm for solving pressure correction equation (finest level)",
//    tuple<std::string>("Umfpack","KLU","ILU","Jacobi","Gauss-Seidel","symmetric Gauss-Seidel","Jacobi stand-alone","Gauss-Seidel stand-alone","symmetric Gauss-Seidel stand-alone"),
//    tuple<int>(0,1,2,3,4,5,6,7,8),
//    &list);

  // switch order of blocks in BGS2x2 preconditioner
  setStringToIntegralParameter<int>(
    "BGS2X2_FLIPORDER","block0_block1_order","BGS2x2 flip order parameter",
    tuple<std::string>("block0_block1_order","block1_block0_order"),
    tuple<int>(0,1),
    &list);

  // damping parameter for BGS2X2
  DoubleParameter("BGS2X2_GLOBAL_DAMPING",1.,"damping parameter for BGS2X2 preconditioner",&list);
  DoubleParameter("BGS2X2_BLOCK1_DAMPING",1.,"damping parameter for BGS2X2 preconditioner block1",&list);
  DoubleParameter("BGS2X2_BLOCK2_DAMPING",1.,"damping parameter for BGS2X2 preconditioner block2",&list);

  // parameters for permutation of linear systems
  {
    Teuchos::Tuple<std::string,5> name;
    Teuchos::Tuple<int,5>  number;
    name[0] = "none";         number[0] = INPAR::SOLVER::Permutation_none;
    name[1] = "Algebraic";    number[1] = INPAR::SOLVER::Permutation_algebraic;
    name[2] = "algebraic";    number[2] = INPAR::SOLVER::Permutation_algebraic;
    name[3] = "Local";        number[3] = INPAR::SOLVER::Permutation_local;
    name[4] = "local";        number[4] = INPAR::SOLVER::Permutation_local;

    setStringToIntegralParameter<int>(
      "PERMUTE_SYSTEM", "none",
      "allow linear solver to permute linear system to improve properties of linear system for iterative methods.",
      name,
      number,
      &list
      );
  }
  DoubleParameter("NON_DIAGDOMINANCE_RATIO",1.,"matrix rows with diagEntry/maxEntry<nonDiagDominanceRatio are marked to be significantly non-diagonal dominant (default: 1.0 = mark all non-diagonal dominant rows)",&list);

  // verbosity flag (for Belos)
  IntParameter("VERBOSITY",0,"verbosity level (0=no output,... 10=extreme), for Belos only",&list);

  // the only one stratimikos specific parameter
  StringParameter("STRATIMIKOS_XMLFILE","none",
                  "xml file for stratimikos parameters",
                  &list);

  // user-given name of solver block (just for beauty)
  StringParameter("NAME","No_name",
                  "User specified name for solver block",
                  &list);

  // damping parameter for SIMPLE
  DoubleParameter("SIMPLE_DAMPING",1.,"damping parameter for SIMPLE preconditioner",&list);



  // Parameters for AMGnxn Preconditioner
  StringParameter("AMGNXN_TYPE","AMG(BGS)", "Name of the pre-built preconditioner to be used. If set to\"XML\" the preconditioner is defined using a xml file", &list);
  StringParameter("AMGNXN_XML_FILE","none", "xml file defining the AMGnxn preconditioner", &list);


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
  DoubleParameter("SIZERATIOSCALE", 0.9, "This is a safety factor to scale theoretical optimal step size, should be lower than 1 and must be larger than 0", &list);

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
void DRT::INPUT::SetValidRandomFieldParameters(Teuchos::ParameterList& list)
{
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  // define some tuples that are often used to account for different writing of certain key words
  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

   // safety measure to prevent using random fields that were not properly setup or meant to be used at all
  setStringToIntegralParameter<int>("ACTIVE","No",
                                    "Do we want to use this field ?  ",
                                    yesnotuple,yesnovalue,&list);

  // Parameters to simulate random fields
  IntParameter("RANDOM_FIELD_DIMENSION",3,"Dimension of Random Field 2 or 3",&list);
  IntParameter("SIZE_PER_DIM",512,"Number of points per dimension",&list);
  DoubleParameter("CORRLENGTH",30,"Correlation length of Random Field",&list);
  IntParameter("NUM_COS_TERMS",64,"Number of terms in geometric row ",&list);
  setStringToIntegralParameter<int>("SPECTRAL_MATCHING","Yes",
                                   "Perform spectral matching with psd 1 yes ",
                                   yesnotuple,yesnovalue,&list);
   DoubleParameter("SIGMA",1.0,"sigma of random field",&list);
   DoubleParameter("MEAN",0.0,"Mean value of random field",&list);
   setStringToIntegralParameter<int>("CORRSTRUCT","gaussian",
                               "Correlation structure of random field",
                               tuple<std::string>("Gaussian","gaussian"),
                               tuple<int>(
                                   INPAR::MLMC::corr_gaussian,INPAR::MLMC::corr_gaussian),
                               &list);
   setStringToIntegralParameter<int>("MARGINALPDF","gaussian","Target marginal probability distribution function",
                                     tuple<std::string>("Gaussian","gaussian",
                                                         "Beta", "beta",
                                                         "lognormal", "Lognormal"),
                                     tuple<int>(
                                         INPAR::MLMC::pdf_gaussian,INPAR::MLMC::pdf_gaussian,
                                         INPAR::MLMC::pdf_beta,INPAR::MLMC::pdf_beta,
                                         INPAR::MLMC::pdf_lognormal,INPAR::MLMC::pdf_lognormal),
                                  &list);
   DoubleParameter("NONGAUSSPARAM1",0,"First parameter for non-gaussian pdf",&list);
   DoubleParameter("NONGAUSSPARAM2",0,"Second parameter for non-gaussian pdf",&list);
   DoubleParameter("KAPPA_U",6.283185307,"CUTOFF WAVE NUMBER FOR PSD",&list);
   setStringToIntegralParameter<int>("CALC_METHOD","FFT",
                                 "Calculation method for the random field",
                                 tuple<std::string>("FFT","COS","FOURIER"),
                                 tuple<int>(
                                     INPAR::MLMC::calc_m_fft,INPAR::MLMC::calc_m_cos,INPAR::MLMC::calc_m_fourier),
                                 &list);

   DoubleParameter("PERIODICITY_FOURIER",1.0,"Periodic Length of Random Field in case of Fourier Series Expansion",&list);
   // For testing
   setStringToIntegralParameter<int>("USEDETVALUE","No",
                                      "Instead of doing proper MC simulation use DETVALUE for the stochastic parameter",
                                      yesnotuple,yesnovalue,&list);
   DoubleParameter("CONTBLENDVALUE",1.5,"Use this values for parameter continuation",&list);

   IntParameter("FOURIER_TRUNCATION_THRESHOLD",130,"Truncation threshold for fourier series expansion",&list);

   // define cutoff values
   setStringToIntegralParameter<int>("BOUNDED","No",
                                    "Cutoff Randomfield to prevent unrealistically low or high values",
                                    yesnotuple,yesnovalue,&list);
   DoubleParameter("LOWERBOUND",1.5,"Lower cutoff value",&list);
   DoubleParameter("UPPERBOUND",1.5,"Uower cutoff value",&list);
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

