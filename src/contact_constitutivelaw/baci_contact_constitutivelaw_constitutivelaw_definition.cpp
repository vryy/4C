/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of contact constitutive law definitions

\level 3


*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* headers */

#include "baci_contact_constitutivelaw_constitutivelaw_definition.hpp"

#include "baci_contact_constitutivelaw_bundle.hpp"
#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "baci_io_inputreader.hpp"
#include "baci_io_line_parser.hpp"
#include "baci_utils_exceptions.hpp"


FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/*======================================================================*/
CONTACT::CONSTITUTIVELAW::LawDefinition::LawDefinition(
    std::string name, std::string description, INPAR::CONTACT::ConstitutiveLawType type)
    : coconstlawname_(name), description_(description), coconstlawtype_(type)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LawDefinition::AddComponent(Teuchos::RCP<INPUT::LineComponent> c)
{
  inputline_.push_back(c);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LawDefinition::Read(const GLOBAL::Problem& problem,
    INPUT::DatFileReader& reader, Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Bundle> bundle)
{
  std::string name = "--CONTACT CONSTITUTIVE LAWS";
  std::vector<const char*> section = reader.Section(name);

  if (section.size() > 0)
  {
    for (auto& i : section)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(i));

      // add trailing white space to stringstream "condline" to avoid deletion of stringstream upon
      // reading the last entry inside This is required since the material parameters can be
      // specified in an arbitrary order in the input file. So it might happen that the last entry
      // is extracted before all of the previous ones are.
      condline->seekp(0, condline->end);
      *condline << " ";

      IO::LineParser parser("While reading 'CONTACT CONSTITUTIVE LAWS' section: ");

      parser.Consume(*condline, "LAW");
      const int id = parser.Read<int>(*condline);
      const std::string name = parser.Read<std::string>(*condline);

      // Remove the parts that were already read.
      condline->str(condline->str().erase(0, (size_t)condline->tellg()));

      if (name == coconstlawname_)
      {
        if (id <= -1) dserror("Illegal negative ID provided");

        // check if material ID is already in use
        if (bundle->Find(id) != -1)
          dserror("More than one contact constitutivelaw with 'Law %d'", id);

        // the read-in contact constitutive law line
        Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container = Teuchos::rcp(
            new CONTACT::CONSTITUTIVELAW::Container(id, coconstlawtype_, coconstlawname_));
        // fill the latter

        for (const auto& j : inputline_) condline = j->Read(Name(), condline, *container);

        // current material input line contains bad elements
        if (condline->str().find_first_not_of(' ') != std::string::npos)
          dserror(
              "Specification of material '%s' contains the following unknown, redundant, or "
              "incorrect elements: '%s'",
              coconstlawname_.c_str(), condline->str().c_str());

        // put contact constitutive law in map of laws
        bundle->Insert(id, container);
      }
    }
  }
}


std::ostream& CONTACT::CONSTITUTIVELAW::LawDefinition::Print(
    std::ostream& stream, const DRT::Discretization* dis)
{
  // a string holding the comment indicating symbols for DAT input file
  const std::string comment = "//";

  // the descriptive lines (comments)
  stream << comment << std::endl;
  stream << comment << " " << description_ << std::endl;
  for (const auto& i : inputline_)
  {
    std::ostringstream desc;
    i->Describe(desc);
    if (desc.str().size() > 0) stream << comment << desc.str() << std::endl;
  }

  // the default line
  stream << comment << "LAW 0   " << coconstlawname_ << "   ";
  for (const auto& i : inputline_)
  {
    i->DefaultLine(stream);
    stream << " ";
  }

  stream << "\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::AppendCoConstLawComponentDefinition(
    std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& list,
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> def)
{
  // test if material was defined with same name or type
  std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>::const_iterator m;
  for (m = list.begin(); m != list.end(); ++m)
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> mmd = *m;

    if (mmd->Type() == def->Type())
      dserror(
          "Trying to define two contact constitutive laws with the same type '%d'\n"
          "Please revise your definitions of valid contact constitutive laws",
          mmd->Type());

    if (mmd->Name() == def->Name())
      dserror(
          "Trying to define two contact constitutive laws with the same name '%s'\n"
          "Please revise your definitions of valid contact constitutive laws",
          mmd->Name().c_str());
  }

  // no coincidence found
  if (m == list.end())
    list.push_back(def);
  else
    dserror("Trouble in determining coincidences of contact constitutive definitions");
}

FOUR_C_NAMESPACE_CLOSE
