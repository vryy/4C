// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_constitutivelaw_constitutivelaw_definition.hpp"

#include "4C_contact_constitutivelaw_bundle.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_io_inputreader.hpp"
#include "4C_io_value_parser.hpp"
#include "4C_utils_exceptions.hpp"


FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/*======================================================================*/
CONTACT::CONSTITUTIVELAW::LawDefinition::LawDefinition(
    std::string name, std::string description, Inpar::CONTACT::ConstitutiveLawType type)
    : coconstlawname_(name), description_(description), coconstlawtype_(type)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LawDefinition::add_component(Teuchos::RCP<Input::LineComponent> c)
{
  inputline_.push_back(c);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LawDefinition::read(const Global::Problem& problem,
    Core::IO::DatFileReader& reader, CONTACT::CONSTITUTIVELAW::Bundle& bundle)
{
  std::string name = "--CONTACT CONSTITUTIVE LAWS";
  std::vector<const char*> section = reader.section(name);

  if (section.size() > 0)
  {
    for (auto& i : section)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::make_rcp<std::stringstream>(i);

      // add trailing white space to stringstream "condline" to avoid deletion of stringstream upon
      // reading the last entry inside This is required since the material parameters can be
      // specified in an arbitrary order in the input file. So it might happen that the last entry
      // is extracted before all of the previous ones are.
      condline->seekp(0, condline->end);
      *condline << " ";

      Core::IO::ValueParser parser(
          *condline, "While reading 'CONTACT CONSTITUTIVE LAWS' section: ");

      parser.consume("LAW");
      const int id = parser.read<int>();
      const std::string name = parser.read<std::string>();

      // Remove the parts that were already read.
      condline->str(condline->str().erase(0, (size_t)condline->tellg()));

      if (name == coconstlawname_)
      {
        if (id <= -1) FOUR_C_THROW("Illegal negative ID provided");

        // check if material ID is already in use
        if (bundle.find(id) != -1)
          FOUR_C_THROW("More than one contact constitutivelaw with 'Law %d'", id);

        // the read-in contact constitutive law line
        Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container =
            Teuchos::make_rcp<CONTACT::CONSTITUTIVELAW::Container>(
                id, coconstlawtype_, coconstlawname_);
        // fill the latter

        for (const auto& j : inputline_)
          condline = j->read(LawDefinition::name(), condline, *container);

        // current material input line contains bad elements
        if (condline->str().find_first_not_of(' ') != std::string::npos)
          FOUR_C_THROW(
              "Specification of material '%s' contains the following unknown, redundant, or "
              "incorrect elements: '%s'",
              coconstlawname_.c_str(), condline->str().c_str());

        // put contact constitutive law in map of laws
        bundle.insert(id, container);
      }
    }
  }
}


std::ostream& CONTACT::CONSTITUTIVELAW::LawDefinition::print(
    std::ostream& stream, const Core::FE::Discretization* dis)
{
  // a string holding the comment indicating symbols for DAT input file
  const std::string comment = "//";

  // the descriptive lines (comments)
  stream << comment << std::endl;
  stream << comment << " " << description_ << std::endl;
  for (const auto& i : inputline_)
  {
    std::ostringstream desc;
    i->describe(desc);
    if (desc.str().size() > 0) stream << comment << desc.str() << std::endl;
  }

  // the default line
  stream << comment << "LAW 0   " << coconstlawname_ << "   ";
  for (const auto& i : inputline_)
  {
    i->default_line(stream);
    stream << " ";
  }

  stream << "\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::append_co_const_law_component_definition(
    std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& list,
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> def)
{
  // test if material was defined with same name or type
  std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>::const_iterator m;
  for (m = list.begin(); m != list.end(); ++m)
  {
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> mmd = *m;

    if (mmd->type() == def->type())
      FOUR_C_THROW(
          "Trying to define two contact constitutive laws with the same type '%d'\n"
          "Please revise your definitions of valid contact constitutive laws",
          mmd->type());

    if (mmd->name() == def->name())
      FOUR_C_THROW(
          "Trying to define two contact constitutive laws with the same name '%s'\n"
          "Please revise your definitions of valid contact constitutive laws",
          mmd->name().c_str());
  }

  // no coincidence found
  if (m == list.end())
    list.push_back(def);
  else
    FOUR_C_THROW("Trouble in determining coincidences of contact constitutive definitions");
}

FOUR_C_NAMESPACE_CLOSE
