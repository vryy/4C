/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of material definitions

\level 0


*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_mat_materialdefinition.hpp"

#include "4C_io_line_parser.hpp"
#include "4C_mat_par_bundle.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN



/*======================================================================*/
/*======================================================================*/
MAT::MaterialDefinition::MaterialDefinition(
    std::string materialname, std::string description, CORE::Materials::MaterialType mattype)
    : materialname_(std::move(materialname)),
      description_(std::move(description)),
      mattype_(mattype)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::MaterialDefinition::AddComponent(const Teuchos::RCP<INPUT::LineComponent>& c)
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<std::pair<int, CORE::IO::InputParameterContainer>> MAT::MaterialDefinition::Read(
    CORE::IO::DatFileReader& reader)
{
  std::string name = "--MATERIALS";
  std::vector<const char*> section = reader.Section(name);

  std::vector<std::pair<int, CORE::IO::InputParameterContainer>> found_materials;
  if (!section.empty())
  {
    for (std::vector<const char*>::const_iterator i = section.begin(); i != section.end(); ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      // add trailing white space to stringstream "condline" to avoid deletion of stringstream upon
      // reading the last entry inside This is required since the material parameters can be
      // specified in an arbitrary order in the input file. So it might happen that the last entry
      // is extracted before all of the previous ones are.
      condline->seekp(0, condline->end);
      *condline << " ";

      CORE::IO::LineParser parser("While reading 'MATERIALS' section: ");

      parser.Consume(*condline, "MAT");
      const int matid = parser.Read<int>(*condline);
      const std::string name = parser.Read<std::string>(*condline);

      // Remove the parts that were already read.
      condline->str(condline->str().erase(0, (size_t)condline->tellg()));

      if (name == materialname_)
      {
        if (matid <= -1) FOUR_C_THROW("Illegal negative ID provided");

        CORE::IO::InputParameterContainer input_data;
        for (auto& j : inputline_) condline = j->Read(Name(), condline, input_data);

        // current material input line contains bad elements
        if (condline->str().find_first_not_of(' ') != std::string::npos)
          FOUR_C_THROW(
              "Specification of material '%s' contains the following unknown, redundant, or "
              "incorrect elements: '%s'",
              materialname_.c_str(), condline->str().c_str());

        found_materials.emplace_back(matid, std::move(input_data));
      }
    }
  }

  return found_materials;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& MAT::MaterialDefinition::Print(std::ostream& stream, const DRT::Discretization* dis)
{
  // a string holding the comment indicating symbols for DAT input file
  const std::string comment = "//";

  // the descriptive lines (comments)
  stream << comment << std::endl;
  stream << comment << " " << description_ << std::endl;
  for (auto& i : inputline_)
  {
    std::ostringstream desc;
    i->Describe(desc);
    if (!desc.str().empty()) stream << comment << desc.str() << std::endl;
  }

  // the default line
  stream << comment << "MAT 0   " << materialname_ << "   ";
  for (auto& i : inputline_)
  {
    i->DefaultLine(stream);
    stream << " ";
  }

  stream << "\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::AppendMaterialDefinition(std::vector<Teuchos::RCP<MaterialDefinition>>& matlist,
    const Teuchos::RCP<MaterialDefinition>& mat)
{
  // test if material was defined with same name or type
  std::vector<Teuchos::RCP<MAT::MaterialDefinition>>::const_iterator m;
  for (m = matlist.begin(); m != matlist.end(); ++m)
  {
    Teuchos::RCP<MAT::MaterialDefinition> mmd = *m;

    if (mmd->Type() == mat->Type())
      FOUR_C_THROW(
          "Trying to define two materials with the same type '%d'\n"
          "Please revise your definitions of valid materials",
          mmd->Type());

    if (mmd->Name() == mat->Name())
      FOUR_C_THROW(
          "Trying to define two materials with the same name '%s'\n"
          "Please revise your definitions of valid materials",
          mmd->Name().c_str());
  }

  // no coincidence found
  if (m == matlist.end())
    matlist.push_back(mat);
  else
    FOUR_C_THROW("Trouble in determining coincidences of material definitions");
}

FOUR_C_NAMESPACE_CLOSE
