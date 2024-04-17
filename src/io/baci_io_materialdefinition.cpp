/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of material definitions

\level 0


*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_io_materialdefinition.hpp"

#include "baci_io_line_parser.hpp"
#include "baci_mat_material.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN



/*======================================================================*/
/*======================================================================*/
INPUT::MaterialDefinition::MaterialDefinition(
    std::string materialname, std::string description, INPAR::MAT::MaterialType mattype)
    : materialname_(std::move(materialname)),
      description_(std::move(description)),
      mattype_(mattype)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void INPUT::MaterialDefinition::AddComponent(const Teuchos::RCP<INPUT::LineComponent>& c)
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void INPUT::MaterialDefinition::Read(
    DatFileReader& reader, const Teuchos::RCP<MAT::PAR::Bundle>& mmap)
{
  std::string name = "--MATERIALS";
  std::vector<const char*> section = reader.Section(name);

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

      IO::LineParser parser("While reading 'MATERIALS' section: ");

      parser.Consume(*condline, "MAT");
      const int matid = parser.Read<int>(*condline);
      const std::string name = parser.Read<std::string>(*condline);

      // Remove the parts that were already read.
      condline->str(condline->str().erase(0, (size_t)condline->tellg()));

      if (name == materialname_)
      {
        if (matid <= -1) dserror("Illegal negative ID provided");

        // check if material ID is already in use
        if (mmap->Find(matid) != -1) dserror("More than one material with 'MAT %d'", matid);

        // the read-in material line
        Teuchos::RCP<MAT::PAR::Material> material =
            Teuchos::rcp(new MAT::PAR::Material(matid, mattype_, materialname_));
        // fill the latter

        for (auto& j : inputline_) condline = j->Read(Name(), condline, *material);

        // current material input line contains bad elements
        if (condline->str().find_first_not_of(' ') != std::string::npos)
          dserror(
              "Specification of material '%s' contains the following unknown, redundant, or "
              "incorrect elements: '%s'",
              materialname_.c_str(), condline->str().c_str());

        // put material in map of materials
        mmap->Insert(matid, material);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& INPUT::MaterialDefinition::Print(std::ostream& stream, const DRT::Discretization* dis)
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
void INPUT::AppendMaterialDefinition(std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>& matlist,
    const Teuchos::RCP<INPUT::MaterialDefinition>& mat)
{
  // test if material was defined with same name or type
  std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>::const_iterator m;
  for (m = matlist.begin(); m != matlist.end(); ++m)
  {
    Teuchos::RCP<INPUT::MaterialDefinition> mmd = *m;

    if (mmd->Type() == mat->Type())
      dserror(
          "Trying to define two materials with the same type '%d'\n"
          "Please revise your definitions of valid materials",
          mmd->Type());

    if (mmd->Name() == mat->Name())
      dserror(
          "Trying to define two materials with the same name '%s'\n"
          "Please revise your definitions of valid materials",
          mmd->Name().c_str());
  }

  // no coincidence found
  if (m == matlist.end())
    matlist.push_back(mat);
  else
    dserror("Trouble in determining coincidences of material definitions");
}

FOUR_C_NAMESPACE_CLOSE
