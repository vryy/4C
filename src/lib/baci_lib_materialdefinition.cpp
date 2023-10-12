/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of material definitions

\level 0


*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_lib_materialdefinition.H"

#include "baci_lib_globalproblem.H"
#include "baci_mat_material.H"

#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>



/*======================================================================*/
/*======================================================================*/
DRT::INPUT::MaterialDefinition::MaterialDefinition(
    std::string materialname, std::string description, INPAR::MAT::MaterialType mattype)
    : materialname_(std::move(materialname)),
      description_(std::move(description)),
      mattype_(mattype)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::MaterialDefinition::AddComponent(const Teuchos::RCP<::INPUT::LineComponent>& c)
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::MaterialDefinition::Read(
    const Problem& problem, DatFileReader& reader, const Teuchos::RCP<MAT::PAR::Bundle>& mmap)
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

      // read header from stringstream and delete it afterwards
      std::string mat, number, name;
      *condline >> mat >> number >> name;
      condline->str(condline->str().erase(0, (size_t)condline->tellg()));

      if (not(*condline) or mat != "MAT") dserror("invalid material line in '%s'", name.c_str());

      if (name == materialname_)
      {
        int matid = -1;
        {
          char* ptr;
          matid = std::strtol(number.c_str(), &ptr, 10);
          if (ptr == number.c_str())
            dserror("failed to read material object number '%s'", number.c_str());
        }

        if (matid <= -1)
          dserror("Either- failed to convert material ID -or- illegal negative ID provided");

        // what was read
        // std::cout << "PE=" << reader.Comm()->MyPID()
        //         << " MAT=" << matid
        //          << " MaterialType=" << mattype_
        //          << " Name=" << materialname_
        //          << std::endl;

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
std::ostream& DRT::INPUT::MaterialDefinition::Print(std::ostream& stream, const Discretization* dis)
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
void DRT::INPUT::AppendMaterialDefinition(
    std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& matlist,
    const Teuchos::RCP<DRT::INPUT::MaterialDefinition>& mat)
{
  // test if material was defined with same name or type
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>::const_iterator m;
  for (m = matlist.begin(); m != matlist.end(); ++m)
  {
    Teuchos::RCP<DRT::INPUT::MaterialDefinition> mmd = *m;

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
