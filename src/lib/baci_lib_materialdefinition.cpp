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
#include "baci_lib_utils_cond_and_mat_definition.H"
#include "baci_mat_material.H"

#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::StringMaterialComponent::StringMaterialComponent(
    std::string name, std::string defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional), defaultvalue_(std::move(defaultvalue))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringMaterialComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringMaterialComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  stream << *cond.Get<std::string>(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringMaterialComponent::Describe(std::ostream& stream) {}


/*---------------------------------------------------------------------------*
 | Read string parameter value from material line of input file   fang 08/14 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::StringMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  // initialize string parameter value to be read
  std::string str = defaultvalue_;

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read string parameter value in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract string parameter value from stringstream "condline"
    *condline >> str;

    // return error in case the extraction was not successful
    if (str.empty())
      dserror("Value of parameter '%s' for material '%s' not properly specified in input file!",
          Name().c_str(), section_name.c_str());

    // remove string parameter value from stringstream "condline"
    condline->str(condline->str().erase((size_t)condline->tellg() - str.size(), str.size()));

    // reset current position in stringstream "condline"
    condline->seekg(position);
  }

  // add double parameter value to material parameter list
  material.Add(Name(), str);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::SeparatorMaterialComponent::SeparatorMaterialComponent(
    std::string separator, std::string description, bool optional)
    : ::INPUT::LineComponent("*SEPARATOR*", optional),
      separator_(std::move(separator)),
      description_(std::move(description))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorMaterialComponent::DefaultLine(std::ostream& stream)
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorMaterialComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  stream << separator_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorMaterialComponent::Describe(std::ostream& stream)
{
  stream << "    " << std::setw(15) << std::left << separator_ << std::setw(15) << std::left
         << (optional_ ? "(optional)" : "") << description_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<std::string> DRT::INPUT::SeparatorMaterialComponent::WriteReadTheDocs() const
{
  std::vector<std::string> tablerow;

  tablerow.push_back(separator_);
  tablerow.push_back((optional_ ? "yes" : ""));
  std::string descriptionstring = "";
  tablerow.push_back(description_);
  return tablerow;
}



/*-----------------------------------------------------------------------------*
 | Find material parameter label in material line of input file     fang 08/14 |
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::SeparatorMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  // try to find material parameter label "separator_" (with leading and trailing white spaces for
  // uniqueness) in stringstream "condline"
  size_t position = condline->str().find(" " + separator_ + " ");

  // case: material parameter label "separator_" not found
  if (position == std::string::npos)
  {
    if (optional_)
      // move stringstream position to end of "condline" in case of optional material parameter
      condline->seekg(0, condline->end);
    else
      // return error in case a required material parameter is not specified
      dserror("Required parameter '%s' for material '%s' not specified in input file!",
          separator_.c_str(), section_name.c_str());
  }
  // case: found material parameter label "separator_"
  else
  {
    // care for leading white space in search string ("position" should indicate position of first
    // actual character of parameter label)
    position++;

    // remove material parameter label "separator_" from stringstream "condline"
    condline->str(condline->str().erase(position, separator_.size()));

    // set current position in stringstream "condline" in front of value associated with material
    // parameter label "separator_"
    condline->seekg((std::streampos)position);
  }

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntMaterialComponent::IntMaterialComponent(
    std::string name, const int defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional), defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntMaterialComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntMaterialComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  stream << cond.GetInt(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntMaterialComponent::Describe(std::ostream& stream) {}


/*----------------------------------------------------------------------------*
 | Read integer parameter value from material line of input file   fang 08/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  // initialize integer parameter value to be read
  int nnumber = defaultvalue_;

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read integer parameter value in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract integer vector component as string
    std::string snumber;
    *condline >> snumber;

    nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<int>(
        snumber, Name(), section_name, 1, optional_);

    // remove parameter value from stringstream "condline"
    condline->str(
        condline->str().erase((size_t)condline->tellg() - snumber.size(), snumber.size()));

    // reset current position in stringstream "condline"
    condline->seekg(position);
  }

  // add int parameter value to material parameter list
  material.Add(Name(), nnumber);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorMaterialComponent::IntVectorMaterialComponent(
    std::string name, int length, const int defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional),
      length_(length),
      lengthname_("*UNDEFINED*"),
      defaultvalue_(defaultvalue)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorMaterialComponent::IntVectorMaterialComponent(
    std::string name, std::string lengthname, const int defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional),
      length_(-1),
      lengthname_(std::move(lengthname)),
      defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorMaterialComponent::DefaultLine(std::ostream& stream)
{
  for (int i = 0; i < length_; ++i) stream << defaultvalue_ << " ";
  if (length_ == 0) stream << "<int vec> ;";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorMaterialComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  const auto* v = cond.Get<std::vector<int>>(Name());
  for (int i : *v)
  {
    //    stream << (*v)[i]+1 << " ";  // ??? : this is used in
    //    DRT::INPUT::IntVectorConditionComponent::Print
    stream << i << " ";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorMaterialComponent::Describe(std::ostream& stream) {}


/*-----------------------------------------------------------------------------*
 | Read integer parameter vector from material line of input file   fang 08/14 |
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntVectorMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  if (lengthname_ != "*UNDEFINED*") length_ = material.GetInt(lengthname_);
  if (length_ == -1) dserror("Trouble to get length of integer vector material component.");

  // initialize integer parameter vector to be read
  std::vector<int> nnumbers(length_, defaultvalue_);

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read integer parameter vector in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract integer parameter vector from stringstream "condline"
    for (auto& current_nnumber : nnumbers)
    {
      // extract integer vector component as string
      std::string snumber;
      *condline >> snumber;

      current_nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<int>(
          snumber, Name(), section_name, length_, optional_);

      // remove parameter value from stringstream "condline"
      condline->str(
          condline->str().erase((size_t)condline->tellg() - snumber.size(), snumber.size()));

      // reset current position in stringstream "condline"
      condline->seekg(position);
    }
  }

  // add int parameter vector to material parameter list
  material.Add(Name(), nnumbers);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealMaterialComponent::RealMaterialComponent(
    std::string name, const double defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional), defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealMaterialComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealMaterialComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  stream << cond.GetDouble(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealMaterialComponent::Describe(std::ostream& stream) {}


/*---------------------------------------------------------------------------*
 | Read double parameter value from material line of input file   fang 08/14 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  // initialize double parameter value to be read
  double nnumber = defaultvalue_;

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read double parameter value in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // read string from stream (need to handle doubles and "none" arguments)
    std::string snumber;
    *condline >> snumber;

    nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<double>(
        snumber, Name(), section_name, 1, optional_);

    // remove parameter value from stringstream "condline"
    condline->str(
        condline->str().erase((size_t)condline->tellg() - snumber.size(), snumber.size()));

    // reset current position in stringstream "condline"
    condline->seekg(position);
  }

  // add double parameter value to material parameter list
  material.Add(Name(), nnumber);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorMaterialComponent::RealVectorMaterialComponent(
    std::string name, int length, const double defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional),
      length_(length),
      lengthname_("*UNDEFINED*"),
      defaultvalue_(defaultvalue)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorMaterialComponent::RealVectorMaterialComponent(
    std::string name, std::string lengthname, const double defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional),
      length_(-1),
      lengthname_(std::move(lengthname)),
      defaultvalue_(defaultvalue)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorMaterialComponent::DefaultLine(std::ostream& stream)
{
  for (int i = 0; i < length_; ++i) stream << defaultvalue_ << " ";
  if (length_ == 0) stream << "<real vec> ;";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorMaterialComponent::Print(
    std::ostream& stream, const DRT::Container& cond)
{
  const auto* v = cond.Get<std::vector<double>>(Name());
  for (double i : *v) stream << i << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorMaterialComponent::Describe(std::ostream& stream) {}


/*----------------------------------------------------------------------------*
 | Read double parameter vector from material line of input file   fang 08/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealVectorMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  if (lengthname_ != "*UNDEFINED*")
    length_ = material.GetInt(lengthname_);
  else if (length_ < 0)
    dserror("Trouble to get length of real vector material component.");

  // initialize double parameter vector to be read
  std::vector<double> nnumbers(length_, defaultvalue_);

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read double parameter vector in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract double parameter vector from stringstream "condline"
    for (auto& current_nnumber : nnumbers)
    {
      // extract double vector component as string
      std::string snumber;
      *condline >> snumber;

      current_nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<double>(
          snumber, Name(), section_name, length_, optional_);

      // remove parameter value from stringstream "condline"
      condline->str(
          condline->str().erase((size_t)condline->tellg() - snumber.size(), snumber.size()));

      // reset current position in stringstream "condline"
      condline->seekg(position);
    }
  }

  // add double parameter vector to material parameter list
  material.Add(Name(), nnumbers);

  return condline;
}


/*======================================================================*/
/*======================================================================*/
const std::string DRT::INPUT::BoolMaterialComponent::lineTrue_ = "Yes";
const std::string DRT::INPUT::BoolMaterialComponent::lineFalse_ = "No";
DRT::INPUT::BoolMaterialComponent::BoolMaterialComponent(
    std::string name, const bool defaultvalue, bool optional)
    : ::INPUT::LineComponent(std::move(name), optional), defaultvalue_(defaultvalue)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::DefaultLine(std::ostream& stream)
{
  PrintYesNo(stream, defaultvalue_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  const bool value = (bool)cond.GetInt(Name());
  PrintYesNo(stream, value);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::PrintYesNo(std::ostream& stream, const bool value) const
{
  if (value)
    stream << lineTrue_;
  else
    stream << lineFalse_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::Describe(std::ostream& stream) {}


/*----------------------------------------------------------------------------*
 | Read boolean parameter value from material line of input file   fang 08/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::BoolMaterialComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& material)
{
  // initialize boolean parameter value to be read
  bool boolean = defaultvalue_;

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read boolean parameter value in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract boolean parameter value from stringstream "condline" as string
    std::string sboolean;
    *condline >> sboolean;

    // try to convert to bool
    if (sboolean == "Yes" or sboolean == "YES" or sboolean == "yes" or sboolean == "True" or
        sboolean == "TRUE" or sboolean == "true")
      boolean = true;
    else if (sboolean == "No" or sboolean == "NO" or sboolean == "no" or sboolean == "False" or
             sboolean == "FALSE" or sboolean == "false")
      boolean = false;
    else
      // return error in case the conversion was not successful
      dserror("Value of parameter '%s' for material '%s' not properly specified in input file!",
          Name().c_str(), section_name.c_str());

    // remove boolean parameter value from stringstream "condline"
    condline->str(
        condline->str().erase((size_t)condline->tellg() - sboolean.size(), sboolean.size()));

    // reset current position in stringstream "condline"
    condline->seekg(position);
  }

  // add boolean parameter value to material parameter list
  material.Add(Name(), boolean);

  return condline;
}


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
