/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of contact constitutive law definitions

\level 3

\maintainer Nora Hagmeyer

*/
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* headers */

#include "constitutivelaw_definition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_inputreader.H"
#include "contact_constitutivelaw_bundle.H"
#include "../drt_lib/drt_colors.H"
#include "coconstlaw_parameter.H"

/*======================================================================*/
/*======================================================================*/
CONTACT::CONSTITUTIVELAW::LawDefinition::LawDefinition(
    std::string name, std::string description, INPAR::CONTACT::ConstitutiveLawType type)
    : coconstlawname_(name), description_(description), coconstlawtype_(type)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LawDefinition::AddComponent(
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::CoConstLawComponent> c)
{
  inputline_.push_back(c);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::LawDefinition::Read(const DRT::Problem& problem,
    DRT::INPUT::DatFileReader& reader, Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Bundle> bundle)
{
  std::string name = "--CONTACT CONSTITUTIVE LAWS";
  std::vector<const char*> section = reader.Section(name);

  if (section.size() > 0)
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
      std::string law, number, name;
      *condline >> law >> number >> name;
      condline->str(condline->str().erase(0, (size_t)condline->tellg()));

      if (not(*condline) or law != "LAW")
        dserror("invalid contact constitutive line in '%s'", name.c_str());

      if (name == coconstlawname_)
      {
        int id = -1;
        {
          char* ptr;
          id = std::strtol(number.c_str(), &ptr, 10);
          if (ptr == number.c_str())
            dserror("failed to read contact constitutive law object number '%s'", number.c_str());
        }

        if (id <= -1)
          dserror(
              "Either- failed to convert contact constitutive law ID -or- illegal negative ID "
              "provided");

        // check if material ID is already in use
        if (bundle->Find(id) != -1)
          dserror("More than one contact constitutivelaw with 'Law %d'", id);

        // the read-in contact constitutive law line
        Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container = Teuchos::rcp(
            new CONTACT::CONSTITUTIVELAW::Container(id, coconstlawtype_, coconstlawname_));
        // fill the latter

        for (unsigned j = 0; j < inputline_.size(); ++j)
          condline = inputline_[j]->Read(this, condline, container);

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

  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::CoConstLawComponent::CoConstLawComponent(std::string name, bool optional)
    : optional_(optional), name_(name)
{
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::StringCoConstLawComponent::StringCoConstLawComponent(
    std::string name, std::string defaultvalue, bool optional)
    : CoConstLawComponent(name, optional), defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::StringCoConstLawComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::StringCoConstLawComponent::Print(
    std::ostream& stream, const CONTACT::CONSTITUTIVELAW::Container* cond)
{
  stream << *cond->Get<std::string>(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::StringCoConstLawComponent::Describe(std::ostream& stream) {}

/*---------------------------------------------------------------------------*/

Teuchos::RCP<std::stringstream> CONTACT::CONSTITUTIVELAW::StringCoConstLawComponent::Read(
    CONTACT::CONSTITUTIVELAW::LawDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container)
{
  // initialize string parameter value to be read
  std::string str = defaultvalue_;

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read string parameter value in case the associated parameter label appears in
  // line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract string parameter value from stringstream "condline"
    *condline >> str;

    // return error in case the extraction was not successful
    if (str.empty())
      dserror(
          "Value of parameter '%s' for contact constitutive law '%s' not properly specified in "
          "input file!",
          Name().c_str(), def->Name().c_str());

    // remove string parameter value from stringstream "condline"
    condline->str(condline->str().erase((size_t)condline->tellg() - str.size(), str.size()));

    // reset current position in stringstream "condline"
    condline->seekg(position);
  }

  // add double parameter value to material parameter list
  container->Add(Name(), str);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::SeparatorCoConstLawComponent::SeparatorCoConstLawComponent(
    std::string separator, std::string description, bool optional)
    : CoConstLawComponent("*SEPARATOR*", optional), separator_(separator), description_(description)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::SeparatorCoConstLawComponent::DefaultLine(std::ostream& stream)
{
  stream << separator_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::SeparatorCoConstLawComponent::Print(
    std::ostream& stream, const CONTACT::CONSTITUTIVELAW::Container* cond)
{
  stream << separator_;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::SeparatorCoConstLawComponent::Describe(std::ostream& stream)
{
  stream << "    " << std::setw(15) << std::left << separator_ << std::setw(15) << std::left
         << (optional_ ? "(optional)" : "") << description_;
}


/*-----------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> CONTACT::CONSTITUTIVELAW::SeparatorCoConstLawComponent::Read(
    CONTACT::CONSTITUTIVELAW::LawDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container)
{
  // try to find parameter label "separator_" (with leading and trailing white spaces for
  // uniqueness) in stringstream "condline"
  size_t position = condline->str().find(" " + separator_ + " ");

  // case: material parameter label "separator_" not found
  if (position == std::string::npos)
  {
    if (optional_)
      // move stringstream position to end of "condline" in case of optional parameter
      condline->seekg(0, condline->end);
    else
      // return error in case a required parameter is not specified
      dserror(
          "Required parameter '%s' for contact constitutive law '%s' not specified in input file!",
          separator_.c_str(), def->Name().c_str());
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
CONTACT::CONSTITUTIVELAW::RealCoConstLawComponent::RealCoConstLawComponent(
    std::string name, const double defaultvalue, bool optional)
    : CoConstLawComponent(name, optional), defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::RealCoConstLawComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::RealCoConstLawComponent::Print(
    std::ostream& stream, const CONTACT::CONSTITUTIVELAW::Container* cond)
{
  stream << cond->GetDouble(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::RealCoConstLawComponent::Describe(std::ostream& stream) {}


/*---------------------------------------------------------------------------*
 | Read double parameter value from material line of input file   fang 08/14 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> CONTACT::CONSTITUTIVELAW::RealCoConstLawComponent::Read(
    CONTACT::CONSTITUTIVELAW::LawDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container)
{
  // initialize double parameter value to be read
  double number = defaultvalue_;

  // get current position in stringstream "condline"
  std::streampos position = condline->tellg();

  // only try to read double parameter value in case the associated parameter label appears in
  // material line of input file
  if ((size_t)position != condline->str().size())
  {
    // extract double parameter value from stringstream "condline" as string
    std::string snumber;
    *condline >> snumber;

    // try to convert to double
    char* check;
    number = strtod(snumber.c_str(), &check);

    // return error in case the conversion was not successful
    if (snumber == check)
      dserror(
          "Value of parameter '%s' for contact constitutive law '%s' not properly specified in "
          "input file!",
          Name().c_str(), def->Name().c_str());

    // remove double parameter value from stringstream "condline"
    condline->str(
        condline->str().erase((size_t)condline->tellg() - snumber.size(), snumber.size()));

    // reset current position in stringstream "condline"
    condline->seekg(position);
  }

  // add double parameter value to contact constitutive law parameter list
  container->Add(Name(), number);

  return condline;
}
/*======================================================================*/
/*======================================================================*/

std::ostream& CONTACT::CONSTITUTIVELAW::LawDefinition::Print(
    std::ostream& stream, const DRT::Discretization* dis, const bool color)
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


  // a string holding the comment indicating symbols for DAT input file
  const std::string comment = "//";

  // the descriptive lines (comments)
  stream << blue2light << comment << std::endl;
  stream << blue2light << comment << " " << magentalight << description_ << std::endl;
  for (unsigned i = 0; i < inputline_.size(); ++i)
  {
    std::ostringstream desc;
    inputline_[i]->Describe(desc);
    if (desc.str().size() > 0) stream << comment << desc.str() << std::endl;
  }

  // the default line
  stream << blue2light << comment << "LAW 0   " << magentalight << coconstlawname_ << "   ";
  for (unsigned i = 0; i < inputline_.size(); ++i)
  {
    inputline_[i]->DefaultLine(stream);
    stream << " ";
  }

  stream << endcolor << "\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::AppendCoConstLawDefinition(
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
