/*----------------------------------------------------------------------*/
/*!
\file drt_materialdefinition.cpp

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */


/*----------------------------------------------------------------------*/
/* headers */
#include <algorithm>
#include <iterator>

#include "drt_materialdefinition.H"
#include "drt_colors.H"
#include "drt_globalproblem.H"

#include "../drt_mat/material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::MaterialComponent::MaterialComponent(
  std::string name,
  bool optional
  )
: optional_(optional),
  name_(name)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::MaterialComponent::PushBack(
  std::string token,
  Teuchos::RCP<std::stringstream> stream
  )
{
  Teuchos::RCP<std::stringstream> out = Teuchos::rcp(new std::stringstream());
  (*out) << token << " ";
  std::copy(std::istream_iterator<std::string>(*stream),
            std::istream_iterator<std::string>(),
            std::ostream_iterator<std::string>(*out," "));
  return out;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::StringMaterialComponent::StringMaterialComponent(
  std::string name,
  std::string defaultvalue,
  bool optional
  )
: MaterialComponent(name,optional),
  defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  stream << *cond->Get<std::string>(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringMaterialComponent::Describe(
  std::ostream& stream
  )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::StringMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  std::string value;
  (*condline) >> value;

  if ( (value=="") or (optional_) )
    value = defaultvalue_;

  material->Add(Name(),value);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::SeparatorMaterialComponent::SeparatorMaterialComponent(
  std::string separator,
  std::string description,
  bool optional
  )
: MaterialComponent("*SEPARATOR*",optional),
  separator_(separator),
  description_(description)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::SeparatorMaterialComponent::SeparatorMaterialComponent(
  std::string separator,
  bool optional
  )
: MaterialComponent("*SEPARATOR*",optional),
  separator_(separator),
  description_("")
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  stream << separator_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorMaterialComponent::Describe(
  std::ostream& stream
  )
{
  stream << "    "
         << std::setw(15) << std::left << separator_
         << std::setw(15) << std::left << (optional_ ? "(optional)" : "")
         << description_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::SeparatorMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  std::string sep;
  (*condline) >> sep;
  if ( (sep == "") and (optional_) )
    ;
  else if (sep != separator_)
    dserror("word '%s' expected but found '%s' while reading '%s'",
            separator_.c_str(),sep.c_str(),def->Name().c_str());
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntMaterialComponent::IntMaterialComponent(
  std::string name,
  const int defaultvalue,
  bool optional
  )
: MaterialComponent(name,optional),
  defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  stream << cond->GetInt(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntMaterialComponent::Describe(
  std::ostream& stream
  )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  std::string number;
  (*condline) >> number;

  int n;
  if ( (number=="") and (optional_) )
  {
    n = defaultvalue_;
  }
  else
  {
    char* ptr;
    n = strtol(number.c_str(),&ptr,10);
    if (ptr==number.c_str())
      dserror("failed to read number '%s' while reading variable '%s' in '%s'",
              number.c_str(),Name().c_str(),def->Name().c_str());
  }
  material->Add(Name(),n);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorMaterialComponent::IntVectorMaterialComponent(
  std::string name,
  int length,
  const int defaultvalue,
  bool optional
)
: MaterialComponent(name,optional),
  length_(length),
  lengthname_("*UNDEFINED*"),
  defaultvalue_(defaultvalue)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorMaterialComponent::IntVectorMaterialComponent(
  std::string name,
  std::string lengthname,
  const int defaultvalue,
  bool optional
)
: MaterialComponent(name,optional),
  length_(-1),
  lengthname_(lengthname),
  defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  for (int i=0; i<length_; ++i)
    stream << defaultvalue_ << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  const std::vector<int>* v = cond->Get<std::vector<int> >(Name());
  for (unsigned i=0; i<v->size(); ++i)
  {
//    stream << (*v)[i]+1 << " ";  // ??? : this is used in DRT::INPUT::IntVectorConditionComponent::Print
    stream << (*v)[i] << " ";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorMaterialComponent::Describe(
  std::ostream& stream
  )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntVectorMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  if (lengthname_ != "*UNDEFINED*")
    length_ = material->GetInt(lengthname_);
  else
    dserror("Trouble to get length of int vector material component.");

  std::vector<int> numbers(length_,defaultvalue_);

  for (int i=0; i<length_; ++i)
  {
    std::string number;
    (*condline) >> number;

    int n;
    {
      char* ptr;
      n = strtol(number.c_str(),&ptr,10);
      if (ptr==number.c_str())
      {
        if (optional_ and i==0)
        {
          // failed to read the numbers, fall back to default values
          condline = PushBack(number,condline);
          break;
        }
        dserror("failed to read number '%s' while reading variable '%s' in '%s'",
                number.c_str(),Name().c_str(),def->Name().c_str());
      }
    }
    numbers[i] = n;
  }
  material->Add(Name(),numbers);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealMaterialComponent::RealMaterialComponent(
  std::string name,
  const double defaultvalue,
  bool optional
  )
  : MaterialComponent(name,optional),
    defaultvalue_(defaultvalue)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  stream << cond->GetDouble(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealMaterialComponent::Describe(
  std::ostream& stream
  )
{
}




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  double number = 0;
  (*condline) >> number;
  material->Add(Name(),number);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorMaterialComponent::RealVectorMaterialComponent(
  std::string name,
  int length,
  const double defaultvalue,
  bool optional
)
: MaterialComponent(name,optional),
  length_(length),
  lengthname_("*UNDEFINED*"),
  defaultvalue_(defaultvalue)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorMaterialComponent::RealVectorMaterialComponent(
  std::string name,
  std::string lengthname,
  const double defaultvalue,
  bool optional
)
: MaterialComponent(name,optional),
  length_(-1),
  lengthname_(lengthname),
  defaultvalue_(defaultvalue)
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  for (int i=0; i<length_; ++i)
    stream << defaultvalue_ << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  const std::vector<double>* v = cond->Get<std::vector<double> >(Name());
  for (unsigned i=0; i<v->size(); ++i)
    stream << (*v)[i] << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorMaterialComponent::Describe(
  std::ostream& stream
  )
{
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealVectorMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  if (lengthname_ != "*UNDEFINED*")
    length_ = material->GetInt(lengthname_);
  else
    dserror("Trouble to get length of real vector material component.");

  std::vector<double> numbers(length_,defaultvalue_);

  for (int i=0; i<length_; ++i)
  {
    double number = 0;
    (*condline) >> number;
    numbers[i] = number;
  }
  material->Add(Name(),numbers);
  return condline;
}


/*======================================================================*/
/*======================================================================*/
const std::string DRT::INPUT::BoolMaterialComponent::lineTrue_ = "Yes";
const std::string DRT::INPUT::BoolMaterialComponent::lineFalse_ = "No";
DRT::INPUT::BoolMaterialComponent::BoolMaterialComponent(
  std::string name,
  const bool defaultvalue,
  bool optional
  )
: MaterialComponent(name,optional),
  defaultvalue_(defaultvalue)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::DefaultLine(
  std::ostream& stream
  )
{
  PrintYesNo(stream,defaultvalue_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::Print(
  std::ostream& stream,
  const MAT::PAR::Material* cond
  )
{
  const bool value = (bool)cond->GetInt(Name());
  PrintYesNo(stream,value);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::PrintYesNo(
  std::ostream& stream,
  const bool value
  ) const
{
  if (value)
    stream << lineTrue_;
  else
    stream << lineFalse_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::BoolMaterialComponent::Describe(
  std::ostream& stream
  )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::BoolMaterialComponent::Read(
  DRT::INPUT::MaterialDefinition* def,
  Teuchos::RCP<std::stringstream> condline,
  Teuchos::RCP<MAT::PAR::Material> material
  )
{
  std::string sval;
  (*condline) >> sval;

  int ival = -1;

  // added the treatment of the case, if the paramter is optional and therefore not existent in condline
  // Furthermore expanded the existing case to allow the input Yes, YES, yes, TRUE, true, True etc ...
  //																																														nagler 2012

  // Setting ival for optional input
  if (("" == sval) and (optional_))
	{
  	if (false == defaultvalue_)
  	{
  		ival = 0;;
  	}
  	else if (true == defaultvalue_)
  	{
  		ival = 1;
  	}
  	else
  	{
  		dserror("Failed to interpret given default value '%s'", defaultvalue_);
  	}
	}

  // setting ival according to given value
  else if (("Yes" == sval) or ("YES" == sval) or ("yes" == sval) or ("True" == sval) or ("TRUE" == sval) or ("true" == sval))
    ival = 1;
  else if (("No" == sval) or ("NO" == sval) or ("no" == sval) or ("False" == sval) or ("FALSE" == sval) or ("false" == sval))
    ival = 0;
  else
    dserror("Failed to read Boolean value '%s' while reading Boolean variable '%s' in '%s'",
            sval.c_str(),Name().c_str(),def->Name().c_str());

  material->Add(Name(),ival);
  return condline;
}


/*======================================================================*/
/*======================================================================*/
DRT::INPUT::MaterialDefinition::MaterialDefinition(
  std::string materialname,
  std::string description,
  INPAR::MAT::MaterialType mattype
  )
: materialname_(materialname),
  description_(description),
  mattype_(mattype)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::MaterialDefinition::AddComponent(
  Teuchos::RCP<MaterialComponent> c
  )
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::MaterialDefinition::Read(
  const Problem& problem,
  DatFileReader& reader,
  Teuchos::RCP<MAT::PAR::Bundle> mmap
  )
{
  std::string name = "--MATERIALS";
  std::vector<const char*> section = reader.Section(name);

  if (section.size() > 0)
  {

    for (std::vector<const char*>::const_iterator i=section.begin();
         i!=section.end();
         ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      std::string mat;
      std::string number;
      std::string name;
      (*condline) >> mat >> number >> name;
      if (not (*condline) or mat!="MAT")
        dserror("invalid material line in '%s'",name.c_str());

      if (name == materialname_)
      {
        int matid = -1;
        {
          char* ptr;
          matid = std::strtol(number.c_str(),&ptr,10);
          if (ptr == number.c_str())
            dserror("failed to read material object number '%s'",
                    number.c_str());
        }

        if (matid <= -1)
          dserror("Either- failed to convert material ID -or- illegal negative ID provided");

        // what was read
        //std::cout << "PE=" << reader.Comm()->MyPID()
        //         << " MAT=" << matid
        //          << " MaterialType=" << mattype_
        //          << " Name=" << materialname_
        //          << std::endl;

        // check if material ID is already in use
        if (mmap->Find(matid) != -1)
          dserror("More than one material with 'MAT %d'", matid);

        // the read-in material line
        Teuchos::RCP<MAT::PAR::Material> material = Teuchos::rcp(new MAT::PAR::Material(matid,mattype_,materialname_));
        // fill the latter

        for (unsigned j=0; j<inputline_.size(); ++j)
        {
          condline = inputline_[j]->Read(this,condline,material);
        }

        // put material in map of materials
        mmap->Insert(matid,material);

      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& DRT::INPUT::MaterialDefinition::Print(
  std::ostream& stream,
  const Discretization* dis,
  const bool color
  )
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
  for (unsigned i=0; i<inputline_.size(); ++i)
  {
    std::ostringstream desc;
    inputline_[i]->Describe(desc);
    if (desc.str().size() > 0)
      stream << comment << desc.str() << std::endl;
  }

  // the default line
  stream << blue2light << comment << "MAT 0   " << magentalight << materialname_ << "   ";
  for (unsigned i=0; i<inputline_.size(); ++i)
  {
    inputline_[i]->DefaultLine(stream);
    stream << " ";
  }

  stream << endcolor << "\n";

  return stream;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::AppendMaterialDefinition(
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >& matlist,
  Teuchos::RCP<DRT::INPUT::MaterialDefinition> mat
  )
{
  // test if material was defined with same name or type
  std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition> >::const_iterator m;
  for (m=matlist.begin(); m!=matlist.end(); ++m)
  {
    Teuchos::RCP<DRT::INPUT::MaterialDefinition> mmd = *m;

    if (mmd->Type() == mat->Type())
      dserror("Trying to define two materials with the same type '%d'\n"
              "Please revise your definitions of valid materials", mmd->Type());

    if (mmd->Name() == mat->Name())
      dserror("Trying to define two materials with the same name '%s'\n"
              "Please revise your definitions of valid materials", mmd->Name().c_str());
  }

  // no coincidence found
  if (m==matlist.end())
    matlist.push_back(mat);
  else
    dserror("Trouble in determining coincidences of material definitions");
}

