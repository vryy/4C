/*---------------------------------------------------------------------*/
/*!

\brief Implementation of base class for all conditions

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/


#include <algorithm>
#include <iterator>

#include "drt_conditiondefinition.H"
#include "drt_colors.H"
#include "drt_globalproblem.H"
#include "drt_discret.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::ConditionComponent::ConditionComponent(std::string name) : name_(name) {}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::ConditionComponent::PushBack(
    std::string token, Teuchos::RCP<std::stringstream> stream)
{
  Teuchos::RCP<std::stringstream> out = Teuchos::rcp(new std::stringstream());
  (*out) << token << " ";
  std::copy(std::istream_iterator<std::string>(*stream), std::istream_iterator<std::string>(),
      std::ostream_iterator<std::string>(*out, " "));
  return out;
}

/* -----------------------------------------------------------------------------------------------*
 | Class StringConditionComponent                                                       ehrl 09/12|
 * -----------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | StringConditionComponent::Constructor()                              |
 *----------------------------------------------------------------------*/
DRT::INPUT::StringConditionComponent::StringConditionComponent(std::string name,
    std::string defaultvalue, const Teuchos::Array<std::string>& datfilevalues,
    const Teuchos::Array<std::string>& stringcondvalues, bool optional)
    : ConditionComponent(name),
      defaultvalue_(defaultvalue),
      datfilevalues_(datfilevalues),
      stringcondvalues_(stringcondvalues),
      intcondvalues_(Teuchos::tuple<int>(-1)),
      stringtostring_(true)
{
  if (std::find(datfilevalues_.begin(), datfilevalues_.end(), defaultvalue_) ==
      datfilevalues_.end())
  {
    dserror("invalid default value '%s'", defaultvalue_.c_str());
  }
  if (datfilevalues_.size() != stringcondvalues_.size())
  {
    dserror("dat file values must match condition values");
  }
}


/*----------------------------------------------------------------------*
| StringConditionComponent::Constructor()                     ehrl 09/12|
*----------------------------------------------------------------------*/
DRT::INPUT::StringConditionComponent::StringConditionComponent(std::string name,
    std::string defaultvalue, const Teuchos::Array<std::string>& datfilevalues,
    const Teuchos::Array<int>& intcondvalues, bool optional)
    : ConditionComponent(name),
      defaultvalue_(defaultvalue),
      datfilevalues_(datfilevalues),
      stringcondvalues_(Teuchos::tuple<std::string>("notdefined")),
      intcondvalues_(intcondvalues),
      stringtostring_(false)
{
  if (std::find(datfilevalues_.begin(), datfilevalues_.end(), defaultvalue_) ==
      datfilevalues_.end())
  {
    dserror("invalid default value '%s'", defaultvalue_.c_str());
  }
  if (datfilevalues_.size() != intcondvalues_.size())
  {
    dserror("dat file values must match condition values");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringConditionComponent::DefaultLine(std::ostream& stream)
{
  stream << defaultvalue_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::StringConditionComponent::Print(std::ostream& stream, const Condition* cond)
{
  stream << *cond->Get<std::string>(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::StringConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  std::string value;
  (*condline) >> value;

  if (value == "") value = defaultvalue_;

  Teuchos::Array<std::string>::iterator i =
      std::find(datfilevalues_.begin(), datfilevalues_.end(), value);
  if (i == datfilevalues_.end())
  {
    dserror("unrecognized std::string '%s' while reading variable '%s' in '%s'", value.c_str(),
        Name().c_str(), def->SectionName().c_str());
  }
  unsigned pos = &*i - &datfilevalues_[0];
  // choose, if we have an array based on std::string or int
  if (stringtostring_)
    condition->Add(Name(), stringcondvalues_[pos]);
  else
    condition->Add(Name(), intcondvalues_[pos]);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::SeparatorConditionComponent::SeparatorConditionComponent(
    std::string separator, bool optional)
    : ConditionComponent("*SEPARATOR*"), separator_(separator), optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorConditionComponent::DefaultLine(std::ostream& stream)
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorConditionComponent::Print(
    std::ostream& stream, const DRT::Condition* cond)
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::SeparatorConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  std::string sep;
  (*condline) >> sep;

  if (sep != separator_)
  {
    if (sep == "" && optional_)
    {
      // not given, fall back to default line
      condline = PushBack(separator_, condline);
      // recursive call -- fix from 01/18 by hiermeier
      Read(def, condline, condition);
    }
    else
    {
      dserror("word '%s' expected but found '%s' while reading '%s'", separator_.c_str(),
          sep.c_str(), def->SectionName().c_str());
    }
  }
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntConditionComponent::IntConditionComponent(
    std::string name, bool fortranstyle, bool noneallowed)
    : ConditionComponent(name), fortranstyle_(fortranstyle), noneallowed_(noneallowed)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntConditionComponent::DefaultLine(std::ostream& stream)
{
  if (noneallowed_)
    stream << "none";
  else
    stream << 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  int n = cond->GetInt(Name());
  if (noneallowed_ and n == -1)
    stream << "none ";
  else
  {
    if (fortranstyle_) n += 1;
    stream << n;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  std::string number;
  (*condline) >> number;

  int n;
  if (noneallowed_ and number == "none")
  {
    n = -1;
  }
  else
  {
    char* ptr;
    n = strtol(number.c_str(), &ptr, 10);
    if (ptr == number.c_str())
      dserror("failed to read number '%s' while reading variable '%s' in '%s'", number.c_str(),
          Name().c_str(), def->SectionName().c_str());
  }
  if (fortranstyle_)
  {
    if (not noneallowed_ or n != -1) n -= 1;
  }
  condition->Add(Name(), n);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorConditionComponent::IntVectorConditionComponent(
    std::string name, int length, bool fortranstyle, bool noneallowed, bool optional)
    : ConditionComponent(name),
      length_(length),
      fortranstyle_(fortranstyle),
      noneallowed_(noneallowed),
      optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  if (noneallowed_)
  {
    for (int i = 0; i < length_; ++i) stream << "none ";
  }
  else
  {
    if (fortranstyle_)
    {
      for (int i = 0; i < length_; ++i) stream << -1 << " ";
    }
    else
    {
      for (int i = 0; i < length_; ++i) stream << 0 << " ";
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::Print(
    std::ostream& stream, const DRT::Condition* cond)
{
  const std::vector<int>* v = cond->Get<std::vector<int>>(Name());
  for (unsigned i = 0; i < v->size(); ++i)
  {
    if (noneallowed_ and (*v)[i] == -1)
      stream << "none ";
    else
    {
      if (fortranstyle_)
        stream << (*v)[i] + 1 << " ";
      else
        stream << (*v)[i] + 1 << " ";
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntVectorConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  // Added this as fortranstyle input was not initialized correctly if optional_ was true.
  int initialize_value = 0;
  if (fortranstyle_) initialize_value = -1;

  std::vector<int> numbers(length_, initialize_value);

  for (int i = 0; i < length_; ++i)
  {
    std::string number;
    (*condline) >> number;

    int n;
    if (noneallowed_ and number == "none")
    {
      n = -1;
    }
    else
    {
      char* ptr;
      n = strtol(number.c_str(), &ptr, 10);
      if (ptr == number.c_str())
      {
        if (optional_ and i == 0)
        {
          // failed to read the numbers, fall back to default values
          condline =
              PushBack("", condline);  // This line has been changed to incorporate optional flag!!
          // condline = PushBack(number,condline); //Old implementation -> Not working!!!
          break;
        }
        dserror(
            "Expected %i input parameters for variable '%s' in '%s'\n"
            "or \n"
            "Failed to read number '%s' while reading variable '%s' in '%s'",
            length_, Name().c_str(), def->SectionName().c_str(), number.c_str(), Name().c_str(),
            def->SectionName().c_str());
      }
    }
    if (fortranstyle_)
    {
      if (not noneallowed_ or n != -1) n -= 1;
    }
    numbers[i] = n;
  }
  condition->Add(Name(), numbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::SetLength(int length)
{
  length_ = length;
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealConditionComponent::RealConditionComponent(std::string name)
    : ConditionComponent(name)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealConditionComponent::DefaultLine(std::ostream& stream) { stream << "0.0"; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealConditionComponent::Print(std::ostream& stream, const DRT::Condition* cond)
{
  stream << cond->GetDouble(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  double number = 0;
  (*condline) >> number;
  condition->Add(Name(), number);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorConditionComponent::RealVectorConditionComponent(
    std::string name, int length, bool optional)
    : ConditionComponent(name), length_(length), optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  for (int i = 0; i < length_; ++i) stream << "0.0 ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::Print(
    std::ostream& stream, const DRT::Condition* cond)
{
  const std::vector<double>* v = cond->Get<std::vector<double>>(Name());
  for (unsigned i = 0; i < v->size(); ++i) stream << (*v)[i] << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealVectorConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  std::vector<double> numbers(length_, 0.0);

  for (int i = 0; i < length_; ++i)
  {
    std::string number;
    (*condline) >> number;

    char* ptr;
    double n = 0.0;
    n = strtod(number.c_str(), &ptr);
    if (ptr == number.c_str())
    {
      if (optional_ and i == 0)
      {
        // failed to read the numbers, fall back to default values
        condline =
            PushBack("", condline);  // This line has been changed to incorporate optional flag!!
        // condline = PushBack(number,condline); //Old implementation -> Not working!!!
        break;
      }
      dserror(
          "Expected %i input parameters for variable '%s' in '%s'\n"
          "or \n"
          "Failed to read number '%s' while reading variable '%s' in '%s'",
          length_, Name().c_str(), def->SectionName().c_str(), number.c_str(), Name().c_str(),
          def->SectionName().c_str());
    }

    numbers[i] = n;
  }
  condition->Add(Name(), numbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::SetLength(int length)
{
  length_ = length;
  return;
}

DRT::INPUT::DirichletNeumannBundle::DirichletNeumannBundle(std::string name,
    Teuchos::RCP<IntConditionComponent> intcomp,
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> intvectsepcomp,
    std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp,
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> realvectsepcomp,
    std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp)
    : ConditionComponent(name),
      intcomp_(intcomp),
      intvectsepcomp_(intvectsepcomp),
      intvectcomp_(intvectcomp),
      realvectsepcomp_(realvectsepcomp),
      realvectcomp_(realvectcomp){};

void DRT::INPUT::DirichletNeumannBundle::DefaultLine(std::ostream& stream)
{
  intcomp_->DefaultLine(stream);
  stream << "  ";
  intvectsepcomp_[0]->DefaultLine(stream);
  stream << " ";
  intvectcomp_[0]->DefaultLine(stream);
  stream << " ";
  realvectsepcomp_[0]->DefaultLine(stream);
  stream << " ";
  realvectcomp_[0]->DefaultLine(stream);
  stream << " ";
  intvectsepcomp_[1]->DefaultLine(stream);
  stream << " ";
  intvectcomp_[1]->DefaultLine(stream);
  stream << " ";
}

void DRT::INPUT::DirichletNeumannBundle::Print(std::ostream& stream, const DRT::Condition* cond)
{
  intcomp_->Print(stream, cond);
  stream << "  ";
  intvectsepcomp_[0]->Print(stream, cond);
  stream << " ";
  intvectcomp_[0]->Print(stream, cond);
  stream << " ";
  realvectsepcomp_[0]->Print(stream, cond);
  stream << " ";
  realvectcomp_[0]->Print(stream, cond);
  stream << " ";
  intvectsepcomp_[1]->Print(stream, cond);
  stream << " ";
  intvectcomp_[1]->Print(stream, cond);
  stream << " ";
};

Teuchos::RCP<std::stringstream> DRT::INPUT::DirichletNeumannBundle::Read(ConditionDefinition* def,
    Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition)
{
  intcomp_->Read(def, condline, condition);
  int length = condition->GetInt(intcomp_->Name());

  intvectsepcomp_[0]->Read(def, condline, condition);
  intvectcomp_[0]->SetLength(length);
  intvectcomp_[0]->Read(def, condline, condition);
  realvectsepcomp_[0]->Read(def, condline, condition);
  realvectcomp_[0]->SetLength(length);
  realvectcomp_[0]->Read(def, condline, condition);
  intvectsepcomp_[1]->Read(def, condline, condition);
  intvectcomp_[1]->SetLength(length);
  intvectcomp_[1]->Read(def, condline, condition);

  return condline;
}

/* -----------------------------------------------------------------------------------------------*
 | Class IntRealBundle                                                                  ehrl 09/12|
 * -----------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | IntRealBundle::Constructor()                               ehrl 09/12|
 *----------------------------------------------------------------------*/
// TODO: More general formulation possible (A. Ehrl)
//      DirichletNeumannBundle can be replaced by a more general IntRealBundle
DRT::INPUT::IntRealBundle::IntRealBundle(std::string name,
    Teuchos::RCP<IntConditionComponent> intcomp,
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp,
    std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp,
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp,
    std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp)
    : ConditionComponent(name),
      intcomp_(intcomp),
      intsepcomp_(intsepcomp),
      intvectcomp_(intvectcomp),
      realsepcomp_(realsepcomp),
      realvectcomp_(realvectcomp){};

/*----------------------------------------------------------------------*
 | CondCompBundle::DefaultLine()                              ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntRealBundle::DefaultLine(std::ostream& stream)
{
  intcomp_->DefaultLine(stream);
  stream << " ";

  // handling of different int vectors including an optional separator
  for (unsigned int i = 0; i < intvectcomp_.size(); ++i)
  {
    if (intsepcomp_[i] != Teuchos::null)
    {
      intsepcomp_[i]->DefaultLine(stream);
      stream << " ";
    }
    intvectcomp_[i]->DefaultLine(stream);
  }

  // handling of different real vectors including an optional separator
  for (unsigned int i = 0; i < realvectcomp_.size(); ++i)
  {
    if (realsepcomp_[i] != Teuchos::null)
    {
      realsepcomp_[i]->DefaultLine(stream);
      stream << " ";
    }
    realvectcomp_[i]->DefaultLine(stream);
  }

  return;
}

/*----------------------------------------------------------------------*
 | CondCompBundle::Print()                                    ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntRealBundle::Print(std::ostream& stream, const DRT::Condition* cond)
{
  intcomp_->Print(stream, cond);
  stream << " ";

  // handling of different int vectors including an optional separator
  for (unsigned int i = 0; i < intvectcomp_.size(); ++i)
  {
    if (intsepcomp_[i] != Teuchos::null)
    {
      intsepcomp_[i]->Print(stream, cond);
      stream << " ";
    }
    intvectcomp_[i]->Print(stream, cond);
  }

  // handling of different real vectors including an optional separator
  for (unsigned int i = 0; i < realvectcomp_.size(); ++i)
  {
    if (realsepcomp_[i] != Teuchos::null)
    {
      realsepcomp_[i]->Print(stream, cond);
      stream << " ";
    }
    realvectcomp_[i]->Print(stream, cond);
  }

  return;
};

/*----------------------------------------------------------------------*
 | CondCompBundle::Read()                                     ehrl 09/12|
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntRealBundle::Read(ConditionDefinition* def,
    Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition)
{
  // Read general length of int and real vectors
  // All vectors have the same length specified by this value
  intcomp_->Read(def, condline, condition);

  // get length based on stored Name()
  int length = condition->GetInt(intcomp_->Name());

  // handling of different int vectors including an optional separator
  for (unsigned int i = 0; i < intvectcomp_.size(); ++i)
  {
    if (intsepcomp_[i] != Teuchos::null) intsepcomp_[i]->Read(def, condline, condition);

    // set length of the i-th vector component
    intvectcomp_[i]->SetLength(length);
    intvectcomp_[i]->Read(def, condline, condition);
  }

  // handling of different real vectors including an optional separator
  for (unsigned int i = 0; i < realvectcomp_.size(); ++i)
  {
    if (realsepcomp_[i] != Teuchos::null) realsepcomp_[i]->Read(def, condline, condition);

    // set length of the i-th vector component
    realvectcomp_[i]->SetLength(length);
    realvectcomp_[i]->Read(def, condline, condition);
  }

  return condline;
}

/* -----------------------------------------------------------------------------------------------*
 | Class CondCompBundle                                                                 ehrl 09/12|
 * -----------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | CondCompBundle::Constructor()                              ehrl 09/12|
 *----------------------------------------------------------------------*/
DRT::INPUT::CondCompBundle::CondCompBundle(
    std::string name, std::vector<Teuchos::RCP<ConditionComponent>> condcomp, int model)
    : ConditionComponent(name), condcomp_(condcomp), model_(model)
{
}


/*----------------------------------------------------------------------*
 | CondCompBundle::DefaultLine()                              ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::CondCompBundle::DefaultLine(std::ostream& stream)
{
  // default line of selected condition component bundle
  for (unsigned int i = 0; i < condcomp_.size(); ++i)
  {
    condcomp_[i]->DefaultLine(stream);
    stream << " ";
  }

  return;
}

/*----------------------------------------------------------------------*
 | CondCompBundle::Print()                                    ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::CondCompBundle::Print(std::ostream& stream, const DRT::Condition* cond)
{
  // printing of selected condition component bundle
  for (unsigned int i = 0; i < condcomp_.size(); ++i)
  {
    condcomp_[i]->Print(stream, cond);
    stream << " ";
  }

  return;
}

/*----------------------------------------------------------------------*
 | CondCompBundle::Read()                                     ehrl 09/12|
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::CondCompBundle::Read(ConditionDefinition* def,
    Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition)
{
  // reading of selected condition component bundle
  for (unsigned int i = 0; i < condcomp_.size(); ++i) condcomp_[i]->Read(def, condline, condition);

  return condline;
}

/* -----------------------------------------------------------------------------------------------*
 | Class CondCompBundleSelector                                                         ehrl 09/12|
 * -----------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | CondCompBundleSelector::Constructor()                      ehrl 09/12|
 *----------------------------------------------------------------------*/
DRT::INPUT::CondCompBundleSelector::CondCompBundleSelector(std::string name,
    Teuchos::RCP<StringConditionComponent> stringcomp,
    std::vector<Teuchos::RCP<CondCompBundle>> condcomp)
    : ConditionComponent(name), stringcomp_(stringcomp), condcomp_(condcomp)
{
}

/*----------------------------------------------------------------------*
 | CondCompBundleSelector::DefaultLine()                      ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::CondCompBundleSelector::DefaultLine(std::ostream& stream)
{
  // Attention: default value defined here may not be identical to the printed condition component
  // bundle
  stringcomp_->DefaultLine(stream);
  stream << " ";

  // compare default std::string with std::vector<std::string> of model types
  std::ostringstream str;
  stringcomp_->DefaultLine(str);
  std::string defaultvalue = str.str();

  for (unsigned int ii = 0; ii < condcomp_.size(); ++ii)
  {
    if (defaultvalue.compare((condcomp_[ii])->Name()) == 0)
    {
      // print default condition component bundle (default bundle)
      condcomp_[ii]->DefaultLine(stream);
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | CondCompBundleSelector::Print()                            ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::CondCompBundleSelector::Print(std::ostream& stream, const DRT::Condition* cond)
{
  // Attention: default value defined here may not be identical to the printed condition component
  // bundle
  stringcomp_->Print(stream, cond);
  stream << " ";

  // compare default std::string with std::vector<std::string> of model types
  std::ostringstream str;
  stringcomp_->DefaultLine(str);
  std::string defaultvalue = str.str();

  for (unsigned int ii = 0; ii < condcomp_.size(); ++ii)
  {
    if (defaultvalue.compare((condcomp_[ii])->Name()) == 0)
    {
      // print default condition component bundle (default bundle)
      condcomp_[ii]->DefaultLine(stream);
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | CondCompBundleSelector::Read()                             ehrl 09/12|
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::CondCompBundleSelector::Read(ConditionDefinition* def,
    Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition)
{
  stringcomp_->Read(def, condline, condition);
  // get model (number is associated with a enum)
  const int model = condition->GetInt(stringcomp_->Name());

  // check if model defined in the condition match model defined in CondCompBundle
  // safety check, if models in condcomp_ are ordered in the same way as the enum defined by you
  if (model != condcomp_[model]->Model())
    dserror(
        "The model defined in your dat-file does not match the model type stored for the "
        "CondCompBundle.\n"
        "Probably, the order of the CondCompBundle in std::vector<CondCompBundle> does not match \n"
        "the model order defined in the enum!!");

  // read associated parameters
  condcomp_[model]->Read(def, condline, condition);

  return condline;
}

/* -----------------------------------------------------------------------------------------------*
 | Class ConditionDefinition                                                                      |
 * -----------------------------------------------------------------------------------------------*/

DRT::INPUT::ConditionDefinition::ConditionDefinition(std::string sectionname,
    std::string conditionname, std::string description, DRT::Condition::ConditionType condtype,
    bool buildgeometry, DRT::Condition::GeometryType gtype)
    : sectionname_(sectionname),
      conditionname_(conditionname),
      description_(description),
      condtype_(condtype),
      buildgeometry_(buildgeometry),
      gtype_(gtype)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::ConditionDefinition::AddComponent(Teuchos::RCP<ConditionComponent> c)
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::ConditionDefinition::Read(const Problem& problem, DatFileReader& reader,
    std::multimap<int, Teuchos::RCP<DRT::Condition>>& cmap)
{
  std::string name = "--";
  name += sectionname_;
  std::vector<const char*> section = reader.Section(name);

  if (section.size() > 0)
  {
    std::stringstream line(section[0]);
    std::string dobj;
    int condcount = -1;
    line >> dobj;
    line >> condcount;

    if (condcount < 0)
    {
      dserror("condcount<0 in section %s\n", sectionname_.c_str());
    }

    bool success = false;
    switch (gtype_)
    {
      case DRT::Condition::Point:
        success = dobj == "DPOINT";
        break;
      case DRT::Condition::Line:
        success = dobj == "DLINE";
        break;
      case DRT::Condition::Surface:
        success = dobj == "DSURF";
        break;
      case DRT::Condition::Volume:
        success = dobj == "DVOL";
        break;
      case DRT::Condition::Particle:
        success = dobj == "DPARTICLE";
        break;
      default:
        dserror("geometry type unspecified");
        break;
    }

    if (not success)
    {
      dserror(
          "expected design object type but got '%s' in '%s'", dobj.c_str(), sectionname_.c_str());
    }

    if (condcount != static_cast<int>(section.size() - 1))
    {
      dserror("Got %d condition lines but expect %d in '%s'", section.size() - 1, condcount,
          sectionname_.c_str());
    }

    for (std::vector<const char*>::iterator i = section.begin() + 1; i != section.end(); ++i)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(*i));

      std::string E;
      std::string number;
      std::string minus;

      (*condline) >> E >> number >> minus;
      if (not(*condline) or E != "E" or minus != "-")
        dserror("invalid condition line in '%s'", sectionname_.c_str());

      int dobjid = -1;
      {
        char* ptr;
        dobjid = strtol(number.c_str(), &ptr, 10);
        if (ptr == number.c_str())
          dserror("failed to read design object number '%s' in '%s'", number.c_str(),
              sectionname_.c_str());
        dobjid -= 1;
      }

      // cout << "condition " << dobjid << " of type " << condtype_ << " on " << gtype_ << "\n";

      Teuchos::RCP<DRT::Condition> condition =
          Teuchos::rcp(new DRT::Condition(dobjid, condtype_, buildgeometry_, gtype_));

      for (unsigned j = 0; j < inputline_.size(); ++j)
      {
        condline = inputline_[j]->Read(this, condline, condition);
      }

      //------------------------------- put condition in map of conditions
      cmap.insert(std::pair<int, Teuchos::RCP<DRT::Condition>>(dobjid, condition));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& DRT::INPUT::ConditionDefinition::Print(
    std::ostream& stream, const Discretization* dis, bool color)
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

  unsigned l = sectionname_.length();
  stream << redlight << "--";
  for (int i = 0; i < std::max<int>(65 - l, 0); ++i) stream << '-';
  stream << greenlight << sectionname_ << endcolor << '\n';

  std::string name;
  switch (gtype_)
  {
    case DRT::Condition::Point:
      name = "DPOINT";
      break;
    case DRT::Condition::Line:
      name = "DLINE";
      break;
    case DRT::Condition::Surface:
      name = "DSURF";
      break;
    case DRT::Condition::Volume:
      name = "DVOL";
      break;
    case DRT::Condition::Particle:
      name = "DPARTICLE";
      break;
    default:
      dserror("geometry type unspecified");
      break;
  }

  int count = 0;
  if (dis != NULL)
  {
    std::vector<Condition*> conds;
    dis->GetCondition(conditionname_, conds);
    for (unsigned c = 0; c < conds.size(); ++c)
    {
      if (conds[c]->GType() == gtype_)
      {
        count += 1;
      }
    }
  }

  stream << bluelight << name << endcolor;
  l = name.length();
  for (int i = 0; i < std::max<int>(31 - l, 0); ++i) stream << ' ';
  stream << ' ' << yellowlight << count << endcolor << '\n';

  stream << blue2light << "//" << magentalight << "E num - ";
  for (unsigned i = 0; i < inputline_.size(); ++i)
  {
    inputline_[i]->DefaultLine(stream);
    stream << " ";
  }

  stream << endcolor << "\n";

  if (dis != NULL)
  {
    std::vector<Condition*> conds;
    dis->GetCondition(conditionname_, conds);

    for (unsigned c = 0; c < conds.size(); ++c)
    {
      if (conds[c]->GType() == gtype_)
      {
        stream << "E " << conds[c]->Id() << " - ";
        for (unsigned i = 0; i < inputline_.size(); ++i)
        {
          inputline_[i]->Print(stream, conds[c]);
          stream << " ";
        }
        stream << "\n";
      }
    }
  }

  return stream;
}
