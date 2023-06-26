/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of base class for all conditions

\level 0


*/
/*---------------------------------------------------------------------*/


#include <algorithm>
#include <iterator>
#include <utility>
#include "lib_conditiondefinition.H"
#include "lib_globalproblem.H"
#include "lib_discret.H"
#include "lib_utils_cond_and_mat_definition.H"
#include "create_rtdfiles_utils.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::ConditionComponent::ConditionComponent(std::string name) : name_(std::move(name)) {}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::ConditionComponent::PushBack(
    const std::string& token, const Teuchos::RCP<std::stringstream>& stream)
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
    : ConditionComponent(std::move(name)),
      defaultvalue_(std::move(defaultvalue)),
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
    : ConditionComponent(std::move(name)),
      defaultvalue_(std::move(defaultvalue)),
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

std::string DRT::INPUT::StringConditionComponent::WriteReadTheDocs() { return "<" + Name() + ">"; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::StringConditionComponent::GetOptions()
{
  return datfilevalues_;
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

  if (value.empty()) value = defaultvalue_;

  auto i = std::find(datfilevalues_.begin(), datfilevalues_.end(), value);
  if (i == datfilevalues_.end())
  {
    std::stringstream error_message;
    error_message << "\n \nUnrecognized std::string '" << value.c_str()
                  << "' while reading variable '" << Name().c_str() << "' in '"
                  << def->SectionName().c_str() << "'.\n"
                  << "Possible values are: \n";

    for (const auto& datfilevalue : datfilevalues_)
    {
      error_message << std::string(23, ' ') << datfilevalue << "\n";
    }

    dserror(error_message.str());
  }
  const unsigned pos = std::distance(datfilevalues_.begin(), i);
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
    : ConditionComponent("*SEPARATOR*"), separator_(std::move(separator)), optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorConditionComponent::DefaultLine(std::ostream& stream)
{
  stream << separator_;
}

std::string DRT::INPUT::SeparatorConditionComponent::WriteReadTheDocs() { return separator_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::SeparatorConditionComponent::Print(
    std::ostream& stream, const DRT::Condition* cond)
{
  stream << separator_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::SeparatorConditionComponent::GetOptions()
{
  return Teuchos::Array<std::string>();
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
    if (sep.empty() && optional_)
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
    std::string name, bool fortranstyle, bool noneallowed, bool optional)
    : ConditionComponent(std::move(name)),
      fortranstyle_(fortranstyle),
      noneallowed_(noneallowed),
      optional_(optional)
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

std::string DRT::INPUT::IntConditionComponent::WriteReadTheDocs()
{
  return (noneallowed_) ? "none" : "0";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::IntConditionComponent::GetOptions()
{
  return Teuchos::Array<std::string>();
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
  // read string from stream (need to handle ints and "none" arguments)
  std::string snumber;
  (*condline) >> snumber;

  int nnumber = 0;

  // in case 'none' is allowed as an input value
  if ((noneallowed_ and snumber == "none"))
  {
    nnumber = -1;
  }
  // in case the parameter is optional and no value is given
  else if (optional_ and snumber.empty())
  {
    condline = PushBack("", condline);
  }
  // all other cases
  else
  {
    nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<int>(
        snumber, Name(), def->SectionName(), 1, optional_);
  }
  if (fortranstyle_)
  {
    if (not noneallowed_ or nnumber != -1) nnumber -= 1;
  }

  condition->Add(Name(), nnumber);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorConditionComponent::IntVectorConditionComponent(
    std::string name, int length, bool fortranstyle, bool noneallowed, bool optional)
    : ConditionComponent(std::move(name)),
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

std::string DRT::INPUT::IntVectorConditionComponent::WriteReadTheDocs()
{
  std::string parameterstring = "<int vec";
  if (noneallowed_) parameterstring += " [incl none]";
  parameterstring += ":" + Name() + "> ";
  return parameterstring;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::IntVectorConditionComponent::GetOptions()
{
  return Teuchos::Array<std::string>();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::Print(
    std::ostream& stream, const DRT::Condition* cond)
{
  const auto* v = cond->Get<std::vector<int>>(Name());
  for (int i : *v)
  {
    if (noneallowed_ and i == -1)
      stream << "none ";
    else
      stream << i + 1 << " ";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntVectorConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  // in order to initialize fortran style input correctly if optional_ is true
  int initialize_value = 0;
  if (fortranstyle_) initialize_value = -1;

  std::vector<int> nnumbers(length_, initialize_value);

  for (auto& current_nnumber : nnumbers)
  {
    // read string from stream (need to handle ints and "none" arguments)
    std::string snumber;
    (*condline) >> snumber;

    // in case 'none' is allowed as an input value
    if (noneallowed_ and snumber == "none")
    {
      current_nnumber = -1;
    }
    // in case the parameter is optional and no value is given
    else if (optional_ and snumber.empty())
    {
      condline = PushBack("", condline);
      break;
    }
    // all other cases
    else
    {
      current_nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<int>(
          snumber, Name(), def->SectionName(), length_, optional_);
    }

    if (fortranstyle_)
    {
      if (not noneallowed_ or current_nnumber != -1) current_nnumber -= 1;
    }
  }

  condition->Add(Name(), nnumbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::SetLength(int length) { length_ = length; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealConditionComponent::RealConditionComponent(std::string name)
    : ConditionComponent(std::move(name))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealConditionComponent::DefaultLine(std::ostream& stream) { stream << "0.0"; }

std::string DRT::INPUT::RealConditionComponent::WriteReadTheDocs() { return "0.0"; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::RealConditionComponent::GetOptions()
{
  return Teuchos::Array<std::string>();
}


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
  // read string from stream
  std::string snumber;
  (*condline) >> snumber;

  double nnumber = 0.0;

  // no optional_ parameter for RealConditionComponent, parameter is always optional
  // hence no check for optional but only check if no value is given
  if (snumber.empty())
  {
    condline = PushBack("", condline);
  }
  // all other cases
  else
  {
    nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<double>(
        snumber, Name(), def->SectionName(), 1, true);
  }

  condition->Add(Name(), nnumber);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorConditionComponent::RealVectorConditionComponent(
    std::string name, int length, bool optional)
    : ConditionComponent(std::move(name)), length_(length), optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  for (int i = 0; i < length_; ++i) stream << "0.0 ";
}

std::string DRT::INPUT::RealVectorConditionComponent::WriteReadTheDocs()
{
  return "<real vec:" + Name() + "> ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::RealVectorConditionComponent::GetOptions()
{
  return Teuchos::Array<std::string>();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::Print(
    std::ostream& stream, const DRT::Condition* cond)
{
  const auto* v = cond->Get<std::vector<double>>(Name());
  for (double i : *v) stream << i << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealVectorConditionComponent::Read(
    DRT::INPUT::ConditionDefinition* def, Teuchos::RCP<std::stringstream> condline,
    Teuchos::RCP<DRT::Condition> condition)
{
  std::vector<double> nnumbers(length_, 0.0);

  for (auto& current_nnumber : nnumbers)
  {
    // read string from stream
    std::string snumber;
    (*condline) >> snumber;

    // in case the parameter is optional and no value is
    // given
    if (optional_ and snumber.empty())
    {
      condline = PushBack("", condline);
      break;
    }
    // all other cases
    else
    {
      current_nnumber = DRT::UTILS::ConvertAndValidateStringToNumber<double>(
          snumber, Name(), def->SectionName(), length_, optional_);
    }
  }

  condition->Add(Name(), nnumbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::SetLength(int length) { length_ = length; }

DRT::INPUT::DirichletNeumannBundle::DirichletNeumannBundle(std::string name,
    Teuchos::RCP<IntConditionComponent> intcomp,
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> intvectsepcomp,
    std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp,
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> realvectsepcomp,
    std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp)
    : ConditionComponent(std::move(name)),
      intcomp_(std::move(intcomp)),
      intvectsepcomp_(std::move(intvectsepcomp)),
      intvectcomp_(std::move(intvectcomp)),
      realvectsepcomp_(std::move(realvectsepcomp)),
      realvectcomp_(std::move(realvectcomp))
{
}

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

std::string DRT::INPUT::DirichletNeumannBundle::WriteReadTheDocs()
{
  std::string parameterstring = "";
  // numdof int
  parameterstring += intcomp_->WriteReadTheDocs() + " ";
  // ONOFF
  parameterstring += intvectsepcomp_[0]->WriteReadTheDocs() + " ";
  // onoff vector
  parameterstring += intvectcomp_[0]->WriteReadTheDocs() + " ";
  // VAL
  parameterstring += realvectsepcomp_[0]->WriteReadTheDocs() + " ";
  // val vector
  parameterstring += realvectcomp_[0]->WriteReadTheDocs() + " ";
  // FUNCT
  parameterstring += intvectsepcomp_[1]->WriteReadTheDocs() + " ";
  // funct vector
  parameterstring += intvectcomp_[1]->WriteReadTheDocs() + " ";
  return parameterstring;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::DirichletNeumannBundle::GetOptions()
{
  return Teuchos::Array<std::string>();
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
}

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
    : ConditionComponent(std::move(name)),
      intcomp_(std::move(intcomp)),
      intsepcomp_(std::move(intsepcomp)),
      intvectcomp_(std::move(intvectcomp)),
      realsepcomp_(std::move(realsepcomp)),
      realvectcomp_(std::move(realvectcomp))
{
}

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
}

std::string DRT::INPUT::IntRealBundle::WriteReadTheDocs()
{
  std::ostringstream parameterstream;
  DRT::INPUT::IntRealBundle::DefaultLine(parameterstream);
  return parameterstream.str();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::IntRealBundle::GetOptions()
{
  return Teuchos::Array<std::string>();
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
}

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
    : ConditionComponent(std::move(name)), condcomp_(std::move(condcomp)), model_(model)
{
}


/*----------------------------------------------------------------------*
 | CondCompBundle::DefaultLine()                              ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::CondCompBundle::DefaultLine(std::ostream& stream)
{
  // default line of selected condition component bundle
  for (auto& i : condcomp_)
  {
    i->DefaultLine(stream);
    stream << " ";
  }
}

std::string DRT::INPUT::CondCompBundle::WriteReadTheDocs()
{
  std::string parameterstring;
  for (auto& component : condcomp_)
  {
    parameterstring += component->WriteReadTheDocs() + " ";
  }
  return parameterstring;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::CondCompBundle::GetOptions()
{
  return Teuchos::Array<std::string>();
}


/*----------------------------------------------------------------------*
 | CondCompBundle::Print()                                    ehrl 09/12|
 *----------------------------------------------------------------------*/
void DRT::INPUT::CondCompBundle::Print(std::ostream& stream, const DRT::Condition* cond)
{
  // printing of selected condition component bundle
  for (auto& i : condcomp_)
  {
    i->Print(stream, cond);
    stream << " ";
  }
}

/*----------------------------------------------------------------------*
 | CondCompBundle::Read()                                     ehrl 09/12|
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::CondCompBundle::Read(ConditionDefinition* def,
    Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition)
{
  // reading of selected condition component bundle
  for (auto& i : condcomp_) i->Read(def, condline, condition);

  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::CondCompBundleSelector::CondCompBundleSelector(std::string name_condition_components,
    const std::vector<Teuchos::RCP<CondCompBundle>>& condcomp)
    : ConditionComponent(name_condition_components),  // + "_selector_internal"),
      stringcomp_(),
      condcomp_()
{
  Teuchos::Array<std::string> names;
  Teuchos::Array<int> models;

  for (const auto& component : condcomp)
  {
    if (condcomp_.find(component->Model()) != condcomp_.end())
      dserror("condition model number is not unique.");

    condcomp_.emplace(component->Model(), component);
    names.push_back(component->Name());
    models.push_back(component->Model());
  }

  stringcomp_ = Teuchos::rcp(
      new StringConditionComponent(std::move(name_condition_components), names[0], names, models));
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

  for (auto& ii : condcomp_)
  {
    if (defaultvalue == ii.second->Name())
    {
      // print default condition component bundle (default bundle)
      ii.second->DefaultLine(stream);
      break;
    }
  }
}
/*----------------------------------------------------------------------*
| CondCompBundleSelector::DefaultLines()                      ische 03/21|
*----------------------------------------------------------------------*/
std::vector<std::string> DRT::INPUT::CondCompBundleSelector::WriteReadTheDocsLines()
{
  std::vector<std::string> condCompStrings;
  for (auto& componentBundle : condcomp_)
  {
    std::string condCompString(componentBundle.second->Name() + " ");
    // print default condition component bundle (default bundle)
    condCompString += componentBundle.second->WriteReadTheDocs();
    condCompStrings.push_back(condCompString);
  }
  return condCompStrings;
}

std::string DRT::INPUT::CondCompBundleSelector::WriteReadTheDocs()
{
  std::string parameterstring = "";
  parameterstring += stringcomp_->WriteReadTheDocs() + " [further parameters]";
  return parameterstring;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::Array<std::string> DRT::INPUT::CondCompBundleSelector::GetOptions()
{
  return stringcomp_->GetOptions();
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

  for (auto& ii : condcomp_)
  {
    if (defaultvalue == ii.second()->Name())
    {
      // print default condition component bundle (default bundle)
      ii.second->DefaultLine(stream);
      break;
    }
  }
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
    : sectionname_(std::move(sectionname)),
      conditionname_(std::move(conditionname)),
      description_(std::move(description)),
      condtype_(condtype),
      buildgeometry_(buildgeometry),
      gtype_(gtype)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::ConditionDefinition::AddComponent(const Teuchos::RCP<ConditionComponent>& c)
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

  if (!section.empty())
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

    for (auto i = section.begin() + 1; i != section.end(); ++i)
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
          dserror("Failed to read design object number '%s' in '%s'", number.c_str(),
              sectionname_.c_str());
        if (ptr[0])  // check for remaining characters that were not read
        {
          dserror(
              "Failed to read design object number '%s' in '%s'. Unable to read %s, so the "
              "specified number format is probably not supported. The design object number has to "
              "be an integer.",
              number.c_str(), sectionname_.c_str(), ptr, number.c_str());
        }

        dobjid -= 1;
      }

      // cout << "condition " << dobjid << " of type " << condtype_ << " on " << gtype_ << "\n";

      Teuchos::RCP<DRT::Condition> condition =
          Teuchos::rcp(new DRT::Condition(dobjid, condtype_, buildgeometry_, gtype_));

      for (auto& j : inputline_)
      {
        condline = j->Read(this, condline, condition);
      }

      //------------------------------- put condition in map of conditions
      cmap.insert(std::pair<int, Teuchos::RCP<DRT::Condition>>(dobjid, condition));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& DRT::INPUT::ConditionDefinition::Print(
    std::ostream& stream, const Discretization* dis)
{
  unsigned l = sectionname_.length();
  stream << "--";
  for (int i = 0; i < std::max<int>(65 - l, 0); ++i) stream << '-';
  stream << sectionname_ << '\n';

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
    default:
      dserror("geometry type unspecified");
      break;
  }

  int count = 0;
  if (dis != nullptr)
  {
    std::vector<Condition*> conds;
    dis->GetCondition(conditionname_, conds);
    for (auto& cond : conds)
    {
      if (cond->GType() == gtype_)
      {
        count += 1;
      }
    }
  }

  stream << name;
  l = name.length();
  for (int i = 0; i < std::max<int>(31 - l, 0); ++i) stream << ' ';
  stream << ' ' << count << '\n';

  stream << "//"
         << "E num - ";
  for (auto& i : inputline_)
  {
    i->DefaultLine(stream);
    stream << " ";
  }

  stream << "\n";

  if (dis != nullptr)
  {
    std::vector<Condition*> conds;
    dis->GetCondition(conditionname_, conds);

    for (auto& cond : conds)
    {
      if (cond->GType() == gtype_)
      {
        stream << "E " << cond->Id() << " - ";
        for (auto& i : inputline_)
        {
          i->Print(stream, cond);
          stream << " ";
        }
        stream << "\n";
      }
    }
  }

  return stream;
}
