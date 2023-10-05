/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of base class for all conditions

\level 0


*/
/*---------------------------------------------------------------------*/


#include "baci_lib_conditiondefinition.H"

#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_linecomponent.H"
#include "baci_lib_utils_cond_and_mat_definition.H"
#include "baci_lib_utils_reader.H"

#include <algorithm>
#include <iterator>
#include <utility>

/* -----------------------------------------------------------------------------------------------*
 | Class StringConditionComponent                                                       ehrl 09/12|
 * -----------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | StringConditionComponent::Constructor()                              |
 *----------------------------------------------------------------------*/
DRT::INPUT::StringConditionComponent::StringConditionComponent(std::string name,
    std::string defaultvalue, const Teuchos::Array<std::string>& datfilevalues,
    const Teuchos::Array<std::string>& stringcondvalues, bool optional)
    : ::INPUT::LineComponent(std::move(name)),
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
    : ::INPUT::LineComponent(std::move(name)),
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
void DRT::INPUT::StringConditionComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  stream << *cond.Get<std::string>(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::StringConditionComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
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
                  << section_name.c_str() << "'.\n"
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
    condition.Add(Name(), stringcondvalues_[pos]);
  else
    condition.Add(Name(), intcondvalues_[pos]);

  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::SeparatorConditionComponent::SeparatorConditionComponent(
    std::string separator, bool optional)
    : ::INPUT::LineComponent("*SEPARATOR*"), separator_(std::move(separator)), optional_(optional)
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
    std::ostream& stream, const DRT::Container& cond)
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
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
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
      Read(section_name, condline, condition);
    }
    else
    {
      dserror("word '%s' expected but found '%s' while reading '%s'", separator_.c_str(),
          sep.c_str(), section_name.c_str());
    }
  }
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntConditionComponent::IntConditionComponent(
    std::string name, bool fortranstyle, bool noneallowed, bool optional)
    : ::INPUT::LineComponent(std::move(name)),
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
void DRT::INPUT::IntConditionComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  int n = cond.GetInt(Name());
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
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
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
        snumber, Name(), section_name, 1, optional_);
  }
  if (fortranstyle_)
  {
    if (not noneallowed_ or nnumber != -1) nnumber -= 1;
  }

  condition.Add(Name(), nnumber);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::IntVectorConditionComponent::IntVectorConditionComponent(
    std::string name, int length, bool fortranstyle, bool noneallowed, bool optional)
    : ::INPUT::LineComponent(std::move(name)),
      length_(length),
      fortranstyle_(fortranstyle),
      noneallowed_(noneallowed),
      optional_(optional)
{
}


DRT::INPUT::IntVectorConditionComponent::IntVectorConditionComponent(std::string name,
    LengthDefinition length_from_component, bool fortranstyle, bool noneallowed, bool optional)
    : ::INPUT::LineComponent(std::move(name)),
      length_(std::move(length_from_component)),
      fortranstyle_(fortranstyle),
      noneallowed_(noneallowed),
      optional_(optional)
{
}


namespace
{
  struct DefaultLengthVisitor
  {
    int operator()(int length) { return length; }
    int operator()(const DRT::INPUT::LengthDefinition& length) { return 1; }
  };
}  // namespace


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  const int default_length = std::visit(DefaultLengthVisitor{}, length_);
  const char* default_value = std::invoke(
      [&]()
      {
        if (noneallowed_) return "none ";
        if (fortranstyle_)
          return "-1 ";
        else
          return "0 ";
      });

  for (int i = 0; i < default_length; ++i) stream << default_value;
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
    std::ostream& stream, const DRT::Container& cond)
{
  const auto* v = cond.Get<std::vector<int>>(Name());
  for (int i : *v)
  {
    if (noneallowed_ and i == -1)
      stream << "none ";
    else
      stream << i + 1 << " ";
  }
}


namespace
{
  struct LengthVisitor
  {
    int operator()(int length) { return length; }
    int operator()(const DRT::INPUT::LengthDefinition& length) { return length(condition); }

    const DRT::Container& condition;
  };
}  // namespace
   //

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::IntVectorConditionComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
{
  // in order to initialize fortran style input correctly if optional_ is true
  const int initialize_value = fortranstyle_ ? -1 : 0;
  const int dynamic_length = std::visit(LengthVisitor{condition}, length_);
  std::vector<int> nnumbers(dynamic_length, initialize_value);

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
          snumber, Name(), section_name, dynamic_length, optional_);
    }

    if (fortranstyle_)
    {
      if (not noneallowed_ or current_nnumber != -1) current_nnumber -= 1;
    }
  }

  condition.Add(Name(), nnumbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::IntVectorConditionComponent::SetLength(int length) { length_ = length; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealConditionComponent::RealConditionComponent(std::string name)
    : ::INPUT::LineComponent(std::move(name))
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
void DRT::INPUT::RealConditionComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  stream << cond.GetDouble(Name());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealConditionComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
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
        snumber, Name(), section_name, 1, true);
  }

  condition.Add(Name(), nnumber);
  return condline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::INPUT::RealVectorConditionComponent::RealVectorConditionComponent(
    std::string name, int length, bool optional)
    : ::INPUT::LineComponent(std::move(name)), length_(length), optional_(optional)
{
}



DRT::INPUT::RealVectorConditionComponent::RealVectorConditionComponent(
    std::string name, LengthDefinition length, bool optional)
    : ::INPUT::LineComponent(std::move(name)), length_(std::move(length)), optional_(optional)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::DefaultLine(std::ostream& stream)
{
  const int default_length = std::visit(DefaultLengthVisitor{}, length_);
  for (int i = 0; i < default_length; ++i) stream << "0.0 ";
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
    std::ostream& stream, const DRT::Container& cond)
{
  const auto* v = cond.Get<std::vector<double>>(Name());
  for (double i : *v) stream << i << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::stringstream> DRT::INPUT::RealVectorConditionComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
{
  const int dynamic_length = std::visit(LengthVisitor{condition}, length_);
  std::vector<double> nnumbers(dynamic_length, 0.0);

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
          snumber, Name(), section_name, dynamic_length, optional_);
    }
  }

  condition.Add(Name(), nnumbers);
  return condline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::INPUT::RealVectorConditionComponent::SetLength(int length) { length_ = length; }



DRT::INPUT::SwitchConditionComponent::SwitchConditionComponent(std::string name,
    const KeyType& default_key,
    std::map<KeyType, std::pair<std::string, std::vector<Teuchos::RCP<::INPUT::LineComponent>>>>
        choices)
    : ::INPUT::LineComponent(std::move(name)),
      default_key_(default_key),
      choices_(std::move(choices))
{
  Teuchos::Array<int> keys;
  keys.reserve(choices_.size());
  Teuchos::Array<std::string> names_for_keys;
  names_for_keys.reserve(choices_.size());

  for (const auto& [key, choice] : choices_)
  {
    keys.push_back(key);
    names_for_keys.push_back(choice.first);
  }

  component_for_key_ = std::make_unique<StringConditionComponent>(
      Name(), choices_[default_key_].first, names_for_keys, keys);
}



void DRT::INPUT::SwitchConditionComponent::DefaultLine(std::ostream& stream)
{
  component_for_key_->DefaultLine(stream);
  stream << " ";

  for (const auto& component : choices_[default_key_].second)
  {
    component->DefaultLine(stream);
    stream << " ";
  }
}



std::vector<std::string> DRT::INPUT::SwitchConditionComponent::WriteReadTheDocsLines()
{
  std::vector<std::string> all_choices_as_rtd;
  std::transform(choices_.begin(), choices_.end(), std::back_inserter(all_choices_as_rtd),
      [this](const auto& key_components)
      {
        const auto& [key, components] = key_components;

        std::stringstream stream;
        stream << choices_[key].first << " ";
        for (const auto& c : components.second) stream << c->WriteReadTheDocs() << " ";
        return stream.str();
      });

  return all_choices_as_rtd;
}


std::string DRT::INPUT::SwitchConditionComponent::WriteReadTheDocs()
{
  return component_for_key_->WriteReadTheDocs() + " [further parameters]";
}



Teuchos::Array<std::string> DRT::INPUT::SwitchConditionComponent::GetOptions()
{
  return component_for_key_->GetOptions();
}



void DRT::INPUT::SwitchConditionComponent::Print(std::ostream& stream, const DRT::Container& cond)
{
  component_for_key_->Print(stream, cond);
  stream << " ";

  const KeyType selected_key = static_cast<KeyType>(cond.GetInt(component_for_key_->Name()));

  dsassert(choices_.count(selected_key) == 1, "Internal error.");
  for (const auto& component : choices_[selected_key].second)
  {
    component->Print(stream, cond);
    stream << " ";
  }
}



Teuchos::RCP<std::stringstream> DRT::INPUT::SwitchConditionComponent::Read(
    const std::string& section_name, Teuchos::RCP<std::stringstream> condline,
    DRT::Container& condition)
{
  component_for_key_->Read(section_name, condline, condition);
  const KeyType key = static_cast<KeyType>(condition.GetInt(component_for_key_->Name()));

  dsassert(choices_.count(key) == 1, "Internal error.");

  for (const auto& component : choices_[key].second)
  {
    component->Read(section_name, condline, condition);
  }

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
void DRT::INPUT::ConditionDefinition::AddComponent(const Teuchos::RCP<::INPUT::LineComponent>& c)
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
              "specified number format is probably not supported. The design object number has "
              "to "
              "be an integer.",
              number.c_str(), sectionname_.c_str(), ptr, number.c_str());
        }

        dobjid -= 1;
      }

      Teuchos::RCP<DRT::Condition> condition =
          Teuchos::rcp(new DRT::Condition(dobjid, condtype_, buildgeometry_, gtype_));

      for (auto& j : inputline_)
      {
        condline = j->Read(SectionName(), condline, *condition);
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
          i->Print(stream, *cond);
          stream << " ";
        }
        stream << "\n";
      }
    }
  }

  return stream;
}
