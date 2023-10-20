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

#include <algorithm>
#include <iterator>
#include <utility>



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

      // add trailing white space to stringstream "condline" to avoid deletion of stringstream upon
      // reading the last entry inside This is required since the material parameters can be
      // specified in an arbitrary order in the input file. So it might happen that the last entry
      // is extracted before all of the previous ones are.
      condline->seekp(0, condline->end);
      *condline << " ";
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
