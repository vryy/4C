/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of base class for all conditions

\level 0


*/
/*---------------------------------------------------------------------*/


#include "4C_lib_conditiondefinition.hpp"

#include "4C_global_data.hpp"
#include "4C_io_line_parser.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_lib_discret.hpp"
#include "4C_utils_demangle.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <iterator>
#include <utility>

FOUR_C_NAMESPACE_OPEN



/* -----------------------------------------------------------------------------------------------*
 | Class ConditionDefinition                                                                      |
 * -----------------------------------------------------------------------------------------------*/

INPUT::ConditionDefinition::ConditionDefinition(std::string sectionname, std::string conditionname,
    std::string description, CORE::Conditions::ConditionType condtype, bool buildgeometry,
    CORE::Conditions::GeometryType gtype)
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
void INPUT::ConditionDefinition::AddComponent(const Teuchos::RCP<INPUT::LineComponent>& c)
{
  inputline_.push_back(c);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void INPUT::ConditionDefinition::Read(const GLOBAL::Problem& problem, DatFileReader& reader,
    std::multimap<int, Teuchos::RCP<DRT::Condition>>& cmap)
{
  std::vector<const char*> section = reader.Section("--" + sectionname_);

  if (section.empty()) return;

  IO::LineParser parser("While reading condition section '" + sectionname_ + "': ");

  // First we read a header for the current section: It needs to start with the
  // geometry type followed by the number of lines:
  //
  // ("DPOINT" | "DLINE" | "DSURF" | "DVOL" ) <number>
  {
    std::stringstream line(section[0]);

    const std::string expected_geometry_type = std::invoke(
        [this]()
        {
          switch (gtype_)
          {
            case CORE::Conditions::geometry_type_point:
              return "DPOINT";
            case CORE::Conditions::geometry_type_line:
              return "DLINE";
            case CORE::Conditions::geometry_type_surface:
              return "DSURF";
            case CORE::Conditions::geometry_type_volume:
              return "DVOL";
            default:
              FOUR_C_THROW("Geometry type unspecified");
          }
        });

    parser.Consume(line, expected_geometry_type);
    const int condition_count = parser.Read<int>(line);

    if (condition_count != static_cast<int>(section.size() - 1))
    {
      FOUR_C_THROW("Got %d condition lines but expected %d in section '%s'", section.size() - 1,
          condition_count, sectionname_.c_str());
    }
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

    parser.Consume(*condline, "E");
    // Read a one-based condition number but convert it to zero-based for internal use.
    const int dobjid = parser.Read<int>(*condline) - 1;
    parser.Consume(*condline, "-");

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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& INPUT::ConditionDefinition::Print(
    std::ostream& stream, const DRT::Discretization* dis)
{
  unsigned l = sectionname_.length();
  stream << "--";
  for (int i = 0; i < std::max<int>(65 - l, 0); ++i) stream << '-';
  stream << sectionname_ << '\n';

  std::string name;
  switch (gtype_)
  {
    case CORE::Conditions::geometry_type_point:
      name = "DPOINT";
      break;
    case CORE::Conditions::geometry_type_line:
      name = "DLINE";
      break;
    case CORE::Conditions::geometry_type_surface:
      name = "DSURF";
      break;
    case CORE::Conditions::geometry_type_volume:
      name = "DVOL";
      break;
    default:
      FOUR_C_THROW("geometry type unspecified");
      break;
  }

  int count = 0;
  if (dis != nullptr)
  {
    std::vector<DRT::Condition*> conds;
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
    std::vector<DRT::Condition*> conds;
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

FOUR_C_NAMESPACE_CLOSE
