/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for all conditions

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_CONDITIONDEFINITION_HPP
#define FOUR_C_LIB_CONDITIONDEFINITION_HPP

#include "4C_config.hpp"

#include "4C_inpar_container.hpp"
#include "4C_io_inputreader.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_lib_condition.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace GLOBAL
{
  class Problem;
}

namespace INPUT
{

  //--------------------------------------------------------------

  /// definition of a valid condition in 4C input
  /*!

    A ConditionDefinition is the definition of a condition dat file
    section. This definition includes the knowledge what this section looks
    like, how to read it and how to write it. In particular given a
    ConditionDefinition object it is possible to (a) write an empty dat file
    section that describes this condition, (b) read a dat file and create
    DRT::Condition objects for each line in this section and (c) write the dat
    file section filled with all corresponding conditions from a given
    DRT::Discretization.

    So this is quite sophisticated internal stuff here. If you want to
    introduce a new condition to 4C, all you have to do is add an
    appropriate definition in ValidConditions(). This will take care of the
    reading part and you will get your DRT::Discretization filled with proper
    DRT::Condition objects.

    \author u.kue
    \date 01/08
   */
  class ConditionDefinition
  {
   public:
    /// construction of a condition definition
    /*!
      \param sectionname name of dat file section
      \param conditionname name of conditions in DRT::Discretization
      \param description description of condition type
      \param condtype type of conditions to be build
      \param buildgeometry whether we need conditions elements
      \param gtype type of geometry the condition lives on
     */
    ConditionDefinition(std::string sectionname, std::string conditionname, std::string description,
        DRT::Condition::ConditionType condtype, bool buildgeometry,
        DRT::Condition::GeometryType gtype);

    /// add a concrete component to the condition line definition
    /*!
      Add new components to the input line. One at a time. Form left to
      right. The order is important! On reading we try and read component
      after component.
     */
    void AddComponent(const Teuchos::RCP<INPUT::LineComponent>& c);

    /// read all conditions from my input file section
    /*!
      \param problem (i) global problem instance that manages the input
      \param reader (i) the actual dat file reader that has access to the dat file
      \param cmap (o) the conditions we read here
     */
    void Read(const GLOBAL::Problem& problem, DatFileReader& reader,
        std::multimap<int, Teuchos::RCP<DRT::Condition>>& cmap);

    /// print my dat file section and possible conditions from the Discretization
    std::ostream& Print(std::ostream& stream, const DRT::Discretization* dis = nullptr);

    /// name of my section in input file
    std::string SectionName() const { return sectionname_; }

    /// my condition name
    std::string Name() const { return conditionname_; }

    /// my condition description
    std::string Description() const { return description_; }

    /// my condition inputline
    std::vector<Teuchos::RCP<INPUT::LineComponent>> Inputline() const { return inputline_; }

    /// my GeometryType
    DRT::Condition::GeometryType GeometryType() const { return gtype_; }

   private:
    std::string sectionname_;
    std::string conditionname_;
    std::string description_;
    DRT::Condition::ConditionType condtype_;
    bool buildgeometry_;
    DRT::Condition::GeometryType gtype_;

    ///
    std::vector<Teuchos::RCP<INPUT::LineComponent>> inputline_;
  };

}  // namespace INPUT


FOUR_C_NAMESPACE_CLOSE

#endif
