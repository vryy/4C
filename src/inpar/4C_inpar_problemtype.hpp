/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_PROBLEMTYPE_HPP
#define FOUR_C_INPAR_PROBLEMTYPE_HPP

#include "4C_config.hpp"

#include "4C_legacy_enum_definitions_problem_type.hpp"
#include "4C_utils_parameter_list.hpp"

#include <map>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace INPAR
{
  namespace PROBLEMTYPE
  {
    /*! \brief Define valid parameters
     *
     * @param[in/out] list Parameter list to be filled with valid parameters and their defaults
     */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// create map of problem name and problem type enum
    std::map<std::string, CORE::ProblemType> StringToProblemTypeMap();

    /// return problem type enum for a given problem name
    CORE::ProblemType StringToProblemType(std::string name);


  }  // namespace PROBLEMTYPE
}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
