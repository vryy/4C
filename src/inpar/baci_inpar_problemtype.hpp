/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_PROBLEMTYPE_HPP
#define FOUR_C_INPAR_PROBLEMTYPE_HPP

#include "baci_config.hpp"

#include "baci_global_data_enums.hpp"
#include "baci_utils_parameter_list.hpp"

#include <map>
#include <string>

BACI_NAMESPACE_OPEN

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
    std::map<std::string, GLOBAL::ProblemType> StringToProblemTypeMap();

    /// return problem type enum for a given problem name
    GLOBAL::ProblemType StringToProblemType(std::string name);


  }  // namespace PROBLEMTYPE
}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
