/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_PROBLEMTYPE_H
#define INPAR_PROBLEMTYPE_H

#include <string>
#include <map>

#include "inpar_parameterlist_utils.H"

#include "../drt_lib/drt_globalproblem_enums.H"

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
    std::map<std::string, ProblemType> StringToProblemTypeMap();

    /// return problem type enum for a given problem name
    ProblemType StringToProblemType(std::string name);

    /// create map of shape function names and shape function type enum
    std::map<std::string, ShapeFunctionType> StringToShapeFunctionTypeMap();

    /// return shape function type enum for a given shape function name
    ShapeFunctionType StringToShapeFunctionType(std::string name);

    /// return shape function name for a given shape function type
    std::string ShapeFunctionTypeToString(ShapeFunctionType shapefunctiontype);
  }  // namespace PROBLEMTYPE
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_PROBLEMTYPE_H
