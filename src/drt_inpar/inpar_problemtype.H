/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
\maintainer Martin Kronbichler
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
    /*! \brief Type of shape functions used in spatial discretitzation
     *
     * Assigning explicit integer values here for historic reasons. This should go away, once we
     * consequently use this enum throughout the entire code.
     */
    enum ShapeFunctionType
    {
      shapefunction_undefined = -1,  ///< Undefined
      shapefunction_polynomial = 1,  ///< Polynomial shape functions
      shapefunction_nurbs = 0,       ///< NURBS shape functions
      shapefunction_meshfree = 2,    ///< Meshfree
      shapefunction_hdg = 3          ///< Hybridizable Discontinuous Galerkin
    };

    /*! \brief Define valid parameters
     *
     * @param[in/out] list Parameter list to be filled with valid parameters and their defaults
     */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// create map of problem name and problem type enum
    std::map<std::string, PROBLEM_TYP> StringToProblemTypeMap();

    /// return problem type enum for a given problem name
    PROBLEM_TYP StringToProblemType(std::string name);

  }  // namespace PROBLEMTYPE
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_PROBLEMTYPE_H
