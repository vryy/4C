/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for geometric search strategy

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_INPAR_GEOMETRIC_SEARCH_HPP
#define BACI_INPAR_GEOMETRIC_SEARCH_HPP

#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

namespace INPAR::GEOMETRICSEARCH
{
  //! set the parameters for the geometric search strategy
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace INPAR::GEOMETRICSEARCH

BACI_NAMESPACE_CLOSE

#endif
