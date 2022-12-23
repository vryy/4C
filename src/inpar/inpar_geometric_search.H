/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for geometric search strategy

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef INPAR_GEOMETRIC_SEARCH_H
#define INPAR_GEOMETRIC_SEARCH_H

#include "inpar_parameterlist_utils.H"

namespace INPAR::GEOMETRICSEARCH
{
  //! set the parameters for the geometric search strategy
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace INPAR::GEOMETRICSEARCH

#endif
