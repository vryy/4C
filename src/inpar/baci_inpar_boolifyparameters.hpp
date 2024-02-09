/*----------------------------------------------------------------------*/
/*! \file

\brief Auxiliar routine to boolify integral Yes/No data


\level 0
*/

#ifndef BACI_INPAR_BOOLIFYPARAMETERS_HPP
#define BACI_INPAR_BOOLIFYPARAMETERS_HPP

#include "baci_config.hpp"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN


namespace INPUT
{
  /// Auxiliar routine to boolify integral Yes/No data
  ///
  /// Parameters consisting integral values for Yes/No tuples
  /// are removed and replaced by a bool parameter holding true/false
  /// respectively.
  void BoolifyValidInputParameters(Teuchos::ParameterList& list  ///< the valid input parameter list
  );

}  // namespace INPUT

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // INPAR_BOOLIFYPARAMETERS_H
