/*
 * inpar_loca_continuation.H
 *
 *  Created on: Nov 20, 2015
 *      Author: hiermeier
 */

#ifndef SRC_DRT_INPAR_INPAR_LOCA_CONTINUATION_H_
#define SRC_DRT_INPAR_INPAR_LOCA_CONTINUATION_H_

#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

namespace INPAR
{
  namespace LOCA
  {
    //! Set valid loca continuation parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    void SetValidPredictorParams(Teuchos::ParameterList& predictor);

    void SetValidBorderedSolverMethod(Teuchos::ParameterList& list);
  }  // namespace LOCA
}  // namespace INPAR


#endif /* SRC_DRT_INPAR_INPAR_LOCA_CONTINUATION_H_ */
