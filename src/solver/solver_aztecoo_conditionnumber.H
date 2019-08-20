/*-----------------------------------------------------------*/
/*! \file

\maintainer  Martin Kronbichler

\brief Condition number via AztecOO

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_SOLVER_AZTECOO_CONDITIONNUMBER_H_
#define SRC_SOLVER_SOLVER_AZTECOO_CONDITIONNUMBER_H_

#include <AztecOO_ConditionNumber.h>

namespace LINALG
{
  class AztecOOConditionNumber : public ::AztecOOConditionNumber
  {
   public:
    AztecOOConditionNumber() : ::AztecOOConditionNumber(){};
    virtual ~AztecOOConditionNumber() = default;

    /// \brief Access the aztec solver status
    /**
     *  See az_aztec_defs.h for a comprehensive overview of possible az_key
     *  values.
     *
     *  \author hiermeier \date 04/18 */
    double getStatus(int az_key) const;
  };
}  // namespace LINALG


#endif /* SRC_SOLVER_SOLVER_AZTECOO_CONDITIONNUMBER_H_ */
