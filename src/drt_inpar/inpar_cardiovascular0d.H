/*----------------------------------------------------------------------*/
/*!
\file inpar_cardiovascular0d.H

\brief Input parameters for 0d cardiovascular-structure coupling

\level 2

\maintainer Marc Hirschvogel
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_CARDIOVASCULAR0D_H_
#define INPAR_CARDIOVASCULAR0D_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

namespace INPAR
{

namespace CARDIOVASCULAR0D
{

  /// possible 0D cardiovascular-structural solvers
  enum Cardvasc0DSolveAlgo
  {
    cardvasc0dsolve_direct, ///< build monolithic 0D cardiovascular-structural system
    cardvasc0dsolve_simple  ///< use simple preconditioner for iterative solve
  };

  enum Cardvasc0DAtrimModel
  {
    elastance_0d,
    structure_3d
  };

  /// set the 0Dcardiovascular parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  /// set specific 0Dcardiovascular conditions
  void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

}// namespace CARDIOVASCULAR0D
}// namespace INPAR
/*----------------------------------------------------------------------*/
#endif /* INPAR_CARDIOVASCULAR0D_H_ */
