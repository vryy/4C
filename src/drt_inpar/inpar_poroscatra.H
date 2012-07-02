/*
 * inpar_poroscatra.H
 *
 *  Created on: Mar 9, 2012
 *      Author: urrecha
 */
#ifndef INPAR_POROSCATRA_H_
#define INPAR_POROSCATRA_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for poroelasticity
namespace PORO_SCATRA
{

  /// Type of coupling strategy for poro scatra problems
  enum SolutionSchemeOverFields
  {
    OneWay
 //   SequStagg,
 //   IterStagg,
 //   Monolithic
  };


}  // namespace PORO_SCATRA

}  // namespace INPAR

/*----------------------------------------------------------------------*/



#endif /* INPAR_POROSCATRA_H_ */
