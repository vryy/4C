/*!------------------------------------------------------------------------------------------------*
\file inpar_poroelast.H

\brief

<pre>
Maintainer: Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef INPAR_POROELAST_H_
#define INPAR_POROELAST_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for poroelasticity
namespace POROELAST
{

  /// Type of coupling strategy for poroelasticity problems
  enum SolutionSchemeOverFields
  {
//    OneWay,
 //   SequStagg,
 //   IterStagg,
    Partitioned,
    Monolithic,
    Monolithic_structuresplit,
    Monolithic_fluidsplit
  };

  /// @name Solution technique and related

  /// type of norm to check for convergence
  enum ConvNorm
  {
    convnorm_abs  ///< absolute norm
  };

  /// type of norm to check for convergence
  enum BinaryOp
  {
    bop_and,  ///<  and
    bop_or	  ///<  or
  };

  /// type of initial field for poroelasticity problem
  enum InitialField
  {
  //  initfield_zero_field,
    initfield_field_by_function
 //   initfield_field_by_condition
  };

  //@}


}  // namespace POROELAST

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_POROELAST_H_ */


