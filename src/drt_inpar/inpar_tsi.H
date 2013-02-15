/*----------------------------------------------------------------------*/
/*!
\file inpar_tsi.H

\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 11/09 |
 *----------------------------------------------------------------------*/
#ifndef INPAR_TSI_H
#define INPAR_TSI_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for TSI
namespace TSI
{

  /// Type of coupling strategy for TSI problems
  enum SolutionSchemeOverFields
  {
    OneWay,
    SequStagg,
    IterStagg,
    IterStaggAitken,
    IterStaggAitkenIrons,
    IterStaggFixedRel,
    Monolithic
  };

  /// Type of coupling strategy for BACI-INCA coupling (TFSI)
  enum BaciIncaCoupling
  {
    undefinded,
    conforming,
    mortar_mortar_dual,
    mortar_mortar_std,
    proj_mortar_std,
    proj_LSI
  };

  /// @name Solution technique and related

  /// type of norm to check for convergence
  enum ConvNorm
  {
    convnorm_abs,  ///< absolute norm
    convnorm_rel,  ///< scaled norm of the single field residual with the length of the vector
    convnorm_reliter0,  ///< relative norm of TSI problem with inital TSI rhs
    convnorm_mix   ///< mixed absolute-relative norm
  };

  /// type of norm to check for convergence
  enum BinaryOp
  {
    bop_and  ///<  and
  };

  //@}


}  // namespace TSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_TSI_H
