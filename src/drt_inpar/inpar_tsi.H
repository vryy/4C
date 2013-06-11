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

//! Input parameters for TSI
namespace TSI
{

  //! Type of coupling strategy for TSI problems
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

  //! Type of coupling strategy for BACI-INCA coupling (TFSI)
  enum BaciIncaCoupling
  {
    undefinded,
    conforming,
    mortar_mortar_dual,
    mortar_mortar_std,
    proj_mortar_std,
    proj_mortar_dual,
    proj_RBFI
  };

  //! @name Solution technique and related stuff

  //! type of norm to check for convergence
  enum ConvNorm
  {
    convnorm_abs,  //!< absolute norm
    convnorm_rel,  //!< relative norm of TSI problem with inital TSI rhs
    convnorm_mix   //!< mixed absolute-relative norm
  };

  //! type of norm to check for convergence
  enum BinaryOp
  {
    bop_and,  //!< and
    bop_or,  //!< or
    bop_coupl_or_singl,  //!< either TSI problem or single field problems converged
    bop_coupl_and_singl,  //!< either TSI problem or single field problems converged
    bop_and_singl,  //!< and in single field problems
    bop_or_singl  //!< or in single field problems
  };

  //@}

  //! @name General
  //@{

  //! type of vector norm used for error/residual vectors
  enum VectorNorm
  {
    norm_vague = 0,  //!< undetermined norm
    norm_l1,  //!< L1/linear norm
    norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
    norm_l2,  //!< L2/Euclidean norm
    norm_rms,  //!< root mean square (RMS) norm
    norm_inf  //!< Maximum/infinity norm
  };

  //! map enum term to string
  static inline const std::string VectorNormString(
    const enum VectorNorm norm  //!< input enum term
    )
  {
    switch (norm)
    {
    case INPAR::TSI::norm_vague:
      return "Vague";
      break;
    case INPAR::TSI::norm_l1:
      return "L1";
      break;
    case INPAR::TSI::norm_l1_scaled:
      return "L1_scaled";
      break;
    case INPAR::TSI::norm_l2:
      return "L2";
      break;
    case INPAR::TSI::norm_rms:
      return "Rms";
      break;
    case INPAR::TSI::norm_inf:
      return "Inf";
      break;
    default:
      dserror("Cannot make string to vector norm %d", norm);
      return "";
    }
  }

  //@}


}  // namespace TSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_TSI_H
