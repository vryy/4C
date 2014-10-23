/*!------------------------------------------------------------------------------------------------*
\file inpar_poroelast.H

\brief

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef INPAR_POROELAST_H_
#define INPAR_POROELAST_H_

#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_dserror.H"

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
    Monolithic_fluidsplit,
    Monolithic_nopenetrationsplit
  };

  /// type of time discretization for continuity equation
  enum TimeDisTypeConti
  {
//    OneWay,
 //   SequStagg,
 //   IterStagg,
    porosity,
    pressure
  };

  /// @name Solution technique and related

  /// type of norm to check for convergence
  enum ConvNorm
  {
    convnorm_abs_global,       ///< absolute norm of global solution vectors
    convnorm_abs_singlefields ///< absolute norm of single field solution vectors
//    convnorm_rel_global,       ///< absolute norm of global solution vectors
//    convnorm_rel_singlefields  ///< absolute norm of single field solution vectors
  };

  /// type of norm to be calculated
  enum VectorNorm
  {
    norm_l1,  //!< L1/linear norm
    norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
    norm_l2,  //!< L2/Euclidean norm
    norm_rms,  //!< root mean square (RMS) norm
    norm_inf  //!< Maximum/infinity norm
  };

  /// type of norm to check for convergence
  enum BinaryOp
  {
    bop_and,  ///<  and
    bop_or    ///<  or
  };

  /// type of initial field for poroelasticity problem
  enum InitialField
  {
  //  initfield_zero_field,
    initfield_field_by_function
 //   initfield_field_by_condition
  };

  //! map enum term to std::string
  static inline const std::string VectorNormString(
    const enum VectorNorm norm  //!< input enum term
    )
  {
    switch (norm)
    {
    case INPAR::POROELAST::norm_l1 :
      return "L1";
      break;
    case INPAR::POROELAST::norm_l1_scaled :
      return "L1_scaled";
      break;
    case INPAR::POROELAST::norm_l2 :
      return "L2";
      break;
    case INPAR::POROELAST::norm_rms :
      return "Rms";
      break;
    case INPAR::POROELAST::norm_inf :
      return "Inf";
      break;
    default:
      dserror("Cannot make std::string to vector norm %d", norm);
      return "";
    }
  }

  //@}


}  // namespace POROELAST

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_POROELAST_H_ */


