/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for poro elasticity

\level 2

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
      Monolithic_nopenetrationsplit,
      Monolithic_meshtying
    };

    /// flag for control, in which equation the transient terms are included
    enum TransientEquationsOfPoroFluid
    {
      transient_none,
      transient_momentum_only,
      transient_continuity_only,
      transient_all
    };

    /// @name Solution technique and related

    /// type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_undefined,
      convnorm_abs_global,       ///< absolute norm of global solution vectors
      convnorm_abs_singlefields  ///< absolute norm of single field solution vectors
      //    convnorm_rel_global,       ///< absolute norm of global solution vectors
      //    convnorm_rel_singlefields  ///< absolute norm of single field solution vectors
    };

    /// type of norm to be calculated
    enum VectorNorm
    {
      norm_undefined,
      norm_l1,         //!< L1/linear norm
      norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
      norm_l2,         //!< L2/Euclidean norm
      norm_rms,        //!< root mean square (RMS) norm
      norm_inf         //!< Maximum/infinity norm
    };

    /// type of norm to check for convergence
    enum BinaryOp
    {
      bop_undefined,
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

    /// equilibration methods for global system of equations
    enum EquilibrationMethods
    {
      equilibration_none,
      equilibration_rows
    };

    //! map enum term to std::string
    static inline const std::string VectorNormString(
        const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case INPAR::POROELAST::norm_l1:
          return "L1";
          break;
        case INPAR::POROELAST::norm_l1_scaled:
          return "L1_scaled";
          break;
        case INPAR::POROELAST::norm_l2:
          return "L2";
          break;
        case INPAR::POROELAST::norm_rms:
          return "Rms";
          break;
        case INPAR::POROELAST::norm_inf:
          return "Inf";
          break;
        default:
          dserror("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }

    //@}


    /// set the poroelast parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace POROELAST

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_POROELAST_H_ */
