/*----------------------------------------------------------------------*/
/*!
\file inpar_thermo.H

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef CCADISCRET
#ifndef INPAR_THERMO_H
#define INPAR_THERMO_H

#ifdef __cplusplus /* stuff which is c++ specific */

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// General input parameters for thermal problems
  namespace THR
  {

    /// @name Time integration
    //@{

    /// Type of time integrator including statics
    enum DynamicType
    {
      dyna_undefined,                  ///< undefined integrator (sth like a default)
      dyna_statics,                    ///< static analysis
      dyna_onesteptheta,               ///< one-step-theta time integrator (implicit)
      dyna_gemm                        ///< generalised energy-momementum method
    };

    /// Map time integrator to string
    static inline const std::string DynamicTypeString
    (
      const enum DynamicType name  ///< enum to convert
    )
    {
      switch (name)
      {
      case dyna_undefined :
        return "Undefined";
        break;
      case dyna_statics :
        return "Statics";
        break;
      case dyna_onesteptheta :
        return "OneStepTheta";
        break;
      case dyna_gemm :
        return "GEMM";
        break;
      default :
        dserror("Cannot make string for time integrator %d", name);
        return "";
      }
    }

    //@}

    /// @name Solution technique and related
    //@{

    /// type of solution techniques
    enum NonlinSolTech
    {
      soltech_vague,  ///< undefined
      soltech_newtonfull  ///< full Newton-Raphson iteration
    };

    /// Map solution technique enum to string
    static inline const std::string NonlinSolTechString
    (
      const enum NonlinSolTech name  ///< enum to convert
    )
    {
      switch (name)
      {
      case soltech_vague :
        return "vague";
        break;
      case soltech_newtonfull :
        return "fullnewton";
        break;
       default :
        dserror("Cannot make string for solution technique %d", name);
        return "";
      }
    }

    /// Type of predictor
    enum PredEnum
    {
      pred_vague,  ///< undetermined
      pred_consttemp,  ///< constant temperatures
      pred_consttemprate, ///< constant temperatures and rates
      pred_tangtemp  ///< linearised solution obeying DBC temperature via tangent
                    ///< T_{n+1}^{<0>} = T_{n} + Ktang_{n,eff}^{-1} . (- Ktang_{n} . (T_{n+1}^{DBC} - T_{n}))
                    ///< This looks hilarious, but remember Ktan_{n,eff}^{-1} is not
                    ///< the inverse of Ktan_{n} due to the application of the
                    ///< Dirichlet BCs (i.e. the reduction to the test space).
    };

    /// Map predictor enum term to string
    static inline const std::string PredEnumString
    (
      const PredEnum name  ///< identifier
    )
    {
      switch (name)
      {
      case pred_vague :
        return "Vague";
        break;
      case pred_consttemp :
        return "ConstTemp";
        break;
      case pred_consttemprate :
              return "ConstTempRate";
              break;
      case pred_tangtemp :
        return "TangTemp";
        break;
      default :
        dserror("Cannot make string for predictor %d", name);
        return "";
      }
    }

     /// type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_abs,         ///< absolute norm
      convnorm_rel,         ///< relative norm
      convnorm_mix         ///< mixed absolute-relative norm
    };

    /// type of norm to check for convergence
    enum BinaryOp
    {
      bop_or,         ///<  or
      bop_and         ///<  and
    };

    //@}

     /// @name Output
    //@{

    /// Type of thermal flux output
    /// (this enum represents the input file parameter THERM_HEATFLUX) CHECK IT!
    enum HeatFluxType
    {
      heatflux_none,    ///< no heatflux output
      heatflux_current,  ///< output of heatflux in current configuration
      heatflux_initial      ///< output of heat flux in initial configuration
    };

    /// Type of thermal gradient output
     /// (this enum represents the input file parameter THERM_TEMPGRAD) CHECK IT!
     enum TempGradType
     {
       tempgrad_none,    ///< no thermal gradient output
       tempgrad_current,  ///< output of thermal gradient in current configuration
       tempgrad_initial      ///< output of thermal gradient in initial configuration
     };

     //@}

     /// @name General
    //@{

    /// type of vector norm used for error/residual vectors
    enum VectorNorm
    {
      norm_vague = 0,  ///< undetermined norm
      norm_l1,         ///< L1/linear norm
      norm_l2,         ///< L2/Euclidean norm
      norm_rms,        ///< root mean square (RMS) norm
      norm_inf         ///< Maximum/infinity norm
    };

    /// map enum term to string
    static inline const std::string VectorNormString
    (
      const enum VectorNorm norm  ///< input enum term
    )
    {
      switch (norm)
      {
      case INPAR::THR::norm_vague:
        return "Vague";
        break;
      case INPAR::THR::norm_l1:
        return "L1";
        break;
      case INPAR::THR::norm_l2:
        return "L2";
        break;
      case INPAR::THR::norm_rms:
        return "Rms";
        break;
      case INPAR::THR::norm_inf:
        return "Inf";
        break;
      default:
        dserror("Cannot make string to vector norm %d", norm);
        return "";
      }
    }

    //@}

  }  // namespace THR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifdef __cplusplus
#endif  // #ifndef INPAR_THERMO_H
#endif  // #ifdef CCADISCRET
