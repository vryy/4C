/*----------------------------------------------------------------------*/
/*!
\file inpar_structure.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef CCADISCRET
#ifndef INPAR_STRUCTURE_H
#define INPAR_STRUCTURE_H

#ifdef __cplusplus /* stuff which is c++ specific */

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// General input parameters for structural problems
  namespace STR
  {

    /// @name Time integration
    //@{

    /// Type of time integrator including statics
    enum DynamicType
    {
      dyna_gen_alfa,                   ///< generalised-alpha time integrator
      dyna_gen_alfa_statics,           ///< static analysis with gen_alfa
      dyna_centr_diff,                 ///< central differences (explicit)
      dyna_Gen_EMM,                    ///< generalised energy-momentum method
      dyna_statics,                    ///< static analysis
      dyna_genalpha,                   ///< generalised-alpha time integrator, 
                                       ///< new style (implicit)
      dyna_onesteptheta,               ///< one-step-theta time integrator (implicit)
      dyna_gemm,                       ///< generalised energy-momementum method
      dyna_ab2,                        ///< Adams-Bashforth 2nd order (explicit)
      dyna_euma,                       ///< Euler-Maruyama scheme (stochastic, explicit)
      dyna_euimsto                     ///< Euler implicit scheme for SDE's (stochastic, implicit)
    };

    /// Map time integrator to string
    static inline const std::string DynamicTypeString
    (
      const enum DynamicType name  ///< enum to convert
    )
    {
      switch (name)
      {
      case dyna_gen_alfa :
        return "Gen_Alfa";
        break;
      case dyna_gen_alfa_statics :
        return "Static";
      case dyna_centr_diff :
        return "Centr_Diff";
        break;
      case dyna_Gen_EMM :
        return "Gen_EMM";
        break;
      case dyna_statics :
        return "Statics";
        break;
      case dyna_genalpha :
        return "GenAlpha";
        break;
      case dyna_onesteptheta :
        return "OneStepTheta";
        break;
      case dyna_gemm :
        return "GEMM";
        break;
      case dyna_ab2 :
        return "AdamsBashforth2" ;
        break;
      default :
        dserror("Cannot make string for time integrator %d", name);
        return "";
      }
    }

    /// Type of (global) damping
    enum DampKind
    {
      damp_none = 0,            ///< damping off
      damp_rayleigh,            ///< globally applied Rayleigh damping
      damp_material,            ///< element-wise applied damping
                                ///< using element velocities
      damp_brownianmotion       ///< damping due to brownian motion
    };

    /// Mid-average type of internal forces for generalised-alpha-like
    /// time integration schemes
    ///
    /// \author bborn
    /// \date 11/08
    enum MidAverageEnum
    {
      midavg_vague = 0,  ///< undefined mid-averaging type
      midavg_imrlike,  ///< alphaf-mid-averaging is done IMR-like, i.e. 
                       ///< \f$F_{int,m}\f$ 
                       ///< \f$= F_{int}(D_m)\f$
                       ///< \f$= F_{int}( (1-\alpha_f)*D_{n+1} + \alpha_f*D_n )\f$
                       ///< (IMR means implicit mid-point rule.)
      midavg_trlike    ///< alphaf-mid-averaging is done TR-like, i.e.
                       ///< \f$F_{int,m}\f$
                       ///< \f$= (1-\alpha_f)*F_{int,n+1} + \alpha_f*F_{int,n}\f$
                       ///<  \f$= (1-\alpha_f)*F_{int}(D_{n+1}) + \alpha_f*F_{int}(D_n)\f$
                       ///<  (TR means trapezoidal rule.)
    };

    //@}

    /// @name Solution technique and related
    //@{

    /// type of solution techniques
    enum NonlinSolTech
    {
      soltech_vague,  ///< undefined
      soltech_newtonfull,  ///< full Newton-Raphson iteration
      soltech_newtonls,  ///< line search Newton-Raphson
      soltech_newtonmod,  ///< modified Newton-Raphson iteration
      soltech_newtonuzawalin,  ///< linear Uzawa iteration for 
                               ///< constraint system
      soltech_newtonuzawanonlin,  ///< non-linear Uzawa iteration 
                                  ///< for constraint system
      soltech_nlncg,  ///< nonlinear CG iteration using Nox
      soltech_ptc,  ///< pseudo transient continuation nonlinear iteration
      soltech_noxnewtonlinesearch,  ///< Line search Newton
                                    ///< utilising NOX
      soltech_noxgeneral  ///< non-linear solution with NOX
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
      case soltech_newtonls :
        return "lsnewton";
        break;
      case soltech_newtonmod :
        return "modnewton";
        break;
      case soltech_newtonuzawalin :
        return "newtonlinuzawa";
        break;
      case soltech_newtonuzawanonlin :
        return "augmentedlagrange";
        break;
      case soltech_noxnewtonlinesearch :
        return "NoxNewtonLineSearch";
        break;
      case soltech_noxgeneral :
        return "noxgeneral";
        break;
      case soltech_nlncg :
        return "nlncg";
        break;
      case soltech_ptc :
        return "ptc";
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
      pred_constdis,  ///< constant displacements
      pred_constdisvelacc,  ///< constant displacements, 
                           ///< velocities and accelerations
      pred_tangdis  ///< linearised solution obeying DBC displacements via tangent
                    ///< D_{n+1}^{<0>} = D_{n} + Ktang_{n,eff}^{-1} . (- Ktang_{n} . (D_{n+1}^{DBC} - D_{n}))
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
      case pred_constdis :
        return "ConstDis";
        break;
      case pred_constdisvelacc :
        return "ConstDisVelAcc";
        break;
      case pred_tangdis :
        return "TangDis";
        break;
      default :
        dserror("Cannot make string for predictor %d", name);
        return "";
      }
    }

    /// kind of iteration convergence checks
    /// \author lw \date 12/07
    enum ConvCheck
    {
      convcheck_absres_or_absdis,  ///< absolute norms of residual forces
                                   ///< OR iterative displacement increments
      convcheck_absres_and_absdis,  ///< absolute norms of residual forces
                                    ///< AND iterative displacement increments
      convcheck_relres_or_absdis,  ///< relative norm of residual forces
                                   ///< OR absolute norm if iterative
                                   ///< displacement increments
      convcheck_relres_and_absdis,  ///< relative norm of residual forces
                                    ///< AND absolute norm if iterative
                                    ///< displacement increments
      convcheck_relres_or_reldis,  ///< relative norms of residual forces
                                   ///< OR iterative displacement increments
      convcheck_relres_and_reldis,  ///< relative norms of residual forces
                                    ///< AND iterative displacement increments
      convcheck_mixres_or_mixdis,  ///< mixed absolute-relative norms of residual forces
                                    ///< OR mixed absolute-relative norm of residual displacements
      convcheck_mixres_and_mixdis,  ///< mixed absolute-relative norms of residual forces
                                    ///< AND mixed absolute-relative norm of residual displacements
      convcheck_vague  ///< undetermined
    };

    /// map convergence check to enum term
    static inline const std::string ConvCheckString
    (
      const enum ConvCheck name  ///< enum term
    )
    {
      switch (name)
      {
      case convcheck_absres_or_absdis :
        return "AbsRes_Or_AbsDis";
        break;
      case convcheck_absres_and_absdis :
        return "AbsRes_And_AbsDis";
        break;
      case convcheck_relres_or_absdis :
        return "RelRes_Or_AbsDis";
        break;
      case convcheck_relres_and_absdis :
        return "RelRes_And_AbsDis";
        break;
      case convcheck_relres_or_reldis :
        return "RelRes_Or_RelDis";
        break;
      case convcheck_relres_and_reldis :
        return "RelRes_And_RelDis";
        break;
      case convcheck_mixres_or_mixdis :
        return "MixRes_Or_MixDis";
        break;
      case convcheck_mixres_and_mixdis :
        return "MixRes_And_MixDis";
        break;
      default :
        dserror("Cannot make string for convergence check %d", name);
        return "";
      }
    }

    //@}

    /// @name Constraints (global, geometric)
    //@{

    /// possible constraint solvers
    enum ConSolveAlgo
    {
      consolve_direct,  ///< ???
      consolve_iterative  ///< ???
    };

    //@}

    /// @name Output
    //@{

    /// Type of structural stress output
    /// (this enum represents the input file parameter STRUCT_STRESS)
    enum StressType
    {
      stress_none,    ///< no stress output
      stress_cauchy,  ///< output of Cauchy stresses
      stress_2pk      ///< output of 2nd Piola-Kirchhoff stresses
    };

    /// Type of structural strain output
    /// (this enum represents the input file parameter STRUCT_STRAIN)
    enum StrainType
    {
      strain_none,    ///< no strain output
      strain_ea,      ///< output of Euler-Almansi strains
      strain_gl       ///< output of Green-Lagrange strains
    };

    //@}

    /// @name Time adapivity
    //@{

    /// type of adaptivity in time
    enum TimAdaKind
    {
      timada_kind_none,           ///< no time adaptivity
      timada_kind_zienxie,        ///< Zienkiewicz-Xie indicator
      timada_kind_ab2             ///< Adams-Bahsforth2 indicator
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
      case INPAR::STR::norm_vague:
        return "Vague";
        break;
      case INPAR::STR::norm_l1:
        return "L1";
        break;
      case INPAR::STR::norm_l2:
        return "L2";
        break;
      case INPAR::STR::norm_rms:
        return "Rms";
        break;
      case INPAR::STR::norm_inf:
        return "Inf";
        break;
      default:
        dserror("Cannot make string to vector norm %d", norm);
        return "";
      }
    }

    //@}

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifdef __cplusplus
#endif  // #ifndef INPAR_STRUCTURE_H
#endif  // #ifdef CCADISCRET
