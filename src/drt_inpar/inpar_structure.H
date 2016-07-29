/*----------------------------------------------------------------------*/
/*!
\file inpar_structure.H

\brief structural dynamics input parameters

\level 1

<pre>
\maintainer Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef INPAR_STRUCTURE_H
#define INPAR_STRUCTURE_H


#include "inpar_parameterlist_utils.H"
#include <boost/algorithm/string/predicate.hpp>

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// General input parameters for structural problems
  namespace STR
  {

    /// @name Time integration
    //@{
    /// Active element technologies
    enum EleTech
    {
      eletech_eas,        ///< EAS detected
      eletech_fbar,       ///< FBAR detected
      eletech_plasticity, ///< plasticity detected
      eletech_pressure    ///< additional pressure degree of freedom
    };

    /// Map element technology to string
    static inline const std::string EleTechString
    (
      const enum EleTech name  ///< enum to convert
    )
    {
      switch (name)
      {
        case eletech_eas:
          return "EAS";
          break;
        case eletech_fbar:
          return "FBar";
          break;
        case eletech_plasticity:
          return "plasticity";
          break;
        case eletech_pressure:
          return "pressure";
          break;
        default :
          dserror("Cannot make std::string for element technology %d", name);
          break;
      }
      return "";
    }

    /// Type of the used models/conditions
    /// necessary for the model evaluator
    enum ModelType
    {
      model_structure,          ///< evaluate the structural model
      model_contact,            ///< evaluate the contact model
      model_meshtying,          ///< evaluate the meshtying model
      model_cardiovascular0d,   ///< evaluate the cardiovascular0d model
      model_springdashpot,      ///< evaluate the springdashpot model
      model_lag_pen_constraint, ///< evaluate the Lagrange or/and penalty enforced constraint models
      model_vague               ///< undefined model type
    };

    /// Map model type to string
    static inline const std::string ModelTypeString
    (
      const enum ModelType name  ///< enum to convert
    )
    {
      switch (name)
      {
      case model_structure :
        return "Structure";
        break;
      case model_contact :
        return "Contact";
        break;
      case model_meshtying :
        return "Meshtying";
        break;
      case model_cardiovascular0d :
        return "Cardiovascular0D";
        break;
      case model_springdashpot :
        return "SpringDashpot";
        break;
      case model_lag_pen_constraint :
        return "Lag-Pen-Constraint";
        break;
      default :
        dserror("Cannot make std::string for model type %d", name);
        return "";
      }
    };

    //! Map model std::string to enum
    inline enum ModelType String2ModelType(const std::string& name)
    {
      ModelType type = model_vague;
      if (boost::iequals(name,"structure"))
        type = model_structure;
      else if (boost::iequals(name,"contact"))
        type = model_contact;
      else if (boost::iequals(name,"meshtying"))
        type = model_meshtying;
      else if (boost::iequals(name,"cardiovascular0d"))
        type = model_cardiovascular0d;
      else if (boost::iequals(name,"springdashpot"))
        type = model_springdashpot;
      else if (boost::iequals(name,"lag-pen-constraint"))
        type = model_lag_pen_constraint;

      return type;
    };

    enum IntegrationStrategy
    {
      int_old,            // old structure (FixMe deprecated and should be deleted, as soon as the clean-up is finished!)
      int_standard,       // standard time integration (implicit or explicit)
      int_loca            // path following with LOCA
    };

    /// Type of time integrator including statics
    enum DynamicType
    {
      dyna_statics,                    ///< static analysis
      dyna_genalpha,                   ///< generalised-alpha time integrator (implicit)
      dyna_onesteptheta,               ///< one-step-theta time integrator (implicit)
      dyna_gemm,                       ///< generalised energy-momementum method (implicit)
      dyna_expleuler,                  ///< forward Euler (explicit)
      dyna_centrdiff,                  ///< central differences (explicit)
      dyna_ab2,                        ///< Adams-Bashforth 2nd order (explicit)
      dyna_euma,                       ///< Euler-Maruyama scheme (stochastic, explicit)
      dyna_euimsto,                    ///< Euler implicit scheme for SDE's (stochastic, implicit)
      dyna_statmech                    ///< Statistical Mechanics Time Integration
    };

    /// Map time integrator to std::string
    static inline const std::string DynamicTypeString
    (
      const enum DynamicType name  ///< enum to convert
    )
    {
      switch (name)
      {
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
      case dyna_expleuler :
        return "ExplEuler";
        break;
      case dyna_centrdiff :
        return "CentrDiff";
        break;
      case dyna_ab2 :
        return "AdamsBashforth2" ;
        break;
      case dyna_statmech :
        return "StatMech" ;
        break;
      default :
        dserror("Cannot make std::string for time integrator %d", name);
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

    /// Map mid-averaging to std::string
    static inline const std::string MidAverageString
    (
      const enum MidAverageEnum name  ///< enum to convert
    )
    {
      switch (name)
      {
      case midavg_vague :
        return "Vague";
        break;
      case midavg_imrlike :
        return "IMR-like";
      case midavg_trlike :
        return "TR-like";
        break;
      default :
        dserror("Cannot make std::string for time integrator %d", name);
        return "";
      }
    }

    //@}

    /// @name Solution technique and related
    //@{

    /// type of control
    enum ControlType
    {
      control_load,  ///< load control
      control_disp,  ///< disp control
      control_arc1,  ///< arclength RIKS control
      control_arc2   ///< arclength Crisfield control
    };

    /// have inertia forces to be linearized?
    enum MassLin
    {
      ml_none,      ///< constant mass matrix
      ml_standard,  ///< nonlinear inertia terms, translational DoFs
      ml_rotations  ///< nonlinear inertia terms, rotational DoFs
    };

    /// type of solution techniques
    enum NonlinSolTech
    {
      soltech_vague,              ///< undefined
      soltech_newtonfull,         ///< full Newton-Raphson iteration
      soltech_newtonls,           ///< line search Newton-Raphson
      soltech_newtonopp,          ///< oppositely converging Newton-Raphson
      soltech_newtonmod,          ///< modified Newton-Raphson iteration
      soltech_newtonuzawalin,     ///< linear Uzawa iteration for
                                  ///< constraint system
      soltech_newtonuzawanonlin,  ///< non-linear Uzawa iteration
                                  ///< for constraint system
      soltech_ptc,                ///< pseudo transient continuation nonlinear iteration
      soltech_noxnewtonlinesearch,///< Line search Newton
                                  ///< utilizing NOX
      soltech_noxgeneral,         ///< non-linear solution with NOX
      soltech_nox_nln,            ///< non-linear solution with NOX (new)
      soltech_nlnsol              ///< nonlinear solver toolbox
    };

    /// Map solution technique enum to std::string
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
      case soltech_newtonopp :
         return "oppnewton";
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
      case soltech_nox_nln :
        return "nox_nln";
        break;
      case soltech_ptc :
        return "ptc";
        break;
      case soltech_nlnsol :
        return "NLNSOL";
        break;
      default :
        dserror("Cannot make std::string for solution technique %d", name);
        return "";
      }
    }

    /// Handling of non-converged nonlinear solver
    enum DivContAct
    {
      divcont_stop,               ///< abort simulation
      divcont_continue,           ///< continue nevertheless
      divcont_repeat_step,        ///< repeat time step
      divcont_halve_step,         ///< halve time step and carry on with simulation
      divcont_adapt_step,         ///< adapt (halve or double) time step and carry on with simulation
      divcont_rand_adapt_step,    ///< adapt randomly time step and carry on with simulation
      divcont_rand_adapt_step_ele_err,    ///< adapt randomly time step and carry on with simulation, including acceptance of element errors in form of negative Jacobian determinant
      divcont_repeat_simulation,  ///< repeat the whole simulation
      divcont_adapt_penaltycontact,  ///< slightly adapt the penalty contact parameter if timestep doesn't converge
    };

    /// Map  enum to string
    static inline const std::string DivContActString
    (
      const enum DivContAct name  ///< enum to convert
    )
    {
      switch (name)
      {
      case divcont_stop :
        return "stop";
        break;
      case divcont_continue :
        return "continue";
        break;
      case  divcont_repeat_step :
        return "repeat_step";
        break;
      case  divcont_halve_step :
        return "halve_step";
        break;
      case  divcont_adapt_step :
        return "adapt_step";
        break;
      case  divcont_rand_adapt_step :
        return "rand_adapt_step";
        break;
      case  divcont_rand_adapt_step_ele_err :
        return "rand_adapt_step_ele_err";
        break;
      case divcont_repeat_simulation :
        return "repeat_simulation";
        break;
      case divcont_adapt_penaltycontact :
        return "adapt_penaltycontact";
        break;
      default :
        dserror("Cannot make string for solution div cont technique %d", name);
        return "";
      }
    }

    /// Handling of non-converged nonlinear solver
    enum ConvergenceStatus
    {
      conv_success = 0,         ///< converged successfully
      conv_nonlin_fail = 1,     ///< nonlinear solution procedure failed
      conv_lin_fail = 2,        ///< linear system failed
      conv_ele_fail = 3,        ///< failure in element in form of negative Jac. det.
      conv_fail_repeat = 4      ///< nonlinear solver failed, repeat step according to divercont action set in input file
    };


    /// Type of predictor
    enum PredEnum
    {
      pred_vague,  ///< undetermined
      pred_constdis,  ///< constant displacements, consistent velocities and accelerations
      pred_constvel,  ///< constant velocity, extrapolated displacements, consistent accelerations
      pred_constacc,  ///< constant acceleration, extrapolated displacements and velocities
      pred_constdisvelacc,  ///< constant displacements,
                           ///< velocities and accelerations
      pred_tangdis, ///< linearised solution obeying DBC displacements via tangent
                    ///< D_{n+1}^{<0>} = D_{n} + Ktang_{n,eff}^{-1} . (- Ktang_{n} . (D_{n+1}^{DBC} - D_{n}))
                    ///< This looks hilarious, but remember Ktan_{n,eff}^{-1} is not
                    ///< the inverse of Ktan_{n} due to the application of the
                    ///< Dirichlet BCs (i.e. the reduction to the test space).
      pred_constdispres,  ///< constant displacements and pressure
      pred_constdisvelaccpres  ///< constant displacements,
                           ///< velocities and accelerations
                           ///< and pressures
    };

    /// Map predictor enum term to std::string
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
      case pred_constvel :
        return "ConstVel";
        break;
      case pred_constacc :
        return "ConstAcc";
        break;
      case pred_constdisvelacc :
        return "ConstDisVelAcc";
        break;
      case pred_tangdis :
        return "TangDis";
        break;
      case pred_constdispres :
        return "ConstDisPres";
        break;
      case pred_constdisvelaccpres :
        return "ConstDisVelAccPres";
        break;
      default :
        dserror("Cannot make std::string for predictor %d", name);
        return "";
      }
    }

    /// type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_abs,   ///< absolute norm
      convnorm_rel,   ///< relative norm
      convnorm_mix    ///< mixed absolute-relative norm
    };

    /// type of norm to check for convergence
    enum BinaryOp
    {
      bop_or,         ///<  or
      bop_and         ///<  and
    };

    /// type of norm to check for convergence
    enum PreStress
    {
      prestress_none, ///<  none
      prestress_mulf, ///<  mulf
      prestress_id    ///<  id
    };

    /// STC scaling for thin shell structures
    enum STC_Scale
    {
      stc_none=0, ///< no scaling
      stc_curr,   ///< Non-symmetric STC
      stc_currsym ///< Symmetric STC
    };

    /// map convergence check to enum term
    static inline const std::string STCString
    (
      const enum STC_Scale name  ///< enum term
    )
    {
      switch (name)
      {
      case stc_none :
        return "stc_none";
        break;
      case stc_curr :
        return "stc_curr";
        break;
      case stc_currsym :
        return "stc_currsym";
        break;
      default :
        dserror("Cannot make std::string for stc method %d", name);
        return "";
      }
    }

    //@}

    /// @name Constraints (global, geometric)
    //@{

    /// possible constraint solvers
    enum ConSolveAlgo
    {
      consolve_direct, ///< build monolythic system for Lagrange multipliers
      consolve_uzawa,  ///< solve linear system iteratively (partitioned)
      consolve_simple  ///< use simple preconditioner for iterative solve
    };

    //@}

    /// @name Output
    //@{

    /// Type of structural stress output
    /// (this enum represents the input file parameter STRUCT_STRESS and
    /// STRUCT_COUPLING_STRESS)
    enum StressType
    {
      stress_none,    ///< no stress output
      stress_cauchy,  ///< output of Cauchy stresses
      stress_2pk      ///< output of 2nd Piola-Kirchhoff stresses
    };

    /// Type of structural strain output
    /// (this enum represents the input file parameter STRUCT_STRAIN and
    /// STRUCT_PLASTIC_STRAIN)
    enum StrainType
    {
      strain_none,    ///< no strain output
      strain_ea,      ///< output of Euler-Almansi strains
      strain_gl,      ///< output of Green-Lagrange strains
      strain_log      ///< output of Logarithmic (or Hencky) strains
    };

    //@}

    /// @name Time adapivity
    //@{

    /// type of adaptivity in time
    enum TimAdaKind
    {
      timada_kind_none, ///< no time adaptivity
      timada_kind_zienxie, ///< Zienkiewicz-Xie indicator
      timada_kind_ab2, ///< Adams-Bahsforth2 indicator
      timada_kind_expleuler, ///< Explicit Euler indicator
      timada_kind_centraldiff ///< Central difference indicator
    };

    //@}

    /// @name General
    //@{

    /// type of vector norm used for error/residual vectors
    enum VectorNorm
    {
      norm_vague = 0,  ///< undetermined norm
      norm_l1,         ///< L1/linear norm: \f$\vert x\vert_1 = \sum_{i=1}^N |x_i|\f$
      norm_l2,         ///< L2/Euclidean norm: \f$\vert x\vert_2 = \sum_{i=1}^N x_i^2\f$
      norm_rms,        ///< root mean square (RMS) norm: \f$\vert x\vert_{rms} = \left(\sum_{i=1}^N x_i^2\right) / \sqrt{N}\f$
      norm_inf         ///< Maximum/infinity norm: \f$\vert x\vert_\infty = \max{x_i}\f$
    };

    /// map enum term to std::string
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
        dserror("Cannot make std::string to vector norm %d", norm);
        return "";
      }
    }

    /// kinematic description
    enum KinemType
    {
      kinem_vague = 0,      ///< undetermined kinematics
      kinem_linear,         ///< linear kinematics
      kinem_nonlinearTotLag ///< nonlinear kinematics Total Lagrange
    };
    //@}

    /// set the structure parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set structure-specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STRUCTURE_H
