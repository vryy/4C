/*----------------------------------------------------------------------*/
/*!
\file inpar_contact.H

\brief Input parameters for contact

\level 2

\maintainer Alexander Popp

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_CONTACT_H
#define INPAR_CONTACT_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for structural CONTACT and meshtying
  namespace CONTACT
  {
    /// Type of contact friction law
    /// (this enum represents the input file parameter FRICTION)
    enum FrictionType
    {
      friction_none=0,    ///< no friction
      friction_stick,   ///< perfect stick
      friction_tresca,  ///< Tresca friction law
      friction_coulomb  ///< Coulomb friction law
    };

    /// Type of contact adhesion law
    /// (this enum represents the input file parameter ADHESION)
    enum AdhesionType
    {
      adhesion_none,    ///< no adhesion
      adhesion_bound    ///< fix bounded adhesion
    };

    /// Type of employed solving strategy
    /// (this enum represents the input file parameter STRATEGY)
    enum SolvingStrategy
    {
      solution_vague,           ///< no solving strategy defined
      solution_lagmult,         ///< method of Lagrange multipliers
      solution_penalty,         ///< penalty approach
      solution_uzawa,           ///< Uzawa augmented Lagrange approach
      solution_combo,           ///< combination of multiple augmented solving strategies
      solution_augmented,       ///< augmented Lagrange approach
      solution_steepest_ascent, ///< steepest ascent approach (based on the augmented formulation)
      solution_std_lagrange,    ///< lagrange strategy (based on the augmented formulation)
      solution_xcontact,        ///< eXtended contact approach
      solution_nitsche,         ///< Nitsche contact solution method
      solution_ehl              ///< method for elasto-hydrodynamic lubrication
    };

    inline std::string SolvingStrategy2String( enum SolvingStrategy stype )
    {

      switch ( stype )
      {
        case solution_vague:
          return "solution_vague";
        case solution_lagmult:
          return "solution_lagmult";
        case solution_penalty:
          return "solution_penalty";
        case solution_uzawa:
          return "solution_uzawa";
        case solution_combo:
          return "solution_combo";
        case solution_augmented:
          return "solution_augmented";
        case solution_steepest_ascent:
          return "solution_steepest_ascent";
        case solution_std_lagrange:
          return "solution_std_lagrange";
        case solution_xcontact:
          return "solution_xcontact";
        case solution_nitsche:
          return "solution_nitsche";
        default:
          return "INVALID SolvingStrategy";
      }

    }

    enum SwitchingStrategy
    {
      switch_preasymptotic    /**< switch between two stratgies.
                               *   One for the pre- and one for the asymptotic
                               *   solution phase */
    };

    /// Type of linear system setup and solution
    /// (this enum represents the input file parameter SYSTEM)
    enum SystemType
    {
      system_none,        ///< no system defined
      system_condensed,   ///< condensed system
      system_condensed_lagmult,   ///< system with condensed lagrange multiplier (differs just in meshtying case)
      system_saddlepoint  ///< saddle point system
    };

    /// Type of energy and momentum output
    /// (this enum represents the input file parameter EMOUTPUT)
    enum EmOutputType
    {
      output_none,    ///< no output
      output_screen,  ///< print to screen
      output_file,    ///< print to file
      output_both     ///< print to screen and file
    };

    /// Type of analytical solution for error norm computation
    /// (this enum represents the input file parameter ERROR_NORMS)
    /// (more details on analytical solutions, see drt_contact/contact_analytical.cpp)
    enum ErrorNorms
    {
      errornorms_none,       ///< no error norm computation
      errornorms_zero,       ///< error norms for zero analytical solution
      errornorms_bending,    ///< error norms for beam bending problem
      errornorms_sphere,     ///< error norms for pressurized sphere problem
      errornorms_thicksphere,///< error norms for thick pressurized sphere problem
      errornorms_infiniteplate ///< error norms for infinite plate with a circular hole
    };

    /// Type of formulation of constraint equations
    /// (this enum represents the input file parameter CONSTRAINT_DIRECTIONS)
    enum ConstraintDirection
    {
      constr_vague,          ///< no constraint directions defined
      constr_ntt,            ///< local normal and tangential coordinates
      constr_xyz             ///< global Cartesian coordinates
    };

    enum Regularization
    {
      reg_none,              ///< no regularization is applied
      reg_tanh               ///< regularization with tanh smoothing is applied
    };

    /// Local definition of problemtype to avoid use of drt_globalproblem.H
    enum Problemtype
    {
      structure,             ///< structural contact problem
      tsi,                   ///< coupled TSI problem with contact
      structalewear,         ///< wear problem including ALE shape changes
      poro,                  ///< poroelasticity problem with contact
      ehl,                   ///< elasto-hydrodymanic lubrication
      other                  ///< other problemtypes
    };

    /// weighting in Nitsche contact
    enum NitscheWeighting
    {
      NitWgt_slave,
      NitWgt_master,
      NitWgt_harmonic,
      NitWgt_phyiscal
    };

    /// Constraint enfrocement method method for thermal conduction and frictional dissipation
    enum NitscheThermoMethod
    {
      NitThr_substitution,
      NitThr_nitsche
    };

    /// Assemble strategy for the augmented Lagrangian framework
    enum AssembleStrategy : int
    {
      assemble_none,
      assemble_node_based   ///< assemble based on nodal data containers
    };

    /// convert assemble strategy to string
    inline std::string AssembleStrategy2String( enum AssembleStrategy assemble_type )
    {
      switch ( assemble_type )
      {
        case assemble_none:
          return "assemble_none";
        case assemble_node_based:
          return "assemble_node_based";
        default:
          return "INVALID assemble strategy";
      }
    }

    /// Variational approach for the augmented Lagrangian framework
    enum VariationalApproach : int
    {
      var_unknown,
      var_complete,   ///< complete variation
      var_incomplete  ///< incomplete variation
    };

    /// convert variational approach to string
    inline std::string VariationalApproach2String( enum VariationalApproach vartype )
    {
      switch ( vartype )
      {
        case var_unknown:
          return "var_unknown";
        case var_complete:
          return "var_complete";
        case var_incomplete:
          return "var_incomplete";
        default:
          return "INVALID variational approach";
      }
    }

    enum class FDCheck : char
    {
      off,
      global,
      gauss_point
    };

    /// Penalty update types
    enum class PenaltyUpdate : char
    {
      vague,                    ///< no update strategy defined (default)
      lm_gap_ratio,             ///< use the ratio between Lagrange multiplier and weighted gap
      complementarity,          ///< enforce the complementarity condition
      sufficient_lin_reduction, ///< enforce sufficient reduction of the linear infeasibility model
      sufficient_angle,         ///< enforce sufficient angle
      none                      ///< perform no penalty parameter update
    };

    /// Penalty update type to string
    inline std::string PenaltyUpdate2String( const enum PenaltyUpdate putype )
    {
      switch ( putype )
      {
        case PenaltyUpdate::vague:
          return "PenaltyUpdate::vague";
        case PenaltyUpdate::lm_gap_ratio:
          return "PenaltyUpdate::lm_gap_ratio";
        case PenaltyUpdate::complementarity:
          return "PenaltyUpdate::complementarity";
        case PenaltyUpdate::sufficient_lin_reduction:
          return "PenaltyUpdate::sufficient_lin_reduction";
        case PenaltyUpdate::sufficient_angle:
          return "PenaltyUpdate::sufficient_angle";
        case PenaltyUpdate::none:
          return "PenaltyUpdate::none";
        default:
          return "INVALID penalty update type";
      }
    }

    enum class PlotMode : int
    {
      off,
      write_single_iteration_of_step,
      write_each_iteration_of_step,
      write_last_iteration_of_step
    };

    enum class PlotType : int
    {
      vague = -2,
      scalar = -1,
      line  = 0,
      surface = 1,
      vector_field_2d = 2
    };

    inline std::string PlotType2String( const enum PlotType ptype )
    {
      switch ( ptype )
      {
        case PlotType::vague:
          return "vague";
        case PlotType::scalar:
          return "scalar";
        case PlotType::line:
          return "line";
        case PlotType::surface:
          return "surface";
        case PlotType::vector_field_2d:
          return "vectorfield2d";
        default:
          return "invalid";
      }
    }

    enum class PlotFuncName : char
    {
      vague,
      lagrangian,
      infeasibility,
      energy,
      energy_gradient,
      weighted_gap,
      weighted_gap_gradient,
      weighted_gap_mod_gradient,
      weighted_gap_gradient_error,
      weighted_gap_gradient_nodal_jacobian_error,
      weighted_gap_gradient_nodal_ma_proj_error
    };

    inline std::string PlotFuncName2String( const enum PlotFuncName pfunc )
    {
      switch ( pfunc )
      {
        case PlotFuncName::vague:
          return "func_vague";
        case PlotFuncName::lagrangian:
          return "func_lagrangian";
        case PlotFuncName::infeasibility:
          return "func_infeasibility";
        case PlotFuncName::energy:
          return "func_energy";
        case PlotFuncName::energy_gradient:
          return "func_energy_gradient";
        case PlotFuncName::weighted_gap:
          return "func_wgap";
        case PlotFuncName::weighted_gap_gradient:
          return "func_wgap_grad";
        case PlotFuncName::weighted_gap_mod_gradient:
          return "func_wgap_mod_grad";
        case PlotFuncName::weighted_gap_gradient_error:
          return "weighted_gap_gradient_error";
        case PlotFuncName::weighted_gap_gradient_nodal_jacobian_error:
          return "weighted_gap_gradient_nodal_jacobian_error";
        case PlotFuncName::weighted_gap_gradient_nodal_ma_proj_error:
          return "weighted_gap_gradient_nodal_ma_proj_error";
        default:
          return "INVALID plot func name";
      }
    }

    enum class PlotReferenceType : char
    {
      vague,
      current_solution,
      previous_solution
    };

    enum class PlotSupportType : char
    {
      vague,
      step_length,
      characteristic_element_length,
      position_angle
    };

    enum class PlotDirection : char
    {
      vague,
      current_search_direction,
      read_from_file,
      zero
    };

    /// direction split for the surface PlotType
    enum class PlotDirectionSplit : char
    {
      vague,                              ///< undefined split
      displacement_lagrange_multiplier,   ///< split into LM and displ. values
      slave_master_displacements          ///< split into slave/master displ.
    };

    enum class PlotFileFormat : char
    {
      vague,
      matlab,
      pgfplot
    };

    inline std::string PlotFileFormat2String( const enum PlotFileFormat pformat )
    {
      switch ( pformat )
      {
        case PlotFileFormat::matlab:
          return "mtl";
        case PlotFileFormat::pgfplot:
          return "pgf";
        default:
          return "invalid";
      }
    }

    /// set the contact parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
