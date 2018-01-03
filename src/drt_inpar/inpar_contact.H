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
          return "UNSUPPORTED SolvingStrategy";
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
          return "UNSUPPORTED assemble strategy";
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
          return "UNSUPPORTED variational approach";
      }
    }

    /// set the contact parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
