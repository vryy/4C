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
      friction_none,    ///< no friction
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
      solution_vague,    ///< no solving strategy defined
      solution_lagmult,  ///< method of Lagrange multipliers
      solution_penalty,  ///< penalty approach
      solution_uzawa,    ///< Uzawa augmented Lagrange approach
      solution_augmented,///< augmented Lagrange approach
      solution_nitsche   ///< Nitsche contact solution method
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
      other                  ///< other problemtypes
    };

    /// weighting in Nitsche contact
    enum NitscheWeighting
    {
      NitWgt_slave,
      NitWgt_master,
      NitWgt_harmonic,
    };

    /// set the contact parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
