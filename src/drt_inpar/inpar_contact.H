/*----------------------------------------------------------------------*/
/*!
\file inpar_contact.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_CONTACT_H
#define INPAR_CONTACT_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for structural CONTACT and meshtying
  namespace CONTACT
  {

    /// Type of contact or meshtying application
    /// (this enum represents the input file parameter APPLICATION)
    enum ApplicationType
    {
      app_none,             ///< no meshtying or contact
      app_mortarcontact,    ///< mortar contact
      app_mortarmeshtying,  ///< mortar meshtying
      app_beamcontact       ///< beam contact
    };

    /// Type of contact friction law
    /// (this enum represents the input file parameter FRICTION)
    enum FrictionType
    {
      friction_none,    ///< no friction
      friction_stick,   ///< perfect stick
      friction_tresca,  ///< Tresca friction law
      friction_coulomb  ///< Coulomb friction law
    };

    /// Type of employed solving strategy
    /// (this enum represents the input file parameter STRATEGY)
    enum SolvingStrategy
    {
      solution_lagmult,  ///< method of lagrange multipliers
      solution_penalty,  ///< penalty approach
      solution_auglag    ///< augmented lagrange approach
    };
    
    /// Type of linear system setup and solution
    /// (this enum represents the input file parameter SYSTEM)
    enum SystemType
    {
      system_condensed,  ///< condensed system
      system_spcoupled,  ///< saddle point system, coupled solution
      system_spsimpler   ///< saddle point system, SIMPLE(R) preconditioner
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

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
#endif  // #ifdef CCADISCRET
