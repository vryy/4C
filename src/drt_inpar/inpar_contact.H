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

  /// Input parameters for structural contact
  namespace CONTACT
  {

    /// Type of structural contact problem
    /// (this enum represents the input file parameter CONTACT)
    enum ContactType
    {
      contact_none,  ///< no contact
      contact_normal,  ///< frictionless contact
      contact_frictional,  ///< frictional contact
      contact_meshtying,  ///< tied contact (mesh tying)
      contact_beams  ///< beam contact
    };

    /// Type of contact friction law
    /// (this enum represents the input file parameter FRICTION)
    enum ContactFrictionType
    {
      friction_none,  ///< no friction
      friction_stick,  ///< perfect stick
      friction_tresca,  ///< Tresca friction law
      friction_coulomb  ///< Coulomb friction law
    };

    /// Type of employed solving strategy
    /// (this enum represents the input file parameter STRATEGY)
    enum SolvingStrategy
    {
      solution_lagmult,  ///< method of lagrange multipliers
      solution_penalty,  ///< penalty approach
      solution_auglag  ///< augmented lagrange approach
    };

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
#endif  // #ifdef CCADISCRET
