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
      contact_none,  ///< ???
      contact_normal,  ///< ???
      contact_frictional,  ///< ???
      contact_meshtying  ///< ???
    };

    /// Type of contact friction law
    /// (this enum represents the input file parameter FRICTION)
    enum ContactFrictionType
    {
      friction_none,  ///< ???
      friction_stick,  ///< ???
      friction_tresca,  ///< ???
      friction_coulomb  ///< ???
    };

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
#endif  // #ifdef CCADISCRET
