/*----------------------------------------------------------------------*/
/*!
\file inpar_statmech.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_STATMECH_H
#define INPAR_STATMECH_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace STATMECH
  {

    /// Type of thermal bath affecting elements
    /// (this enum represents the input file parameter THERMALBATH)
    enum ThermalBathType
    {
      thermalbath_none,  ///< ???
      thermalbath_uniform,  ///< ???
      thermalbath_shearflow  ///< ???
    };

    /// (this enum represents the input file parameter STAT_OUTPUT)
    enum StatOutput
    {
      statout_none,  ///< ???
      statout_endtoendlog,  ///< ???
      statout_endtoendergodicity,  ///< ???
      statout_viscoelasticity  ///< ???
    };

  }  // namespace STATMECH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATMECH_H
#endif  // #ifdef CCADISCRET
