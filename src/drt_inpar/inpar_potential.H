/*----------------------------------------------------------------------*/
/*!
\file inpar_potential.H

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_POTENTIAL_H
#define INPAR_POTENTIAL_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for structural contact
  namespace POTENTIAL
  {
    /// Type of effects which should be considered
    /// (this enum represents the input file parameter POTENTIAL)
    enum PotentialType
    {
      potential_surface,          ///< only surface effects are considered
      potential_volume,           ///< only volume effects are considered
      potential_surfacevolume,    ///< surface and volume effects are considered
      potential_surface_fsi,      ///< only surface and fluid effects are considered
      potential_volume_fsi,       ///< only volume and fluid effects are considered
      potential_surfacevolume_fsi ///< surface, volume and fluid effects are considered
    };

    /// Type of approximation
    /// (this enum represents the input file parameter POTENTIAL)
    enum ApproximationType
    {
      approximation_none,       ///< no approximation
      approximation_surface,    ///< surface approximation
      approximation_point       ///< point approximation
    };
  }  // namespace POTENTIAL

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_POTENTIAL_H
#endif  // #ifdef CCADISCRET

