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
      thermalbath_none,  ///< no theraml bath around structure (ordinary simulation)
      thermalbath_uniform,  ///< uniform thermal bath around structures
      thermalbath_shearflow  ///< structures immersed into shear flow
    };

    /// (this enum represents the input file parameter STAT_OUTPUT)
    enum StatOutput
    {
      statout_none,  ///< no special output
      statout_endtoendlog,  ///< measuring end-to-end length of a single filament
      statout_endtoendergodicity,  ///< measuring ergodicity of single filament simulations
      statout_viscoelasticity, ///< measuring viscoelastic properties of filament networks
      statout_gmsh///< writing special output in order to generate a Gmsh movie
    };

  }  // namespace STATMECH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATMECH_H
#endif  // #ifdef CCADISCRET
