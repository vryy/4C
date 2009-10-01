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
      statout_endtoendlog,  ///< measuring end-to-end length of a single filament at time steps with exponentially increasing distance
      statout_endtoendconst,  ///< measuring end-to-end length of a single filament at time steps with constant distance
      statout_anisotropic,  ///< measuring anisotropic diffusion of filament
      statout_viscoelasticity, ///< measuring viscoelastic properties of filament networks
      statout_gmsh ///< writing special output in order to generate a Gmsh movie
    };
    
    /// (this enum represents the input file parameter STAT_OUTPUT)
    enum FrictionModel
    {
      frictionmodel_none,                  ///< no friction model
      frictionmodel_isotropiclumped,       ///< isotropic, lumped damping matrix
      frictionmodel_isotropicconsistent,   ///< isotropic, consistent damping matrix
      frictionmodel_anisotropicconsistent  ///< anisotropic, consistent damping matrix
    };

  }  // namespace STATMECH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATMECH_H
#endif  // #ifdef CCADISCRET
