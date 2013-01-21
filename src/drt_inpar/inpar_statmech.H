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
#ifndef INPAR_STATMECH_H
#define INPAR_STATMECH_H

#include "inpar_parameterlist_utils.H"

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
      statout_orientationcorrelation,  ///< measuring correlation of orientation between different elements
      statout_anisotropic,  ///< measuring anisotropic diffusion of filament
      statout_viscoelasticity, ///< measuring viscoelastic properties of filament networks
      statout_gmsh, ///< writing special output in order to generate a Gmsh movie
      statout_structanaly, ///< writing special output in order to determine network structure
      statout_octree, ///< writing output for octree visualization and checks
      statout_loom, ///< writing output for loom setup: a) nearest neighbor b) attractive forces c) coverage
      statout_loomelnrg ///< writing output for filament elastic energy
    };

    /// (this enum represents the input file parameter STAT_OUTPUT)
    enum FrictionModel
    {
      frictionmodel_none,                  ///< no friction model
      frictionmodel_isotropiclumped,       ///< isotropic, lumped damping matrix
      frictionmodel_isotropicconsistent,   ///< isotropic, consistent damping matrix
      frictionmodel_anisotropicconsistent  ///< anisotropic, consistent damping matrix
    };

    /// (this enum represents the input file parameter DBCTYPE)
    enum DBCType
    {
      dbctype_none,        ///< no Dirichlet Boundary conditions
      dbctype_std,          ///< standard Dirichlet Boundary conditions
      dbctype_shearfixed,  ///< special DBC application in case of rheological measurements with fixed set of nodes
      dbctype_sheartrans,  ///< special DBC application in case of rheological measurements with transient set of nodes
      dbctype_pinnodes,  ///< pin down individual nodes when a linker is attached
    };
  }  // namespace STATMECH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATMECH_H
