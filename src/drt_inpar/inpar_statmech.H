/*-----------------------------------------------------------*/
/*!
\file inpar_statmech.H

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger

\level 2

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_STATMECH_H
#define INPAR_STATMECH_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace STATMECH
  {

    /// (this enum represents the input file parameter SPECIALOUTPUT)
    enum StatOutput
    {
      statout_none,  ///< no special output
      statout_endtoendlog,  ///< measuring end-to-end length of a single filament at time steps with exponentially increasing distance
      statout_endtoendconst,  ///< measuring end-to-end length of a single filament at time steps with constant distance
      statout_orientationcorrelation,  ///< measuring correlation of orientation between different elements
      statout_anisotropic,  ///< measuring anisotropic diffusion of filament
      statout_viscoelasticity, ///< measuring viscoelastic properties of filament networks
      statout_networkcreep, ///< measuring displacements upon force application
      statout_networkrelax, ///< measuring stress relaxation of a network
      statout_networkdispfield, ///< output network nodal displacement field
      statout_structanaly, ///< writing special output in order to determine network structure
      statout_octree, ///< writing output for octree visualization and checks
      statout_loom, ///< writing output for loom setup: a) nearest neighbor b) attractive forces c) coverage
      statout_loomelnrg, ///< writing output for filament elastic energy
      statout_motassay, ///< writing output for motility assay (active linker setup)
      statout_linkerlength, //< writing output of spring linkers
      statout_deltatheta, // < writing output of change in inclusive angle of filaments when connected to linkers
      // Membrane cases
      statout_vesceqshapes
    };

    /// (this enum represents the input file parameter DBCTYPE)
    enum DBCType
    {
      dbctype_none,        ///< no Dirichlet Boundary conditions
      dbctype_std,          ///< standard Dirichlet Boundary conditions
      dbctype_shearfixed,  ///< special DBC application in case of rheological measurements with fixed set of nodes
      dbctype_shearfixeddel,  ///< special DBC application in case of rheological measurements with fixed set of nodes, delete cut elements
      dbctype_sheartrans,  ///< special DBC application in case of rheological measurements with transient set of nodes
      dbctype_pinnodes,    ///< pin down individual nodes when a linker is attached
      dbctype_affineshear,  ///< apply affine shear displacement to all nodes
      dbctype_affinesheardel,  ///< apply affine shear displacement to all nodes, delete cut elements
      dbctype_movablesupport1d,  ///< constrain all displacement except for one direction
    };

    /// (this enum represents the input file parameter NBCTYPE)
    enum NBCType
    {
      nbctype_std,             ///< standard Neumann Boundary conditions defined by input file
      nbctype_constcreep,      ///< constant Neumann Boundary condition
      nbctype_randompointforce ///< application of time-dependent Neumann term at a random node
    };

    /// (this enum represents the input file parameter BINDINGSITESEARCH)
    enum BSSearchType
    {
      bsstype_volpart, ///< volume partitioning
      bsstype_binning, ///< binning
      bsstype_octree   ///< octree
    };

    /// (this enum represents the input file parameter NETWORKTYPE
    enum NetworkType
    {
      networktype_std, /// < standard actin network
      networktype_casimir, /// < casimir force type setup
      networktype_motassay /// < motility assay network setup
    };

    /// (this enum represents the input file parameter LINKERMODEL
    enum LinkerModel
    {
      linkermodel_none, /// < no linkers at all
      linkermodel_std, /// < standard linkers with chemical association and dissociation rate
      linkermodel_stdintpol, /// <standard linker with interpolated binding position on filament
      linkermodel_bellseq, ///< like standard linker but a force dependent dissociation rate
      linkermodel_bellseqintpol, //< linker with force dependent dissociation rate and interpolated binding spot positions
      linkermodel_active, ///< like standard linker but with a force dependent dissociation rate and active contraction and extension
      linkermodel_activeintpol, ///< active linker with interpolated binding spot positions
      linkermodel_myosinthick ///< active linker myosin mini filament
    };

    /// (this enum represents the input file parameter FILAMENTMODEL
    enum FilamentModel
    {
      filamentmodel_std, /// < standard filament type
      filamentmodel_helical /// < filament with helical orientation of binding spots
    };

    /// (this enum represents the input file parameter FILAMENTMODEL
    enum AreaPenaltyType
    {
      areapenalty_local,   // Local area penalty
      areapenalty_global,  // Global area penalty
      areapenalty_none     // No area penalty
    };

    /// set the statmech parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace STATMECH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATMECH_H
