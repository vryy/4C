/*----------------------------------------------------------------------*/
/*!
\file inpar_mortar.H

\brief input parameters for mortar methods

\level 2

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_MORTAR_H
#define INPAR_MORTAR_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }  // namespace INPUT
}  // namespace DRT

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace MORTAR
  {
    /// Type of employed set of Lagrange multiplier shape functions
    /// (this enum represents the input file parameter LM_SHAPEFCN)
    enum ShapeFcn
    {
      shape_undefined,       ///< undefined
      shape_standard,        ///< standard shape functions
      shape_dual,            ///< dual shape functions
      shape_petrovgalerkin,  ///< Petrov-Galerkin approach
      shape_none             ///< for all methods w/o Lagrange multiplier interpolation
    };

    /// Type of Lagrange multiplier interpolation for quadratic FE case
    /// (this enum represents the input file parameter LM_QUADRATIC)
    enum LagMultQuad
    {
      lagmult_undefined,  ///< undefined
      lagmult_quad,       ///< quadratic interpolation
      lagmult_pwlin,      ///< piecewise linear interpolation
      lagmult_lin,        ///< linear interpolation
      lagmult_const       ///< element-wise constant interpolation (only for quadratic FE)
    };

    /// Type of mortar coupling search algorithm
    /// (this enum represents the input file parameter SEARCH_ALGORITHM)
    enum SearchAlgorithm
    {
      search_bfele,      ///< brute force element-based
      search_binarytree  ///< binary tree element based
    };

    /// Type of binary tree update
    /// (this enum represents the input file parameter BINARYTREE_UPDATETYPE)
    enum BinaryTreeUpdateType
    {
      BOTTOM_UP,  ///< indicates a bottom up update of binary tree
      TOP_DOWN    ///< indicates a top down update of binary tree
    };

    /// Type of mesh relocation
    /// (this enum represents the input file parameter MESH_RELOCATION)
    enum MeshRelocation
    {
      relocation_initial,   ///< only initial mesh relocation
      relocation_timestep,  ///< mesh relocation in every time step, but no initial mesh relocation
      relocation_none       ///< no mesh relocation
    };

    /// Type of redundancy in interface storage
    /// (this enum represents the input file parameter REDUNDANT_STORAGE)
    enum RedundantStorage
    {
      redundant_all,     ///< redundant storage of master AND slave surfaces
      redundant_master,  ///< redundant storage of master surface
      redundant_none     ///< no redundant storage
    };

    /// Type of parallel interface evaluation
    /// (this enum represents the input file parameter PARALLEL_STRATEGY)
    enum ParallelStrategy
    {
      ghosting_redundant,  ///< redundant storage of master or master&slave surfaces
      roundrobinevaluate,  ///< round robin loop for entire interface evaluation
      roundrobinghost,     ///< round robin loop for setup of ghosting only
      binningstrategy      ///< binning strategy
    };

    /// Type of meshtying/contact algorithm
    /// (this enum represents the input file parameter ALGORITHM)
    enum AlgorithmType
    {
      algorithm_mortar,  ///< mortar algorithm (segment-to-segment)
      algorithm_nts,     ///< node-to-segment algorithm
      algorithm_gpts,    ///< gp-to-segment algorithm
      algorithm_lts,     ///< line-to-segment algorithm
      algorithm_ltl,     ///< line-to-segment algorithm
      algorithm_ntl,     ///< node-to-line algorithm (coming soon...)
      algorithm_stl      ///< segment-to-line algorithm
    };

    /// Type of parallel redistribution algorithm
    /// (this enum represents the input file parameter PARALLEL_REDIST)
    enum ParRedist
    {
      parredist_none,    ///< no redistribution
      parredist_static,  ///< static redistribution (at t=0 and after restart)
      parredist_dynamic  ///< dynamic redistribution
    };

    /// Type of integration procedure
    /// (this enum represents the input file parameter INTTYPE)
    enum IntType
    {
      inttype_segments,    ///< segmentation of mortar interface
      inttype_elements,    ///< fast, elementwise integration
      inttype_elements_BS  ///< fast, elementwise intgration with boundary segmentation
    };

    /// Type of triangulation for segment-base d integration
    /// (this enum represents the input file parameter TRIANGULATION)
    enum Triangulation
    {
      triangulation_center,   ///< simpler center-based triangulation (see e.g. Popp et al. 2010)
      triangulation_delaunay  ///< delaunay triangulation
    };

    /// Determining, on which quadrature points biorthogonality is enforced
    enum ConsistentDualType
    {
      consistent_none,      ///< always use element GP (fastest option)
      consistent_boundary,  ///< use triangulation GPs only in partially integrated elements
      consistent_all,       ///< use triangulation GPs for all elements
    };

    enum class DBCHandling : int
    {
      do_nothing,
      remove_dbc_nodes_from_slave_side
    };

    /// set the mortar parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific mortar conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace MORTAR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MORTAR_H
