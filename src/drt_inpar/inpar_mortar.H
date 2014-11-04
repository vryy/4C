/*----------------------------------------------------------------------*/
/*!
\file inpar_mortar.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_MORTAR_H
#define INPAR_MORTAR_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for MORTAR coupling
  namespace MORTAR
  {

    /// Type of employed set of Lagrange multiplier shape functions
    /// (this enum represents the input file parameter LM_SHAPEFCN)
    enum ShapeFcn
    {
      shape_undefined,      ///< undefined
      shape_standard,       ///< standard shape functions
      shape_dual,           ///< dual shape functions
      shape_petrovgalerkin  ///< Petrov-Galerkin approach
    };

    /// Type of Lagrange multiplier interpolation for quadratic FE case
    /// (this enum represents the input file parameter LM_QUADRATIC)
    enum LagMultQuad
    {
      lagmult_undefined,   ///< undefined
      lagmult_quad,        ///< quadratic interpolation
      lagmult_pwlin,       ///< piecewise linear interpolation
      lagmult_lin          ///< linear interpolation
    };

    /// Type of mortar coupling search algorithm
    /// (this enum represents the input file parameter SEARCH_ALGORITHM)
    enum SearchAlgorithm
    {
      search_bfele,       ///< brute force element-based
      search_binarytree   ///< binary tree element based
    };

    /// Type of mesh relocation
    /// (this enum represents the input file parameter MESH_RELOCATION)
    enum MeshRelocation
    {
      relocation_initial,  ///< only initial mesh relocation
      relocation_timestep, ///< mesh relocation in every time step, but no initial mesh relocation
      relocation_none      ///< no mesh relocation
    };

    /// Type of redudancy in interface storage
    /// (this enum represents the input file parameter REDUNDANT_STORAGE)
    enum RedundantStorage
    {
      redundant_all,      ///< redundant storage of master AND slave surfaces
      redundant_master,   ///< redundant storage of master surface
      redundant_none      ///< no redundant storage
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
      algorithm_mortar,    ///< mortar algorithm
      algorithm_nts        ///< node-to-segment algorithm
    };

    /// Type of parallel redistribution algorithm
    /// (this enum represents the input file parameter PARALLEL_REDIST)
    enum ParRedist
    {
      parredist_none,      ///< no redistribution
      parredist_static,    ///< static redistribution (at t=0 and after restart)
      parredist_dynamic    ///< dynamic redistribution
    };

    /// Type of integration procedure
    /// (this enum represents the input file parameter INTTYPE)
    enum IntType
    {
      inttype_segments,    ///< segmentation of mortar interface
      inttype_elements,    ///< fast, elementwise integration
      inttype_elements_BS  ///< fast, elementwise intgration with boundary segmentation
    };
  }  // namespace MORTAR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MORTAR_H
