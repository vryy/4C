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
#ifdef CCADISCRET
#ifndef INPAR_MORTAR_H
#define INPAR_MORTAR_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for MORTAR coupling
  namespace MORTAR
  {
    
    /// Type of employed set of Lagrange multiplier shape functions
    /// (this enum represents the input file parameter SHAPEFCN)
    enum ShapeFcn
    {
      shape_undefined, ///< undefined
      shape_standard,  ///< standard shape functions
      shape_dual       ///< dual shape functions
    };
    
    /// Type of Lagrange multiplier interpolation for 3D quadratic case
    /// (this enum represents the input file parameter LAGMULT_QUAD3D)
    enum LagMultQuad3D
    {
      lagmult_undefined,   ///< undefined (default for cases other than 3D quadratic)
      lagmult_quad_quad,   ///< quadratic interpolation, quadratic testing
      lagmult_quad_pwlin,  ///< quadratic interpolation, piecewise linear testing
      lagmult_quad_lin,    ///< quadratic interpolation, linear testing
      lagmult_pwlin_pwlin, ///< piecewise linear interpolation, piecewise linear testing
      lagmult_lin_lin      ///< linear interpolation, linear testing
    };
        
    /// Type of mortar coupling search algorithm
    /// (this enum represents the input file parameter SEARCH_ALGORITHM)
    enum SearchAlgorithm
    {
      search_bfnode,      ///< brute force node-based
      search_bfele,       ///< brute force element-based
      search_binarytree   ///< binary tree element based
    };

  }  // namespace MORTAR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MORTAR_H
#endif  // #ifdef CCADISCRET
