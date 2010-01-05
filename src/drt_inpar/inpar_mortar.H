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

  /// Input parameters for mortar coupling
  namespace MORTAR
  {
    
    /// Type of employed set of Lagrange multiplier shape functions
    /// (this enum represents the input file parameter SHAPEFCN)
    enum ShapeFcn
    {
      shape_standard,  ///< standard shape functions
      shape_dual  ///< dual shape functions
    };
    
    /// Type of mortar coupling search algorithm
    /// (this enum represents the input file parameter SEARCH_ALGORITHM)
    enum SearchAlgorithm
    {
      search_bfnode,  ///< brute force node-based
      search_bfele,  ///< brute force element-based
      search_binarytree  ///< binary tree element based
    };

  }  // namespace MORTAR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MORTAR_H
#endif  // #ifdef CCADISCRET
