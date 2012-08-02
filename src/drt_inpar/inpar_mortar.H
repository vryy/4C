/*----------------------------------------------------------------------*/
/*!
\file inpar_mortar.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

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
    /// (this enum represents the input file parameter SHAPEFCN)
    enum ShapeFcn
    {
      shape_undefined,      ///< undefined
      shape_standard,       ///< standard shape functions
      shape_dual,           ///< dual shape functions
      shape_petrovgalerkin  ///< Petrov-Galerkin approach
    };
    
    /// Type of Lagrange multiplier interpolation for quadratic FE case
    /// (this enum represents the input file parameter LAGMULT_QUAD)
    enum LagMultQuad
    {
      lagmult_undefined,   ///< undefined (default)
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
      search_bfele,       ///< brute force element-based
      search_binarytree   ///< binary tree element based
    };

    /// Type of parallel redistribution algorithm
    /// (this enum represents the input file parameter PARALLEL_REDIST)
    enum ParRedist
    {
      parredist_none,      ///< no redistribution
      parredist_static,    ///< static redistribution (at t=0 and after restart)
      parredist_dynamic    ///< dynamic redistribution
    };

    /// Type of analytical solution for error norm computation
    /// (this enum represents the input file parameter ERROR_NORMS)
    /// (more details on analytical solutions, see drt_mortar/mortar_analytical.cpp)
    enum ErrorNorms
    {
      errornorms_none,       ///< no error norm computation
      errornorms_zero,       ///< error norms for zero analytical solution
      errornorms_bending,    ///< error norms for beam bending problem
      errornorms_sphere,     ///< error norms for pressurized sphere problem
      errornorms_thicksphere ///< error norms for thick pressurized sphere problem
    };
  }  // namespace MORTAR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MORTAR_H
