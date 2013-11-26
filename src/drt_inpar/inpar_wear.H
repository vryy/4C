/*----------------------------------------------------------------------*/
/*!
\file inpar_wear.H

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

Questions? Contact Prof. Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_WEAR_H
#define INPAR_WEAR_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  /// Input parameters for structural wear
  namespace CONTACT
  {
    /// Type of contact wear law
    /// (this enum represents the input file parameter WEAR)
    enum WearLaw
    {
      wear_none,     ///< no wear
      wear_archard   ///< Archard wear law
    };

    /// Definition of contact wear surface
    /// (this enum represents the input file parameter BOTH_SIDED_WEAR)
    enum WearSide
    {
      wear_slave,     ///< wear on slave side
      wear_both_map  ///< map slave-side wear to master side
    };

    /// Definition of contact wear algorithm
    /// (this enum represents the input file parameter WEARTYPE)
    enum WearType
    {
      wear_impl,
      wear_expl,
      wear_discr
    };

    /// Definition of wear shape functions
    /// (this enum represents the input file parameter WEAR_SHAPEFCN)
    enum WearShape
    {
      wear_shape_dual,
      wear_shape_standard
    };

    /// Definition of wear-ALE coupling algorithm
    /// (this enum represents the input file parameter WEAR_COUPALGO)
    enum WearCoupAlgo
    {
      wear_stagg,
      wear_iterstagg,
      wear_monolithic
    };
  }  // namespace CONTACT
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_WEAR_H
