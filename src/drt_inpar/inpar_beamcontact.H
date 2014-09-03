/*----------------------------------------------------------------------*/
/*!
\file inpar_beamcontact.H

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
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_BEAMCONTACT_H
#define INPAR_BEAMCONTACT_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for beam contact
  namespace BEAMCONTACT
  {
    /// Application of a smoothed tangent field
    /// (this enum represents the input file parameter BEAMS_SMOOTHING)
    enum Smoothing
    {
      bsm_none,         ///< no smoothing
      bsm_cpp           ///< smoothing only for the closest point projection; element evaluation without smoothing"
    };

    enum Damping
    {
      bd_no,           ///< no damping force
      bd_yes           ///< application of a contact damping force
    };

    /// Application of a smoothed tangent field
    /// (this enum represents the input file parameter BEAMS_SMOOTHING)
    enum PenaltyLaw
    {
      pl_lp,          ///< linear penalty law
      pl_qp,          ///< quadratic penalty law
      pl_lnqp,        ///< linear penalty law with quadratic regularization for negative gaps
      pl_lpqp,        ///< linear penalty law with quadratic regularization for positive gaps
      pl_lpcp,        ///< linear penalty law with cubic regularization for positive gaps
      pl_lpdqp,       ///< linear penalty law with double quadratic regularization for positive gaps
      pl_lpep         ///< linear penalty law with exponential regularization for positive gaps
    };

    /// Beam Contact Octree  and Bounding Box Type
    /// (this enum represents the input file parameter BEAMS_OCTREEBBOX)
    enum OctreeType
    {
      boct_none,       ///< no boundig box -> no octree
      boct_aabb,       ///< axis aligned bounding boxes
      boct_cobb,       ///< cylindrical oriented bounding boxes
      boct_spbb        ///< spherical bounding boxes
    };

  }  // namespace BEAMCONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_BEAMCONTACT_H
