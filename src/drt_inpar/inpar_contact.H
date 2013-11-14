/*----------------------------------------------------------------------*/
/*!
\file inpar_contact.H

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
#ifndef INPAR_CONTACT_H
#define INPAR_CONTACT_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for structural CONTACT and meshtying
  namespace CONTACT
  {

    /// Type of contact or meshtying application
    /// (this enum represents the input file parameter APPLICATION)
    enum ApplicationType
    {
      app_none,             ///< no meshtying or contact
      app_mortarcontact,    ///< mortar contact
      app_mortarmeshtying,  ///< mortar meshtying
      app_beamcontact       ///< beam contact
    };

    /// Type of contact friction law
    /// (this enum represents the input file parameter FRICTION)
    enum FrictionType
    {
      friction_none,    ///< no friction
      friction_stick,   ///< perfect stick
      friction_tresca,  ///< Tresca friction law
      friction_coulomb  ///< Coulomb friction law
    };

    /// Type of employed solving strategy
    /// (this enum represents the input file parameter STRATEGY)
    enum SolvingStrategy
    {
      solution_lagmult,  ///< method of lagrange multipliers
      solution_penalty,  ///< penalty approach
      solution_auglag    ///< augmented lagrange approach
    };
    
    /// Type of linear system setup and solution
    /// (this enum represents the input file parameter SYSTEM)
    enum SystemType
    {
      system_condensed,  ///< condensed system
      system_spcoupled,  ///< saddle point system, coupled solution
      system_spsimpler   ///< saddle point system, SIMPLE(R) preconditioner
    };
    
    /// Type of energy and momentum output
    /// (this enum represents the input file parameter EMOUTPUT)
    enum EmOutputType
    {
      output_none,    ///< no output
      output_screen,  ///< print to screen
      output_file,    ///< print to file
      output_both     ///< print to screen and file
    };

    /// Type of analytical solution for error norm computation
    /// (this enum represents the input file parameter ERROR_NORMS)
    /// (more details on analytical solutions, see drt_contact/contact_analytical.cpp)
    enum ErrorNorms
    {
      errornorms_none,       ///< no error norm computation
      errornorms_zero,       ///< error norms for zero analytical solution
      errornorms_bending,    ///< error norms for beam bending problem
      errornorms_sphere,     ///< error norms for pressurized sphere problem
      errornorms_thicksphere ///< error norms for thick pressurized sphere problem
    };

    /// Application of a smoothed tangent field
    /// (this enum represents the input file parameter BEAMS_SMOOTHING)
    enum Smoothing
    {
      bsm_none,         ///< no smoothing
      bsm_smoothed,     ///< smoothing for all pairs
      bsm_partially,    ///< smoothing only for pairs which are not found with "no smoothing"
      bsm_cpp           ///< smoothing only for the closest point projection; element evaluation without smoothing"
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

    // problemtype to avoid use of drt_globalproblem.H
    enum Problemtype
    {
      structure,
      tsi,
      structalewear,
      other
    };

  }  // namespace CONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CONTACT_H
