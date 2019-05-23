/*----------------------------------------------------------------------*/
/*!

\brief input parameter definitions for beam contact

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_BEAMCONTACT_H
#define INPAR_BEAMCONTACT_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace BEAMCONTACT
  {
    /// Type of employed solving strategy
    /// (this enum represents the input file parameter BEAMS_STRATEGY)
    enum Strategy
    {
      bstr_none,     ///< no beam contact
      bstr_penalty,  ///< penalty method
      bstr_gmshonly  ///< misuse of beam contact module for GMSH output
    };

    /// Type of employed solving strategy
    /// (this enum represents the input file parameter BEAMS_STRATEGY)
    enum Modelevaluator
    {
      bstr_old,      ///<  use old beaminteraction model
      bstr_standard  ///<  use new beamcontact submodel
    };

    /// Application of a smoothed tangent field
    /// (this enum represents the input file parameter BEAMS_SMOOTHING)
    enum Smoothing
    {
      bsm_none,  ///< no smoothing
      bsm_cpp    ///< smoothing only for the closest point projection; element evaluation without
                 ///< smoothing"
    };

    /// Application of a contact damping force
    /// (this enum represents the input file parameter BEAMS_DAMPING)
    enum Damping
    {
      bd_no,  ///< no damping force
      bd_yes  ///< application of a contact damping force
    };

    /// Application of a smoothed tangent field
    /// (this enum represents the input file parameter BEAMS_SMOOTHING)
    enum PenaltyLaw
    {
      pl_lp,     ///< linear penalty law
      pl_qp,     ///< quadratic penalty law
      pl_lnqp,   ///< linear penalty law with quadratic regularization for negative gaps
      pl_lpqp,   ///< linear penalty law with quadratic regularization for positive gaps
      pl_lpcp,   ///< linear penalty law with cubic regularization for positive gaps
      pl_lpdqp,  ///< linear penalty law with double quadratic regularization for positive gaps
      pl_lpep    ///< linear penalty law with exponential regularization for positive gaps
    };

    /// Beam Contact Octree and Bounding Box Type
    /// (this enum represents the input file parameter BEAMS_OCTREEBBOX)
    enum OctreeType
    {
      boct_none,  ///< no boundig box -> no octree
      boct_aabb,  ///< axis aligned bounding boxes
      boct_cobb,  ///< cylindrical oriented bounding boxes
      boct_spbb   ///< spherical bounding boxes
    };

    /// data format for written numeric data via vtk
    enum OutputDataFormat
    {
      binary,
      ascii,
      vague
    };

    /// set the beam contact parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace BEAMCONTACT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_BEAMCONTACT_H
