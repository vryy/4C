/*-----------------------------------------------------------------------------------------------*/
/*!
\file beaminteraction_geometry_utils.H

\brief utility functions for geometric problems associated with beam-to-? interactions

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMINTERACTION_GEO_UTILS_H
#define BEAMINTERACTION_GEO_UTILS_H

#include "../linalg/linalg_fixedsizematrix.H"

#include "../drt_lib/drt_element.H"


namespace BEAMINTERACTION
{
namespace GEO
{
  // point-to-curve projection: solve minimal distance problem
  // convergence criteria for local Newton's method
  const unsigned int POINT_TO_CURVE_PROJECTION_MAX_NUM_ITER = 50;
  const double POINT_TO_CURVE_PROJECTION_TOLERANCE_RESIDUUM = 1.0e-10;
  const double POINT_TO_CURVE_PROJECTION_TOLERANCE_INCREMENT = 1.0e-10;
  // threshold values for sanity checks
  const double POINT_TO_CURVE_PROJECTION_IDENTICAL_POINTS_TOLERANCE = 1.0e-12;
  const double POINT_TO_CURVE_PROJECTION_NONUNIQUE_MINIMAL_DISTANCE_TOLERANCE = 1.0e-12;

  /** \brief solves minimal distance problem to find the closest point on a 3D spatial curve
   *         (i.e. its curve parameter value) relative to a given point
   *         a.k.a 'unilateral' closest-point projection
   *
   *  \author grill, meier
   *  \date 10/17, 01/14 */
  template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
  bool PointToCurveProjection(
      LINALG::TMatrix<T,3,1> const& r_slave,
      T& xi_master,
      double const& xi_master_initial_guess,
      const LINALG::TMatrix<T,3*numnodes*numnodalvalues,1>& master_centerline_dof_values,
      const DRT::Element::DiscretizationType& master_distype,
      double master_ele_ref_length);

  /** \brief evaluates residual of orthogonality condition for so-called unilateral closest-point
   *         projection, i.e. a point-to-curve projection
   *
   *  \author grill, meier
   *  \date 10/17, 10/14 */
  template<typename T>
  void EvaluatePointToCurveOrthogonalityCondition(
      T& f,
      const LINALG::TMatrix<T,3,1>& delta_r,
      const double norm_delta_r,
      const LINALG::TMatrix<T,3,1>& r_xi_master);

  /** \brief evaluates Jacobian of orthogonality condition for so-called unilateral closest-point
   *         projection, i.e. a point-to-curve projection
   *
   *  \author grill, meier
   *  \date 10/17, 10/14 */
  template<typename T>
  bool EvaluateLinearizationPointToCurveOrthogonalityCondition(
      T& df,
      const LINALG::TMatrix<T,3,1>& delta_r,
      const double norm_delta_r,
      const LINALG::TMatrix<T,3,1>& r_xi_master,
      const LINALG::TMatrix<T,3,1>& r_xixi_master);

  /** \brief compute linearization of parameter coordinate on master if determined by a
   *         point-to-curve projection
   *
   *  \author grill, meier
   *  \date 10/17, 10/14 */
  template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
  void CalcLinearizationPointToCurveProjectionParameterCoordMaster(
      LINALG::TMatrix<T,1,3*numnodes*numnodalvalues>& lin_xi_master_slaveDofs,
      LINALG::TMatrix<T,1,3*numnodes*numnodalvalues>& lin_xi_master_masterDofs,
      const LINALG::TMatrix<T,3,1>& delta_r,
      const LINALG::TMatrix<T,3,1>& r_xi_master,
      const LINALG::TMatrix<T,3,1>& r_xixi_master,
      const LINALG::TMatrix<double,3,3*numnodes*numnodalvalues>& N_slave,
      const LINALG::TMatrix<T,3,3*numnodes*numnodalvalues>& N_master,
      const LINALG::TMatrix<T,3,3*numnodes*numnodalvalues>& N_xi_master);

  /** \brief point-to-curve projection:
   *         partial derivatives of the parameter coordinate on master xi_master with respect to
   *         centerline position of slave point, master point and centerline tangent of master
   *
   *  \author grill
   *  \date 10/17 */
  template<typename T>
  void CalcPointToCurveProjectionParameterCoordMasterPartialDerivs(
      LINALG::TMatrix<T,1,3>& xi_master_partial_r_slave,
      LINALG::TMatrix<T,1,3>& xi_master_partial_r_master,
      LINALG::TMatrix<T,1,3>& xi_master_partial_r_xi_master,
      const LINALG::TMatrix<T,3,1>& delta_r,
      const LINALG::TMatrix<T,3,1>& r_xi_master,
      const LINALG::TMatrix<T,3,1>& r_xixi_master);

  /** \brief point-to-curve projection:
   *         partial derivative of the orthogonality condition with respect to parameter coordinate
   *         on master xi_master
   *
   *  \author grill
   *  \date 10/17 */
  template<typename T>
  void CalcPTCProjectionOrthogonalityConditionPartialDerivParameterCoordMaster(
      T& orthogon_condition_partial_xi_master,
      const LINALG::TMatrix<T,3,1>& delta_r,
      const LINALG::TMatrix<T,3,1>& r_xi_master,
      const LINALG::TMatrix<T,3,1>& r_xixi_master);

  /** \brief point-to-curve projection:
   *         partial derivative of the orthogonality condition with respect to centerline position
   *         on slave
   *
   *  \author grill
   *  \date 10/17 */
  template<typename T>
  void CalcPTCProjectionOrthogonalityConditionPartialDerivClPosSlave(
      LINALG::TMatrix<T,1,3>& orthogon_condition_partial_r_slave,
      const LINALG::TMatrix<T,3,1>& r_xi_master);

  /** \brief point-to-curve projection:
   *         partial derivative of the orthogonality condition with respect to centerline position
   *         on master
   *
   *  \author grill
   *  \date 10/17 */
  template<typename T>
  void CalcPTCProjectionOrthogonalityConditionPartialDerivClPosMaster(
      LINALG::TMatrix<T,1,3>& orthogon_condition_partial_r_master,
      const LINALG::TMatrix<T,3,1>& r_xi_master);

  /** \brief point-to-curve projection:
   *         partial derivative of the orthogonality condition with respect to centerline tangent
   *         on master
   *
   *  \author grill
   *  \date 10/17 */
  template<typename T>
  void CalcPTCProjectionOrthogonalityConditionPartialDerivClTangentMaster(
      LINALG::TMatrix<T,1,3>& orthogon_condition_partial_r_xi_master,
      const LINALG::TMatrix<T,3,1>& delta_r);

  /** \brief calculate angle enclosed by two vectors a and b
   *
   *  \author grill, meier
   *  \date 10/17, 10/14 */
  template <typename T>
  void CalcEnclosedAngle(T& angle, T& cosine_angle, const LINALG::TMatrix<T,3,1>& a,
      const LINALG::TMatrix<T,3,1>& b);

}
}

#endif
