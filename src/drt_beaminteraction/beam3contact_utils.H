/*----------------------------------------------------------------------------*/
/*!
\file beam3contact_utils.H

\brief A set of utility functions for beam contact

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/
#ifndef BEAM3CONTACT_UTILS_H
#define BEAM3CONTACT_UTILS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../linalg/linalg_fixedsizematrix.H"
#include "Sacado.hpp"

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/largerotations.H"

typedef Sacado::Fad::DFad<double> FAD;

// Todo change to BEAMINTERACTION namespace
namespace BEAMCONTACT
{
  /*!
  \brief Calculate angle encompassed by two lines: returns an angle \in [0;pi/2]
  */
  double CalcAngle(LINALG::TMatrix<double, 3, 1> a, LINALG::TMatrix<double, 3, 1> b);

  /*!
  \brief Get closest distance between the endpoints of two lines
  */
  template <typename type>
  type GetClosestEndpointDist(LINALG::TMatrix<type, 3, 1> r1_a, LINALG::TMatrix<type, 3, 1> r1_b,
      LINALG::TMatrix<type, 3, 1> r2_a, LINALG::TMatrix<type, 3, 1> r2_b);

  /*!
  \brief Set primary displacement DoFs for automatic differentiation with Sacado
  */
  template <int numnodes, int numnodalvalues>
  void SetFADDispDofs(LINALG::TMatrix<FAD, 3 * numnodes * numnodalvalues, 1>& ele1pos_,
      LINALG::TMatrix<FAD, 3 * numnodes * numnodalvalues, 1>& ele2pos_)
  {
    // The 2*3*numnodes*numnodalvalues primary DoFs are the components of the nodal positions /
    // tangents. The two (+2) additional degrees of freedom represent the dependency on the
    // parameter coordinates xi and eta, which is necessary in beam contact.
    for (int i = 0; i < 3 * numnodes * numnodalvalues; i++)
      ele1pos_(i).diff(i, 2 * 3 * numnodes * numnodalvalues + 2);

    for (int i = 0; i < 3 * numnodes * numnodalvalues; i++)
      ele2pos_(i).diff(3 * numnodes * numnodalvalues + i, 2 * 3 * numnodes * numnodalvalues + 2);

    return;
  }

  /*!
\brief BTS-Contact: Set primary displacement DoFs for automatic differentiation with Sacado
*/
  template <int numnodessol, int numnodes, int numnodalvalues>
  void SetFADDispDofs(LINALG::TMatrix<FAD, 3 * numnodes * numnodalvalues, 1>& ele1pos,
      LINALG::TMatrix<FAD, 3 * numnodessol, 1>& ele2pos, const int& numvar)
  {
    for (int i = 0; i < 3 * numnodes * numnodalvalues; i++)
      ele1pos(i).diff(i, 3 * numnodes * numnodalvalues + 3 * numnodessol + numvar);

    for (int i = 0; i < 3 * numnodessol; i++)
      ele2pos(i).diff(3 * numnodes * numnodalvalues + i,
          3 * numnodes * numnodalvalues + 3 * numnodessol + numvar);

    return;
  }

  /*!
  \brief Set primary parameter coordinate DoFs for automatic differentiation with Sacado
  */
  template <int numnodes, int numnodalvalues>
  void SetFADParCoordDofs(FAD& xi, FAD& eta)
  {
    // The 2*3*numnodes*numnodalvalues primary DoFs are the components of the nodal positions /
    // tangents. The two (+2) additional degrees of freedom represent the dependency on the
    // parameter coordinates xi and eta, which is necessary in beam contact.
    xi.diff((2 * 3 * numnodes * numnodalvalues + 1) - 1, 2 * 3 * numnodes * numnodalvalues + 2);
    eta.diff((2 * 3 * numnodes * numnodalvalues + 2) - 1, 2 * 3 * numnodes * numnodalvalues + 2);

    return;
  }

  /*!
  \brief BTS-Contact: Set primary parameter coordinate DoFs for automatic differentiation with
  Sacado
  */
  template <int numnodessol, int numnodes, int numnodalvalues>
  void SetFADParCoordDofs(FAD& xi1, FAD& xi2, FAD& eta)
  {
    xi1.diff(3 * numnodes * numnodalvalues + 3 * numnodessol,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 3);
    xi2.diff(3 * numnodes * numnodalvalues + 3 * numnodessol + 1,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 3);
    eta.diff(3 * numnodes * numnodalvalues + 3 * numnodessol + 2,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 3);

    return;
  }

  /*!
  \brief BTS-Contact: Set primary parameter coordinate DoFs for automatic differentiation with
  Sacado
  */
  template <int numnodessol, int numnodes, int numnodalvalues>
  void SetFADParCoordDofs(FAD& xi1, FAD& xi2, FAD& eta_a, FAD& eta_b)
  {
    xi1.diff(3 * numnodes * numnodalvalues + 3 * numnodessol,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 4);
    xi2.diff(3 * numnodes * numnodalvalues + 3 * numnodessol + 1,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 4);
    eta_a.diff(3 * numnodes * numnodalvalues + 3 * numnodessol + 2,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 4);
    eta_b.diff(3 * numnodes * numnodalvalues + 3 * numnodessol + 3,
        3 * numnodes * numnodalvalues + 3 * numnodessol + 4);

    return;
  }

  /*!
  \brief Check, if current node belongs to a beam element
  */
  bool BeamNode(const DRT::Node& node);

  /*!
  \brief Check, if current node is used for centerline interpolation of a beam element
  */
  bool BeamCenterlineNode(const DRT::Node& node);

  /*!
  \brief Check, if current node belongs to a rigid sphere element
  */
  bool RigidsphereNode(const DRT::Node& node);

  /*!
  \brief Check, if current node is a beam element
  */
  bool BeamElement(const DRT::Element& element);

  /*!
  \brief Check, if current node is a rigid sphere element
  */
  bool RigidsphereElement(const DRT::Element& element);

  /*!
  \brief Check, if current node is a solid contact element
  */
  bool SolidContactElement(const DRT::Element& element);

  /*!
  \brief Check, if current node is a solid meshtying element
  */
  bool SolidMeshtyingElement(const DRT::Element& element);

  /*
  \brief Check, if two elements share a node -> neighbor elements
  */
  bool ElementsShareNode(const DRT::Element& element1, const DRT::Element& element2);

  /*
  \brief Calculate beam radius
  */
  double CalcEleRadius(const DRT::Element* ele);

  /*
  \brief Intersect two parallel cylinders
  */
  bool IntersectParallelCylinders(LINALG::TMatrix<double, 3, 1>& r1_a,
      LINALG::TMatrix<double, 3, 1>& r1_b, LINALG::TMatrix<double, 3, 1>& r2_a,
      LINALG::TMatrix<double, 3, 1>& r2_b, double& distancelimit);

  /*
  \brief Intersect two non-parallel, arbitrary oriented cylinders
  */
  bool IntersectArbitraryCylinders(LINALG::TMatrix<double, 3, 1>& r1_a,
      LINALG::TMatrix<double, 3, 1>& r1_b, LINALG::TMatrix<double, 3, 1>& r2_a,
      LINALG::TMatrix<double, 3, 1>& r2_b, double& distancelimit,
      std::pair<double, double>&
          closestpoints,  // The closest point are only set, if we have detected an intersection at
                          // a valid closest point with eta1_seg, eta2_seg \in [-1.0;1.0]
      bool etaset);       // bool to check, if the closest point coordinates have been set or not

  /*
  \brief Calculate closest distance of a point and a line
  */
  double CalcPointLineDist(LINALG::TMatrix<double, 3, 1>& rline_a,
      LINALG::TMatrix<double, 3, 1>& rline_b, LINALG::TMatrix<double, 3, 1>& rp, double& eta);

  /*
  \brief Determine inpute parameter representing the additive searchbox increment
  */
  double DetermineSearchboxInc(Teuchos::ParameterList& beamcontactparams);

  /*
  \brief Check if a given double lies within a prescribed interval (enlarged by the tolerance
  XIETATOL)
  */
  inline bool WithinInterval(double& testpoint, double& leftbound, double& rightbound)
  {
    // The tolerance XIETATOL makes the test more conservative, i.e. the testpoint
    // is assumed to be within the interval even if it actually is slightly outside
    if (testpoint > leftbound - XIETATOL and testpoint < rightbound + XIETATOL)
      return true;
    else
      return false;
  }

  /*
  \brief Get interval-id out of numberofintervals intervals, in which the given point lies
  */
  inline int GetIntervalId(double& point, int numberofintervals, bool leftbound)
  {
    int interval_id = 0;
    double unrounded_id = 0.0;

    // With the following formula we would get the exact interval Ids when inserting the coordinate
    // of the left bound of the interval. By inserting any double value lying within the interval,
    // we get as result a double value that is larger than the sought-after interval ID but smaller
    // than the next highter ID, i.e. we have to round down the solution.
    unrounded_id = (point + 1.0) / 2.0 * numberofintervals;
    interval_id = floor(unrounded_id);

    // Size Check: If the size of the interval becomes smaller than RELSEGMENTTOL times the interval
    // length we simply shift the point to the next higher/lower interval
    double segmenttol = RELSEGMENTTOL * numberofintervals / 2.0;

    // The necessary shifting procedure depends on the fact, whether we are searching for the left
    // bound
    if (leftbound)
    {
      if (fabs(1 + interval_id - unrounded_id) < segmenttol)
        dserror("Such small segmented integration intervals are not possible so far!");
    }
    // or for the right bound of the integration interval
    else
    {
      if (fabs(interval_id - unrounded_id) < segmenttol)
        dserror("Such small segmented integration intervals are not possible so far!");
    }

    if (interval_id < 0) dserror("Interval-ID can't be negative!");

    return interval_id;
  }

  /*
  \brief Get segment-id out of numberofsegments segmenets, in which the given point lies
  */
  inline int GetSegmentId(double& point, int numberofsegments)
  {
    int segment_id = 0;
    double unrounded_id = 0.0;

    // The method is similar to GetIntervalId() above. However, since it is purely used as output
    // quantity and it does does not influence algorithmic quantities such as the integration
    // interval length, we don't need a Size Check such as the method GetIntervalId() above. By
    // inserting any double value lying within the segment, we get as result a double value that is
    // larger than the sought-after segment ID but smaller than the next highter ID, i.e. we have to
    // round down the solution.
    unrounded_id = (point + 1.0) / 2.0 * numberofsegments;
    segment_id = floor(unrounded_id);

    if (segment_id < 0) dserror("Segment-ID can't be negative!");

    return segment_id;
  }

}  // namespace BEAMCONTACT

#endif
