/*!-----------------------------------------------------------------------------------------------*
\file direct_divergence_refplane.H

\brief Construct reference plane for direct divergence method when used in global
coordinate system

\level 2

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifndef DIRECT_DIVERGENCE_REFPLANE_H_
#define DIRECT_DIVERGENCE_REFPLANE_H_

#include <vector>
#include "cut_element.H"

namespace GEO {
namespace CUT {

class Point;

// When direct divergence is used in local coordinate system, a reference plane can be easily chosen
// to ensure that all the integration points are completely within the background element. However,
// ALE methods requires the integration to be constructed in global coordinate system in which
// obtaining a proper reference plane is not trivial.
// This class handles all possible choices of choosing a reference plane in such a way that all
// integration points are ensured to be within the background element

class DirectDivergenceGlobalRefplane {

public:
  DirectDivergenceGlobalRefplane( Element* elem )
  :elem1_(elem)
  {
  }

  ~DirectDivergenceGlobalRefplane(){};

  /*!
  \Compute the reference plane for this element
   */
  std::vector<double> GetReferencePlane();

  /*!
  \brief Get the reference points that are used to define the reference plane. This is used in
  gmsh output of volume cells
   */
  std::vector<Point*> GetReferencePointGmsh(){ return refPtsGmsh_; }

private:

  /*!
  \brief Compute reference plane based on the diagonals of the element
   */
  bool DiagonalBasedRef( std::vector<double>& RefPlaneEqn );

  /*!
  \brief Compute reference plane based on the sides of the element
   */
  bool SideBasedRef( std::vector<double>& RefPlaneEqn );

  /*!
  \brief Returns true if all the projected points are within the element
   */
  bool isAllProjectedCornersInsideEle( std::vector<double>& RefPlaneEqn );

  /*!
  \brief Scale the given equation of plane
   */
  void scaleEquationOfPlane( std::vector<double>& RefPlaneEqn );

  /*!
  \brief Class used to sort maps in descending order
   */
  struct compareClass
  {
    bool operator()(const double &left, const double &right) const
    {
        return left > right;
    }
  };

  //! background element that contains the volumecell
  Element* elem1_;

  //! Points that define the reference plane used for this volumecell
  std::vector<Point*> refPtsGmsh_;
};

} /* namespace CUT */
} /* namespace GEO */
#endif /* DIRECT_DIVERGENCE_REFPLANE_H_ */
