/*!----------------------------------------------------------------------
\file mortar_projector.H

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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_PROJECTOR_H
#define MORTAR_PROJECTOR_H

#include "../linalg/linalg_fixedsizematrix.H"


/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations
class MortarNode;
class MortarElement;

/*!
\brief A class to perform projections of nodes onto opposing elements

\author popp (popp@lnm.mw.tum.de)
*/

class MortarProjector
{
public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as projections are
  performed in parallel by individual processes.

  \param dim (in): Global problem dimension

  */
  MortarProjector(int dim);

  /*!
  \brief Destructor

  */
  virtual ~MortarProjector() {}

  //! @name Access methods

  /*!
  \brief Get problem dimension

  */
  virtual const int& Dim() {return dim_;};

  //@}

  //! @name 2D and 3D projection methods

  /*!
  \brief Project a node onto an element along the node's normal

  Used to project a slave side node onto an element of the master side

  This method will compute the coordinates of a projection of a node in
  the local coordinate system of an element. The projection point will
  not neccesarily fall inside the element. However, if the projection
  point is far outside the segment's boundaries, problems with the
  internal nonlinear iteration might occur and a warning is issued when
  convergence can not be achieved in a limited number of iterations.

  \param node (in): Slave node to project
  \param ele (in) : Master element to project on
  \param xi (out) : Local coordinates of projection on element

  */
  bool ProjectNodalNormal(MORTAR::MortarNode& node,
                          MORTAR::MortarElement& ele,
                          double* xi);

  /*!
  \brief Project a node onto an element along the interpolated
         outward normal field of the element

  Used to project a master side node onto an element of the slave side

  This method will compute the coordinates of a projection of a node in
  the local coordinate system of an element. The projection point will
  not neccesarily fall inside the element. However, if the projection
  point is far outside the segment's boundaries, problems with the
  internal nonlinear iteration might occur and a warning is issued when
  convergence can not be achieved in a limited number of iterations.

  \param node (in): Master node to project
  \param ele (in) : Slave element to project on
  \param xi (out) : Local coordinates of projection on element

  */
  bool ProjectElementNormal(MORTAR::MortarNode& node,
                            MORTAR::MortarElement& ele,
                            double* xi);

  /*!
  \brief Project a node onto an element along the interpolated
         outward normal field of the element (3D)

  Used to project a master side node onto an element of the slave side (3D)

  This method will compute the coordinates of a projection of a node in
  the local coordinate system of an element. The projection point will
  not neccesarily fall inside the element. However, if the projection
  point is far outside the segment's boundaries, problems with the
  internal nonlinear iteration might occur.

  \param node (in): Master node to project
  \param ele (in) : Slave element to project on
  \param xi (out) : Local coordinates of projection on element
  \param par (out): Projection parameter alpha

  */
  bool ProjectElementNormal3D(MORTAR::MortarNode& node,
                              MORTAR::MortarElement& ele,
                              double* xi, double& par);

  /*!
  \brief Project a Gauss point onto an element along GP normal

  Used to project a slave side GP onto an element of the master side

  This method will compute the coordinates of a projection of a Gauss
  point in the local coordinate system of an element.

  \param gpele (in): Slave element containing GP to project
  \param gpeta (in): Local coordinates of GP on gpele
  \param ele (in)  : Master element to project on
  \param xi (out)  : Local coordinates of projection on master element

  */
  bool ProjectGaussPoint(MORTAR::MortarElement& gpele,
                         const double* gpeta,
                         MORTAR::MortarElement& ele,
                         double* xi);

  /*!
  \brief Project a Gauss point onto an element along GP normal (3D)

  Used to project a slave side GP onto an element of the master side

  This method will compute the coordinates of a projection of a Gauss
  point in the local coordinate system of an element.

  \param gpele (in): Slave element containing GP to project
  \param gpeta (in): Local coordinates of GP on gpele
  \param ele (in)  : Master element to project on
  \param xi (out)  : Local coordinates of projection on master element
  \param par (out ): Projection parameter alpha

  */
  bool ProjectGaussPoint3D(MORTAR::MortarElement& gpele,
                           const double* gpeta,
                           MORTAR::MortarElement& ele,
                           double* xi, double& par);

  /*!
  \brief Check intersection of projection normal with warped element to project on.

  Create aux-plane for every ele-node with adjacent element edges.
  If all intersection points of the projection normal and the aux-planes
  are located inside the corresponding ele-edges, then the gp lie on the element.
  --> return false for projection check

  \param ele(in)    : element to project on
  \param ngp(in)    : gp-normal to project along
  \param globgp(in) : global gp coordinates

  */
  bool CheckProjection4AUXPLANE(MORTAR::MortarElement& ele,
                       double* ngp, double* globgp);

  /*!
  \brief Project a Gauss point onto an element along AuxPlane normal (3D)

  Used to project an AuxPlane GP onto an element of the slave or master side

  This method will compute the coordinates of a projection of a Gauss
  point in the local coordinate system of an element.

  \param globgp(in): Gauss point to project, given in global coords
  \param auxn(in)  : Normal of AuxPlane along which to project
  \param ele (in)  : Slave / master element to project on
  \param xi (out)  : Local coordinates of projection on element
  \param par (out ): Projection parameter alpha

  */
  bool ProjectGaussPointAuxn3D(const double* globgp,
                               const double* auxn,
                               MORTAR::MortarElement& ele,
                               double* xi, double& par);

  /*!
  \brief Evaluate F for nodal normal projection

  */
  double EvaluateFNodalNormal(MORTAR::MortarNode& node,
                              MORTAR::MortarElement& ele,
                              const double* eta);

  /*!
  \brief Evaluate GradF for nodal normal projection

  */
  double EvaluateGradFNodalNormal(MORTAR::MortarNode& node,
                                  MORTAR::MortarElement& ele,
                                  const double* eta);

  /*!
  \brief Evaluate F for element normal projection

  */
  double EvaluateFElementNormal(MORTAR::MortarNode& node,
                                MORTAR::MortarElement& ele,
                                const double* eta);

  /*!
  \brief Evaluate GradF for element normal projection

  */
  double EvaluateGradFElementNormal(MORTAR::MortarNode& node,
                                    MORTAR::MortarElement& ele,
                                    const double* eta);

  /*!
  \brief Evaluate F for element normal projection (3D)

  */
  bool EvaluateFElementNormal3D(double* f,
                                MORTAR::MortarNode& node,
                                MORTAR::MortarElement& ele,
                                const double* eta,
                                const double& alpha);

  /*!
  \brief Evaluate GradF for element normal projection (3D)

  */
  bool EvaluateGradFElementNormal3D(LINALG::Matrix<3,3>& fgrad,
                              MORTAR::MortarNode& node,
                              MORTAR::MortarElement& ele,
                              const double* eta,
                              const double& alpha);

  /*!
  \brief Evaluate F for Gauss point projection

  */
  double EvaluateFGaussPoint(const double* gpx,
                             const double* gpn,
                             MORTAR::MortarElement& ele,
                             const double* eta);

  /*!
  \brief Evaluate GradF for Gauss point projection

  */
  double EvaluateGradFGaussPoint(const double* gpn,
                                 MORTAR::MortarElement& ele,
                                 const double* eta);

  /*!
  \brief Evaluate F for Gauss point projection (3D)

  */
  bool EvaluateFGaussPoint3D(double* f,
                             const double* gpx,
                             const double* gpn,
                             MORTAR::MortarElement& ele,
                             const double* eta,
                             const double& alpha);

  /*!
  \brief Evaluate GradF for Gauss point projection (3D)

  */
  bool EvaluateGradFGaussPoint3D(LINALG::Matrix<3,3>& fgrad,
                                 const double* gpx,
                                 const double* gpn,
                                 MORTAR::MortarElement& ele,
                                 const double* eta,
                                 const double& alpha);

  /*!
  \brief Evaluate F for AuxPlane Gauss point projection (3D)

  */
  bool EvaluateFGaussPointAuxn3D(double* f,
                                 const double* globgp,
                                 const double* auxn,
                                 MORTAR::MortarElement& ele,
                                 const double* eta,
                                 const double& alpha);

  /*!
  \brief Evaluate GradF for AuxPlane Gauss point projection (3D)

  */
  bool EvaluateGradFGaussPointAuxn3D(LINALG::Matrix<3,3>& fgrad,
                                     const double* globgp,
                                     const double* auxn,
                                     MORTAR::MortarElement& ele,
                                     const double* eta,
                                     const double& alpha);

  //@}

protected:

  // don't want = operator and cctor
  MortarProjector operator = (const MortarProjector& old);
  MortarProjector(const MortarProjector& old);

  int       dim_;          // dimension of problem (2D or 3D)

}; // class MortarProjector
}  // namespace MORTAR


#endif  // #ifndef MORTAR_PROJECTOR_H
