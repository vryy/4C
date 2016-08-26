/*!----------------------------------------------------------------------
\file mortar_projector.H

\brief header for projector functions

\level 1

\maintainer Philipp Farah, Alexander Seitz

*----------------------------------------------------------------------*/
#ifndef MORTAR_PROJECTOR_H
#define MORTAR_PROJECTOR_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

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

class MortarProjector : DRT::SingletonDestruction
{
public:

  MortarProjector() {};

  virtual ~MortarProjector() {}

  //! @name Access methods

  /// 1. Internal implementation class
  static MortarProjector* Impl(MortarElement& ele);

  /// 2. Internal implementation class
  static MortarProjector* Impl(MortarElement& sele, MortarElement& mele);

  //! @name virtual functions
  virtual bool ProjectNodalNormal(MORTAR::MortarNode& node,
                          MORTAR::MortarElement& ele,
                          double* xi) = 0;

  virtual bool ProjectElementNormal(MORTAR::MortarNode& node,
                            MORTAR::MortarElement& ele,
                            double* xi) = 0;

  virtual bool ProjectGaussPoint2D(MORTAR::MortarElement& gpele,
                         const double* gpeta,
                         MORTAR::MortarElement& ele,
                         double* xi) = 0;

  virtual bool ProjectGaussPoint3D(MORTAR::MortarElement& gpele,
                           const double* gpeta,
                           MORTAR::MortarElement& ele,
                           double* xi, double& par) = 0;

  virtual bool ProjectGaussPointAuxn3D(const double* globgp,
                               const double* auxn,
                               MORTAR::MortarElement& ele,
                               double* xi, double& par) = 0;

  virtual bool ProjectGaussPointHermit(MORTAR::MortarElement& gpele,
                                       const double* gpeta, MORTAR::MortarElement& ele,
                                       double* xi) = 0;

  virtual bool ProjectSNodeByMNormal(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist) = 0;


  virtual bool ProjectSNodeByMNodalNormalLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin) = 0;

  virtual bool ProjectSNodeByMNormalLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin) = 0;
}; // class MortarProjector

/*!
\brief A class to perform projections of nodes onto opposing elements

\author farah (farah@lnm.mw.tum.de)
*/
template<DRT::Element::DiscretizationType distype>
class MortarProjectorCalc : public MortarProjector
{

public:

  // constructor
  MortarProjectorCalc();

  /// called upon destruction
  virtual void Done();

  /// Singleton access method
  static MortarProjectorCalc<distype> * Instance(bool create=true);

  //! n_: number of element nodes
  static const int n_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions ("+1" due to considering only interface elements)
  static const int ndim_ = DRT::UTILS::DisTypeToDim<distype>::dim + 1;

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
  \brief Project a Gauss point onto an element along GP normal

  Used to project a slave side GP onto an element of the master side

  This method will compute the coordinates of a projection of a Gauss
  point in the local coordinate system of an element.

  \param gpele (in): Slave element containing GP to project
  \param gpeta (in): Local coordinates of GP on gpele
  \param ele (in)  : Master element to project on
  \param xi (out)  : Local coordinates of projection on master element

  */
  bool ProjectGaussPoint2D(MORTAR::MortarElement& gpele,
                         const double* gpeta,
                         MORTAR::MortarElement& ele,
                         double* xi)
  {
    dserror("Called ele-based projection for segment-based integration!!!");
    return false;
  };

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
                           double* xi, double& par)
  {
    dserror("Called ele-based projection for segment-based integration!!!");
    return false;
  };

  /*!
  \brief Project a Gauss point onto an element along GP normal for Hermit

  Used to project a slave side GP onto an element of the master side

  This method will compute the coordinates of a projection of a Gauss
  point in the local coordinate system of an element.

  \param gpele (in): Slave element containing GP to project
  \param gpeta (in): Local coordinates of GP on gpele
  \param ele (in)  : Master element to project on
  \param xi (out)  : Local coordinates of projection on master element

  */
  bool ProjectGaussPointHermit(MORTAR::MortarElement& gpele,
                               const double* gpeta, MORTAR::MortarElement& ele,
                               double* xi)
  {
    dserror("Called ele-based projection for segment-based integration!!!");
    return false;
  };

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

  //TODO explanation
  bool ProjectSNodeByMNormal(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist);

  //TODO explanation
  bool ProjectSNodeByMNodalNormalLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  //TODO explanation
  bool ProjectSNodeByMNormalLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

protected:

  bool ProjectSNodeByMNormal3D(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist);

  bool ProjectSNodeByMNormal3DLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  bool ProjectSNodeByMNormal2D(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist);

  bool ProjectSNodeByMNodalNormal2DLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  bool ProjectSNodeByMNodalNormal3DLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  bool ProjectSNodeByMNormal2DLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

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
};

/*!
\brief A class to perform element based projections of nodes onto opposing elements

\author farah (farah@lnm.mw.tum.de)
*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
class MortarProjectorCalc_EleBased : public MortarProjector
{

public:

  // constructor
  MortarProjectorCalc_EleBased();

  /// called upon destruction
  virtual void Done();

  /// Singleton access method
  static MortarProjectorCalc_EleBased<distypeS, distypeM> * Instance(bool create=true);

  //! ns_: number of slave element nodes
  static const int ns_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;

  //! nm_: number of master element nodes
  static const int nm_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

  //! number of space dimensions ("+1" due to considering only interface elements)
  static const int ndim_ = DRT::UTILS::DisTypeToDim<distypeS>::dim + 1;

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
                          double* xi)
  {
    dserror("Called segment-based projection for element-based integration!!!");
    return false;
  };

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
                            double* xi)
  {
    dserror("Called segment-based projection for element-based integration!!!");
    return false;
  };

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
  bool ProjectGaussPoint2D(MORTAR::MortarElement& gpele,
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
                               double* xi, double& par)
  {
    dserror("Called Aux.-plane projection for element-based integration!!!");
    return false;
  };

  /*!
  \brief Project a Gauss point onto an element along GP normal for Hermit

  Used to project a slave side GP onto an element of the master side

  This method will compute the coordinates of a projection of a Gauss
  point in the local coordinate system of an element.

  \param gpele (in): Slave element containing GP to project
  \param gpeta (in): Local coordinates of GP on gpele
  \param ele (in)  : Master element to project on
  \param xi (out)  : Local coordinates of projection on master element

  */
  bool ProjectGaussPointHermit(MORTAR::MortarElement& gpele,
                               const double* gpeta, MORTAR::MortarElement& ele,
                               double* xi);

  bool ProjectSNodeByMNormal(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist){dserror("ERROR");return false;};

  bool ProjectSNodeByMNodalNormalLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin){dserror("ERROR");return false;};

  bool ProjectSNodeByMNormalLin(MORTAR::MortarNode& snode,
                         MORTAR::MortarElement& mele,
                         double* xi,
                         double* normal,
                         double& dist,
                         std::vector<GEN::pairedvector<int,double> >& normaltolineLin){dserror("ERROR");return false;};

protected:

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
  \brief Evaluate F for Gauss point projection

  */
  double EvaluateFGaussPoint2D(const double* gpx,
                             const double* gpn,
                             MORTAR::MortarElement& ele,
                             const double* eta);

  /*!
  \brief Evaluate GradF for Gauss point projection

  */
  double EvaluateGradFGaussPoint2D(const double* gpn,
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
  \brief Evaluate F for Gauss point projection (2D Hermit)

  */
  double EvaluateFGaussPointHermit(const double* gpx,
                                   const double* gpn,
                                   MORTAR::MortarElement& ele,
                                   const double* eta);

  /*!
  \brief Evaluate GradF for Gauss point projection (2D Hermit)

  */
  double EvaluateGradFGaussPointHermit(const double* gpn,
                                       MORTAR::MortarElement& ele,
                                       const double* eta);
};
}  // namespace MORTAR


#endif  // #ifndef MORTAR_PROJECTOR_H
