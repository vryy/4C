/*!----------------------------------------------------------------------
\file drt_utils_shapevalues_hdg.cpp
\brief

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_UTILS_SHAPEVALUES_HDG_H
#define DRT_UTILS_SHAPEVALUES_HDG_H

#include "drt_utils_polynomial.H"
#include "drt_utils_gausspoints.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

namespace DRT
{
  namespace UTILS
  {

    /// Helper class for evaluating HDG polynomials, geometry, etc.
    /*!

      \author kronbichler
      \date 05/14
    */
    template <DRT::Element::DiscretizationType distype>
    class ShapeValues
    {
    public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const unsigned int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      ///! number of faces on element
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

      ///! number of nodes on faces
      static const unsigned int nfn_    = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

      /*!
      \brief Initialize the data structure for shape functions, set element-independent data

      @param degree  Degree of the polynomial
      @param completepoly  Sets whether to use a polynomial of complete degree or tensor degree (on hex); on tet, only complete degree makes sense
      @param quadratureDegree  Sets the degree of polynomials the Gauss formula should integrate exactly; typically 2*degree
       */
      ShapeValues (const unsigned int degree,
                   const bool         completepoly,
                   const unsigned int quadratureDegree);

      /*!
      \brief Compute element-dependent data, like gradients in real cell, integration weights, for element dofs
      */
      void Evaluate              (const DRT::Element &ele);

      /// polynomial degree
      unsigned int degree_;

      /// underlying polynomial space for element interior, created in constructor
      Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_> >  polySpace_;

      /// scalar dofs per element
      unsigned int ndofs_;

      /// quadrature rule
      Teuchos::RCP<DRT::UTILS::GaussPoints> quadrature_;

      /// complete poly
      bool usescompletepoly_;

      /// number of integration points
      unsigned int nqpoints_;

      LINALG::SerialDenseMatrix     funct;   /// values of mapping shape functions on all quadrature points
      LINALG::Matrix<nsd_,nen_>     deriv;   /// gradients of mapping shape functions
      LINALG::SerialDenseMatrix     xyzreal; /// coordinates of all quadrature points in real space
      LINALG::SerialDenseMatrix     nodexyzunit; /// coordinates of all node (support) points of Lagrange basis functions in unit coordinates (all points at cell center for Legendre-type polynomials)
      LINALG::SerialDenseMatrix     nodexyzreal; /// coordinates of all node (support) points of Lagrange basis functions in real space (all points at cell center for Legendre-type polynomials)

      LINALG::SerialDenseMatrix     shfunct; /// evaluated HDG shape functions on all quadrature points
      LINALG::SerialDenseVector     shfunctAvg;/// average of shfunctF on cell
      LINALG::SerialDenseMatrix     shderiv; /// evaluated HDG shape function gradients in unit coordinates
      LINALG::SerialDenseMatrix     shderxy; /// evaluated HDG shape function gradients in real coordinates

      LINALG::Matrix<nsd_,1>        xsi;     /// quadrature points
      LINALG::Matrix<nsd_,nsd_>     xjm;     /// Jacobi matrix of transformation
      LINALG::Matrix<nsd_,nsd_>     xji;     /// inverse of Jacobi matrix of transformation
      LINALG::Matrix<nsd_,nen_>     xyze;    /// element nodes
      LINALG::SerialDenseVector     jfac;    /// Jacobian determinant times quadrature weight
    };

    /// Helper class for evaluating HDG polynomials, geometry, etc.
    /*!

      \author schoeder
      \date 06/14
    */
    template <DRT::Element::DiscretizationType distype>
    class ShapeValuesFace
    {
    public:
      //! nen_: number of element nodes (T. Hughes: The Finite Element Method)
      static const unsigned int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      ///! number of nodes on faces
      static const unsigned int nfn_    = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

      //! number of space dimensions
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      ///! number of faces on element
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

      /*!
      \brief Constructor. Does nothing because of p-adaptivity where EvaluateFace needs to create all the data fields
       */
      ShapeValuesFace();

      /*!
      \brief Compute element-dependent data on faces, like integration weights, normal vectors, correctly oriented trace variables
      */
      void EvaluateFace(const unsigned int degree,
                        const bool         completepoly,
                        const unsigned int quadratureDegree,
                        const DRT::Element &ele,
                        const unsigned int  face,
                        const DRT::UTILS::ShapeValues<distype> shapes);

      /*!
      \brief Consider the orientation of faces for face degrees of freedom
       */
      void AdjustFaceOrientation (const DRT::Element &ele,
                                  const unsigned int  face);

      /// polynomial degree
      unsigned int degree_;

      /// underlying polynomial space for faces
      Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_-1> > polySpace_;

      /// scalar dofs per face
      unsigned int nfdofs_;

      /// quadrature rule
      Teuchos::RCP<DRT::UTILS::GaussPoints> quadrature_;

      /// number of integration points on face
      unsigned int nqpoints_;

      LINALG::SerialDenseMatrix     funct;  /// values of mapping shape functions on all face quadrature points
      LINALG::Matrix<nsd_-1,nfn_>   deriv;  /// gradients of mapping shape functions on face
      LINALG::Matrix<nsd_-1,nsd_-1> metricTensor; /// metric tensor on face
      LINALG::Matrix<nsd_,1>        normal;  /// normal vector

      LINALG::SerialDenseMatrix     xyzreal;/// coordinates of face quadrature points in real space
      LINALG::SerialDenseMatrix     nodexyzunit; /// coordinates of all node (support) points of face Lagrange basis functions in unit coordinates (invalid for Legendre-type polynomials)
      LINALG::SerialDenseMatrix     nodexyzreal; /// coordinates of all node (support) points of face Lagrange basis functions in real space (invalid for Legendre-type polynomials)

      LINALG::SerialDenseMatrix     shfunct;/// evaluated shape functions for HDG face polynomials, permuted to account for face orientation
      LINALG::SerialDenseMatrix     shfunctNoPermute;/// evaluated shape functions for HDG face polynomials in natural ordering
      std::vector<LINALG::SerialDenseMatrix> shfunctI;/// evaluated shape functions on face for interior HDG polynomials
      LINALG::SerialDenseMatrix     normals; /// normal vectors on a single face for all quadrature points

      LINALG::Matrix<nsd_-1,1>      xsi;    /// face quadrature points
      LINALG::Matrix<nsd_,nfn_>     xyze;   /// face nodes
      LINALG::SerialDenseVector     jfac;   /// face Jacobian determinant times quadrature weight

      std::vector<std::vector<int> > faceNodeOrder; /// numbering of nodes belonging to faces

    private:
      LINALG::SerialDenseMatrix     faceQPoints; /// Quadrature points projected onto face
      LINALG::SerialDenseVector     faceValues; /// Evaluated basis functions on face
    };
  }
}

#endif

