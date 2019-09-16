/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of arbitrary lagrange polynomials in HDG context

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

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
    /*!
    \brief helper data type to evaluate only nonzero element shape functions on face

     Author: schoeder 10/14
    */
    struct ShapeValuesInteriorOnFace
    {
      // shape function
      void Shape(int numrows, int numcols)
      {
        matrix_.LightShape(numrows, numcols);
        isNonzero_.resize(numrows);
        return;
      }

      // for analog usage as a pure matrix
      double &operator()(const int row, const int col) { return matrix_(row, col); }

      // for analog usage as a pure matrix
      double operator()(const int row, const int col) const { return matrix_(row, col); }

      bool NonzeroOnFace(const int index) const { return isNonzero_[index]; }

      LINALG::SerialDenseMatrix matrix_;
      std::vector<bool> isNonzero_;

      // assignment operator
      ShapeValuesInteriorOnFace &operator=(ShapeValuesInteriorOnFace &other)
      {
        matrix_ = other.matrix_;
        isNonzero_ = other.isNonzero_;
        return *this;
      }
    };

    /*!
     \brief helper holding the parameters, which we need to construct a the ShapeValuesFace
    */
    struct ShapeValuesFaceParams
    {
      explicit ShapeValuesFaceParams(
          unsigned int degree, bool completepoly, unsigned int quadraturedegree)
          : degree_(degree),
            completepoly_(completepoly),
            quadraturedegree_(quadraturedegree),
            face_(0)
      {
      }

      unsigned int degree_;
      bool completepoly_;
      unsigned int quadraturedegree_;
      unsigned int face_;

      /// convert the data structure to integers to simplify comparisons
      /// (rather than implementing operator < with many if statements)
      std::size_t ToInt() const
      {
        dsassert(degree_ < 32, "Not implemented");
        dsassert(quadraturedegree_ < 64, "Not implemented");
        dsassert(face_ < 8, "Not implemented");
        // Simply encode the various integer values into a single number, all of
        // them take small numbers so this is easily possible. Beware about this when
        // adding further information, tough (might overflow when larger than 4 billion).
        return std::size_t(degree_) + 32 * completepoly_ + 64 * quadraturedegree_ + 4096 * face_;
      }
    };


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
      static const unsigned int nen_ =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      ///! number of faces on element
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

      ///! number of nodes on faces
      static const unsigned int nfn_ = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

      /*!
      \brief Initialize the data structure for shape functions, set element-independent data

      @param degree  Degree of the polynomial
      @param completepoly  Sets whether to use a polynomial of complete degree or tensor degree (on
      hex); on tet, only complete degree makes sense
      @param quadratureDegree  Sets the degree of polynomials the Gauss formula should integrate
      exactly; typically 2*degree
       */
      ShapeValues(
          const unsigned int degree, const bool completepoly, const unsigned int quadratureDegree);

      /*!
      \brief Compute element-dependent data, like gradients in real cell, integration weights, for
      element dofs
      */
      void Evaluate(const DRT::Element &ele, const std::vector<double> &aleDis = {});

      /// polynomial degree
      unsigned int degree_;

      /// underlying polynomial space for element interior, created in constructor
      Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_>> polySpace_;

      /// scalar dofs per element
      unsigned int ndofs_;

      /// quadrature rule
      Teuchos::RCP<DRT::UTILS::GaussPoints> quadrature_;

      /// complete poly
      bool usescompletepoly_;

      /// number of integration points
      unsigned int nqpoints_;

      LINALG::SerialDenseMatrix
          funct;  /// values of mapping shape functions on all quadrature points
      LINALG::Matrix<nsd_, nen_>
          deriv;                        /// gradients of mapping shape functions in unit coordinates
      LINALG::SerialDenseMatrix derxy;  /// gradients of mapping shape functions in real coordinates
                                        /// on all quadrature points
      LINALG::SerialDenseMatrix xyzreal;  /// coordinates of all quadrature points in real space
      LINALG::SerialDenseMatrix
          nodexyzunit;  /// coordinates of all node (support) points of Lagrange basis functions in
                        /// unit coordinates (all points at cell center for Legendre-type
                        /// polynomials)
      LINALG::SerialDenseMatrix
          nodexyzreal;  /// coordinates of all node (support) points of Lagrange basis functions in
                        /// real space (all points at cell center for Legendre-type polynomials)

      LINALG::SerialDenseMatrix shfunct;  /// evaluated HDG shape functions on all quadrature points
      LINALG::SerialDenseVector shfunctAvg;  /// average of shfunctF on cell
      LINALG::SerialDenseMatrix
          shderiv;  /// evaluated HDG shape function gradients in unit coordinates
      LINALG::SerialDenseMatrix
          shderxy;  /// evaluated HDG shape function gradients in real coordinates

      LINALG::Matrix<nsd_, 1> xsi;      /// quadrature points
      LINALG::Matrix<nsd_, nsd_> xjm;   /// Jacobi matrix of transformation
      LINALG::Matrix<nsd_, nsd_> xji;   /// inverse of Jacobi matrix of transformation
      LINALG::Matrix<nsd_, nen_> xyze;  /// element nodes
      LINALG::SerialDenseVector jfac;   /// Jacobian determinant times quadrature weight
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
      static const unsigned int nen_ =
          DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      ///! number of nodes on faces
      static const unsigned int nfn_ = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

      //! number of space dimensions
      static const unsigned int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      ///! number of faces on element
      static const unsigned int nfaces_ = DRT::UTILS::DisTypeToNumFaces<distype>::numFaces;

      /*!
      \brief Constructor which does the things which do not need to be redone for the same
      parameters
      */
      ShapeValuesFace(ShapeValuesFaceParams params);

      /*!
      \brief Compute element-dependent data on faces, like integration weights, normal vectors,
      correctly oriented trace variables
      */
      void EvaluateFace(
          const DRT::Element &ele, const unsigned int face, const std::vector<double> &aleDis = {});

      /*!
      \brief Consider the orientation of faces for face degrees of freedom
       */
      void AdjustFaceOrientation(const DRT::Element &ele, const unsigned int face);

      /// Parameters underlying this structure, necessary for evaluating element basis functions on
      /// faces
      ShapeValuesFaceParams params_;

      /// polynomial degree
      unsigned int degree_;

      /// underlying polynomial space for faces
      Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_ - 1>> polySpace_;

      /// scalar dofs per face
      unsigned int nfdofs_;

      /// quadrature rule
      Teuchos::RCP<DRT::UTILS::GaussPoints> quadrature_;

      /// number of integration points on face
      unsigned int nqpoints_;

      LINALG::SerialDenseMatrix
          funct;  /// values of mapping shape functions on all face quadrature points
      LINALG::Matrix<nsd_ - 1, nfn_> deriv;  /// gradients of mapping shape functions on face
      LINALG::Matrix<nsd_ - 1, nsd_ - 1> metricTensor;  /// metric tensor on face
      LINALG::Matrix<nsd_, 1> normal;                   /// normal vector
      LINALG::Matrix<nsd_, nsd_ - 1> tangent;  /// Face reference frame wrt real one (face -> real)

      LINALG::SerialDenseMatrix xyzreal;  /// coordinates of face quadrature points in real space
      LINALG::SerialDenseMatrix
          nodexyzunit;  /// coordinates of all node (support) points of face Lagrange basis
                        /// functions in unit coordinates (invalid for Legendre-type polynomials)
      LINALG::SerialDenseMatrix
          nodexyzreal;  /// coordinates of all node (support) points of face Lagrange basis
                        /// functions in real space (invalid for Legendre-type polynomials)

      LINALG::SerialDenseMatrix shfunct;  /// evaluated shape functions for HDG face polynomials,
                                          /// permuted to account for face orientation
      LINALG::SerialDenseMatrix shfunctNoPermute;  /// evaluated shape functions for HDG face
                                                   /// polynomials in natural ordering
      ShapeValuesInteriorOnFace
          shfunctI;  /// evaluated shape functions on face for interior HDG polynomials
      LINALG::SerialDenseMatrix
          normals;  /// normal vectors on a single face for all quadrature points

      LINALG::Matrix<nsd_ - 1, 1> xsi;  /// face quadrature points
      LINALG::Matrix<nsd_, nfn_> xyze;  /// face nodes
      LINALG::SerialDenseVector jfac;   /// face Jacobian determinant times quadrature weight

      std::vector<std::vector<int>> faceNodeOrder;  /// numbering of nodes belonging to faces

     private:
      LINALG::SerialDenseVector faceValues;  /// Evaluated basis functions on face

      /*!
      \brief Computes the face reference system considering the ordering of the master element

      The face reference system created is orthonormal and therefore the inverse of the
      transformationm matrix is just the transpose of the matrix (T)^{-1} = (T)'.

      \note
      The routine is only used in 3D so far and therefore tested only in this case.

      \Author: Berardocco
       */
      void ComputeFaceReferenceSystem(const DRT::Element &ele, const unsigned int face);
    };

    /*!
    \brief Cache for the face-based shape values, so we do not need to calculate again

     Author: schoeder 08/14
    */
    template <DRT::Element::DiscretizationType distype>
    class ShapeValuesFaceCache : public SingletonDestruction
    {
     public:
      /// return instance
      static ShapeValuesFaceCache<distype> &Instance();

      /// cleanup
      virtual void Done();

      /// give pointer to corresponding shape values face
      Teuchos::RCP<ShapeValuesFace<distype>> Create(ShapeValuesFaceParams params);

     private:
      /// instance
      static ShapeValuesFaceCache<distype> *instance_;

      /// the actual cache; contains all pairs
      std::map<std::size_t, Teuchos::RCP<ShapeValuesFace<distype>>> svf_cache_;
    };

    /*!
    \brief Cache for the interior values of shape velocities

     Author: kronbichler 09/14
    */
    template <DRT::Element::DiscretizationType distype>
    class ShapeValuesInteriorOnFaceCache : public SingletonDestruction
    {
     public:
      /// return instance
      static ShapeValuesInteriorOnFaceCache &Instance();

      /// cleanup
      virtual void Done();

      /// give pointer to corresponding shape values face
      Teuchos::RCP<ShapeValuesInteriorOnFace> Create(ShapeValuesFaceParams params);

     private:
      /// instance
      static ShapeValuesInteriorOnFaceCache *instance_;

      /// the actual cache; contains all pairs
      std::map<std::size_t, Teuchos::RCP<ShapeValuesInteriorOnFace>> cache_;
    };

  }  // namespace UTILS
}  // namespace DRT

#endif
