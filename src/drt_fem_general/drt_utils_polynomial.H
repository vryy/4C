/*----------------------------------------------------------------------*/
/*! \file

 \brief Generic polynomials for HDG methods in 1D, 2D, 3D

\level 2

 */

#ifndef DRT_UTILS_POLYNOMIAL
#define DRT_UTILS_POLYNOMIAL

#include "drt_utils_local_connectivity_matrices.H"

#include "Epetra_SerialDenseSolver.h"

namespace DRT
{
  namespace UTILS
  {
    /*!
     \brief helper holding the parameters, which we need to know, to construct a polynomial space
     */
    struct PolynomialSpaceParams
    {
      explicit PolynomialSpaceParams(
          DRT::Element::DiscretizationType distype, unsigned int degree, bool completeSpace)
          : distype_(distype), degree_(degree), completeSpace_(completeSpace)
      {
      }

      DRT::Element::DiscretizationType distype_;
      unsigned int degree_;
      bool completeSpace_;

      bool operator<(const PolynomialSpaceParams &otherparams) const
      {
        if (degree_ < otherparams.degree_)
          return true;
        else if (degree_ == otherparams.degree_)
        {
          if ((int)distype_ < (int)(otherparams.distype_))
            return true;
          else if ((int)distype_ == (int)(otherparams.distype_))
          {
            if (completeSpace_ < otherparams.completeSpace_)
              return true;
            else
              return false;
          }
          else
            return false;
        }
        else
          return false;
      }
    };

    /*!
     \brief A class of 1D Lagrange polynomials on a given number of points
     */
    class LagrangePolynomial
    {
     public:
      /*!
       \brief Constructor.
       */
      LagrangePolynomial(const std::vector<double> &supportPoints, const double nodePoint)
          : supportPoints_(supportPoints), nodePoint_(nodePoint)
      {
        long double weight = 1.;
        for (unsigned int i = 0; i < supportPoints.size(); ++i)
          weight *= nodePoint - supportPoints[i];
        dsassert((std::abs(weight) > std::numeric_limits<double>::min() &&
                     std::abs(weight) < std::numeric_limits<double>::max()),
            "Error in evaluation of polynomial");
        weight_ = 1. / weight;
      }

      /*!
       \brief Copy constructor
       */
      LagrangePolynomial(const LagrangePolynomial &other)
          : supportPoints_(other.supportPoints_),
            nodePoint_(other.nodePoint_),
            weight_(other.weight_)
      {
      }

      /*!
       \brief Assignment operator
       */
      LagrangePolynomial &operator=(const LagrangePolynomial &other)
      {
        supportPoints_ = other.supportPoints_;
        nodePoint_ = other.nodePoint_;
        weight_ = other.weight_;
        return *this;
      }

      /*!
       \brief Evaluates the polynomial on a given point on the unit interval [-1, 1].
       */
      double Evaluate(const double point) const
      {
        double value = weight_;
        for (unsigned int i = 0; i < supportPoints_.size(); ++i) value *= point - supportPoints_[i];
        return value;
      }

      /*!
       \brief Evaluates the polynomial and its derivatives (up to the order given by
              the length of values) on the given point
       */
      template <typename M>
      void Evaluate(const double point, M &derivatives) const
      {
        dsassert(derivatives.Columns() == 1, "Only column vectors supported");
        if (derivatives.Rows() == 1)
        {
          derivatives(0) = Evaluate(point);
          return;
        }

        // compute the value and derivatives by expanding the derivatives
        derivatives(0) = weight_;
        for (unsigned int k = 1; k < derivatives.Rows(); ++k) derivatives(k) = 0;
        for (unsigned int i = 0; i < supportPoints_.size(); ++i)
        {
          const double v = point - supportPoints_[i];
          for (int k = derivatives.Rows() - 1; k > 0; --k)
            derivatives(k) = v * derivatives(k) + derivatives(k - 1);
          derivatives(0) *= v;
        }
        double faculty = 1;
        for (unsigned int k = 2; k < derivatives.Rows(); ++k)
        {
          faculty *= static_cast<double>(k);
          derivatives(k) *= faculty;
        }
      }

      double NodePoint() const { return nodePoint_; }

     private:
      std::vector<double> supportPoints_;
      double nodePoint_;
      double weight_;
    };



    /*!
     \brief A class of polynomials based on monomial coefficients

     This class takes a vector of coefficients of the monomials and evaluates a
     polynomial from these coefficients. The evaluation is done by the Horner scheme.
     */
    class Polynomial
    {
     public:
      /*!
       \brief Constructor.
       */
      Polynomial(const std::vector<double> &coefficients) : coefficients_(coefficients) {}

      /*!
       \brief Copy constructor
       */
      Polynomial(const Polynomial &other) : coefficients_(other.coefficients_) {}

      /*!
\brief Assignment operator
*/
      Polynomial &operator=(const Polynomial &other)
      {
        coefficients_ = other.coefficients_;
        return *this;
      }

      /*!
       \brief Evaluates the polynomial on a given point
       */
      double Evaluate(const double point) const
      {
        // classical Horner scheme
        const int length = coefficients_.size() - 1;
        double value = coefficients_[length];
        for (int i = length - 1; i >= 0; --i) value = value * point + coefficients_[i];
        return value;
      }

      /*!
       \brief Evaluates the polynomial and its derivatives (up to the order given by
           the length of values) on the given point
       */
      template <typename M>
      void Evaluate(const double point, M &derivatives) const
      {
        dsassert(derivatives.Columns() == 1, "Only column vectors supported");
        if (derivatives.Rows() == 1)
        {
          derivatives(0) = Evaluate(point);
          return;
        }

        // compute the value and derivatives by Horner scheme
        std::vector<double> temp(coefficients_);
        const int length = coefficients_.size();
        const int nderivs = derivatives.Rows();
        const int maxder = std::min(length, nderivs);
        double kfaculty = 1.;
        // skip derivatives that are necessarily zero
        for (int k = 0; k < maxder; ++k)
        {
          for (int i = length - 2; i >= k; --i) temp[i] += point * temp[i + 1];
          derivatives(k) = kfaculty * temp[k];
          kfaculty *= static_cast<double>(k + 1);
        }
        for (int k = maxder; k < nderivs; ++k) derivatives(k) = 0;
      }

      // Here, we do not know about nodes, so put the center of the 1d element
      double NodePoint() const { return 0.; }


     private:
      std::vector<double> coefficients_;
    };



    /*!
     \brief generates complete Lagrange basis in 1D for [-1,1]^d elements
     */
    std::vector<LagrangePolynomial> generateLagrangeBasis1D(const unsigned int degree);



    /*!
     \brief generates complete Legendre basis in 1D

     Legendre polynomials are orthogonal on the unit interval [-1, 1]. As opposed
     to the usual mathematical definition, we also scale the polynomials such that
     they are actually orthonormal on the unit interval.
     */
    std::vector<Polynomial> generateLegendreBasis1D(const unsigned int degree);



    /*!
    \brief Base class for polynomial spaces in nsd_ dimensions used by HDG
    */
    template <int nsd_>
    class PolynomialSpaceBase
    {
     public:
      virtual ~PolynomialSpaceBase() {}

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      virtual std::size_t Size() const = 0;

      /*
       \brief Evaluates the values of all polynomials on the given point
       */
      virtual void Evaluate(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseVector &values) const = 0;

      /*
       \brief Evaluates the values of all polynomials on the given point
       */
      virtual void Evaluate_deriv1(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const = 0;

      /*
       \brief Evaluates the first derivative of all polynomials on the given point
       */
      virtual void Evaluate_deriv2(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const = 0;

      /*
       \brief Creates an array with coordinates of the nodes supporting the polynomials.
       */
      virtual void FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const = 0;
    };



    /*!
     \brief A class of tensor product polynomials expanded from 1D polynomials

     The 1D polynomial class must feature two Evaluate functions, one taking just
     the point as an argument (evaluating the polynomial's value) and one taking
     the point and an array with as many elements as derivatives are requested.

     Base class for LagrangeBasis
     */
    template <int nsd_, class POLY>
    class PolynomialSpaceTensor : public PolynomialSpaceBase<nsd_>
    {
     public:
      /*
       \brief Constructor from a vector of one-dimensional polynomials
       */
      PolynomialSpaceTensor(const std::vector<POLY> polySpace1d) : polySpace1d_(polySpace1d)
      {
        // renumbering is included in case we want to use these polynomials also
        // for standard finite elements
        renumbering_.resize(Size());

        for (unsigned int i = 0; i < renumbering_.size(); ++i) renumbering_[i] = i;
      }

      virtual ~PolynomialSpaceTensor() {}

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      static std::size_t Size(const std::size_t degree)
      {
        std::size_t size = degree + 1;
        for (unsigned int d = 1; d < nsd_; ++d) size *= degree + 1;
        return size;
      }

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      virtual std::size_t Size() const { return Size(polySpace1d_.size() - 1); }

      /*
       \brief Evaluates the values of the whole polynomial space in the given point
       */
      virtual void Evaluate(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseVector &values) const;

      /*
       \brief Evaluates the first derivative of the whole polynomial space in the given point
       */
      virtual void Evaluate_deriv1(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const;

      /*
       \brief Evaluates the second derivative of the whole polynomial space in the given point
       */
      virtual void Evaluate_deriv2(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const;

      /*
       \brief Evaluates the second derivative of the whole polynomial space in the given point
       */
      template <typename M>
      void Evaluate_deriv2(const LINALG::Matrix<nsd_, 1> &point, M &derivatives) const
      {
      }

      /*
       \brief Convert from a index within the polynomial space to the tensor indices in the
       individual dimensions
       */
      LINALG::Matrix<nsd_, 1, unsigned int> getIndices(const unsigned int index) const
      {
        dsassert(index < Size(), "Access out of range");
        LINALG::Matrix<nsd_, 1, unsigned int> indices;
        const unsigned int npoly = polySpace1d_.size();
        switch (nsd_)
        {
          case 1:
            indices(0) = index;
            break;
          case 2:
            indices(0) = index % npoly;
            indices(1) = index / npoly;
            break;
          case 3:
            indices(0) = index % npoly;
            indices(1) = (index / npoly) % npoly;
            indices(2) = index / (npoly * npoly);
            break;
          default:
            dserror("Invalid dimension");
            break;
        }
        return indices;
      }

      /*
       \brief Creates an array with coordinates of the nodes supporting the polynomials.
       */
      virtual void FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const;

     private:
      std::vector<POLY> polySpace1d_;
      std::vector<int> renumbering_;
    };



    /*!
     \brief A class of polynomials of complete degree expanded from 1D polynomials by truncated
     tensor products

     The 1D polynomial class must feature two Evaluate functions, one taking just
     the point as an argument (evaluating the polynomial's value) and one taking
     the point and an array with as many elements as derivatives are requested.

     Base class for LagrangeBasis
     */
    template <int nsd_, class POLY>
    class PolynomialSpaceComplete : public PolynomialSpaceBase<nsd_>
    {
     public:
      /*
       \brief Constructor from a vector of one-dimensional polynomials
       */
      PolynomialSpaceComplete(const std::vector<POLY> polySpace1d) : polySpace1d_(polySpace1d)
      {
        renumbering_.resize(Size());
        for (unsigned int i = 0; i < renumbering_.size(); ++i) renumbering_[i] = i;
      }

      virtual ~PolynomialSpaceComplete() {}

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      static std::size_t Size(const std::size_t degree)
      {
        std::size_t size = degree + 1;
        for (unsigned int d = 1; d < nsd_; ++d)
        {
          size *= degree + 1 + d;
          size /= (d + 1);  // This integer division is always without remainder
        }
        return size;
      }

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      virtual std::size_t Size() const { return Size(polySpace1d_.size() - 1); }

      /*
       \brief Evaluates the values of the whole polynomial space in the given point
       */
      virtual void Evaluate(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseVector &values) const;

      /*
       \brief Evaluates the first derivative of the whole polynomial space in the given point
       */
      virtual void Evaluate_deriv1(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const;

      /*
       \brief Evaluates the second derivative of the whole polynomial space in the given point
       */
      virtual void Evaluate_deriv2(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const;

      /*
       \brief Creates an array with coordinates of the nodes supporting the polynomials.
       */
      void FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const;

     private:
      std::vector<POLY> polySpace1d_;
      std::vector<int> renumbering_;
    };



    /*!
     \brief Class of tensor product Lagrange polynomials based on 1D Gauss-Lobatto
     support points (gives standard Q1 and Q2 basis functions and good polynomials
     for higher order). This class is usually not directly called in user code, as
     its functionality can be accessed through PolynomialSpace<nsd_>.
    */
    template <int nsd_>
    class LagrangeBasis : public PolynomialSpaceTensor<nsd_, LagrangePolynomial>
    {
     public:
      /*
      \brief Constructor from a vector of one-dimensional polynomials
       */
      LagrangeBasis(const unsigned int degree)
          : PolynomialSpaceTensor<nsd_, LagrangePolynomial>(generateLagrangeBasis1D(degree))
      {
      }

      ~LagrangeBasis() {}
    };



    /*!
    \brief Class of Legendre polynomials that are pairwise orthogonal. This class
     is usually not directly called in user code, as its functionality can be
     accessed through PolynomialSpace<nsd_>.
    */
    template <int nsd_>
    class LegendreBasis : public PolynomialSpaceComplete<nsd_, Polynomial>
    {
     public:
      /*
     \brief Constructor from a vector of one-dimensional polynomials
       */
      LegendreBasis(const unsigned int degree)
          : PolynomialSpaceComplete<nsd_, Polynomial>(generateLegendreBasis1D(degree))
      {
      }
    };



    /*!
    \brief Lagrange basis for triangular/tetrahedral elements.

     This class constructs a Lagrangian basis for tetrahedral elements of arbitrary degree
     by a transformation of Legendre polynomials onto so-called Fekete nodes, which are
     a generalization of Gauss-Lobatto points to higher dimensions on triangles. For the
     transformation, a Vandermonde matrix is used to map the nodeless orthonormal Legendre
     polynomials onto the desired node distribution.

     The transformation with the Vandermonde matrix is described in
     R. Sevilla, S. Fernandez-Mendez, A. Huerta, NURBS-Enhanced Finite Element Method,
     Arch. Comput. Methods Eng. (2011), 441-484

     Author: kronbichler 08/14
     */
    template <int nsd_>
    class LagrangeBasisTet : public PolynomialSpaceBase<nsd_>
    {
     public:
      LagrangeBasisTet(const unsigned int degree) : legendre_(degree)
      {
        FillFeketePoints(degree);
        ComputeVandermondeMatrices(degree);
      }

      virtual ~LagrangeBasisTet() {}

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      static std::size_t Size(const std::size_t degree)
      {
        std::size_t size = degree + 1;
        for (unsigned int d = 1; d < nsd_; ++d)
        {
          size *= degree + 1 + d;
          size /= (d + 1);  // This integer division is always without remainder
        }
        return size;
      }

      /*
       \brief Return the number of polynomials (over all dimensions)
       */
      virtual std::size_t Size() const { return legendre_.Size(); }

      /*
       \brief Evaluates the values of the whole polynomial space in the given point
       */
      virtual void Evaluate(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseVector &values) const;

      /*
       \brief Evaluates the first derivative of the whole polynomial space in the given point
       */
      virtual void Evaluate_deriv1(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const;

      /*
       \brief Evaluates the second derivative of the whole polynomial space in the given point
       */
      virtual void Evaluate_deriv2(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const;

      /*
       \brief Creates an array with coordinates of the nodes supporting the polynomials.
       */
      virtual void FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const;

     private:
      void FillFeketePoints(const unsigned int degree);
      void ComputeVandermondeMatrices(const unsigned int degree);

      Epetra_SerialDenseMatrix vandermonde_;
      mutable Epetra_SerialDenseSolver vandermondeFactor_;
      mutable Epetra_SerialDenseMatrix evaluateVec_;
      Epetra_SerialDenseMatrix feketePoints_;
      LegendreBasis<nsd_> legendre_;
    };



    /*!
     \brief General polynomial class that encapsulates the various polynomial shapes

     For QUAD/HEX elements, this class allows to select between the usual tensor
     polynomial space LagrangeBasisHEX and the complete polynomial space
     LegendreBasis. For TRI/TET elements, Legendre polynomials are always used.
     We might want to implement a Lagrange basis for triangles at some point but
     for that we need a special transform to be able to use the above truncated
     tensor product.
     */
    template <int nsd_>
    class PolynomialSpace
    {
     public:
      PolynomialSpace(const DRT::Element::DiscretizationType distype, const unsigned int degree,
          const bool completeSpace)
          : polyspace_((getNumberOfElementFaces(distype) == 1 + nsd_ && nsd_ > 1)
                           ? static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_> *>(
                                 new DRT::UTILS::LagrangeBasisTet<nsd_>(degree))
                           : completeSpace ? static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_> *>(
                                                 new DRT::UTILS::LegendreBasis<nsd_>(degree))
                                           : static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_> *>(
                                                 new DRT::UTILS::LagrangeBasis<nsd_>(degree)))
      {
        if (nsd_ != getDimension(distype))
          dserror("Dimension of shape does not match template argument nsd_ in PolynomialSpace");
      }

      PolynomialSpace(PolynomialSpaceParams params)
          : polyspace_((getNumberOfElementFaces(params.distype_) == 1 + nsd_ && nsd_ > 1)
                           ? static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_> *>(
                                 new DRT::UTILS::LagrangeBasisTet<nsd_>(params.degree_))
                           : params.completeSpace_
                                 ? static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_> *>(
                                       new DRT::UTILS::LegendreBasis<nsd_>(params.degree_))
                                 : static_cast<DRT::UTILS::PolynomialSpaceBase<nsd_> *>(
                                       new DRT::UTILS::LagrangeBasis<nsd_>(params.degree_)))
      {
      }

      /*
     \brief Return the number of polynomials (over all dimensions)
       */
      std::size_t Size() const { return polyspace_->Size(); }

      /*
     \brief Evaluates the values of all polynomials on the given point
       */
      void Evaluate(const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseVector &values) const
      {
        polyspace_->Evaluate(point, values);
      }

      /*
     \brief Evaluates the values of all polynomials on the given point
       */
      void Evaluate_deriv1(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const
      {
        polyspace_->Evaluate_deriv1(point, derivatives);
      }

      /*
     \brief Evaluates the first derivative of all polynomials on the given point
       */
      void Evaluate_deriv2(
          const LINALG::Matrix<nsd_, 1> &point, Epetra_SerialDenseMatrix &derivatives) const
      {
        polyspace_->Evaluate_deriv2(point, derivatives);
      }

      /*
     \brief Creates an array with coordinates of the nodes supporting the polynomials.

     For Lagrange polynomials, these are the node points of the respective Gauss-Lobatto
     quadrature formula. For Legendre polynomials which are non-nodal, all points will be
     the element center coordinates.

     The output matrix is resized to the correct dimension.
       */
      void FillUnitNodePoints(LINALG::SerialDenseMatrix &matrix) const
      {
        polyspace_->FillUnitNodePoints(matrix);
      }

     private:
      Teuchos::RCP<PolynomialSpaceBase<nsd_>> polyspace_;
    };

    /*!
    \brief Cache for polynomial space so we do not need to calculate again

     In analogy to GaussPointCache
     Author: schoeder 06/14
     */
    template <int nsd_>
    class PolynomialSpaceCache : public DRT::SingletonDestruction
    {
     public:
      static PolynomialSpaceCache<nsd_> &Instance();

      /// cleanup
      virtual void Done();

      Teuchos::RCP<PolynomialSpace<nsd_>> Create(PolynomialSpaceParams params);

     private:
      static PolynomialSpaceCache<nsd_> *instance_;

      /// cache of already created polynomial spaces
      std::map<PolynomialSpaceParams, Teuchos::RCP<PolynomialSpace<nsd_>>> ps_cache_;
    };

    /*!
     \brief Returns the size of the polynomial space.

     Note: This class must always be synchronized with PolynomialSpace above!
     */
    inline int getBasisSize(
        const DRT::Element::DiscretizationType distype, const int degree, const bool completeSpace)
    {
      const int dim = getDimension(distype);
      const int nfaces = getNumberOfElementFaces(distype);
      switch (dim)
      {
        case 3:
          if (completeSpace || (nfaces == dim + 1))
            return LegendreBasis<3>::Size(degree);
          else
            return LagrangeBasis<3>::Size(degree);
          break;
        case 2:
          if (completeSpace || (nfaces == dim + 1))
            return LegendreBasis<2>::Size(degree);
          else
            return LagrangeBasis<2>::Size(degree);
          break;
        case 1:
          return degree + 1;
          break;
        default:
          dserror("Invalid dimension");
          return -1;
      }
    }

  }  // namespace UTILS
}  // namespace DRT

#endif  // #ifdef DRT_UTILS_POLYNOMIAL
