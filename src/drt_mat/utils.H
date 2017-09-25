/*----------------------------------------------------------------------*/
/*!
\file utils.H

\brief Utility file for the material evaluation
\level 3

\maintainer Michael Hiermeier

\date Sep 20, 2017

*/
/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_MAT_UTILS_H_
#define SRC_DRT_MAT_UTILS_H_

#include "../linalg/linalg_fixedsizematrix.H"

namespace MAT {

  /// Voigt notation types
  enum class Notation {
    stress,
    strain
  };

  /*--------------------------------------------------------------------------*/
  /// Index maps
  struct IMap
  {
    IMap() = default;

    /// map second order tensor indices to perturbed Voigt vector notation
    static const unsigned second_[3][3];

    /// map perturbed Voigt matrix notation for the fourth order elasticity
    /// tensor to the fourth order indices
    static const unsigned fourth_[6][6][4];
  };

  /*--------------------------------------------------------------------------*/
  /** \brief Utility routines for the perturbed Voigt tensor notation
   *
   *  \author hiermeier \date 09/17 */
  template < Notation type >
  class VoigtUtils
  {
  public:
    /// constructor
    VoigtUtils() = default;

    /// destructor
    ~VoigtUtils() = default;

    /** \brief compute power of a symmetric 3x3 matrix in perturbed Voigt notation
     *
     *  \f[
     *  [vtensorpow]_{AB} = [vtensor^{pow}]_{AB}
     *  \f]
     *
     *  \param[in] pow          positive integer exponent
     *  \param[in] vtensor      input tensor in Voigt <type> notation
     *  \param[in] vtensor_pow  result, i.e. input tensor to the given power
     *
     *  \author hiermeier \date 09/17 */
    static void PowerOfSymmetricTensor(
        const unsigned pow,
        const LINALG::Matrix<6,1>& vtensor,
        LINALG::Matrix<6,1>& vtensor_pow );

    /** \brief Compute the inverse tensor in perturbed Voigt notation
     *
     *  \param[in]  vtensor      tensor in Voigt <type> notation
     *  \param[out] vtensor_inv  inverse tensor in Voigt <type> notation
     *
     *  \author hiermeier \date 09/17 */
    static void InverseTensor(
        const LINALG::Matrix<6,1>& vtensor,
        LINALG::Matrix<6,1>& vtensor_inv );

    /** \brief Compute the product of a tensor in perturbed Voigt notation
     *  and a vector
     *
     *  \f$ [vecres]_{A} = vtensor_{AB} vec^{B} \f$
     *
     *  \param[in]  vtensor      tensor in Voigt <type> notation
     *  \param[out] vtensor_inv  inverser tensor in Voigt <type> notation
     *
     *  \author hiermeier \date 09/17 */
    static void MultiplyTensorVector(
        const LINALG::Matrix<6,1>& vtensor,
        const LINALG::Matrix<3,1>& vec,
        LINALG::Matrix<3,1>& vec_res );

    /** \brief Compute the symmetric outer product of two vectors
     *
     *  \f$ [abba]_{AB} = [veca]_{A} [vecb]_{B} + [veca]_{B} [vecb]_{A} \f$
     *
     *  \param[in]  vec_a  first vector
     *  \param[in]  vec_b  second vector
     *  \param[out] ab_ba  symmetric outer product of the two input vectors
     *                     in the Voigt <type> notation
     *
     *  \author hiermeier \date 09/17 */
    static void SymmetricOuterProduct(
        const LINALG::Matrix<3,1>& vec_a,
        const LINALG::Matrix<3,1>& vec_b,
        LINALG::Matrix<6,1>& ab_ba );

    /** \brief scale off diagonal values
     *
     *  \note This function changes the values only if the strain notation is used.
     *
     *  \param[out] tensor  scale the off-diagonal values of this tensor
     *
     *  \author hiermeier \date 09/17 */
    static void ScaleOffDiagonalVals( LINALG::Matrix<6,1>& tensor );

    /** \brief unscale off diagonal values
     *
     *  \note This function changes the values only if the strain notation is used.
     *
     *  \param[out] tensor  unscale the off-diagonal values of this tensor
     *
     *  \author hiermeier \date 09/17 */
    static void UnscaleOffDiagonalVals( LINALG::Matrix<6,1>& strain );

    /** \brief unscale factors for the perturbed Voigt strain notation
     *
     *  These factors are meaningful if the strain convention is followed. */
    static const double unscale_fac_[6];

    /** \brief scale factors for the perturbed Voigt stress notation
     *
     *  These factors are meaningful if the strain convention is followed. */
    static const double scale_fac_[6];
  };

  /// @name type defs for easier use
  /// @{

  /// Utility class for the perturbed Voigt stress notation
  typedef VoigtUtils<Notation::stress> VStressUtils;

  /// Utility class for the perturbed Voigt strain notation
  typedef VoigtUtils<Notation::strain> VStrainUtils;

  /// @}

  /// access scaling factors
  template< Notation type > inline double ScaleFactor( unsigned i );
  template <> inline double ScaleFactor<Notation::stress>( unsigned i )
  {
    return 1.0;
  };
  template <> inline double ScaleFactor<Notation::strain>( unsigned i )
  {
    return VStrainUtils::scale_fac_[i];
  };

  /// access unscaling factors
  template< Notation type > inline double UnscaleFactor( unsigned i );
  template <> inline double UnscaleFactor<Notation::stress>( unsigned i )
  {
    return 1.0;
  };
  template <> inline double UnscaleFactor<Notation::strain>( unsigned i )
  {
    return VStrainUtils::unscale_fac_[i];
  };

  /** \brief convert tensor from perturbed Voigt strain notation to
   *  perturbed Voigt stress notation
   *
   * \author hiermeier \date 09/17 */
  inline void ConvertStrainToStressNotation( LINALG::Matrix<6,1>& strain )
  {
    VStrainUtils::UnscaleOffDiagonalVals( strain );
  }

  /** \brief convert tensor from perturbed Voigt stress notation to
   *  perturbed Voigt strain notation
   *
   * \author hiermeier \date 09/17 */
  inline void ConvertStressToStrainNotation( LINALG::Matrix<6,1>& strain )
  {
    VStrainUtils::ScaleOffDiagonalVals( strain );
  }
}


#endif /* SRC_DRT_MAT_UTILS_H_ */
