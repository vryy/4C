/*! \file
\brief Voigt notation definition and utilities
\level 1
\maintainer Sebastian Proell
*/

#include "voigt_notation.H"
#include <Sacado.hpp>

using NotationType = UTILS::VOIGT::NotationType;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void UTILS::VOIGT::Matrix3x3to9x1(LINALG::Matrix<3, 3> const& in, LINALG::Matrix<9, 1>& out)
{
  for (int i = 0; i < 3; ++i) out(i) = in(i, i);
  out(3) = in(0, 1);
  out(4) = in(1, 2);
  out(5) = in(0, 2);
  out(6) = in(1, 0);
  out(7) = in(2, 1);
  out(8) = in(2, 0);
}

template <NotationType rows_notation, NotationType cols_notation>
void UTILS::VOIGT::FourthOrderIdentityMatrix(LINALG::Matrix<6, 6>& id)
{
  id.Clear();

  for (unsigned int i = 0; i < 3; ++i) id(i, i) = 1.0;

  for (unsigned int i = 3; i < 6; ++i)
    id(i, i) =
        0.5 * VoigtUtils<rows_notation>::ScaleFactor(i) * VoigtUtils<cols_notation>::ScaleFactor(i);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::SymmetricOuterProduct(const LINALG::Matrix<3, 1>& vec_a,
    const LINALG::Matrix<3, 1>& vec_b, LINALG::Matrix<6, 1>& ab_ba)
{
  std::fill(ab_ba.A(), ab_ba.A() + 6, 0.0);

  LINALG::Matrix<3, 3> outer_product;
  outer_product.MultiplyNT(vec_a, vec_b);

  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = i; j < 3; ++j)
      ab_ba(IndexMappings::SymToVoigt6(i, j)) += outer_product(i, j) + outer_product(j, i);

  // scale off-diagonal values
  ScaleOffDiagonalVals(ab_ba);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::MultiplyTensorVector(
    const LINALG::Matrix<6, 1>& strain, const LINALG::Matrix<3, 1>& vec, LINALG::Matrix<3, 1>& res)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
    {
      const double fac = UnscaleFactor(IndexMappings::SymToVoigt6(i, j));
      res(i, 0) += strain(IndexMappings::SymToVoigt6(i, j)) * fac * vec(j, 0);
    }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::PowerOfSymmetricTensor(
    const unsigned pow, const LINALG::Matrix<6, 1>& strain, LINALG::Matrix<6, 1>& strain_pow)
{
  std::copy(strain.A(), strain.A() + 6, strain_pow.A());

  if (pow > 1)
  {
    // unscale the off-diagonal values
    UnscaleOffDiagonalVals(strain_pow);

    LINALG::Matrix<6, 1> prod(false);

    using vmap = IndexMappings;
    for (unsigned p = 1; p < pow; ++p)
    {
      std::fill(prod.A(), prod.A() + 6, 0.0);

      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = i; j < 3; ++j)
          for (unsigned k = 0; k < 3; ++k)
            prod(vmap::SymToVoigt6(i, j), 0) += strain_pow(vmap::SymToVoigt6(i, k), 0) *
                                                unscale_fac_[vmap::SymToVoigt6(k, j)] *
                                                strain(vmap::SymToVoigt6(k, j), 0);

      std::copy(prod.A(), prod.A() + 6, strain_pow.A());
    }

    // scale the off-diagonal values again
    ScaleOffDiagonalVals(strain_pow);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::InverseTensor(
    const LINALG::Matrix<6, 1>& tens, LINALG::Matrix<6, 1>& tens_inv)
{
  double det = Determinant(tens);
  tens_inv(0) = (tens(1) * tens(2) - UnscaleFactor(4) * UnscaleFactor(4) * tens(4) * tens(4)) /
                det * ScaleFactor(0);
  tens_inv(1) = (tens(0) * tens(2) - UnscaleFactor(5) * UnscaleFactor(5) * tens(5) * tens(5)) /
                det * ScaleFactor(1);
  tens_inv(2) = (tens(0) * tens(1) - UnscaleFactor(3) * UnscaleFactor(3) * tens(3) * tens(3)) /
                det * ScaleFactor(2);
  tens_inv(3) = (UnscaleFactor(5) * UnscaleFactor(4) * tens(5) * tens(4) -
                    UnscaleFactor(3) * UnscaleFactor(2) * tens(3) * tens(2)) /
                det * ScaleFactor(3);
  tens_inv(4) = (UnscaleFactor(3) * UnscaleFactor(5) * tens(3) * tens(5) -
                    UnscaleFactor(0) * UnscaleFactor(4) * tens(0) * tens(4)) /
                det * ScaleFactor(4);
  tens_inv(5) = (UnscaleFactor(3) * UnscaleFactor(4) * tens(3) * tens(4) -
                    UnscaleFactor(5) * UnscaleFactor(1) * tens(5) * tens(1)) /
                det * ScaleFactor(5);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::ToStressLike(
    const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<6, 1>& vtensor_out)
{
  for (unsigned i = 0; i < 6; ++i) vtensor_out(i) = UnscaleFactor(i) * vtensor_in(i);
}

template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::ToStrainLike(
    const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<6, 1>& vtensor_out)
{
  for (unsigned i = 0; i < 6; ++i)
    vtensor_out(i) = UnscaleFactor(i) * vtensor_in(i) * VStrainUtils::ScaleFactor(i);
}

template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::ToMatrix(
    const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<3, 3>& tensor_out)
{
  for (int i = 0; i < 3; ++i) tensor_out(i, i) = vtensor_in(i);
  tensor_out(0, 1) = tensor_out(1, 0) = UnscaleFactor(3) * vtensor_in(3);
  tensor_out(1, 2) = tensor_out(2, 1) = UnscaleFactor(4) * vtensor_in(4);
  tensor_out(0, 2) = tensor_out(2, 0) = UnscaleFactor(5) * vtensor_in(5);
}

template <NotationType type>
template <typename T>
void UTILS::VOIGT::VoigtUtils<type>::MatrixToVector(
    const LINALG::Matrix<3, 3, T>& tensor_in, LINALG::Matrix<6, 1, T>& vtensor_out)
{
  for (int i = 0; i < 3; ++i) vtensor_out(i) = tensor_in(i, i);
  vtensor_out(3) = 0.5 * ScaleFactor(3) * (tensor_in(0, 1) + tensor_in(1, 0));
  vtensor_out(4) = 0.5 * ScaleFactor(4) * (tensor_in(1, 2) + tensor_in(2, 1));
  vtensor_out(5) = 0.5 * ScaleFactor(5) * (tensor_in(0, 2) + tensor_in(2, 0));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::ScaleOffDiagonalVals(LINALG::Matrix<6, 1>& strain)
{
  for (unsigned i = 3; i < 6; ++i) strain(i, 0) *= ScaleFactor(i);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <NotationType type>
void UTILS::VOIGT::VoigtUtils<type>::UnscaleOffDiagonalVals(LINALG::Matrix<6, 1>& strain)
{
  for (unsigned i = 3; i < 6; ++i) strain(i, 0) *= UnscaleFactor(i);
}

template <>
const double UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::strain>::unscale_fac_[6] = {
    1.0, 1.0, 1.0, 0.5, 0.5, 0.5};
template <>
const double UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::stress>::unscale_fac_[6] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
template <>
const double UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::strain>::scale_fac_[6] = {
    1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
template <>
const double UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::stress>::scale_fac_[6] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// explicit template declarations
template class UTILS::VOIGT::VoigtUtils<NotationType::strain>;
template class UTILS::VOIGT::VoigtUtils<NotationType::stress>;

template void UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::strain>::MatrixToVector<double>(
    LINALG::Matrix<3, 3, double> const& in, LINALG::Matrix<6, 1, double>& out);
template void UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::stress>::MatrixToVector<double>(
    LINALG::Matrix<3, 3, double> const& in, LINALG::Matrix<6, 1, double>& out);

using FAD = Sacado::Fad::DFad<double>;
template void UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::strain>::MatrixToVector<FAD>(
    LINALG::Matrix<3, 3, FAD> const& in, LINALG::Matrix<6, 1, FAD>& out);
template void UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::stress>::MatrixToVector<FAD>(
    LINALG::Matrix<3, 3, FAD> const& in, LINALG::Matrix<6, 1, FAD>& out);

template void UTILS::VOIGT::FourthOrderIdentityMatrix<NotationType::stress, NotationType::stress>(
    LINALG::Matrix<6, 6>& id);
template void UTILS::VOIGT::FourthOrderIdentityMatrix<NotationType::stress, NotationType::strain>(
    LINALG::Matrix<6, 6>& id);
