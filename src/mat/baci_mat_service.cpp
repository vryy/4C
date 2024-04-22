/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_mat_service.hpp"

#include "baci_linalg_four_tensor.hpp"
#include "baci_linalg_utils_densematrix_eigen.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mixture_constituent_remodelfiber_material_exponential.hpp"
#include "baci_mixture_constituent_remodelfiber_material_exponential_active.hpp"
#include "baci_mixture_growth_strategy_anisotropic.hpp"
#include "baci_mixture_growth_strategy_isotropic.hpp"
#include "baci_mixture_growth_strategy_stiffness.hpp"
#include "baci_mixture_prestress_strategy_constant.hpp"
#include "baci_mixture_prestress_strategy_isocyl.hpp"
#include "baci_mixture_prestress_strategy_iterative.hpp"
#include "baci_mixture_rule_function.hpp"
#include "baci_mixture_rule_growthremodel.hpp"
#include "baci_mixture_rule_map.hpp"
#include "baci_mixture_rule_simple.hpp"

#include <Sacado.hpp>

using FAD = Sacado::Fad::DFad<double>;

FOUR_C_NAMESPACE_OPEN

template <typename T>
void MAT::AddtoCmatHolzapfelProduct(
    CORE::LINALG::Matrix<6, 6, T>& cmat, const CORE::LINALG::Matrix<6, 1, T>& invc, const T scalar)
{
#ifdef BACI_DEBUG
  if (cmat.numRows() != 6 or cmat.numCols() != 6 or invc.numRows() != 6)
    dserror("Wrong dimensions in function AddtoCmatHolzapfelProduct");
#endif

  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0, 0) += scalar * invc(0) * invc(0);
  cmat(0, 1) += scalar * invc(3) * invc(3);
  cmat(0, 2) += scalar * invc(5) * invc(5);
  cmat(0, 3) += scalar * invc(0) * invc(3);
  cmat(0, 4) += scalar * invc(3) * invc(5);
  cmat(0, 5) += scalar * invc(0) * invc(5);

  cmat(1, 0) += scalar * invc(3) * invc(3);
  cmat(1, 1) += scalar * invc(1) * invc(1);
  cmat(1, 2) += scalar * invc(4) * invc(4);
  cmat(1, 3) += scalar * invc(3) * invc(1);
  cmat(1, 4) += scalar * invc(1) * invc(4);
  cmat(1, 5) += scalar * invc(3) * invc(4);

  cmat(2, 0) += scalar * invc(5) * invc(5);
  cmat(2, 1) += scalar * invc(4) * invc(4);
  cmat(2, 2) += scalar * invc(2) * invc(2);
  cmat(2, 3) += scalar * invc(5) * invc(4);
  cmat(2, 4) += scalar * invc(4) * invc(2);
  cmat(2, 5) += scalar * invc(5) * invc(2);

  cmat(3, 0) += scalar * invc(0) * invc(3);
  cmat(3, 1) += scalar * invc(3) * invc(1);
  cmat(3, 2) += scalar * invc(5) * invc(4);
  cmat(3, 3) += scalar * 0.5 * (invc(0) * invc(1) + invc(3) * invc(3));
  cmat(3, 4) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(3, 5) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));

  cmat(4, 0) += scalar * invc(3) * invc(5);
  cmat(4, 1) += scalar * invc(1) * invc(4);
  cmat(4, 2) += scalar * invc(4) * invc(2);
  cmat(4, 3) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(4, 4) += scalar * 0.5 * (invc(1) * invc(2) + invc(4) * invc(4));
  cmat(4, 5) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));

  cmat(5, 0) += scalar * invc(0) * invc(5);
  cmat(5, 1) += scalar * invc(3) * invc(4);
  cmat(5, 2) += scalar * invc(5) * invc(2);
  cmat(5, 3) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));
  cmat(5, 4) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));
  cmat(5, 5) += scalar * 0.5 * (invc(0) * invc(2) + invc(5) * invc(5));
}

void MAT::ElastSymTensorMultiply(CORE::LINALG::Matrix<6, 6>& C, const double ScalarAB,
    const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B,
    const double ScalarThis)
{
  // everything in Voigt-Notation
  CORE::LINALG::Matrix<6, 1> AVoigt;
  CORE::LINALG::Matrix<6, 1> BVoigt;

  AVoigt(0, 0) = A(0, 0);
  AVoigt(1, 0) = A(1, 1);
  AVoigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3, 0) = A(1, 0);
  AVoigt(4, 0) = A(2, 1);
  AVoigt(5, 0) = A(2, 0);

  BVoigt(0, 0) = B(0, 0);
  BVoigt(1, 0) = B(1, 1);
  BVoigt(2, 0) = B(2, 2);
  BVoigt(3, 0) = B(1, 0);
  BVoigt(4, 0) = B(2, 1);
  BVoigt(5, 0) = B(2, 0);

  C.MultiplyNT(ScalarAB, AVoigt, BVoigt, ScalarThis);
}

void MAT::ElastSymTensorMultiplyAddSym(CORE::LINALG::Matrix<6, 6>& C, const double ScalarAB,
    const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B,
    const double ScalarThis)
{
#ifdef BACI_DEBUG
  // check sizes
  if (A.numRows() != A.numCols() || B.numRows() != B.numCols() || A.numRows() != 3 ||
      B.numRows() != 3)
  {
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.numRows() != C.numCols() || C.numRows() != 6) dserror("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  CORE::LINALG::Matrix<6, 1> AVoigt;
  CORE::LINALG::Matrix<6, 1> BVoigt;

  AVoigt(0, 0) = A(0, 0);
  AVoigt(1, 0) = A(1, 1);
  AVoigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  AVoigt(3, 0) = A(1, 0);
  AVoigt(4, 0) = A(2, 1);
  AVoigt(5, 0) = A(2, 0);

  BVoigt(0, 0) = B(0, 0);
  BVoigt(1, 0) = B(1, 1);
  BVoigt(2, 0) = B(2, 2);
  BVoigt(3, 0) = B(1, 0);
  BVoigt(4, 0) = B(2, 1);
  BVoigt(5, 0) = B(2, 0);

  C.MultiplyNT(ScalarAB, AVoigt, BVoigt, ScalarThis);
  C.MultiplyNT(ScalarAB, BVoigt, AVoigt, 1.0);
}

void MAT::ElastSymTensor_o_Multiply(CORE::LINALG::Matrix<6, 6>& C, const double ScalarAB,
    const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B,
    const double ScalarThis)
{
  const double ScalarABhalf = ScalarAB * 0.5;
  C(0, 0) = ScalarThis * C(0, 0) + ScalarAB * (A(0, 0) * B(0, 0));                          // C1111
  C(0, 1) = ScalarThis * C(0, 1) + ScalarAB * (A(0, 1) * B(0, 1));                          // C1122
  C(0, 2) = ScalarThis * C(0, 2) + ScalarAB * (A(0, 2) * B(0, 2));                          // C1133
  C(0, 3) = ScalarThis * C(0, 3) + ScalarABhalf * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));  // C1112
  C(0, 4) = ScalarThis * C(0, 4) + ScalarABhalf * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));  // C1123
  C(0, 5) = ScalarThis * C(0, 5) + ScalarABhalf * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));  // C1113

  C(1, 0) = ScalarThis * C(1, 0) + ScalarAB * (A(1, 0) * B(1, 0));                          // C2211
  C(1, 1) = ScalarThis * C(1, 1) + ScalarAB * (A(1, 1) * B(1, 1));                          // C2222
  C(1, 2) = ScalarThis * C(1, 2) + ScalarAB * (A(1, 2) * B(1, 2));                          // C2233
  C(1, 3) = ScalarThis * C(1, 3) + ScalarABhalf * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));  // C2212
  C(1, 4) = ScalarThis * C(1, 4) + ScalarABhalf * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));  // C2223
  C(1, 5) = ScalarThis * C(1, 5) + ScalarABhalf * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));  // C2213

  C(2, 0) = ScalarThis * C(2, 0) + ScalarAB * (A(2, 0) * B(2, 0));                          // C3311
  C(2, 1) = ScalarThis * C(2, 1) + ScalarAB * (A(2, 1) * B(2, 1));                          // C3322
  C(2, 2) = ScalarThis * C(2, 2) + ScalarAB * (A(2, 2) * B(2, 2));                          // C3333
  C(2, 3) = ScalarThis * C(2, 3) + ScalarABhalf * (A(2, 1) * B(2, 1) + A(2, 1) * B(2, 0));  // C3312
  C(2, 4) = ScalarThis * C(2, 4) + ScalarABhalf * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));  // C3323
  C(2, 5) = ScalarThis * C(2, 5) + ScalarABhalf * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));  // C3313

  C(3, 0) = ScalarThis * C(3, 0) + ScalarAB * (A(0, 0) * B(1, 0));                          // C1211
  C(3, 1) = ScalarThis * C(3, 1) + ScalarAB * (A(0, 1) * B(1, 1));                          // C1222
  C(3, 2) = ScalarThis * C(3, 2) + ScalarAB * (A(0, 2) * B(1, 2));                          // C1233
  C(3, 3) = ScalarThis * C(3, 3) + ScalarABhalf * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));  // C1212
  C(3, 4) = ScalarThis * C(3, 4) + ScalarABhalf * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));  // C1223
  C(3, 5) = ScalarThis * C(3, 5) + ScalarABhalf * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));  // C1213

  C(4, 0) = ScalarThis * C(4, 0) + ScalarAB * (A(1, 0) * B(2, 0));                          // C2311
  C(4, 1) = ScalarThis * C(4, 1) + ScalarAB * (A(1, 1) * B(2, 1));                          // C2322
  C(4, 2) = ScalarThis * C(4, 2) + ScalarAB * (A(1, 2) * B(2, 2));                          // C2333
  C(4, 3) = ScalarThis * C(4, 3) + ScalarABhalf * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));  // C2312
  C(4, 4) = ScalarThis * C(4, 4) + ScalarABhalf * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));  // C2323
  C(4, 5) = ScalarThis * C(4, 5) + ScalarABhalf * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));  // C2313

  C(5, 0) = ScalarThis * C(5, 0) + ScalarAB * (A(0, 0) * B(2, 0));                          // C1311
  C(5, 1) = ScalarThis * C(5, 1) + ScalarAB * (A(0, 1) * B(2, 1));                          // C1322
  C(5, 2) = ScalarThis * C(5, 2) + ScalarAB * (A(0, 2) * B(2, 2));                          // C1333
  C(5, 3) = ScalarThis * C(5, 3) + ScalarABhalf * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));  // C1312
  C(5, 4) = ScalarThis * C(5, 4) + ScalarABhalf * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));  // C1323
  C(5, 5) = ScalarThis * C(5, 5) + ScalarABhalf * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313
}

void MAT::VolumetrifyAndIsochorify(CORE::LINALG::Matrix<6, 1>* pk2vol,
    CORE::LINALG::Matrix<6, 6>* cvol, CORE::LINALG::Matrix<6, 1>* pk2iso,
    CORE::LINALG::Matrix<6, 6>* ciso, const CORE::LINALG::Matrix<6, 1>& gl,
    const CORE::LINALG::Matrix<6, 1>& pk2, const CORE::LINALG::Matrix<6, 6>& cmat)
{
  // useful call?
#ifdef BACI_DEBUG
  if ((pk2vol == nullptr) and (cvol == nullptr) and (pk2iso == nullptr) and (ciso == nullptr))
    dserror("Useful call? Apparently you do not want to compute anything");
#endif

  // right Cauchy--Green tensor
  // REMARK: stored in _strain_-like 6-Voigt vector
  CORE::LINALG::Matrix<6, 1> rcg(gl);
  rcg.Scale(2.0);
  for (int i = 0; i < 3; i++) rcg(i) += 1.0;

  // third invariant (determinant) of right Cauchy--Green strains
  const double rcg3rd = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                        0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                        0.25 * rcg(0) * rcg(4) * rcg(4);

  // inverse right Cauchy--Green tensor C^{-1}
  // REMARK: stored in as _stress_ 6-Voigt vector
  CORE::LINALG::Matrix<6, 1> icg(false);
  icg(0) = (rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4)) / rcg3rd;        // (C^{-1})^{11}
  icg(1) = (rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5)) / rcg3rd;        // (C^{-1})^{22}
  icg(2) = (rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3)) / rcg3rd;        // (C^{-1})^{33}
  icg(3) = (0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2)) / rcg3rd;  // (C^{-1})^{12}
  icg(4) = (0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4)) / rcg3rd;  // (C^{-1})^{23}
  icg(5) = (0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1)) / rcg3rd;  // (C^{-1})^{31}

  // double contraction of 2nd Piola--Kirchhoff stress and right Cauchy--Green strain,
  // i.e. in index notation S^{AB} C_{AB}
  // REMARK: equal to S^T C, because S is stress-like and C is strain-like 6-Voigt vector
  const double pk2rcg = pk2.Dot(rcg);

  // stress splitting
  {
    // volumetric 2nd Piola--Kirchhoff stress
    CORE::LINALG::Matrix<6, 1> pk2vol_(false);
    if (pk2vol != nullptr) pk2vol_.SetView(*pk2vol);
    pk2vol_.Update(pk2rcg / 3.0, icg);

    // isochoric 2nd Piola--Kirchhoff stress
    // S^{AB}_iso = S^{AB} - S^{AB}_{vol}
    if (pk2iso != nullptr) pk2iso->Update(1.0, pk2, -1.0, pk2vol_);
  }

  // elasticity tensor splitting
  {
    // 'linearised' 2nd Piola--Kirchhoff stress
    // S^{CD}_lin = S^{CD} + 1/2 C_{AB} C^{ABCD}
    CORE::LINALG::Matrix<6, 1> pk2lin(pk2);
    pk2lin.MultiplyTN(0.5, cmat, rcg, 1.0);  // transpose on purpose

    // volumetric part of constitutive tensor
    // C^{ABCD}_vol = 2/3 (C^{-1})^{AB} S^{CD}_lin
    //              - 2/3 (S^{EF} C_{EF}) ( 1/2 (
    //                (C^{-1})^{AC} (C^{-1})^{BD} + (C^{-1})^{AD} (C^{-1})^{BC}
    //              ) )
    CORE::LINALG::Matrix<6, 6> cvol_(false);
    if (cvol != nullptr) cvol_.SetView(*cvol);
    cvol_.MultiplyNT(2.0 / 3.0, icg, pk2lin);
    AddtoCmatHolzapfelProduct(cvol_, icg, -2.0 / 3.0 * pk2rcg);

    // isochoric part of constitutive tensor
    // C^{ABCD}_iso = C^{ABCD} - C^{ABCD}_vol
    if (ciso != nullptr) ciso->Update(1.0, cmat, -1.0, cvol_);
  }
}

void MAT::AddToCmatDerivTensorSquare(CORE::LINALG::Matrix<6, 6>& C, double ScalarDX2,
    CORE::LINALG::Matrix<3, 3> X, double ScalarThis)
{
  C(0, 0) = ScalarThis * C(0, 0) + ScalarDX2 * 2. * X(0, 0);  // C1111
  C(0, 1) = ScalarThis * C(0, 1);                             // C1122
  C(0, 2) = ScalarThis * C(0, 2);                             // C1133
  C(0, 3) = ScalarThis * C(0, 3) + ScalarDX2 * X(0, 1);       // C1112
  C(0, 4) = ScalarThis * C(0, 4);                             // C1123
  C(0, 5) = ScalarThis * C(0, 5) + ScalarDX2 * X(0, 2);       // C1113

  C(1, 0) = ScalarThis * C(1, 0);                             // C2211
  C(1, 1) = ScalarThis * C(1, 1) + ScalarDX2 * 2. * X(1, 1);  // C2222
  C(1, 2) = ScalarThis * C(1, 2);                             // C2233
  C(1, 3) = ScalarThis * C(1, 3) + ScalarDX2 * X(0, 1);       // C2212
  C(1, 4) = ScalarThis * C(1, 4) + ScalarDX2 * X(1, 2);       // C2223
  C(1, 5) = ScalarThis * C(1, 5);                             // C2213

  C(2, 0) = ScalarThis * C(2, 0);                             // C3311
  C(2, 1) = ScalarThis * C(2, 1);                             // C3322
  C(2, 2) = ScalarThis * C(2, 2) + ScalarDX2 * 2. * X(2, 2);  // C3333
  C(2, 3) = ScalarThis * C(2, 3);                             // C3312
  C(2, 4) = ScalarThis * C(2, 4) + ScalarDX2 * X(1, 2);       // C3323
  C(2, 5) = ScalarThis * C(2, 5) + ScalarDX2 * X(0, 2);       // C3313

  C(3, 0) = ScalarThis * C(3, 0) + ScalarDX2 * X(0, 1);                    // C1211
  C(3, 1) = ScalarThis * C(3, 1) + ScalarDX2 * X(0, 1);                    // C1222
  C(3, 2) = ScalarThis * C(3, 2);                                          // C1233
  C(3, 3) = ScalarThis * C(3, 3) + ScalarDX2 * 0.5 * (X(0, 0) + X(1, 1));  // C1212
  C(3, 4) = ScalarThis * C(3, 4) + ScalarDX2 * 0.5 * X(0, 2);              // C1223
  C(3, 5) = ScalarThis * C(3, 5) + ScalarDX2 * 0.5 * X(1, 2);              // C1213

  C(4, 0) = ScalarThis * C(4, 0);                                          // C2311
  C(4, 1) = ScalarThis * C(4, 1) + ScalarDX2 * X(1, 2);                    // C2322
  C(4, 2) = ScalarThis * C(4, 2) + ScalarDX2 * X(1, 2);                    // C2333
  C(4, 3) = ScalarThis * C(4, 3) + ScalarDX2 * 0.5 * X(0, 2);              // C2312
  C(4, 4) = ScalarThis * C(4, 4) + ScalarDX2 * 0.5 * (X(1, 1) + X(2, 2));  // C2323
  C(4, 5) = ScalarThis * C(4, 5) + ScalarDX2 * 0.5 * X(0, 1);              // C2313

  C(5, 0) = ScalarThis * C(5, 0) + ScalarDX2 * X(0, 2);                    // C1311
  C(5, 1) = ScalarThis * C(5, 1);                                          // C1322
  C(5, 2) = ScalarThis * C(5, 2) + ScalarDX2 * X(0, 2);                    // C1333
  C(5, 3) = ScalarThis * C(5, 3) + ScalarDX2 * 0.5 * X(1, 2);              // C1312
  C(5, 4) = ScalarThis * C(5, 4) + ScalarDX2 * 0.5 * X(0, 1);              // C1323
  C(5, 5) = ScalarThis * C(5, 5) + ScalarDX2 * 0.5 * (X(2, 2) + X(0, 0));  // C1313
}

void MAT::AddSymmetricHolzapfelProduct(CORE::LINALG::Matrix<6, 6>& X,
    const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B, const double fac)
{
  X(0, 0) += 4 * fac * A(0, 0) * B(0, 0);
  X(0, 3) += fac * (2 * A(0, 0) * B(1, 0) + 2 * A(1, 0) * B(0, 0));
  X(0, 5) += fac * (2 * A(0, 0) * B(2, 0) + 2 * A(2, 0) * B(0, 0));
  X(0, 1) += 4 * fac * A(1, 0) * B(1, 0);
  X(0, 4) += fac * (2 * A(1, 0) * B(2, 0) + 2 * A(2, 0) * B(1, 0));
  X(0, 2) += 4 * fac * A(2, 0) * B(2, 0);

  X(3, 0) += fac * (2 * A(0, 0) * B(0, 1) + 2 * A(0, 1) * B(0, 0));
  X(3, 3) += fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1) + A(1, 1) * B(0, 0) + A(0, 1) * B(1, 0));
  X(3, 5) += fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1) + A(2, 1) * B(0, 0) + A(0, 1) * B(2, 0));
  X(3, 1) += fac * (2 * A(1, 0) * B(1, 1) + 2 * A(1, 1) * B(1, 0));
  X(3, 4) += fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1) + A(2, 1) * B(1, 0) + A(1, 1) * B(2, 0));
  X(3, 2) += fac * (2 * A(2, 0) * B(2, 1) + 2 * A(2, 1) * B(2, 0));

  X(5, 0) += fac * (2 * A(0, 0) * B(0, 2) + 2 * A(0, 2) * B(0, 0));
  X(5, 3) += fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2) + A(1, 2) * B(0, 0) + A(0, 2) * B(1, 0));
  X(5, 5) += fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2) + A(2, 2) * B(0, 0) + A(0, 2) * B(2, 0));
  X(5, 1) += fac * (2 * A(1, 0) * B(1, 2) + 2 * A(1, 2) * B(1, 0));
  X(5, 4) += fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2) + A(2, 2) * B(1, 0) + A(1, 2) * B(2, 0));
  X(5, 2) += fac * (2 * A(2, 0) * B(2, 2) + 2 * A(2, 2) * B(2, 0));

  X(1, 0) += 4 * fac * A(0, 1) * B(0, 1);
  X(1, 3) += fac * (2 * A(0, 1) * B(1, 1) + 2 * A(1, 1) * B(0, 1));
  X(1, 5) += fac * (2 * A(0, 1) * B(2, 1) + 2 * A(2, 1) * B(0, 1));
  X(1, 1) += 4 * fac * A(1, 1) * B(1, 1);
  X(1, 4) += fac * (2 * A(1, 1) * B(2, 1) + 2 * A(2, 1) * B(1, 1));
  X(1, 2) += 4 * fac * A(2, 1) * B(2, 1);

  X(4, 0) += fac * (2 * A(0, 1) * B(0, 2) + 2 * A(0, 2) * B(0, 1));
  X(4, 3) += fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2) + A(1, 2) * B(0, 1) + A(0, 2) * B(1, 1));
  X(4, 5) += fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2) + A(2, 2) * B(0, 1) + A(0, 2) * B(2, 1));
  X(4, 1) += fac * (2 * A(1, 1) * B(1, 2) + 2 * A(1, 2) * B(1, 1));
  X(4, 4) += fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(1, 1) + A(1, 2) * B(2, 1));
  X(4, 2) += fac * (2 * A(2, 1) * B(2, 2) + 2 * A(2, 2) * B(2, 1));

  X(2, 0) += 4 * fac * A(0, 2) * B(0, 2);
  X(2, 3) += fac * (2 * A(0, 2) * B(1, 2) + 2 * A(1, 2) * B(0, 2));
  X(2, 5) += fac * (2 * A(0, 2) * B(2, 2) + 2 * A(2, 2) * B(0, 2));
  X(2, 1) += 4 * fac * A(1, 2) * B(1, 2);
  X(2, 4) += fac * (2 * A(1, 2) * B(2, 2) + 2 * A(2, 2) * B(1, 2));
  X(2, 2) += 4 * fac * A(2, 2) * B(2, 2);
}

template <typename T>
void MAT::AddRightNonSymmetricHolzapfelProduct(CORE::LINALG::Matrix<6, 9, T>& out,
    CORE::LINALG::Matrix<3, 3, T> const& A, CORE::LINALG::Matrix<3, 3, T> const& B, T const fac)
{
  out(0, 0) += 2 * fac * A(0, 0) * B(0, 0);
  out(0, 3) += 2 * fac * A(0, 0) * B(0, 1);
  out(0, 5) += 2 * fac * A(0, 0) * B(0, 2);
  out(0, 6) += 2 * fac * A(0, 1) * B(0, 0);
  out(0, 1) += 2 * fac * A(0, 1) * B(0, 1);
  out(0, 4) += 2 * fac * A(0, 1) * B(0, 2);
  out(0, 8) += 2 * fac * A(0, 2) * B(0, 0);
  out(0, 7) += 2 * fac * A(0, 2) * B(0, 1);
  out(0, 2) += 2 * fac * A(0, 2) * B(0, 2);

  out(1, 0) += 2 * fac * A(1, 0) * B(1, 0);
  out(1, 3) += 2 * fac * A(1, 0) * B(1, 1);
  out(1, 5) += 2 * fac * A(1, 0) * B(1, 2);
  out(1, 6) += 2 * fac * A(1, 1) * B(1, 0);
  out(1, 1) += 2 * fac * A(1, 1) * B(1, 1);
  out(1, 4) += 2 * fac * A(1, 1) * B(1, 2);
  out(1, 8) += 2 * fac * A(1, 2) * B(1, 0);
  out(1, 7) += 2 * fac * A(1, 2) * B(1, 1);
  out(1, 2) += 2 * fac * A(1, 2) * B(1, 2);

  out(2, 0) += 2 * fac * A(2, 0) * B(2, 0);
  out(2, 3) += 2 * fac * A(2, 0) * B(2, 1);
  out(2, 5) += 2 * fac * A(2, 0) * B(2, 2);
  out(2, 6) += 2 * fac * A(2, 1) * B(2, 0);
  out(2, 1) += 2 * fac * A(2, 1) * B(2, 1);
  out(2, 4) += 2 * fac * A(2, 1) * B(2, 2);
  out(2, 8) += 2 * fac * A(2, 2) * B(2, 0);
  out(2, 7) += 2 * fac * A(2, 2) * B(2, 1);
  out(2, 2) += 2 * fac * A(2, 2) * B(2, 2);

  out(3, 0) += fac * (A(0, 0) * B(1, 0) + A(1, 0) * B(0, 0));
  out(3, 3) += fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1));
  out(3, 5) += fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2));
  out(3, 6) += fac * (A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0));
  out(3, 1) += fac * (A(0, 1) * B(1, 1) + A(1, 1) * B(0, 1));
  out(3, 4) += fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2));
  out(3, 8) += fac * (A(0, 2) * B(1, 0) + A(1, 2) * B(0, 0));
  out(3, 7) += fac * (A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1));
  out(3, 2) += fac * (A(0, 2) * B(1, 2) + A(1, 2) * B(0, 2));

  out(4, 0) += fac * (A(1, 0) * B(2, 0) + A(2, 0) * B(1, 0));
  out(4, 3) += fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1));
  out(4, 5) += fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2));
  out(4, 6) += fac * (A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0));
  out(4, 1) += fac * (A(1, 1) * B(2, 1) + A(2, 1) * B(1, 1));
  out(4, 4) += fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2));
  out(4, 8) += fac * (A(1, 2) * B(2, 0) + A(2, 2) * B(1, 0));
  out(4, 7) += fac * (A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1));
  out(4, 2) += fac * (A(1, 2) * B(2, 2) + A(2, 2) * B(1, 2));

  out(5, 0) += fac * (A(0, 0) * B(2, 0) + A(2, 0) * B(0, 0));
  out(5, 3) += fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1));
  out(5, 5) += fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2));
  out(5, 6) += fac * (A(0, 1) * B(2, 0) + A(2, 1) * B(0, 0));
  out(5, 1) += fac * (A(0, 1) * B(2, 1) + A(2, 1) * B(0, 1));
  out(5, 4) += fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2));
  out(5, 8) += fac * (A(0, 2) * B(2, 0) + A(2, 2) * B(0, 0));
  out(5, 7) += fac * (A(0, 2) * B(2, 1) + A(2, 2) * B(0, 1));
  out(5, 2) += fac * (A(0, 2) * B(2, 2) + A(2, 2) * B(0, 2));
}

template <typename T>
void MAT::AddRightNonSymmetricHolzapfelProductStrainLike(CORE::LINALG::Matrix<6, 9, T>& out,
    CORE::LINALG::Matrix<3, 3, T> const& A, CORE::LINALG::Matrix<3, 3, T> const& B, T const fac)
{
  out(0, 0) += 2 * fac * A(0, 0) * B(0, 0);
  out(0, 3) += 2 * fac * A(0, 0) * B(0, 1);
  out(0, 5) += 2 * fac * A(0, 0) * B(0, 2);
  out(0, 6) += 2 * fac * A(0, 1) * B(0, 0);
  out(0, 1) += 2 * fac * A(0, 1) * B(0, 1);
  out(0, 4) += 2 * fac * A(0, 1) * B(0, 2);
  out(0, 8) += 2 * fac * A(0, 2) * B(0, 0);
  out(0, 7) += 2 * fac * A(0, 2) * B(0, 1);
  out(0, 2) += 2 * fac * A(0, 2) * B(0, 2);

  out(1, 0) += 2 * fac * A(1, 0) * B(1, 0);
  out(1, 3) += 2 * fac * A(1, 0) * B(1, 1);
  out(1, 5) += 2 * fac * A(1, 0) * B(1, 2);
  out(1, 6) += 2 * fac * A(1, 1) * B(1, 0);
  out(1, 1) += 2 * fac * A(1, 1) * B(1, 1);
  out(1, 4) += 2 * fac * A(1, 1) * B(1, 2);
  out(1, 8) += 2 * fac * A(1, 2) * B(1, 0);
  out(1, 7) += 2 * fac * A(1, 2) * B(1, 1);
  out(1, 2) += 2 * fac * A(1, 2) * B(1, 2);

  out(2, 0) += 2 * fac * A(2, 0) * B(2, 0);
  out(2, 3) += 2 * fac * A(2, 0) * B(2, 1);
  out(2, 5) += 2 * fac * A(2, 0) * B(2, 2);
  out(2, 6) += 2 * fac * A(2, 1) * B(2, 0);
  out(2, 1) += 2 * fac * A(2, 1) * B(2, 1);
  out(2, 4) += 2 * fac * A(2, 1) * B(2, 2);
  out(2, 8) += 2 * fac * A(2, 2) * B(2, 0);
  out(2, 7) += 2 * fac * A(2, 2) * B(2, 1);
  out(2, 2) += 2 * fac * A(2, 2) * B(2, 2);

  out(3, 0) += 2 * fac * (A(0, 0) * B(1, 0) + A(1, 0) * B(0, 0));
  out(3, 3) += 2 * fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1));
  out(3, 5) += 2 * fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2));
  out(3, 6) += 2 * fac * (A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0));
  out(3, 1) += 2 * fac * (A(0, 1) * B(1, 1) + A(1, 1) * B(0, 1));
  out(3, 4) += 2 * fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2));
  out(3, 8) += 2 * fac * (A(0, 2) * B(1, 0) + A(1, 2) * B(0, 0));
  out(3, 7) += 2 * fac * (A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1));
  out(3, 2) += 2 * fac * (A(0, 2) * B(1, 2) + A(1, 2) * B(0, 2));

  out(4, 0) += 2 * fac * (A(1, 0) * B(2, 0) + A(2, 0) * B(1, 0));
  out(4, 3) += 2 * fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1));
  out(4, 5) += 2 * fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2));
  out(4, 6) += 2 * fac * (A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0));
  out(4, 1) += 2 * fac * (A(1, 1) * B(2, 1) + A(2, 1) * B(1, 1));
  out(4, 4) += 2 * fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2));
  out(4, 8) += 2 * fac * (A(1, 2) * B(2, 0) + A(2, 2) * B(1, 0));
  out(4, 7) += 2 * fac * (A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1));
  out(4, 2) += 2 * fac * (A(1, 2) * B(2, 2) + A(2, 2) * B(1, 2));

  out(5, 0) += 2 * fac * (A(0, 0) * B(2, 0) + A(2, 0) * B(0, 0));
  out(5, 3) += 2 * fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1));
  out(5, 5) += 2 * fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2));
  out(5, 6) += 2 * fac * (A(0, 1) * B(2, 0) + A(2, 1) * B(0, 0));
  out(5, 1) += 2 * fac * (A(0, 1) * B(2, 1) + A(2, 1) * B(0, 1));
  out(5, 4) += 2 * fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2));
  out(5, 8) += 2 * fac * (A(0, 2) * B(2, 0) + A(2, 2) * B(0, 0));
  out(5, 7) += 2 * fac * (A(0, 2) * B(2, 1) + A(2, 2) * B(0, 1));
  out(5, 2) += 2 * fac * (A(0, 2) * B(2, 2) + A(2, 2) * B(0, 2));
}

void MAT::AddLeftNonSymmetricHolzapfelProduct(CORE::LINALG::Matrix<9, 6>& out,
    CORE::LINALG::Matrix<3, 3> const& A, CORE::LINALG::Matrix<3, 3> const& B, double const fac)
{
  out(0, 0) += fac * 2 * A(0, 0) * B(0, 0);
  out(0, 3) += fac * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));
  out(0, 5) += fac * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));
  out(0, 1) += fac * 2 * A(0, 1) * B(0, 1);
  out(0, 4) += fac * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));
  out(0, 2) += fac * 2 * A(0, 2) * B(0, 2);

  out(3, 0) += fac * 2 * A(0, 0) * B(1, 0);
  out(3, 3) += fac * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));
  out(3, 5) += fac * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));
  out(3, 1) += fac * 2 * A(0, 1) * B(1, 1);
  out(3, 4) += fac * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));
  out(3, 2) += fac * 2 * A(0, 2) * B(1, 2);

  out(5, 0) += fac * 2 * A(0, 0) * B(2, 0);
  out(5, 3) += fac * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));
  out(5, 5) += fac * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));
  out(5, 1) += fac * 2 * A(0, 1) * B(2, 1);
  out(5, 4) += fac * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));
  out(5, 2) += fac * 2 * A(0, 2) * B(2, 2);

  out(6, 0) += fac * 2 * A(1, 0) * B(0, 0);
  out(6, 3) += fac * (A(1, 0) * B(0, 1) + A(1, 1) * B(0, 0));
  out(6, 5) += fac * (A(1, 0) * B(0, 2) + A(1, 2) * B(0, 0));
  out(6, 1) += fac * 2 * A(1, 1) * B(0, 1);
  out(6, 4) += fac * (A(1, 1) * B(0, 2) + A(1, 2) * B(0, 1));
  out(6, 2) += fac * 2 * A(1, 2) * B(0, 2);

  out(1, 0) += fac * 2 * A(1, 0) * B(1, 0);
  out(1, 3) += fac * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));
  out(1, 5) += fac * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));
  out(1, 1) += fac * 2 * A(1, 1) * B(1, 1);
  out(1, 4) += fac * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));
  out(1, 2) += fac * 2 * A(1, 2) * B(1, 2);

  out(4, 0) += fac * 2 * A(1, 0) * B(2, 0);
  out(4, 3) += fac * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));
  out(4, 5) += fac * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));
  out(4, 1) += fac * 2 * A(1, 1) * B(2, 1);
  out(4, 4) += fac * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));
  out(4, 2) += fac * 2 * A(1, 2) * B(2, 2);

  out(8, 0) += fac * 2 * A(2, 0) * B(0, 0);
  out(8, 3) += fac * (A(2, 0) * B(0, 1) + A(2, 1) * B(0, 0));
  out(8, 5) += fac * (A(2, 0) * B(0, 2) + A(2, 2) * B(0, 0));
  out(8, 1) += fac * 2 * A(2, 1) * B(0, 1);
  out(8, 4) += fac * (A(2, 1) * B(0, 2) + A(2, 2) * B(0, 1));
  out(8, 2) += fac * 2 * A(2, 2) * B(0, 2);

  out(7, 0) += fac * 2 * A(2, 0) * B(1, 0);
  out(7, 3) += fac * (A(2, 0) * B(1, 1) + A(2, 1) * B(1, 0));
  out(7, 5) += fac * (A(2, 0) * B(1, 2) + A(2, 2) * B(1, 0));
  out(7, 1) += fac * 2 * A(2, 1) * B(1, 1);
  out(7, 4) += fac * (A(2, 1) * B(1, 2) + A(2, 2) * B(1, 1));
  out(7, 2) += fac * 2 * A(2, 2) * B(1, 2);

  out(2, 0) += fac * 2 * A(2, 0) * B(2, 0);
  out(2, 3) += fac * (A(2, 0) * B(2, 1) + A(2, 1) * B(2, 0));
  out(2, 5) += fac * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));
  out(2, 1) += fac * 2 * A(2, 1) * B(2, 1);
  out(2, 4) += fac * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));
  out(2, 2) += fac * 2 * A(2, 2) * B(2, 2);
}

void MAT::AddNonSymmetricProduct(double const& fac, CORE::LINALG::Matrix<3, 3> const& A,
    CORE::LINALG::Matrix<3, 3> const& B, CORE::LINALG::Matrix<9, 9>& out)
{
  out(0, 0) += fac * A(0, 0) * B(0, 0);
  out(0, 3) += fac * A(0, 0) * B(1, 0);
  out(0, 5) += fac * A(0, 0) * B(2, 0);
  out(0, 6) += fac * A(0, 1) * B(0, 0);
  out(0, 1) += fac * A(0, 1) * B(1, 0);
  out(0, 4) += fac * A(0, 1) * B(2, 0);
  out(0, 8) += fac * A(0, 2) * B(0, 0);
  out(0, 7) += fac * A(0, 2) * B(1, 0);
  out(0, 2) += fac * A(0, 2) * B(2, 0);

  out(3, 0) += fac * A(0, 0) * B(0, 1);
  out(3, 3) += fac * A(0, 0) * B(1, 1);
  out(3, 5) += fac * A(0, 0) * B(2, 1);
  out(3, 6) += fac * A(0, 1) * B(0, 1);
  out(3, 1) += fac * A(0, 1) * B(1, 1);
  out(3, 4) += fac * A(0, 1) * B(2, 1);
  out(3, 8) += fac * A(0, 2) * B(0, 1);
  out(3, 7) += fac * A(0, 2) * B(1, 1);
  out(3, 2) += fac * A(0, 2) * B(2, 1);

  out(5, 0) += fac * A(0, 0) * B(0, 2);
  out(5, 3) += fac * A(0, 0) * B(1, 2);
  out(5, 5) += fac * A(0, 0) * B(2, 2);
  out(5, 6) += fac * A(0, 1) * B(0, 2);
  out(5, 1) += fac * A(0, 1) * B(1, 2);
  out(5, 4) += fac * A(0, 1) * B(2, 2);
  out(5, 8) += fac * A(0, 2) * B(0, 2);
  out(5, 7) += fac * A(0, 2) * B(1, 2);
  out(5, 2) += fac * A(0, 2) * B(2, 2);

  out(6, 0) += fac * A(1, 0) * B(0, 0);
  out(6, 3) += fac * A(1, 0) * B(1, 0);
  out(6, 5) += fac * A(1, 0) * B(2, 0);
  out(6, 6) += fac * A(1, 1) * B(0, 0);
  out(6, 1) += fac * A(1, 1) * B(1, 0);
  out(6, 4) += fac * A(1, 1) * B(2, 0);
  out(6, 8) += fac * A(1, 2) * B(0, 0);
  out(6, 7) += fac * A(1, 2) * B(1, 0);
  out(6, 2) += fac * A(1, 2) * B(2, 0);

  out(1, 0) += fac * A(1, 0) * B(0, 1);
  out(1, 3) += fac * A(1, 0) * B(1, 1);
  out(1, 5) += fac * A(1, 0) * B(2, 1);
  out(1, 6) += fac * A(1, 1) * B(0, 1);
  out(1, 1) += fac * A(1, 1) * B(1, 1);
  out(1, 4) += fac * A(1, 1) * B(2, 1);
  out(1, 8) += fac * A(1, 2) * B(0, 1);
  out(1, 7) += fac * A(1, 2) * B(1, 1);
  out(1, 2) += fac * A(1, 2) * B(2, 1);

  out(4, 0) += fac * A(1, 0) * B(0, 2);
  out(4, 3) += fac * A(1, 0) * B(1, 2);
  out(4, 5) += fac * A(1, 0) * B(2, 2);
  out(4, 6) += fac * A(1, 1) * B(0, 2);
  out(4, 1) += fac * A(1, 1) * B(1, 2);
  out(4, 4) += fac * A(1, 1) * B(2, 2);
  out(4, 8) += fac * A(1, 2) * B(0, 2);
  out(4, 7) += fac * A(1, 2) * B(1, 2);
  out(4, 2) += fac * A(1, 2) * B(2, 2);

  out(8, 0) += fac * A(2, 0) * B(0, 0);
  out(8, 3) += fac * A(2, 0) * B(1, 0);
  out(8, 5) += fac * A(2, 0) * B(2, 0);
  out(8, 6) += fac * A(2, 1) * B(0, 0);
  out(8, 1) += fac * A(2, 1) * B(1, 0);
  out(8, 4) += fac * A(2, 1) * B(2, 0);
  out(8, 8) += fac * A(2, 2) * B(0, 0);
  out(8, 7) += fac * A(2, 2) * B(1, 0);
  out(8, 2) += fac * A(2, 2) * B(2, 0);

  out(7, 0) += fac * A(2, 0) * B(0, 1);
  out(7, 3) += fac * A(2, 0) * B(1, 1);
  out(7, 5) += fac * A(2, 0) * B(2, 1);
  out(7, 6) += fac * A(2, 1) * B(0, 1);
  out(7, 1) += fac * A(2, 1) * B(1, 1);
  out(7, 4) += fac * A(2, 1) * B(2, 1);
  out(7, 8) += fac * A(2, 2) * B(0, 1);
  out(7, 7) += fac * A(2, 2) * B(1, 1);
  out(7, 2) += fac * A(2, 2) * B(2, 1);

  out(2, 0) += fac * A(2, 0) * B(0, 2);
  out(2, 3) += fac * A(2, 0) * B(1, 2);
  out(2, 5) += fac * A(2, 0) * B(2, 2);
  out(2, 6) += fac * A(2, 1) * B(0, 2);
  out(2, 1) += fac * A(2, 1) * B(1, 2);
  out(2, 4) += fac * A(2, 1) * B(2, 2);
  out(2, 8) += fac * A(2, 2) * B(0, 2);
  out(2, 7) += fac * A(2, 2) * B(1, 2);
  out(2, 2) += fac * A(2, 2) * B(2, 2);
}

void MAT::AddDerivInvABInvBProduct(double const& fac, const CORE::LINALG::Matrix<6, 1>& invA,
    const CORE::LINALG::Matrix<6, 1>& invABinvA, CORE::LINALG::Matrix<6, 6>& out)
{
  out(0, 0) -= 2.0 * fac * (invA(0) * invABinvA(0));
  out(0, 1) -= 2.0 * fac * (invA(3) * invABinvA(3));
  out(0, 2) -= 2.0 * fac * (invA(5) * invABinvA(5));
  out(0, 3) -= fac * (invA(0) * invABinvA(3) + invA(3) * invABinvA(0));
  out(0, 4) -= fac * (invA(3) * invABinvA(5) + invA(5) * invABinvA(3));
  out(0, 5) -= fac * (invA(0) * invABinvA(5) + invA(5) * invABinvA(0));

  out(1, 0) -= 2.0 * fac * (invA(3) * invABinvA(3));
  out(1, 1) -= 2.0 * fac * (invA(1) * invABinvA(1));
  out(1, 2) -= 2.0 * fac * (invA(4) * invABinvA(4));
  out(1, 3) -= fac * (invA(3) * invABinvA(1) + invA(1) * invABinvA(3));
  out(1, 4) -= fac * (invA(1) * invABinvA(4) + invA(4) * invABinvA(1));
  out(1, 5) -= fac * (invA(3) * invABinvA(4) + invA(4) * invABinvA(3));

  out(2, 0) -= 2.0 * fac * (invA(5) * invABinvA(5));
  out(2, 1) -= 2.0 * fac * (invA(4) * invABinvA(4));
  out(2, 2) -= 2.0 * fac * (invA(2) * invABinvA(2));
  out(2, 3) -= fac * (invA(5) * invABinvA(4) + invA(4) * invABinvA(5));
  out(2, 4) -= fac * (invA(4) * invABinvA(2) + invA(2) * invABinvA(4));
  out(2, 5) -= fac * (invA(5) * invABinvA(2) + invA(2) * invABinvA(5));

  out(3, 0) -= fac * (invA(0) * invABinvA(3) + invA(3) * invABinvA(0));
  out(3, 1) -= fac * (invA(3) * invABinvA(1) + invA(1) * invABinvA(3));
  out(3, 2) -= fac * (invA(5) * invABinvA(4) + invA(4) * invABinvA(5));
  out(3, 3) -=
      fac * (0.5 * (invA(0) * invABinvA(1) + invA(1) * invABinvA(0)) + invA(3) * invABinvA(3));
  out(3, 4) -= 0.5 * fac *
               (invA(3) * invABinvA(4) + invA(5) * invABinvA(1) + invA(1) * invABinvA(5) +
                   invA(4) * invABinvA(3));
  out(3, 5) -= 0.5 * fac *
               (invA(0) * invABinvA(4) + invA(5) * invABinvA(3) + invA(3) * invABinvA(5) +
                   invA(4) * invABinvA(0));

  out(4, 0) -= fac * (invA(3) * invABinvA(5) + invA(5) * invABinvA(3));
  out(4, 1) -= fac * (invA(1) * invABinvA(4) + invA(4) * invABinvA(1));
  out(4, 2) -= fac * (invA(4) * invABinvA(2) + invA(2) * invABinvA(4));
  out(4, 3) -= 0.5 * fac *
               (invA(3) * invABinvA(4) + invA(1) * invABinvA(5) + invA(5) * invABinvA(1) +
                   invA(4) * invABinvA(3));
  out(4, 4) -=
      fac * (0.5 * (invA(1) * invABinvA(2) + invA(2) * invABinvA(1)) + invA(4) * invABinvA(4));
  out(4, 5) -= 0.5 * fac *
               (invA(3) * invABinvA(2) + invA(4) * invABinvA(5) + invA(5) * invABinvA(4) +
                   invA(2) * invABinvA(3));

  out(5, 0) -= fac * (invA(0) * invABinvA(5) + invA(5) * invABinvA(0));
  out(5, 1) -= fac * (invA(3) * invABinvA(4) + invA(4) * invABinvA(3));
  out(5, 2) -= fac * (invA(5) * invABinvA(2) + invA(2) * invABinvA(5));
  out(5, 3) -= 0.5 * fac *
               (invA(0) * invABinvA(4) + invA(3) * invABinvA(5) + invA(5) * invABinvA(3) +
                   invA(4) * invABinvA(0));
  out(5, 4) -= 0.5 * fac *
               (invA(3) * invABinvA(2) + invA(5) * invABinvA(4) + invA(4) * invABinvA(5) +
                   invA(2) * invABinvA(3));
  out(5, 5) -=
      fac * (0.5 * (invA(0) * invABinvA(2) + invA(2) * invABinvA(0)) + invA(5) * invABinvA(5));
}

void MAT::InvariantsModified(
    CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<3, 1>& prinv)
{
  // 1st invariant, trace
  modinv(0) = prinv(0) * std::pow(prinv(2), -1. / 3.);
  // 2nd invariant
  modinv(1) = prinv(1) * std::pow(prinv(2), -2. / 3.);
  // J
  modinv(2) = std::pow(prinv(2), 1. / 2.);
}

void MAT::StretchesPrincipal(CORE::LINALG::Matrix<3, 1>& prstr, CORE::LINALG::Matrix<3, 3>& prdir,
    const CORE::LINALG::Matrix<6, 1>& rcg)
{
  // create right Cauchy-Green 2-tensor
  CORE::LINALG::Matrix<3, 3> rcgt(false);
  rcgt(0, 0) = rcg(0);
  rcgt(1, 1) = rcg(1);
  rcgt(2, 2) = rcg(2);
  rcgt(0, 1) = rcgt(1, 0) = 0.5 * rcg(3);
  rcgt(1, 2) = rcgt(2, 1) = 0.5 * rcg(4);
  rcgt(2, 0) = rcgt(0, 2) = 0.5 * rcg(5);

  // eigenvalue decomposition
  CORE::LINALG::Matrix<3, 3> prstr2;  // squared principal stretches
  CORE::LINALG::SYEV(rcgt, prstr2, prdir);

  // THE principal stretches
  for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));
}

void MAT::StretchesModified(
    CORE::LINALG::Matrix<3, 1>& modstr, const CORE::LINALG::Matrix<3, 1>& prstr)
{
  // determinant of deformation gradient
  const double detdefgrad = prstr(0) * prstr(1) * prstr(2);

  // determine modified principal stretches
  modstr.Update(std::pow(detdefgrad, -1.0 / 3.0), prstr);
}

template <class T>
T* MAT::CreateMaterialParameterInstance(Teuchos::RCP<MAT::PAR::Material> curmat)
{
  if (curmat->Parameter() == nullptr)
  {
    curmat->SetParameter(new T(curmat));
  }
  auto* params = dynamic_cast<T*>(curmat->Parameter());
  return params;
}

template <int dim>
void MAT::ClearFourTensor(CORE::LINALG::FourTensor<dim>& fourTensor)
{
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          fourTensor(i, j, k, l) = 0.0;
        }
      }
    }
  }
}

template <int dim>
void MAT::MultiplyFourTensorMatrix(CORE::LINALG::FourTensor<dim>& fourTensorResult,
    const CORE::LINALG::FourTensor<dim>& fourTensor, const CORE::LINALG::Matrix<dim, dim>& matrix,
    const bool clearResultTensor)
{
  if (clearResultTensor) ClearFourTensor(fourTensorResult);
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          for (int m = 0; m < dim; ++m)
          {  // C^ijkl = A^ijkm * B_m^l
            fourTensorResult(i, j, k, l) += fourTensor(i, j, k, m) * matrix(m, l);
          }
        }
      }
    }
  }
}

template <int dim>
void MAT::MultiplyMatrixFourTensor(CORE::LINALG::FourTensor<dim>& fourTensorResult,
    const CORE::LINALG::Matrix<dim, dim>& matrix, const CORE::LINALG::FourTensor<dim>& fourTensor,
    const bool clearResultTensor)
{
  if (clearResultTensor) ClearFourTensor(fourTensorResult);
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          for (int m = 0; m < dim; ++m)
          {
            // C^ijkl = B^i_m * A^mjkl
            fourTensorResult(i, j, k, l) += matrix(i, m) * fourTensor(m, j, k, l);
          }
        }
      }
    }
  }
}

template <int dim>
void MAT::MultiplyMatrixFourTensorBySecondIndex(CORE::LINALG::FourTensor<dim>& fourTensorResult,
    const CORE::LINALG::Matrix<dim, dim>& matrix, const CORE::LINALG::FourTensor<dim>& fourTensor,
    const bool clearResultTensor)
{
  if (clearResultTensor) ClearFourTensor(fourTensorResult);
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          for (int m = 0; m < dim; ++m)
          {
            // C^ijkl = B_m^j * A^imkl
            fourTensorResult(i, j, k, l) += matrix(m, j) * fourTensor(i, m, k, l);
          }
        }
      }
    }
  }
}

template <int dim>
void MAT::MultiplyFourTensorFourTensor(CORE::LINALG::FourTensor<dim>& fourTensorResult,
    const CORE::LINALG::FourTensor<dim>& fourTensor1,
    const CORE::LINALG::FourTensor<dim>& fourTensor2, const bool clearResultTensor)
{
  if (clearResultTensor) ClearFourTensor(fourTensorResult);
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          // C^ijkl = A^ij_ab * B^abkl
          for (int a = 0; a < dim; ++a)
            for (int b = 0; b < dim; ++b)
              fourTensorResult(i, j, k, l) += fourTensor1(i, j, a, b) * fourTensor2(a, b, k, l);
        }
      }
    }
  }
}

template <int dim>
CORE::LINALG::Matrix<6, 6> MAT::PullBackFourTensor(
    const CORE::LINALG::Matrix<dim, dim>& defgr, const CORE::LINALG::Matrix<6, 6>& cMatVoigt)
{
#ifdef BACI_DEBUG
  if (dim != 3) dserror("Current implementation only valid for dim = 3.");
#endif

  CORE::LINALG::FourTensor<dim> cMatTensor(true);
  SetupFourTensor(cMatTensor, cMatVoigt);

  // We can use the fact that cMatResultVoigt(i,j,k,l)=cMatResultVoigt(k,l,i,j) if we have a
  // hyper-elastic material
  CORE::LINALG::Matrix<6, 6> cMatResultVoigt(true);

  cMatResultVoigt(0, 0) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 0, 0, 0);
  cMatResultVoigt(0, 1) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 0, 1, 1);
  cMatResultVoigt(0, 2) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 0, 2, 2);
  cMatResultVoigt(0, 3) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 0, 0, 1);
  cMatResultVoigt(0, 4) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 0, 1, 2);
  cMatResultVoigt(0, 5) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 0, 0, 2);
  cMatResultVoigt(1, 0) = cMatResultVoigt(0, 1);
  cMatResultVoigt(1, 1) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 1, 1, 1);
  cMatResultVoigt(1, 2) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 1, 2, 2);
  cMatResultVoigt(1, 3) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 1, 0, 1);
  cMatResultVoigt(1, 4) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 1, 1, 2);
  cMatResultVoigt(1, 5) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 1, 0, 2);
  cMatResultVoigt(2, 0) = cMatResultVoigt(0, 2);
  cMatResultVoigt(2, 1) = cMatResultVoigt(1, 2);
  cMatResultVoigt(2, 2) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 2, 2, 2, 2);
  cMatResultVoigt(2, 3) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 2, 2, 0, 1);
  cMatResultVoigt(2, 4) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 2, 2, 1, 2);
  cMatResultVoigt(2, 5) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 2, 2, 0, 2);
  cMatResultVoigt(3, 0) = cMatResultVoigt(0, 3);
  cMatResultVoigt(3, 1) = cMatResultVoigt(1, 3);
  cMatResultVoigt(3, 2) = cMatResultVoigt(2, 3);
  cMatResultVoigt(3, 3) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 1, 0, 1);
  cMatResultVoigt(3, 4) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 1, 1, 2);
  cMatResultVoigt(3, 5) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 1, 0, 2);
  cMatResultVoigt(4, 0) = cMatResultVoigt(0, 4);
  cMatResultVoigt(4, 1) = cMatResultVoigt(1, 4);
  cMatResultVoigt(4, 2) = cMatResultVoigt(2, 4);
  cMatResultVoigt(4, 3) = cMatResultVoigt(3, 4);
  cMatResultVoigt(4, 4) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 2, 1, 2);
  cMatResultVoigt(4, 5) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 1, 2, 0, 2);
  cMatResultVoigt(5, 0) = cMatResultVoigt(0, 5);
  cMatResultVoigt(5, 1) = cMatResultVoigt(1, 5);
  cMatResultVoigt(5, 2) = cMatResultVoigt(2, 5);
  cMatResultVoigt(5, 3) = cMatResultVoigt(3, 5);
  cMatResultVoigt(5, 4) = cMatResultVoigt(4, 5);
  cMatResultVoigt(5, 5) = PullBackFourTensorijkl<dim>(defgr, cMatTensor, 0, 2, 0, 2);

  return cMatResultVoigt;
}

template <int dim>
double MAT::PullBackFourTensorijkl(const CORE::LINALG::Matrix<dim, dim>& defgr,
    const CORE::LINALG::FourTensor<dim>& fourTensor, const int i, const int j, const int k,
    const int l)
{
  double cMatResult_ijkl(0.0);

  for (int A = 0; A < dim; ++A)
  {
    for (int B = 0; B < dim; ++B)
    {
      for (int C = 0; C < dim; ++C)
      {
        for (int D = 0; D < dim; ++D)
        {
          cMatResult_ijkl +=
              defgr(i, A) * defgr(j, B) * defgr(k, C) * defgr(l, D) * fourTensor(A, B, C, D);
        }
      }
    }
  }

  return cMatResult_ijkl;
}

template <int dim>
void MAT::SetupFourTensor(
    CORE::LINALG::FourTensor<dim>& fourTensor, const CORE::LINALG::Matrix<6, 6>& matrixVoigt)
{
#ifdef BACI_DEBUG
  if (dim != 3) dserror("Current implementation only valid for dim = 3.");
#endif
  // Setup 4-Tensor from 6x6 Voigt matrix (which has to be the representative of a 4 tensor with at
  // least minor symmetries)
  fourTensor(0, 0, 0, 0) = matrixVoigt(0, 0);  // C1111
  fourTensor(0, 0, 1, 1) = matrixVoigt(0, 1);  // C1122
  fourTensor(0, 0, 2, 2) = matrixVoigt(0, 2);  // C1133
  fourTensor(0, 0, 0, 1) = matrixVoigt(0, 3);  // C1112
  fourTensor(0, 0, 1, 0) = matrixVoigt(0, 3);  // C1121
  fourTensor(0, 0, 1, 2) = matrixVoigt(0, 4);  // C1123
  fourTensor(0, 0, 2, 1) = matrixVoigt(0, 4);  // C1132
  fourTensor(0, 0, 0, 2) = matrixVoigt(0, 5);  // C1113
  fourTensor(0, 0, 2, 0) = matrixVoigt(0, 5);  // C1131

  fourTensor(1, 1, 0, 0) = matrixVoigt(1, 0);  // C2211
  fourTensor(1, 1, 1, 1) = matrixVoigt(1, 1);  // C2222
  fourTensor(1, 1, 2, 2) = matrixVoigt(1, 2);  // C2233
  fourTensor(1, 1, 0, 1) = matrixVoigt(1, 3);  // C2212
  fourTensor(1, 1, 1, 0) = matrixVoigt(1, 3);  // C2221
  fourTensor(1, 1, 1, 2) = matrixVoigt(1, 4);  // C2223
  fourTensor(1, 1, 2, 1) = matrixVoigt(1, 4);  // C2232
  fourTensor(1, 1, 0, 2) = matrixVoigt(1, 5);  // C2213
  fourTensor(1, 1, 2, 0) = matrixVoigt(1, 5);  // C2231

  fourTensor(2, 2, 0, 0) = matrixVoigt(2, 0);  // C3311
  fourTensor(2, 2, 1, 1) = matrixVoigt(2, 1);  // C3322
  fourTensor(2, 2, 2, 2) = matrixVoigt(2, 2);  // C3333
  fourTensor(2, 2, 0, 1) = matrixVoigt(2, 3);  // C3312
  fourTensor(2, 2, 1, 0) = matrixVoigt(2, 3);  // C3321
  fourTensor(2, 2, 1, 2) = matrixVoigt(2, 4);  // C3323
  fourTensor(2, 2, 2, 1) = matrixVoigt(2, 4);  // C3332
  fourTensor(2, 2, 0, 2) = matrixVoigt(2, 5);  // C3313
  fourTensor(2, 2, 2, 0) = matrixVoigt(2, 5);  // C3331

  fourTensor(0, 1, 0, 0) = matrixVoigt(3, 0);
  fourTensor(1, 0, 0, 0) = matrixVoigt(3, 0);  // C1211 = C2111
  fourTensor(0, 1, 1, 1) = matrixVoigt(3, 1);
  fourTensor(1, 0, 1, 1) = matrixVoigt(3, 1);  // C1222 = C2122
  fourTensor(0, 1, 2, 2) = matrixVoigt(3, 2);
  fourTensor(1, 0, 2, 2) = matrixVoigt(3, 2);  // C1233 = C2133
  fourTensor(0, 1, 0, 1) = matrixVoigt(3, 3);
  fourTensor(1, 0, 0, 1) = matrixVoigt(3, 3);  // C1212 = C2112
  fourTensor(0, 1, 1, 0) = matrixVoigt(3, 3);
  fourTensor(1, 0, 1, 0) = matrixVoigt(3, 3);  // C1221 = C2121
  fourTensor(0, 1, 1, 2) = matrixVoigt(3, 4);
  fourTensor(1, 0, 1, 2) = matrixVoigt(3, 4);  // C1223 = C2123
  fourTensor(0, 1, 2, 1) = matrixVoigt(3, 4);
  fourTensor(1, 0, 2, 1) = matrixVoigt(3, 4);  // C1232 = C2132
  fourTensor(0, 1, 0, 2) = matrixVoigt(3, 5);
  fourTensor(1, 0, 0, 2) = matrixVoigt(3, 5);  // C1213 = C2113
  fourTensor(0, 1, 2, 0) = matrixVoigt(3, 5);
  fourTensor(1, 0, 2, 0) = matrixVoigt(3, 5);  // C1231 = C2131

  fourTensor(1, 2, 0, 0) = matrixVoigt(4, 0);
  fourTensor(2, 1, 0, 0) = matrixVoigt(4, 0);  // C2311 = C3211
  fourTensor(1, 2, 1, 1) = matrixVoigt(4, 1);
  fourTensor(2, 1, 1, 1) = matrixVoigt(4, 1);  // C2322 = C3222
  fourTensor(1, 2, 2, 2) = matrixVoigt(4, 2);
  fourTensor(2, 1, 2, 2) = matrixVoigt(4, 2);  // C2333 = C3233
  fourTensor(1, 2, 0, 1) = matrixVoigt(4, 3);
  fourTensor(2, 1, 0, 1) = matrixVoigt(4, 3);  // C2312 = C3212
  fourTensor(1, 2, 1, 0) = matrixVoigt(4, 3);
  fourTensor(2, 1, 1, 0) = matrixVoigt(4, 3);  // C2321 = C3221
  fourTensor(1, 2, 1, 2) = matrixVoigt(4, 4);
  fourTensor(2, 1, 1, 2) = matrixVoigt(4, 4);  // C2323 = C3223
  fourTensor(1, 2, 2, 1) = matrixVoigt(4, 4);
  fourTensor(2, 1, 2, 1) = matrixVoigt(4, 4);  // C2332 = C3232
  fourTensor(1, 2, 0, 2) = matrixVoigt(4, 5);
  fourTensor(2, 1, 0, 2) = matrixVoigt(4, 5);  // C2313 = C3213
  fourTensor(1, 2, 2, 0) = matrixVoigt(4, 5);
  fourTensor(2, 1, 2, 0) = matrixVoigt(4, 5);  // C2331 = C3231

  fourTensor(0, 2, 0, 0) = matrixVoigt(5, 0);
  fourTensor(2, 0, 0, 0) = matrixVoigt(5, 0);  // C1311 = C3111
  fourTensor(0, 2, 1, 1) = matrixVoigt(5, 1);
  fourTensor(2, 0, 1, 1) = matrixVoigt(5, 1);  // C1322 = C3122
  fourTensor(0, 2, 2, 2) = matrixVoigt(5, 2);
  fourTensor(2, 0, 2, 2) = matrixVoigt(5, 2);  // C1333 = C3133
  fourTensor(0, 2, 0, 1) = matrixVoigt(5, 3);
  fourTensor(2, 0, 0, 1) = matrixVoigt(5, 3);  // C1312 = C3112
  fourTensor(0, 2, 1, 0) = matrixVoigt(5, 3);
  fourTensor(2, 0, 1, 0) = matrixVoigt(5, 3);  // C1321 = C3121
  fourTensor(0, 2, 1, 2) = matrixVoigt(5, 4);
  fourTensor(2, 0, 1, 2) = matrixVoigt(5, 4);  // C1323 = C3123
  fourTensor(0, 2, 2, 1) = matrixVoigt(5, 4);
  fourTensor(2, 0, 2, 1) = matrixVoigt(5, 4);  // C1332 = C3132
  fourTensor(0, 2, 0, 2) = matrixVoigt(5, 5);
  fourTensor(2, 0, 0, 2) = matrixVoigt(5, 5);  // C1313 = C3113
  fourTensor(0, 2, 2, 0) = matrixVoigt(5, 5);
  fourTensor(2, 0, 2, 0) = matrixVoigt(5, 5);  // C1331 = C3131
}

template <int dim>
void MAT::Setup6x6VoigtMatrix(
    CORE::LINALG::Matrix<6, 6>& matrixVoigt, const CORE::LINALG::FourTensor<dim>& fourTensor)
{
#ifdef BACI_DEBUG
  if (dim != 3) dserror("Current implementation only valid for dim = 3.");
#endif

  // Setup 6x6 Voigt matrix from 4-Tensor
  matrixVoigt(0, 0) = fourTensor(0, 0, 0, 0);                                   // C1111
  matrixVoigt(0, 1) = fourTensor(0, 0, 1, 1);                                   // C1122
  matrixVoigt(0, 2) = fourTensor(0, 0, 2, 2);                                   // C1133
  matrixVoigt(0, 3) = 0.5 * (fourTensor(0, 0, 0, 1) + fourTensor(0, 0, 1, 0));  // 0.5*(C1112+C1121)
  matrixVoigt(0, 4) = 0.5 * (fourTensor(0, 0, 1, 2) + fourTensor(0, 0, 2, 1));  // 0.5*(C1123+C1132)
  matrixVoigt(0, 5) = 0.5 * (fourTensor(0, 0, 0, 2) + fourTensor(0, 0, 2, 0));  // 0.5*(C1113+C1131)

  matrixVoigt(1, 0) = fourTensor(1, 1, 0, 0);                                   // C2211
  matrixVoigt(1, 1) = fourTensor(1, 1, 1, 1);                                   // C2222
  matrixVoigt(1, 2) = fourTensor(1, 1, 2, 2);                                   // C2233
  matrixVoigt(1, 3) = 0.5 * (fourTensor(1, 1, 0, 1) + fourTensor(1, 1, 1, 0));  // 0.5*(C2212+C2221)
  matrixVoigt(1, 4) = 0.5 * (fourTensor(1, 1, 1, 2) + fourTensor(1, 1, 2, 1));  // 0.5*(C2223+C2232)
  matrixVoigt(1, 5) = 0.5 * (fourTensor(1, 1, 0, 2) + fourTensor(1, 1, 2, 0));  // 0.5*(C2213+C2231)

  matrixVoigt(2, 0) = fourTensor(2, 2, 0, 0);                                   // C3311
  matrixVoigt(2, 1) = fourTensor(2, 2, 1, 1);                                   // C3322
  matrixVoigt(2, 2) = fourTensor(2, 2, 2, 2);                                   // C3333
  matrixVoigt(2, 3) = 0.5 * (fourTensor(2, 2, 0, 1) + fourTensor(2, 2, 1, 0));  // 0.5*(C3312+C3321)
  matrixVoigt(2, 4) = 0.5 * (fourTensor(2, 2, 1, 2) + fourTensor(2, 2, 2, 1));  // 0.5*(C3323+C3332)
  matrixVoigt(2, 5) = 0.5 * (fourTensor(2, 2, 0, 2) + fourTensor(2, 2, 2, 0));  // 0.5*(C3313+C3331)

  matrixVoigt(3, 0) = 0.5 * (fourTensor(0, 1, 0, 0) + fourTensor(1, 0, 0, 0));  // 0.5*(C1211+C2111)
  matrixVoigt(3, 1) = 0.5 * (fourTensor(0, 1, 1, 1) + fourTensor(1, 0, 1, 1));  // 0.5*(C1222+C2122)
  matrixVoigt(3, 2) = 0.5 * (fourTensor(0, 1, 2, 2) + fourTensor(1, 0, 2, 2));  // 0.5*(C1233+C2133)
  matrixVoigt(3, 3) =
      0.25 * (fourTensor(0, 1, 0, 1) + fourTensor(1, 0, 0, 1) + fourTensor(0, 1, 1, 0) +
                 fourTensor(1, 0, 1, 0));  // 0.5*(C1212+C2112+C1221+C2121)
  matrixVoigt(3, 4) =
      0.25 * (fourTensor(0, 1, 1, 2) + fourTensor(1, 0, 1, 2) + fourTensor(0, 1, 2, 1) +
                 fourTensor(1, 0, 2, 1));  // 0.5*(C1223+C2123+C1232+C2132)
  matrixVoigt(3, 5) =
      0.25 * (fourTensor(0, 1, 0, 2) + fourTensor(1, 0, 0, 2) + fourTensor(0, 1, 2, 0) +
                 fourTensor(1, 0, 2, 0));  // 0.5*(C1213+C2113+C1231+C2131)

  matrixVoigt(4, 0) = 0.5 * (fourTensor(1, 2, 0, 0) + fourTensor(2, 1, 0, 0));  // 0.5*(C2311+C3211)
  matrixVoigt(4, 1) = 0.5 * (fourTensor(1, 2, 1, 1) + fourTensor(2, 1, 1, 1));  // 0.5*(C2322+C3222)
  matrixVoigt(4, 2) = 0.5 * (fourTensor(1, 2, 2, 2) + fourTensor(2, 1, 2, 2));  // 0.5*(C2333+C3233)
  matrixVoigt(4, 3) =
      0.25 * (fourTensor(1, 2, 0, 1) + fourTensor(2, 1, 0, 1) + fourTensor(1, 2, 1, 0) +
                 fourTensor(2, 1, 1, 0));  // 0.5*(C2312+C3212+C2321+C3221)
  matrixVoigt(4, 4) =
      0.25 * (fourTensor(1, 2, 1, 2) + fourTensor(2, 1, 1, 2) + fourTensor(1, 2, 2, 1) +
                 fourTensor(2, 1, 2, 1));  // 0.5*(C2323+C3223+C2332+C3232)
  matrixVoigt(4, 5) =
      0.25 * (fourTensor(1, 2, 0, 2) + fourTensor(2, 1, 0, 2) + fourTensor(1, 2, 2, 0) +
                 fourTensor(2, 1, 2, 0));  // 0.5*(C2313+C3213+C2331+C3231)

  matrixVoigt(5, 0) = 0.5 * (fourTensor(0, 2, 0, 0) + fourTensor(2, 0, 0, 0));  // 0.5*(C1311+C3111)
  matrixVoigt(5, 1) = 0.5 * (fourTensor(0, 2, 1, 1) + fourTensor(2, 0, 1, 1));  // 0.5*(C1322+C3122)
  matrixVoigt(5, 2) = 0.5 * (fourTensor(0, 2, 2, 2) + fourTensor(2, 0, 2, 2));  // 0.5*(C1333+C3133)
  matrixVoigt(5, 3) =
      0.25 * (fourTensor(0, 2, 0, 1) + fourTensor(2, 0, 0, 1) + fourTensor(0, 2, 1, 0) +
                 fourTensor(2, 0, 1, 0));  // 0.5*(C1312+C3112+C1321+C3121)
  matrixVoigt(5, 4) =
      0.25 * (fourTensor(0, 2, 1, 2) + fourTensor(2, 0, 1, 2) + fourTensor(0, 2, 2, 1) +
                 fourTensor(2, 0, 2, 1));  // 0.5*(C1323+C3123+C1332+C3132)
  matrixVoigt(5, 5) =
      0.25 * (fourTensor(0, 2, 0, 2) + fourTensor(2, 0, 0, 2) + fourTensor(0, 2, 2, 0) +
                 fourTensor(2, 0, 2, 0));  // 0.5*(C1313+C3113+C1331+C3131)
}

template <int dim>
void MAT::TransposeFourTensor12(
    CORE::LINALG::FourTensor<dim>& resultTensor, const CORE::LINALG::FourTensor<dim>& inputTensor)
{
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          resultTensor(i, j, k, l) = inputTensor(j, i, k, l);
        }
      }
    }
  }
}

void MAT::AddDyadicProductMatrixMatrix(CORE::LINALG::FourTensor<3>& fourTensorResult,
    const CORE::LINALG::Matrix<3, 3>& matrixA, const CORE::LINALG::Matrix<3, 3>& matrixB)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          fourTensorResult(i, j, k, l) += matrixA(i, j) * matrixB(k, l);
}

void MAT::AddContractionMatrixFourTensor(CORE::LINALG::Matrix<3, 3>& matrixResult,
    const CORE::LINALG::Matrix<3, 3>& matrix, const CORE::LINALG::FourTensor<3>& fourTensor)
{
  for (unsigned k = 0; k < 3; ++k)
    for (unsigned l = 0; l < 3; ++l)
      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
          matrixResult(k, l) += matrix(i, j) * fourTensor(i, j, k, l);
}
double MAT::ContractMatrixMatrix(
    const CORE::LINALG::Matrix<3, 3>& matrixA, const CORE::LINALG::Matrix<3, 3>& matrixB)
{
  double scalarContraction = 0.0;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) scalarContraction += matrixA(i, j) * matrixB(i, j);

  return scalarContraction;
}

// explicit instantiation of template functions
template void MAT::AddRightNonSymmetricHolzapfelProduct<double>(CORE::LINALG::Matrix<6, 9, double>&,
    CORE::LINALG::Matrix<3, 3, double> const&, CORE::LINALG::Matrix<3, 3, double> const&,
    double const);
template void MAT::AddRightNonSymmetricHolzapfelProduct<FAD>(CORE::LINALG::Matrix<6, 9, FAD>&,
    CORE::LINALG::Matrix<3, 3, FAD> const&, CORE::LINALG::Matrix<3, 3, FAD> const&, FAD const);
template void MAT::AddRightNonSymmetricHolzapfelProductStrainLike<double>(
    CORE::LINALG::Matrix<6, 9, double>& out, CORE::LINALG::Matrix<3, 3, double> const& A,
    CORE::LINALG::Matrix<3, 3, double> const& B, double const fac);
template void MAT::AddRightNonSymmetricHolzapfelProductStrainLike<FAD>(
    CORE::LINALG::Matrix<6, 9, FAD>& out, CORE::LINALG::Matrix<3, 3, FAD> const& A,
    CORE::LINALG::Matrix<3, 3, FAD> const& B, FAD const fac);
template void MAT::AddtoCmatHolzapfelProduct<double>(CORE::LINALG::Matrix<6, 6, double>&,
    const CORE::LINALG::Matrix<6, 1, double>&, const double scalar);
template void MAT::AddtoCmatHolzapfelProduct<FAD>(
    CORE::LINALG::Matrix<6, 6, FAD>&, const CORE::LINALG::Matrix<6, 1, FAD>&, const FAD scalar);

template MIXTURE::PAR::IsotropicCylinderPrestressStrategy* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::ConstantPrestressStrategy* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::IterativePrestressStrategy* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::FunctionMixtureRule* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::MapMixtureRule* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::SimpleMixtureRule* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::GrowthRemodelMixtureRule* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::IsotropicGrowthStrategy* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::AnisotropicGrowthStrategy* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::StiffnessGrowthStrategy* MAT::CreateMaterialParameterInstance(
    Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::RemodelFiberMaterialExponential<double>*
MAT::CreateMaterialParameterInstance(Teuchos::RCP<MAT::PAR::Material> curmat);
template MIXTURE::PAR::RemodelFiberMaterialExponentialActive<double>*
MAT::CreateMaterialParameterInstance(Teuchos::RCP<MAT::PAR::Material> curmat);

template void MAT::ClearFourTensor<3>(CORE::LINALG::FourTensor<3>& fourTensor);

template void MAT::MultiplyFourTensorMatrix<3>(CORE::LINALG::FourTensor<3>& fourTensorResult,
    const CORE::LINALG::FourTensor<3>& fourTensor, const CORE::LINALG::Matrix<3, 3>& matrix,
    const bool clearResultTensor);

template void MAT::MultiplyMatrixFourTensor<3>(CORE::LINALG::FourTensor<3>& fourTensorResult,
    const CORE::LINALG::Matrix<3, 3>& matrix, const CORE::LINALG::FourTensor<3>& fourTensor,
    const bool clearResultTensor);

template void MAT::MultiplyMatrixFourTensorBySecondIndex<3>(
    CORE::LINALG::FourTensor<3>& fourTensorResult, const CORE::LINALG::Matrix<3, 3>& matrix,
    const CORE::LINALG::FourTensor<3>& fourTensor, const bool clearResultTensor);

template void MAT::MultiplyFourTensorFourTensor<3>(CORE::LINALG::FourTensor<3>& fourTensorResult,
    const CORE::LINALG::FourTensor<3>& fourTensor1, const CORE::LINALG::FourTensor<3>& fourTensor2,
    const bool clearResultTensor);

template CORE::LINALG::Matrix<6, 6> MAT::PullBackFourTensor<3>(
    const CORE::LINALG::Matrix<3, 3>& defgr, const CORE::LINALG::Matrix<6, 6>& cMatVoigt);

template double MAT::PullBackFourTensorijkl<3>(const CORE::LINALG::Matrix<3, 3>& defgr,
    const CORE::LINALG::FourTensor<3>& fourTensor, const int i, const int j, const int k,
    const int l);

template void MAT::SetupFourTensor<3>(
    CORE::LINALG::FourTensor<3>& fourTensor, const CORE::LINALG::Matrix<6, 6>& matrixVoigt);

template void MAT::Setup6x6VoigtMatrix<3>(
    CORE::LINALG::Matrix<6, 6>& matrixVoigt, const CORE::LINALG::FourTensor<3>& fourTensor);

template void MAT::TransposeFourTensor12<3>(
    CORE::LINALG::FourTensor<3>& resultTensor, const CORE::LINALG::FourTensor<3>& inputTensor);
FOUR_C_NAMESPACE_CLOSE
