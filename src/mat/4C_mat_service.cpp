/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for complex materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_service.hpp"

#include "4C_linalg_four_tensor.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential_active.hpp"
#include "4C_mixture_growth_strategy_anisotropic.hpp"
#include "4C_mixture_growth_strategy_isotropic.hpp"
#include "4C_mixture_growth_strategy_stiffness.hpp"
#include "4C_mixture_prestress_strategy_constant.hpp"
#include "4C_mixture_prestress_strategy_isocyl.hpp"
#include "4C_mixture_prestress_strategy_iterative.hpp"
#include "4C_mixture_rule_function.hpp"
#include "4C_mixture_rule_growthremodel.hpp"
#include "4C_mixture_rule_map.hpp"
#include "4C_mixture_rule_simple.hpp"

#include <Sacado.hpp>

using FAD = Sacado::Fad::DFad<double>;

FOUR_C_NAMESPACE_OPEN

template <typename T>
void Mat::add_holzapfel_product(
    Core::LinAlg::Matrix<6, 6, T>& cmat, const Core::LinAlg::Matrix<6, 1, T>& invc, const T scalar)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (cmat.num_rows() != 6 or cmat.num_cols() != 6 or invc.num_rows() != 6)
    FOUR_C_THROW("Wrong dimensions in function add_holzapfel_product");
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

double Mat::second_invariant_of_deviatoric_stress(const Core::LinAlg::Matrix<3, 3>& stress)
{
  const double p = 1.0 / 3 * (stress(0, 0) + stress(1, 1) + stress(2, 2));
  const double s11 = stress(0, 0) - p;
  const double s22 = stress(1, 1) - p;
  const double s33 = stress(2, 2) - p;
  const double s12 = stress(0, 1);
  const double s23 = stress(1, 2);
  const double s13 = stress(0, 2);
  const double J2 =
      0.5 * (s11 * s11 + s22 * s22 + s33 * s33 + 2 * (s12 * s12 + s23 * s23 + s13 * s13));
  return J2;
}

void Mat::add_elasticity_tensor_product(Core::LinAlg::Matrix<6, 6>& C, const double scalar_AB,
    const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B,
    const double scalar_this)
{
  // everything in Voigt-Notation
  Core::LinAlg::Matrix<6, 1> A_voigt;
  Core::LinAlg::Matrix<6, 1> B_voigt;

  A_voigt(0, 0) = A(0, 0);
  A_voigt(1, 0) = A(1, 1);
  A_voigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  A_voigt(3, 0) = A(1, 0);
  A_voigt(4, 0) = A(2, 1);
  A_voigt(5, 0) = A(2, 0);

  B_voigt(0, 0) = B(0, 0);
  B_voigt(1, 0) = B(1, 1);
  B_voigt(2, 0) = B(2, 2);
  B_voigt(3, 0) = B(1, 0);
  B_voigt(4, 0) = B(2, 1);
  B_voigt(5, 0) = B(2, 0);

  C.multiply_nt(scalar_AB, A_voigt, B_voigt, scalar_this);
}

void Mat::add_symmetric_elasticity_tensor_product(Core::LinAlg::Matrix<6, 6>& C,
    const double scalar_AB, const Core::LinAlg::Matrix<3, 3>& A,
    const Core::LinAlg::Matrix<3, 3>& B, const double scalar_this)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check sizes
  if (A.num_rows() != A.num_cols() || B.num_rows() != B.num_cols() || A.num_rows() != 3 ||
      B.num_rows() != 3)
  {
    FOUR_C_THROW("2nd order tensors must be 3 by 3");
  }
  if (C.num_rows() != C.num_cols() || C.num_rows() != 6)
    FOUR_C_THROW("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  Core::LinAlg::Matrix<6, 1> A_voigt;
  Core::LinAlg::Matrix<6, 1> B_voigt;

  A_voigt(0, 0) = A(0, 0);
  A_voigt(1, 0) = A(1, 1);
  A_voigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  A_voigt(3, 0) = A(1, 0);
  A_voigt(4, 0) = A(2, 1);
  A_voigt(5, 0) = A(2, 0);

  B_voigt(0, 0) = B(0, 0);
  B_voigt(1, 0) = B(1, 1);
  B_voigt(2, 0) = B(2, 2);
  B_voigt(3, 0) = B(1, 0);
  B_voigt(4, 0) = B(2, 1);
  B_voigt(5, 0) = B(2, 0);

  C.multiply_nt(scalar_AB, A_voigt, B_voigt, scalar_this);
  C.multiply_nt(scalar_AB, B_voigt, A_voigt, 1.0);
}

void Mat::add_kronecker_tensor_product(Core::LinAlg::Matrix<6, 6>& C, const double scalar_AB,
    const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B,
    const double scalar_this)
{
  const double scalar_AB_half = scalar_AB * 0.5;
  C(0, 0) = scalar_this * C(0, 0) + scalar_AB * (A(0, 0) * B(0, 0));  // C1111
  C(0, 1) = scalar_this * C(0, 1) + scalar_AB * (A(0, 1) * B(0, 1));  // C1122
  C(0, 2) = scalar_this * C(0, 2) + scalar_AB * (A(0, 2) * B(0, 2));  // C1133
  C(0, 3) =
      scalar_this * C(0, 3) + scalar_AB_half * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));  // C1112
  C(0, 4) =
      scalar_this * C(0, 4) + scalar_AB_half * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));  // C1123
  C(0, 5) =
      scalar_this * C(0, 5) + scalar_AB_half * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));  // C1113

  C(1, 0) = scalar_this * C(1, 0) + scalar_AB * (A(1, 0) * B(1, 0));  // C2211
  C(1, 1) = scalar_this * C(1, 1) + scalar_AB * (A(1, 1) * B(1, 1));  // C2222
  C(1, 2) = scalar_this * C(1, 2) + scalar_AB * (A(1, 2) * B(1, 2));  // C2233
  C(1, 3) =
      scalar_this * C(1, 3) + scalar_AB_half * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));  // C2212
  C(1, 4) =
      scalar_this * C(1, 4) + scalar_AB_half * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));  // C2223
  C(1, 5) =
      scalar_this * C(1, 5) + scalar_AB_half * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));  // C2213

  C(2, 0) = scalar_this * C(2, 0) + scalar_AB * (A(2, 0) * B(2, 0));  // C3311
  C(2, 1) = scalar_this * C(2, 1) + scalar_AB * (A(2, 1) * B(2, 1));  // C3322
  C(2, 2) = scalar_this * C(2, 2) + scalar_AB * (A(2, 2) * B(2, 2));  // C3333
  C(2, 3) =
      scalar_this * C(2, 3) + scalar_AB_half * (A(2, 1) * B(2, 1) + A(2, 1) * B(2, 0));  // C3312
  C(2, 4) =
      scalar_this * C(2, 4) + scalar_AB_half * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));  // C3323
  C(2, 5) =
      scalar_this * C(2, 5) + scalar_AB_half * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));  // C3313

  C(3, 0) = scalar_this * C(3, 0) + scalar_AB * (A(0, 0) * B(1, 0));  // C1211
  C(3, 1) = scalar_this * C(3, 1) + scalar_AB * (A(0, 1) * B(1, 1));  // C1222
  C(3, 2) = scalar_this * C(3, 2) + scalar_AB * (A(0, 2) * B(1, 2));  // C1233
  C(3, 3) =
      scalar_this * C(3, 3) + scalar_AB_half * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));  // C1212
  C(3, 4) =
      scalar_this * C(3, 4) + scalar_AB_half * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));  // C1223
  C(3, 5) =
      scalar_this * C(3, 5) + scalar_AB_half * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));  // C1213

  C(4, 0) = scalar_this * C(4, 0) + scalar_AB * (A(1, 0) * B(2, 0));  // C2311
  C(4, 1) = scalar_this * C(4, 1) + scalar_AB * (A(1, 1) * B(2, 1));  // C2322
  C(4, 2) = scalar_this * C(4, 2) + scalar_AB * (A(1, 2) * B(2, 2));  // C2333
  C(4, 3) =
      scalar_this * C(4, 3) + scalar_AB_half * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));  // C2312
  C(4, 4) =
      scalar_this * C(4, 4) + scalar_AB_half * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));  // C2323
  C(4, 5) =
      scalar_this * C(4, 5) + scalar_AB_half * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));  // C2313

  C(5, 0) = scalar_this * C(5, 0) + scalar_AB * (A(0, 0) * B(2, 0));  // C1311
  C(5, 1) = scalar_this * C(5, 1) + scalar_AB * (A(0, 1) * B(2, 1));  // C1322
  C(5, 2) = scalar_this * C(5, 2) + scalar_AB * (A(0, 2) * B(2, 2));  // C1333
  C(5, 3) =
      scalar_this * C(5, 3) + scalar_AB_half * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));  // C1312
  C(5, 4) =
      scalar_this * C(5, 4) + scalar_AB_half * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));  // C1323
  C(5, 5) =
      scalar_this * C(5, 5) + scalar_AB_half * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313
}

void Mat::volumetrify_and_isochorify(Core::LinAlg::Matrix<6, 1>* pk2vol,
    Core::LinAlg::Matrix<6, 6>* cvol, Core::LinAlg::Matrix<6, 1>* pk2iso,
    Core::LinAlg::Matrix<6, 6>* ciso, const Core::LinAlg::Matrix<6, 1>& gl,
    const Core::LinAlg::Matrix<6, 1>& pk2, const Core::LinAlg::Matrix<6, 6>& cmat)
{
  // useful call?
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if ((pk2vol == nullptr) and (cvol == nullptr) and (pk2iso == nullptr) and (ciso == nullptr))
    FOUR_C_THROW("Useful call? Apparently you do not want to compute anything");
#endif

  // right Cauchy--Green tensor
  // REMARK: stored in _strain_-like 6-Voigt vector
  Core::LinAlg::Matrix<6, 1> rcg(gl);
  rcg.scale(2.0);
  for (int i = 0; i < 3; i++) rcg(i) += 1.0;

  // third invariant (determinant) of right Cauchy--Green strains
  const double rcg3rd = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                        0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                        0.25 * rcg(0) * rcg(4) * rcg(4);

  // inverse right Cauchy--Green tensor C^{-1}
  // REMARK: stored in as _stress_ 6-Voigt vector
  Core::LinAlg::Matrix<6, 1> icg(false);
  icg(0) = (rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4)) / rcg3rd;        // (C^{-1})^{11}
  icg(1) = (rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5)) / rcg3rd;        // (C^{-1})^{22}
  icg(2) = (rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3)) / rcg3rd;        // (C^{-1})^{33}
  icg(3) = (0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2)) / rcg3rd;  // (C^{-1})^{12}
  icg(4) = (0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4)) / rcg3rd;  // (C^{-1})^{23}
  icg(5) = (0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1)) / rcg3rd;  // (C^{-1})^{31}

  // double contraction of 2nd Piola--Kirchhoff stress and right Cauchy--Green strain,
  // i.e. in index notation S^{AB} C_{AB}
  // REMARK: equal to S^T C, because S is stress-like and C is strain-like 6-Voigt vector
  const double pk2rcg = pk2.dot(rcg);

  // stress splitting
  {
    // volumetric 2nd Piola--Kirchhoff stress
    Core::LinAlg::Matrix<6, 1> pk2vol_tmp(false);
    if (pk2vol != nullptr) pk2vol_tmp.set_view(*pk2vol);
    pk2vol_tmp.update(pk2rcg / 3.0, icg);

    // isochoric 2nd Piola--Kirchhoff stress
    // S^{AB}_iso = S^{AB} - S^{AB}_{vol}
    if (pk2iso != nullptr) pk2iso->update(1.0, pk2, -1.0, pk2vol_tmp);
  }

  // elasticity tensor splitting
  {
    // 'linearised' 2nd Piola--Kirchhoff stress
    // S^{CD}_lin = S^{CD} + 1/2 C_{AB} C^{ABCD}
    Core::LinAlg::Matrix<6, 1> pk2lin(pk2);
    pk2lin.multiply_tn(0.5, cmat, rcg, 1.0);  // transpose on purpose

    // volumetric part of constitutive tensor
    // C^{ABCD}_vol = 2/3 (C^{-1})^{AB} S^{CD}_lin
    //              - 2/3 (S^{EF} C_{EF}) ( 1/2 (
    //                (C^{-1})^{AC} (C^{-1})^{BD} + (C^{-1})^{AD} (C^{-1})^{BC}
    //              ) )
    Core::LinAlg::Matrix<6, 6> cvol_tmp(false);
    if (cvol != nullptr) cvol_tmp.set_view(*cvol);
    cvol_tmp.multiply_nt(2.0 / 3.0, icg, pk2lin);
    add_holzapfel_product(cvol_tmp, icg, -2.0 / 3.0 * pk2rcg);

    // isochoric part of constitutive tensor
    // C^{ABCD}_iso = C^{ABCD} - C^{ABCD}_vol
    if (ciso != nullptr) ciso->update(1.0, cmat, -1.0, cvol_tmp);
  }
}

void Mat::add_derivative_of_squared_tensor(Core::LinAlg::Matrix<6, 6>& C, double scalar_squared_dx,
    Core::LinAlg::Matrix<3, 3> X, double scalar_this)
{
  C(0, 0) = scalar_this * C(0, 0) + scalar_squared_dx * 2. * X(0, 0);  // C1111
  C(0, 1) = scalar_this * C(0, 1);                                     // C1122
  C(0, 2) = scalar_this * C(0, 2);                                     // C1133
  C(0, 3) = scalar_this * C(0, 3) + scalar_squared_dx * X(0, 1);       // C1112
  C(0, 4) = scalar_this * C(0, 4);                                     // C1123
  C(0, 5) = scalar_this * C(0, 5) + scalar_squared_dx * X(0, 2);       // C1113

  C(1, 0) = scalar_this * C(1, 0);                                     // C2211
  C(1, 1) = scalar_this * C(1, 1) + scalar_squared_dx * 2. * X(1, 1);  // C2222
  C(1, 2) = scalar_this * C(1, 2);                                     // C2233
  C(1, 3) = scalar_this * C(1, 3) + scalar_squared_dx * X(0, 1);       // C2212
  C(1, 4) = scalar_this * C(1, 4) + scalar_squared_dx * X(1, 2);       // C2223
  C(1, 5) = scalar_this * C(1, 5);                                     // C2213

  C(2, 0) = scalar_this * C(2, 0);                                     // C3311
  C(2, 1) = scalar_this * C(2, 1);                                     // C3322
  C(2, 2) = scalar_this * C(2, 2) + scalar_squared_dx * 2. * X(2, 2);  // C3333
  C(2, 3) = scalar_this * C(2, 3);                                     // C3312
  C(2, 4) = scalar_this * C(2, 4) + scalar_squared_dx * X(1, 2);       // C3323
  C(2, 5) = scalar_this * C(2, 5) + scalar_squared_dx * X(0, 2);       // C3313

  C(3, 0) = scalar_this * C(3, 0) + scalar_squared_dx * X(0, 1);                    // C1211
  C(3, 1) = scalar_this * C(3, 1) + scalar_squared_dx * X(0, 1);                    // C1222
  C(3, 2) = scalar_this * C(3, 2);                                                  // C1233
  C(3, 3) = scalar_this * C(3, 3) + scalar_squared_dx * 0.5 * (X(0, 0) + X(1, 1));  // C1212
  C(3, 4) = scalar_this * C(3, 4) + scalar_squared_dx * 0.5 * X(0, 2);              // C1223
  C(3, 5) = scalar_this * C(3, 5) + scalar_squared_dx * 0.5 * X(1, 2);              // C1213

  C(4, 0) = scalar_this * C(4, 0);                                                  // C2311
  C(4, 1) = scalar_this * C(4, 1) + scalar_squared_dx * X(1, 2);                    // C2322
  C(4, 2) = scalar_this * C(4, 2) + scalar_squared_dx * X(1, 2);                    // C2333
  C(4, 3) = scalar_this * C(4, 3) + scalar_squared_dx * 0.5 * X(0, 2);              // C2312
  C(4, 4) = scalar_this * C(4, 4) + scalar_squared_dx * 0.5 * (X(1, 1) + X(2, 2));  // C2323
  C(4, 5) = scalar_this * C(4, 5) + scalar_squared_dx * 0.5 * X(0, 1);              // C2313

  C(5, 0) = scalar_this * C(5, 0) + scalar_squared_dx * X(0, 2);                    // C1311
  C(5, 1) = scalar_this * C(5, 1);                                                  // C1322
  C(5, 2) = scalar_this * C(5, 2) + scalar_squared_dx * X(0, 2);                    // C1333
  C(5, 3) = scalar_this * C(5, 3) + scalar_squared_dx * 0.5 * X(1, 2);              // C1312
  C(5, 4) = scalar_this * C(5, 4) + scalar_squared_dx * 0.5 * X(0, 1);              // C1323
  C(5, 5) = scalar_this * C(5, 5) + scalar_squared_dx * 0.5 * (X(2, 2) + X(0, 0));  // C1313
}

void Mat::add_symmetric_holzapfel_product(Core::LinAlg::Matrix<6, 6>& X,
    const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B, const double fac)
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
void Mat::add_right_non_symmetric_holzapfel_product(Core::LinAlg::Matrix<6, 9, T>& out,
    Core::LinAlg::Matrix<3, 3, T> const& A, Core::LinAlg::Matrix<3, 3, T> const& B, T const fac)
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
void Mat::add_right_non_symmetric_holzapfel_product_strain_like(Core::LinAlg::Matrix<6, 9, T>& out,
    Core::LinAlg::Matrix<3, 3, T> const& A, Core::LinAlg::Matrix<3, 3, T> const& B, T const fac)
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

void Mat::add_left_non_symmetric_holzapfel_product(Core::LinAlg::Matrix<9, 6>& out,
    Core::LinAlg::Matrix<3, 3> const& A, Core::LinAlg::Matrix<3, 3> const& B, double const fac)
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

void Mat::add_non_symmetric_product(double const& fac, Core::LinAlg::Matrix<3, 3> const& A,
    Core::LinAlg::Matrix<3, 3> const& B, Core::LinAlg::Matrix<9, 9>& out)
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

void Mat::add_derivative_of_inva_b_inva_product(double const& fac,
    const Core::LinAlg::Matrix<6, 1>& invA, const Core::LinAlg::Matrix<6, 1>& invABinvA,
    Core::LinAlg::Matrix<6, 6>& out)
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

void Mat::invariants_principal(
    Core::LinAlg::Matrix<3, 1>& prinv, const Core::LinAlg::Matrix<3, 3>& tens)
{
  // 1st invariant, trace tens
  prinv(0) = tens(0, 0) + tens(1, 1) + tens(2, 2);

  // 2nd invariant, 0.5( (trace(tens))^2 - trace(tens^2))
  prinv(1) = tens(0, 0) * tens(1, 1) + tens(1, 1) * tens(2, 2) + tens(0, 0) * tens(2, 2) -
             tens(0, 1) * tens(1, 0) - tens(1, 2) * tens(2, 1) - tens(0, 2) * tens(2, 0);

  // 3rd invariant, determinant tens
  prinv(2) = tens.determinant();
}

void Mat::invariants_modified(
    Core::LinAlg::Matrix<3, 1>& modinv, const Core::LinAlg::Matrix<3, 1>& prinv)
{
  // 1st invariant, trace
  modinv(0) = prinv(0) * std::pow(prinv(2), -1. / 3.);
  // 2nd invariant
  modinv(1) = prinv(1) * std::pow(prinv(2), -2. / 3.);
  // J
  modinv(2) = std::pow(prinv(2), 1. / 2.);
}

void Mat::stretches_principal(Core::LinAlg::Matrix<3, 1>& prstr, Core::LinAlg::Matrix<3, 3>& prdir,
    const Core::LinAlg::Matrix<6, 1>& rcg)
{
  // create right Cauchy-Green 2-tensor
  Core::LinAlg::Matrix<3, 3> rcgt(false);
  rcgt(0, 0) = rcg(0);
  rcgt(1, 1) = rcg(1);
  rcgt(2, 2) = rcg(2);
  rcgt(0, 1) = rcgt(1, 0) = 0.5 * rcg(3);
  rcgt(1, 2) = rcgt(2, 1) = 0.5 * rcg(4);
  rcgt(2, 0) = rcgt(0, 2) = 0.5 * rcg(5);

  // eigenvalue decomposition
  Core::LinAlg::Matrix<3, 3> prstr2;  // squared principal stretches
  Core::LinAlg::SYEV(rcgt, prstr2, prdir);

  // THE principal stretches
  for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));
}

void Mat::stretches_modified(
    Core::LinAlg::Matrix<3, 1>& modstr, const Core::LinAlg::Matrix<3, 1>& prstr)
{
  // determinant of deformation gradient
  const double detdefgrad = prstr(0) * prstr(1) * prstr(2);

  // determine modified principal stretches
  modstr.update(std::pow(detdefgrad, -1.0 / 3.0), prstr);
}

template <int dim>
void Mat::clear_four_tensor(Core::LinAlg::FourTensor<dim>& four_tensor)
{
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          four_tensor(i, j, k, l) = 0.0;
        }
      }
    }
  }
}

template <int dim>
void Mat::multiply_four_tensor_matrix(Core::LinAlg::FourTensor<dim>& four_tensor_result,
    const Core::LinAlg::FourTensor<dim>& four_tensor, const Core::LinAlg::Matrix<dim, dim>& matrix,
    const bool clear_result_tensor)
{
  if (clear_result_tensor) clear_four_tensor(four_tensor_result);
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
            four_tensor_result(i, j, k, l) += four_tensor(i, j, k, m) * matrix(m, l);
          }
        }
      }
    }
  }
}

template <int dim>
void Mat::multiply_matrix_four_tensor(Core::LinAlg::FourTensor<dim>& four_tensor_result,
    const Core::LinAlg::Matrix<dim, dim>& matrix, const Core::LinAlg::FourTensor<dim>& four_tensor,
    const bool clear_result_tensor)
{
  if (clear_result_tensor) clear_four_tensor(four_tensor_result);
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
            four_tensor_result(i, j, k, l) += matrix(i, m) * four_tensor(m, j, k, l);
          }
        }
      }
    }
  }
}

template <int dim>
void Mat::multiply_matrix_four_tensor_by_second_index(
    Core::LinAlg::FourTensor<dim>& four_tensor_result, const Core::LinAlg::Matrix<dim, dim>& matrix,
    const Core::LinAlg::FourTensor<dim>& four_tensor, const bool clear_result_tensor)
{
  if (clear_result_tensor) clear_four_tensor(four_tensor_result);
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
            four_tensor_result(i, j, k, l) += matrix(m, j) * four_tensor(i, m, k, l);
          }
        }
      }
    }
  }
}

template <int dim>
void Mat::multiply_four_tensor_four_tensor(Core::LinAlg::FourTensor<dim>& four_tensor_result,
    const Core::LinAlg::FourTensor<dim>& four_tensor_1,
    const Core::LinAlg::FourTensor<dim>& four_tensor_2, const bool clear_result_tensor)
{
  if (clear_result_tensor) clear_four_tensor(four_tensor_result);
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
              four_tensor_result(i, j, k, l) +=
                  four_tensor_1(i, j, a, b) * four_tensor_2(a, b, k, l);
        }
      }
    }
  }
}

template <int dim>
Core::LinAlg::Matrix<6, 6> Mat::pull_back_four_tensor(
    const Core::LinAlg::Matrix<dim, dim>& defgrd, const Core::LinAlg::Matrix<6, 6>& cmat_voigt)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (dim != 3) FOUR_C_THROW("Current implementation only valid for dim = 3.");
#endif

  Core::LinAlg::FourTensor<dim> cmat_tensor(true);
  setup_four_tensor_from_6x6_voigt_matrix(cmat_tensor, cmat_voigt);

  // We can use the fact that cmat_result_voigt(i,j,k,l)=cmat_result_voigt(k,l,i,j) if we have a
  // hyper-elastic material
  Core::LinAlg::Matrix<6, 6> cmat_result_voigt(true);

  cmat_result_voigt(0, 0) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 0, 0, 0);
  cmat_result_voigt(0, 1) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 0, 1, 1);
  cmat_result_voigt(0, 2) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 0, 2, 2);
  cmat_result_voigt(0, 3) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 0, 0, 1);
  cmat_result_voigt(0, 4) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 0, 1, 2);
  cmat_result_voigt(0, 5) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 0, 0, 2);
  cmat_result_voigt(1, 0) = cmat_result_voigt(0, 1);
  cmat_result_voigt(1, 1) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 1, 1, 1);
  cmat_result_voigt(1, 2) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 1, 2, 2);
  cmat_result_voigt(1, 3) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 1, 0, 1);
  cmat_result_voigt(1, 4) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 1, 1, 2);
  cmat_result_voigt(1, 5) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 1, 0, 2);
  cmat_result_voigt(2, 0) = cmat_result_voigt(0, 2);
  cmat_result_voigt(2, 1) = cmat_result_voigt(1, 2);
  cmat_result_voigt(2, 2) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 2, 2, 2, 2);
  cmat_result_voigt(2, 3) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 2, 2, 0, 1);
  cmat_result_voigt(2, 4) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 2, 2, 1, 2);
  cmat_result_voigt(2, 5) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 2, 2, 0, 2);
  cmat_result_voigt(3, 0) = cmat_result_voigt(0, 3);
  cmat_result_voigt(3, 1) = cmat_result_voigt(1, 3);
  cmat_result_voigt(3, 2) = cmat_result_voigt(2, 3);
  cmat_result_voigt(3, 3) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 1, 0, 1);
  cmat_result_voigt(3, 4) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 1, 1, 2);
  cmat_result_voigt(3, 5) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 1, 0, 2);
  cmat_result_voigt(4, 0) = cmat_result_voigt(0, 4);
  cmat_result_voigt(4, 1) = cmat_result_voigt(1, 4);
  cmat_result_voigt(4, 2) = cmat_result_voigt(2, 4);
  cmat_result_voigt(4, 3) = cmat_result_voigt(3, 4);
  cmat_result_voigt(4, 4) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 2, 1, 2);
  cmat_result_voigt(4, 5) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 1, 2, 0, 2);
  cmat_result_voigt(5, 0) = cmat_result_voigt(0, 5);
  cmat_result_voigt(5, 1) = cmat_result_voigt(1, 5);
  cmat_result_voigt(5, 2) = cmat_result_voigt(2, 5);
  cmat_result_voigt(5, 3) = cmat_result_voigt(3, 5);
  cmat_result_voigt(5, 4) = cmat_result_voigt(4, 5);
  cmat_result_voigt(5, 5) = get_pull_back_four_tensor_entry<dim>(defgrd, cmat_tensor, 0, 2, 0, 2);

  return cmat_result_voigt;
}

template <int dim>
double Mat::get_pull_back_four_tensor_entry(const Core::LinAlg::Matrix<dim, dim>& defgrd,
    const Core::LinAlg::FourTensor<dim>& four_tensor, const int i, const int j, const int k,
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
              defgrd(i, A) * defgrd(j, B) * defgrd(k, C) * defgrd(l, D) * four_tensor(A, B, C, D);
        }
      }
    }
  }

  return cMatResult_ijkl;
}

template <int dim>
void Mat::setup_four_tensor_from_6x6_voigt_matrix(
    Core::LinAlg::FourTensor<dim>& four_tensor, const Core::LinAlg::Matrix<6, 6>& matrix_voigt)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (dim != 3) FOUR_C_THROW("Current implementation only valid for dim = 3.");
#endif
  // Setup 4-Tensor from 6x6 Voigt matrix (which has to be the representative of a 4 tensor with at
  // least minor symmetries)
  four_tensor(0, 0, 0, 0) = matrix_voigt(0, 0);  // C1111
  four_tensor(0, 0, 1, 1) = matrix_voigt(0, 1);  // C1122
  four_tensor(0, 0, 2, 2) = matrix_voigt(0, 2);  // C1133
  four_tensor(0, 0, 0, 1) = matrix_voigt(0, 3);  // C1112
  four_tensor(0, 0, 1, 0) = matrix_voigt(0, 3);  // C1121
  four_tensor(0, 0, 1, 2) = matrix_voigt(0, 4);  // C1123
  four_tensor(0, 0, 2, 1) = matrix_voigt(0, 4);  // C1132
  four_tensor(0, 0, 0, 2) = matrix_voigt(0, 5);  // C1113
  four_tensor(0, 0, 2, 0) = matrix_voigt(0, 5);  // C1131

  four_tensor(1, 1, 0, 0) = matrix_voigt(1, 0);  // C2211
  four_tensor(1, 1, 1, 1) = matrix_voigt(1, 1);  // C2222
  four_tensor(1, 1, 2, 2) = matrix_voigt(1, 2);  // C2233
  four_tensor(1, 1, 0, 1) = matrix_voigt(1, 3);  // C2212
  four_tensor(1, 1, 1, 0) = matrix_voigt(1, 3);  // C2221
  four_tensor(1, 1, 1, 2) = matrix_voigt(1, 4);  // C2223
  four_tensor(1, 1, 2, 1) = matrix_voigt(1, 4);  // C2232
  four_tensor(1, 1, 0, 2) = matrix_voigt(1, 5);  // C2213
  four_tensor(1, 1, 2, 0) = matrix_voigt(1, 5);  // C2231

  four_tensor(2, 2, 0, 0) = matrix_voigt(2, 0);  // C3311
  four_tensor(2, 2, 1, 1) = matrix_voigt(2, 1);  // C3322
  four_tensor(2, 2, 2, 2) = matrix_voigt(2, 2);  // C3333
  four_tensor(2, 2, 0, 1) = matrix_voigt(2, 3);  // C3312
  four_tensor(2, 2, 1, 0) = matrix_voigt(2, 3);  // C3321
  four_tensor(2, 2, 1, 2) = matrix_voigt(2, 4);  // C3323
  four_tensor(2, 2, 2, 1) = matrix_voigt(2, 4);  // C3332
  four_tensor(2, 2, 0, 2) = matrix_voigt(2, 5);  // C3313
  four_tensor(2, 2, 2, 0) = matrix_voigt(2, 5);  // C3331

  four_tensor(0, 1, 0, 0) = matrix_voigt(3, 0);
  four_tensor(1, 0, 0, 0) = matrix_voigt(3, 0);  // C1211 = C2111
  four_tensor(0, 1, 1, 1) = matrix_voigt(3, 1);
  four_tensor(1, 0, 1, 1) = matrix_voigt(3, 1);  // C1222 = C2122
  four_tensor(0, 1, 2, 2) = matrix_voigt(3, 2);
  four_tensor(1, 0, 2, 2) = matrix_voigt(3, 2);  // C1233 = C2133
  four_tensor(0, 1, 0, 1) = matrix_voigt(3, 3);
  four_tensor(1, 0, 0, 1) = matrix_voigt(3, 3);  // C1212 = C2112
  four_tensor(0, 1, 1, 0) = matrix_voigt(3, 3);
  four_tensor(1, 0, 1, 0) = matrix_voigt(3, 3);  // C1221 = C2121
  four_tensor(0, 1, 1, 2) = matrix_voigt(3, 4);
  four_tensor(1, 0, 1, 2) = matrix_voigt(3, 4);  // C1223 = C2123
  four_tensor(0, 1, 2, 1) = matrix_voigt(3, 4);
  four_tensor(1, 0, 2, 1) = matrix_voigt(3, 4);  // C1232 = C2132
  four_tensor(0, 1, 0, 2) = matrix_voigt(3, 5);
  four_tensor(1, 0, 0, 2) = matrix_voigt(3, 5);  // C1213 = C2113
  four_tensor(0, 1, 2, 0) = matrix_voigt(3, 5);
  four_tensor(1, 0, 2, 0) = matrix_voigt(3, 5);  // C1231 = C2131

  four_tensor(1, 2, 0, 0) = matrix_voigt(4, 0);
  four_tensor(2, 1, 0, 0) = matrix_voigt(4, 0);  // C2311 = C3211
  four_tensor(1, 2, 1, 1) = matrix_voigt(4, 1);
  four_tensor(2, 1, 1, 1) = matrix_voigt(4, 1);  // C2322 = C3222
  four_tensor(1, 2, 2, 2) = matrix_voigt(4, 2);
  four_tensor(2, 1, 2, 2) = matrix_voigt(4, 2);  // C2333 = C3233
  four_tensor(1, 2, 0, 1) = matrix_voigt(4, 3);
  four_tensor(2, 1, 0, 1) = matrix_voigt(4, 3);  // C2312 = C3212
  four_tensor(1, 2, 1, 0) = matrix_voigt(4, 3);
  four_tensor(2, 1, 1, 0) = matrix_voigt(4, 3);  // C2321 = C3221
  four_tensor(1, 2, 1, 2) = matrix_voigt(4, 4);
  four_tensor(2, 1, 1, 2) = matrix_voigt(4, 4);  // C2323 = C3223
  four_tensor(1, 2, 2, 1) = matrix_voigt(4, 4);
  four_tensor(2, 1, 2, 1) = matrix_voigt(4, 4);  // C2332 = C3232
  four_tensor(1, 2, 0, 2) = matrix_voigt(4, 5);
  four_tensor(2, 1, 0, 2) = matrix_voigt(4, 5);  // C2313 = C3213
  four_tensor(1, 2, 2, 0) = matrix_voigt(4, 5);
  four_tensor(2, 1, 2, 0) = matrix_voigt(4, 5);  // C2331 = C3231

  four_tensor(0, 2, 0, 0) = matrix_voigt(5, 0);
  four_tensor(2, 0, 0, 0) = matrix_voigt(5, 0);  // C1311 = C3111
  four_tensor(0, 2, 1, 1) = matrix_voigt(5, 1);
  four_tensor(2, 0, 1, 1) = matrix_voigt(5, 1);  // C1322 = C3122
  four_tensor(0, 2, 2, 2) = matrix_voigt(5, 2);
  four_tensor(2, 0, 2, 2) = matrix_voigt(5, 2);  // C1333 = C3133
  four_tensor(0, 2, 0, 1) = matrix_voigt(5, 3);
  four_tensor(2, 0, 0, 1) = matrix_voigt(5, 3);  // C1312 = C3112
  four_tensor(0, 2, 1, 0) = matrix_voigt(5, 3);
  four_tensor(2, 0, 1, 0) = matrix_voigt(5, 3);  // C1321 = C3121
  four_tensor(0, 2, 1, 2) = matrix_voigt(5, 4);
  four_tensor(2, 0, 1, 2) = matrix_voigt(5, 4);  // C1323 = C3123
  four_tensor(0, 2, 2, 1) = matrix_voigt(5, 4);
  four_tensor(2, 0, 2, 1) = matrix_voigt(5, 4);  // C1332 = C3132
  four_tensor(0, 2, 0, 2) = matrix_voigt(5, 5);
  four_tensor(2, 0, 0, 2) = matrix_voigt(5, 5);  // C1313 = C3113
  four_tensor(0, 2, 2, 0) = matrix_voigt(5, 5);
  four_tensor(2, 0, 2, 0) = matrix_voigt(5, 5);  // C1331 = C3131
}

template <int dim>
void Mat::setup_6x6_voigt_matrix_from_four_tensor(
    Core::LinAlg::Matrix<6, 6>& matrix_voigt, const Core::LinAlg::FourTensor<dim>& four_tensor)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (dim != 3) FOUR_C_THROW("Current implementation only valid for dim = 3.");
#endif

  // Setup 6x6 Voigt matrix from 4-Tensor
  matrix_voigt(0, 0) = four_tensor(0, 0, 0, 0);  // C1111
  matrix_voigt(0, 1) = four_tensor(0, 0, 1, 1);  // C1122
  matrix_voigt(0, 2) = four_tensor(0, 0, 2, 2);  // C1133
  matrix_voigt(0, 3) =
      0.5 * (four_tensor(0, 0, 0, 1) + four_tensor(0, 0, 1, 0));  // 0.5*(C1112+C1121)
  matrix_voigt(0, 4) =
      0.5 * (four_tensor(0, 0, 1, 2) + four_tensor(0, 0, 2, 1));  // 0.5*(C1123+C1132)
  matrix_voigt(0, 5) =
      0.5 * (four_tensor(0, 0, 0, 2) + four_tensor(0, 0, 2, 0));  // 0.5*(C1113+C1131)

  matrix_voigt(1, 0) = four_tensor(1, 1, 0, 0);  // C2211
  matrix_voigt(1, 1) = four_tensor(1, 1, 1, 1);  // C2222
  matrix_voigt(1, 2) = four_tensor(1, 1, 2, 2);  // C2233
  matrix_voigt(1, 3) =
      0.5 * (four_tensor(1, 1, 0, 1) + four_tensor(1, 1, 1, 0));  // 0.5*(C2212+C2221)
  matrix_voigt(1, 4) =
      0.5 * (four_tensor(1, 1, 1, 2) + four_tensor(1, 1, 2, 1));  // 0.5*(C2223+C2232)
  matrix_voigt(1, 5) =
      0.5 * (four_tensor(1, 1, 0, 2) + four_tensor(1, 1, 2, 0));  // 0.5*(C2213+C2231)

  matrix_voigt(2, 0) = four_tensor(2, 2, 0, 0);  // C3311
  matrix_voigt(2, 1) = four_tensor(2, 2, 1, 1);  // C3322
  matrix_voigt(2, 2) = four_tensor(2, 2, 2, 2);  // C3333
  matrix_voigt(2, 3) =
      0.5 * (four_tensor(2, 2, 0, 1) + four_tensor(2, 2, 1, 0));  // 0.5*(C3312+C3321)
  matrix_voigt(2, 4) =
      0.5 * (four_tensor(2, 2, 1, 2) + four_tensor(2, 2, 2, 1));  // 0.5*(C3323+C3332)
  matrix_voigt(2, 5) =
      0.5 * (four_tensor(2, 2, 0, 2) + four_tensor(2, 2, 2, 0));  // 0.5*(C3313+C3331)

  matrix_voigt(3, 0) =
      0.5 * (four_tensor(0, 1, 0, 0) + four_tensor(1, 0, 0, 0));  // 0.5*(C1211+C2111)
  matrix_voigt(3, 1) =
      0.5 * (four_tensor(0, 1, 1, 1) + four_tensor(1, 0, 1, 1));  // 0.5*(C1222+C2122)
  matrix_voigt(3, 2) =
      0.5 * (four_tensor(0, 1, 2, 2) + four_tensor(1, 0, 2, 2));  // 0.5*(C1233+C2133)
  matrix_voigt(3, 3) =
      0.25 * (four_tensor(0, 1, 0, 1) + four_tensor(1, 0, 0, 1) + four_tensor(0, 1, 1, 0) +
                 four_tensor(1, 0, 1, 0));  // 0.5*(C1212+C2112+C1221+C2121)
  matrix_voigt(3, 4) =
      0.25 * (four_tensor(0, 1, 1, 2) + four_tensor(1, 0, 1, 2) + four_tensor(0, 1, 2, 1) +
                 four_tensor(1, 0, 2, 1));  // 0.5*(C1223+C2123+C1232+C2132)
  matrix_voigt(3, 5) =
      0.25 * (four_tensor(0, 1, 0, 2) + four_tensor(1, 0, 0, 2) + four_tensor(0, 1, 2, 0) +
                 four_tensor(1, 0, 2, 0));  // 0.5*(C1213+C2113+C1231+C2131)

  matrix_voigt(4, 0) =
      0.5 * (four_tensor(1, 2, 0, 0) + four_tensor(2, 1, 0, 0));  // 0.5*(C2311+C3211)
  matrix_voigt(4, 1) =
      0.5 * (four_tensor(1, 2, 1, 1) + four_tensor(2, 1, 1, 1));  // 0.5*(C2322+C3222)
  matrix_voigt(4, 2) =
      0.5 * (four_tensor(1, 2, 2, 2) + four_tensor(2, 1, 2, 2));  // 0.5*(C2333+C3233)
  matrix_voigt(4, 3) =
      0.25 * (four_tensor(1, 2, 0, 1) + four_tensor(2, 1, 0, 1) + four_tensor(1, 2, 1, 0) +
                 four_tensor(2, 1, 1, 0));  // 0.5*(C2312+C3212+C2321+C3221)
  matrix_voigt(4, 4) =
      0.25 * (four_tensor(1, 2, 1, 2) + four_tensor(2, 1, 1, 2) + four_tensor(1, 2, 2, 1) +
                 four_tensor(2, 1, 2, 1));  // 0.5*(C2323+C3223+C2332+C3232)
  matrix_voigt(4, 5) =
      0.25 * (four_tensor(1, 2, 0, 2) + four_tensor(2, 1, 0, 2) + four_tensor(1, 2, 2, 0) +
                 four_tensor(2, 1, 2, 0));  // 0.5*(C2313+C3213+C2331+C3231)

  matrix_voigt(5, 0) =
      0.5 * (four_tensor(0, 2, 0, 0) + four_tensor(2, 0, 0, 0));  // 0.5*(C1311+C3111)
  matrix_voigt(5, 1) =
      0.5 * (four_tensor(0, 2, 1, 1) + four_tensor(2, 0, 1, 1));  // 0.5*(C1322+C3122)
  matrix_voigt(5, 2) =
      0.5 * (four_tensor(0, 2, 2, 2) + four_tensor(2, 0, 2, 2));  // 0.5*(C1333+C3133)
  matrix_voigt(5, 3) =
      0.25 * (four_tensor(0, 2, 0, 1) + four_tensor(2, 0, 0, 1) + four_tensor(0, 2, 1, 0) +
                 four_tensor(2, 0, 1, 0));  // 0.5*(C1312+C3112+C1321+C3121)
  matrix_voigt(5, 4) =
      0.25 * (four_tensor(0, 2, 1, 2) + four_tensor(2, 0, 1, 2) + four_tensor(0, 2, 2, 1) +
                 four_tensor(2, 0, 2, 1));  // 0.5*(C1323+C3123+C1332+C3132)
  matrix_voigt(5, 5) =
      0.25 * (four_tensor(0, 2, 0, 2) + four_tensor(2, 0, 0, 2) + four_tensor(0, 2, 2, 0) +
                 four_tensor(2, 0, 2, 0));  // 0.5*(C1313+C3113+C1331+C3131)
}

void Mat::add_dyadic_product_matrix_matrix(Core::LinAlg::FourTensor<3>& four_tensor_result,
    const Core::LinAlg::Matrix<3, 3>& matrix_A, const Core::LinAlg::Matrix<3, 3>& matrix_B)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          four_tensor_result(i, j, k, l) += matrix_A(i, j) * matrix_B(k, l);
}

void Mat::add_dyadic_product_matrix_matrix(Core::LinAlg::FourTensor<3>& four_tensor_result,
    const double scale, const Core::LinAlg::Matrix<3, 3>& matrix_A,
    const Core::LinAlg::Matrix<3, 3>& matrix_B)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          four_tensor_result(i, j, k, l) += scale * matrix_A(i, j) * matrix_B(k, l);
}

void Mat::add_contraction_matrix_four_tensor(Core::LinAlg::Matrix<3, 3>& matrix_result,
    const Core::LinAlg::Matrix<3, 3>& matrix, const Core::LinAlg::FourTensor<3>& four_tensor)
{
  for (unsigned k = 0; k < 3; ++k)
    for (unsigned l = 0; l < 3; ++l)
      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
          matrix_result(k, l) += matrix(i, j) * four_tensor(i, j, k, l);
}

void Mat::add_contraction_matrix_four_tensor(Core::LinAlg::Matrix<3, 3>& matrix_result,
    const double scale, const Core::LinAlg::FourTensor<3>& four_tensor,
    const Core::LinAlg::Matrix<3, 3>& matrix)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          matrix_result(i, j) += scale * four_tensor(i, j, k, l) * matrix(k, l);
}

void Mat::calculate_linear_isotropic_elastic_tensor(Core::LinAlg::FourTensor<3>& elasticity_tensor,
    const double youngs_modulus, const double poisson_ratio)
{
  const double lambda =
      poisson_ratio * youngs_modulus / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
  const double mu = youngs_modulus / (2 * (1 + poisson_ratio));

  const auto eye = [](int i, int j) { return i == j ? 1.0 : 0.0; };

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
          elasticity_tensor(i, j, k, l) =
              lambda * eye(i, j) * eye(k, l) + mu * (eye(i, k) * eye(j, l) + eye(i, l) * eye(j, k));
}

void Mat::calculate_deviatoric_projection_tensor(
    Core::LinAlg::FourTensor<3>& four_tensor, const double scale)
{
  const auto eye = [](int i, int j) { return i == j ? 1.0 : 0.0; };

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int k = 0; k < 3; ++k)
      {
        for (unsigned int l = 0; l < 3; ++l)
          four_tensor(i, j, k, l) =
              scale * (0.5 * eye(i, k) * eye(j, l) + 0.5 * eye(i, l) * eye(j, k) -
                          1.0 / 3 * eye(i, j) * eye(k, l));
      }
    }
  }
}

double Mat::contract_matrix_matrix(
    const Core::LinAlg::Matrix<3, 3>& matrix_A, const Core::LinAlg::Matrix<3, 3>& matrix_B)
{
  double scalarContraction = 0.0;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) scalarContraction += matrix_A(i, j) * matrix_B(i, j);

  return scalarContraction;
}

void Mat::four_tensor_to_matrix(const Core::LinAlg::FourTensor<3>& T, Core::LinAlg::Matrix<6, 6>& A)
{
  A(0, 0) = T(0, 0, 0, 0);  // xx-xx
  A(0, 1) = T(0, 0, 1, 1);  // xx-yy
  A(0, 2) = T(0, 0, 2, 2);  // xx-zz
  A(0, 3) = T(0, 0, 0, 1);  // xx-xy
  A(0, 4) = T(0, 0, 1, 2);  // xx-yz
  A(0, 5) = T(0, 0, 0, 2);  // xx-xz

  A(1, 0) = T(1, 1, 0, 0);
  A(1, 1) = T(1, 1, 1, 1);
  A(1, 2) = T(1, 1, 2, 2);
  A(1, 3) = T(1, 1, 0, 1);
  A(1, 4) = T(1, 1, 1, 2);
  A(1, 5) = T(1, 1, 0, 2);

  A(2, 0) = T(2, 2, 0, 0);
  A(2, 1) = T(2, 2, 1, 1);
  A(2, 2) = T(2, 2, 2, 2);
  A(2, 3) = T(2, 2, 0, 1);
  A(2, 4) = T(2, 2, 1, 2);
  A(2, 5) = T(2, 2, 0, 2);

  A(3, 0) = T(0, 1, 0, 0);
  A(3, 1) = T(0, 1, 1, 1);
  A(3, 2) = T(0, 1, 2, 2);
  A(3, 3) = T(0, 1, 0, 1);
  A(3, 4) = T(0, 1, 1, 2);
  A(3, 5) = T(0, 1, 0, 2);

  A(4, 0) = T(1, 2, 0, 0);
  A(4, 1) = T(1, 2, 1, 1);
  A(4, 2) = T(1, 2, 2, 2);
  A(4, 3) = T(1, 2, 0, 1);
  A(4, 4) = T(1, 2, 1, 2);
  A(4, 5) = T(1, 2, 0, 2);

  A(5, 0) = T(0, 2, 0, 0);
  A(5, 1) = T(0, 2, 1, 1);
  A(5, 2) = T(0, 2, 2, 2);
  A(5, 3) = T(0, 2, 0, 1);
  A(5, 4) = T(0, 2, 1, 2);
  A(5, 5) = T(0, 2, 0, 2);
}

// explicit instantiation of template functions
template void Mat::add_right_non_symmetric_holzapfel_product<double>(
    Core::LinAlg::Matrix<6, 9, double>&, Core::LinAlg::Matrix<3, 3, double> const&,
    Core::LinAlg::Matrix<3, 3, double> const&, double const);
template void Mat::add_right_non_symmetric_holzapfel_product<FAD>(Core::LinAlg::Matrix<6, 9, FAD>&,
    Core::LinAlg::Matrix<3, 3, FAD> const&, Core::LinAlg::Matrix<3, 3, FAD> const&, FAD const);
template void Mat::add_right_non_symmetric_holzapfel_product_strain_like<double>(
    Core::LinAlg::Matrix<6, 9, double>& out, Core::LinAlg::Matrix<3, 3, double> const& A,
    Core::LinAlg::Matrix<3, 3, double> const& B, double const fac);
template void Mat::add_right_non_symmetric_holzapfel_product_strain_like<FAD>(
    Core::LinAlg::Matrix<6, 9, FAD>& out, Core::LinAlg::Matrix<3, 3, FAD> const& A,
    Core::LinAlg::Matrix<3, 3, FAD> const& B, FAD const fac);
template void Mat::add_holzapfel_product<double>(Core::LinAlg::Matrix<6, 6, double>&,
    const Core::LinAlg::Matrix<6, 1, double>&, const double scalar);
template void Mat::add_holzapfel_product<FAD>(
    Core::LinAlg::Matrix<6, 6, FAD>&, const Core::LinAlg::Matrix<6, 1, FAD>&, const FAD scalar);

template void Mat::clear_four_tensor<3>(Core::LinAlg::FourTensor<3>& four_tensor);

template void Mat::multiply_four_tensor_matrix<3>(Core::LinAlg::FourTensor<3>& four_tensor_result,
    const Core::LinAlg::FourTensor<3>& four_tensor, const Core::LinAlg::Matrix<3, 3>& matrix,
    const bool clear_result_tensor);

template void Mat::multiply_matrix_four_tensor<3>(Core::LinAlg::FourTensor<3>& four_tensor_result,
    const Core::LinAlg::Matrix<3, 3>& matrix, const Core::LinAlg::FourTensor<3>& four_tensor,
    const bool clear_result_tensor);

template void Mat::multiply_matrix_four_tensor_by_second_index<3>(
    Core::LinAlg::FourTensor<3>& four_tensor_result, const Core::LinAlg::Matrix<3, 3>& matrix,
    const Core::LinAlg::FourTensor<3>& four_tensor, const bool clear_result_tensor);

template void Mat::multiply_four_tensor_four_tensor<3>(
    Core::LinAlg::FourTensor<3>& four_tensor_result,
    const Core::LinAlg::FourTensor<3>& four_tensor_1,
    const Core::LinAlg::FourTensor<3>& four_tensor_2, const bool clear_result_tensor);

template Core::LinAlg::Matrix<6, 6> Mat::pull_back_four_tensor<3>(
    const Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<6, 6>& cmat_voigt);

template double Mat::get_pull_back_four_tensor_entry<3>(const Core::LinAlg::Matrix<3, 3>& defgrd,
    const Core::LinAlg::FourTensor<3>& four_tensor, const int i, const int j, const int k,
    const int l);

template void Mat::setup_four_tensor_from_6x6_voigt_matrix<3>(
    Core::LinAlg::FourTensor<3>& four_tensor, const Core::LinAlg::Matrix<6, 6>& matrix_voigt);

template void Mat::setup_6x6_voigt_matrix_from_four_tensor<3>(
    Core::LinAlg::Matrix<6, 6>& matrix_voigt, const Core::LinAlg::FourTensor<3>& four_tensor);

FOUR_C_NAMESPACE_CLOSE
