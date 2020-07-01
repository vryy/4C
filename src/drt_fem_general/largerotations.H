/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of utility functions for large rotations

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef LARGEROTATIONS_DEFINES_H
#define LARGEROTATIONS_DEFINES_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../headers/FAD_utils.H"

/************************************************************************/
/* Large rotations functions                                            */
/************************************************************************/
namespace LARGEROTATIONS
{
  //! computing spin matrix S from rotation angle vector theta, Crisfield, Vol. 2, equation (16.8)
  // Fixme: remove this function and all non-templates in this file; should be obsolete!
  void computespin(LINALG::Matrix<3, 3>& S, const LINALG::Matrix<3, 1>& theta);

  //! this is a templated version of the method computespin as declared above in order to also
  //! handle the datatype FAD which is needed for automatic differentiation
  template <typename T>
  void computespin(LINALG::Matrix<3, 3, T>& S, const LINALG::Matrix<3, 1, T>& theta)
  {
    // function based on Crisfield Vol. 2, Section 16 (16.8)
    S(0, 0) = 0.0;
    S(0, 1) = -theta(2);
    S(0, 2) = theta(1);
    S(1, 0) = theta(2);
    S(1, 1) = 0.0;
    S(1, 2) = -theta(0);
    S(2, 0) = -theta(1);
    S(2, 1) = theta(0);
    S(2, 2) = 0.0;
  }

  //! computing rotation matrix R from quaternion q, Crisfield, Vol. 2, equation (16.70)
  void quaterniontotriad(const LINALG::Matrix<4, 1>& q, LINALG::Matrix<3, 3>& R);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <typename T>
  void quaterniontotriad(const LINALG::Matrix<4, 1, T>& q, LINALG::Matrix<3, 3, T>& R)
  {
    // separate storage of vector part of q
    LINALG::Matrix<3, 1, T> qvec;
    for (int i = 0; i < 3; i++) qvec(i) = q(i);

    // setting R to third summand of equation (16.70)
    LARGEROTATIONS::computespin(R, qvec);
    R.Scale(2 * q(3));

    // adding second summand of equation (16.70)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) R(i, j) += 2 * q(i) * q(j);

    // adding diagonal entries according to first summand of equation (16.70)
    R(0, 0) = 1 - 2 * (q(1) * q(1) + q(2) * q(2));
    R(1, 1) = 1 - 2 * (q(0) * q(0) + q(2) * q(2));
    R(2, 2) = 1 - 2 * (q(0) * q(0) + q(1) * q(1));

    return;
  }  // quaterniontotriad

  //! computing quaternion q from rotation angle theta, Crisfield, Vol. 2, equation (16.67)
  void angletoquaternion(const LINALG::Matrix<3, 1>& theta, LINALG::Matrix<4, 1>& q);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <typename T>
  void angletoquaternion(const LINALG::Matrix<3, 1, T>& theta, LINALG::Matrix<4, 1, T>& q)
  {
    // absolute value of rotation angle theta
    T abs_theta = 0.0;
    for (int i = 0; i < 3; i++) abs_theta += theta(i) * theta(i);

    abs_theta = pow(abs_theta, 0.5);

    // computing quaterion for rotation by angle theta, Crisfield, Vol. 2, equation (16.67)
    if (abs_theta > 1e-12)
    {
      q(0) = theta(0) * sin(abs_theta / 2) / abs_theta;
      q(1) = theta(1) * sin(abs_theta / 2) / abs_theta;
      q(2) = theta(2) * sin(abs_theta / 2) / abs_theta;
      q(3) = cos(abs_theta / 2);
    }
    else  // altered due to FAD - with theta->0 the equations above simplify in the manner below
          // NOTE: the equations below are chosen in a way that the linearisation (with FAD) is the
          // same of the equations above for the case theta->0
    {
      q(0) = 0.5 * theta(0);
      q(1) = 0.5 * theta(1);
      q(2) = 0.5 * theta(2);
      q(3) = cos(abs_theta / 2);
    }

    return;
  }  // angletoquaternion

  //! computing rotation matrix R from angle theta, Crisfield Vol. 2, equation (16.22)
  void angletotriad(const LINALG::Matrix<3, 1>& theta, LINALG::Matrix<3, 3>& R);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <typename T>
  void angletotriad(const LINALG::Matrix<3, 1, T>& theta, const LINALG::Matrix<3, 3, T>& triad_ref,
      LINALG::Matrix<3, 3, T>& triad)
  {
    LINALG::Matrix<3, 3, T> rotation_matrix;
    LINALG::Matrix<4, 1, T> quaternion;

    rotation_matrix.Clear();
    quaternion.Clear();
    triad.Clear();

    LARGEROTATIONS::angletoquaternion(theta, quaternion);
    LARGEROTATIONS::quaterniontotriad(quaternion, rotation_matrix);
    triad.Multiply(triad_ref, rotation_matrix);  // ???
  }

  //! computing rotation angle vector theta from quaternion q, Crisfield, Vol. 2, equation (16.67)
  void quaterniontoangle(const LINALG::Matrix<4, 1>& q, LINALG::Matrix<3, 1>& theta);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <typename T>
  void quaterniontoangle(const LINALG::Matrix<4, 1, T>& q, LINALG::Matrix<3, 1, T>& theta)
  {
    /*the following function computes from a quaternion q an angle theta within [-PI; PI]; such an
     * interval is imperative for the use of the resulting angle together with formulae like
     * Crisfield, Vol. 2, equation (16.90); note that these formulae comprise not only trigonometric
     * functions, but rather the angle theta directly. Hence they are not 2*PI-invariant !!! */

    // if the rotation angle is pi we have q(3) == 0 and the rotation angle vector can be computed
    // by
    if (q(3) == 0)
    {
      // note that with q(3) == 0 the first three elements of q represent the unit direction vector
      // of the angle according to Crisfield, Vol. 2, equation (16.67)
      for (int i = 0; i < 3; i++) theta(i) = q(i) * M_PI;
    }
    else
    {
      // otherwise the angle can be computed from a quaternion via Crisfield, Vol. 2, eq. (16.79)
      LINALG::Matrix<3, 1, T> omega;
      T abs_omega = 0.0;
      for (int i = 0; i < 3; i++)
      {
        omega(i) = q(i) * 2 / q(3);
        abs_omega += omega(i) * omega(i);
      }

      abs_omega = pow(abs_omega, 0.5);
      T tanhalf = abs_omega / 2;
      T thetaabs = atan(tanhalf) * 2;

      // if the rotation angle is zero we return a zero rotation angle vector at once
      if (abs_omega < 1e-12)
      {
        for (int i = 0; i < 3; i++)
          theta(i) =
              2 * q(i);  // altered due to FAD, the equations simplify in that way for theta->0
      }
      else
        for (int i = 0; i < 3; i++) theta(i) = thetaabs * omega(i) / abs_omega;
    }

    return;
  }

  //! computing from quaternion q the Rodrigues parameters omega, Crisfield, Vol. 2, equation
  //! (16.79)
  void quaterniontorodrigues(const LINALG::Matrix<4, 1>& q, LINALG::Matrix<3, 1>& omega);

  //! computing quaternion q from a rotation matrix R, Crisfield, Vol. 2, equation (16.74) - (16.78)
  void triadtoquaternion(const LINALG::Matrix<3, 3>& R, LINALG::Matrix<4, 1>& q);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <typename T>
  void triadtoquaternion(const LINALG::Matrix<3, 3, T>& R, LINALG::Matrix<4, 1, T>& q)
  {
    T trace = R(0, 0) + R(1, 1) + R(2, 2);

    if (trace > R(0, 0) && trace > R(1, 1) && trace > R(2, 2))
    {
      q(3) = 0.5 * pow(1 + trace, 0.5);
      /*note: if trace is greater than each element on diagonal, all diagonal elements are positive
       *and hence also the trace is positive and thus q(3) > 0 so that division by q(3) is allowed*/
      q(0) = (R(2, 1) - R(1, 2)) / (4 * q(3));
      q(1) = (R(0, 2) - R(2, 0)) / (4 * q(3));
      q(2) = (R(1, 0) - R(0, 1)) / (4 * q(3));
    }
    else
    {
      for (int i = 0; i < 3; i++)
      {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;

        if (R(i, i) >= R(j, j) && R(i, i) >= R(k, k))
        {
          // equation (16.78a)
          q(i) = pow(0.5 * R(i, i) + 0.25 * (1 - trace), 0.5);

          // equation (16.78b)
          q(3) = 0.25 * (R(k, j) - R(j, k)) / q(i);

          // equation (16.78c)
          q(j) = 0.25 * (R(j, i) + R(i, j)) / q(i);
          q(k) = 0.25 * (R(k, i) + R(i, k)) / q(i);
        }
      }
    }

    return;
  }

  //! computing the rotation angle theta rotating the unit direction vector d1 into the unit
  //! direction vector d2
  void directionstoangle(
      const LINALG::Matrix<3, 1>& d1, const LINALG::Matrix<3, 1>& d2, LINALG::Matrix<3, 1>& theta);

  //! computing the triad R rotating the unit direction vector d1 into the unit direction vector d2
  void directionstotriad(
      const LINALG::Matrix<3, 1>& d1, const LINALG::Matrix<3, 1>& d2, LINALG::Matrix<3, 3>& R);

  //! quaternion product q12 = q2*q1, Crisfield, Vol. 2, equation (16.71)
  //  void quaternionproduct(const LINALG::Matrix<4,1>& q1,const LINALG::Matrix<4,1>&
  //  q2,LINALG::Matrix<4,1>& q12);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <class T1, class T2, class T3>
  void quaternionproduct(const T1& q1, const T2& q2, T3& q12)
  {
    dsassert(q12.M() == 4 and q12.N() == 1 and q1.M() == 4 and q1.N() == 1 and q2.M() == 4 and
                 q2.N() == 1,
        "size mismatch: expected 4x1 vector for quaternion");

    q12(0) = q2(3) * q1(0) + q1(3) * q2(0) + q2(1) * q1(2) - q1(1) * q2(2);
    q12(1) = q2(3) * q1(1) + q1(3) * q2(1) + q2(2) * q1(0) - q1(2) * q2(0);
    q12(2) = q2(3) * q1(2) + q1(3) * q2(2) + q2(0) * q1(1) - q1(0) * q2(1);
    q12(3) = q2(3) * q1(3) - q2(2) * q1(2) - q2(1) * q1(1) - q2(0) * q1(0);
  }

  //! compute inverse quaternion q^{-1} of quaternion q
  LINALG::Matrix<4, 1> inversequaternion(const LINALG::Matrix<4, 1>& q);

  //! this is a templated version of the method declared above in order to also handle
  //! the datatype FAD which is needed for automatic differentiation
  template <typename T>
  LINALG::Matrix<4, 1, T> inversequaternion(const LINALG::Matrix<4, 1, T>& q)
  {
    // square norm ||q||^2 of quaternion q
    T qnormsq = 0.0;

    for (int i = 0; i < 4; i++)
    {
      qnormsq += q(i) * q(i);
    }

    // declaration of variable for inverse quaternion
    LINALG::Matrix<4, 1, T> qinv;

    // inverse quaternion q^(-1) = [-q0, -q1, -q2, q3] / ||q||^2;
    for (int i = 0; i < 3; i++) qinv(i) = -q(i) / qnormsq;

    qinv(3) = q(3) / qnormsq;

    return qinv;
  }

  //! matrix T(\theta) from Jelenic 1999, eq. (2.5), equivalent to matrix H^{-1} in Crisfield,
  //! (16.94)
  LINALG::Matrix<3, 3> Tmatrix(LINALG::Matrix<3, 1> theta);

  //! matrix T(\theta) from Jelenic 1999, eq. (2.5), equivalent to matrix H^{-1} in Crisfield,
  //! (16.94)
  template <typename T>
  LINALG::Matrix<3, 3, T> Tmatrix(LINALG::Matrix<3, 1, T> theta)
  {
    LINALG::Matrix<3, 3, T> result(true);
    T theta_abs = pow(theta(0) * theta(0) + theta(1) * theta(1) + theta(2) * theta(2), 0.5);

    // in case of theta_abs == 0 the following computation has problems with singularities
    if (theta_abs > 1e-8)
    {
      computespin(result, theta);
      result.Scale(-0.5);

      T theta_abs_half = theta_abs / 2.0;

      for (int i = 0; i < 3; i++) result(i, i) += theta_abs / (2.0 * tan(theta_abs_half));

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          result(i, j) += theta(i) * theta(j) * (1.0 - theta_abs / (2.0 * tan(theta_abs_half))) /
                          (theta_abs * theta_abs);
    }
    // based on the small angle approximation tan(x)=x, we get: T = I - 0.5*S(theta)
    else
    {
      LARGEROTATIONS::computespin(result, theta);
      result.Scale(-0.5);
      for (int j = 0; j < 3; j++) result(j, j) += 1.0;
    }
    return result;
  }  // LARGEROTATIONS::Tmatrix

  //! matrix T(\theta)^{-1} from Jelenic 1999, eq. (2.5)
  LINALG::Matrix<3, 3> Tinvmatrix(LINALG::Matrix<3, 1> theta);

  //! this is a templated version of the method Tinvmatrix as declared in LARGEROTATIONS in order to
  //! also handle the datatype FAD which is needed for automatic differentiation
  template <typename T>
  LINALG::Matrix<3, 3, T> Tinvmatrix(LINALG::Matrix<3, 1, T> theta)
  {
    LINALG::Matrix<3, 3, T> result(true);
    T theta_abs = pow(theta(0) * theta(0) + theta(1) * theta(1) + theta(2) * theta(2), 0.5);

    // in case of theta_abs == 0 the following computation has problems with ill-conditioning /
    // singularities
    if (theta_abs > 1e-8)
    {
      // ultimate term in eq. (2.5)
      LARGEROTATIONS::computespin(result, theta);
      result.Scale((1.0 - cos(theta_abs)) / (theta_abs * theta_abs));

      // penultimate term in eq. (2.5)
      for (int i = 0; i < 3; i++) result(i, i) += sin(theta_abs) / (theta_abs);

      // first term on the right side in eq. (2.5)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          result(i, j) +=
              theta(i) * theta(j) * (1.0 - sin(theta_abs) / (theta_abs)) / (theta_abs * theta_abs);
    }
    // in case of theta_abs == 0 H(theta) is the identity matrix and hence also Hinv
    else
    {
      // based on the small angle approximation sin(x)=x and 1-cos(x)=x^2/2 we get: Tinv = I +
      // 0.5*S(theta) -> Tinv'=0.5*S(theta')
      LARGEROTATIONS::computespin(result, theta);
      result.Scale(0.5);
      for (int j = 0; j < 3; j++) result(j, j) += 1.0;
    }

    return result;
  }  // Tinvmatrix

  //! compute dT^{-1}(\theta)/dx according to the two-lined equation below (3.19) on page 152 of
  //! Jelenic 1999
  void computedTinvdx(const LINALG::Matrix<3, 1>& Psil, const LINALG::Matrix<3, 1>& Psilprime,
      LINALG::Matrix<3, 3>& dTinvdx);

  //! compute dT^{-1}(\theta)/dx according to the two-lined equation below (3.19) on page 152 of
  //! Jelenic 1999
  template <typename T>
  void computedTinvdx(const LINALG::Matrix<3, 1, T>& Psil, const LINALG::Matrix<3, 1, T>& Psilprime,
      LINALG::Matrix<3, 3, T>& dTinvdx)
  {
    // auxiliary matrix for storing intermediate results
    LINALG::Matrix<3, 3, T> auxmatrix;

    // norm of \Psi^l:
    T normPsil = pow(Psil(0) * Psil(0) + Psil(1) * Psil(1) + Psil(2) * Psil(2), 0.5);

    // for relative rotations smaller then 1e-12 we use the limit for Psil -> 0 according to the
    // comment above NOTE 4 on page 152, Jelenic 1999
    if (normPsil < 1e-8)
    {
      computespin(dTinvdx, Psilprime);
      dTinvdx.Scale(0.5);
    }
    else
    {
      // scalarproduct \Psi^{l,t} \cdot \Psi^{l,'}
      T scalarproductPsilPsilprime = 0.0;
      for (int i = 0; i < 3; i++) scalarproductPsilPsilprime += Psil(i) * Psilprime(i);

      // spin matrices of Psil and Psilprime
      LINALG::Matrix<3, 3, T> spinPsil;
      LINALG::Matrix<3, 3, T> spinPsilprime;
      computespin(spinPsil, Psil);
      computespin(spinPsilprime, Psilprime);

      // third summand
      dTinvdx.Multiply(spinPsilprime, spinPsil);
      auxmatrix.Multiply(spinPsil, spinPsilprime);
      dTinvdx += auxmatrix;
      dTinvdx.Scale((1.0 - sin(normPsil) / normPsil) / (normPsil * normPsil));

      // first summand
      auxmatrix.PutScalar(0.0);
      auxmatrix += spinPsil;
      auxmatrix.Scale(scalarproductPsilPsilprime *
                      (normPsil * sin(normPsil) - 2 * (1.0 - cos(normPsil))) /
                      (normPsil * normPsil * normPsil * normPsil));
      dTinvdx += auxmatrix;

      // second summand
      auxmatrix.PutScalar(0.0);
      auxmatrix += spinPsilprime;
      auxmatrix.Scale((1.0 - cos(normPsil)) / (normPsil * normPsil));
      dTinvdx += auxmatrix;

      // fourth summand
      auxmatrix.Multiply(spinPsil, spinPsil);
      auxmatrix.Scale(scalarproductPsilPsilprime *
                      (3.0 * sin(normPsil) - normPsil * (2.0 + cos(normPsil))) /
                      (normPsil * normPsil * normPsil * normPsil * normPsil));
      dTinvdx += auxmatrix;
    }

    return;
  }  // LARGEROTATIONS::computedTinvdx

  //! Transformation from node number according to Crisfield 1999 to storage position applied in
  //! BACI
  unsigned int NumberingTrafo(const unsigned int j, const unsigned int numnode);

  //! Rotate an arbitrary triad around its first base vector (tangent)
  template <typename T>
  void RotateTriad(
      const LINALG::Matrix<3, 3, T>& triad, const T& alpha, LINALG::Matrix<3, 3, T>& triad_rot)
  {
    for (int i = 0; i < 3; i++)
    {
      triad_rot(i, 0) = triad(i, 0);
      triad_rot(i, 1) = triad(i, 1) * cos(alpha) + triad(i, 2) * sin(alpha);
      triad_rot(i, 2) = triad(i, 2) * cos(alpha) - triad(i, 1) * sin(alpha);
    }

    return;
  }

  //! Calculate the SR mapping for a given reference system triad_ref and a given tangent vector r_s
  template <typename T>
  void CalculateSRTriads(const LINALG::Matrix<3, 1, T>& r_s,
      const LINALG::Matrix<3, 3, T>& triad_ref, LINALG::Matrix<3, 3, T>& triad)
  {
    // In this calculation, r_s does not necessarily have to be a unit vector.
    T temp_scalar1 = 0.0;
    T temp_scalar2 = 0.0;
    T temp_scalar3 = 0.0;
    T fac_n0 = 0.0;
    T fac_b0 = 0.0;
    T abs_r_s = FADUTILS::Norm<T>(r_s);

    for (int i = 0; i < 3; i++)
    {
      temp_scalar1 += triad_ref(i, 1) * r_s(i);
      temp_scalar2 += triad_ref(i, 2) * r_s(i);
      temp_scalar3 += triad_ref(i, 0) * r_s(i);
    }

    // This line is necessary in order to avoid a floating point exception and subsequent abortion
    // of the simulation in case the denominator equals zero. In such a case, convergence is still
    // unlikely. However, since the simulation can continue, functionalities such as
    // DIVERCONT==adapt_step will still work!
    if (std::fabs(abs_r_s + temp_scalar3) < 1.0e-10)
    {
      // Todo @grill: think about using runtime error here
      std::cout << "Warning: SR-Denominator scaled!" << std::endl;
      temp_scalar3 = 0.99 * temp_scalar3;
    }

    fac_n0 = temp_scalar1 / (abs_r_s + temp_scalar3);
    fac_b0 = temp_scalar2 / (abs_r_s + temp_scalar3);

    for (int i = 0; i < 3; i++)
    {
      triad(i, 0) = r_s(i) / abs_r_s;
      triad(i, 1) = triad_ref(i, 1) - fac_n0 * (r_s(i) / abs_r_s + triad_ref(i, 0));
      triad(i, 2) = triad_ref(i, 2) - fac_b0 * (r_s(i) / abs_r_s + triad_ref(i, 0));
    }

    return;
  }

  //! compute the relative angle theta between triad_ref and triad so that exp(theta) =
  //! triad_ref^T*triad -> inversion of the right translation triad=triad_ref*exp(theta)
  template <typename T>
  void triadtoangleright(LINALG::Matrix<3, 1, T>& theta, const LINALG::Matrix<3, 3, T>& triad_ref,
      const LINALG::Matrix<3, 3, T>& triad)
  {
    LINALG::Matrix<3, 3, T> rotation_matrix;
    LINALG::Matrix<4, 1, T> quaternion;

    rotation_matrix.Clear();
    quaternion.Clear();

    rotation_matrix.MultiplyTN(triad_ref, triad);
    LARGEROTATIONS::triadtoquaternion(rotation_matrix, quaternion);
    LARGEROTATIONS::quaterniontoangle(quaternion, theta);
  }

  //! compute the relative angle theta between triad_ref and triad so that exp(theta) =
  //! triad*triad_ref^T -> inversion of the right translation triad=exp(theta)*triad_ref
  template <typename T>
  void triadtoangleleft(LINALG::Matrix<3, 1, T>& theta, const LINALG::Matrix<3, 3, T>& triad_ref,
      const LINALG::Matrix<3, 3, T>& triad)
  {
    LINALG::Matrix<3, 3, T> rotation_matrix;
    LINALG::Matrix<4, 1, T> quaternion;

    rotation_matrix.Clear();
    quaternion.Clear();

    rotation_matrix.MultiplyNT(triad, triad_ref);
    LARGEROTATIONS::triadtoquaternion(rotation_matrix, quaternion);
    LARGEROTATIONS::quaterniontoangle(quaternion, theta);
  }


}  // namespace LARGEROTATIONS

#endif  // #ifndef LARGEROTATIONS_DEFINES_H
