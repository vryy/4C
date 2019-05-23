/*----------------------------------------------------------------------------*/
/*!

\brief A set of utility functions for large rotations

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "largerotations.H"
#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 |computes from a quaternion q rodrigues parameters omega (public)cyron02/09|
 *---------------------------------------------------------------------------*/
void LARGEROTATIONS::quaterniontorodrigues(
    const LINALG::Matrix<4, 1>& q, LINALG::Matrix<3, 1>& omega)
{
  /*the Rodrigues parameters are defined only for angles whose absolute valued is smaller than PI,
   * i.e. for which the fourth component of the quaternion is unequal zero; if this is not satisfied
   * for the quaternion passed into this method an error is thrown*/
  if (q(3) == 0)
    dserror("cannot compute Rodrigues parameters for angles with absolute valued PI !!!");

  // in any case except for the one dealt with above the angle can be computed from a quaternion via
  // Crisfield, Vol. 2, eq. (16.79)
  for (int i = 0; i < 3; i++) omega(i) = q(i) * 2 / q(3);

  return;
}  // LARGEROTATIONS::quaterniontorodrigues

/*----------------------------------------------------------------------*
 |computes from a quaternion q the related angle theta (public)cyron10/08|
 *----------------------------------------------------------------------*/
void LARGEROTATIONS::quaterniontoangle(const LINALG::Matrix<4, 1>& q, LINALG::Matrix<3, 1>& theta)
{
  /*the following function computes from a quaternion q an angle theta within ]-PI; PI]; such an
   * interval is imperative for the use of the resulting angle together with formulae like
   * Crisfield, Vol. 2, equation (16.90); note that these formulae comprise not only trigonometric
   * functions, but rather the angle theta directly. Hence they are not 2*PI-invariant !!! */

  // if the rotation angle is pi we have q(3) == 0 and the rotation angle vector can be computed by
  if (q(3) == 0)
  {
    // note that with q(3) == 0 the first three elements of q represent the unit direction vector of
    // the angle according to Crisfield, Vol. 2, equation (16.67)
    for (int i = 0; i < 3; i++) theta(i) = q(i) * M_PI;
  }
  else
  {
    // otherwise the angle can be computed from a quaternion via Crisfield, Vol. 2, eq. (16.79)
    LINALG::Matrix<3, 1> omega;
    for (int i = 0; i < 3; i++) omega(i) = q(i) * 2 / q(3);

    double tanhalf = omega.Norm2() / 2;
    double thetaabs = atan(tanhalf) * 2;

    // if the rotation angle is zero we return a zero rotation angle vector at once
    if (omega.Norm2() == 0)
    {
      for (int i = 0; i < 3; i++) theta(i) = 0;
    }
    else
      for (int i = 0; i < 3; i++) theta(i) = thetaabs * omega(i) / omega.Norm2();
  }

  return;
}  // LARGEROTATIONS::quaterniontoangle()

/*---------------------------------------------------------------------------*
 |computes a spin matrix out of a rotation vector        (public)cyron02/09|
 *---------------------------------------------------------------------------*/
void LARGEROTATIONS::computespin(LINALG::Matrix<3, 3>& S, const LINALG::Matrix<3, 1>& theta)
{
  // function based on Crisfield Vol. 2, Section 16 (16.8)
  S(0, 0) = 0;
  S(0, 1) = -theta(2);
  S(0, 2) = theta(1);
  S(1, 0) = theta(2);
  S(1, 1) = 0;
  S(1, 2) = -theta(0);
  S(2, 0) = -theta(1);
  S(2, 1) = theta(0);
  S(2, 2) = 0;

  return;
}  // LARGEROTATIONS::computespin

/*----------------------------------------------------------------------*
 |computes a rotation matrix R from a quaternion q            |
 |cf. Crisfield, Vol. 2, equation (16.70)         (public)cyron10/08|
 *----------------------------------------------------------------------*/
void LARGEROTATIONS::quaterniontotriad(const LINALG::Matrix<4, 1>& q, LINALG::Matrix<3, 3>& R)
{
  // separate storage of vector part of q
  LINALG::Matrix<3, 1> qvec;
  for (int i = 0; i < 3; i++) qvec(i) = q(i);

  // setting R to third summand of equation (16.70)
  computespin(R, qvec);
  R.Scale(2 * q(3));

  // adding second summand of equation (16.70)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) R(i, j) += 2 * q(i) * q(j);

  // adding diagonal entries according to first summand of equation (16.70)
  R(0, 0) = 1 - 2 * (q(1) * q(1) + q(2) * q(2));
  R(1, 1) = 1 - 2 * (q(0) * q(0) + q(2) * q(2));
  R(2, 2) = 1 - 2 * (q(0) * q(0) + q(1) * q(1));

  return;
}  // LARGEROTATIONS::quaterniontotriad

/*---------------------------------------------------------------------------*
 |computes a quaternion from an angle vector           (public)cyron02/09|
 *---------------------------------------------------------------------------*/
void LARGEROTATIONS::angletoquaternion(const LINALG::Matrix<3, 1>& theta, LINALG::Matrix<4, 1>& q)
{
  // absolute value of rotation angle theta
  double abs_theta = theta.Norm2();

  // computing quaterion for rotation by angle theta, Crisfield, Vol. 2, equation (16.67)
  if (abs_theta > 0)
  {
    q(0) = theta(0) * sin(abs_theta / 2) / abs_theta;
    q(1) = theta(1) * sin(abs_theta / 2) / abs_theta;
    q(2) = theta(2) * sin(abs_theta / 2) / abs_theta;
    q(3) = cos(abs_theta / 2);
  }
  else
  {
    q.PutScalar(0.0);
    q(3) = 1;
  }

  return;
}  // LARGEROTATIONS::angletoquaternion

/*----------------------------------------------------------------------*
 | Compute rotation matrix R from angle theta                cyron 11/10|
 *----------------------------------------------------------------------*/
void LARGEROTATIONS::angletotriad(const LINALG::Matrix<3, 1>& theta, LINALG::Matrix<3, 3>& R)
{
  // compute spin matrix according to Crisfield Vol. 2, equation (16.8)
  LINALG::Matrix<3, 3> spin;
  computespin(spin, theta);

  // nompute norm of theta
  double theta_abs = theta.Norm2();

  // build an identity matrix
  LINALG::Matrix<3, 3> identity;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      if (i == j)
        identity(i, j) = 1.0;
      else
        identity(i, j) = 0.0;
    }

  // square of spin matrix
  LINALG::Matrix<3, 3> spin2;
  spin2.Multiply(spin, spin);

  // compute rotation matrix according to Crisfield Vol. 2, equation (16.22)
  if (theta_abs > 1.0e-12)
  {
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        R(i, j) = identity(i, j) + spin(i, j) * (sin(theta_abs)) / theta_abs +
                  (1 - (cos(theta_abs))) / (theta_abs * theta_abs) * spin2(i, j);
      }
    }
  }
  else
  {
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        R(i, j) = identity(i, j);
      }
    }
  }

  return;
}

/*---------------------------------------------------------------------------*
 |computes a quaternion q from a rotation matrix R; all operations are      |
 |performed according to Crisfield, Vol. 2, section 16.10 and the there      |
 |described Spurrier's algorithm                   (public)cyron02/09|
 *---------------------------------------------------------------------------*/
void LARGEROTATIONS::triadtoquaternion(const LINALG::Matrix<3, 3>& R, LINALG::Matrix<4, 1>& q)
{
  double trace = R(0, 0) + R(1, 1) + R(2, 2);

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
        q(i) = std::pow(0.5 * R(i, i) + 0.25 * (1 - trace), 0.5);

        // equation (16.78b)
        q(3) = 0.25 * (R(k, j) - R(j, k)) / q(i);

        // equation (16.78c)
        q(j) = 0.25 * (R(j, i) + R(i, j)) / q(i);
        q(k) = 0.25 * (R(k, i) + R(i, k)) / q(i);
      }
    }
  }
  return;
}  // LARGEROTATIONS::TriadToQuaternion

/*---------------------------------------------------------------------------*
 |matrix T(\theta) from Jelenic 1999, eq. (2.5), equivalent to matrix H^(-1) |
 |in Crisfield, Vol. 2, equation (16.93)                 (public) cyron 04/10|
 *---------------------------------------------------------------------------*/
LINALG::Matrix<3, 3> LARGEROTATIONS::Tmatrix(LINALG::Matrix<3, 1> theta)
{
  LINALG::Matrix<3, 3> result(true);
  double theta_abs = sqrt(theta(0) * theta(0) + theta(1) * theta(1) + theta(2) * theta(2));

  // in case of theta_abs == 0 the following computation has problems with singularities
  if (theta_abs > 1.0e-8)
  {
    computespin(result, theta);
    result.Scale(-0.5);

    double theta_abs_half = theta_abs / 2.0;

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

/*---------------------------------------------------------------------------*
 |matrix T(\theta)^{-1} from Jelenic 1999, eq. (2.5)     (public) cyron 04/10|
 *---------------------------------------------------------------------------*/
LINALG::Matrix<3, 3> LARGEROTATIONS::Tinvmatrix(LINALG::Matrix<3, 1> theta)
{
  LINALG::Matrix<3, 3> result;
  double rootarg = theta(0) * theta(0) + theta(1) * theta(1) + theta(2) * theta(2);
  double theta_abs = sqrt(rootarg);

  // in case of theta_abs == 0 the following computation has problems with ill-conditioning /
  // singularities
  if (theta_abs > 1.0e-8)
  {
    // ultimate term in eq. (2.5)
    computespin(result, theta);
    result.Scale((1 - cos(theta_abs)) / (theta_abs * theta_abs));

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
    result.PutScalar(0.0);
    LARGEROTATIONS::computespin(result, theta);
    result.Scale(0.5);
    for (int j = 0; j < 3; j++) result(j, j) += 1.0;
  }

  return result;
}  // LARGEROTATIONS::Tinvmatrix

/*----------------------------------------------------------------------------------------------------------------------*
 |compute d(T^{-1})/dx according to the two-lined equation below (3.19) on page 152 of Jelenic 1999
 | | cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
void LARGEROTATIONS::computedTinvdx(const LINALG::Matrix<3, 1>& Psil,
    const LINALG::Matrix<3, 1>& Psilprime, LINALG::Matrix<3, 3>& dTinvdx)
{
  // auxiliary matrix for storing intermediate results
  LINALG::Matrix<3, 3> auxmatrix;

  // norm of \Psi^l:
  double normPsil = Psil.Norm2();

  // for relative rotations smaller then 1e-12 we use the limit for Psil -> 0 according to the
  // comment above NOTE 4 on page 152, Jelenic 1999
  if (Psil.Norm2() < 1e-8)
  {
    computespin(dTinvdx, Psilprime);
    dTinvdx.Scale(0.5);
  }
  else
  {
    // scalarproduct \Psi^{l,t} \cdot \Psi^{l,'}
    double scalarproductPsilPsilprime = 0;
    for (int i = 0; i < 3; i++) scalarproductPsilPsilprime += Psil(i) * Psilprime(i);

    // spin matrices of Psil and Psilprime
    LINALG::Matrix<3, 3> spinPsil;
    LINALG::Matrix<3, 3> spinPsilprime;
    computespin(spinPsil, Psil);
    computespin(spinPsilprime, Psilprime);

    // third summand
    dTinvdx.Multiply(spinPsilprime, spinPsil);
    auxmatrix.Multiply(spinPsil, spinPsilprime);
    dTinvdx += auxmatrix;
    dTinvdx.Scale((1 - sin(normPsil) / normPsil) / (normPsil * normPsil));

    // first summand
    auxmatrix.PutScalar(0);
    auxmatrix += spinPsil;
    auxmatrix.Scale(scalarproductPsilPsilprime *
                    (normPsil * sin(normPsil) - 2 * (1 - cos(normPsil))) /
                    (normPsil * normPsil * normPsil * normPsil));
    dTinvdx += auxmatrix;

    // second summand
    auxmatrix.PutScalar(0);
    auxmatrix += spinPsilprime;
    auxmatrix.Scale((1 - cos(normPsil)) / (normPsil * normPsil));
    dTinvdx += auxmatrix;

    // fourth summand
    auxmatrix.Multiply(spinPsil, spinPsil);
    auxmatrix.Scale(scalarproductPsilPsilprime *
                    (3 * sin(normPsil) - normPsil * (2 + cos(normPsil))) /
                    (normPsil * normPsil * normPsil * normPsil * normPsil));
    dTinvdx += auxmatrix;
  }

  return;
}  // LARGEROTATIONS::computedTinvdx

/*-----------------------------------------------------------------------------------*
 |computes inverse quaternion q^{-1} for input quaternion q      (public)cyron02/09|
 *-----------------------------------------------------------------------------------*/
LINALG::Matrix<4, 1> LARGEROTATIONS::inversequaternion(const LINALG::Matrix<4, 1>& q)
{
  // square norm ||q||^2 of quaternion q
  double qnormsq = q.Norm2() * q.Norm2();

  // declaration of variable for inverse quaternion
  LINALG::Matrix<4, 1> qinv;

  // inverse quaternion q^(-1) = [-q0, -q1, -q2, q3] / ||q||^2;
  for (int i = 0; i < 3; i++) qinv(i) = -q(i) / qnormsq;

  qinv(3) = q(3) / qnormsq;

  return qinv;

}  // LARGEROTATIONS::inversequaternion

/*---------------------------------------------------------------------------------------------------*
 |quaternion product q12 = q2*q1, Crisfield, Vol. 2, equation (16.71); | |explanation: if q1 and q2
 correspond to the rotation matrices R1 and R2, respectively, the compound| |rotation R12 = R2*R1
 corresponds to the compound quaternion q12 = q2*q1          (public)cyron02/09|
 *---------------------------------------------------------------------------------------------------*/
// void LARGEROTATIONS::quaternionproduct(const LINALG::Matrix<4,1>& q1,const LINALG::Matrix<4,1>&
// q2,LINALG::Matrix<4,1>& q12)
//{
//  q12(0) = q2(3)*q1(0) + q1(3)*q2(0) + q2(1)*q1(2) - q1(1)*q2(2);
//  q12(1) = q2(3)*q1(1) + q1(3)*q2(1) + q2(2)*q1(0) - q1(2)*q2(0);
//  q12(2) = q2(3)*q1(2) + q1(3)*q2(2) + q2(0)*q1(1) - q1(0)*q2(1);
//  q12(3) = q2(3)*q1(3) - q2(2)*q1(2) - q2(1)*q1(1) - q2(0)*q1(0);
//} //LARGEROTATIONS::quaternionproduct

/*---------------------------------------------------------------------------------------------------*
 |computing the rotation angle theta which rotates the given unit direction vector d1 into the given
 | |unit direction vector d2; the unit direction vector of theta is the axis about which the
 rotation  | |from d1 to d2 is carried out and the absolute value of theta is the absolute value of
 the rotation | |angle (public) cyron07/10|
 *---------------------------------------------------------------------------------------------------*/
void LARGEROTATIONS::directionstoangle(
    const LINALG::Matrix<3, 1>& d1, const LINALG::Matrix<3, 1>& d2, LINALG::Matrix<3, 1>& theta)
{
  // check whether d1 and d2 are really unit direction vectors
  if (fabs(d1.Norm2() - 1) > 1e-8 || fabs(d2.Norm2() - 1) > 1e-8)
    dserror("d1 or d2 is not a unit direction vector!");

  /*theta is orthgonal to the plane defined by d1 and d2. Thus its direction follows from a
   *crossproduct between these two vectors*/
  theta(0) = d1(1) * d2(2) - d1(2) * d2(1);
  theta(1) = d1(2) * d2(0) - d1(0) * d2(2);
  theta(2) = d1(0) * d2(1) - d1(1) * d2(0);

  // compute scalar product between d1 and d2
  double scalarprod = 0.0;
  for (int j = 0; j < 3; j++) scalarprod += d1(j) * d2(j);

  // length of rotation angle vector theta given by the angle between d1 and d2
  theta.Scale(acos(scalarprod) / theta.Norm2());

}  // LARGEROTATIONS::directionstoangle

/*---------------------------------------------------------------------------------------------------*
 |computing the triad R which rotates the given unit direction vector d1 into the given unit
 direction| |vector d2 (public) cyron07/10|
 *---------------------------------------------------------------------------------------------------*/
void LARGEROTATIONS::directionstotriad(
    const LINALG::Matrix<3, 1>& d1, const LINALG::Matrix<3, 1>& d2, LINALG::Matrix<3, 3>& R)
{
  // computing rotation angle corresponding to triad R
  LINALG::Matrix<3, 1> theta;
  directionstoangle(d1, d2, theta);

  // computing rotation quaternion R corresponding to theta
  LINALG::Matrix<4, 1> q;
  angletoquaternion(theta, q);

  // computing rotation matrix R corresponding to theta
  quaterniontotriad(q, R);

}  // LARGEROTATIONS::directionstotriad

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
unsigned int LARGEROTATIONS::NumberingTrafo(const unsigned int j, const unsigned int numnode)
{
  // Node numbering j=1,...,NumNode() according to Crisfield 1999:
  // LIN2  1---2
  // LIN3  1---2---3
  // LIN4  1---2---3---4
  // LIN5  1---2---3---4---5

  // Storage position i=0,...,NumNode()-1 of nodal quantities applied in BACI:
  // LIN2  (1,2)
  // LIN3  (1,3,2)
  // LIN4  (1,4,2,3)
  // LIN5  (1,5,2,3,4)

  // Initialization
  unsigned int i = 0;

  // Transformation
  if (j == 1)
    i = 0;
  else if (j == numnode)
    i = 1;
  else
    i = j;

  return i;
}
