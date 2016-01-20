/*!----------------------------------------------------------------------
\file largerotations.H

\brief A set of utility functions for large rotations

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Christian J. Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifndef LARGEROTATIONS_DEFINES_H
#define LARGEROTATIONS_DEFINES_H

#include "../linalg/linalg_fixedsizematrix.H"


/************************************************************************/
/* Large rotations functions                                            */
/************************************************************************/
/*!
\brief LARGEROTATIONS: namespace for operations with finite rotations according
to Crisfield 2003, Non-linear Finite Element Analysis of Solids and Structures,
Volume 2, Chapter 16 and according to  Jelenic 1999, Geometrically exact 3D beam
theory: implementation of a strain-invariant finite element for statics and dynamics
(Computational Methods in Applied Mechanics and Engineering 171, p. 141 - 171)
*/
namespace LARGEROTATIONS
{
  //!computing spin matrix S from rotation angle vector theta, Crisfield, Vol. 2, equation (16.8)
  void computespin(LINALG::Matrix<3,3>& S, const LINALG::Matrix<3,1>& theta);

  //!this is a templated version of the method computespin as declared above in order to also handle
  //!the datatype FAD which is needed for automatic differentiation
  template<typename T>
  void computespin(LINALG::TMatrix<T,3,3>& S, const LINALG::TMatrix<T,3,1>& theta)
  {
    //function based on Crisfield Vol. 2, Section 16 (16.8)
    S(0,0) = 0;
    S(0,1) = -theta(2);
    S(0,2) = theta(1);
    S(1,0) = theta(2);
    S(1,1) = 0;
    S(1,2) = -theta(0);
    S(2,0) = -theta(1);
    S(2,1) = theta(0);
    S(2,2) = 0;

    return;
  }

  //!computing rotation matrix R from quaternion q, Crisfield, Vol. 2, equation (16.70)
  void quaterniontotriad(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,3>& R);

  //!computing quaternion q from rotation angle theta, Crisfield, Vol. 2, equation (16.67)
  void angletoquaternion(const LINALG::Matrix<3,1>& theta, LINALG::Matrix<4,1>& q);

  //!computing rotation matrix R from angle theta, Crisfield Vol. 2, equation (16.22)
  void angletotriad(const LINALG::Matrix<3,1>& theta, LINALG::Matrix<3,3>& R);

  //!computing rotation angle vector theta from quaternion q, Crisfield, Vol. 2, equation (16.67)
  void quaterniontoangle(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& theta);

  //!computing from quaternion q the Rodrigues parameters omega, Crisfield, Vol. 2, equation (16.79)
  void quaterniontorodrigues(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& omega);

  //!computing quaternion q from a rotation matrix R, Crisfield, Vol. 2, equation (16.74) - (16.78)
  void triadtoquaternion(const LINALG::Matrix<3,3>& R, LINALG::Matrix<4,1>& q);

  //!computing the rotation angle theta rotating the unit direction vector d1 into the unit direction vector d2
  void directionstoangle(const LINALG::Matrix<3,1>& d1,const LINALG::Matrix<3,1>& d2, LINALG::Matrix<3,1>& theta);

  //!computing the triad R rotating the unit direction vector d1 into the unit direction vector d2
  void directionstotriad(const LINALG::Matrix<3,1>& d1,const LINALG::Matrix<3,1>& d2, LINALG::Matrix<3,3>& R);

  //!quaternion product q12 = q2*q1, Crisfield, Vol. 2, equation (16.71)
  void quaternionproduct(const LINALG::Matrix<4,1>& q1,const LINALG::Matrix<4,1>& q2,LINALG::Matrix<4,1>& q12);

  //!compute inverse quaternion q^{-1} of quaternion q
  LINALG::Matrix<4,1> inversequaternion(const LINALG::Matrix<4,1>& q);

  //!matrix T(\theta) from Jelenic 1999, eq. (2.5), equivalent to matrix H^{-1} in Crisfield, (16.94)
  LINALG::Matrix<3,3> Tmatrix(LINALG::Matrix<3,1> theta);

  //!matrix T(\theta)^{-1} from Jelenic 1999, eq. (2.5)
  LINALG::Matrix<3,3> Tinvmatrix(LINALG::Matrix<3,1> theta);

  //!compute dT^{-1}(\theta)/dx according to the two-lined equation below (3.19) on page 152 of Jelenic 1999
  void computedTinvdx(const LINALG::Matrix<3,1>& Psil, const LINALG::Matrix<3,1>& Psilprime, LINALG::Matrix<3,3>& dTinvdx);

  //!Transformation from node number according to Crisfield 1999 to storage position applied in BACI
  int NumberingTrafo(const int j, const int numnode);

} // LARGEROTATIONS

#endif  // #ifndef LARGEROTATIONS_DEFINES_H
