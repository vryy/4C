/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for coordinate systems

\level 3

*----------------------------------------------------------------------*/

#ifndef SRC_GEOMETRY_COORDINATE_SYSTEM_UTILS_H_
#define SRC_GEOMETRY_COORDINATE_SYSTEM_UTILS_H_

#include "linalg_fixedsizematrix.H"

namespace GEO
{
  /*!
   * create a right-handed orthornormal basis from a unit vector based on
   * https://doi.org/10.1080/10867651.1999.10487513
   * @param unitvec[i]    : unit vector based on which basis will be created (third basis vector)
   * @param basisvec_1[o] : first basis vector
   * @param basisvec_2[o] : second basis vector
   */
  void BuildOrthonormalBasisFromUnitVector(const LINALG::Matrix<3, 1>& unitvec,
      LINALG::Matrix<3, 1>& basisvec_1, LINALG::Matrix<3, 1>& basisvec_2)
  {
    if (std::abs(unitvec.Norm2() - 1.0) > 1.0e-14) dserror("given vector not normalized!");

    if ((std::abs(unitvec(0)) <= std::abs(unitvec(1))) and
        (std::abs(unitvec(0)) <= std::abs(unitvec(2))))
    {
      basisvec_1(0) = 0.0;
      basisvec_1(1) = -unitvec(2);
      basisvec_1(2) = unitvec(1);
    }
    else if (std::abs(unitvec(1)) <= std::abs(unitvec(2)))
    {
      basisvec_1(0) = -unitvec(2);
      basisvec_1(1) = 0.0;
      basisvec_1(2) = unitvec(0);
    }
    else
    {
      basisvec_1(0) = -unitvec(1);
      basisvec_1(1) = unitvec(0);
      basisvec_1(2) = 0.0;
    }

    basisvec_1.Scale(1.0 / basisvec_1.Norm2());

    basisvec_2.CrossProduct(unitvec, basisvec_1);
  }
}  // namespace GEO


#endif /* SRC_GEOMETRY_COORDINATE_SYSTEM_UTILS_H_ */
