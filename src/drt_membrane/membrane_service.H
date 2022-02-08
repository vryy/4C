/*----------------------------------------------------------------------*/
/*! \file
\brief Various service routines related to membranes

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef MEMBRANE_SERVICE_H
#define MEMBRANE_SERVICE_H

#include "../linalg/linalg_fixedsizematrix.H"

namespace MEMBRANE
{
  /*!
   * @brief Get stress like voigt notation stress tensor from 3D stress matrix assuming plane stress
   * in the x-y plane
   *
   * @param stress (in) : Stress in Matrix notation
   * @param planeStressStressLike (out) : Stress in voigt notation assuming plane stress in the x-y
   * plane
   */
  inline void LocalPlaneStressToStressLikeVoigt(
      const LINALG::Matrix<3, 3>& stress, LINALG::Matrix<3, 1>& planeStressStressLike)
  {
    planeStressStressLike(0) = stress(0, 0);
    planeStressStressLike(1) = stress(1, 1);
    planeStressStressLike(2) = 0.5 * (stress(0, 1) + stress(1, 0));
  }

  /*!
   * @brief Subpart of the linearization assuming plane stress in the x-y plane
   *
   * @param cmat (in) : Full linearization
   * @param cmatred (out) : Reduced linearization assuming plane stress in the x-y plane
   */
  inline void LocalFourthTensorPlaneStressToStressLikeVoigt(
      const LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<3, 3>& cmatred)
  {
    cmatred(0, 0) = cmat(0, 0);
    cmatred(0, 1) = cmat(0, 1);
    cmatred(0, 2) = cmat(0, 3);
    cmatred(1, 0) = cmat(1, 0);
    cmatred(1, 1) = cmat(1, 1);
    cmatred(1, 2) = cmat(1, 3);
    cmatred(2, 0) = cmat(3, 0);
    cmatred(2, 1) = cmat(3, 1);
    cmatred(2, 2) = cmat(3, 3);
  }

}  // namespace MEMBRANE

#endif
