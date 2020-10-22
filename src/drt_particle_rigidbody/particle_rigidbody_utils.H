/*---------------------------------------------------------------------------*/
/*! \file
\brief utils for rigid bodies
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RIGIDBODY_UTILS_H
#define PARTICLE_RIGIDBODY_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <cmath>

/*---------------------------------------------------------------------------*
 | declarations                                                              |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  namespace UTILS
  {
    //! @name collection of quaternion operations
    //@{

    /*!
     * \brief clear quaternion
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    template <class T>
    inline void quaternion_clear(T* q)
    {
      q[0] = 0.0;
      q[1] = 0.0;
      q[2] = 0.0;
      q[3] = 1.0;
    }

    /*!
     * \brief set quaternion
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    template <class T>
    inline void quaternion_set(T* q1, const T* q2)
    {
      q1[0] = q2[0];
      q1[1] = q2[1];
      q1[2] = q2[2];
      q1[3] = q2[3];
    }

    /*!
     * \brief invert quaternion
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    template <class T>
    inline void quaternion_invert(T* q1, const T* q2)
    {
      q1[0] = -q2[0];
      q1[1] = -q2[1];
      q1[2] = -q2[2];
      q1[3] = q2[3];
    }

    /*!
     * \brief quaternion product
     *
     * q12 = q2 * q1, Crisfield, Vol. 2, equation (16.71)
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    template <class T>
    inline void quaternion_product(T* q12, const T* q2, const T* q1)
    {
      q12[0] = q2[3] * q1[0] + q1[3] * q2[0] + q2[1] * q1[2] - q1[1] * q2[2];
      q12[1] = q2[3] * q1[1] + q1[3] * q2[1] + q2[2] * q1[0] - q1[2] * q2[0];
      q12[2] = q2[3] * q1[2] + q1[3] * q2[2] + q2[0] * q1[1] - q1[0] * q2[1];
      q12[3] = q2[3] * q1[3] - q2[2] * q1[2] - q2[1] * q1[1] - q2[0] * q1[0];
    }

    /*!
     * \brief get quaternion from angle
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    template <class T>
    inline void quaternion_from_angle(T* q, const T* phi)
    {
      double absphi = std::sqrt(phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);

      if (absphi > 1E-14)
      {
        double fac = std::sin(0.5 * absphi) / absphi;

        q[0] = fac * phi[0];
        q[1] = fac * phi[1];
        q[2] = fac * phi[2];
        q[3] = std::cos(0.5 * absphi);
      }
      else
        quaternion_clear(q);
    }

    /*!
     * \brief rotate vector with given quaternion
     *
     * Note that the three dimensional vectors v and w are considered as quaternions with a real
     * coordinate equal to zero and that the inverse of a quaterion q = [q1 q2 q3 q4] evaluates to
     * q^-1 = [-q1 -q2 -q3 q4]
     *
     * [w 0] = q * [v 0] * q^-1
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    template <class T>
    inline void quaternion_rotate_vector(T* w, const T* q, const T* v)
    {
      double qv[4];
      qv[0] = +q[3] * v[0] + q[1] * v[2] - v[1] * q[2];
      qv[1] = +q[3] * v[1] + q[2] * v[0] - v[2] * q[0];
      qv[2] = +q[3] * v[2] + q[0] * v[1] - v[0] * q[1];
      qv[3] = -q[2] * v[2] - q[1] * v[1] - q[0] * v[0];

      w[0] = -qv[3] * q[0] + q[3] * qv[0] - qv[1] * q[2] + q[1] * qv[2];
      w[1] = -qv[3] * q[1] + q[3] * qv[1] - qv[2] * q[0] + q[2] * qv[0];
      w[2] = -qv[3] * q[2] + q[3] * qv[2] - qv[0] * q[1] + q[0] * qv[1];
    }

    //@}

  }  // namespace UTILS
}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
