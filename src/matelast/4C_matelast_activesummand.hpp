/*----------------------------------------------------------------------*/
/*! \file
\brief Summand of the ElastHyper Toolbox with active behavior

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATELAST_ACTIVESUMMAND_HPP
#define FOUR_C_MATELAST_ACTIVESUMMAND_HPP
#include "4C_config.hpp"

#include "4C_matelast_summand.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    /*!
     * \brief This is a pure abstract extension of the Summand class to be used for active
     * materials.
     */
    class ActiveSummand : public Summand
    {
     public:
      /*!
       * \brief adds the active part of the summand to the stress and stiffness matrix
       *
       * \note This is ONLY the ACTIVE PART of the stress response
       *
       * \param CM Right Cauchy-Green deformation tensor in Tensor notation
       * \param cmat Material stiffness matrix
       * \param stress 2nd PK-stress
       * \param eleGID global element id
       */
      virtual void add_active_stress_cmat_aniso(
          const CORE::LINALG::Matrix<3, 3>& CM,  ///< right Cauchy Green tensor
          CORE::LINALG::Matrix<6, 6>& cmat,      ///< material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,    ///< 2nd PK-stress
          int gp,                                ///< Gauss point
          int eleGID) const = 0;                 ///< element GID

      /*!
       * \brief Evaluates the first derivative of active fiber potential w.r.t. the active fiber
       * stretch
       *
       * \return double First derivative of active fiber potential w.r.t. the active fiber stretch
       */
      virtual double get_derivative_aniso_active() const = 0;
    };
  }  // namespace ELASTIC

}  // namespace MAT


FOUR_C_NAMESPACE_CLOSE

#endif