/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains a summand of the ElastHyper Toolbox with active behavior

\level 3

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/
#ifndef ELAST_ACTIVESUMMAND_H_
#define ELAST_ACTIVESUMMAND_H_
#include "elast_summand.H"

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
      ~ActiveSummand() override = default;
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
      virtual void AddActiveStressCmatAniso(
          const LINALG::Matrix<3, 3>& CM,  ///< rigtht Cauchy Green tensor
          LINALG::Matrix<6, 6>& cmat,      ///< material stiffness matrix
          LINALG::Matrix<6, 1>& stress,    ///< 2nd PK-stress
          int gp,                          ///< Gauss point
          int eleGID) const = 0;           ///< element GID


      /*!
       * \brief Evaluates the first derivative of active fiber potential w.r.t. the active fiber
       * stretch
       *
       * \return double First derivative of active fiber potential w.r.t. the active fiber stretch
       */
      virtual double GetDerivativeAnisoActive() const = 0;
    };
  }  // namespace ELASTIC

}  // namespace MAT


#endif  // ELAST_ACTIVESUMMAND_H_