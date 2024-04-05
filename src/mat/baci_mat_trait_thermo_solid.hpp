/*! \file
\brief Interface for every material that can evaluate coupled thermo-solid material laws

\level 3

*/

#include "baci_config.hpp"

#include "baci_mat_trait_solid.hpp"
#include "baci_mat_trait_thermo.hpp"

#ifndef FOUR_C_MAT_TRAIT_THERMO_SOLID_HPP
#define FOUR_C_MAT_TRAIT_THERMO_SOLID_HPP

BACI_NAMESPACE_OPEN
namespace MAT
{
  namespace TRAIT
  {
    class ThermoSolid : public Thermo, public Solid
    {
     public:
      /*!
       * Set current quantities for this material
       *
       * The quantities are used for evaluation and possibly in CommitCurrentState()
       * @param defgrd
       * @param glstrain
       * @param temperature
       * @param gp
       *
       */
      virtual void Reinit(const CORE::LINALG::Matrix<3, 3>* defgrd,
          const CORE::LINALG::Matrix<6, 1>* glstrain, double temperature, unsigned gp) = 0;

      //! @name Coupled derivatives
      //! @{

      /// @brief  get derivative of 2nd PK stress wrt temperature
      ///
      /// this term arises for coupled thermo-mechanical materials
      virtual void GetdSdT(CORE::LINALG::Matrix<6, 1>* dS_dT) = 0;

      //! @}

      /*!
       * Return stress-temperature modulus and thermal derivative for coupled thermomechanics
       *
       * @param stm tensor to be filled with stress-temperature moduli
       * @param stm_deriv tensor to be filled with derivatives
       */
      virtual void StressTemperatureModulusAndDeriv(
          CORE::LINALG::Matrix<6, 1>& stm, CORE::LINALG::Matrix<6, 1>& stm_dT) = 0;
    };
  }  // namespace TRAIT
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif