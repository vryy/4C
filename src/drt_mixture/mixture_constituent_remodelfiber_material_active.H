/*----------------------------------------------------------------------*/
/*! \file
\brief Base class for active remodel fiber materials.
\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_ACTIVE_H
#define BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_ACTIVE_H

#include "mixture_constituent_remodelfiber_material.H"

namespace MIXTURE

{
  template <typename T>
  class RemodelFiberMaterialActive
  {
   public:
    /*!
     * @brief Evaluates the first derivative of the strain energy function w.r.t. the active stretch
     *
     * @return T
     */
    [[nodiscard]] virtual T GetFirstDerivativeActive() const = 0;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_ACTIVE_H