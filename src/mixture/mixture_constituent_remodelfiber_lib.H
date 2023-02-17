/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of helper functions for the remodel fiber constituent
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef MIXTURE_CONSTITUENT_REMODELFIBER_LIB_H
#define MIXTURE_CONSTITUENT_REMODELFIBER_LIB_H

namespace MIXTURE::PAR
{
  template <typename T>
  class RemodelFiberMaterial;

  /*!
   * @brief Create a remodel fiber material from the material id given in the input file
   *
   * @param matid material id given in the input file
   * @return const MIXTURE::PAR::RemodelFiberMaterial<double>*
   */
  [[nodiscard]] const MIXTURE::PAR::RemodelFiberMaterial<double>* FiberMaterialFactory(int matid);
}  // namespace MIXTURE::PAR

#endif  // MIXTURE_CONSTITUENT_REMODELFIBER_LIB_H
