/*! \file
\brief Interface for every material monolithic solid-scalar material (ssi, tsi)

\level 3

*/

#ifndef FOUR_C_MAT_MONOLITHIC_SOLID_SCALAR_MATERIAL_HPP
#define FOUR_C_MAT_MONOLITHIC_SOLID_SCALAR_MATERIAL_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class MonolithicSolidScalarMaterial
  {
   public:
    virtual ~MonolithicSolidScalarMaterial() = default;

    /*!
     * @brief Evaluates the added derivatives of the stress w.r.t. all scalars
     *
     * @param defgrad (in) : Deformation gradient
     * @param glstrain (in) : Green-Lagrange strain
     * @param params (in) : ParameterList for additional parameters
     * @param gp (in) : Gauss points
     * @param eleGID (in) : global element id
     * @return std::vector<std::optional<Core::LinAlg::Matrix<6, 1>>>
     */
    virtual Core::LinAlg::Matrix<6, 1> evaluate_d_stress_d_scalar(
        const Core::LinAlg::Matrix<3, 3>& defgrad, const Core::LinAlg::Matrix<6, 1>& glstrain,
        Teuchos::ParameterList& params, int gp, int eleGID) = 0;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif