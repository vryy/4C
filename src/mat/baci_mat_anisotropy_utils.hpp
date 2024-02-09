/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for the Anisotropy classes

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MAT_ANISOTROPY_UTILS_HPP
#define BACI_MAT_ANISOTROPY_UTILS_HPP

#include "baci_config.hpp"

#include "baci_io_linedefinition.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

#include <Teuchos_RCPDecl.hpp>

#include <string>

BACI_NAMESPACE_OPEN

// forward declaration
namespace CORE::COMM
{
  class PackBuffer;
}
namespace INPUT
{
  class LineDefinition;
}

namespace MAT
{
  // forward declaration
  namespace ELASTIC
  {
    class StructuralTensorStrategyBase;
  }
  /*!
   * Reads a fiber with a specification from the input file definition
   *
   * @param linedef (in) : Input line definition
   * @param specifier (in) : Identifier of the fiber
   * @param fiber_vector (out) : Fiber vector
   */
  void ReadAnisotropyFiber(INPUT::LineDefinition* linedef, std::string specifier,
      CORE::LINALG::Matrix<3, 1>& fiber_vector);

  /*!
   * \brief Compute structural tensors of a 2D vector of fibers with the structural tensor
   * strategy
   *
   * \tparam T Output type of the structural tensor (either matrix notation or stress-like Voigt
   * notation)
   *
   * \tparam numfib number of fibers
   *
   * \param fibers 2D vector of fibers (3x1 matrices)
   * \param structural_tensor 2D vector of structural tensors (3x3 or 6x1 matrices)
   * \param strategy Reference to the structural tensor strategy
   */
  template <typename T, unsigned int numfib>
  void ComputeStructuralTensors(std::vector<std::array<CORE::LINALG::Matrix<3, 1>, numfib>>& fibers,
      std::vector<std::array<T, numfib>>& structural_tensor,
      const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy);

  /*!
   * \brief Pack 2D vector of fibers and structural tensors
   *
   * \tparam T Type of the fiber (3x1, 6x1 or 3x3 matrices)
   * \param buffer buffer where to pack the data
   * \param vct vector
   */
  template <typename T>
  void PackFiberVector(CORE::COMM::PackBuffer& buffer, const std::vector<std::vector<T>>& vct);

  /*!
   * \brief Pack 2D vector of fibers and structural tensors
   *
   * \tparam T Type of the fiber (3x1, 6x1 or 3x3 matrices)
   * \tparam numfib number of fibers
   * \param buffer buffer where to pack the data
   * \param vct vector
   */
  template <typename T, unsigned int numfib>
  void PackFiberArray(
      CORE::COMM::PackBuffer& buffer, const std::vector<std::array<T, numfib>>& vct);

  /*!
   * \brief Unpack 2D vector of fibers and structural tensors
   *
   * \tparam T Type of the fiber (3x1, 6x1 or 3x3 matrices)
   * \param position Position where to start to unpack the data
   * \param data data where to unpack the data from
   * \param vct destination 2D array
   */
  template <typename T>
  void UnpackFiberVector(std::vector<char>::size_type& position, const std::vector<char>& data,
      std::vector<std::vector<T>>& vct);

  /*!
   * \brief Unpack 2D vector of fibers and structural tensors
   *
   * \tparam T Type of the fiber (3x1, 6x1 or 3x3 matrices)
   * \tparam numfib number of fibers
   * \param position Position where to start to unpack the data
   * \param data data where to unpack the data from
   * \param vct destination 2D array
   */
  template <typename T, unsigned int numfib>
  void UnpackFiberArray(std::vector<char>::size_type& position, const std::vector<char>& data,
      std::vector<std::array<T, numfib>>& vct);
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_ANISOTROPY_UTILS_H
