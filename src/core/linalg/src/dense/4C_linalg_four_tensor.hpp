/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of general 4C 4-tensor

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_FOUR_TENSOR_HPP
#define FOUR_C_LINALG_FOUR_TENSOR_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <array>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  /*!
   * @class FourTensor
   *
   * @brief Basic implementation of a 4-tensor that is templated on the problem dimension.
   *
   * @tparam dim  dimension of the problem
   */
  template <int dim>
  class FourTensor
  {
   public:
    //! constructor that initializes a 4-tensor object with only zero entries if setZero is true
    explicit FourTensor(const bool setZero = true)
    {
      if (setZero) four_tensor_ = {{{{}}}};
    };

    //! @name Access methods
    //@{

    /*!
     * @brief return mutable 4-tensor object
     *
     * @return 4-tensor
     */
    inline std::array<std::array<std::array<std::array<double, dim>, dim>, dim>, dim>& Get()
    {
      return four_tensor_;
    }

    /*!
     * @brief return const 4-tensor object
     *
     * @return 4-tensor
     */
    inline const std::array<std::array<std::array<std::array<double, dim>, dim>, dim>, dim>&
    GetConst() const
    {
      return four_tensor_;
    }

    /*!
     * @brief returns writeable element C_{i1,i2,i3,i4} from 4-tensor C
     *
     * @param[in] i1 index of first basis vector
     * @param[in] i2 index of second basis vector
     * @param[in] i3 index of third basis vector
     * @param[in] i4 index of fourth basis vector
     * @return C_{i1,i2,i3,i4}
     */
    inline double& operator()(int i1, int i2, int i3, int i4);

    /*!
     * @brief returns const element C_{i1,i2,i3,i4} from 4-tensor C
     *
     * @param[in] i1 index of first basis vector
     * @param[in] i2 index of second basis vector
     * @param[in] i3 index of third basis vector
     * @param[in] i4 index of fourth basis vector
     * @return C_{i1,i2,i3,i4}
     */
    inline const double& operator()(int i1, int i2, int i3, int i4) const;
    //@}

   private:
    // 4-tensor
    std::array<std::array<std::array<std::array<double, dim>, dim>, dim>, dim> four_tensor_;
  };

  template <int dim>
  inline double& FourTensor<dim>::operator()(int i1, int i2, int i3, int i4)
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (i1 >= dim or i2 >= dim or i3 >= dim or i4 >= dim)
      FOUR_C_THROW("Indices %i,%i,%i,%i out of range in FourTensor<%i>.", i1, i2, i3, i4, dim);
#endif
    return Get()[i1][i2][i3][i4];
  }

  template <int dim>
  inline const double& FourTensor<dim>::operator()(int i1, int i2, int i3, int i4) const
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (i1 >= dim or i2 >= dim or i3 >= dim or i4 >= dim)
      FOUR_C_THROW("Indices %i,%i,%i,%i out of range in FourTensor<%i>.", i1, i2, i3, i4, dim);
#endif
    return GetConst()[i1][i2][i3][i4];
  }
}  // namespace CORE::LINALG

FOUR_C_NAMESPACE_CLOSE

#endif
