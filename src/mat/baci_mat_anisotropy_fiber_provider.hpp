/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of the a pure virtual class of a fiber provider

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ANISOTROPY_FIBER_PROVIDER_HPP
#define FOUR_C_MAT_ANISOTROPY_FIBER_PROVIDER_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace MAT
{
  /*!
   * @brief Pure abstract class that defines the interface of a fiber holder
   */
  class FiberProvider
  {
   public:
    virtual ~FiberProvider() = default;

    /// @name Getter methods for the fibers
    //@{
    /**
     * \brief Returns the i-th fiber vector at the Integration point
     *
     * \note Use gp=#GPDEFAULT if element fibers are used
     *
     * @param gp (in) : Id of the integration point (use #GPDEFAULT for Element fibers)
     * @param i (in) : Id of the fiber
     * @return Reference to the vector of the fiber
     */
    virtual const CORE::LINALG::Matrix<3, 1>& GetFiber(int gp, int i) const = 0;

    /**
     * \brief Returns the i-th structural tensor at the Integration point in stress-like Voigt
     * notation
     *
     * \note Use gp=#GPDEFAULT if element fibers are used
     *
     * @param gp (in) : Id of the integration point (use #GPDEFAULT for Element fibers)
     * @param i (in) : Id of the fiber
     * @return Martix of the structural tensor in stress-like Voigt notation
     */
    virtual const CORE::LINALG::Matrix<6, 1>& GetStructuralTensor_stress(int gp, int i) const = 0;

    /**
     * \brief Returns the i-th structural tensor at the Integration point in tensor notation
     *
     * \note Use gp=#GPDEFAULT if element fibers are used
     *
     * @param gp (in) : Id of the integration point (use #GPDEFAULT for Element fibers)
     * @param i (in) : Id of the fiber
     * @return Reference to Matrix of the structural tensor in tensor notation
     */
    virtual const CORE::LINALG::Matrix<3, 3>& GetStructuralTensor(int gp, int i) const = 0;
    //@}
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
