/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an provider of an Mat::FiberAnisotropyExtension

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ANISOTROPY_EXTENSION_PROVIDER_HPP
#define FOUR_C_MAT_ANISOTROPY_EXTENSION_PROVIDER_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration
  template <unsigned int numfib>
  class FiberAnisotropyExtension;

  /*!
   * \brief An pure abstract interface for an FiberAnisotropyExtension provider
   *
   * \tparam numfib Number of fibers
   */
  template <unsigned int numfib>
  class FiberAnisotropyExtensionProvider
  {
   public:
    /*!
     * \brief Default destructor
     */
    virtual ~FiberAnisotropyExtensionProvider() = default;

    /*!
     * \brief Returns the reference to an FiberAnisotropyExtension
     *
     * \return FiberAnisotropyExtension& Reference to a fiber anisotropy extension
     */
    virtual FiberAnisotropyExtension<numfib>& get_fiber_anisotropy_extension() = 0;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
