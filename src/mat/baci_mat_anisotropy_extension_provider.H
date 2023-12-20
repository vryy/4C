/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an provider of an MAT::FiberAnisotropyExtension

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MAT_ANISOTROPY_EXTENSION_PROVIDER_H
#define BACI_MAT_ANISOTROPY_EXTENSION_PROVIDER_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

namespace MAT
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
    virtual FiberAnisotropyExtension<numfib>& GetFiberAnisotropyExtension() = 0;
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MAT_ANISOTROPY_EXTENSION_PROVIDER_H
