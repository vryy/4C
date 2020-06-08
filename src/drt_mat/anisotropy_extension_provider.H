/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of an provider of an MAT::FiberAnisotropyExtension

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FIBER_ANISOTROPY_EXTENSION_PROVIDER_H
#define BACI_FIBER_ANISOTROPY_EXTENSION_PROVIDER_H

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

#endif  // BACI_FIBER_ANISOTROPY_EXTENSION_PROVIDER_H
