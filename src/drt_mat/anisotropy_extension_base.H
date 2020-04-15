/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a base anisotropy extension to be used by anisotropic materials with
@MAT::Anisotropy

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_BASE_ANISOTROPY_EXTENSION_H
#define BACI_BASE_ANISOTROPY_EXTENSION_H

#include <Teuchos_RCPDecl.hpp>

// forward declarations
namespace DRT
{
  class PackBuffer;
}
namespace MAT
{
  // forward declaration
  class Anisotropy;

  class BaseAnisotropyExtension
  {
    // Anisotropy is a friend to create back reference
    friend class Anisotropy;

   public:
    /// If element fibers are used, they are stored at the beginning of the list
    static constexpr int GPDEFAULT = 0;

    virtual ~BaseAnisotropyExtension() = default;

    ///@name Packing and Unpacking
    /// @{

    /*!
     * \brief Pack all data for parallel distribution and restart
     *
     * \param data
     */
    virtual void PackAnisotropy(DRT::PackBuffer& data) const = 0;

    /*!
     * \brief Unpack all data from parallel distribution or restart
     *
     * \param data whole data array
     * \param position position of the current reader
     */
    virtual void UnpackAnisotropy(
        const std::vector<char>& data, std::vector<char>::size_type& position) = 0;
    /// @}

    /*!
     * \brief This method will be called by MAT::Anisotropy if element and Gauss point fibers are
     * available
     */
    virtual void OnGlobalDataInitialized() = 0;

   protected:
    /*!
     * \brief Returns the reference to the anisotropy
     *
     * \return Teuchos::RCP<Anisotropy>& Reference to the anisotropy
     */
    Teuchos::RCP<Anisotropy>& GetAnisotropy() { return anisotropy_; }
    /*!
     * \brief Returns the reference to the anisotropy
     *
     * \return Teuchos::RCP<Anisotropy>& Reference to the anisotropy
     */
    const Teuchos::RCP<Anisotropy>& GetAnisotropy() const { return anisotropy_; }

   private:
    /*!
     * \brief This method will be called by MAT::Anisotropy to notify that element information is
     * available.
     */
    virtual void OnGlobalElementDataInitialized() = 0;


    /*!
     * \brief This method will be called by MAT::Anisotropy to notify that Gauss point information
     * is available.
     */
    virtual void OnGlobalGPDataInitialized() = 0;

    /// \name Private methods called by the friend class MAT::Anisotropy
    /// \{
    /*!
     * \brief Set the anisotropy. This method will only be used by Anisotropy itself to give the
     * extension access to all anisotropy information.
     *
     * \param anisotropy
     */
    void SetAnisotropy(Anisotropy& anisotropy) { anisotropy_ = Teuchos::rcpFromRef(anisotropy); }

    /// Reference to Anisotropy
    Teuchos::RCP<Anisotropy> anisotropy_;
  };
}  // namespace MAT
#endif  // BACI_BASE_ANISOTROPY_EXTENSION_H