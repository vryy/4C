/*----------------------------------------------------------------------*/
/*! \file

\brief

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ANISOTROPY_DEFAULT_EXTENSION_H
#define BACI_ANISOTROPY_DEFAULT_EXTENSION_H

#include "anisotropy_extension.H"
namespace MAT
{
  /*!
   * \brief Default anisotropy extension
   *
   * This is an anisotropy extension that contains the functionality of most of the anisotropic
   * materials. Fibers are either initialized using the FIBERX notation, via a cylinder coordinate
   * system or in x direction.
   */
  class DefaultAnisotropyExtension : public FiberAnisotropyExtension
  {
   public:
    //! @name Initialization modes of the fibers
    /// @{
    //! Fibers defined in material on element basis
    static constexpr int INIT_MODE_ELEMENT_EXTERNAL = 0;
    //! Fibers defined in dat file on element basis
    static constexpr int INIT_MODE_ELEMENT_FIBERS = 1;
    //! Fibers defined in material on Gauss point basis
    static constexpr int INIT_MODE_NODAL_EXTERNAL = 4;
    //! Fibers defined in dat file on Gauss point basis
    static constexpr int INIT_MODE_NODAL_FIBERS = 3;
    /// @}

    /*!
     * \brief Constructor of the default anisotropy extension
     *
     * \param init_mode initialization mode. Use one of the following:
     * DefaultAnisotropyExtension::INIT_MODE_ELEMENT_EXTERNAL,
     * DefaultAnisotropyExtension::INIT_MODE_ELEMENT_FIBERS,
     * DefaultAnisotropyExtension::INIT_MODE_NODAL_EXTERNAL,
     * DefaultAnisotropyExtension::INIT_MODE_NODAL_FIBERS
     *
     * \param gamma angle
     *
     * \param adapt_angle flag whether the angle is subject to growth and remodeling
     *
     * \param stucturalTensorStrategy Reference to the sturctural tensor strategy
     */
    DefaultAnisotropyExtension(int init_mode, double gamma, bool adapt_angle,
        const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy);

    ~DefaultAnisotropyExtension() override = default;

    /*!
     * \brief Initializes the Element fibers
     *
     * \return true if the material is parametrized so that element fibers should be used
     * \return false otherwise
     */
    bool DoElementFiberInitialization() override;

    /*!
     * \brief Initializes Gauss point fibers
     *
     * \return true if the materials is parametrized so that Gauss point fibers should be used
     * \return false otherwise
     */
    bool DoGPFiberInitialization() override;

    /*!
     * \brief Set Fiber vectors by a new angle gamma in the current configuration
     *
     * \note this method is here for backwards compatibility
     *
     * \param newgamma New angle
     * \param locsys local coordinate system
     * \param defgrd deformation gradient
     */
    void SetFiberVecs(
        double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd);

    /*!
     * \brief Set the new element fibers directly
     *
     * \param fibervec unit vector pointing in fiber direction
     */
    void SetFiberVecs(const LINALG::Matrix<3, 1>& fibervec);

   protected:
    /*!
     * \brief The single fiber is aligned in z-direction
     */
    virtual void DoExternalFiberInitialization();

   private:
    /// Initialization mode
    const int init_mode_;

    /// angle of the fiber
    const double gamma_;

    /// flag whether the angle should be adapted with growth and remodeling (for comoatibility
    /// reasons)
    const bool adapt_angle_;
  };
}  // namespace MAT

#endif  // BACI_ANISOTROPY_DEFAULT_EXTENSION_H
