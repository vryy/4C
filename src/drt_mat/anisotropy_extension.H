/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of anisotropy extension to be used by anisotropic materials with @MAT::Anisotropy

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ANISOTROPY_EXTENSION_H
#define BACI_ANISOTROPY_EXTENSION_H

#include <functional>
#include <utility>
#include <boost/optional.hpp>
#include "../linalg/linalg_fixedsizematrix.H"
#include "anisotropy.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "anisotropy_extension_base.H"
#include "anisotropy_fiber_provider.H"

namespace MAT
{
  /*!
   * @brief definition, which kind of fibers should be used (Element fibers or nodal (aka Gauss
   * point) fibers)
   */
  enum class FiberLocation
  {
    /// Undefined fiber location
    None,
    /// Fibers are constant per element
    ElementFibers,
    /// Fibers are defined on the GP
    GPFibers
  };

  /*!
   * \brief A fiber extension to be used with MAT::Anisotropy
   *
   * The anisotropy extension is the actual provider of the fibers for the material law.
   * MAT::Anisotropy initializes the anisotropy information (fibers and coordinate systems) and
   * notifies the MAT::FiberAnisotropyExtension to build their needed fibers and structural tensors.
   * Every anisotropic material should have an MAT::FiberAnisotropyExtension and should register it
   * to the Anisotropy framework.
   *
   * After the setup, the fibers and structural tensors (or any other needed quantity) is computed
   * and stored. Anisotropy extension can handle element and Gauss-point fibers
   *
   * \tparam numfib Number of fibers
   */
  template <unsigned int numfib>
  class FiberAnisotropyExtension : public BaseAnisotropyExtension, public FiberProvider
  {
    // Anisotropy is a friend to create back reference
    friend class Anisotropy;

   public:
    //! @name Tensors needed for the evaluation
    /// @{
    /// The material evaluation needs the fiber vectors
    static constexpr std::uint_fast8_t FIBER_VECTORS{1U << 0U};
    /// The material evaluation needs the structural tensors
    static constexpr std::uint_fast8_t STRUCTURAL_TENSOR{1U << 1U};
    /// The material evaluation needs the structural tensors in stress like voigt notation
    static constexpr std::uint_fast8_t STRUCTURAL_TENSOR_STRESS{1U << 2U};
    /// @}

    /*!
     * \brief Create an anisotropy extension with a structural tensor strategy
     *
     * \param stucturalTensorStrategy
     */
    explicit FiberAnisotropyExtension(
        const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy);

    /*!
     * \brief Create an anisotropy extension without a structural tensor strategy. Computing
     * structural tensors may not work.
     */
    FiberAnisotropyExtension();

    ///@name Packing and Unpacking
    /// @{

    /*!
     * \brief Pack all data for parallel distribution and restart
     *
     * \param data
     */
    void PackAnisotropy(DRT::PackBuffer& data) const override;

    /*!
     * \brief Unpack all data from parallel distribution or restart
     *
     * \param data whole data array
     * \param position position of the current reader
     */
    void UnpackAnisotropy(
        const std::vector<char>& data, std::vector<char>::size_type& position) override;
    /// @}

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
    const LINALG::Matrix<3, 1>& GetFiber(int gp, int i) const override;

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
    const LINALG::Matrix<6, 1>& GetStructuralTensor_stress(int gp, int i) const override;

    /**
     * \brief Returns the i-th structural tensor at the Integration point in tensor notation
     *
     * \note Use gp=#GPDEFAULT if element fibers are used
     *
     * @param gp (in) : Id of the integration point (use #GPDEFAULT for Element fibers)
     * @param i (in) : Id of the fiber
     * @return Reference to Matrix of the structural tensor in tensor notation
     */
    const LINALG::Matrix<3, 3>& GetStructuralTensor(int gp, int i) const override;
    //@}

    /*!
     * \brief Needed structural tensors should be registered here before setup so that only needed
     * structural tensors are computed
     *
     * Use FiberAnisotropyExtension::FIBER_VECTORS, FiberAnisotropyExtension::STRUCTURAL_TENSOR or
     * FiberAnisotropyExtension::STRUCTURAL_TENSOR_STRESS
     *
     * They can be combined with binary or, hence
     * FiberAnisotropyExtension::FIBER_VECTORS|FiberAnisotropyExtension::STRUCTURAL_TENSOR
     *
     * \param tensor_flags Flags of the needed structural tensors
     */
    void RegisterNeededTensors(std::uint_fast8_t tensor_flags) { tensor_flags_ |= tensor_flags; }

    /*!
     * \brief Returns the Id of the fiber to be returned at the Gauss point. If Element fibers are
     * used, this returns FiberAnisotropyExtension::GPDEFAULT. In case of Gauss point fibers, it
     * returns the Gauss point
     *
     * \param params Container that has stored the current Gauss point
     *
     * \return int
     */
    int GetVirtualGaussPoint(const Teuchos::ParameterList& params) const;

    /*!
     * \brief Returns the Id of the fiber to be returned at the Gauss point. If Element fibers are
     * used, this returns FiberAnisotropyExtension::GPDEFAULT. In case of Gauss point fibers, it
     * returns the Gauss point
     *
     * \param gp Current Gauss point
     *
     * \return int
     */
    int GetVirtualGaussPoint(int gp) const;

    /*!
     * \brief Returns the number of fibers per element. If Element fiber are used, it returns 1. In
     * case of Gauss point fibers, it is the number of Gauss points
     *
     * \return int
     */
    int GetFibersPerElement() const;

   protected:
    /*!
     * \brief Set all fibers at the Gauss point gp
     *
     * \param gp Gauss point. Use FiberAnisotropyExtension::GPDEFAULT in case of element fibers
     * \param fibers Vector of all fibers
     */
    void SetFibers(int gp, const std::array<LINALG::Matrix<3, 1>, numfib>& fibers);

    /*!
     * \brief Set all fibers of the element
     *
     * \param fibers The first index are the Gauss points, the second index the fibers. In case of
     * element fiebers, the first vector should only contain one element.
     */
    void SetFibers(const std::vector<std::array<LINALG::Matrix<3, 1>, numfib>>& fibers);

    /*!
     * \brief Method that compute all structural tensors. Should be executed after a change of the
     * fibers.
     */
    void ComputeNeededStructuralTensors();

    /*!
     * \brief Method that initializes element fibers.
     *
     * This method should only return true if element fibers are used.
     *
     * \return true if the fibers are initialized
     * \return false if the fibers are not initialized
     */
    virtual bool DoElementFiberInitialization() { return false; }

    /*!
     * \brief Method that initialized Gauss point fibers.
     *
     * This method should only return true if Gauss point fibers are used.
     *
     * \return true if the fibers are initialized
     * \return false if the fibers are not initialized
     */
    virtual bool DoGPFiberInitialization() { return false; }

    /*!
     * \brief Method that will be called of the fibers are initialized.
     */
    virtual void OnFibersInitialized()
    {
      // do nothing in the default case
    }

    /*!
     * \brief Set the flag where fibers lie (Gausspoint or Element fibers)
     *
     * \param location
     */
    void SetFiberLocation(FiberLocation location);

    /*!
     * \brief This method will be called by MAT::Anisotropy if element and Gauss point fibers are
     * available
     */
    void OnGlobalDataInitialized() override {}

   private:
    /*!
     * \brief This method will be called by MAT::Anisotropy to notify that element information is
     * available.
     */
    void OnGlobalElementDataInitialized() override
    {
      const bool initialized = DoElementFiberInitialization();
      if (initialized) OnFibersInitialized();
    }

    /*!
     * \brief This method will be called by MAT::Anisotropy to notify that Gauss point information
     * is available.
     */
    void OnGlobalGPDataInitialized() override
    {
      const bool initialized = DoGPFiberInitialization();
      if (initialized) OnFibersInitialized();
    }
    /// \}

    /*!
     * \brief Method that computes structural tensors from all given fibers
     */
    void ComputeStructuralTensors();

    /*!
     * \brief Method that computes all structural tensors in stress like Voigt notation from all
     * given fibers.
     */
    void ComputeStructuralTensors_stress();

    /// Indication of the fiber location
    FiberLocation fiberLocation_ = FiberLocation::None;

    /// Tensors needed for the evaluation
    std::uint_fast8_t tensor_flags_{};

    /**
     * Fibers of the element. The first index is for the Gauss points, the second index is for
     * the fiber id
     */
    std::vector<std::array<LINALG::Matrix<3, 1>, numfib>> fibers_;

    /**
     * Structural tensors of the fibers in stress like Voigt notation. The ordering is the same as
     * in #fibers_
     */
    std::vector<std::array<LINALG::Matrix<6, 1>, numfib>> fiberStructuralTensors_stress_;

    /**
     * Structural tensors of the fibers. The ordering is the same as in #fibers_
     */
    std::vector<std::array<LINALG::Matrix<3, 3>, numfib>> fiberStructuralTensors_;

    /// Structural tensor strategy
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase> structuralTensorStrategy_ =
        Teuchos::null;
  };
}  // namespace MAT

#endif  // BACI_ANISOTROPY_EXTENSION_H
