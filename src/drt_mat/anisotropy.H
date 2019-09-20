/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of common functionality for anisotropic materials

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ANISOTROPY_H
#define BACI_ANISOTROPY_H

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_pack_buffer.H"

namespace MAT
{
  // forward declaration of structural tensor strategy
  namespace ELASTIC
  {
    class StructuralTensorStrategyBase;
  }

  /*!
   * Class that handles the initialization of the anisotropic parts of a material model.
   *
   * \note: This is an optional class. There are anisotropic materials that do not use this
   * interface.
   */
  class Anisotropy
  {
   public:
    //! Numerical tolerance used to check for unit vectors
    static constexpr float TOLERANCE = 1e-9;

    //! @name Initialization modes of the fibers
    //@{
    //! External fiber initialization (Material specific)
    static constexpr int INIT_MODE_EXTERNAL = 0;

    //! Fibers are defined in the DAT-file on the elements
    static constexpr int INIT_MODE_ELEMENT_FIBERS = 1;

    //! Fibers are defined in the DAT file on the nodes
    static constexpr int INIT_MODE_NODAL_FIBERS = 3;
    //@}

    /*!
     * For element fibers, the position of the fibers at the GP are all at the beginning of the list
     */
    static constexpr int GPDEFAULT = 0;

    /*!
     * Fiber definition mode
     */
    enum DefinitionMode
    {
      undefined,  //!< Unknown definition mode
      fiberi,     //!< Fibers are defined with the FIBERi notation in the input file
      radaxicir,  //!< Fibers are defined with RAD-AXI-CIR notation in the input file
      external    //!< Fibers are defined external (likely by the material)
    };

    /*!
     * Minimal constructor
     *
     * \note When using this constructor, the initialization mode and the structural tensor strategy
     * have to be set with #Initialize(int, const
     * Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>&)
     *
     * @param number_fibers Number of fibers needed for the material
     */
    explicit Anisotropy(int number_fibers);

    /*!
     * Constructor to directly set the material parameters
     *
     * @param number_fibers (in) : Number of fibers needed for the material
     * @param init (in) : Initialization mode
     * @param structuralTensorStrategy (in) : Strategy of the structural tensor
     */
    Anisotropy(int number_fibers, int init,
        const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy);

    ///@name Packing and Unpacking
    //@{

    /*!
     * Pack all data for parallel distribution
     *
     * @param data (in/out) : data object
     */
    virtual void PackAnisotropy(DRT::PackBuffer& data) const;

    /*!
     * Unpack all data from another processor
     *
     * @param data (in) : data object
     * @param position (in/out) : current position in the data
     */
    virtual void UnpackAnisotropy(
        const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /*!
     * Initializes the fiber types
     *
     * @param init (in) : initialization mode of the fibers (0: no fibers given, 1: Element fibers,
     * 3: Nodal fibers)
     * @param structuralTensorStrategy (in) : Strategy to build the structural tensor
     */
    virtual void Initialize(int init,
        const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy);

    /*!
     * This method should be called as soon as the number of Gauss points is known.
     * @param numgp Number of Gauss points
     */
    void SetNumberOfGaussPoints(int numgp);

    /*!
     * Reads the line definition of an element to get the fibers defined on the element
     *
     * @param linedef (in) : Input line of the corresponding element
     */
    virtual void ReadAnisotropyFromElement(DRT::INPUT::LineDefinition* linedef);

    /*!
     * This method extracts the Gauss-point fibers written by the elements into the ParameterList
     * and stores them internally. This method should be called during the PostSetup. This method
     * will only check for Gauss-point fibers if the initialization mode is #INIT_MODE_NODAL_FIBERS.
     *
     * @param params Container that hold the Gauss-point fibers.
     */
    virtual void ReadAnisotropyFromParameterList(Teuchos::ParameterList& params);

    /**
     * Sets the fiber vectors at the Gauss point level
     *
     * @param gp (in) : Id of the Gauss point
     * @param fibers (in) : Vector of fibers
     */
    virtual void SetFibers(int gp, const std::vector<LINALG::Matrix<3, 1>>& fibers);

    /**
     * Method to check the initialization status of the fibers
     * @return true if the fibers are already initialized, otherwise false
     */
    virtual bool FibersInitialized() { return fibers_initialized_; }

    /**
     * Do the fiber initialization. Needed to be overriden by the base class
     */
    virtual void DoFiberInitialization()
    {
      // empty in base class. Needed for external fiber initialization
      dserror(
          "You need to override this method to implement the external fiber initialization. Call "
          "SetFibers(...) from here to set the fibers.");
    }

    /**
     * \brief A notifier method that is called after the initialization of the fibers.
     *
     * Could be used to setup structural tensors or other quantities that only need to be computed
     * once per fiber initialization
     */
    virtual void OnFibersInitialized()
    {
      // empty in base class. Override if needed.
    }

    /**
     * \brief Setup fibers by the given local coordinate system
     *
     * This method will be called after the local coordinate system was read from the RAD-ACI-CIR
     * notation.
     *
     * The columns of #locsys matrix are:
     *  1.: radial unit vector
     *  2.: axial unit vector
     *  3.: circular unit vector
     *
     * @param locsys Local coordinate system
     */
    virtual void SetupFiberByCosy(LINALG::Matrix<3, 3>& locsys)
    {
      dserror(
          "This material does not support to set the fiber by the coordinate system. The material "
          "should override SetupFiberByCosy(...) for the setup.");
    }

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
    LINALG::Matrix<3, 1>& GetFiber(int gp = GPDEFAULT, int i = 0);

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
    LINALG::Matrix<6, 1>& GetStructuralTensor_stress(int gp = GPDEFAULT, int i = 0);

    /**
     * \brief Returns the i-th structural tensor at the Integration point in tensor notation
     *
     * \note Use gp=#GPDEFAULT if element fibers are used
     *
     * @param gp (in) : Id of the integration point (use #GPDEFAULT for Element fibers)
     * @param i (in) : Id of the fiber
     * @return Reference to Matrix of the structural tensor in tensor notation
     */
    LINALG::Matrix<3, 3>& GetStructuralTensor(int gp = GPDEFAULT, int i = 0);
    //@}

   protected:
    /*!
     * Reads a fiber with a specification from the input file definition
     *
     * @param linedef (in) : Input line definition
     * @param specifier (in) : Identifier of the fiber
     * @param fiber_vector (out) : Fiber vector
     */
    static void ReadAnisotropyFiber(DRT::INPUT::LineDefinition* linedef, std::string specifier,
        LINALG::Matrix<3, 1>& fiber_vector);

    /*!
     * Returns the number of fibers per element.
     *
     * If the fibers are constant per element, this is
     * 1. In case of nodal fibers, this is the number of integration points.
     *
     * @return Number of fibers per element
     */
    int GetFibersPerElement();

    /*!
     * Returns the current Gauss point if nodal fibers are used, otherwise returns #GPDEFAULT
     *
     * @param params (in) : Parameter list that has to contain the Integration point
     * @return Current Gauss point id if nodal fibers are used, otherwise #GPDEFAULT
     */
    int GetGPId(Teuchos::ParameterList& params);

   private:
    template <typename T>
    void PackFiberVector(DRT::PackBuffer& buffer, const std::vector<std::vector<T>>& vct) const;

    template <typename T>
    void UnpackFiberVector(std::vector<char>::size_type& position, const std::vector<char>& data,
        std::vector<std::vector<T>>& vct) const;

    /*!
     * Compute and store all structural tensors in stress like Voigt notation for all given fibers
     * internally.
     */
    void ComputeStructuralTensors_stress();

    /*!
     * Compute and store all structural tensors in tensor notation for all given fibers internally.
     */
    void ComputeStructuralTensors();

    /// Number of Gauss points
    int numgp_ = 0;

    /// Total number of fibers
    const unsigned number_fibers_ = 0;

    /// Initialization mode of the fibers
    int init_mode_ = 0;

    /**
     * Fibers of the element. The first index is for the Gauss points, the second index is for the
     * fiber id
     */
    std::vector<std::vector<LINALG::Matrix<3, 1>>> fibers_;

    /**
     * Structural tensors of the fibers in stress like Voigt notation. The ordering is the same as
     * in #fibers_
     */
    std::vector<std::vector<LINALG::Matrix<6, 1>>> structuralTensors_stress_;

    /**
     * Structural tensors of the fibers. The ordering is the same as in #fibers_
     */
    std::vector<std::vector<LINALG::Matrix<3, 3>>> structuralTensors_;

    /// Flag whether the fibers were set
    bool fibers_initialized_;

    /// mode of the definition of fibers
    DefinitionMode definitionMode_;

    /// Structural tensor strategy
    Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase> structuralTensorStrategy_ = Teuchos::null;
  };
}  // namespace MAT

#endif  // BACI_ANISOTROPY_H
