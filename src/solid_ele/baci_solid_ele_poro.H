/*! \file

\brief The declaration for the solid-poro element

 This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

#ifndef BACI_SOLID_ELE_PORO_H
#define BACI_SOLID_ELE_PORO_H

#include "baci_inpar_poro.H"
#include "baci_inpar_scatra.H"
#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_solid_ele_poro_calc_interface.H"
#include "baci_structure_new_elements_paramsinterface.H"

#include <memory>


namespace MAT
{
  class StructPoro;
  class FluidPoroMultiPhase;
}  // namespace MAT

namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidEleCalcInterface;

  class SolidPoroType : public DRT::ElementType
  {
   public:
    /// setup the dat file input line definitions for this type of element
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        override;

    /// create an element from a dat file specifier
    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
        const int id, const int owner) override;

    /// create an empty element
    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    // create and unpack element from data
    DRT::ParObject* Create(const std::vector<char>& data) override;

    std::string Name() const override { return "SolidPoroType"; }

    void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
        DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static SolidPoroType& Instance();

   private:
    static SolidPoroType instance_;

  };  // class SolidPoroType


  class SolidPoro : public DRT::Element
  {
    friend class SolidPoroType;

   public:
    //! @name Constructors and destructors and related methods
    //!@{

    /*!
    \brief Standard Constructor

    \param id : A unique global id
    \param owner : elements owner
    */
    SolidPoro(int id, int owner);


    //! destructor
    ~SolidPoro() = default;

    //! Copy Constructor
    SolidPoro(const SolidPoro& other);

    //! copy assignment operator
    SolidPoro& operator=(const SolidPoro& other);

    //! move constructor
    SolidPoro(SolidPoro&& other) noexcept = default;

    //! move assignment operator
    SolidPoro& operator=(SolidPoro&& other) noexcept = default;

    //!@}

    /*!
    \brief Deep copy
    */
    DRT::Element* Clone() const override;


    /*!
    \brief Return unique ParObject id
    */
    [[nodiscard]] int UniqueParObjectId() const override
    {
      return SolidPoroType::Instance().UniqueParObjectId();
    };

    /*!
    \brief Return number of lines to this element
    */
    [[nodiscard]] int NumLine() const override;

    /*!
    \brief Return number of surfaces to this element
    */
    [[nodiscard]] int NumSurface() const override;

    /*!
    \brief Return number of volumes to this element
    */
    [[nodiscard]] int NumVolume() const override;

    /*!
    \brief Get vector of Teuchos::RCPs to the lines of this element

    */
    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    /*!
    \brief Get vector of Teuchos::RCPs to the surfaces of this element

    */
    std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

    /*!
    \brief Get vector of Teuchos::RCPs to the volumes of this element

    */
    std::vector<Teuchos::RCP<DRT::Element>> Volumes() override;

    /*!
    \brief Get number of degrees of freedom of a certain node
    */
    [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override
    {
      return 3;
    }  // todo fix this for 2D elements


    /*!
    \brief Get number of degrees of freedom per element
           (implements pure virtual DRT::Element)

    The element decides how many element degrees of freedom it has.
    It can redecide along the way of a simulation.

    \note Element degrees of freedom mentioned here are dofs that are visible
          at the level of the total system of equations. Purely internal
          element dofs that are condensed internally should NOT be considered.
    */
    [[nodiscard]] int NumDofPerElement() const override { return 0; }

    /*!
\   brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this element

            */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Unpack(const std::vector<char>& data) override;

    /*!
    \brief Set kinematic type of element
    */
    void SetKinematicType(INPAR::STR::KinemType kintype) { kintype_ = kintype; }

    /*!
    \brief Get shape type of element
    */
    [[nodiscard]] DRT::Element::DiscretizationType Shape() const override { return distype_; };

    //! return ElementType instance
    [[nodiscard]] DRT::ElementType& ElementType() const override
    {
      return SolidPoroType::Instance();
    }


    /*!
    \brief Read input for this element
    */
    bool ReadElement(const std::string& eletype, const std::string& distype,
        DRT::INPUT::LineDefinition* linedef) override;

    /*!
    \brief Evaluate an element

    \param params (in/out): ParameterList for communication between control routine
                            and elements
    \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                            the controling method does not epxect the element to fill
                            this matrix.
    \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                            the controling method does not epxect the element to fill
                            this matrix.
    \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                            the controlling method does not epxect the element
                            to fill this vector
    \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                            the controlling method does not epxect the element
                            to fill this vector
    \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                            the controlling method does not epxect the element
                            to fill this vector
    \return 0 if successful, negative otherwise
    */
    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
        CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseVector& elevec2,
        CORE::LINALG::SerialDenseVector& elevec3) override;

    /*!
    \brief Evaluate Neumann boundary condition

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): reference to the underlying discretization
    \param condition (in)     : condition to be evaluated
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

    \return 0 if successful, negative otherwise
    */
    int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

    //! read anisotropic permeability directions in the element definition
    void ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(
        DRT::INPUT::LineDefinition* linedef);

    //! read nodal anisotropic permeability scaling coefficients in the element definition
    void ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(
        DRT::INPUT::LineDefinition* linedef);

    //! return anisotropic permeability directions (used for cloning)
    [[nodiscard]] const std::vector<std::vector<double>>& GetAnisotropicPermeabilityDirections()
        const
    {
      return anisotropic_permeability_directions_;
    }

    //! return scaling coefficients for anisotropic permeability (used for cloning)
    [[nodiscard]] const std::vector<std::vector<double>>& GetAnisotropicPermeabilityNodalCoeffs()
        const
    {
      return anisotropic_permeability_nodal_coeffs_;
    }

    /** \brief set the parameter interface ptr for the solid-poro elements
     *
     *  \param p (in): Parameter list coming from the time integrator.
     *
     *  \author hiermeier
     *  \date 04/16 */
    void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;


    //@}

    /** \brief get access to the parameter interface pointer
     *
     *  \author hiermeier
     *  \date 04/16 */
    Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr() override
    {
      return interface_ptr_;
    }

    /** \brief returns true if the parameter interface is defined and initialized, otherwise false
     *
     *  \author hiermeier
     *  \date 04/16 */
    [[nodiscard]] inline bool IsParamsInterface() const override
    {
      return (not interface_ptr_.is_null());
    }

    /** \brief get access to the interface
     *
     *  \author hiermeier
     *  \date 04/16 */
    [[nodiscard]] inline STR::ELEMENTS::ParamsInterface& ParamsInterface() const
    {
      if (not IsParamsInterface()) dserror("The interface ptr is not set!");
      return *interface_ptr_;
    }

    [[nodiscard]] MAT::StructPoro& StructPoroMaterial(int nummat = 0) const;

    [[nodiscard]] MAT::FluidPoroMultiPhase& FluidPoroMultiMaterial(int nummat = 1) const;

    [[nodiscard]] MAT::So3Material& SolidPoroMaterial(int nummat = 0) const;


    INPAR::PORO::PoroType GetElePoroType() { return porotype_; }
    std::set<INPAR::STR::EleTech> GetEleTech() { return eletech_; }
    INPAR::STR::KinemType GetEleKinematicType() { return kintype_; }
    ::STR::ELEMENTS::EasType GetEAStype() { return eastype_; }
    INPAR::SCATRA::ImplType GetImplType() { return impltype_; }
    void VisNames(std::map<std::string, int>& names) override;
    bool VisData(const std::string& name, std::vector<double>& data) override;

    /*!
    \brief Anything that needs to be done after the standard Evaluate.
    Empty function in the base class that may be overloaded in derived elements.
    */
    virtual int PostEvaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
        CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3);

   private:
    //! discretization type
    DRT::Element::DiscretizationType distype_;

    //! kinematic type
    INPAR::STR::KinemType kintype_;

    //! element techonology
    std::set<INPAR::STR::EleTech> eletech_;

    //! specify EAS type in case  contains EAS
    ::STR::ELEMENTS::EasType eastype_;

    //! scalar transport implementation type (physics)
    INPAR::PORO::PoroType porotype_;

    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype_;

    //! directions for anisotropic permeability
    std::vector<std::vector<double>> anisotropic_permeability_directions_;

    //! scaling coefficients for nodal anisotropic permeability
    std::vector<std::vector<double>> anisotropic_permeability_nodal_coeffs_;

    /** \brief interface ptr
     *
     *  data exchange between the element and the time integrator. */
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! element calculation interface
    std::shared_ptr<SolidEleCalcInterface> solid_interface_ = nullptr;

    //! element calculation interface
    std::shared_ptr<SolidPoroEleCalcInterface> solidporo_interface_ = nullptr;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;


  };  // class Solid

}  // namespace DRT::ELEMENTS

#endif  // SO_PORO_ELE_H
