/*! \file

\brief Declaration of the solid-poro element

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
#include "baci_solid_ele_factory.H"
#include "baci_solid_ele_poro_factory.H"
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
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    DRT::ParObject* Create(const std::vector<char>& data) override;

    [[nodiscard]] std::string Name() const override { return "SolidPoroType"; }

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

    \param id    (in): A globally unique element id
    \param owner (in): owner processor of the element
    */
    SolidPoro(int id, int owner);

    //!@}

    [[nodiscard]] DRT::Element* Clone() const override;

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return SolidPoroType::Instance().UniqueParObjectId();
    };

    [[nodiscard]] int NumLine() const override;

    [[nodiscard]] int NumSurface() const override;

    [[nodiscard]] int NumVolume() const override;

    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

    [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override { return 3; }

    [[nodiscard]] int NumDofPerElement() const override { return 0; }

    void Pack(DRT::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    void SetKinematicType(INPAR::STR::KinemType kintype) { kintype_ = kintype; }

    [[nodiscard]] CORE::FE::CellType Shape() const override { return distype_; };

    [[nodiscard]] DRT::ElementType& ElementType() const override
    {
      return SolidPoroType::Instance();
    }

    bool ReadElement(const std::string& eletype, const std::string& distype,
        DRT::INPUT::LineDefinition* linedef) override;

    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
        CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseVector& elevec2,
        CORE::LINALG::SerialDenseVector& elevec3) override;

    int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

    void ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(
        DRT::INPUT::LineDefinition* linedef);

    void ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(
        DRT::INPUT::LineDefinition* linedef);

    [[nodiscard]] const std::vector<std::vector<double>>& GetAnisotropicPermeabilityDirections()
        const
    {
      return anisotropic_permeability_directions_;
    }

    [[nodiscard]] const std::vector<std::vector<double>>& GetAnisotropicPermeabilityNodalCoeffs()
        const
    {
      return anisotropic_permeability_nodal_coeffs_;
    }

    void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;

    Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr() override
    {
      return interface_ptr_;
    }

    [[nodiscard]] inline bool IsParamsInterface() const override
    {
      return (not interface_ptr_.is_null());
    }

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

   private:
    //! discretization type
    CORE::FE::CellType distype_;

    //! kinematic type
    INPAR::STR::KinemType kintype_;

    //! element techonology
    std::set<INPAR::STR::EleTech> eletech_;

    //! specify EAS type in case contains enhanced assumed strains
    ::STR::ELEMENTS::EasType eastype_;

    //! porosity implementation type (physics)
    INPAR::PORO::PoroType porotype_;

    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype_;

    //! directions for anisotropic permeability
    std::vector<std::vector<double>> anisotropic_permeability_directions_;

    //! scaling coefficients for nodal anisotropic permeability
    std::vector<std::vector<double>> anisotropic_permeability_nodal_coeffs_;

    //! interface pointer for data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! element calculation holding one of the implemented variants
    SolidCalcVariant solid_calc_variant_;

    //! poro element calculation holding one of the implemented variants
    SolidPoroCalcVariant solidporo_calc_variant_;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

  };  // class Solid

}  // namespace DRT::ELEMENTS

#endif  // SO_PORO_ELE_H
