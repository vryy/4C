/*! \file

\brief Declaration of the solid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

#ifndef BACI_SOLID_ELE_H
#define BACI_SOLID_ELE_H

#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_solid_ele_calc_eas.H"
#include "baci_solid_ele_factory.H"
#include "baci_structure_new_elements_paramsinterface.H"

#include <memory>

namespace MAT
{
  class So3Material;
}
namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidEleCalcInterface;

  class SolidType : public DRT::ElementType
  {
   public:
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    DRT::ParObject* Create(const std::vector<char>& data) override;

    [[nodiscard]] std::string Name() const override { return "SolidType"; }

    void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
        DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static SolidType& Instance();

   private:
    static SolidType instance_;

  };  // class SolidType

  class Solid : public DRT::Element
  {
    friend class SolidType;

   public:
    //! @name Constructors and destructors and related methods
    //!@{

    /*!
    \brief Standard Constructor

    \param id    (in): A globally unique element id
    \param owner (in): owner processor of the element
    */
    Solid(int id, int owner);

    //!@}

    [[nodiscard]] DRT::Element* Clone() const override;

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return SolidType::Instance().UniqueParObjectId();
    };

    void Pack(DRT::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    [[nodiscard]] DRT::ElementType& ElementType() const override { return SolidType::Instance(); }

    [[nodiscard]] DRT::Element::DiscretizationType Shape() const override { return distype_; };

    //! Get kinematic type of element
    [[nodiscard]] INPAR::STR::KinemType KinematicType() const { return kintype_; }

    void SetKinematicType(INPAR::STR::KinemType kintype) { kintype_ = kintype; }

    [[nodiscard]] virtual Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat = 0) const;

    [[nodiscard]] int NumLine() const override;

    [[nodiscard]] int NumSurface() const override;

    [[nodiscard]] int NumVolume() const override;

    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

    std::vector<Teuchos::RCP<DRT::Element>> Volumes() override;

    [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override { return 3; }

    [[nodiscard]] int NumDofPerElement() const override { return 0; }

    bool ReadElement(const std::string& eletype, const std::string& distype,
        DRT::INPUT::LineDefinition* linedef) override;

    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
        CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseVector& elevec2,
        CORE::LINALG::SerialDenseVector& elevec3) override;

    int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

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

    void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;

    void SetEASType(STR::ELEMENTS::EasType type) { eastype_ = type; }

    [[nodiscard]] STR::ELEMENTS::EasType GetEASType() const { return eastype_; }

    [[nodiscard]] INPAR::STR::KinemType GetKinemType() const { return kintype_; }

    [[nodiscard]] const std::set<INPAR::STR::EleTech>& GetEleTech() const { return eletech_; }

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

   private:
    //! discretization type
    DRT::Element::DiscretizationType distype_ = DRT::Element::dis_none;

    //! kinematic type
    INPAR::STR::KinemType kintype_ = INPAR::STR::kinem_vague;

    //! element technology
    std::set<INPAR::STR::EleTech> eletech_;

    //! specify EAS type in case contains enhanced assumed strains
    STR::ELEMENTS::EasType eastype_ = STR::ELEMENTS::EasType::eastype_undefined;

    //! interface pointer for data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! element calculation holding one of the implemented variants
    SolidCalcVariant solid_calc_variant_;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

  };  // class Solid

}  // namespace DRT::ELEMENTS

#endif  // BACI_SOLID_ELE_H
