/*! \file

\brief Declaration of a solid-scatra coupling element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

#ifndef BACI_SOLID_SCATRA_3D_ELE_HPP
#define BACI_SOLID_SCATRA_3D_ELE_HPP

#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_solid_3D_ele_calc_eas.hpp"
#include "baci_solid_scatra_3D_ele_factory.hpp"
#include "baci_structure_new_elements_paramsinterface.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  class So3Material;
}
namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidScatraType : public DRT::ElementType
  {
   public:
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions) override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string elecelltype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

    [[nodiscard]] std::string Name() const override { return "SolidScatraType"; }

    void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
        DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static SolidScatraType& Instance();

   private:
    static SolidScatraType instance_;

  };  // class SolidType

  class SolidScatra : public DRT::Element
  {
    friend class SolidScatraType;

   public:
    SolidScatra(int id, int owner);

    [[nodiscard]] DRT::Element* Clone() const override;

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return SolidScatraType::Instance().UniqueParObjectId();
    }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    [[nodiscard]] DRT::ElementType& ElementType() const override
    {
      return SolidScatraType::Instance();
    }

    [[nodiscard]] CORE::FE::CellType Shape() const override { return celltype_; }

    [[nodiscard]] virtual MAT::So3Material& SolidMaterial(int nummat = 0) const;

    [[nodiscard]] int NumLine() const override;

    [[nodiscard]] int NumSurface() const override;

    [[nodiscard]] int NumVolume() const override;

    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

    [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override { return 3; }

    [[nodiscard]] int NumDofPerElement() const override { return 0; }

    bool ReadElement(const std::string& eletype, const std::string& celltype,
        INPUT::LineDefinition* linedef) override;

    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
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

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

    /// return SCATRA::ImplType
    [[nodiscard]] INPAR::SCATRA::ImplType ImplType() const { return properties_.impltype; }

   private:
    //! cell type
    CORE::FE::CellType celltype_ = CORE::FE::CellType::dis_none;

    //! solid-scatra properties
    SolidScatraElementProperties properties_{};

    //! interface pointer for data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! solid element calculation holding one of the implemented variants
    SolidScatraCalcVariant solid_scatra_calc_variant_;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

  };  // class SolidScatra

}  // namespace DRT::ELEMENTS


BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_SCATRA_3D_ELE_HPP
