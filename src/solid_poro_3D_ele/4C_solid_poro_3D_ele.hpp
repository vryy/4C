/*! \file

\brief Declaration of the solid-poro element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_HPP

#include "4C_config.hpp"

#include "4C_inpar_poro.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_lib_element.hpp"
#include "4C_lib_elementtype.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_solid_3D_ele_factory.hpp"
#include "4C_solid_poro_3D_ele_factory.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


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
    void setup_element_definition(
        std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions) override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string elecelltype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

    [[nodiscard]] std::string Name() const override { return "SolidPoroType"; }

    void nodal_block_information(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

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

    [[nodiscard]] int num_dof_per_element() const override { return 0; }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    [[nodiscard]] CORE::FE::CellType Shape() const override { return celltype_; };

    [[nodiscard]] DRT::ElementType& ElementType() const override
    {
      return SolidPoroType::Instance();
    }

    bool ReadElement(const std::string& eletype, const std::string& celltype,
        INPUT::LineDefinition* linedef) override;

    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
        CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseVector& elevec2,
        CORE::LINALG::SerialDenseVector& elevec3) override;

    int evaluate_neumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        CORE::Conditions::Condition& condition, std::vector<int>& lm,
        CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

    void set_params_interface_ptr(const Teuchos::ParameterList& p) override;

    Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr() override
    {
      return interface_ptr_;
    }

    [[nodiscard]] inline bool IsParamsInterface() const override
    {
      return (not interface_ptr_.is_null());
    }

    [[nodiscard]] inline STR::ELEMENTS::ParamsInterface& params_interface() const
    {
      if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
      return *interface_ptr_;
    }

    [[nodiscard]] MAT::StructPoro& StructPoroMaterial(int nummat = 0) const;

    [[nodiscard]] MAT::FluidPoroMultiPhase& fluid_poro_multi_material(int nummat = 1) const;

    [[nodiscard]] MAT::So3Material& SolidPoroMaterial(int nummat = 0) const;

    INPAR::PORO::PoroType GetElePoroType() { return poro_ele_property_.porotype; }

    [[nodiscard]] bool HaveEAS() const
    {
      return solid_ele_property_.element_technology == ElementTechnology::eas_full ||
             solid_ele_property_.element_technology == ElementTechnology::eas_mild;
    }

    INPAR::STR::KinemType GetEleKinematicType() { return solid_ele_property_.kintype; }

    INPAR::SCATRA::ImplType GetImplType() { return poro_ele_property_.impltype; }

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

   private:
    //! cell type
    CORE::FE::CellType celltype_{CORE::FE::CellType::dis_none};

    //! solid element properties
    SolidElementProperties solid_ele_property_{};

    //! additional poro element properties
    SolidPoroElementProperties poro_ele_property_{};

    //! interface pointer for data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! element calculation holding one of the implemented variants
    SolidCalcVariant solid_calc_variant_;

    //! poro element calculation holding one of the implemented variants
    SolidPoroCalcVariant solidporo_calc_variant_;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

  };  // class SolidPoro

}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
