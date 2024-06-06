/*! \file

\brief Declaration of a solid-scatra coupling element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_solid_3D_ele_calc_eas.hpp"
#include "4C_solid_scatra_3D_ele_calc_lib_nitsche.hpp"
#include "4C_solid_scatra_3D_ele_factory.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class So3Material;
}
namespace Discret::ELEMENTS
{
  // forward declaration
  class SolidScatraType : public Core::Elements::ElementType
  {
   public:
    void setup_element_definition(
        std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions) override;

    Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
        const std::string elecelltype, const int id, const int owner) override;

    Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

    [[nodiscard]] std::string Name() const override { return "SolidScatraType"; }

    void nodal_block_information(
        Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
        Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static SolidScatraType& Instance();

   private:
    static SolidScatraType instance_;

  };  // class SolidType

  class SolidScatra : public Core::Elements::Element
  {
    friend class SolidScatraType;

   public:
    SolidScatra(int id, int owner);

    [[nodiscard]] Core::Elements::Element* Clone() const override;

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return SolidScatraType::Instance().UniqueParObjectId();
    }

    void Pack(Core::Communication::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    [[nodiscard]] Core::Elements::ElementType& ElementType() const override
    {
      return SolidScatraType::Instance();
    }

    [[nodiscard]] Core::FE::CellType Shape() const override { return celltype_; }

    [[nodiscard]] virtual Mat::So3Material& SolidMaterial(int nummat = 0) const;

    [[nodiscard]] int NumLine() const override;

    [[nodiscard]] int NumSurface() const override;

    [[nodiscard]] int NumVolume() const override;

    std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

    std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

    [[nodiscard]] int NumDofPerNode(const Core::Nodes::Node& node) const override { return 3; }

    [[nodiscard]] int num_dof_per_element() const override { return 0; }

    bool ReadElement(const std::string& eletype, const std::string& celltype,
        Input::LineDefinition* linedef) override;

    int Evaluate(Teuchos::ParameterList& params, Discret::Discretization& discretization,
        Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
        Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
        Core::LinAlg::SerialDenseVector& elevec2,
        Core::LinAlg::SerialDenseVector& elevec3) override;

    int evaluate_neumann(Teuchos::ParameterList& params, Discret::Discretization& discretization,
        Core::Conditions::Condition& condition, std::vector<int>& lm,
        Core::LinAlg::SerialDenseVector& elevec1,
        Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

    Teuchos::RCP<Core::Elements::ParamsInterface> ParamsInterfacePtr() override
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

    void set_params_interface_ptr(const Teuchos::ParameterList& p) override;

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

    /// return ScaTra::ImplType
    [[nodiscard]] Inpar::ScaTra::ImplType ImplType() const { return properties_.impltype; }

    /*!
     * @brief Returns the Cauchy stress in the direction @p dir at @p xi with normal @p n
     *
     * @param disp Nodal displacements of the element
     * @param scalars Scalars at the nodes of the element
     * @param xi
     * @param n
     * @param dir
     * @param linearizations [in/out] : Struct holding the linearizations that are possible for
     * evaluation
     * @return double
     *
     * @note @p scalars is an optional since it might not be set in the very initial call of the
     * stucture. Once the structure does not evaluate itself after setup, this optional parameter
     * can be made mandatory.
     */
    double GetCauchyNDirAtXi(const std::vector<double>& disp,
        const std::optional<std::vector<double>>& scalars, const Core::LinAlg::Matrix<3, 1>& xi,
        const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
        SolidScatraCauchyNDirLinearizations<3>& linearizations);

   private:
    //! cell type
    Core::FE::CellType celltype_ = Core::FE::CellType::dis_none;

    //! solid-scatra properties
    SolidScatraElementProperties properties_{};

    //! interface pointer for data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! solid element calculation holding one of the implemented variants
    SolidScatraCalcVariant solid_scatra_calc_variant_;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

  };  // class SolidScatra

}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
