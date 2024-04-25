/*! \file

\brief Declaration of the solid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_HPP
#define FOUR_C_SOLID_3D_ELE_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_lib_element.hpp"
#include "4C_lib_elementtype.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_solid_3D_ele_calc_eas.hpp"
#include "4C_solid_3D_ele_factory.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

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
        std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions) override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string elecelltype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

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

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    [[nodiscard]] DRT::ElementType& ElementType() const override { return SolidType::Instance(); }

    [[nodiscard]] CORE::FE::CellType Shape() const override { return celltype_; };

    void SetKinematicType(INPAR::STR::KinemType kintype) { solid_ele_property_.kintype = kintype; }

    [[nodiscard]] virtual Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat = 0) const;

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
      if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
      return *interface_ptr_;
    }

    void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;

    [[nodiscard]] bool HaveEAS() const
    {
      return solid_ele_property_.element_technology == ElementTechnology::eas_full ||
             solid_ele_property_.element_technology == ElementTechnology::eas_mild;
    }

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

    /*!
     * @brief Evaluate the Cauchy stress at @p xi with the normal vector @p n in the direction @p
     * dir and compute the linearizations w.r.t. all input parameters (disp, xi, n and dir)
     *
     * @tparam dim (in) : Dimension of the element
     * @param disp (in) : Nodal displacements of the element
     * @param xi (in) : Position of the point in the parameter space of the element
     * @param n (in) : Normal vector
     * @param dir (in) : Direction of the Cauchy stress
     * @return CauchyNDirAndLinearization<dim> (out) : A struct containing the Cauchy stress
     * component and its linearizations
     */
    template <int dim>
    CauchyNDirAndLinearization<dim> GetCauchyNDirAndDerivativesAtXi(const std::vector<double>& disp,
        const CORE::LINALG::Matrix<dim, 1>& xi, const CORE::LINALG::Matrix<dim, 1>& n,
        const CORE::LINALG::Matrix<dim, 1>& dir);

   private:
    //! cell type
    CORE::FE::CellType celltype_ = CORE::FE::CellType::dis_none;

    //! solid element properties
    SolidElementProperties solid_ele_property_{};

    //! interface pointer for data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! element calculation holding one of the implemented variants
    SolidCalcVariant solid_calc_variant_;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

  };  // class Solid

}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
