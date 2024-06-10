/*! \file

\brief A 2D shell element with ScaTra functionality

\level 3
*/

#ifndef FOUR_C_SHELL7P_ELE_SCATRA_HPP
#define FOUR_C_SHELL7P_ELE_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_shell7p_ele_calc_interface.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  class Shell7pEleCalcInterface;
  class Shell7pLine;

  class Shell7pScatraType : public Core::Elements::ElementType
  {
   public:
    void setup_element_definition(
        std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions) override;

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

    Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
        const std::string eledistype, const int id, const int owner) override;

    Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

    [[nodiscard]] std::string Name() const override { return "Shell7pScatraType"; }

    int Initialize(Core::FE::Discretization& dis) override;

    void nodal_block_information(
        Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
        Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static Shell7pScatraType& Instance();

   private:
    static Shell7pScatraType instance_;
  };

  class Shell7pScatra : public Core::Elements::Element
  {
    //! @name Friends
    friend class Shell7pScatraType;
    friend class Shell7pLine;

   public:
    //! @name Constructors and destructors related methods
    //! @{
    /*!
     * \brief Standard Constructor
     * @param id (in) : A unique global id
     * @param owner (in) : elements owner
     */
    Shell7pScatra(int id, int owner) : Core::Elements::Element(id, owner){};

    //! copy Constructor
    Shell7pScatra(const Shell7pScatra& other);


    //! copy assignment operator
    Shell7pScatra& operator=(const Shell7pScatra& other);

    //! move constructor
    Shell7pScatra(Shell7pScatra&& other) noexcept = default;

    //! move assignment operator
    Shell7pScatra& operator=(Shell7pScatra&& other) noexcept = default;
    //! @}

    [[nodiscard]] Core::Elements::Element* Clone() const override;

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return Shell7pScatraType::Instance().UniqueParObjectId();
    }

    [[nodiscard]] int NumLine() const override;

    [[nodiscard]] int NumSurface() const override;

    std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

    std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

    [[nodiscard]] int NumDofPerNode(const Core::Nodes::Node& node) const override { return 6; }

    [[nodiscard]] int num_dof_per_element() const override { return 0; }

    void Pack(Core::Communication::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    void Print(std::ostream& os) const override;


    [[nodiscard]] Core::Elements::ElementType& ElementType() const override
    {
      return Shell7pScatraType::Instance();
    }

    [[nodiscard]] Core::FE::CellType Shape() const override { return distype_; };

    bool ReadElement(const std::string& eletype, const std::string& distype,
        Input::LineDefinition* linedef) override;

    //! @name Evaluation
    //! @{
    int Evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
        Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
        Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
        Core::LinAlg::SerialDenseVector& elevec2,
        Core::LinAlg::SerialDenseVector& elevec3) override;

    int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
        Core::Conditions::Condition& condition, std::vector<int>& la,
        Core::LinAlg::SerialDenseVector& elevec1,
        Core::LinAlg::SerialDenseMatrix* elemat1) override;

    //@}

    //! @name Query methods
    //! @{
    [[nodiscard]] inline bool IsParamsInterface() const override
    {
      return (not interface_ptr_.is_null());
    }

    [[nodiscard]] inline STR::ELEMENTS::ParamsInterface& str_params_interface() const
    {
      if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
      return *interface_ptr_;
    }

    void set_params_interface_ptr(const Teuchos::ParameterList& p) override;
    //! @}

    [[nodiscard]] const std::set<Inpar::STR::EleTech>& GetEleTech() const { return eletech_; }

    [[nodiscard]] Teuchos::RCP<Mat::So3Material> SolidMaterial(int nummat = 0) const;

    [[nodiscard]] const Inpar::ScaTra::ImplType& ImplType() const { return impltype_; };

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

    [[nodiscard]] const double& GetThickness() const { return thickness_; }

    [[nodiscard]] const Core::LinAlg::SerialDenseMatrix& GetDirectors() const
    {
      return nodal_directors_;
    }

    inline void set_all_nodal_directors(const Core::LinAlg::SerialDenseMatrix& nodal_directors)
    {
      nodal_directors_ = nodal_directors;
    }

    inline void SetNodalDirector(const int& node_id, const std::vector<double>& director)
    {
      nodal_directors_(node_id, 0) = director[0];
      nodal_directors_(node_id, 1) = director[1];
      nodal_directors_(node_id, 2) = director[2];
    }

   private:
    //! discretization type
    Core::FE::CellType distype_ = Core::FE::CellType::dis_none;

    //! interface ptr, data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_ = Teuchos::null;

    //! element technology
    std::set<Inpar::STR::EleTech> eletech_ = {};

    //! shell thickness in reference frame
    double thickness_ = 0.0;

    //! nodal director
    Core::LinAlg::SerialDenseMatrix nodal_directors_ = {};

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

    //! shell calculation interface
    std::shared_ptr<Shell7pEleCalcInterface> shell_interface_ = nullptr;

    //! scalar transport implementation type (physics)
    Inpar::ScaTra::ImplType impltype_ = Inpar::ScaTra::ImplType::impltype_undefined;
  };

}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
