/*! \file

\brief A 2D shell element with ScaTra functionality

\level 3
*/

#ifndef FOUR_C_SHELL7P_ELE_SCATRA_HPP
#define FOUR_C_SHELL7P_ELE_SCATRA_HPP

#include "baci_config.hpp"

#include "baci_inpar_scatra.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_lib_utils.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_shell7p_ele_calc_interface.hpp"
#include "baci_structure_new_elements_paramsinterface.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <memory>

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  class Shell7pEleCalcInterface;
  class Shell7pLine;

  class Shell7pScatraType : public DRT::ElementType
  {
   public:
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions) override;

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
        const int id, const int owner) override;

    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    [[nodiscard]] std::string Name() const override { return "Shell7pScatraType"; }

    int Initialize(DRT::Discretization& dis) override;

    void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
        DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static Shell7pScatraType& Instance();

   private:
    static Shell7pScatraType instance_;
  };

  class Shell7pScatra : public DRT::Element
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
    Shell7pScatra(int id, int owner) : DRT::Element(id, owner){};

    //! copy Constructor
    Shell7pScatra(const Shell7pScatra& other);


    //! copy assignment operator
    Shell7pScatra& operator=(const Shell7pScatra& other);

    //! move constructor
    Shell7pScatra(Shell7pScatra&& other) noexcept = default;

    //! move assignment operator
    Shell7pScatra& operator=(Shell7pScatra&& other) noexcept = default;
    //! @}

    [[nodiscard]] DRT::Element* Clone() const override;

    [[nodiscard]] int UniqueParObjectId() const override
    {
      return Shell7pScatraType::Instance().UniqueParObjectId();
    }

    [[nodiscard]] int NumLine() const override;

    [[nodiscard]] int NumSurface() const override;

    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

    [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override { return 6; }

    [[nodiscard]] int NumDofPerElement() const override { return 0; }

    void Pack(CORE::COMM::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    void Print(std::ostream& os) const override;


    [[nodiscard]] DRT::ElementType& ElementType() const override
    {
      return Shell7pScatraType::Instance();
    }

    [[nodiscard]] CORE::FE::CellType Shape() const override { return distype_; };

    bool ReadElement(const std::string& eletype, const std::string& distype,
        INPUT::LineDefinition* linedef) override;

    //! @name Evaluation
    //! @{
    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
        CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseVector& elevec2,
        CORE::LINALG::SerialDenseVector& elevec3) override;

    int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& la, CORE::LINALG::SerialDenseVector& elevec1,
        CORE::LINALG::SerialDenseMatrix* elemat1) override;

    //@}

    //! @name Query methods
    //! @{
    [[nodiscard]] inline bool IsParamsInterface() const override
    {
      return (not interface_ptr_.is_null());
    }

    [[nodiscard]] inline STR::ELEMENTS::ParamsInterface& StrParamsInterface() const
    {
      if (not IsParamsInterface()) dserror("The interface ptr is not set!");
      return *interface_ptr_;
    }

    void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;
    //! @}

    [[nodiscard]] const std::set<INPAR::STR::EleTech>& GetEleTech() const { return eletech_; }

    [[nodiscard]] Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat = 0) const;

    [[nodiscard]] const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data) override;

    [[nodiscard]] const double& GetThickness() const { return thickness_; }

    [[nodiscard]] const CORE::LINALG::SerialDenseMatrix& GetDirectors() const
    {
      return nodal_directors_;
    }

    inline void SetAllNodalDirectors(const CORE::LINALG::SerialDenseMatrix& nodal_directors)
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
    CORE::FE::CellType distype_ = CORE::FE::CellType::dis_none;

    //! interface ptr, data exchange between the element and the time integrator.
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_ = Teuchos::null;

    //! element technology
    std::set<INPAR::STR::EleTech> eletech_ = {};

    //! shell thickness in reference frame
    double thickness_ = 0.0;

    //! nodal director
    CORE::LINALG::SerialDenseMatrix nodal_directors_ = {};

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_ = false;

    //! shell calculation interface
    std::shared_ptr<Shell7pEleCalcInterface> shell_interface_ = nullptr;

    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype_ = INPAR::SCATRA::ImplType::impltype_undefined;
  };

}  // namespace DRT::ELEMENTS


BACI_NAMESPACE_CLOSE

#endif  // BACI_SHELL7P_ELE_SCATRA_H
