/*! \file

\brief The 7-parameter shell element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

In addition, it contains the interface between element call and gauss point
loop for the evaluation.

\level 3
*/


#ifndef FOUR_C_SHELL7P_ELE_HPP
#define FOUR_C_SHELL7P_ELE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_shell7p_ele_calc_interface.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Mat
{
  class So3Material;
}  // namespace Mat

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declaration
    class Shell7pEleCalcInterface;
    class Shell7pLine;


    class Shell7pType : public Core::Elements::ElementType
    {
     public:
      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      int Initialize(Core::FE::Discretization& dis) override;

      [[nodiscard]] std::string Name() const override { return "Shell7pType"; }

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      static Shell7pType& Instance();

     private:
      static Shell7pType instance_;

    };  // class Shell7pType

    /*!
      *  \brief  The 7-parameter shell element
      *
      *  A nonlinear shell element based on the 7-parameter shell formulation with relative
      *  displacements as parameterization of rotations. The formulation is based on [1]. This
      *  formulation is adopted by a scaled-director approach, see [2], to remove the
      *  extra-ill conditioning of the resulting matrices that appears with three-dimensional shell
      *  formulations as compared to formulations that neglect the thickness change of the shell.
      *
      *  The shell element can be used with tri3, tri6, quad4, quad8 and quad9 elements. The
      *  7-parameter formulation has three translational and three relative-translational degrees
      *  of freedom with an additional linear thickness change introduced via a hybrid formulation.
      *
      *  The current formulation is able to remedy different types of locking for quadrilateral
      *  elements due to the EAS and ANS method. The EAS input parameters corresponds to membrane
      *  locking, bending locking, locking due to thickness changes, locking due to transverse
      *  shear strain (constant and linear). The individual mechanism to remedy locking can be
      *  turned off/on separately. The combined EAS/ANS element can be practically considered
      *  locking free. Currently, the number of gauss points in thickness direction are prescribed
      *  by two integration points. Otherwise, the shell element would suffer from a
      *  nonlinear poisson stiffening affect, because the implemented formulation assumes that
      *  the strains are linear in thickness direction.
      *
      *  The element supports the use of 3D Solid Materials, however only St. Venant Kirchhoff,
      *  Coupled NeoHooke, IsoNeoHooke were tested now.
      *
      *  In this shell formulation the strains and stresses are not sorted as usual in the
      *  voigt notation, instead they are sorted as follows: alpha = {alpha_11 alpha_12 alpha_13
      *  alpha_22 alpha_23 alpha_33 beta_11 beta_12 beta_13 beta_22 beta_23 beta_33 }
      *

      *  References:
      *  [1] Theorie und Numerik einer dreidimensionalen Schalenformulierung (1999), M.Bischoff,
      *      PhD-Thesis
      *  [2] Effiziente Loesungsstrategien in der nichtlinearen Schalenmechanik (2004), M.Gee,
      *      PhD-Thesis
      *
      */
    class Shell7p : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class Shell7pType;
      friend class Shell7pLine;

      //! @name Constructors and destructors related methods
      //! @{
      /*!
      \brief Standard Constructor

      @param id    (in): A globally unique element id
      @param owner (in): owner processor of the element
      */
      Shell7p(int id, int owner) : Core::Elements::Element(id, owner){};


      //! copy constructor
      Shell7p(const Discret::ELEMENTS::Shell7p& other);

      //! copy assignment operator
      Shell7p& operator=(const Shell7p& other);

      //! move constructor
      Shell7p(Shell7p&& other) noexcept = default;

      //! move assignment operator
      Shell7p& operator=(Shell7p&& other) noexcept = default;
      //! @}

      [[nodiscard]] Core::Elements::Element* Clone() const override;

      [[nodiscard]] int UniqueParObjectId() const override
      {
        return Shell7pType::Instance().UniqueParObjectId();
      };

      void pack(Core::Communication::PackBuffer& data) const override;

      void unpack(const std::vector<char>& data) override;

      [[nodiscard]] Core::Elements::ElementType& ElementType() const override
      {
        return Shell7pType::Instance();
      }

      [[nodiscard]] Core::FE::CellType Shape() const override { return distype_; };

      [[nodiscard]] int NumLine() const override;

      [[nodiscard]] int NumSurface() const override;

      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

      [[nodiscard]] int NumDofPerNode(const Core::Nodes::Node& node) const override { return 6; }

      [[nodiscard]] int num_dof_per_element() const override { return 0; }

      //! @name Evaluation
      //! @{
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& dof_index_array, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& dof_index_array,
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

      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          Input::LineDefinition* linedef) override;

      [[nodiscard]] const std::set<Inpar::STR::EleTech>& GetEleTech() const { return eletech_; }

      [[nodiscard]] Teuchos::RCP<Mat::So3Material> SolidMaterial(int nummat = 0) const;

      void Print(std::ostream& os) const override;

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

      //! nodal director matrix
      Core::LinAlg::SerialDenseMatrix nodal_directors_{};

      //! flag, whether the post setup of materials is already called
      bool material_post_setup_ = false;

      //! shell calculation interface
      std::shared_ptr<Shell7pEleCalcInterface> shell_interface_ = nullptr;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif