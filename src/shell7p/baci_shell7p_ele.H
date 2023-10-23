/*! \file

\brief The 7-parameter shell element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

In addition, it contains the interface between element call and gauss point
loop for the evaluation.

\level 3
*/


#ifndef BACI_SHELL7P_ELE_H
#define BACI_SHELL7P_ELE_H

#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_utils.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_shell7p_ele_calc_interface.H"
#include "baci_structure_new_elements_paramsinterface.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <memory>


// forward declaration
namespace MAT
{
  class So3Material;
}  // namespace MAT

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declaration
    class Shell7pEleCalcInterface;
    class Shell7pLine;


    class Shell7pType : public DRT::ElementType
    {
     public:
      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      DRT::ParObject* Create(const std::vector<char>& data) override;

      int Initialize(DRT::Discretization& dis) override;

      [[nodiscard]] std::string Name() const override { return "Shell7pType"; }

      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

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
    class Shell7p : public DRT::Element
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
      Shell7p(int id, int owner)
          : DRT::Element(id, owner),
            distype_(DRT::Element::DiscretizationType::dis_none),
            interface_ptr_(Teuchos::null),
            thickness_(0.0),
            material_post_setup_(false)
      {
        nodal_directors_(0, 0);
        shell_interface_ = nullptr;
      };

      //! destructor
      ~Shell7p() override = default;

      //! copy constructor
      Shell7p(const DRT::ELEMENTS::Shell7p& other);

      //! copy assignment operator
      Shell7p& operator=(const Shell7p& other);

      //! move constructor
      Shell7p(Shell7p&& other) noexcept = default;

      //! move assignment operator
      Shell7p& operator=(Shell7p&& other) noexcept = default;
      //! @}

      [[nodiscard]] DRT::Element* Clone() const override;

      [[nodiscard]] int UniqueParObjectId() const override
      {
        return Shell7pType::Instance().UniqueParObjectId();
      };

      void Pack(DRT::PackBuffer& data) const override;

      void Unpack(const std::vector<char>& data) override;

      [[nodiscard]] DRT::ElementType& ElementType() const override
      {
        return Shell7pType::Instance();
      }

      [[nodiscard]] DRT::Element::DiscretizationType Shape() const override { return distype_; };

      [[nodiscard]] int NumLine() const override;

      [[nodiscard]] int NumSurface() const override;

      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override { return 6; }

      [[nodiscard]] int NumDofPerElement() const override { return 0; }

      //! @name Evaluation
      //! @{
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& dof_index_array,
          CORE::LINALG::SerialDenseVector& elevec1,
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

      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          DRT::INPUT::LineDefinition* linedef) override;

      [[nodiscard]] const std::set<INPAR::STR::EleTech>& GetEleTech() const { return eletech_; }

      [[nodiscard]] Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat = 0) const;

      void Print(std::ostream& os) const override;

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
      DRT::Element::DiscretizationType distype_;

      //! interface ptr, data exchange between the element and the time integrator.
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      //! element technology
      std::set<INPAR::STR::EleTech> eletech_;

      //! shell thickness in reference frame
      double thickness_;

      //! nodal director matrix
      CORE::LINALG::SerialDenseMatrix nodal_directors_;

      //! flag, whether the post setup of materials is already called
      bool material_post_setup_ = false;

      //! shell calculation interface
      std::shared_ptr<Shell7pEleCalcInterface> shell_interface_ = nullptr;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif