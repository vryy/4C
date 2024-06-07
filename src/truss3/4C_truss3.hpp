/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/
#ifndef FOUR_C_TRUSS3_HPP
#define FOUR_C_TRUSS3_HPP

#include "4C_config.hpp"

#include "4C_fem_general_elementtype.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN
using FAD = Sacado::Fad::DFad<double>;

namespace STR
{
  namespace ELEMENTS
  {
    class ParamsInterface;
  }
}  // namespace STR

namespace Discret
{
  namespace ELEMENTS
  {
    class Truss3Type : public Core::Elements::ElementType
    {
     public:
      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Discret::Discretization& dis) override;

      static Truss3Type& Instance();

      std::string Name() const override { return "Truss3Type"; }

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static Truss3Type instance_;
    };

    /*!
     \brief three dimensional total Lagrange truss element

     */
    class Truss3 : public Core::Elements::Element
    {
     public:
      //! @name Friends
      friend class Truss3Type;

      /*!
       \brief Standard Constructor

       \param id    (in): A globally unique element id
       \param owner (in): owner processor of the element
       */
      Truss3(int id, int owner);

      /*!
       \brief Copy Constructor

       Makes a deep copy of a Element
       */
      Truss3(const Truss3& old);

      Core::Elements::Element* Clone() const override;

      //! prepare elemental specific geometric values
      //! \param[in] ele_state              elemental states (depending on the instantiated element)
      //! \param[out] curr_nodal_coords     nodal coordinates
      //! \param[out] dcurr_nodal_coords_du deriv. of nodal coordinates w.r.t. global displacement
      //! \param[out] dN_dx                 derivative of shape functions
      void prep_calc_internal_force_stiff_tot_lag(
          const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::Matrix<6, 1>& curr_nodal_coords,
          Core::LinAlg::Matrix<6, 6>& dcurr_nodal_coords_du, Core::LinAlg::Matrix<6, 1>& dN_dx);

      //! \brief calculate internal force vector and stiffness matrix based on absolute nodal
      //! positions (using kinematic type tr3_totlag)
      //!
      //! \param[in] ele_state    elemental states (depending on the instantiated element)
      //! \param[out] forcevec    element force vector
      //! \param[out] stiffmat    element stiffness matrix
      virtual void calc_internal_force_stiff_tot_lag(
          const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::SerialDenseVector& forcevec, Core::LinAlg::SerialDenseMatrix& stiffmat);

      //! calcluate stresses at Gauss point
      //! \param[in] params      parameter list
      //! \param[in] ele_state   elemental states (depending on the instantiated element)
      virtual void CalcGPStresses(Teuchos::ParameterList& params,
          const std::map<std::string, std::vector<double>>& ele_state);

      Core::Elements::ElementType& ElementType() const override { return Truss3Type::Instance(); }

      int Evaluate(Teuchos::ParameterList& params, Discret::Discretization& discretization,
          LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      int evaluate_neumann(Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //! get internal (elastic) energy of element
      double GetInternalEnergy() const { return eint_; };

      inline bool IsParamsInterface() const override { return (not interface_ptr_.is_null()); }

      //! cross section area
      double CrossSection() const { return crosssec_; }

      //! Return the current length of the truss from @p curr_nodal_coords
      double CurrLength(const Core::LinAlg::Matrix<6, 1>& curr_nodal_coords) const
      {
        return curr_nodal_coords.Norm2() * M_SQRT1_2;
      }

      //! Return the squared value of the current length of the truss from @p curr_nodal_coords
      double CurrLength2(const Core::LinAlg::Matrix<6, 1>& curr_nodal_coords) const
      {
        return CurrLength(curr_nodal_coords) * CurrLength(curr_nodal_coords);
      }

      //! derivative of current length w.r.t. nodal coordinate (entry @p col) from @p
      //! curr_nodal_coords
      double dCurrLengthdu(const Core::LinAlg::Matrix<6, 1>& curr_nodal_coords, const int col) const
      {
        return curr_nodal_coords(col) / curr_nodal_coords.Norm2() * M_SQRT1_2;
      }

      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      // TODO: remove once truss3 element is fixed and no longer expects more dofs (6) than it can
      // inherently handle (3)...
      void LocationVector(
          const Discretization& dis, LocationArray& la, bool doDirichlet) const override;

      int num_dof_per_element() const override { return 0; }

      int NumDofPerNode(const Core::Nodes::Node& node) const override
      {
        /*note: this is not necessarily the number of DOF assigned to this node by the
         *discretization finally, but only the number of DOF requested for this node by this
         *element; the discretization will finally assign the maximal number of DOF to this node
         *requested by any element connected to this node*/
        return 3;
      }

      int NumLine() const override { return 1; }

      void Pack(Core::Communication::PackBuffer& data) const override;

      Teuchos::RCP<Core::Elements::ParamsInterface> ParamsInterfacePtr() override;

      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //! scale truss reference length
      void scale_reference_length(double scalefac);

      //! set cross section area of this element
      void SetCrossSec(const double& crosssec);

      void set_params_interface_ptr(const Teuchos::ParameterList& p) override;

      //! \brief sets reference coordinates X_ and reference length lrefe_ for elements added to
      //! the discretization
      //!
      //! \param xrefe     nodal coordinates in reference frame
      void set_up_reference_geometry(const std::vector<double>& xrefe);

      Core::FE::CellType Shape() const override;

      int UniqueParObjectId() const override { return Truss3Type::Instance().UniqueParObjectId(); }

      void Unpack(const std::vector<char>& data) override;

      //! coordinates of nodes in reference configuration
      const Core::LinAlg::Matrix<6, 1>& X() const { return x_; }

     protected:
      //! kind of integration to be performed
      enum IntegrationType
      {
        gaussexactintegration,
        gaussunderintegration,
        lobattointegration
      };

      //! get access to the parameter interface
      inline STR::ELEMENTS::ParamsInterface& params_interface()
      {
        if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
        return *interface_ptr_;
      }

      //! extract elemental quantities from nodal quantities
      //!
      //! \param[in] la              location array
      //! \param[in] discretization  discretization
      //! \param[in] params          parameter list
      //! \param[out] ele_state      elemental states (depending on the instantiated element)
      virtual void extract_elemental_variables(LocationArray& la,
          const Discret::Discretization& discretization, const Teuchos::ParameterList& params,
          std::map<std::string, std::vector<double>>& ele_state);

      //! determine Gauss rule from required type of integration
      Core::FE::GaussRule1D my_gauss_rule(int nnode, IntegrationType integrationtype);

      //! calculation of elastic energy
      //!
      //! \param ele_state [in]   elemental states (depending on the instantiated element)
      //! \param params
      //! \param intenergy
      virtual void energy(const std::map<std::string, std::vector<double>>& ele_state,
          Teuchos::ParameterList& params, Core::LinAlg::SerialDenseVector& intenergy);

      //! cross section area
      double crosssec_;

      //! internal energy of element
      double eint_;

      //! length in reference configuration
      double lrefe_;

      //! gaussrule_ will be initialized automatically to a 2 point integration rule
      Core::FE::GaussRule1D gaussrule_;

     private:
      //! possible kinematic types
      enum class KinematicType
      {
        tr3_totlag,
        tr3_engstrain
      };

      //! lump mass matrix
      void lump_mass(Core::LinAlg::SerialDenseMatrix* emass);

      //! calculation of nonlinear stiffness and mass matrix switching between total lagrange
      //! and enginerring strains
      //!
      //! \param[in] ele_state     elemental states (depending on the instantiated element)
      //! \param[out] stiffmatrix  elemental sitffness matrix
      //! \param[out] massmatrix   elemental mass matrix
      //! \param[out] force        elemental force vector
      void nln_stiff_mass(const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::SerialDenseMatrix* stiffmatrix, Core::LinAlg::SerialDenseMatrix* massmatrix,
          Core::LinAlg::SerialDenseVector* force);

      //! \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain
      //! measure.
      //!
      //! Unlike the fully nonlinear implementation of #t3_nlnstiffmass_totlag, this function uses
      //! \f$\varepsilon=\Delta d / d\f$ as strain measure.
      //!
      //! \param[in] ele_state            elemental states (depending on the instantiated element)
      //! \param[out] DummyStiffMatrix    elemental sitffness matrix
      //! \param[out] massmatrix          elemental mass matrix
      //! \param[out] DummyForce          elemental force vector
      void nln_stiff_mass_eng_str(const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::SerialDenseMatrix& DummyStiffMatrix,
          Core::LinAlg::SerialDenseMatrix* massmatrix, Core::LinAlg::SerialDenseVector& DummyForce);

      //! calculation of nonlinear stiffness and mass matrix
      //!
      //! \param[in] ele_state           elemental states (depending on the instantiated element)
      //! \param[out] DummyStiffMatrix   elemental sitffness matrix
      //! \param[out] massmatrix         elemental mass matrix
      //! \param[out] DummyForce         elemental force vector
      void nln_stiff_mass_tot_lag(const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::SerialDenseMatrix& DummyStiffMatrix,
          Core::LinAlg::SerialDenseMatrix* massmatrix, Core::LinAlg::SerialDenseVector& DummyForce);

      //! reference tangent position
      Core::LinAlg::Matrix<1, 3> diff_disp_ref_;

      //!  data exchange between the element and the time integrator.
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      //! variable saving whether element has already been initialized (then isinit_ == true)
      bool isinit_;

      //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
      std::vector<double> jacobimass_;

      //! vector holding value of Jacobi determinant jacobi at nodes
      std::vector<double> jacobinode_;

      //! Kinematic type
      KinematicType kintype_;

      //! material type
      int material_;

      //! reference node position
      Core::LinAlg::Matrix<6, 1> x_;

      // don't want = operator
      Truss3& operator=(const Truss3& old);
    };

    // << operator
    std::ostream& operator<<(std::ostream& os, const Core::Elements::Element& ele);

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
