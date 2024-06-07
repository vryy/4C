/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly Compressible fluid element based on the HDG method

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_HDG_WEAK_COMP_HPP
#define FOUR_C_FLUID_ELE_HDG_WEAK_COMP_HPP

#include "4C_config.hpp"

#include "4C_fem_general_dg_element.hpp"
#include "4C_fem_general_utils_polynomial.hpp"
#include "4C_fluid_ele.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Input
{
  class LineDefinition;
}

namespace Discret
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidHDGWeakCompType : public FluidType
    {
     public:
      std::string Name() const override { return "FluidHDGWeakCompType"; }

      static FluidHDGWeakCompType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      virtual void ComputeNullSpace(Discret::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static FluidHDGWeakCompType instance_;
    };


    /*!
    \brief HDG weakly compressible fluid element
    */
    class FluidHDGWeakComp : public Fluid, public Core::Elements::DgElement
    {
     public:
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidHDGWeakComp(int id,  ///< A unique global id
          int owner             ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidHDGWeakComp(const FluidHDGWeakComp& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        return FluidHDGWeakCompType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Unpack(const std::vector<char>& data) override;

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Access methods
      /*!
      \brief Get number of degrees of freedom per node, zero for the primary
      and the secondary dof set and equal to the given number for the tertiary dof set
      */
      int NumDofPerNode(const Core::Nodes::Node&) const override { return 0; }

      /*!
      \brief Returns the number of dofs per node for the ALE displacements
       */
      int num_dof_per_node_auxiliary() const override
      {
        return (Core::FE::getDimension(distype_) + 0);
        // maybe + 1 to mimic the standard fluid dofset with velocity + pressure
        // this is maybe needed for the fluid-ALE coupling
      }

      /*!
      \brief Get number of degrees of freedom per face

      */
      int num_dof_per_face(const unsigned face) const override
      {
        return (1 + Core::FE::getDimension(distype_)) * NumDofPerComponent(face);
      }

      /*!
      \brief Get number of dofs per component per face
      */
      int NumDofPerComponent(const unsigned face) const override
      {
        return Core::FE::getBasisSize(
            Core::FE::getEleFaceShapeType(distype_), this->Degree(), completepol_);
      }

      /*!
      \brief Get number of degrees of freedom per element, zero for the primary dof set
      and equal to the given number for the secondary dof set
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Returns the number of dofs per element for the interior variables
       */
      int num_dof_per_element_auxiliary() const override
      {
        const int nsd_ = Core::FE::getDimension(distype_);
        const int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;
        return (msd_ + nsd_ + 1.0) * Core::FE::getBasisSize(distype_, degree_, completepol_);
      }

      /*!
       \brief Returns the degree of the element
       */
      int Degree() const override { return degree_; }

      /*!
       \brief Returns the degree of the element
       */
      int uses_complete_polynomial_space() const { return completepol_; }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, that is, call the element routines to evaluate fluid
      element matrices and vectors or evaluate errors, statistics or updates etc. directly.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, Discret::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      //@}

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override
      {
        return FluidHDGWeakCompType::Instance();
      }

     private:
      // don't want = operator
      FluidHDGWeakComp& operator=(const FluidHDGWeakComp& old);

      // stores the degree of the element
      unsigned char degree_;

      // stores the polynomial type (tensor product or complete polynomial)
      bool completepol_;
    };  // class Fluid


  }  // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
