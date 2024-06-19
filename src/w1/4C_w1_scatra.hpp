/*----------------------------------------------------------------------------*/
/*! \file
\brief a 2D solid-wall element with ScaTra coupling.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_SCATRA_HPP
#define FOUR_C_W1_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_w1.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Wall1ScatraType : public Wall1Type
    {
     public:
      std::string Name() const override { return "Wall1ScatraType"; }

      static Wall1ScatraType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static Wall1ScatraType instance_;
    };

    class Wall1Scatra : public Wall1
    {
     public:
      /// @name Constructors and destructors and related methods

      /// Standard Constructor
      Wall1Scatra(int id,  ///< A unique global id
          int owner);

      /// Copy Constructor
      ///
      /// Makes a deep copy of a Element
      Wall1Scatra(const Wall1Scatra& old);

      /// Deep copy this instance of Wall1 and return pointer to the copy
      ///
      /// The Clone() method is used from the virtual base class Element in cases
      /// where the type of the derived class is unknown and a copy-ctor is needed
      Core::Elements::Element* Clone() const override;

      /// Return unique ParObject id
      ///
      /// every class implementing ParObject needs a unique id defined at the
      /// top of this file.
      int UniqueParObjectId() const override
      {
        return Wall1ScatraType::Instance().UniqueParObjectId();
      }

      /// Pack this class so it can be communicated
      ///
      /// \ref pack and \ref unpack are used to communicate this element
      void pack(Core::Communication::PackBuffer& data) const override;

      /// Unpack data from a char vector into this class
      ///
      /// \ref pack and \ref unpack are used to communicate this element
      void unpack(const std::vector<char>& data) override;

      //@}

      virtual int NumDofPerNode(
          const unsigned nds, const Core::Nodes::Node& node, const std::string disname) const
      {
        if (nds != 0) return 1;
        return Wall1::NumDofPerNode(node);
      };

      /// Print this element
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override
      {
        return Wall1ScatraType::Instance();
      }

      //@}

      /// @name Input and Creation
      //@{

      /// Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      /// @name Evaluation
      //@{

      void pre_evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,   ///< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la  ///< location array for de-assembly
      );

      /// Evaluate an element
      ///
      /// Evaluate Wall1 element stiffness, mass, internal forces etc
      ///
      /// \return 0 if successful, negative otherwise
      int evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,  ///< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la,  ///< location array for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  ///< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   ///< vector to be filled by element
          ) override;

      /// @name params
      /// return ScaTra::ImplType
      const Inpar::ScaTra::ImplType& ImplType() const { return impltype_; };

     private:
      //@{
      //! scalar transport implementation type (physics)
      Inpar::ScaTra::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      Wall1Scatra& operator=(const Wall1Scatra& old);

      int my_evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,  ///< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la,  ///< location array for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  ///< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   ///< vector to be filled by element
      );
    };


  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
