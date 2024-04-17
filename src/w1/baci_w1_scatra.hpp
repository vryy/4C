/*----------------------------------------------------------------------------*/
/*! \file
\brief a 2D solid-wall element with ScaTra coupling.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_SCATRA_HPP
#define FOUR_C_W1_SCATRA_HPP

#include "baci_config.hpp"

#include "baci_inpar_scatra.hpp"
#include "baci_w1.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class Wall1ScatraType : public Wall1Type
    {
     public:
      std::string Name() const override { return "Wall1ScatraType"; }

      static Wall1ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Wall1ScatraType instance_;
    };

    class Wall1_Scatra : public Wall1
    {
     public:
      /// @name Constructors and destructors and related methods

      /// Standard Constructor
      Wall1_Scatra(int id,  ///< A unique global id
          int owner);

      /// Copy Constructor
      ///
      /// Makes a deep copy of a Element
      Wall1_Scatra(const Wall1_Scatra& old);

      /// Deep copy this instance of Wall1 and return pointer to the copy
      ///
      /// The Clone() method is used from the virtual base class Element in cases
      /// where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const override;

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
      /// \ref Pack and \ref Unpack are used to communicate this element
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /// Unpack data from a char vector into this class
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;

      //@}

      virtual int NumDofPerNode(
          const unsigned nds, const DRT::Node& node, const std::string disname) const
      {
        if (nds != 0) return 1;
        return Wall1::NumDofPerNode(node);
      };

      /// Print this element
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return Wall1ScatraType::Instance(); }

      //@}

      /// @name Input and Creation
      //@{

      /// Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

      //@}

      /// @name Evaluation
      //@{

      void PreEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       ///< location array for de-assembly
      );

      /// Evaluate an element
      ///
      /// Evaluate Wall1 element stiffness, mass, internal forces etc
      ///
      /// \return 0 if successful, negative otherwise
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  ///< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   ///< vector to be filled by element
          ) override;

      /// @name params
      /// return SCATRA::ImplType
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

     private:
      //@{
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      Wall1_Scatra& operator=(const Wall1_Scatra& old);

      int MyEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  ///< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   ///< vector to be filled by element
      );
    };


  }  // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
