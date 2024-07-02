/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element including scatra functionality

 \level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SO3_PORO_SCATRA_HPP
#define FOUR_C_SO3_PORO_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_so3_poro.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 3 dimensional solid element with modifications for porous media
    including scatra functionality

    A structural 3 dimensional solid displacement element for large deformations
    and (near)-incompressibility.

    */
    template <class So3Ele, Core::FE::CellType distype>
    class So3PoroScatra : public So3Poro<So3Ele, distype>
    {
      using my = So3Poro<So3Ele, distype>;
      using my::numnod_;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So3PoroScatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So3PoroScatra(const So3PoroScatra& old);

      //@}

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override;

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //! @name Access methods

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate element stiffness, mass, internal forces, etc.

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,  //!< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la,  //!< location array for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  //!< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  //!< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  //!< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   //!< vector to be filled by element
          ) override;

      void pre_evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,   //!< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la  //!< location array for de-assembly
          ) override;

      //!@}

      //! @name Input and Creation
      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          Input::LineDefinition* linedef) override;

      /// @name params
      /// return ScaTra::ImplType
      const Inpar::ScaTra::ImplType& ImplType() const { return impltype_; };

     private:
      //! scalar transport implementation type (physics)
      Inpar::ScaTra::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      So3PoroScatra& operator=(const So3PoroScatra& old);

    };  // class So3_Poro_Scatra
  }     // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
