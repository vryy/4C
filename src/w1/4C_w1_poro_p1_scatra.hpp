/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D wall element for solid-part of porous medium using p1 (mixed) approach including scatra
functionality

\level 2


*/
/*---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*!

 \brief a 2D wall element for solid-part of porous medium using p1 (mixed) approach including scatra
functionality

 \level 2

*----------------------------------------------------------------------*/


#ifndef FOUR_C_W1_PORO_P1_SCATRA_HPP
#define FOUR_C_W1_PORO_P1_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_w1_poro_p1.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media using
    p1 (mixed) approach including scatra functionality

    */
    template <Core::FE::CellType distype>
    class Wall1PoroP1Scatra : public Wall1PoroP1<distype>
    {
      typedef Discret::ELEMENTS::Wall1PoroP1<distype> my;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1PoroP1Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1PoroP1Scatra(const Wall1PoroP1Scatra& old);


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

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override;

      //@}

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

     protected:
      //! don't want = operator
      Wall1PoroP1Scatra& operator=(const Wall1PoroP1Scatra& old);
    };
  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
