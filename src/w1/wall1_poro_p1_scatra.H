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


#ifndef WALL1_PORO_P1_SCATRA_H_
#define WALL1_PORO_P1_SCATRA_H_

#include "wall1_poro_p1.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media using
    p1 (mixed) approach including scatra functionality

    */
    template <DRT::Element::DiscretizationType distype>
    class Wall1_PoroP1Scatra : public Wall1_PoroP1<distype>
    {
      typedef DRT::ELEMENTS::Wall1_PoroP1<distype> my;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1_PoroP1Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1_PoroP1Scatra(const Wall1_PoroP1Scatra& old);

      /*!
      \brief Destructor

      */
      virtual ~Wall1_PoroP1Scatra(){};

      //@}

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const;

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Unpack(const std::vector<char>& data);

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const;

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          DRT::INPUT::LineDefinition* linedef);

      /// @name params
      /// return SCATRA::ImplType
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

     private:
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;

     protected:
      //! don't want = operator
      Wall1_PoroP1Scatra& operator=(const Wall1_PoroP1Scatra& old);
    };
  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* WALL1_PORO_P1_SCATRA_H_ */
