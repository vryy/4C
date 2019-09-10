/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element including scatra functionality

 \level 2

 \maintainer Christoph Schmidt
             schmidt@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 15251
 *----------------------------------------------------------------------*/


#ifndef SO3_PORO_SCATRA_H_
#define SO3_PORO_SCATRA_H_

#include "so3_poro.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 3 dimensional solid element with modifications for porous media
    including scatra functionality

    A structural 3 dimensional solid displacement element for large deformations
    and (near)-incompressibility.

    */
    template <class so3_ele, DRT::Element::DiscretizationType distype>
    class So3_Poro_Scatra : public So3_Poro<so3_ele, distype>
    {
      typedef So3_Poro<so3_ele, distype> my;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So3_Poro_Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So3_Poro_Scatra(const So3_Poro_Scatra& old);

      /*!
      \brief Destructor

      */
      virtual ~So3_Poro_Scatra(){};

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
      //@}

     protected:
      //! don't want = operator
      So3_Poro_Scatra& operator=(const So3_Poro_Scatra& old);

    };  // class So3_Poro_Scatra


  }  // namespace ELEMENTS
}  // namespace DRT
#endif /* SO3_PORO_SCATRA_H_ */
