/*----------------------------------------------------------------------*/
/*!
 \file so_base.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SO_BASE_H_
#define SO_BASE_H_

#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_element.H"

namespace MAT
{
  class So3Material;
}


namespace DRT
{

namespace ELEMENTS
{

//! A wrapper for structural elements
class So_base : public DRT::Element
{

public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param owner (in): owner processor of the element
  */
  So_base(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_base(const So_base& old);

  //! Destructor
  virtual ~So_base(){};

  //@}

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

  // get the kinematic type from the element
  INPAR::STR::KinemType KinematicType() const
  {
    return kintype_;
  }

  // get the kinematic type from the element
  void SetKinematicType(INPAR::STR::KinemType kintype)
  {
    kintype_=kintype;
  }

  /*!
  \brief Return the material of this element

  Note: The input parameter nummat is not the material number from input file
        as in SetMaterial(int matnum), but the number of the material within
        the vector of materials the element holds

  \param nummat (in): number of requested material
  */
  virtual Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat=0) const;

protected:
  //! kinematic type
  INPAR::STR::KinemType kintype_;

private:
  /*!
  \brief Default Constructor must not be called

  */
  So_base();

};  // class So_base

}
}

#endif /* SO_BASE_H_ */
