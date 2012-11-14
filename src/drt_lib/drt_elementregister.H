/*!----------------------------------------------------------------------
\file drt_elementregister.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef ELEMENTTYPE_H
#define ELEMENTTYPE_H


#include "Teuchos_RefCountPtr.hpp"

#include "drt_parobject.H"
#include "drt_element.H"
#include "drt_parobjectfactory.H"
#include "../linalg/linalg_sparsematrix.H"



namespace DRT
{

class ElementRegisterType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ElementRegisterType"; }

  static ElementRegisterType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ElementRegisterType instance_;
};

// forward declarations

/*!
\brief A virtual element register class

This class is used to initialize the finite elements of a specific type
in a discretization. Usage is as follows:<br>
Every concrete implementation should also have a concrete implementation
of the ElementRegister class, if it needs to have an initialization process
done before calculations on the elements can be performed. The method
DRT::ElementRegister::Initialize performs this initialization and is
called from DRT::Discretization::FillComplete ONCE for every type of
finite element included in a discretization.<br>
If a finite element does NOT need an initialization process, it may use
this base class as its ElementRegister as this base class implements a
dummy Initialize method.

\author gee (gee@lnm.mw.tum.de)
*/
class ElementRegister : public DRT::ParObject
{
public:




  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  ElementRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy

  */
  ElementRegister(const DRT::ElementRegister& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual DRT::ElementRegister* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~ElementRegister();

  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return ElementRegisterType::Instance().UniqueParObjectId(); }

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

  //@}

  //! @name Acess methods

  /*!
  \brief Get type of element
  */
  inline DRT::Element::ElementType Type() const { return etype_; }

  /*!
  \brief Print this element

  This class also implements the ostream << operator.
  */
  virtual void Print(ostream& os) const;

  //@}

  //! @name Construction

  /*!
  \brief Call Element init routines

  This method performs all steps needed to initialize the
  finite elements in the discretization dis when reimplemented
  by a derived element register class.<br>
  It is called once by DRT::Discretization::FillComplete for every
  type of finite element contained in dis.
  All necessary steps to init the finite elements of a certain type
   should be performed by the reimplementation of this method.

  \note The Discretization dis may contain a mixture of various types of
        elements, so care of this should be taken when reimplementing this
        method. The reimplementation should only
         initialize the ONE type of finite element it was written for,
         as other elements will have their own call to their distinct Initialize.

  \param dis (in/out): Discretization with finite elements that have to
                       be initialized.

  \note This class implements a dummy version of this method so finite
  elements that do not need an initialization can use this base class
  as their register class.
  */
  virtual int Initialize(DRT::Discretization& dis);


  /*!
  \brief Call pre-evaluation of elements

  \param dis (in/out): Discretization with finite elements that have to
                       be initialized.
  \param p   (in)    : parameter list of the evaluation call

  \note This class implements a dummy version of this method so finite
  elements that do not need pre-evaluation can use this base class dummy
  */
  virtual inline void PreEvaluate(DRT::Discretization& dis,
                                  Teuchos::ParameterList& p,
                                  RCP<LINALG::SparseOperator> systemmatrix1,
                                  RCP<LINALG::SparseOperator> systemmatrix2,
                                  RCP<Epetra_Vector>          systemvector1,
                                  RCP<Epetra_Vector>          systemvector2,
                                  RCP<Epetra_Vector>          systemvector3)
  { return; }
  //@}



private:

  //! type of element this register class refers to
  DRT::Element::ElementType           etype_;

}; // class ElementRegister
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::ElementRegister& ele);









#endif  // #ifndef ELEMENTTYPE_H
