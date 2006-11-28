/*!----------------------------------------------------------------------
\file design.H

\class DRT::Design

\brief A class to manage a CAD design

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DESIGN_H
#define DESIGN_H

#include "designdiscretization.H"
#include "designelement.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

/*!
\class DesignDiscretization

\brief A class to manage a CAD design

\author gee (gee@lnm.mw.tum.de)
*/
class Design
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  Design(RefCountPtr<Epetra_Comm> comm);

  /*!
  \brief Copy Constructor

  */
  Design(const DRT::Design& old);

  /*!
  \brief Destructor

  */
  virtual ~Design();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Get communicator associated with this class
  */
  virtual inline const Epetra_Comm& Comm() const {return *comm_;}

  /*!
  \brief Test whether all discretizations herein have been called FillComplete()
  */
  virtual bool Filled() const 
  { 
    for (int i=0; i<3; ++i)
      if (!(entity_[i]->Filled())) 
        return false;
    return true;
  }

  /*!
  \brief Return one hierarchy level of entities
  
  Returns a handle to a design discretization entity level. 
  That is, given index = 0, returns the discretization containing lines and nodes, 
  index = 1 returns handle to surfaces, index = 2 returns handle to volumes.

  */
  RefCountPtr<DesignDiscretization> operator [] (const int index) const;

  //@}

protected:

  // don't want = operator
  Design operator = (const Design& old);

  RefCountPtr<Epetra_Comm>          comm_;       // epetra communicator
  /*! 
  entity_[0]: lines and nodes
  entity_[1]: surfaces
  entity_[2]: volumes
  */
  RefCountPtr<DesignDiscretization> entity_[3];  

}; // class Design
} // namespace DRT











#endif  // #ifndef DESIGN_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
