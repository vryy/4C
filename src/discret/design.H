/*!----------------------------------------------------------------------
\file design.H

\class CCADISCRETIZATION::Design

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
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
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
  Design(const CCADISCRETIZATION::Design& old);

  /*!
  \brief Destructor

  */
  virtual ~Design();

  //@}

  // @{ \name Acess methods

  RefCountPtr<DesignDiscretization> operator [] (const int index) const
  { return entity_[index]; }

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
} // namespace CCADISCRETIZATION











#endif  // #ifndef DESIGN_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
