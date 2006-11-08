/*!----------------------------------------------------------------------
\file design.H

\class CCADISCRETIZATION::Design

\brief A class to manage a CAD design description

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

#include "discret.H"
#include "designelement.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class Design

\brief A class to manage a CAD design description

design inherits from Discretization

\author gee (gee@lnm.mw.tum.de)
*/
class Design : public CCADISCRETIZATION::Discretization
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  Design(RefCountPtr<Epetra_Comm> comm);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  Design(const CCADISCRETIZATION::Design& old);

  /*!
  \brief Destructor

  */
  virtual ~Design();

  //@}

  // @{ \name Construction methods

  /*!
  \brief Complete construction of this design class
  
  This method calls the base class FillComplete() that creates maps and
  pointers between the elements and the nodes
  */
  int FillComplete() { return Discretization::FillComplete();}

  /*!
  \brief Complete construction of this design class
  
  This method creates pointers between the design elements stored in 
  the 'this' class and the elements stored in the lower class.
  It is used to create pointers between
  design-volumes <-> design-surfaces
  design-surfaces <-> design-lines
  It also calls the base class FillComplete method to create the element maps
  and set the flag Filled() to indicate that this class is complete.

  \param higher: Discretiztion holding the higher entities
  \param lower:  Discretiztion holding the lower  entities
  
  \note NULL is a valid argument for upper and/or lower

  */
  int FillComplete(CCADISCRETIZATION::Discretization* higher, 
                   CCADISCRETIZATION::Discretization* lower);

  //@}

protected:




}; // class Design
} // namespace CCADISCRETIZATION











#endif  // #ifndef DESIGN_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
