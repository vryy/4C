/*!----------------------------------------------------------------------
\file designdiscretization.H

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
#ifndef DESIGNDISCRETIZATION_H
#define DESIGNDISCRETIZATION_H

#include "discret.H"
#include "designelement.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class DesignDiscretization

\brief A class to manage a CAD design description

design inherits from Discretization

\author gee (gee@lnm.mw.tum.de)
*/
class DesignDiscretization : public CCADISCRETIZATION::Discretization
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  DesignDiscretization(RefCountPtr<Epetra_Comm> comm);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  DesignDiscretization(const CCADISCRETIZATION::DesignDiscretization& old);

  /*!
  \brief Destructor

  */
  virtual ~DesignDiscretization();

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
  int FillComplete(CCADISCRETIZATION::DesignDiscretization* higher, 
                   CCADISCRETIZATION::DesignDiscretization* lower);

  //@}

protected:




}; // class DesignDiscretization
} // namespace CCADISCRETIZATION











#endif  // #ifndef DESIGNDISCRETIZATION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
