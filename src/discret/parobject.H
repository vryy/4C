/*!----------------------------------------------------------------------
\file parobject.H

\class CCADISCRETIZATION::ParObject

\brief A virtual class with functionality to pack, unpack and communicate
       classes in parallel

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef PAROBJECT_H
#define PAROBJECT_H


using namespace std;


/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class ParObject

\brief A virtual class with functionality to pack, unpack and communicate
       classes in parallel

\author gee (gee@lnm.mw.tum.de)
*/
class ParObject
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  ParObject();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  ParObject(const CCADISCRETIZATION::ParObject& old);

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const = 0;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data) = 0;

  /*!
  \brief Destructor

  */
  virtual ~ParObject();

  //@}

  // @{ \name Acess methods

  //@}

protected:

}; // class ParObject
} // namespace CCADISCRETIZATION











#endif  // #ifndef PAROBJECT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
