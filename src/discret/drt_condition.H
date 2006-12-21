/*!----------------------------------------------------------------------
\file drt_condition.H

\class DRT::Condition

\brief A condition of any kind

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef CONDITION_H
#define CONDITION_H


#include "Teuchos_RefCountPtr.hpp"

using namespace std;
using namespace Teuchos;

#include "drt_container.H"
#include "drt_dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations

/*!
\class Condition

\brief A condition of any kind

\author gee (gee@lnm.mw.tum.de)
*/
class Condition : public DRT::Container
{
public:


  /*!
  \brief Type of condition
         
  */
  enum ConditionType 
  {
    none,
    Dirichlet,
    PointNeumann,
    LineNeumann,
    SurfaceNeumann,
    VolumeNeumann
  };


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  Condition(const int id, const ConditionType type);

  /*!
  \brief Constructor with type condition_none
  
  */
  Condition();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Node
  
  */
  Condition(const DRT::Condition& old);

  /*!
  \brief Destructor

  */
  virtual ~Condition();

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this class
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual bool Unpack(const char* data);

  //@}

  // @{ \name Acess methods

  /*!
  \brief Print this element
  */
  inline virtual int Id() const { return id_; }

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Return type of condition
  */
  inline virtual ConditionType Type() const { return type_; }

  //@}

  // @{ \name Construction


  //@}

protected:

  // don't want = operator
  Condition operator = (const Condition& old);
  

protected:

  int           id_;     // unique id of this condition, no second condition of the same type with same id may exist
  ConditionType type_;   // type of this condition

}; // class Condition
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Condition& node);


#endif  // #ifndef CONDITION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
