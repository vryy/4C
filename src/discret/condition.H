/*!----------------------------------------------------------------------
\file condition.H

\class CCADISCRETIZATION::Condition

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

#include "container.H"
#include "dserror.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

// forward declarations

/*!
\class Condition

\brief A condition of any kind

\author gee (gee@lnm.mw.tum.de)
*/
class Condition : public CCADISCRETIZATION::Container
{
public:


  /*!
  \brief Type of condition
         
  */
  enum ConditionType 
  {
    condition_none,
    condition_Dirichlet,
    condition_Neumann
  };


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  Condition(ConditionType type);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Node
  
  */
  Condition(const CCADISCRETIZATION::Condition& old);

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
  virtual void Print(ostream& os) const ;

  /*!
  \brief Print this element
  */
  inline virtual ConditionType Type() const { return type_; }

  //@}

  // @{ \name Construction


  //@}

protected:

  // don't want = operator
  Condition operator = (const Condition& old);
  

protected:

  ConditionType type_;   // type of this condition

}; // class Condition
} // namespace CCADISCRETIZATION


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Condition& node);


#endif  // #ifndef CONDITION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
