/*!----------------------------------------------------------------------
\file shell8.H

\class CCADISCRETIZATION::Shell8

\brief A C++ wrapper for the shell8 element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef SHELL8_H
#define SHELL8_H


#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;

#include "element.H"


/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class Shell8

\brief A C++ wrapper for the shell8 element

\author gee (gee@lnm.mw.tum.de)
*/
class Shell8 : public CCADISCRETIZATION::Element
{
public:
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id : A unique global id
  */
  Shell8(int id);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  Shell8(const CCADISCRETIZATION::Shell8& old);

  /*!
  \brief Deep copy this instance of Shell8 and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  CCADISCRETIZATION::Element* Clone() const;

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data);

  /*!
  \brief Destructor

  */
  virtual ~Shell8();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Print this element
  */
  void Print() const;

  //@}

  // @{ \name Other

  /*!
  \brief Read input for this element
  */
  bool ReadElement();

  //@}


private:

  // don't want = operator
  Shell8 operator = (const Shell8& old);

  // element specific data
  
  enum ForceType
  {
    s8_none,
    s8_xyz,
    s8_rst,
    s8_rst_ortho
  };
  ForceType        forcetype_;  // type of force calculation
  double           thickness_;  // shell thickness
  int              ngp_[3];     // no. gaussian points in each spatial dimension
  int              ngptri_;     // no. gaussian points for triangle elements
  int              nhyb_;       // no. of internal dofs
  int              eas_[5];     // eas parameters
  int              ans_;        // type of ANS: 0=none 1=Q 2=T 3=QT  
  double           sdc_;        // scaled director conditioning
  int              material_;   // number of the material law

}; // class Shell8
} // namespace CCADISCRET











#endif  // #ifndef SHELL8_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
