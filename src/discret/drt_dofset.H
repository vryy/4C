/*!----------------------------------------------------------------------
\file drt_dofset.H

\class DRT::DofSet

\brief A set of degrees of freedom

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DOFSET_H
#define DOFSET_H


using namespace std;
#include "drt_parobject.H"
#include "drt_dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations

/*!
\class DofSet

\brief A set of degrees of freedom

\author gee (gee@lnm.mw.tum.de)
*/
class DofSet : public DRT::ParObject
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  DofSet();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Node
  
  */
  DofSet(const DRT::DofSet& old);

  /*!
  \brief Destructor

  */
  virtual ~DofSet();

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
  \brief Print this class
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Return no. of degrees of freedom
  */
  virtual inline int NumDof() const { return (int)dofs_.size(); }

  /*!
  \brief Return degrees of freedom
  */
  virtual inline const int* Dofs() const 
  { if (NumDof()) return &dofs_[0]; else return NULL;}

  /*!
  \brief Return a degree of freedom
  */
  virtual inline int operator [] (const int index) const 
  { if (NumDof()>index) return dofs_[index]; else return -1;}

  //@}

  // @{ \name Construction

  /*!
  \brief change no. of degrees of freedom
  */
  virtual void SetNumDof(const int size); 

  /*!
  \brief Set degrees of freedom
  */
  virtual void SetDof(const int* dofs, const int size); 

  //@}

private:

  // don't want = operator
  DofSet operator = (const DofSet& old);
  

protected:

  vector<int>   dofs_;  // vector of degrees of freedom

}; // class DofSet
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::DofSet& dofset);


#endif  // #ifndef DOFSET_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
