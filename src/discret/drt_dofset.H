/*!----------------------------------------------------------------------
\file drt_dofset.H

\class DRT::DofSet

\brief A set of degrees of freedom

A set of degrees of freedom to be used in a node or an element or so.
This is a VERY simple class...

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

  Makes a deep copy
  
  */
  DofSet(const DRT::DofSet& old);

  /*!
  \brief Destructor

  */
  virtual ~DofSet();

  /*!
  \brief Return unique ParObject id
  
  see DRT::ParObject
  */
  virtual int UniqueParObjectId() const { return ParObject_DofSet;}

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this class
  
  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Unpack(const vector<char>& data);

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
  \brief Return one local degree of freedom
  
  \param index (in): local index of the degree of freedom to be accessed
  */
  virtual inline int operator [] (const int index) const 
  { 
#ifdef DEBUG
    if (NumDof()>index) return dofs_[index]; else return -1;
#else
    return dofs_[index];
#endif
  }

  //@}

  // @{ \name Construction

  /*!
  \brief change no. of degrees of freedom
  
  \param size (in) : new number of degrees of freedom in this set
  */
  virtual void SetNumDof(const int size); 

  /*!
  \brief Set degrees of freedom
  
  \param dofs (in): degrees of freedom to be stored in this class
  \param size (in): number of degrees of freedom in dofs
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
