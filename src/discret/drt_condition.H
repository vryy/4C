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


#include "Epetra_Comm.h"
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
class Element;

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

  /*!
  \brief Return communicator
  */
  inline const RefCountPtr<Epetra_Comm> Comm() const { return comm_; }

  //@}

  // @{ \name Construction

  /*!
  \brief Set a communicator
  */
  inline void SetComm(RefCountPtr<Epetra_Comm> comm) 
  { comm_ = comm; return; }

  /*!
  \brief Add a geometry description to the condition
  
  A geometry condition can be added to the condition.
  In case the condition refers to lines, surfaces or volumes, a
  geometry description is needed to properly evaluate this condition.
  Such a geometry description is build in \ref DRT::Discretization::BoundaryConditionsGeometry
  and then added to this Condition.
  
  \param geom (in): Map of elements describing the geometry.
                    A deep copy of the map is made and stored.
                    Normally though, these elements are a line, surface or
                    volume description produced and shared with the discretization.
                    Do not mess with the RefCountPtr!
  
  */
  virtual void AddGeometry(map<int,RefCountPtr<DRT::Element> >& geom) 
  { geometry_ = geom; return; }

  /*!
  \brief Delete a geometry description of the condition
  
  */
  virtual void ClearGeometry() { geometry_.clear(); return; }

  /*!
  \brief Get a reference to the geometry description of the condition
  
  */
  virtual map<int,RefCountPtr<DRT::Element> >& Geometry() { return geometry_; }

  //@}

protected:

  // don't want = operator
  Condition operator = (const Condition& old);
  

protected:

  int                                 id_;       // unique id of this condition, no second condition of the same type with same id may exist
  ConditionType                       type_;     // type of this condition
  map<int,RefCountPtr<DRT::Element> > geometry_; // geometry description of this condition
  RefCountPtr<Epetra_Comm>            comm_;     // a communicator 
}; // class Condition
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Condition& node);


#endif  // #ifndef CONDITION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
