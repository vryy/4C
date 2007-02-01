/*!----------------------------------------------------------------------
\file drt_elementsurface.H

\class DRT::ElementSurface

\brief An element representing a surface of a finite element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef ELEMENTSURFACE_H
#define ELEMENTSURFACE_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "drt_element.H"
#include "drt_elementregister.H"


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;
class Shell8Register;

/*!
\class ElementSurface

\brief An element representing a surface of a finite element

\author gee (gee@lnm.mw.tum.de)
*/
class ElementSurface : public DRT::Element
{
public:

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id : A unique global id
  \param owner: Processor owning this line
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  */
  ElementSurface(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  ElementSurface(const ElementSurface& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id
  
  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_ElementSurface;}

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
  virtual ~ElementSurface();

  //@}

  // @{ \name Acess methods

  
  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)
  
  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const {return 0;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)
  
  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.
  
  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Allocate and return a ElementRegister class
  
  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const;

  /*!
  \brief Integrate a Surface Neumann boundary condition
  
  */
  void SurfaceLoad(ParameterList& params,
                   Epetra_SerialDenseVector& elevec1,
                   const vector<double>& disp,
                   const vector<int>& lm);
  //@}


private:

  // don't want = operator
  ElementSurface operator = (const ElementSurface& old);
  

}; // class ElementSurface
} // namespace DRT











#endif  // #ifndef ELEMENTSURFACE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
