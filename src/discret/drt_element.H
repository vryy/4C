/*!----------------------------------------------------------------------
\file drt_element.H

\class DRT::Element

\brief A pure virtual element class

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef ELEMENT_H
#define ELEMENT_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
using namespace std;
using namespace Teuchos;
#include "drt_parobject.H"
#include "drt_condition.H"
#include "drt_dofset.H"


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations
class Discretization;
class Node;
class ElementRegister;
class ElementSurface;
class ElementLine;
class ElementVolume;

/*!
\class Element

\brief A virtual class all elements that are used in DRT have to implement

\author gee (gee@lnm.mw.tum.de)
*/
class Element : public DRT::ParObject
{
public:

  
  /*!
  \brief The Discretization is a friend of Element
  */
  friend class Discretization;
  
  /*!
  \brief Type of element
         
  */
  enum ElementType 
  {
    element_none,
    element_designline,
    element_designsurface,
    element_designvolume,
    element_shell8,
    element_elementvolume,
    element_elementsurface,
    element_elementline
  };


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id   : A globally unique element id
  \param etype: Type of element
  */
  Element(int id, ElementType etype, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  Element(const DRT::Element& old);

  /*!
  \brief Deep copy the derived class and return pointer to it
  
  */
  virtual DRT::Element* Clone() const = 0;

  /*!
  \brief Destructor

  */
  virtual ~Element();

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

  //@}

  // @{ \name Acess methods

  /*!
  \brief Return global id of this element
  */
  inline int Id() const { return id_; }

  /*!
  \brief Return owner of this element
  */
  virtual inline int Owner() const { return owner_; }

  /*!
  \brief Get type of element
  */
  inline ElementType Type() const { return etype_; }

  /*!
  \brief Allocate and return an ElementRegister class
  
  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const = 0;

  /*!
  \brief Return number of nodes of this element
  */
  inline int NumNode() const { return nodeid_.size(); }

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 0;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 0;}

  /*!
  \brief Return number of volumes of this element
  */
  virtual int NumVolume() const {return 0;}

  /*!
  \brief Return id's of nodes adjacent to this element
  */
  inline const int* NodeIds() const { if (nodeid_.size()) return &nodeid_[0]; 
                                      else                return NULL;     }

  /*!
  \brief Get vector of ptrs to nodes
  
  */
  virtual Node** Nodes()
  { if (node_.size()) return &node_[0]; else return NULL; }

  /*!
  \brief Get vector of ptrs to the lines of this element
  
  */
  virtual ElementLine** Lines() { return NULL; }

  /*!
  \brief Get vector of ptrs to the surfaces of this element
  
  */
  virtual ElementSurface** Surfaces() { return NULL; }

  /*!
  \brief Get vector of ptrs to the volumes of this element
  
  */
  virtual ElementVolume** Volumes() { return NULL; }

  /*!
  \brief Get number of degrees of freedom of a certain node
  
  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const = 0;

  /*!
  \brief Get number of degrees of freedom per element
  
  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.
  
  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const = 0;

  /*!
  \brief Return reference to the set of degrees of freedom
  */
  virtual inline DRT::DofSet& Dof() { return dofset_; }

  /*!
  \brief Return reference to the set of degrees of freedom
  */
  virtual inline const DRT::DofSet& Dof() const { return dofset_; }

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const;

  //@}

  // @{ \name Construction

  /*!
  \brief Read input for this element
  
  This class implements a dummy of this method that prints a warning and
  returns false.
  */
  virtual bool ReadElement();

  /*!
  \brief Set ownership
  
  \param owner: Proc owning this node
  */
  virtual inline void SetOwner(const int owner) { owner_ = owner; return; }

  /*!
  \brief Set a list of node ids this element is connected to
  
  Sets the nodal ids of the nodes adjacent to this element and the number
  of nodes
  
  \param nnode : number of nodes
  \param nodes : list of unique global nodal ids
  
  */
  void SetNodeIds(const int nnode, const int* nodes);

  /*!
  \brief Set a condition with a certain name
  
  Store a condition with a certain name in the element. The name need not
  be unique, meaning multiple conditions with the same name can be stored.
  Conditions can then be accessed with the GetCondition methods.
  
  \param name : Name of condition
  \param cond : The Condition class
  
  \note this takes ownership of cond
  
  \warning If a condition with the exact same name already exists, it will
           NOT be overwritten but stored twice in the element
  
  */
  virtual void SetCondition(const string& name,RefCountPtr<Condition> cond)
  { condition_.insert(pair<string,RefCountPtr<Condition> >(name,cond)); return; }

  /*!
  \brief Get all conditions with a certain name
  
  Get all conditions with a certain name. A vector of ptrs to all conditions
  with name name is returned in out. The number of conditions found with name
  name is out.size(). out.size() is 0 if no condition with that name is found.
  
  \param name (in): Name of condition
  \param out  (out): vector of pointers to all conditions with that name
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual void GetCondition(const string& name, vector<DRT::Condition*>& out);

  /*!
  \brief Get a condition with a certain name
  
  Returns the first condition with name name found in the multimap. 
  If multiple conditions with the same name exist, the first condition is
  returned and behaviour is therefore non-deterministic. this method should
  therefore only be used in cases where the user is sure that name is unique.
  
  \param name (in): Name of condition
  
  \return Returns NULL if condition with that name does not exist
  */
  virtual DRT::Condition* GetCondition(const string& name);

  /*!
  \brief Delete all conditions set to this element
  */
  virtual void ClearConditions() { condition_.clear(); return; }

  //@}

  // @{ \name Evaluation
  
  /*!
  \brief Return the location vector of this element
  
  The method computes degrees of freedom this element is using.
  Numbering is as follows:<br>
  First all degrees of freedom of nodes are numbered in nodal order,
  then the element internal degrees of freedom are given if present.<br>
  If your implementation of an Element has to use a different numbering scheme, 
  you are welcome to overload this method.<br>
  Length of the output vector must match number of degrees of freedom
  exactly.<br>
  
  \note The degrees of freedom returned are not neccessarily only nodal dofs.
        Depending on the element implementation, output also includes
        element dofs.
        
  \param lm (out): vector of degrees of freedom adressed by this element
  \param lmdirich (out): vector of zeros and ones indicating which
                         dofs have dirichlet boundary conditions. Ordering
                         matches dofs in lm.
  \param lmowner (out) : vector of proc numbers indicating which dofs are owned
                         by which procs in a dof row map. Ordering
                         matches dofs in lm.
  
  */
  virtual void LocationVector(vector<int>& lm, vector<int>& lmdirich, vector<int>& lmowner);
  
  /*!
  \brief Evaluate an element
  
  This class implements a dummy of this method that prints a warning and
  returns false.
  
  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input, 
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input, 
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input, 
                          the controlling method does not epxect the element
                          to fill this vector 
  \param elevec2 (out)  : vector to be filled by element. If NULL on input, 
                          the controlling method does not epxect the element
                          to fill this vector 
  \param elevec3 (out)  : vector to be filled by element. If NULL on input, 
                          the controlling method does not epxect the element
                          to fill this vector 
  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(ParameterList& params, 
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);
  


  //@}

protected:

  /*!
  \brief Build pointer vector from map of nodes
  
  */
  virtual bool BuildNodalPointers(map<int,RefCountPtr<DRT::Node> >& nodes);

  /*!
  \brief Build pointer vector from vector of nodal pointers
  
  */
  virtual bool BuildNodalPointers(DRT::Node** nodes); 

private:

  int                                      id_;        // a unique global element id
  int                                      owner_;     // owner of this element
  ElementType                              etype_;     // type of element
  vector<int>                              nodeid_;    // list of my nodal ids
  vector<DRT::Node*>                       node_;      // pointers to adjacent nodes
  DRT::DofSet                              dofset_;    // a set of degrees of freedom

  multimap<string,RefCountPtr<Condition> > condition_; // some conditions e.g. BCs

}; // class Element
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);









#endif  // #ifndef ELEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
