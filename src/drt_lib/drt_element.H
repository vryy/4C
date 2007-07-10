/*!----------------------------------------------------------------------
\file drt_element.H

\class DRT::Element

\brief A pure virtual element class

This is the pure virtual base class for all finite elements to be used with
the DRT::Discretization. Every element (and boundary condition) to be used
with the discretization management module has to implement this class.
It implements various basic element methods and
stores basic information such as element to node connectivity.

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

// forward declarations
namespace MAT
{
  class Material;
}

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations
class Discretization;
class Node;
class ElementRegister;

/*!
\class Element

\brief A virtual class all elements that are used in DRT have to implement

This is the pure virtual base class for all finite elements to be used with
the DRT::Discretization. Every element (and boundary condition) to be used
with the discretization management module has to implement this class.
It implements various basic element methods and
stores basic information such as element to node connectivity.

\author gee (gee@lnm.mw.tum.de)
*/
class Element : public DRT::ParObject
{
public:

  //! @name Enums and Friends

  /*!
  \brief Type of element

  */
  enum ElementType
  {
    element_none,
    element_shell8,
    element_shell8line,
    element_wall1,
    element_wall1line,
    element_fluid2,
    element_fluid2line,
    element_condif2,
    element_condif2line,
    element_fluid3,
    element_fluid3surface,
    element_fluid3line,
    element_ale3,
    element_ale3surface,
    element_so_hex8,
    element_soh8surface,
    element_soh8line,
    element_soh8volume,
    element_sosh8,
    element_ale2,
    element_ale2line,
    element_so_tet10,
    element_sotet10surface,
    element_sotet10line,
    element_sotet10volume,
  };

  /*!
   \brief Element shape type
  */
  enum DiscretizationType
  {
    dis_none,                   // unknown dis type
    quad4,                      // 4 noded quadrilateral
    quad8,                      // 8 noded quadrilateral
    quad9,                      // 9 noded quadrilateral
    tri3,                       // 3 noded triangle
    tri6,                       // 6 noded triangle
    hex8,                       // 8 noded hexahedra
    hex20,                      // 20 noded hexahedra
    hex27,                      // 27 noded hexahedra
    tet4,                       // 4 noded tetrahedra
    tet10,                      // 4 noded tetrahedra
    wedge6,                     // 6 noded wedge
    wedge15,                    // 15 noded wedge
    line2,                      // 2 noded line
    line3,                      // 3 noded line
    max_distype                 //  end marker. must be the last entry
  };

  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Element(int id, ElementType etype, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Element(const DRT::Element& old);

  /*!
  \brief Destructor

  */
  virtual ~Element();

  /*!
  \brief Deep copy the derived class and return pointer to it

  This method is sort of a copy constructor for a class derived from Element.
  It allows to copy construct the derived class without knowing what it
  actually is using the base class Element.

  */
  virtual DRT::Element* Clone() const = 0;


  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return ParObject_Element; }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Allocate and return an ElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  The ElementRegister class is then used to hold track of what type of
  elements are in one discretization.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const = 0;

  //@}

  //! @name Query methods

  /*!
  \brief Return global id of this element
  */
  inline int Id() const { return id_; }

  /*!
  \brief Return processor local col map id
  */
  virtual inline int LID() const { return lid_; }

  /*!
  \brief Return owner of this element
  */
  virtual inline int Owner() const { return owner_; }

  /*!
  \brief Get type of element
  */
  inline ElementType Type() const { return etype_; }

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const = 0;

  /*!
  \brief Return number of nodes to this element
  */
  inline int NumNode() const { return nodeid_.size(); }

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 0;}

  /*!
  \brief Return number of surfaces to this element
  */
  virtual int NumSurface() const {return 0;}

  /*!
  \brief Return number of volumes to this element
  */
  virtual int NumVolume() const {return 0;}

  /*!
  \brief Return id's of nodes adjacent to this element
  */
  inline const int* NodeIds() const { if (nodeid_.size()) return &nodeid_[0];
                                      else                return NULL;     }

  /*!
  \brief Get vector of ptrs to nodes

  \warning The pointers to the nodes are build in
           DRT::Discretization::FillComplete. A standalone
           element that has not been added to a discretization
           (or the discretization has not been called FillComplete)
           does not have pointers to nodes. In this case, the method returns
           NULL.
  \return Ptr to pointers to nodes of the element in local nodal ordering.
          Returns NULL if pointers to not exist.
  */
  virtual Node** Nodes()
  { if (node_.size()) return &node_[0]; else return NULL; }

  /*!
  \brief Get vector of ptrs to the lines of this element

  This is a base class dummy routine that always returns NULL.
  The derived element class is expected to allocate and store
  a vector of elements that represent the edges of this element.
  These edges are then used to create and evaluate boundary conditions.

  \note A 1D type of element (e.g. a beam) may return
       a vector of length one pointing to itself. It then does not need
       to be able to spin of an explicit separate line element.
       In this case, the element itself must be able to evaluate the
       Neumann boundary conditions on a line itself.

  */
  virtual Element** Lines() { return NULL; }

  /*!
  \brief Get vector of ptrs to the surfaces of this element

  This is a base class dummy routine that always returns NULL.
  The derived element class is expected to allocate and store
  a vector of elements that represent the surfaces of this element.
  These surfaces are then used to create and evaluate boundary conditions.

  \note A 2D type of element (e.g. a shell, wall, fluid2 etc) may return
       a vector of length one pointing to itself. It then does not need
       to be able to spin of an explicit separate surface element.
       In this case, the element itself must be able to evaluate the
       Neumann boundary conditions on a surface itself.

  */
  virtual Element** Surfaces() { return NULL; }

  /*!
  \brief Get vector of ptrs to the volumes of this element

  This is a base class dummy routine that always returns NULL.
  The derived element class is expected to allocate and store
  a vector of elements that represent the volume of this element.
  Length of this vector would usually be just one.
  This volume element is then used to create and evaluate boundary conditions.

  \note A 3D type of element (e.g. a solid, fluid3, etc) may return
       a vector of length one pointing to itself. It then does not need
       to be able to spin of an explicit separate volume element.
       In this case, the element itself must be able to evaluate the
       Neumann boundary conditions on a volume itself.
  */
  virtual Element** Volumes() { return NULL; }

  /*!
  \brief Get number of degrees of freedom of a certain node

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.<br>
  This method is used by the DRT::Discretization to determine how many
  degrees of freedom should be assigned to each node. The discretization will
  assign a node the max number of dofs the adjacent elements demand by this
  method.

  */
  virtual int NumDofPerNode(const DRT::Node& node) const = 0;

  /*!
  \brief Get number of degrees of freedom per element

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.<br>
  Usually, a standard finite element would not have element degrees of
  freedom and would therefore return zero here.

  \note Element degrees of freedom mentioned here are dofs that are supposed
        to be visible at the global system level. Purely internal
        element dofs that are condensed internally should NOT be considered.
        The DRT::Discretization will use this method to determine how many degrees
        of freedom it should include in the degree of freedom row and column map
        for the global system of equations.
  */
  virtual int NumDofPerElement() const = 0;

  /*!
  \brief Print this element

  Prints basic information about this element to ostream. This method would
  usually be called by the print method of a derived class.
  */
  virtual void Print(ostream& os) const;

  /// return the material of this element
  RefCountPtr<MAT::Material> Material() const { return mat_; }

  //@}

  //! @name Construction

  /*!
  \brief Set global id of this element
  */
  virtual inline void SetId(const int id) { id_ = id; return; }

  /*!
  \brief Read input for this element

  This class implements a dummy of this method that prints a warning and
  returns false. A derived class would read one line from the input file and
  store all necessary information.

  */
  virtual bool ReadElement();

  /*!
    \brief Set processor local col id
    \param lid: processor local col id
   */
  virtual inline void SetLID(int lid) { lid_ = lid; }

  /*!
  \brief Set ownership

  This method is used by the DRT::Discretiztion to change the ownership of
  an element that got communicated from one processor to another.

  \param owner: Proc owning this node

  \warning You should be very careful with changing the ownership
           of an element by hand as this might significantly confuse
           the DRT::Discretization the element is stored in.

  */
  virtual inline void SetOwner(const int owner) { owner_ = owner; return; }

  /*!
  \brief Set Type

  This method is used by the DRT::Discretiztion to change the type of
  an element.

  \param type: elementtype to be set

  */
  virtual inline void SetType(const ElementType type) { etype_ = type; return; }

  /*!
  \brief Set a list of node ids this element is connected to

  Sets the nodal ids of the nodes adjacent to this element and the number
  of nodes. This method is used in the construction phase of a discretization.
  It allows, that elements and nodes are created separately and be combined later
  on the way.

  \param nnode : number of nodes
  \param nodes : list of unique global nodal ids

  */
  void SetNodeIds(const int nnode, const int* nodes);

  /// Set element material
  /*!
    Material numbers are read from the input file. The element stores
    a corresponding material object. These material objects can be
    anything from very simple (just a little calculation) to highly
    sophisticated with history data. The material is packed and
    unpacked along with its element.

    \param matnum : material number from input file
   */
  void SetMaterial(int matnum);

  //@}

  //! @name Conditions

  /*!
  \brief Set a condition with a certain name

  Store a condition with a certain name in the element. The name need not
  be unique, meaning multiple conditions with the same name can be stored.
  Conditions can then be accessed with the GetCondition methods.

  \param name : Name of condition
  \param cond : The Condition class

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

  */
  virtual void GetCondition(const string& name, vector<DRT::Condition*>& out);

  /*!
  \brief Get a condition with a certain name

  Returns the first condition with name name found in the multimap.
  If multiple conditions with the same name exist, the first condition with
  that name is returned and behaviour is therefore non-deterministic.
  This method should therefore only be used in cases where the user is
  positive that name is unique.

  \param name (in): Name of condition

  \return Returns NULL if condition with that name does not exist
  */
  virtual DRT::Condition* GetCondition(const string& name);

  /*!
  \brief Delete all conditions set to this element
  */
  virtual void ClearConditions() { condition_.clear(); return; }

  //@}

  //! @name Evaluation methods

  /*!
  \brief Return the location vector of this element

  The method computes degrees of freedom this element adresses.
  Degree of freedom ordering is as follows:<br>
  First all degrees of freedom of adjacent nodes are numbered in
  local nodal order, then the element internal degrees of freedom are
  given if present.<br>
  If a derived element has to use a different ordering scheme,
  it is welcome to overload this method as the assembly routines actually
  don't care as long as matrices and vectors evaluated by the element
  match the ordering, which is implicitly assumed.<br>
  Length of the output vector matches number of degrees of freedom
  exactly.<br>

  \note The degrees of freedom returned are not neccessarily only nodal dofs.
        Depending on the element implementation, output might also include
        element dofs.

  \param dis (in)      : the discretization this element belongs to
  \param lm (out)      : vector of degrees of freedom adressed by this element
  \param lmdirich (out): vector of zeros and ones indicating which
                         dofs have dirichlet boundary conditions. Ordering
                         matches dofs in lm.
  \param lmowner (out) : vector of proc numbers indicating which dofs are owned
                         by which procs in a dof row map. Ordering
                         matches dofs in lm.

  */
  virtual void LocationVector(const Discretization& dis, vector<int>& lm, vector<int>& lmdirich, vector<int>& lmowner);

  /*!
  \brief Return the location vector of this element

  The method computes degrees of freedom this element adresses.
  Degree of freedom ordering is as follows:<br>
  First all degrees of freedom of adjacent nodes are numbered in
  local nodal order, then the element internal degrees of freedom are
  given if present.<br>
  If a derived element has to use a different ordering scheme,
  it is welcome to overload this method as the assembly routines actually
  don't care as long as matrices and vectors evaluated by the element
  match the ordering, which is implicitly assumed.<br>
  Length of the output vector matches number of degrees of freedom
  exactly.<br>

  \note The degrees of freedom returned are not neccessarily only nodal dofs.
        Depending on the element implementation, output might also include
        element dofs.

  \param dis (in)      : the discretization this element belongs to
  \param lm (out)      : vector of degrees of freedom adressed by this element
  \param lmowner (out) : vector of proc numbers indicating which dofs are owned
                         by which procs in a dof row map. Ordering
                         matches dofs in lm.

  */
  virtual void LocationVector(const Discretization& dis, vector<int>& lm, vector<int>& lmowner);

  /*!
  \brief Evaluate an element

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates element matrices and
  vectors accoring to the command in params.

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elemat2 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elevec1 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec2 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec3 (out)      : vector to be filled by element depending on commands
                              given in params
  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(ParameterList&            params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);

  /*!
  \brief Evaluate a Neumann boundary condition

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


  //@}


  /*!
  \brief Build pointer vector from map of nodes

  The method is used to build the variable node_ in this element. It is called from
  DRT::Discretization in DRT::Discretization::FillComplete() to create the
  pointers from elements to nodes (and nodes to elements)

  \param nodes (in): A map of all nodes of a discretization
  */
  virtual bool BuildNodalPointers(map<int,RefCountPtr<DRT::Node> >& nodes);

  /*!
  \brief Build pointer vector from vector of nodal pointers

  The method is used to build the variable node_ in this element. It can be used
  to explicitly pass the nodal pointers to the element.

  \param nodes (in): Pointer to array of pointers to nodes. The array of pointers
                     is implicitly expected to be of length NumNode() and contain pointers
                     to nodes in the correct element local ordering scheme.
  */
  virtual bool BuildNodalPointers(DRT::Node** nodes);


private:

  //! \brief A unique global element id
  int id_;

  //! local col map id
  int lid_;

  //! \brief owner processor of this element
  int owner_;

  //! \brief Type of element
  ElementType etype_;

  //! \brief List of my nodal ids, length NumNode()
  vector<int> nodeid_;

  //! \brief Pointers to adjacent nodes in element local ordering
  vector<DRT::Node*> node_;

  //! \brief Some conditions e.g. BCs
  multimap<string,RefCountPtr<Condition> > condition_;

  //! material object of element
  RefCountPtr<MAT::Material> mat_;

}; // class Element
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);









#endif  // #ifndef ELEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
