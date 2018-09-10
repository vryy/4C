/*!----------------------------------------------------------------------
\file drt_element.H

\brief Base class of all elements used in BACI
\level 0
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
#ifndef ELEMENT_H
#define ELEMENT_H


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "drt_parobject.H"
#include "drt_condition.H"

// forward declarations
namespace MAT
{
  class Material;
}

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition;
  }

  namespace MESHFREE
  {
    template <typename>
    class MeshfreeBin;
  }

  namespace ELEMENTS
  {
    class ParamsInterface;
  }

  // forward declarations
  class Discretization;
  class DiscretizationHDG;
  class Node;
  class ElementType;

  class FaceElement;

  /*!
  \brief A virtual class all elements that are used in DRT have to implement

  This is the pure virtual base class for all finite elements to be used with
  the DRT::Discretization. Every element (and boundary condition) to be used
  with the discretization management module has to implement this class.
  It implements various basic element methods and
  stores basic information such as element to node connectivity.

  For elements that are attached to another 'parent' element, we use the derived
  class FaceElement. All boundary element types as well as interior face element
  classes are of type FaceElement. A FaceElement behaves as a usual Element most
  of the time, but it does have additional query methods for the parent element
  and other connectivity information.

  \author gee (gee@lnm.mw.tum.de)
  */
  class Element : public DRT::ParObject
  {
   public:
    //! @name Enums and Friends

    template <typename>
    friend class MESHFREE::MeshfreeBin;

    /*!
     \brief Element shape type
    */
    enum DiscretizationType
    {
      dis_none,    ///< unknown dis type
      quad4,       ///< 4 noded quadrilateral
      quad6,       ///< 6 noded quadrilateral (linear in one direction, quadratic in the other)
      quad8,       ///< 8 noded quadrilateral
      quad9,       ///< 9 noded quadrilateral
      tri3,        ///< 3 noded triangle
      tri6,        ///< 6 noded triangle
      hex8,        ///< 8 noded hexahedra
      hex16,       ///< 16 noded hexahedra
      hex18,       ///< 18 noded hexahedra
      hex20,       ///< 20 noded hexahedra
      hex27,       ///< 27 noded hexahedra
      tet4,        ///< 4 noded tetrahedra
      tet10,       ///< 10 noded tetrahedra
      wedge6,      ///< 6 noded wedge
      wedge15,     ///< 15 noded wedge
      pyramid5,    ///< 5 noded pyramid
      line2,       ///< 2 noded line
      line3,       ///< 3 noded line
      line4,       ///< 4 noded line
      line5,       ///< 5 noded line
      line6,       ///< 6 noded line
      point1,      ///< 1 noded point
      nurbs2,      ///< 2 control point first order nurbs line element
      nurbs3,      ///< 3 control point second order nurbs line element
      nurbs4,      ///< 4 control point first order nurbs surface element
      nurbs9,      ///< 9 control point second order nurbs surface element
      nurbs8,      ///< 8 control point first order nurbs volume element
      nurbs27,     ///< 27 control point second order nurbs volume element
      max_distype  ///< end marker. must be the last entry
    };

    //@}


    /// Location data for one dof set
    /*!
      A helper that manages location vectors. Required since there can be an
      arbitrary number of location vectors, matching the number of dofsets in
      the discretization.

      \author u.kue
      \date 12/09
     */
    class LocationData
    {
     public:
      /// clear all vectors
      void Clear()
      {
        lm_.clear();
        lmdirich_.clear();
        lmowner_.clear();
        stride_.clear();
      }

      /// return number of dofs collected
      int Size() const { return lm_.size(); }

      /// global dof numbers of elemental dofs
      std::vector<int> lm_;

      /// Dirichlet on/off flags
      std::vector<int> lmdirich_;

      /// Owner of dof (that is owner of node or element the dof belongs to)
      std::vector<int> lmowner_;

      /// Nodal stride (that is, how many dofs are guaranteed to be contiguous in the system matrix)
      std::vector<int> stride_;
    };


    /// Location data for all dof sets
    /*!
      A helper that manages location vectors. Required since there can be an
      arbitrary number of location vectors, matching the number of dofsets in
      the discretization.

      \author u.kue
      \date 12/09
     */
    class LocationArray
    {
     public:
      /// constructed with number of dofsets in discretization
      explicit LocationArray(int size) : data_(size) {}

      /// clear all location entries
      void Clear()
      {
        for (unsigned i = 0; i < data_.size(); ++i) data_[i].Clear();
      }

      /// access location entry
      LocationData& operator[](int i) { return data_[i]; }

      /// number of location entries, that is number of dofsets in discretization
      int Size() const { return data_.size(); }

     private:
      std::vector<LocationData> data_;
    };


    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id    (in): A globally unique element id
    \param owner (in): owner processor of the element
    */
    Element(int id, int owner);

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
    virtual int UniqueParObjectId() const = 0;

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this element

    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this element

    */
    virtual void Unpack(const std::vector<char>& data);

    /// return ElementType instance
    virtual DRT::ElementType& ElementType() const = 0;

    //@}

    //! @name Query methods

    /*!
    \brief Return global id of this element
    */
    int Id() const { return id_; }

    /*!
    \brief Return processor local col map id
    */
    virtual int LID() const { return lid_; }

    /*!
    \brief Return owner of this element
    */
    virtual int Owner() const { return owner_; }

    /*!
    \brief Get shape type of element
    */
    virtual DiscretizationType Shape() const = 0;

    /*!
    \brief Return number of nodes to this element
    */
    virtual int NumNode() const { return nodeid_.size(); }

    /*!
     * \brief Return number of points forming the geometry of this element
     *
     * For finite elements, this is equal to NumNodes(). Needed for meshfree
     * discretization where the integration cell (i.e. the "element") is not
     * formed by the nodes and many more nodes can have influence.
     */
    virtual int NumPoint() const { return nodeid_.size(); }

    /*!
    \brief Return number of lines to this element
    */
    virtual int NumLine() const { return 0; }

    /*!
    \brief Return number of surfaces to this element
    */
    virtual int NumSurface() const { return 0; }

    /*!
    \brief Return number of volumes to this element
    */
    virtual int NumVolume() const { return 0; }

    /*!
     * \brief Return the number of faces to the element.
     *
     * As opposed to NumLine(), NumSurface(), etc which have fixed dimension,
     * this is always the object one dimension less than the element dimension.
     */
    virtual int NumFace() const;

    /*!
    \brief Return id's of nodes adjacent to this element
    */
    virtual const int* NodeIds() const
    {
      if (nodeid_.size())
        return &nodeid_[0];
      else
        return NULL;
    }

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
    {
      if (node_.size())
        return &node_[0];
      else
        return NULL;
    }

    // mgee: change to this (as const and !const):
    // virtual std::vector<DRT::Node*>& Nodes()
    //{ return node_; }

    /*!
    \brief Get const vector of ptrs to nodes

    \warning The pointers to the nodes are build in
             DRT::Discretization::FillComplete. A standalone
             element that has not been added to a discretization
             (or the discretization has not been called FillComplete)
             does not have pointers to nodes. In this case, the method returns
             NULL.
    \return Ptr to pointers to nodes of the element in local nodal ordering.
            Returns NULL if pointers to not exist.
    */
    virtual const Node* const* Nodes() const
    {
      if (node_.size())
        return (const DRT::Node* const*)(&node_[0]);
      else
        return NULL;
    }

    /*!
     * \brief Return id's of points forming the geometry of this element
     *
     * For finite elements, this is equal to NodeIds(). Needed for meshfree
     * discretization where the integration cell (i.e. the "element") is not
     * formed by the nodes and many more nodes can have influence.
     */
    virtual const int* PointIds() const { return NodeIds(); }

    /*!
     * \brief Get vector of ptrs to points forming the geometry of this element
     *
     * For finite elements, this is equal to Nodes(). Needed for meshfree
     * discretization where the integration cell (i.e. the "element") is not
     * formed by the nodes and many more nodes can have influence.
     */
    virtual Node** Points() { return Nodes(); }

    /*!
     * \brief Get const vector of ptrs to points forming the geometry of this element
     *
     * For finite elements, this is equal to Nodes() const. Needed for meshfree
     * discretization where the integration cell (i.e. the "element") is not
     * formed by the nodes and many more nodes can have influence.
     */
    virtual const Node* const* Points() const { return Nodes(); }


    /*!
    \brief Get vector of Teuchos::RCPs to the lines of this element

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
    virtual std::vector<Teuchos::RCP<DRT::Element>> Lines()
    {
      return std::vector<Teuchos::RCP<DRT::Element>>(0);
    }

    // virtual const Element*const* Lines() const { dserror("unexpected base method called.");
    // return NULL; }

    /*!
    \brief Get vector of Teuchos::RCPs to the surfaces of this element

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
    virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces()
    {
      return std::vector<Teuchos::RCP<DRT::Element>>(0);
    }

    /*!
    \brief Get vector of Teuchos::RCPs to the volumes of this element

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
    virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes()
    {
      return std::vector<Teuchos::RCP<DRT::Element>>(0);
    }

    /*!
    \brief Set whether the given element is a nurbs element
    */
    virtual bool& SetNurbsElement() { return isNurbs_; }


    /*!
    \brief Get whether the given element is a nurbs element
    */
    virtual bool IsNurbsElement() const { return isNurbs_; }

    /*!
    \brief Returns whether the given element actually is a face element with degrees of freedom
    living there
    */
    virtual bool IsFaceElement() const { return false; }

    /*!
    \brief Get vector of Teuchos::RCPs to the faces of this element (as opposed to the Lines or
    Surfaces)
    */
    Teuchos::RCP<DRT::FaceElement>* Faces() { return face_.empty() ? NULL : &face_[0]; }

    /*!
    \brief Get vector of Teuchos::RCPs to the faces of this element (as opposed to the Lines or
    Surfaces)
    */
    Teuchos::RCP<DRT::FaceElement> const* Faces() const { return face_.empty() ? NULL : &face_[0]; }

    /*!
    \brief Get a pointer to the neighboring element behind the given face. Returns 0 if at boundary
    or faces are not created
    */
    DRT::Element* Neighbor(const int face) const;

    /*!
    \brief Construct a face element between this element and the given slave element

    This is a base class routine that constructs a plain element without association to physics,
    containing only the geometry and topology information. A surface/line element is created
    that represents the face between this element (master element) and parent_slave element (in case
    it exists). The derived element class is expected to provide physics information when allocate
    and store an surface/line element that represents the face between this element (master element)
    and parent_slave element.

    This element is e.g. used to create and evaluate edge stabilizations or for HDG discretizations

    \author schott 03/12
    */
    virtual Teuchos::RCP<DRT::Element> CreateFaceElement(
        DRT::Element* parent_slave,            //!< parent slave element
        int nnode,                             //!< number of nodes
        const int* nodeids,                    //!< node ids
        DRT::Node** nodes,                     //!< nodeids
        const int lsurface_master,             //!< local index of surface w.r.t master element
        const int lsurface_slave,              //!< local index of surface w.r.t slave element
        const std::vector<int>& localtrafomap  //!< local trafo map
    )
    {
      Teuchos::RCP<DRT::FaceElement> face;
      return face;
    }


    /*!
    \brief Get nodal connectivity and weights for nodes

    The method is used to build the connectivity between all the nodes adjacent
    to this element and how expensive its evaluation is

    \param edgeweights (out): A Epetra_SerialDenseMatrix containing weights of all connected nodes
    \param nodeweights (out): A Epetra_SerialDenseVector containing weights of all nodes
    */
    virtual void NodalConnectivity(
        Epetra_SerialDenseMatrix& edgeweights, Epetra_SerialDenseVector& nodeweights);

    /*!
    \brief Return value how expensive it is to evaluate this element

    \param double (out): cost to evaluate this element
    */
    virtual double EvaluationCost() { return 10.0; }

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
    virtual int NumDofPerNode(const DRT::Node& node) const
    {
      dserror("not implemented");
      return -1;
    }

    /*!
    \brief Get number of degrees of freedom per face

    The element decides how many element degrees of freedom it has.
    It can redecide along the way of a simulation.

    Usually, a standard finite element would not have face degrees of
    freedom and would therefore return zero here.

    \note Face degrees of freedom mentioned here are dofs that are supposed
          to be visible at the global system level. Purely internal
          element dofs that are condensed internally should NOT be considered.
          The DRT::Discretization will use this method to determine how many degrees
          of freedom it should include in the degree of freedom row and column map
          for the global system of equations.
    */
    virtual int NumDofPerFace(const unsigned /* face */) const { return 0; }

    /*!
    \brief Get number of degrees of freedom per component

    This is especially important for face elements in HDG discretizations: Since
    the trace field can be in terms of a scalar quantity, e.g. the pressure, or in terms
    of a vector quantity, e.g. the velocity, it is important to know, which case actually
    occurs. With this function, the number of scalar dofs may be asked. Standard elements
    return the number of nodes, where the number of dofs per node specifies the number of
    components. HDG face elements overwrite this function and return the number of scalar
    dofs.
    */
    virtual int NumDofPerComponent(const unsigned /* face */) const { return NumNode(); }

    /*!
    \brief Get the degree of an element

    The element can have variable degree, in case of HDG discretizations.
    HDG elements implement this functions, standard elements return the degree
    depending on their discretization type and assuming isoparametric concept.
    */
    virtual int Degree() const;


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
    virtual int NumDofPerElement() const
    {
      dserror("not implemented");
      return -1;
    }

    /*!
    \brief Print this element

    Prints basic information about this element to ostream. This method would
    usually be called by the print method of a derived class.
    */
    virtual void Print(std::ostream& os) const;

    /*!
    \brief Return the material of this element

    Note: The input parameter nummat is not the material number from input file
          as in SetMaterial(int matnum), but the number of the material within
          the vector of materials the element holds

    \param nummat (in): number of requested material
    */
    virtual Teuchos::RCP<MAT::Material> Material(int nummat = 0) const
    {
      dsassert(nummat < (int)mat_.size(), "invalid material number");
      return mat_[nummat];
    }

    /*!
    \brief Check whether the element has only ghost nodes

    Note: If the element has only ghost nodes, the element will not be allowed
          to assemble in any global vectors or matrixes, and, hence, can be
          skiped during Evaluate(). The only reason why it is ghosted on this
          proc is to provide access to its data. This might be necessary for volumetric
          coupling of non conforming meshes.

    \param mypid (in): id of calling processor
    */
    virtual bool HasOnlyGhostNodes(const int mypid) const;

    /*!
    \brief Query names of element data to be visualized using BINIO

    This method is to be overloaded by a derived class.
    The element is supposed to fill the provided map with key names of
    visualization data the element wants to visualize AT THE CENTER
    of the element geometry. The values is supposed to be dimension of the
    data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
    or 9 (nonsym. tensor)

    Example:
    \code
      // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
      names.insert(std::pair<string,int>("StressesXYZ",6));
    \endcode

    \param names (out): On return, the derived class has filled names with
                        key names of data it wants to visualize and with int dimensions
                        of that data.
    */
    virtual void VisNames(std::map<std::string, int>& names) { return; }

    /*!
    \brief Visuzalize the owner of the element using BINIO

    \param names (out): Owner is added to the key names
    */
    virtual void VisOwner(std::map<std::string, int>& names)
    {
      names.insert(std::pair<std::string, int>("Owner", 1));
      // names.insert(std::pair<string,int>("EleGId",1));
      return;
    }

    /*!
    \brief Query data to be visualized using BINIO of a given name

    This method is to be overloaded by a derived method.
    The derived method is supposed to call this base method to visualize the owner of
    the element.
    If the derived method recognizes a supported data name, it shall fill it
    with corresponding data.
    If it does NOT recognizes the name, it shall do nothing.

    \warning The method must not change size of variable data

    \param name (in):   Name of data that is currently processed for visualization
    \param data (out):  data to be filled by element if it recognizes the name
    */
    virtual bool VisData(const std::string& name, std::vector<double>& data)
    {
      if (name == "Owner")
      {
        if ((int)data.size() < 1) dserror("Size mismatch");
        data[0] = Owner();
        return true;
      }
      if (name == "EleGId")
      {
        if ((int)data.size() < 1) dserror("Size mismatch");
        data[0] = Id();
        return true;
      }
      return false;
    }

    //@}

    //! @name Construction

    /*!
    \brief Set global id of this element
    */
    void SetId(const int id) { id_ = id; }

    /*!
    \brief Read input for this element
    */
    virtual bool ReadElement(const std::string& eletype, const std::string& distype,
        DRT::INPUT::LineDefinition* linedef);

    /*!
      \brief Set processor local col id

      \param lid: processor local col id
     */
    void SetLID(int lid) { lid_ = lid; }

    /*!
    \brief Set ownership

    This method is used by the DRT::Discretiztion to change the ownership of
    an element that got communicated from one processor to another.

    \param owner: Proc owning this node

    \warning You should be very careful with changing the ownership
             of an element by hand as this might significantly confuse
             the DRT::Discretization the element is stored in.

    */
    void SetOwner(const int owner) { owner_ = owner; }

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

    /*!
    \brief Set a list of node ids this element is connected to

    Here the node ids are directly taken from an input line.
    */
    void SetNodeIds(const std::string& distype, DRT::INPUT::LineDefinition* linedef);

    /*!
    \brief Set a the face with index faceindex this element is connected to

    Sets the face pointer of the face adjacent to this element, using NumFace() as
    number of faces. This method is used in the construction phase of a discretization.
    It allows, that elements and faces are created separately and be combined later
    on the way.

    Faces, that are set via this method are created as non-owning (weak) RCPs in
    the std::vector face_.

    \param faceindex   : index of the given face
    \param faceelement : face object
    */
    void SetFace(const int faceindex, DRT::FaceElement* faceelement);

    /*!
    \brief Set a the face with index faceindex this element is connected to

    Sets the face pointer of the face adjacent to this element, using NumFace() as
    number of faces.

    Faces, that are set via this method copy the incoming RCP so strong and weak
    RCPs in the std::vector face_ can be created.

    \param faceindex   : index of the given face
    \param faceelement : face object
    */
    void SetFace(const int faceindex, Teuchos::RCP<DRT::FaceElement> faceelement);

    /// Set element material
    /*!
      Material numbers are read from the input file. The element stores
      a corresponding material object. These material objects can be
      anything from very simple (just a little calculation) to highly
      sophisticated with history data. The material is packed and
      unpacked along with its element.

      \param matnum : material number from input file
     */
    virtual void SetMaterial(int matnum);

    /// Add element material
    /*!
      In case of volume coupled problems the element needs information
      from the material of the other field. Therefore, it is possible to
      add pointers on other materials from other elements.

      \param mat: material to be added
      \param nummat (out):  number of materials the element holds
     */
    int AddMaterial(Teuchos::RCP<MAT::Material> mat);

    /// Number of materials of the element
    /*!
      In case of volume coupled problems the element needs information
      from the material of the other field. Therefore, it is possible to
      add pointers on other materials from other elements.

      \param nummat (out):  number of materials the element holds
     */
    int NumMaterial() const { return mat_.size(); };

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
    virtual void SetCondition(const std::string& name, Teuchos::RCP<Condition> cond)
    {
      condition_.insert(std::pair<std::string, Teuchos::RCP<Condition>>(name, cond));
      return;
    }

    /*!
    \brief Get all conditions with a certain name

    Get all conditions with a certain name. A vector of ptrs to all conditions
    with name name is returned in out. The number of conditions found with name
    name is out.size(). out.size() is 0 if no condition with that name is found.

    \param name (in): Name of condition
    \param out  (out): vector of pointers to all conditions with that name

    */
    virtual void GetCondition(const std::string& name, std::vector<DRT::Condition*>& out) const;

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
    virtual DRT::Condition* GetCondition(const std::string& name) const;

    /*!
    \brief Delete all conditions set to this element
    */
    virtual void ClearConditions()
    {
      condition_.clear();
      return;
    }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Return the location vector of this element

    Extended version that features a nodal dof set (nds) array. Each node might
    have multiple sets of (physical) dofs. That is the actual number of dofs at
    a node represent a (small) number of independent physical dofs. This is not
    to be confused with multiple DRT::DofSet objects within one
    DRT::Discretization.

    The case of multiple set of dofs occurs in xfem without enrichments. E.g. a
    node of a cutted fluid element might own two or more velocity-pressure
    pairs. In that case we need to choose one set of physical dofs for each
    node. The nds array gives the dof set number of each node.

    \note The degrees of freedom returned are not neccessarily only nodal dofs.
          Depending on the element implementation, output might also include
          element dofs.

    \param dis (in)      : the discretization this element belongs to
    \param nds           : dof set number of each node
    \param la (out)      : location data for all dofsets of the discretization
    \param doDirichlet (in): whether to get the Dirichlet flags

    */
    virtual void LocationVector(const DRT::Discretization& dis, const std::vector<int>& nds,
        DRT::Element::LocationArray& la, bool doDirichlet) const;

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
    \param la (out)      : location data for all dofsets of the discretization
    \param doDirichlet (in): whether to get the Dirichlet flags

    */
    virtual void LocationVector(
        const Discretization& dis, LocationArray& la, bool doDirichlet) const;


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
    This version is intended to fill the LocationArray with the dofs
    the element will assemble into. In the standard case these dofs are
    the dofs of the element itself. For some special conditions (e.g.
    the weak dirichlet boundary condtion) a surface element will assemble
    into the dofs of a volume element. These elements need to overwrite this
    method.<br>

    \note The degrees of freedom returned are not neccessarily only nodal dofs.
          Depending on the element implementation, output might also include
          element dofs.

    \param dis (in)      : the discretization this element belongs to
    \param la (out)      : location data for all dofsets of the discretization
    \param doDirichlet (in): whether to get the Dirichlet flags
    \param condstring (in): Name of condition to be evaluated
    \param condstring (in):  List of parameters for use at element level
    */
    virtual void LocationVector(const Discretization& dis, LocationArray& la, bool doDirichlet,
        const std::string& condstring, Teuchos::ParameterList& params) const;
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
    virtual void LocationVector(const Discretization& dis, std::vector<int>& lm,
        std::vector<int>& lmdirich, std::vector<int>& lmowner, std::vector<int>& lmstride) const;

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
    virtual void LocationVector(const Discretization& dis, std::vector<int>& lm,
        std::vector<int>& lmowner, std::vector<int>& lmstride) const;

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
    \param la (in)            : location data for all dofsets of the discretization
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
    virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        LocationArray& la, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3);

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
    virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
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
    virtual int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
        Epetra_SerialDenseMatrix* elemat1 = NULL) = 0;


    //@}


    //! @name Public methods to be used by DRT::Discretization only

    /*!
    \brief Build pointer vector from map of nodes

    \warning (public, but to be used by DRT::Discretization ONLY!)

    The method is used to build the variable node_ in this element. It is called from
    DRT::Discretization in DRT::Discretization::FillComplete() to create the
    pointers from elements to nodes (and nodes to elements)

    \param nodes (in): A map of all nodes of a discretization
    */
    virtual bool BuildNodalPointers(std::map<int, Teuchos::RCP<DRT::Node>>& nodes);

    /*!
    \brief Build pointer vector from vector of nodal pointers

    \warning (public, but to be used by DRT::Discretization ONLY!)

    The method is used to build the variable node_ in this element. It can be used
    to explicitly pass the nodal pointers to the element.

    \param nodes (in): Pointer to array of pointers to nodes. The array of pointers
                       is implicitly expected to be of length NumNode() and contain pointers
                       to nodes in the correct element local ordering scheme.
    */
    virtual bool BuildNodalPointers(DRT::Node** nodes);

    /*!
    \brief Build pointer vector from map of elements

    \warning (public, but to be used by DRT::Discretization ONLY!)

    The method is used to build the element connectivity in this element. For standard elements this
    procedure returns true and does nothing. For interface-elements a connection is made between the
    interface and it's left and right element. It is called from DRT::Discretization in
    DRT::Discretization::FillComplete() to create the pointers from elements to elements
    (next to the node-element and nodes-elements connectivity).

    \param elements (in): A map of all elements of a discretization
    */
    virtual bool BuildElementPointers(std::map<int, Teuchos::RCP<DRT::Element>>& elements)
    {
      return true;
    }

    //@}

    /// @name Parameter interface access and management functions
    /// @{
    /** \brief set the interface ptr for the elements
     *
     *  \param p (in): Parameter list coming from the time integrator.
     *
     *  \author hiermeier
     *  \date 04/16 */
    virtual void SetParamsInterfacePtr(const Teuchos::ParameterList& p)
  { /* This is a dummy function. Please implement the function in the derived classes, if necessary. */ };

    /** \brief returns true if the interface is defined and initialized, otherwise false
     *
     *  \author hiermeier
     *  \date 04/16 */
    virtual inline bool IsParamsInterface() const
    { /* dummy implementation */
      return false;
    }

    /** \brief get access to the interface pointer
     *
     *  \author hiermeier
     *  \date 04/16 */
    virtual Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr()
    {
      dserror(
          "This is a dummy function. Please implement the function in the derived classes, if "
          "necessary.");
      return Teuchos::null;
    }
    /// @}

   private:
    /*!
    \brief Default Constructor must not be called

    */
    Element()
    {
      dserror(
          "Default constructor of DRT::Element must not be called. Due to virtual"
          "inheritance from DRT::Element it can be necessary to call the non-default constructor "
          "explicitly.");
    }

    //! \brief A unique global element id
    int id_;

    //! local col map id
    int lid_;

    //! \brief owner processor of this element
    int owner_;

    //! \brief List of my nodal ids, length NumNode()
    std::vector<int> nodeid_;

    //! \brief Pointers to adjacent nodes in element local ordering
    std::vector<DRT::Node*> node_;

    //! \brief List of my faces, length NumFace(). Only filled if face elements are created, when
    //! using DiscretizationFaces
    std::vector<Teuchos::RCP<DRT::FaceElement>> face_;

    //! \brief Some conditions e.g. BCs
    std::multimap<std::string, Teuchos::RCP<Condition>> condition_;

    //! vector of material objects of element
    std::vector<Teuchos::RCP<MAT::Material>> mat_;

    //! is this a nurbs element?
    bool isNurbs_;

  };  // class Element


  /*!
  \brief A virtual class all face elements that are used in DRT have to implement

  This is the base class for all face elements. Every element that is attached to another
  (parent) element needs to implement this class, in particular boundary elements used for
  boundary conditions.

  A face element can be attached to one (for boundary faces) or two
  elements (for interior faces). In the latter case, we distinguish between the so-called
  parent master element and the parent slave element. The parent master element is
  the element that shares the node numbering with the face element. Boundary elements
  have only one parent and the methods ParentElement() and ParentMasterElement()
  return the same element. The parent slave element is the element on the other
  side that usually has a different orientation of nodes.

  \author kronbichler (kronbichler@lnm.mw.tum.de)
  */

  class FaceElement : public Element
  {
   public:
    /*!
    \brief Standard Constructor

    \param id    (in): A globally unique element id
    \param owner (in): owner processor of the element
    */
    FaceElement(int id, int owner);

    /*!
    \brief Copy Constructor
    */
    FaceElement(const FaceElement& old);

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this face element

    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this face element

    */
    virtual void Unpack(const std::vector<char>& data);

    /*!
    \brief Returns whether the given element actually is a face element with degrees of freedom
    living there
    */
    virtual bool IsFaceElement() const { return true; }

    /*!
    \brief Return the parent element id the face element is connected to (for interior faces, the
    parent master element).

    If no parent master has been assigned, -1 is returned (e.g. on non-face discretizations)

    This Id is also available, if the calling processor is not owner of the ParentElement()!
    */
    int ParentElementId() const { return parent_id_; }

    /*!
    \brief Return the parent element the face element is connected to (for interior faces, the
    parent master element).

    If no parent master has been assigned, NULL is returned (e.g. on non-face discretizations)
    */
    DRT::Element* ParentElement() const { return parent_master_; }

    /*!
    \brief Return the master element the face element is connected to

    If no parent master has been assigned, NULL is returned (e.g. on non-face discretizations)
    */
    DRT::Element* ParentMasterElement() const { return parent_master_; }

    /*!
    \brief Return the slave element the face element is connected to

    If no parent slave has been assigned, NULL is returned (non-face discretizations, boundary
    faces)
    */
    DRT::Element* ParentSlaveElement() const { return parent_slave_; }

    /*!
     \brief Get the index of a face element within the parent element
     */
    int FaceParentNumber() const
    {
      dsassert(lface_master_ != -1,
          "Face information has not been filled or this is not a face element");
      return lface_master_;
    }

    /*!
     \brief Get the index of a face element within the parent master element
     */
    int FaceMasterNumber() const
    {
      dsassert(lface_master_ != -1,
          "Face information has not been filled or this is not a face element");
      return lface_master_;
    }

    /*!
     \brief Get the index of a face element within the parent slave element
     */
    int FaceSlaveNumber() const
    {
      dsassert(lface_slave_ != -1,
          "Face information has not been filled or this is not a face element with slave parent");
      return lface_slave_;
    }

    /*!
    \brief return a tranformation for the face's nodes between the local coordinate systems of the
           face w.r.t the master parent element's face's coordinate system and the slave element's
           face's coordinate system

    Only filled for interior faces, otherwise zero length-vector (empty).
    */
    const std::vector<int>& GetLocalTrafoMap() const { return localtrafomap_; }

    /*!
    \brief Set the master element this face element is connected to

    Sets the pointers of the main parent element adjacent to the current face.
    This method is used in the construction phase of a discretization.
    It allows that elements and faces are created separately and be combined later
    on the way.

    \param master: parent element which shares the orientation of faces
    \param lface_master: local number of face within master element
    */
    void SetParentMasterElement(DRT::Element* master, const int lface_master)
    {
      parent_master_ = master;
      lface_master_ = lface_master;
      if (master != NULL) parent_id_ = master->Id();
    }

    /*!
    \brief Set the slave element this face element is connected to

    Sets the pointers of the second parent element adjacent to the current face.
    This method is used in the construction phase of a discretization.
    It allows that elements and faces are created separately and be combined later
    on the way.

    \param slave: second parent element when there already is a master element
    \param lface_slave: local number of face within slave element
    */
    void SetParentSlaveElement(DRT::Element* slave, const int lface_slave)
    {
      parent_slave_ = slave;
      lface_slave_ = lface_slave;
    }

   protected:
    /*!
    \brief Set the slave element this face element is connected to

    Sets the pointers of the second parent element adjacent to the current face.
    This method is used in the construction phase of a discretization.
    It allows that elements and faces are created separately and be combined later
    on the way.

    \param trace: transformation between
    */
    void SetLocalTrafoMap(const std::vector<int>& trafo);

   private:
    //! \brief The master parent element of face element. NULL for usual elements
    DRT::Element* parent_master_;

    //! \brief The master parent element of face element. NULL for usual elements
    DRT::Element* parent_slave_;

    //! The local surface number of this surface w.r.t to the parent_master_ element
    int lface_master_;

    //! The local surface number of this surface w.r.t to the parent_slave_ element
    int lface_slave_;

    /*!
     \brief map for the face's nodes between the local coordinate systems of the face w.r.t the
     master parent element's face's coordinate system and the slave element's face's coordinate
     system
     */
    std::vector<int> localtrafomap_;

    //! The parent element id
    int parent_id_;

    //! Friend declaration
    friend class DRT::DiscretizationHDG;
  };


  /*!
   * \brief translate to string for screen output
   */
  inline std::string DistypeToString(const DRT::Element::DiscretizationType distype)
  {
    std::string s = "";
    switch (distype)
    {
      case DRT::Element::dis_none:
        s = "DIS_NONE";
        break;
      case DRT::Element::quad4:
        s = "QUAD4";
        break;
      case DRT::Element::quad8:
        s = "QUAD8";
        break;
      case DRT::Element::quad9:
        s = "QUAD9";
        break;
      case DRT::Element::tri3:
        s = "TRI3";
        break;
      case DRT::Element::tri6:
        s = "TRI6";
        break;
      case DRT::Element::hex8:
        s = "HEX8";
        break;
      case DRT::Element::hex16:
        s = "HEX16";
        break;
      case DRT::Element::hex18:
        s = "HEX18";
        break;
      case DRT::Element::hex20:
        s = "HEX20";
        break;
      case DRT::Element::hex27:
        s = "HEX27";
        break;
      case DRT::Element::tet4:
        s = "TET4";
        break;
      case DRT::Element::tet10:
        s = "TET10";
        break;
      case DRT::Element::wedge6:
        s = "WEDGE6";
        break;
      case DRT::Element::wedge15:
        s = "WEDGE15";
        break;
      case DRT::Element::pyramid5:
        s = "PYRAMID5";
        break;
      case DRT::Element::line2:
        s = "LINE2";
        break;
      case DRT::Element::line3:
        s = "LINE3";
        break;
      case DRT::Element::point1:
        s = "POINT1";
        break;
      case DRT::Element::nurbs2:
        s = "NURBS2";
        break;
      case DRT::Element::nurbs3:
        s = "NURBS3";
        break;
      case DRT::Element::nurbs4:
        s = "NURBS4";
        break;
      case DRT::Element::nurbs9:
        s = "NURBS9";
        break;
      case DRT::Element::nurbs8:
        s = "NURBS8";
        break;
      case DRT::Element::nurbs27:
        s = "NURBS27";
        break;
      case DRT::Element::max_distype:
        s = "MAX_DISTYPE";
        break;
      default:
        std::cout << distype << std::endl;
        dserror("no string for this distype defined!");
        break;
    };
    return s;
  }

  /// translate shards::key to enum
  DRT::Element::DiscretizationType ShardsKeyToDisType(const unsigned& key);

  /// translate distype string to enum
  DRT::Element::DiscretizationType StringToDistype(std::string name);


}  // namespace DRT


// << operator
std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);



#endif  // #ifndef ELEMENT_H
