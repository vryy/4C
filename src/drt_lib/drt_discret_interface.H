/*---------------------------------------------------------------------*/
/*!

\brief Interface class for DRT::Discretization, the main discretization class

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_DISCRET_INTERFACE_H_
#define SRC_DRT_LIB_DRT_DISCRET_INTERFACE_H_

#include "drt_dserror.H"
#include <Teuchos_RCP.hpp>

// forward declarations ...
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_SerialDenseVector;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace LINALG
{
  class SparseOperator;
  class MapExtractor;
}  // namespace LINALG
namespace IO
{
  class DiscretizationWriter;
}  // namespace IO
namespace DRT
{
  class Node;
  class Element;
  class Condition;

  /** \class DiscretizationInterface
   *
   *  \brief Generic discretization interface class.
   *
   *  This class represents the generic base class of all discretization
   *  classes. It is supposed to be used as an interface to different
   *  discretization objects. One example are function calls where a concrete
   *  discretization or a discretization wrapper can be handed in. See for
   *  example the use in the structural time integration.
   *
   *  \note You can add additional pure virtual functions if you need them.
   *  But do not add any concrete function implementations and do not
   *  add a cpp-file! This class is constructed as a pure virtual interface and
   *  should stay a pure virtual interface! All concrete implementations can be
   *  found in the derived classes.
   *
   *  The only exception of the rule are wrapper functions.
   *
   *  \author hiermeier
   *  \date 10/16 */
  class DiscretizationInterface
  {
   public:
    /// constructor
    DiscretizationInterface(){};

    /// destructor
    virtual ~DiscretizationInterface(){};

    /// Access the communicator object
    virtual const Epetra_Comm& Comm() const = 0;

    /// Get communicator associated with this class
    virtual Teuchos::RCP<IO::DiscretizationWriter> Writer() const = 0;

    /// Get flag indicating whether FillComplete() has been called
    virtual bool Filled() const = 0;

    /// Get name of the discretization
    virtual const std::string& Name() const = 0;

    /** \brief Print this discretization to os (Filled()==true NOT prerequisite)
     *
     *  \note This is a collective call */
    virtual void Print(std::ostream& os) const = 0;

    /** \brief Complete construction of a discretization  (Filled()==true NOT prerequisite)
     *
     *  After adding or deleting nodes or elements or redistributing them in parallel,
     *  or adding/deleting boundary conditions, this method has to be called to (re)construct
     *  pointer topologies.<br>
     *  It builds in this order:<br>
     *  - row map of nodes
     *  - column map of nodes
     *  - row map of elements
     *  - column map of elements
     *  - pointers from elements to nodes
     *  - pointers from nodes to elements
     *  - assigns degrees of freedoms
     *  - map of element register classes
     *  - calls all element register initialize methods
     *  - build geometries of all Dirichlet and Neumann boundary conditions
     *
     *  \param assigndegreesoffreedom (in): if true, resets existing dofsets and performs
     *                                      assigning of degrees of freedoms to nodes and
     *                                      elements.
     *  \param initelements           (in): if true, build element register classes and call
     *                                      Initialize() on each type of finite element present
     *  \param doboundaryconditions   (in): if true, build geometry of boundary conditions
     *                                      present.
     *
     *  \note In order to receive a fully functional discretization, this method must be called
     *        with all parameters set to true (the default). The parameters though can be
     *        used to turn off specific tasks to allow for more flexibility in the
     *        construction of a discretization, where it is known that this method will
     *        be called more than once.
     *
     *  \note Sets Filled()=true */
    virtual int FillComplete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true) = 0;

    /** \brief Get map associated with the distribution of the ownership of nodes
     *         (Filled()==true prerequisite)
     *
     * This map includes all nodes stored on this proc and also owned by this proc.
     * This map is non-ambiguous, meaning that it is a non-overlapping map.
     *
     * \return NULL if Filled() is false. A call to FillComplete() is a prerequisite. */
    virtual const Epetra_Map* NodeRowMap() const = 0;

    /** \brief Get map associated with the distribution of the ownership of elements
     *         (Filled()==true prerequisite)
     *
     *  This map includes all elements stored on this proc and also owned by this proc.
     *  This map is non-ambiguous, meaning that it is a non-overlapping map.
     *
     *  \return NULL if Filled() is false. A call to FillComplete() is a prerequisite. */
    virtual const Epetra_Map* ElementRowMap() const = 0;

    /** \brief Get processor local number of nodes owned by this processor
     *         (Filled()==true prerequisite) */
    virtual int NumMyRowNodes() const = 0;

    /** \brief Query whether a Node with global id gid is stored on this proc
     *         (Filled()==true NOT prerequisite) */
    virtual bool HaveGlobalNode(int gid) const = 0;

    /** \brief Get the node with global row id gid (Filled()==true NOT prerequisite)
     *
     *  Returns the node with global row id gid if node is on this proc.
     *  Will return row or column node, ghosted or  not.
     *  This is an individual call
     *
     *  \return Address of node if node is stored on calling proc */
    virtual DRT::Node* gNode(int gid) const = 0;

    /** \brief Get the element with global id gid (Filled()==true NOT prerequisite)
     *
     *  Returns the element with global row id gid if element is on this proc.
     *  Will return row or column element, ghosted or not.
     *  This is an individual call. Will test for existence of element in
     *  DEBUG version and throw an error message if not. Will crash in non-DEBUG
     *  version if element does not exist on calling processor.
     *
     *  \return Address of element if element is owned by calling proc, returns NULL
     *          otherwise */
    virtual DRT::Element* gElement(int gid) const = 0;

    /** \brief Get the node with local row id lid (Filled()==true prerequisite)
     *  Returns the node with local row index lid.
     *  Will not return any ghosted node.
     *  This is an individual call and Filled()=true is a prerequisite
     *
     *  \return Address of node if node is owned and stored by calling proc */
    virtual DRT::Node* lRowNode(int lid) const = 0;

    /** \brief Get processor local number of elements owned by this processor
     *         (Filled()==true prerequisite) */
    virtual int NumMyRowElements() const = 0;

    /** \brief Get the element with local row id lid (Filled()==true prerequisite)
     *
     *  Returns the element with local row index lid.
     *  Will not return any ghosted element.
     *  This is an individual call and Filled()=true is a prerequisite
     *
     *  \return Address of element if element is owned by calling proc */
    virtual DRT::Element* lRowElement(int lid) const = 0;

    /** \brief Get degree of freedom row map (Filled()==true prerequisite)
     *
     *  Return ptr to degree of freedom row distribution map of this discretization.
     *  If it does not exist yet, build it.
     *
     *  - Filled()==true prerequisite
     *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom())) */
    virtual const Epetra_Map* DofRowMap(unsigned nds = 0) const = 0;

    /** \brief True if degrees of freedom where assigned
     *
     *  Degrees of freedom need to be assigned using AssignDegreesOfFreedom()
     *  before any calculations using the discretization can be made */
    virtual bool HaveDofs() const = 0;

    /// Get the number of DoF sets
    virtual int NumDofSets() const = 0;

    /** \brief Get number of dofs for given node.
     *
     *  Ask the current DofSet for the number of dofs of this node.
     *  There is a variable number of DofSets and one is currently selected.
     *
     *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
     *
     *  \param dofset_id  (in): id of the desired dofset
     *  \param node       (in): the node those number of dofs are requested */
    int NumDof(const Node* node) const
    {
      dsassert(NumDofSets() == 1, "expect just one dof set");
      return NumDof(0, node);
    };
    virtual int NumDof(unsigned dofset_id, const Node* node) const = 0;

    /** \brief Get number of standard (w/o enrichment) dofs for given node.
     *
     *  \remark In the default case, we can use the NumDof() routine. This method
     *          gives you only for enriched nodes a different result!
     *          ( see DRT::DiscretizationXFEM )
     *
     *  \param dofset_id (in): number of dofset
     *  \param node      (in): the node those number of dofs are requested
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumStandardDof(const unsigned& dofset_id, const Node* node) const = 0;

    /** \brief Get the gid of all dofs of a node.
     *
     *  Ask the current DofSet for the gids of the dofs of this node. The
     *  required vector is created and filled on the fly. So better keep it
     *  if you need more than one dof gid.
     *
     *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
     *
     *  \param dofset_id (in) : number of dofset
     *  \param node (in)      : the node */
    std::vector<int> Dof(const Node* node) const
    {
      dsassert(NumDofSets() == 1, "expect just one dof set");
      return Dof(0, node);
    };
    virtual std::vector<int> Dof(unsigned dofset_id, const Node* node) const = 0;

    /** \brief Get the gid of a dof for given node.
     *
     *  Ask the current DofSet for the gid of the dof of this node.
     *  There is a variable number of DofSets and one is currently selected.
     *
     *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
     *
     *  \param dofset_id (in): number of dofset
     *  \param node      (in): the node
     *  \param dof       (in): the node local dof number */
    int Dof(const Node* node, const int ldof) const
    {
      dsassert(NumDofSets() == 1, "expect just one dof set");
      return Dof(0, node, ldof);
    };
    virtual int Dof(unsigned dofset_id, const Node* node, const int ldof) const = 0;

    /** \brief Get the gid of all dofs of a node.
     *
     *  Ask the current DofSet for the gids of the dofs of this node. The
     *  required vector is created and filled on the fly. So better keep it
     *  if you need more than one dof gid.
     *
     *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
     *
     *  \param dofset_id (in): number of dofset
     *  \param node      (in): the node
     *  \param lm    (in/out): lm vector the dofs are appended to */
    void Dof(const Node* node, std::vector<int>& lm) const
    {
      dsassert(NumDofSets() == 1, "expect just one dof set");
      Dof((unsigned)0, node, lm);
    };
    virtual void Dof(unsigned dofset_id, const Node* node, std::vector<int>& lm) const = 0;

    /*! \brief Get the gid of all dofs of a node.
     *
     *  Ask the current DofSet for the GID's of the DoF's of this node. The
     *  required vector is created and filled on the fly. So better keep it
     *  if you need more than one DoF GID.
     *  - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))
     *
     *  Additional input nodal DoF set: If the node contains more than one set
     *  of DoF's, which can be evaluated, the number of the set needs to be
     *  given. Currently only the case for XFEM.
     *
     *  \param dof           (out): vector of dof gids (to be filled)
     *  \param node           (in): the node
     *  \param dofset_id      (in): id of the desired dofset
     *  \param nodaldofset_id (in): id of nodal dofset (currently only for XFEM)
     *  \param element        (in): the element (optional)
     */
    virtual void Dof(std::vector<int>& dof, const Node* node, unsigned dofset_id,
        unsigned nodaldofset_id, const Element* element = NULL) const = 0;

    /** \brief Get all conditions with a certain name
     *
     *  Get all conditions with a certain name. A vector of ptrs to all conditions
     *  with name name is returned in out. The number of conditions found with name
     *  name is out.size(). out.size() is 0 if no condition with that name is found.
     *
     *  \note Conditions attached to the discretization have to be
     *        completely redundant meaning that nodal cloud in the
     *        condition is the same on each processor and spans all
     *        nodes that hold this condition
     *
     *  \param name ( in): Name of condition
     *  \param out  (out): vector of pointers to all conditions with that name
     *
     *  \return Returns out.size()=0 if condition with that name does not exist */
    virtual void GetCondition(const std::string& name, std::vector<DRT::Condition*>& out) const = 0;
    virtual void GetCondition(
        const std::string& name, std::vector<Teuchos::RCP<Condition>>& out) const = 0;

    /** \brief Get a condition with a certain name
     *
     *  Returns the first condition with name name found in the multimap.
     *  If multiple conditions with the same name exist, the first condition is
     *  returned and behavior is therefore non-deterministic. This method should
     *  therefore only be used in cases where the user is sure that name is unique.
     *
     *  \note Conditions attached to the discretization have to be
     *     completely redundant meaning that nodal cloud in the
     *     condition is the same on each processor and spans all
     *     nodes that hold this condition
     *
     *  \param name (in): Name of condition
     *
     *  \return Returns NULL if condition with that name does not exist */
    virtual DRT::Condition* GetCondition(const std::string& name) const = 0;

    /** \brief Set a reference to a data vector
     *
     *  Using this method, a reference to a vector can
     *  be supplied to the discretization. The elements can access
     *  this vector by using the name of that vector.
     *  The method expects state to be either of dof row map or of
     *  dof column map.
     *  If the vector is supplied in DofColMap() a reference to it will be stored.
     *  If the vector is NOT supplied in DofColMap(), but in DofRowMap(),
     *  a vector with column map is allocated and the supplied vector is exported to it.
     *  Everything is stored/referenced using Teuchos::RCP.
     *
     *  \param nds   (in): number of dofset
     *  \param name  (in): Name of data
     *  \param state (in): vector of some data
     *
     *  \note This class will not take ownership or in any way modify the solution vector. */
    void SetState(const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
    {
      SetState(0, name, state);
    }
    virtual void SetState(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state) = 0;

    /** \brief Get a reference to a data vector
     *
     *  Providing a name of a solution state, get a reference to the solution vector.
     *  If a vector under the provided name does not exist, the method will throw
     *  a lethal error message.
     *
     *  \param nds (in): number of dofset
     *  \param name (in): Name of solution state
     *
     *  \return Reference to solution state */
    Teuchos::RCP<const Epetra_Vector> GetState(const std::string& name) const
    {
      return GetState(0, name);
    };
    virtual Teuchos::RCP<const Epetra_Vector> GetState(
        unsigned nds, const std::string& name) const = 0;

    /** \brief Clear solution state references
     *
     *  \param clearalldofsets (in): If FALSE, only the state belonging to the
     *                               own dofset will be deleted.
     *
     *  The method deletes references to any solution data */
    virtual void ClearState(bool clearalldofsets = false) = 0;

    /** \brief Evaluate Neumann boundary conditions
     *
     *  Loop all Neumann conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  "PointNeumann", LineNeumann", "SurfaceNeumann" and "VolumeNeumann".
     *  It takes a current time from the parameter list params named "total time"
     *  and evaluates the appropiate time curves at that time for each
     *  Neumann condition separately. If "total time" is not included
     *  in the parameters, no time curves are used.
     *  Parameters recognized by this method:
     *  \code
     *  params.set("total time",acttime); // current total time
     *  \endcode
     *
     *  \param params        (in): List of parameters
     *  \param systemvector (out): Vector to assemble Neumann BCs to.
     *                             The vector is NOT initialized to zero by
     *                             this method.
     *  \param systemvector (out): Optional matrix to assemble linearization
     *                             of Neumann BCs to. */
    void EvaluateNeumann(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix = Teuchos::null)
    {
      if (systemmatrix.is_null())
        EvaluateNeumann(params, *systemvector);
      else
        EvaluateNeumann(params, *systemvector, systemmatrix.get());
      return;
    }
    virtual void EvaluateNeumann(Teuchos::ParameterList& params, Epetra_Vector& systemvector,
        LINALG::SparseOperator* systemmatrix = NULL) = 0;

    /** \brief Call elements to evaluate
     *
     *  Call element routines to perform integration and return element contributions to
     *  system vectors and matrices. Type of action taken by the elements is
     *  controlled by the params parameter.<br>
     *  Parameters that control element behavior are:<br>
     *  \code
     *  params.set("action","<element_action>"); // <element_action> somthing that elements
     * understand \endcode Other parameters eventually recognized by the elements:<br> \code
     *  params.set("total time",1.23);     // current total time of simulation
     *  params.set("delta time",0.01);     // time increment
     *  \endcode
     *
     *
     *  \param params         (in): Parameter list past to the elements containing
     *                              commands and parameters for the elements and
     *                              containing assembly intructions
     *  \param systemmatrix1 (out): Sparse matrix that may be filled by
     *                              assembly of element contributions.
     *                              May be Teuchos::null on entry.
     *                              Matrix must be systemmatrix1->Filled()==false on input.
     *  \param systemmatrix2 (out): Sparse matrix that may be filled by
     *                              assembly of element contributions.
     *                              May be Teuchos::null on entry.
     *                              Matrix must be systemmatrix2->Filled()==false on input.
     *  \param systemvector1 (out): Distributed vector that may be filled by
     *                              assembly of element contributions.
     *                              May be Teuchos::null on entry.
     *                              Vector will NOT be initialized to zero by
     *                              the underlying assembly methods that add element
     *                              contributions.
     *  \param systemvector2 (out): Distributed vector that may be filled by
     *                              assembly of element contributions.
     *                              May be Teuchos::null on entry.
     *                              Vector will NOT be initialized to zero by
     *                              the underlying assembly methods that add element
     *                              contributions.
     *  \param systemvector3 (out): Distributed vector that may be filled by
     *                              assembly of element contributions.
     *                              May be Teuchos::null on entry.
     *                              Vector will NOT be initialized to zero by
     *                              the underlying assembly methods that add element
     *                              contributions. */
    virtual void Evaluate(Teuchos::ParameterList& params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3) = 0;

    /** \brief Evaluate Dirichlet boundary conditions
     *
     *  Loop all Dirichlet conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  "PointDirichlet", "LineDirichlet", "SurfaceDirichlet" and "VolumeDirichlet".
     *  It takes a current time from the parameter list params named "total time"
     *  and evaluates the appropiate time curves at that time for each
     *  Dirichlet condition separately. If "total time" is not included
     *  in the parameters, no time curves are used.
     *
     *  \note Opposed to the other 'Evaluate' method does this one NOT assembly but
     *        OVERWRITE values in the output vector systemvector. For this reason,
     *        dirichlet boundary conditions are evaluated in the following order:
     *        First "VolumeDirichlet", then "SurfaceDirichlet", then "LineDirichlet"
     *        and finally "PointDirichlet". This way, the lower entity dirichlet BCs override
     *        the higher ones and a point Dirichlet BCs has priority over other dirichlet
     *        BCs in the input file.
     *
     *  Parameters recognized by this method:
     *  \code
     *  params.set("total time",acttime); // current total time
     *  \endcode
     *
     *  \param params           (in): List of parameters
     *  \param systemvector    (out): Vector holding prescribed Dirichlet values
     *  \param systemvectord   (out): Vector holding 1st time derivative of
     *                                prescribed Dirichlet values
     *  \param systemvectordd  (out): Vector holding 2nd time derivative prescribed
     *                                Dirichlet values
     *  \param toggle          (out): Vector containing 1.0 for each Dirichlet dof
     *                                and 0 for everything else
     *  \param dbcmapextractor (out): Map extractor containing maps for the DOFs
     *                                subjected to Dirichlet boundary conditions and
     *                                the remaining/free DOFs */
    virtual void EvaluateDirichlet(Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> systemvector, Teuchos::RCP<Epetra_Vector> systemvectord,
        Teuchos::RCP<Epetra_Vector> systemvectordd, Teuchos::RCP<Epetra_Vector> toggle,
        Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor = Teuchos::null) = 0;

    /** \brief Evaluate a specified condition
     *
     *  Loop all conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  matching the user-provided string condstring.
     *  Calls more general EvaluateCondition method, see below.
     *
     *  \param params        (in): List of parameters for use at element level
     *  \param systemvector (out): Vector to assemble BCs to.(NOT initialized to zero
     *                             by this method)
     *  \param condstring    (in): Name of condition to be evaluated
     *  \param condid        (in): condition ID */
    void EvaluateCondition(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_Vector> systemvector,
        const std::string& condstring, const int condid = -1)
    {
      EvaluateCondition(params, Teuchos::null, Teuchos::null, systemvector, Teuchos::null,
          Teuchos::null, condstring, condid);
    };

    /** \brief Evaluate a specified condition
     *
     *  Loop all conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  matching the user-provided string condstring.
     *  Calls more general EvaluateCondition method, see below.
     *
     *  \param params     (in): List of parameters for use at element level
     *  \param condstring (in): Name of condition to be evaluated
     *  \param condid     (in): condition ID */
    void EvaluateCondition(
        Teuchos::ParameterList& params, const std::string& condstring, const int condid = -1)
    {
      EvaluateCondition(params, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null,
          Teuchos::null, condstring, condid);
    }

    /** \brief Evaluate a specific condition
     *
     *  Loop all conditions attached to the discretization and evaluate them.
     *  This method considers all conditions in condition_ with the names
     *  matching the user-provided string condstring.
     *  It takes a current time from the parameter list params named "total time"
     *  and evaluates the appropriate time curves at that time for each
     *  condition separately. If "total time" is not included
     *  in the parameters, no time curves are used.
     *
     *  \param params         (in): List of parameters for use at element level
     *  \param systemmatrix1 (out): Sparse matrix that may be changed by
     *                              assembly of boundary element contributions.
     *                              May not be Teuchos::null.
     *                              Matrix must be systemmatrix->Filled()==false on input.
     *  \param systemmatrix2 (out): Sparse matrix that may be changed by
     *                              assembly of boundary element contributions.
     *                              May not be Teuchos::null.
     *                              Matrix must be systemmatrix->Filled()==false on input.
     *  \param systemvector1 (out): Vector to assemble BCs to.
     *                              The vector is NOT initialized to zero by this method.
     *  \param systemvector2 (out): Vector to assemble BCs to.
     *                              The vector is NOT initialized to zero by this method.
     *  \param systemvector3 (out): Vector to assemble BCs to.
     *                              The vector is NOT initialized to zero by this method.
     *  \param condstring     (in): Name of condition to be evaluated
     *  \param condid         (in): Condition ID */
    virtual void EvaluateCondition(Teuchos::ParameterList& params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring,
        const int condid = -1) = 0;

    /** \brief Assemble scalar quantities across elements
     *
     *  Every element is only called \b once by its owning processor. (We are
     *  parsing the element row map.) At this call the element return its value(s),
     *  i.e. its scalar(s), and contributes to the global value of the respective
     *  scalar quantity(ies).
     *
     *  Example: strain energy in structures
     *
     *  \param params   (in): parameters
     *  \param scalars (out): scalars */
    virtual void EvaluateScalars(
        Teuchos::ParameterList& params, Teuchos::RCP<Epetra_SerialDenseVector> scalars) = 0;

    /** \brief Compute the nullspace of the discretization
     *
     *  This method looks in the solver parameters whether algebraic multigrid (AMG)
     *  is used as preconditioner. AMG desires the nullspace of the system of equations
     *  which is then computed here if it does not already exist in the parameter list.
     *
     *  \param solveparams (in): List of parameters
     *  \param recompute   (in): force method to recompute the nullspace */
    virtual void ComputeNullSpaceIfNecessary(
        Teuchos::ParameterList& solveparams, bool recompute = false) = 0;

  };  // class DiscretizationInterface
}  // namespace DRT


#endif /* SRC_DRT_LIB_DRT_DISCRET_INTERFACE_H_ */
