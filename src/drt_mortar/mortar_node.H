/*-----------------------------------------------------------------------*/
/*! \file
\brief header for mortar node and data container

\level 2

*/
/*-----------------------------------------------------------------------*/
#ifndef MORTAR_NODE_H
#define MORTAR_NODE_H

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../headers/pairedvector.H"

namespace MORTAR
{
  /*!
  \brief A subclass of DRT::ParObjectType that adds mortar element type specific methods

  */
  class MortarNodeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MortarNodeType"; }

    static MortarNodeType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static MortarNodeType instance_;
  };

  /*!
   \brief A class containing additional data from mortar nodes

   This class contains additional information from mortar nodes which are
   are not needed for contact search and therefore are only available on the
   node's processor (ColMap). The class MortarNodeDataContainer must be declared
   before the MortarNode itself.

   */
  class MortarNodeDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    /*!
     \brief Standard Constructor

     */
    MortarNodeDataContainer();

    //! @}

    //! @name Pack/Unpack
    //! @{

    /*!
     \brief Pack this class so that it can be communicated

     This function packs the data container. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void Pack(DRT::PackBuffer& data) const;

    /*!
     \brief Unpack data from a vector into this class

     This function unpacks the data container. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

    //!@}

    //! @name Access methods (general stuff)
    //! @{

    /*!
     \brief Return current nodal normal (only for slave side!) (length 3)
     */
    double* n() { return n_; }

    /*!
     \brief Return current edge tangent(length 3)
     */
    double* EdgeTangent() { return edgeTangent_; }

    //! @}

    //! @name Access methods for Lagrange multipliers
    //! @{

    /*!
     \brief Return current Lagrange multiplier in step n+1 (only for slave side!) (length 3)
     */
    double* lm() { return lm_; }

    /*!
     \brief Return const ptr to current Lagrange multiplier in step n+1
     (only for slave side!) (length 3)
     */
    const double* lm() const { return lm_; }

    /*!
     \brief Return old Lagrange multiplier from step n (only for slave side!) (length 3)
     */
    double* lmold() { return lmold_; }

    /*!
     \brief Return Lagrange multiplier from last Uzawa step (only for slave side!) (length 3)
     */
    double* lmuzawa() { return lmuzawa_; }

    //! @}

    //! @name Access methods for segment-to-segment data
    //! @{

    /*!
     \brief Return the 'D' map (vector) of this node
     */
    GEN::pairedvector<int, double>& GetD() { return drows_; }

    /*!
     \brief Return the 'M' map of this node
     */
    std::map<int, double>& GetM() { return mrows_; }

    //! @}

    //! @name Access methods for node-to-segment data
    //! @{

    /*!
     \brief Return the 'D' map (vector) of this node for node-to-segment
     */
    GEN::pairedvector<int, double>& GetDnts() { return drows_nts_; }

    /*!
     \brief Return the 'M' map of this node for node-to-segment
     */
    std::map<int, double>& GetMnts() { return mrows_nts_; }

    //! @}

    //! @name Access methods for line-to-segment data
    //! @{

    /*!
     \brief Return the 'D' map (vector) of this node for line-to-segment
     */
    GEN::pairedvector<int, double>& GetDlts() { return drows_lts_; }

    /*!
     \brief Return the 'M' map of this node for line-to-segment
     */
    std::map<int, double>& GetMlts() { return mrows_lts_; }

    //! @}

    //! @name Access methods for line-to-line data
    //! @{

    /*!
     \brief Return the 'D' map (vector) of this node for line-to-line
     */
    GEN::pairedvector<int, double>& GetDltl() { return drows_ltl_; }

    /*!
     \brief Return the 'M' map of this node for line-to-line
     */
    std::map<int, double>& GetMltl() { return mrows_ltl_; }

    //! @}

    /*!
     \brief Return the 'Mmod' map of this node
     */
    std::map<int, double>& GetMmod() { return mmodrows_; }

    /*!
     \brief Return the d matrix scale factor of this node
     */
    double& GetDscale() { return d_nonsmooth_; }

   private:
    //! don't want = operator
    MortarNodeDataContainer operator=(const MortarNodeDataContainer& old);

    //! don't want copy constructor
    MortarNodeDataContainer(const MortarNodeDataContainer& old);

   protected:
    //! @name Geometry-related vectors
    //! @{

    //! Nodal normal for contact methods
    double n_[3];

    //! Current edge tangent
    double edgeTangent_[3];

    //! @}

    //! @name Lagrange multipliers
    //! @{

    //! Current Lagrange multiplier value (n+1)
    double lm_[3];

    //! Old Lagrange multiplier value (last converged state n)
    double lmold_[3];

    //! Uzawa Lagrange multiplier value (last Uzawa step k)
    double lmuzawa_[3];

    //! @}
    double d_nonsmooth_;

    //! @name Segment-to-segment (sts) specifics
    //! @{

    //! Nodal rows of D matrix for segment-to-segment
    GEN::pairedvector<int, double> drows_;

    //! Nodal rows of M matrix for segment-to-segment
    std::map<int, double> mrows_;

    //! Nodal rows of Mmod matrix (segment-to-segment)
    std::map<int, double> mmodrows_;

    //! @}

    //! @name Node-to-segment (nts) specifics
    //! @{

    //! Nodal rows of D matrix for node-to-segment
    GEN::pairedvector<int, double> drows_nts_;

    //! Nodal rows of M matrix for node-to-segment
    std::map<int, double> mrows_nts_;

    //! @}

    //! @name Line-to-segment (lts) specifics
    //! @{

    //! Nodal rows of D matrix for line-to-segment
    GEN::pairedvector<int, double> drows_lts_;

    //! Nodal rows of M matrix for line-to-segment
    std::map<int, double> mrows_lts_;

    //! @}

    //! @name Line-to-line (ltl) specifics
    //! @{

    //! Nodal rows of D matrix for lts
    GEN::pairedvector<int, double> drows_ltl_;

    //! Nodal rows of M matrix for lts
    std::map<int, double> mrows_ltl_;

    //! @}
  };  // class MortarNodeDataContainer

  /*!
   \brief A class for a mortar node derived from DRT::Node

   This class represents a finite element node capable of mortar coupling.

   */
  class MortarNode : public DRT::Node
  {
   public:
    //! @name Enums and Friends
    //! @{

    /*!
     \brief The Discretization is a friend of MortarNode
     */
    friend class DRT::Discretization;

    //! @}

    //! @name Constructors and destructors and related methods

    /*!
     \brief Standard Constructor

     \param id     (in): A globally unique node id
     \param coords (in): vector of nodal coordinates, length 3
     \param owner  (in): Owner of this node.
     \param numdof (in): Number of degrees of freedom to this node
     \param dofs   (in): list of global degrees of freedom
     \param isslave(in): flag indicating whether node is slave or master

     */
    MortarNode(int id, const double* coords, const int owner, const int numdof,
        const std::vector<int>& dofs, const bool isslave);

    /*!
     \brief Copy Constructor

     Makes a deep copy of a MortarNode

     */
    MortarNode(const MORTAR::MortarNode& old);

    /*!
     \brief Deep copy the derived class and return pointer to it

     */
    virtual MORTAR::MortarNode* Clone() const;

    /*!
     \brief Destructor

     */
    virtual ~MortarNode() { return; }

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_lib/drt_parobject.H.

     */
    int UniqueParObjectId() const override
    {
      return MortarNodeType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     \ref Pack and \ref Unpack are used to communicate this node

     */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
     \brief Unpack data from a char vector into this class

     \ref Pack and \ref Unpack are used to communicate this node

     */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! @name Access methods

    /*!
     \brief Print this MortarNode
     */
    virtual void Print(std::ostream& os) const;

    /*!
     \brief Is node on slave or master side of mortar interface
     */
    bool IsSlave() const { return isslave_; }

    /*!
     \brief Modify slave / master status of current node

     This changing of topology becomes necessary for self contact
     simulations, where slave and master status are assigned dynamically

     This belated modification is also necessary to be able to deal with
     boundary nodes on the slave side of the interface. Their status is
     changed to master, thy do NOT carry Lagrange multipliers and their
     neighbors' dual shape fct. are modified!
     */
    bool& SetSlave() { return isslave_; }

    /*!
     \brief Return attached status

     */
    bool IsDetected() { return detected_; }

    /*!
     \brief Return attached status

     */
    bool& SetDetected() { return detected_; }

    /*!
     \brief Is slave node tied or untied
     */
    bool IsTiedSlave() const { return istiedslave_; }

    /*!
     \brief Modify tying status of current slave node

     This changing of status becomes necessary for meshtying simulations
     where the given slave surface only partially overlaps with the master
     surface. Then the flag istiedslave_ needs to be initialized according
     to the actual tying status. True means that the node is participating
     in meshtying and thus carries mortar contributions. False means that
     the node is not involved in meshtying and does not need to carry
     Lagrange multipliers anyway. There is some similarity with an active
     set definition in contact mechanics, yet the set is static here in the
     meshtying case, of course.

     During problem initialization this flag is first set to the same value
     as isslave_, i.e. we assume a complete projection of the slave surface.
     Then, the actual meshtying zone is identified and this flag is adapted
     for each node accordingly.

     */
    bool& SetTiedSlave() { return istiedslave_; }

    /*!
     \brief Is node on boundary of slave side of mortar interface
     */
    bool IsOnBound() const { return isonbound_; }

    /*!
     \brief Set slave side boundary status of current node
     */
    bool& SetBound() { return isonbound_; }

    /*!
     \brief Is this node on a boundary, edge, or corner?
     */
    bool IsOnBoundorCE() const { return ((isonedge_ or isoncorner_) or isonbound_); }

    /*!
     \brief Is this node on an edge or corner?
     */
    bool IsOnCornerEdge() const { return (isonedge_ or isoncorner_); }

    /*!
     \brief Is this nnode on a geometrical edge
     */
    bool IsOnEdge() const { return isonedge_; }

    /*!
     \brief Set node is on geometrical edge
     */
    bool& SetOnEdge() { return isonedge_; }

    /*!
     \brief Is this node on geometrical corner?
     */
    bool IsOnCorner() const { return isoncorner_; }

    /*!
     \brief Is this node on a boundary or a corner?
     */
    bool IsOnCornerorBound() const { return (isoncorner_ or isonbound_); }

    /*!
     \brief Set node is on geometrical corner
     */
    bool& SetOnCorner() { return isoncorner_; }

    /*!
     \brief Return D.B.C. status of this node (true if at least one dof with D.B.C)
     */
    bool IsDbc() const { return isdbc_; }

    /*!
     \brief Set D.B.C. status of current node
     */
    bool& SetDbc() { return isdbc_; }

    /*!
     \brief Get number of degrees of freedom
     */
    int NumDof() const { return numdof_; }

    /*!
     \brief Get predefined degrees of freedom
     */
    const int* Dofs() const { return &dofs_[0]; }

    /*!
     \brief Return current configuration (length 3)
     */
    double* xspatial() { return xspatial_; }

    /*!
     \brief Return const ptr to current configuration (length 3)
     */
    const double* xspatial() const { return xspatial_; }

    /*!
     \brief Return old displacement (length 3)
     */
    double* uold() { return uold_; }

    /*!
     \brief Return projection status of this node (only for slave side!)
     */
    bool& HasProj() { return hasproj_; }

    /*!
     \brief Return segmentation / cell status of this node (only for slave side!)
     */
    bool& HasSegment() { return hassegment_; }

    /// get the upper bound for the number of D-entries
    int GetNumDentries() const { return dentries_; };

    /*!
     \brief Return of data container of this node

     This method returns the data container of this mortar node where additional
     contact specific quantities/information are stored.

     */
    inline MORTAR::MortarNodeDataContainer& MoData() const
    {
      if (modata_ == Teuchos::null) dserror("No mortar data attached. (node-id = %d)", Id());
      return *modata_;
    }

    //@}

    //! @name Evaluation methods

    /*!
     \brief Add a value to the 'D' map of this node

     The 'D' map is later assembled to the D matrix.
     Note that drows_ here is a vector.

     \param row : local dof row id to add to (rowwise)
     \param val : value to be added
     \param col : global dof column id of the value added

     */
    void AddDValue(const int& colnode, const double& val);
    void AddDntsValue(const int& colnode, const double& val);
    void AddDltsValue(const int& colnode, const double& val);
    void AddDltlValue(const int& colnode, const double& val);

    /*!
     \brief Add a value to the 'M' map of this node

     The 'M' map is later assembled to the M matrix.
     Note that mrows_ here is a vector.

     \param row : local dof row id to add to (rowwise)
     \param val : value to be added
     \param col : global dof column id of the value added

     */
    void AddMValue(const int& colnode, const double& val);
    void AddMntsValue(const int& colnode, const double& val);
    void AddMltsValue(const int& colnode, const double& val);
    void AddMltlValue(const int& colnode, const double& val);

    /*!
     \brief Add a value to the 'Mmod' map of this node

     The 'Mmod' map is later assembled to the M matrix.
     Note that mmodrows_ here is a vector.

     \param row : local dof row id to add to (rowwise)
     \param val : value to be added
     \param col : global dof column id of the value added

     */
    void AddMmodValue(const int& colnode, const double& val);

    /*!
     \brief Build nodal normal
     */
    virtual void BuildAveragedNormal();

    /*!
     \brief Find closest node from given node set and return pointer

     This method will compute the distance of the active node to all
     nodes of the given Epetra_Map on the given DRT::Discretization

     \param intdis (in):         Master Node to project
     \param nodesearchmap (in):  Slave Celement to project on
     \param mindist (out):       Distance to closest node

     */
    virtual MORTAR::MortarNode* FindClosestNode(const Teuchos::RCP<DRT::Discretization> intdis,
        const Teuchos::RCP<Epetra_Map> nodesearchmap, double& mindist);

    /*!
     \brief Check if mesh re-initialization for this node was feasible

     This method checks whether mesh distortion due to the applied relocation
     of the current node stays below a certain limit. This check is very
     empirical, of course!

     The mesh distortion status is returned via a boolean parameter.
     (TRUE indicates that mesh distortion is acceptable!)

     \param relocation (in):    relocation distance
     \param limit      (in):    limit in percent of edge length

     */
    virtual bool CheckMeshDistortion(double& relocation, const double& limit);

    /*!
     \brief Initializes the data container of the node

     With this function, the container with contact specific quantities/information
     is initialized.

     */
    virtual void InitializeDataContainer();

    virtual void InitializePoroDataContainer()
    {
      dserror("Shouldn't land here...");
      return;
    };
    virtual void InitializeEhlDataContainer()
    {
      dserror("Shouldn't land here...");
      return;
    };

    /*!
     \brief Resets the data container of the node

     With this function, the container with contact specific quantities/information
     is deleted / reset to Teuchos::null pointer

     */
    virtual void ResetDataContainer();

    // return dirichlet satus at node
    bool* DbcDofs() { return dbcdofs_; }

    /*!
     \brief Return weighting for nurbs control point (node)

     */
    inline double& NurbsW() { return nurbsw_; }
    inline double NurbsW() const { return nurbsw_; }

    //@}

   protected:
    //! true if this node is on slave side of mortar interface
    bool isslave_;

    //! true if this node is on slave side and actually tied
    bool istiedslave_;

    //! true if this node is on slave side boundary
    bool isonbound_;

    //! true if this node is on geometrical edge
    bool isonedge_;

    //! true if this node is on geometrical corner
    bool isoncorner_;

    //! true if this node has a Dirichlet boundary condition applied to at least one of its DOFs
    bool isdbc_;

    //! true at each node that has a Dirichlet boundary condition applied at i-th dof
    bool dbcdofs_[3];

    //! number of degrees of freedom (length of dofs_)
    int numdof_;

    //! degrees of freedom of this node
    std::vector<int> dofs_;

    //! position in current configuration
    double xspatial_[3];

    //! old displacement (last converged state)
    double uold_[3];

    //! true if feasible projection within search set (slave nodes only)
    bool hasproj_;

    //! true if some integration segment / cell attached (slave nodes only)
    bool hassegment_;

    //! true if this node is detected in active(slip) maps (only for master side)
    bool detected_;

    //! max number of d matrix entries for this node
    int dentries_;

    //! additional information of proc's mortar nodes
    Teuchos::RCP<MORTAR::MortarNodeDataContainer> modata_;

    //! @name Nurbs specifics
    //! @{

    //! Nurbs control point weighting
    double nurbsw_;

    //! @}

  };  // class MortarNode
}  // namespace MORTAR


//! << operator
std::ostream& operator<<(std::ostream& os, const MORTAR::MortarNode& cnode);


#endif  // #ifndef MORTAR_NODE_H
