/*!----------------------------------------------------------------------
\file mortar_node.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_NODE_H
#define MORTAR_NODE_H

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../headers/pairedvector.H"

/*!
 \brief MORTAR: namespace of the BACI mortar coupling algorithms

 */
namespace MORTAR
{

class MortarNodeType: public DRT::ParObjectType
{
public:

  std::string Name() const
  {
    return "MortarNodeType";
  }

  static MortarNodeType & Instance()
  {
    return instance_;
  }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

private:

  static MortarNodeType instance_;
};

/*!
 \brief A class containing additional data from mortar nodes

 This class contains additional information from mortar nodes which are
 are not needed for contact search and therefore are only available on the
 node's processor (ColMap). The class MortarNodeDataContainer must be declared
 before the MortarNode itself.

 \author gitterle (gitterle@lnm.mw.tum.de)
 */
class MortarNodeDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
   \brief Standard Constructor

   */
  MortarNodeDataContainer();

  /*!
   \brief Destructor

   */
  virtual ~MortarNodeDataContainer()
  {
    return;
  }

  /*!
   \brief Pack this class so that it can be communicated

   This function packs the datacontainer. This is only called
   when the class has been initialized and the pointer to this
   class exists.

   */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
   \brief Unpack data from a vector into this class

   This function unpacks the datacontainer. This is only called
   when the class has been initialized and the pointer to this
   class exists.

   */
  virtual void Unpack(std::vector<char>::size_type& position,
      const std::vector<char>& data);

  //@}

  //! @name Access methods

  /*!
   \brief Return current nodal normal (only for slave side!) (length 3)
   */
  virtual double* n()
  {
    return n_;
  }

  /*!
   \brief Return current Lagrange mutliplier in step n+1 (only for slave side!) (length 3)
   */
  virtual double* lm()
  {
    return lm_;
  }

  /*!
   \brief Return old Lagrange mutliplier from step n (only for slave side!) (length 3)
   */
  virtual double* lmold()
  {
    return lmold_;
  }

  /*!
   \brief Return Lagrange mutliplier from last Uzawa step (only for slave side!) (length 3)
   */
  virtual double* lmuzawa()
  {
    return lmuzawa_;
  }

  /*!
   \brief Return the 'D' map (vector) of this node
   */
  virtual std::vector<GEN::pairedvector<int, double> >& GetD()
  {
    return drows_;
  }

  /*!
   \brief Return the 'M' map (vector) of this node
   */
  virtual std::vector<std::map<int, double> >& GetM()
  {
    return mrows_;
  }

  /*!
   \brief Return the 'Mmod' map (vector) of this node
   */
  virtual std::vector<std::map<int, double> >& GetMmod()
  {
    return mmodrows_;
  }

  /*!
   \brief Return the scale factor of this node
   */
  virtual double& GetScale()
  {
    return scalefac_;
  }

  //@}

protected:

  // don't want = operator and cctor
  MortarNodeDataContainer operator =(const MortarNodeDataContainer& old);
  MortarNodeDataContainer(const MortarNodeDataContainer& old);

  double n_[3]; // nodal normal for contact methods
  double lm_[3]; // current Lagrange multiplier value (n+1)
  double lmold_[3]; // old Lagrange multiplier value (last converged state n)
  double lmuzawa_[3]; // Uzawa Lagrange multiplier value (last Uzawa step k)

  std::vector<GEN::pairedvector<int, double> > drows_; // nodal rows of D matrix
  std::vector<std::map<int, double> > mrows_; // nodal rows of M matrix
  std::vector<std::map<int, double> > mmodrows_; // nodal rows of Mmod matrix

  double scalefac_; // nodal scaling factor for better conditioning

};
// class MortarNodeDataContainer

/*!
 \brief A class for a mortar node derived from DRT::Node

 This class represents a finite element node capable of mortar coupling.

 \author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
 */
class MortarNode: public DRT::Node
{
public:

  //! @name Enums and Friends

  /*!
   \brief The Discretization is a friend of MortarNode
   */
  friend class DRT::Discretization;

  //@}

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
  virtual ~MortarNode()
  {
    return;
  }

  /*!
   \brief Return unique ParObject id

   every class implementing ParObject needs a unique id defined at the
   top of drt_lib/drt_parobject.H.

   */
  virtual int UniqueParObjectId() const
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
  virtual bool IsSlave() const
  {
    return isslave_;
  }

  /*!
   \brief Modify slave / master status of current node

   This changing of topology becomes necessary for self contact
   simulations, where slave and master status are assigned dynamically

   This belated modification is also necessary to be able to deal with
   boundary nodes on the slave side of the interface. Their status is
   changed to master, thy do NOT carry Lagrange multipliers and their
   neighbors' dual shape fct. are modified!
   */
  virtual bool& SetSlave()
  {
    return isslave_;
  }

  /*!
   \brief Return attached status

   */
  virtual bool IsDetected()
  {
    return detected_;
  }

  /*!
   \brief Return attached status

   */
  virtual bool& SetDetected()
  {
    return detected_;
  }

  /*!
   \brief Is slave node tied or untied
   */
  virtual bool IsTiedSlave() const
  {
    return istiedslave_;
  }

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
  virtual bool& SetTiedSlave()
  {
    return istiedslave_;
  }

  /*!
   \brief Is node on boundary of slave side of mortar interface
   */
  virtual bool IsOnBound() const
  {
    return isonbound_;
  }

  /*!
   \brief Set slave side boundary status of current node
   */
  virtual bool& SetBound()
  {
    return isonbound_;
  }

  /*!
   \brief Return D.B.C. status of this node (true if at least one dof with D.B.C)
   */
  virtual bool IsDbc() const
  {
    return isdbc_;
  }

  /*!
   \brief Set D.B.C. status of current node
   */
  virtual bool& SetDbc()
  {
    return isdbc_;
  }

  /*!
   \brief Get number of degrees of freedom
   */
  virtual int NumDof() const
  {
    return numdof_;
  }

  /*!
   \brief Get predefined degrees of freedom
   */
  virtual const int* Dofs() const
  {
    return &dofs_[0];
  }

  /*!
   \brief Return current configuration (length 3)
   */
  virtual double* xspatial()
  {
    return xspatial_;
  }

  /*!
   \brief Return old displacement (length 3)
   */
  virtual double* uold()
  {
    return uold_;
  }

  /*!
   \brief Return projection status of this node (only for slave side!)
   */
  virtual bool& HasProj()
  {
    return hasproj_;
  }

  /*!
   \brief Return segmentation / cell status of this node (only for slave side!)
   */
  virtual bool& HasSegment()
  {
    return hassegment_;
  }

  /*!
   \brief Return of data container of this node

   This method returns the data container of this mortar node where additional
   contact specific quantities/information are stored.

   */
  inline MORTAR::MortarNodeDataContainer& MoData()
  {
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
  void AddDValue(int& row, int& col, double& val);

  /*!
   \brief Add a value to the 'M' map of this node

   The 'M' map is later assembled to the M matrix.
   Note that mrows_ here is a vector.

   \param row : local dof row id to add to (rowwise)
   \param val : value to be added
   \param col : global dof column id of the value added

   */
  void AddMValue(int& row, int& col, double& val);

  /*!
   \brief Add a value to the 'Mmod' map of this node

   The 'Mmod' map is later assembled to the M matrix.
   Note that mmodrows_ here is a vector.

   \param row : local dof row id to add to (rowwise)
   \param val : value to be added
   \param col : global dof column id of the value added

   */
  void AddMmodValue(int& row, int& col, double& val);

  /*!
   \brief Add a value to the scale factor of this node

   \param val : value to be added

   */
  void AddScValue(double& val);

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
  virtual MORTAR::MortarNode* FindClosestNode(
      const Teuchos::RCP<DRT::Discretization> intdis,
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
  virtual bool CheckMeshDistortion(double& relocation, double& limit);

  /*!
   \brief Initializes the data container of the node

   With this function, the container with contact specific quantities/information
   is initialized.

   */
  virtual void InitializeDataContainer();

  virtual void InitializePoroDataContainer() {dserror("Shouldn't land here..."); return;};

  /*!
   \brief Resets the data container of the node

   With this function, the container with contact specific quantities/information
   is deleted / reset to Teuchos::null pointer

   */
  virtual void ResetDataContainer();

  // return dirichlet satus at node
  bool* DbcDofs()
  {
    return dbcdofs_;
  }

  /*!
   \brief Return weighting for nurbs control point (node)

   */
  virtual inline double& NurbsW()
  {
    return nurbsw_;
  }

  //@}

protected:

  bool              isslave_;        // true if cnode is on slave side of contact problem
  bool              istiedslave_;    // true if cnode is on slave side and actually tied
  bool              isonbound_;      // true if cnode is on slave side boundary
  bool              isdbc_;          // true if cnode has D.B.C.applied to at least one dof
  bool              dbcdofs_[3];      // true at each cnode has D.B.C applied at i-th dof
  int               numdof_;         // number of degrees of freedom (length of dofs_)
  std::vector<int>  dofs_;           // degrees of freedom of this node
  double            xspatial_[3];    // current configuration
  double            uold_[3];        // old displacement (last converged state)
  bool              hasproj_;        // true if feasible projection within search set (slave nodes only)
  bool              hassegment_;     // true if some integration segment / cell attached (slave nodes only)
  bool              detected_;       // true if cnode is detected in active(slip) maps (only for master side)
  int               dentries_;       // max number of d matrix entries for this node

  Teuchos::RCP<MORTAR::MortarNodeDataContainer> modata_;   // additional information of proc's mortar nodes

  // nurbs specific:
  double            nurbsw_;         // nurbs control point weighting

}; // class MortarNode
} // namespace MORTAR


// << operator
std::ostream& operator << (std::ostream& os, const MORTAR::MortarNode& cnode);


#endif  // #ifndef MORTAR_NODE_H
