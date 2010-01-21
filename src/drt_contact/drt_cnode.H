/*!----------------------------------------------------------------------
\file drt_cnode.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CNODE_H
#define CNODE_H


#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "drt_cnodedatacontainer.H"

using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

/*!
\brief A class for a contact node derived from DRT::Node

This class represents a finite element node capable of contact.

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CNode : public DRT::Node
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of CNode
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
  \param initactive (in): flag indicating whether initially set to active

  */
  CNode(int id, const double* coords, const int owner,
        const int numdof, const vector<int>& dofs, const bool isslave,
        const bool initactive);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a CNode

  */
  CNode(const CONTACT::CNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~CNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.

  */
  virtual int UniqueParObjectId() const { return ParObject_CNode;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Print this cnode
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Is Node on slave or master side of contact interface
  */
  inline virtual bool IsSlave() const { return isslave_; }

  /*!
  \brief Modify slave / master status of current node

  This changing of contact topology becomes necessary for self contact
  simulations, where slave and master status are assigned dynamically

  This belated modification is also necessary to be able to deal with
  boundary nodes on the slave side of the interface. Their status is
  changed to master, thy do NOT carry Lagrange multipliers and their
  neighbors' dual shape fct. are modified!
  */
  inline virtual bool& SetSlave() { return isslave_; }

  /*!
  \brief Is Node initialized as active node (only slave nodes)
  */
  inline virtual bool IsInitActive() const
  { if (!IsSlave()) dserror("ERROR: InitActive requested for Master node");
    return initactive_; }

  /*!
  \brief Modify initial active status of slave node

  This belated modification is necessary to be able to use
  the binary search tree for contact initialization in the
  load-controlled quasistatic case (instead of input file
  information Active/Inactive)

  */
  inline virtual bool& SetInitActive()
  { if (!IsSlave()) dserror("ERROR: InitActive requested for Master node");
    return initactive_; }

  /*!
  \brief Is Node on boundary of slave side of contact interface
  */
  inline virtual bool IsOnBound() const { return isonbound_; }

  /*!
  \brief Set slave side boundary status of current node
  */
  inline virtual bool& SetBound() { return isonbound_; }

  /*!
  \brief Return D.B.C. status of this node (true if at least one dof with D.B.C)
  */
  inline virtual bool IsDbc() const { return isdbc_; }

  /*!
  \brief Set D.B.C. status of current node
  */
  inline virtual bool& SetDbc() { return isdbc_; }

  /*!
  \brief Get number of degrees of freedom
  */
  inline virtual int NumDof() const { return numdof_; }

  /*!
  \brief Get predefined degrees of freedom
  */
  inline virtual const int* Dofs() const { return &dofs_[0]; }

  /*!
  \brief Return current configuration (length 3)
  */
  virtual inline double* xspatial() { return xspatial_; }

  /*!
  \brief Return old displacement (length 3)
  */
  virtual inline double* uold() { return uold_; }

  /*!
  \brief Return current closest node's global id (by reference)
  */
  virtual inline int& ClosestNode() { return closestnode_; }

  /*!
  \brief Return projection status of this node (only for slave side!)
  */
  virtual inline bool& HasProj() { return hasproj_; }

  /*!
  \brief Return contact status of this node (active=true)
  */
  virtual inline bool& Active() { return active_; }

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
  \brief Add a value to the 'D' map of this node (Petrov-Galerkin approach)

  The 'D' map is later assembled to the D matrix.
  Note that drows_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddDValuePG(int& row, int& col, double& val);

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
  \brief Add a value to the 'M' map of this node (Petrov-Galerkin approach)

  The 'M' map is later assembled to the M matrix.
  Note that mrows_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddMValuePG(int& row, int& col, double& val);

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
  \brief Add a value to the SNode set of this node

  */
  void AddSNode (int node);

  /*!
  \brief Add a value to the MNode set of this node

  */
  void AddMNode (int node);

  /*!
  \brief Add a value to the weighted gap of this node

  This value is later assembled to the weighted gap vec.
  Note that grow_ here is a scalar.

  \param val : value to be added

  */
  void AddgValue(double& val);

  /*!
  \brief Add a value to the map of Jump derivatives of this node

  Note that derivjump_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddDerivJumpValue(int& row, const int& col, double val);

  /*!
  \brief Add a value to the map of LM derivatives of this node

  The 'DerivZ' map is later assembled to the global DerivZ matrix.
  Note that derivz_ here is a vector.

  Note: This is only calculated when performing a penalty strategy

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddDerivZValue(int& row, const int& col, double val);

  /*!
  \brief Write nodal entries of D and M to Dold and Mold

  At the end of a time step the nodal entries (vector) of the Mortar
  Matrices D and M are stored to the old ones.

  */
  void StoreDMOld();

  /*!
  \brief Write nodal entries of D and M to Dold and Mold (PG approach)

  At the end of a time step the nodal entries (vector) of the Mortar
  Matrices D and M are stored to the old ones.

  */
  void StoreDMOldPG();


  /*!
  \brief Write nodal entries of Penalty tractions to old ones

  At the end of a time step the nodal entries (vector) of the Penalty
  tractions are stored to the old ones.

  */
  void StoreTracOld();

  /*!
  \brief Initializes the data container of the node

  With this function, the container with contact specific quantities/information
  is initialized.

  */
  void InitializeCnodeDataContainer();

  //@}

  //! @name Construction methods

  /*!
  \brief Build nodal normal
  */
  virtual void BuildAveragedNormal();

  /*!
  \brief Build directional derivative of nodal normal + tangents

  This method will be called after having finished the method
  BuildAveragedNormal() for the computation of nodal normals. The result
  (directional derivatives) will be stored in the nodal maps derivn_ and
  will later be assembled to the global system of equations.
  Please note that we also compute the directional derivative of the
  nodal tangent here,and the results will be stored analogously in
  the nodal maps derivtxi_ (2D and 3D) and derivteta_ (only 3D).

  \param elens (in):  Matrix containing normals of adjacent elements
  \param length (in): Length of the nodal averaged normal
  \param ltxi (in):   Length of the nodal tangent txi

  */
  virtual void DerivAveragedNormal(Epetra_SerialDenseMatrix& elens,
                                   double length, double ltxi);

  /*!
  \brief Find closest node from given node set and return global ID

  This method will compute the distance of the active CNode to all
  nodes of the given Epetra_Map on the given DRT::Discretization

  \param intdis (in):         Master Node to project
  \param nodesearchmap (in) : Slave Celement to project on
  \param xi (out) : Local coordinates of projection on Celement

  */
  virtual CONTACT::CNode* FindClosestNode(const RCP<DRT::Discretization> intdis,
                                          const RCP<Epetra_Map> nodesearchmap,
                                          double& mindist);

  /*!
  \brief Return of data container of this node

  This method returns the data container of this node where additional
  contact specific quantities/information are stored.

  */
  virtual inline CONTACT::CNodeDataContainer& GetData() { return *data_; }

  //@}

protected:

  bool            isslave_;        // true if cnode is on slave side of contact problem
  bool            initactive_;     // true if node is initialized as active node
  bool            isonbound_;      // true if cnode is on slave side boundary
  bool            isdbc_;          // true if cnode has D.B.C.applied to at least one dof
  int             numdof_;         // number of degrees of freedom (length of dofs_)
  vector<int>     dofs_;           // degrees of freedom of this node

  double          xspatial_[3];    // current configuration
  double          uold_[3];        // old displacement (last converged state)

  int             closestnode_;    // closest node id on the other side of interface
  bool            hasproj_;        // true if feasible projection within search set (slave nodes only)
  bool            active_;         // true if cnode is in contact (active set strategy)

  RCP<CONTACT::CNodeDataContainer> data_;         // additional information of proc's contact nodes

}; // class CNode
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CNode& cnode);




#endif  // #ifndef CNODE_H
#endif  // #ifdef CCADISCRET
