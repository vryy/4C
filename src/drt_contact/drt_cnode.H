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
  \brief Return current penalty traction (length 3)
  */
  virtual inline double* traction() { return traction_; }
  
  /*!
  \brief Return old penalty traction (length 3)
  */
  virtual inline double* tractionold() { return tractionold_; }

  /*!
  \brief Return current nodal normal (only for slave side!) (length 3)
  */
  virtual inline double* n() { return n_; }

  /*!
  \brief Return current nodal tangent t_xi (only for slave side!) (length 3)
  */
  virtual inline double* txi() { return txi_; }

  /*!
  \brief Return current nodal tangent t_eta (only for slave side!) (length 3)
  */
  virtual inline double* teta() { return teta_; }

  /*!
  \brief Return current Lagrange mutliplier in step n+1 (only for slave side!) (length 3)
  */
  virtual inline double* lm() { return lm_; }

  /*!
  \brief Return old Lagrange mutliplier from step n (only for slave side!) (length 3)
  */
  virtual inline double* lmold() { return lmold_; }
  
  /*!
  \brief Return Lagrange mutliplier from last Uzawa step (only for slave side!) (length 3)
  */
  virtual inline double* lmuzawa() { return lmuzawa_; }

  /*!
  \brief Return summation of incremental jumps per time step (only for slave side!) (length 3)
  */
  virtual inline double* jump() { return jump_; }

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
  \brief Return contact status of last converged state n (active=true)
  */
  virtual inline bool& ActiveOld() { return activeold_; }

  /*!
  \brief Return friction status of this node (slip=true)
  */
  virtual inline bool& Slip() { return slip_; }

  /*!
  \brief Return scaling factor for weighted gap (penalty strategy)
  */
  virtual inline double& Kappa() { return kappa_; }
    
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
  \brief Return the 'D' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetD() { return drows_; }

  /*!
  \brief Return the 'M' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetM() { return mrows_; }

  /*!
  \brief Return the old 'D' map (vector) of this node (last converged state)
  */
  virtual inline vector<map<int,double> >& GetDOld() { return drowsold_; }

  /*!
  \brief Return the old 'M' map (vector) of this node (last converged state)
  */
  virtual inline vector<map<int,double> >& GetMOld() { return mrowsold_; }

  /*!
  \brief Return the 'Mmod' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetMmod() { return mmodrows_; }

  /*!
  \brief Return the 'D' map (vector) of this node (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetDPG() { return drowsPG_; }

  /*!
  \brief Return the 'M' map (vector) of this node (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetMPG() { return mrowsPG_; }

  /*!
  \brief Return the old 'D' map (vector) of this node (last converged state) (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetDOldPG() { return drowsoldPG_; }

  /*!
  \brief Return the old 'M' map (vector) of this node (last converged state) (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetMOldPG() { return mrowsoldPG_; }

  /*!
  \brief Return the map with according master nodes
  */
  virtual inline set<int>& GetSNodes() { return snodes_; }
  
  /*!
  \brief Return the map with according master nodes
  */
  virtual inline set<int>& GetMNodes() { return mnodes_; }

  /*!
  \brief Return the old map with according master nodes
  */
  virtual inline set<int>& GetMNodesOld() { return mnodesold_; }

  /*!
  \brief Return the weighted gap (scalar) of this node
  */
  virtual inline double& Getg() { return grow_; }

  /*!
  \brief Return the 'DerivN' map (vector) of this node

  These maps contain the directional derivatives of the node's
  averaged normal vector with respect to the slave displacements.
  A vector is used because the normal itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivN() { return derivn_; }

  /*!
  \brief Return the 'DerivTxi' map (vector) of this node

  These maps contain the directional derivatives of the node's
  unit tangent vector t_xi with respect to the slave displacements.
  A vector is used because the tangent itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivTxi() { return derivtxi_; }

  /*!
  \brief Return the 'DerivTeta' map (vector) of this node

  These maps contain the directional derivatives of the node's
  unit tangent vector t_eta with respect to the slave displacements.
  A vector is used because the tangent itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivTeta() { return derivteta_; }

  /*!
  \brief Return the 'DerivD' map of this node

  D-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for D_ik.
  
  Note: In earlier versions this was just a simple map accounting for the
  diagonality of D when using dual shape functions. Due to the extension to
  support arbitrary types of shape functions, this is not possible anymore.

  */
  virtual inline map<int,map<int,double> >& GetDerivD() { return derivd_; }

  /*!
  \brief Return the 'DerivM' map of maps of this node

  This map contains the directional derivatives of the node's
  M-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for M_ik.

  */
  virtual inline map<int,map<int,double> >& GetDerivM() { return derivm_; }

  /*!
  \brief Return one specific 'DerivD' map of this node

  This method returns the map of directional derivatives of one
  specific D-matrix D_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivD(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivd_.find(k);
    if (p==derivd_.end())
     dserror("ERROR: GetDerivD: No map entry existing for given index");
    return derivd_[k];
  }
  
  /*!
  \brief Return one specific 'DerivM' map of this node

  This method returns the map of directional derivatives of one
  specific M-matrix M_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivM(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivm_.find(k);
    if (p==derivm_.end())
     dserror("ERROR: GetDerivM: No map entry existing for given index");
    return derivm_[k];
  }

  /*!
  \brief Return the 'DerivD' map of this node (Petrov-Galerkin-approach)

  D-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for D_ik.

  Note: In earlier versions this was just a simple map accounting for the
  diagonality of D when using dual shape functions. Due to the extension to
  support arbitrary types of shape functions, this is not possible anymore.

  */
  virtual inline map<int,map<int,double> >& GetDerivDPG() { return derivdPG_; }

  /*!
  \brief Return the 'DerivM' map of maps of this node (Petrov-Galerkin-approach)

  This map contains the directional derivatives of the node's
  M-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for M_ik.

  */
  virtual inline map<int,map<int,double> >& GetDerivMPG() { return derivmPG_; }

  /*!
  \brief Return one specific 'DerivD' map of this node (Petrov-Galerkin-approach)

  This method returns the map of directional derivatives of one
  specific D-matrix D_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivDPG(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivdPG_.find(k);
    if (p==derivdPG_.end())
     dserror("ERROR: GetDerivDPG: No map entry existing for given index");
    return derivdPG_[k];
  }

  /*!
  \brief Return one specific 'DerivM' map of this node (Petrov-Galerkin-approach)

  This method returns the map of directional derivatives of one
  specific M-matrix M_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivMPG(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivmPG_.find(k);
    if (p==derivmPG_.end())
     dserror("ERROR: GetDerivM: No map entry existing for given index");
    return derivmPG_[k];
  }


  /*!
  \brief Return the 'DerivG' map of this node

  This map contains the directional derivatives of the node's
  weighted gap entry g~ with respect to the slave/master displacements.

  */
  virtual inline map<int,double>& GetDerivG() { return derivg_; }
  
  /*!
  \brief Return the 'DerivJump' map (vector) of this node

  These maps contain the directional derivatives of the node's
  relative movement (jump).
  A vector is used because the jump itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivJump() { return derivjump_; }

  /*!
  \brief Return the 'DerivZ map of this node
  
  This map contains the directional derivatives of the node's
  lagrange multiplier entries with respect to the slave/master displacements.

  Note: This is only calculated when performing a penalty strategy

   */
  virtual inline vector<map<int,double> >& GetDerivZ() { return derivz_; }
    
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

  //@}

protected:

  bool            isslave_;        // true if cnode is on slave side of contact problem
  bool            initactive_;     // true if node is initialized as active node
  bool            isonbound_;      // true if cnode is on slave side boundary
  bool            isdbc_;          // true if cnode has D.B.C.applied to at least one dof
  int             numdof_;         // number of degrees of freedom (length of dofs_)
  vector<int>     dofs_;           // degrees of freedom of this node
  
  double          n_[3];           // nodal normal for contact methods
  double          txi_[3];         // nodal tangent t_xi for contact methods
  double          teta_[3];        // nodal tangent t_eta for contact methods
  double          xspatial_[3];    // current configuration
  double          uold_[3];        // old displacement (last converged state)
  double          lm_[3];          // current Lagrange multiplier value (n+1)
  double          lmold_[3];       // old Lagrange multiplier value (last converged state n)
  double          lmuzawa_[3];     // Uzawa Lagrange multiplier value (last Uzawa step k)
  double          jump_[3];        // summation of incremental jumps per time step
  
  int             closestnode_;    // closest node id on the other side of interface
  bool            hasproj_;        // true if feasible projection within search set (slave nodes only)
  bool            active_;         // true if cnode is in contact (active set strategy)
  bool            activeold_;      // true if cnode was in contact (last converged state n)
  bool            slip_;           // true if cnode slips (contact force reaches friction bound)

  vector<map<int,double> >      drows_;          // nodal rows of D matrix
  vector<map<int,double> >      mrows_;          // nodal rows of M matrix
  vector<map<int,double> >      drowsPG_;      // nodal rows of D matrix with petrov galerkin approach
  vector<map<int,double> >      mrowsPG_;      // nodal rows of M matrix with petrov galerkin approach
  vector<map<int,double> >      drowsold_;       // nodal rows of old D matrix
  vector<map<int,double> >      mrowsold_;       // nodal rows of old M matrix
  vector<map<int,double> >      drowsoldPG_;   // nodal rows of old D matrix with petrov galerkin approach
  vector<map<int,double> >      mrowsoldPG_;   // nodal rows of old M matrix with petrov galerkin approach
  vector<map<int,double> >      mmodrows_;       // nodal rows of Mmod matrix
  set<int>                      snodes_;         // nodal set of according slave nodes
  set<int>                      mnodes_;         // nodal set of according master nodes
  set<int>                      mnodesold_;      // nodal set of old according masternodes
  double                        grow_;           // nodal entry of g vector

  vector<map<int,double> >      derivn_;     // directional derivative of nodal normal
  vector<map<int,double> >      derivtxi_;   // directional derivative of nodal tangent t_xi
  vector<map<int,double> >      derivteta_;  // directional derivative of nodal tangent t_eta
  map<int,map<int,double> >     derivd_;     // directional derivative of nodal D-matrix value
  map<int,map<int,double> >     derivm_;     // directional derivative of nodal M-matrix values
  map<int,map<int,double> >     derivdPG_; // directional derivative of nodal D-matrix value
  map<int,map<int,double> >     derivmPG_; // directional derivative of nodal M-matrix values
  map<int,double>               derivg_;     // directional derivative of nodal weighted gap value
  vector<map<int,double> >      derivjump_;  // directional derivative of nodal weighted jump vector

  // Penalty-related quantities
  double                        traction_[3];     // traction vector of current time step
  double                        tractionold_[3];  // traction vector of previous time step
  double                        kappa_;           // gap-scaling factor
  vector<map<int,double> >      derivz_;          // direction derivative of nodal z-matrix value

}; // class CNode
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CNode& cnode);




#endif  // #ifndef CNODE_H
#endif  // #ifdef CCADISCRET
