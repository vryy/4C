/*!----------------------------------------------------------------------
\file mortar_node.H

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_NODE_H
#define MORTAR_NODE_H


#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

using namespace std;
using namespace Teuchos;


/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

/*!
\brief A class for a mortar node derived from DRT::Node

This class represents a finite element node capable of mortar coupling.

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarNode : public DRT::Node
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
  MortarNode(int id, const double* coords, const int owner,
             const int numdof, const vector<int>& dofs, const bool isslave);

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
  virtual int UniqueParObjectId() const { return ParObject_MortarNode;}

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
  \brief Print this MortarNode
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Is node on slave or master side of mortar interface
  */
  virtual bool IsSlave() const { return isslave_; }

  /*!
  \brief Modify slave / master status of current node

  This changing of topology becomes necessary for self contact
  simulations, where slave and master status are assigned dynamically
  
  This belated modification is also necessary to be able to deal with
  boundary nodes on the slave side of the interface. Their status is
  changed to master, thy do NOT carry Lagrange multipliers and their
  neighbors' dual shape fct. are modified!
  */
  virtual bool& SetSlave() { return isslave_; }

  /*!
  \brief Is node on boundary of slave side of mortar interface
  */
  virtual bool IsOnBound() const { return isonbound_; }

  /*!
  \brief Set slave side boundary status of current node
  */
  virtual bool& SetBound() { return isonbound_; }

  /*!
  \brief Return D.B.C. status of this node (true if at least one dof with D.B.C)
  */
  virtual bool IsDbc() const { return isdbc_; }
  
  /*!
  \brief Set D.B.C. status of current node
  */
  virtual bool& SetDbc() { return isdbc_; }

  /*!
  \brief Get number of degrees of freedom
  */
  virtual int NumDof() const { return numdof_; }

  /*!
  \brief Get predefined degrees of freedom
  */
  virtual const int* Dofs() const { return &dofs_[0]; }

  /*!
  \brief Return current configuration (length 3)
  */
  virtual double* xspatial() { return xspatial_; }

  /*!
  \brief Return old displacement (length 3)
  */
  virtual double* uold() { return uold_; }

  /*!
  \brief Return current nodal normal (only for slave side!) (length 3)
  */
  virtual double* n() { return n_; }
  
  /*!
  \brief Return current Lagrange mutliplier in step n+1 (only for slave side!) (length 3)
  */
  virtual double* lm() { return lm_; }

  /*!
  \brief Return old Lagrange mutliplier from step n (only for slave side!) (length 3)
  */
  virtual double* lmold() { return lmold_; }
  
  /*!
  \brief Return Lagrange mutliplier from last Uzawa step (only for slave side!) (length 3)
  */
  virtual double* lmuzawa() { return lmuzawa_; }

  /*!
  \brief Return current closest node's global id (by reference)
  */
  virtual int& ClosestNode() { return closestnode_; }

  /*!
  \brief Return projection status of this node (only for slave side!)
  */
  virtual bool& HasProj() { return hasproj_; }
    
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
  \brief Return the 'D' map (vector) of this node
  */
  virtual vector<map<int,double> >& GetD() { return drows_; }

  /*!
  \brief Return the 'M' map (vector) of this node
  */
  virtual vector<map<int,double> >& GetM() { return mrows_; }

  /*!
  \brief Return the 'Mmod' map (vector) of this node
  */
  virtual vector<map<int,double> >& GetMmod() { return mmodrows_; }
    
  //@}

  //! @name Evaluation methods
  
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
  virtual MORTAR::MortarNode* FindClosestNode(const RCP<DRT::Discretization> intdis,
                                          const RCP<Epetra_Map> nodesearchmap,
                                          double& mindist);
  
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

  //@}

protected:

  bool            isslave_;        // true if cnode is on slave side of contact problem
  bool            isonbound_;      // true if cnode is on slave side boundary
  bool            isdbc_;          // true if cnode has D.B.C.applied to at least one dof
  int             numdof_;         // number of degrees of freedom (length of dofs_)
  vector<int>     dofs_;           // degrees of freedom of this node
  
  double          n_[3];           // nodal normal for contact methods
  double          xspatial_[3];    // current configuration
  double          uold_[3];        // old displacement (last converged state)
  double          lm_[3];          // current Lagrange multiplier value (n+1)
  double          lmold_[3];       // old Lagrange multiplier value (last converged state n)
  double          lmuzawa_[3];     // Uzawa Lagrange multiplier value (last Uzawa step k)
  
  int             closestnode_;    // closest node id on the other side of interface
  bool            hasproj_;        // true if feasible projection within search set (slave nodes only)

  vector<map<int,double> >      drows_;      // nodal rows of D matrix
  vector<map<int,double> >      mrows_;      // nodal rows of M matrix
  vector<map<int,double> >      mmodrows_;   // nodal rows of Mmod matrix
  
  // Penalty-related quantities 
  double                        kappa_; // gap-scaling factor

}; // class MortarNode
} // namespace MORTAR


// << operator
ostream& operator << (ostream& os, const MORTAR::MortarNode& cnode);


#endif  // #ifndef MORTAR_NODE_H
#endif  // #ifdef CCADISCRET
