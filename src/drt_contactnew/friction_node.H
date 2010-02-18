/*!----------------------------------------------------------------------
\file friction_node.H

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
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef FRICTION_NODE_H
#define FRICTION_NODE_H

#include "../drt_mortar/mortar_node.H"
#include "contact_node.H"

using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact and meshtying algorithms

*/
namespace CONTACT
{

/*!
\brief A class containing additional data from frictional contact nodes

This class contains additional information from frictional contact nodes which are
are not needed for contact search and therefore are only available on the
node's processor (ColMap). The class FriNodeDataContainer must be declared
before the FriNode itself.

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class FriNodeDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  FriNodeDataContainer();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a FriNodeDataContainer

  */
  FriNodeDataContainer(const CONTACT::FriNodeDataContainer& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::FriNodeDataContainer* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~FriNodeDataContainer() { return; }

  /*!
  \brief Pack this class so that it can be communicated

  This function packs the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a vector into this class

  This function unpacks the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Unpack(int& position, const vector<char>& data);
  
  //@}

  //! @name Access methods
  
  /*!
  \brief Return jump per time step (only for slave side!) (length 3)
  */
  virtual inline double* jump() { return jump_; }

  /*!
  \brief Return contact status of last converged state n (active=true)
  */
  virtual inline bool& ActiveOld() { return activeold_; }

  /*!
  \brief Return friction status of this node (slip=true)
  */
  virtual inline bool& Slip() { return slip_; }
  
  /*!
  \brief Return the old 'D' map (vector) of this node (last converged state)
  */
  virtual inline vector<map<int,double> >& GetDOld() { return drowsold_; }

  /*!
  \brief Return the old 'M' map (vector) of this node (last converged state)
  */
  virtual inline vector<map<int,double> >& GetMOld() { return mrowsold_; }
 
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
  \brief Return the 'DerivM' map of maps of this node (Petrov-Galerkin-approach)

  This map contains the directional derivatives of the node's
  M-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for M_ik.

  */
  virtual inline map<int,map<int,double> >& GetDerivMPG() { return derivmPG_; }

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
  \brief Return the 'DerivJump' map (vector) of this node

  These maps contain the directional derivatives of the node's
  relative movement (jump).
  A vector is used because the jump itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivJump() { return derivjump_; }
  
  /*!
  \brief Return current penalty traction (length 3)
  */
  virtual inline double* traction() { return traction_; }

  /*!
  \brief Return old penalty traction (length 3)
  */
  virtual inline double* tractionold() { return tractionold_; }

  //@}

protected:

  double          jump_[3];     // jump per time step
  bool            activeold_;   // true if cnode was in contact (last converged state n)
  bool            slip_;        // true if cnode slips (contact force reaches friction bound)

  vector<map<int,double> >      drowsold_;     // nodal rows of old D matrix
  vector<map<int,double> >      mrowsold_;     // nodal rows of old M matrix
  vector<map<int,double> >      drowsPG_;      // nodal rows of D matrix with petrov galerkin approach
  vector<map<int,double> >      mrowsPG_;      // nodal rows of M matrix with petrov galerkin approach
  vector<map<int,double> >      drowsoldPG_;   // nodal rows of old D matrix with petrov galerkin approach
  vector<map<int,double> >      mrowsoldPG_;   // nodal rows of old M matrix with petrov galerkin approach
  set<int>                      snodes_;       // nodal set of according slave nodes
  set<int>                      mnodes_;       // nodal set of according master nodes
  set<int>                      mnodesold_;    // nodal set of old according masternodes

  map<int,map<int,double> >     derivdPG_; // directional derivative of nodal D-matrix value
  map<int,map<int,double> >     derivmPG_; // directional derivative of nodal M-matrix values
  vector<map<int,double> >      derivjump_;  // directional derivative of nodal weighted jump vector

  // penalty related quantities
  double                        traction_[3];     // traction vector of current time step
  double                        tractionold_[3];  // traction vector of previous time step

}; // class FriNodeDataContainer


/*!
\brief A class for a frictional contact node derived from CONTACT::CoNode

This class represents a finite element node capable of frictional contact.

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class FriNode : public CONTACT::CoNode
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of FriNode
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
  FriNode(int id, const double* coords, const int owner,
         const int numdof, const vector<int>& dofs, const bool isslave,
         const bool initactive);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a FriNode

  */
  FriNode(const CONTACT::FriNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::FriNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~FriNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.

  */
  virtual int UniqueParObjectId() const { return ParObject_FriNode;}

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

  //! @name access methods

  /*!
  \brief Print this cnode
  */
  virtual void Print(ostream& os) const ;
  
  /*!
  \brief Return of data container of this node

  This method returns the data container of this node where additional
  contact specific quantities/information are stored.

  */
  inline CONTACT::FriNodeDataContainer& Data() { return *fridata_; }

  //@}

  //! @name Evaluation methods
 
  /*!
  \brief Add a value to the SNode set of this node

  */
  void AddSNode (int node);

  /*!
  \brief Add a value to the MNode set of this node

  */
  void AddMNode (int node);
  
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
  \brief Add a value to the 'M' map of this node (Petrov-Galerkin approach)

  The 'M' map is later assembled to the M matrix.
  Note that mrows_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddMValuePG(int& row, int& col, double& val);

  /*!
  \brief Add a value to the map of Jump derivatives of this node

  Note that derivjump_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddDerivJumpValue(int& row, const int& col, double val);
  
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
  virtual void InitializeDataContainer();
  
  //@}

protected:
  
  RCP<CONTACT::FriNodeDataContainer> fridata_;       // additional information of proc's contact nodes
  
}; // class FriNode

} // namespace CONTACT

//// << operator
//ostream& operator << (ostream& os, const CONTACT::FriNode& cnode);

#endif  // #ifndef FRICTION_NODE_H
#endif  // #ifdef CCADISCRET
