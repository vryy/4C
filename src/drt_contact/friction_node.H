/*!----------------------------------------------------------------------
\file friction_node.H
\brief A class for a frictional contact node

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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef FRICTION_NODE_H
#define FRICTION_NODE_H

#include "contact_node.H"


/*!
\brief CONTACT: namespace of the baci contact and meshtying algorithms

*/
namespace CONTACT
{

class FriNodeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "FriNodeType"; }

  static FriNodeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static FriNodeType instance_;
};

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
  \brief Destructor

  */
  virtual ~FriNodeDataContainer() { return; }

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
  \brief Return jump per time step (only for slave side!) (length 3)
  */
  virtual inline double* jump() { return jump_; }

  /*!
  \brief Return jump per time step (only for slave side!) (max length 2)
  */
  virtual inline double* jump_var() { return jumpvar_; }
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
  virtual inline std::vector<std::map<int,double> >& GetDOld() { return drowsold_; }

  /*!
  \brief Return the old 'M' map (vector) of this node (last converged state)
  */
  virtual inline std::vector<std::map<int,double> >& GetMOld() { return mrowsold_; }

  /*!
  \brief Return the map with according slave nodes
  */
  virtual inline std::set<int>& GetSNodes() { return snodes_; }

  /*!
  \brief Return the map with according master nodes
  */
  virtual inline std::set<int>& GetMNodes() { return mnodes_; }

  /*!
  \brief Return the old map with according master nodes
  */
  virtual inline std::set<int>& GetMNodesOld() { return mnodesold_; }

   /*!
  \brief Return the 'DerivJump' map (vector) of this node

  These maps contain the directional derivatives of the node's
  relative movement (jump).
  A vector is used because the jump itself is a vector (2 or 3 components).

  */
  virtual inline std::vector<std::map<int,double> >& GetDerivJump() { return derivjump_; }

  /*!
 \brief Return the 'DerivVarJump' map (vector) of this node

 These maps contain the directional derivatives of the node's
 relative -- object variant-- movement (jump). This Jump is already
 multiplied and linearized with the tangent vectors txi and teta.
 Therefore, the vector has 1/2 components instead of 2/3

 */
  virtual inline std::vector<std::map<int,double> >& GetDerivVarJump() { return derivvarjump_; }
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

  // don't want = operator and cctor
  FriNodeDataContainer operator = (const FriNodeDataContainer& old);
  FriNodeDataContainer(const FriNodeDataContainer& old);

  double          jump_[3];     // jump per time step
  double          jumpvar_[2];  // variant jump per time step
  bool            activeold_;   // true if cnode was in contact (last converged state n)
  bool            slip_;        // true if cnode slips (contact force reaches friction bound)

  std::vector<std::map<int,double> > drowsold_;     // nodal rows of old D matrix
  std::vector<std::map<int,double> > mrowsold_;     // nodal rows of old M matrix
  std::set<int>                      snodes_;       // nodal set of according slave nodes
  std::set<int>                      mnodes_;       // nodal set of according master nodes
  std::set<int>                      mnodesold_;    // nodal set of old according masternodes

  std::vector<std::map<int,double> > derivjump_;    // directional derivative of nodal weighted jump vector
  std::vector<std::map<int,double> > derivvarjump_; // directional derivative of obj.-variant. nodal weighted jump vector

  // penalty related quantities
  double                             traction_[3];  // traction vector of current time step
  double                          tractionold_[3];  // traction vector of previous time step

}; // class FriNodeDataContainer

/*!
\brief An additional container for wear and tsi with contact specific data

This class contains additional information to the data container of the 
frictional node. These additional data are needed only for contact 
problems with wear and thermo-structure-interaction problems with contact.

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class FriNodeDataContainerPlus
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  FriNodeDataContainerPlus();

  /*!
  \brief Destructor

  */
  virtual ~FriNodeDataContainerPlus() { return; }

  //@}

  //! @name Access methods

  /*!
  \brief Return the old 'A' map (vector) of this node (last converged state)
  */
  virtual inline std::vector<std::map<int,double> >& GetA() { return arows_; }
  
  /*!
  \brief Return the map with according slave nodes
  */
  virtual inline std::set<int>& GetANodes() { return anodes_; }

  /*!
  \brief Return the weighted wear per node (length 1)
   */
  virtual double& Wear() { return wear_; }
  
  /*!
  \brief Return the delta of the weighted wear per node without wear coefficient (length 1)
   */
   virtual double& DeltaWear() { return deltawear_; }

   /*!
   \brief Return the 'T' map (vector) of this node
   */
   virtual std::vector<std::map<int,double> >& GetT() { return trows_; }

   /*!
   \brief Return the 'E' map (vector) of this node
   */
   virtual std::vector<std::map<int,double> >& GetE() { return erows_; }

   /*!
   \brief Deriv. w.r.t. displ. of E and T matrix entries for
          this node. This matrices are required for the
          discrete wear condition.
   */
   virtual std::map<int,std::map<int,double> >& GetDerivTw() { return derivt_; }
   virtual std::map<int,std::map<int,double> >& GetDerivE() { return derive_; }

   /*!
   \brief Return current discrete wear in step n+1 (only for slave side!) (length 1)
   */
   virtual double* wcurr() { return wcurr_; }
  //@}

protected:

  // don't want = operator and cctor
  FriNodeDataContainerPlus operator = (const FriNodeDataContainerPlus& old);
  FriNodeDataContainerPlus(const FriNodeDataContainerPlus& old);

  // tsi with contact related quantities
  std::vector<std::map<int,double> >      arows_;        // nodal rows of matrix A
  std::set<int>                           anodes_;       // nodal set of according slave nodes
  
  // wear related quantities
  double                                  wear_;         // weighted wear
  double                                  deltawear_;    // deltawear
  double                                  wcurr_[1];        // current discrete wear value (n+1)

  std::vector<std::map<int,double> >      trows_;        // nodal rows of T matrix
  std::vector<std::map<int,double> >      erows_;        // nodal rows of E matrix
  
  std::map<int,std::map<int,double> >     derivt_;        // directional derivative of nodal Tw-matrix values
  std::map<int,std::map<int,double> >     derive_;        // directional derivative of nodal E-matrix values


}; // class FriNodeDataContainerPlus

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
         const int numdof, const std::vector<int>& dofs, const bool isslave,
         const bool initactive, const bool friplus);

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
  virtual int UniqueParObjectId() const { return FriNodeType::Instance().UniqueParObjectId(); }

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

  //! @name access methods

  /*!
  \brief Print this cnode
  */
  virtual void Print(std::ostream& os) const ;

  /*!
  \brief Return mechanical dissipation per node  (length 1)
  */
  virtual double& MechDiss() { return mechdiss_; }
  
  /*!
  \brief Return the 'B' map (vector) of this node 
    */
  virtual inline std::vector<std::map<int,double> >& GetB() { return brows_; }
  
  /*!
  \brief Return the map with according master nodes
  */
  virtual inline std::set<int>& GetBNodes() { return bnodes_; }

  /*!
  \brief Return of data container of this node

  This method returns the data container of this node where additional
  contact specific quantities/information are stored.

  */
  inline CONTACT::FriNodeDataContainer& FriData() { return *fridata_; }

  /*!
  \brief Return of additional data container of this node

  This method returns the additional data container of this node where additional
  tsi and wear specific quantities/information are stored.

  */
  inline CONTACT::FriNodeDataContainerPlus& FriDataPlus() { return *fridataplus_; }
  
  //@}

  //! @name Evaluation methods

  /*!
  \brief Add a value to the SNode set of this node

  */
  void AddSNode (int node);

  /*!
  \brief Add a value to the 'T' map of this node

  The 'T' map is later assembled to the T matrix.
  Note that drows_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddTValue(int& row, int& col, double& val);

  /*!
  \brief Add a value to the 'E' map of this node

  The 'E' map is later assembled to the E matrix.
  Note that drows_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddEValue(int& row, int& col, double& val);

  /*!
  \brief Add a value to the 'WS' map of this node

  The 'WS' map is later assembled to the WS matrix.
  Note that drows_ here is a vector.

  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddWSValue(int& row, int& col, double& val);
  /*!
  \brief Add a value to the MNode set of this node

  */
  void AddMNode (int node);

  /*!
  \brief Add a value to the map of Jump derivatives of this node

  Note that derivjump_ here is a vector.

  \param row : local dof row id to add to (rowise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddDerivJumpValue(int& row, const int& col, double val);
  
  void AddJumpValue(double val, int k);
  /*!
  \brief Add value to mechanical dissipation

  */
  void AddMechDissValue(double& val);

  /*!
  \brief Add a value to the 'A' map of this node

  \param row : local dof row id to add to (rowise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddAValue(int& row, int& col, double& val);
  
  /*!
  \brief Add a value to the ANode set of this node

  */
  void AddANode (int node);

  /*!
  \brief Add a value to the 'B' map of this node

  \param row : local dof row id to add to (rowise)
  \param val : value to be added
  \param col : global dof column id of the value added

  */
  void AddBValue(int& row, int& col, double& val);
  
  /*!
  \brief Add a value to the BNode set of this node

  */
  void AddBNode (int node);

  /*!
  \brief Set deltawear to value

   */
  void AddDeltaWearValue(double& val);
  
  /*!
  \brief Write nodal entries of D and M to Dold and Mold

  At the end of a time step the nodal entries (vector) of the Mortar
  Matrices D and M are stored to the old ones.

  */
  void StoreDMOld();

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

  /*!
  \brief Resets the data container of the node

  With this function, the container with contact specific quantities/information
  is deleted / reset to Teuchos::null pointer

  */
  virtual void ResetDataContainer();

  //@}

protected:
  
  // thermo contact related quantities
  double                              mechdiss_;     // mechanical dissipation
  std::vector<std::map<int,double> >  brows_;        // nodal rows of matrix A
  std::set<int>                       bnodes_;       // nodal set of according slave nodes

  Teuchos::RCP<CONTACT::FriNodeDataContainer>     fridata_;       // additional information of proc's friction nodes
  Teuchos::RCP<CONTACT::FriNodeDataContainerPlus> fridataplus_;   // additional information to proc's data container

  bool friplus_;                                      // bool for tsi || wear
}; // class FriNode

} // namespace CONTACT

//// << operator
//std::ostream& operator << (std::ostream& os, const CONTACT::FriNode& cnode);

#endif  // #ifndef FRICTION_NODE_H
