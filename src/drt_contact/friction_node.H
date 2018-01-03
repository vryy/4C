/*!----------------------------------------------------------------------
\file friction_node.H

\brief A class for a frictional contact node
\level 2
\maintainer Alexander Seitz


*-----------------------------------------------------------------------*/
#ifndef FRICTION_NODE_H
#define FRICTION_NODE_H

#include "contact_node.H"

/*!
 \brief CONTACT: namespace of the BACI contact and meshtying algorithms

 */
namespace CONTACT
{

class FriNodeType: public DRT::ParObjectType
{
public:

  std::string Name() const
  {
    return "FriNodeType";
  }

  static FriNodeType & Instance()
  {
    return instance_;
  }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

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
  virtual ~FriNodeDataContainer()
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
   \brief Return jump per time step (only for slave side!) (length 3)
   */
  virtual inline double* jump()
  {
    return jump_;
  }

  /*!
   \brief Return jump per time step (only for slave side!) (max length 2)
   */
  virtual inline double* jump_var()
  {
    return jumpvar_;
  }

  /*!
   \brief Return friction status of this node (slip=true)
   */
  virtual inline bool& Slip()
  {
    return slip_;
  }

  /*!
   \brief Return the old 'D' map (vector) of this node (last converged state)
   */
  virtual inline GEN::pairedvector<int, double>& GetDOld()
  {
    return drowsold_;
  }

  /*!
   \brief Return the old 'M' map (vector) of this node (last converged state)
   */
  virtual inline std::map<int, double>& GetMOld()
  {
    return mrowsold_;
  }

  /*!
   \brief Return the old 'D' map (vector) of this node (last converged state)
   */
  virtual inline GEN::pairedvector<int, double>& GetDOldLTL()
  {
    return drowsoldLTL_;
  }

  /*!
   \brief Return the old 'M' map (vector) of this node (last converged state)
   */
  virtual inline std::map<int, double>& GetMOldLTL()
  {
    return mrowsoldLTL_;
  }

  /*!
   \brief Return the map with according slave nodes
   */
  virtual inline std::set<int>& GetSNodes()
  {
    return snodes_;
  }

  /*!
   \brief Return the map with according master nodes
   */
  virtual inline std::set<int>& GetMNodes()
  {
    return mnodes_;
  }

  /*!
   \brief Return the old map with according master nodes
   */
  virtual inline std::set<int>& GetMNodesOld()
  {
    return mnodesold_;
  }

  /*!
   \brief Return the 'DerivJump' map (vector) of this node

   These maps contain the directional derivatives of the node's
   relative movement (jump).
   A vector is used because the jump itself is a vector (2 or 3 components).

   */
  virtual inline std::vector<std::map<int, double> >& GetDerivJump()
  {
    return derivjump_;
  }

  /*!
   \brief Return the 'DerivVarJump' map (vector) of this node

   These maps contain the directional derivatives of the node's
   relative -- object variant-- movement (jump). This Jump is already
   multiplied and linearized with the tangent vectors txi and teta.
   Therefore, the vector has 1/2 components instead of 2/3

   */
  virtual inline std::vector<std::map<int, double> >& GetDerivVarJump()
  {
    return derivvarjump_;
  }
  /*!
   \brief Return current penalty traction (length 3)
   */
  virtual inline double* traction()
  {
    return traction_;
  }
  /*!
   \brief Return current penalty traction (length 3)
   */
  virtual inline double* tractionLTL()
  {
    return tractionLTL_;
  }
  /*!
   \brief Return contact status of last converged state n (active=true)
   */
  virtual inline bool& SlipOld()
  {
    return slipold_;
  }

  /*!
   \brief Return old penalty traction (length 3)
   */
  virtual inline double* tractionold()
  {
    return tractionold_;
  }
  /*!
   \brief Return old penalty traction (length 3)
   */
  virtual inline double* tractionoldLTL()
  {
    return tractionoldLTL_;
  }

  //@}

protected:

  // don't want = operator and cctor
  FriNodeDataContainer operator =(const FriNodeDataContainer& old);
  FriNodeDataContainer(const FriNodeDataContainer& old);

  double jump_[3];    // jump per time step
  double jumpvar_[2]; // variant jump per time step
  bool slip_;         // true if cnode slips (contact force reaches friction bound)
  bool slipold_;      // true if cnode slips (contact force reaches friction bound) - last converged timestep

  GEN::pairedvector<int, double> drowsold_; // nodal rows of old D matrix
  std::map<int, double> mrowsold_;          // nodal rows of old M matrix
  GEN::pairedvector<int, double> drowsoldLTL_; // nodal rows of old D matrix
  std::map<int, double> mrowsoldLTL_;          // nodal rows of old M matrix
  std::set<int> snodes_;    // nodal set of according slave nodes
  std::set<int> mnodes_;    // nodal set of according master nodes
  std::set<int> mnodesold_; // nodal set of old according masternodes

  std::vector<std::map<int, double> > derivjump_;    // directional derivative of nodal weighted jump vector
  std::vector<std::map<int, double> > derivvarjump_; // directional derivative of obj.-variant. nodal weighted jump vector

  // penalty related quantities
  double traction_[3];    // traction vector of current time step
  double tractionold_[3]; // traction vector of previous time step

  // penalty related quantities
  double tractionLTL_[3];    // traction vector of current time step
  double tractionoldLTL_[3]; // traction vector of previous time step
};
// class FriNodeDataContainer

/*!
 \brief An additional container for wear and tsi with contact specific data

 This class contains additional information to the data container of the
 frictional node. These additional data are needed only for contact
 problems with wear and thermo-structure-interaction problems with contact.

 \author gitterle (gitterle@lnm.mw.tum.de)
 */
class FriNodeWearDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
   \brief Standard Constructor

   */
  FriNodeWearDataContainer();

  /*!
   \brief Destructor

   */
  virtual ~FriNodeWearDataContainer()
  {
    return;
  }

  //@}

  //! @name Access methods

  virtual void Pack(DRT::PackBuffer& data) const;
  virtual void Unpack(std::vector<char>::size_type& position,
      const std::vector<char>& data);

  /*!
   \brief Return the weighted wear per node (length 1)
   */
  virtual double& WeightedWear()
  {
    return weightedwear_;
  }

  /*!
   \brief Return the delta of the weighted wear per node without wear coefficient (length 1)
   */
  virtual double& DeltaWeightedWear()
  {
    return deltaweightedwear_;
  }

  /*!
   \brief Return the 'T' map (vector) of this node
   */
  virtual std::vector<std::map<int, double> >& GetT()
  {
    return trows_;
  }

  /*!
   \brief Return the 'E' map (vector) of this node
   */
  virtual std::vector<std::map<int, double> >& GetE()
  {
    return erows_;
  }

  /*!
   \brief Deriv. w.r.t. displ. of E and T matrix entries for
   this node. This matrices are required for the
   discrete wear condition.
   */
  virtual std::map<int, std::map<int, double> >& GetDerivTw()
  {
    return derivt_;
  }
  virtual std::map<int, std::map<int, double> >& GetDerivE()
  {
    return derive_;
  }

  /*!
   \brief Return the 'D2' map (vector) of this node
   */
  virtual std::vector<std::map<int, double> >& GetD2()
  {
    return d2rows_;
  }

  /*!
   \brief Return current discrete wear in step n+1 (only for slave side!) (length 1)
   */
  virtual double* wcurr()
  {
    return wcurr_;
  }

  /*!
   \brief Return accumulated wear for different time scales (only for slave side!) (length 1)
   */
  virtual double* waccu()
  {
    return waccu_;
  }

  /*!
   \brief Return old discrete wear in time step i (only for slave side!) (length 1)
   */
  virtual double* wold()
  {
    return wold_;
  }
  //@}

protected:

  // don't want = operator and cctor
  FriNodeWearDataContainer operator =(const FriNodeWearDataContainer& old);
  FriNodeWearDataContainer(const FriNodeWearDataContainer& old);

  // wear related quantities
  double weightedwear_;      // weighted wear
  double deltaweightedwear_; // delta weighted wear
  double wcurr_[1];          // current pv wear value (n+1)
  double wold_[1];           // old pv wear value (i) - for paritioned solution scheme
  double waccu_[1];          // accumulated pv wear value (i) - for different time scales

  std::vector<std::map<int, double> > trows_; // nodal rows of T matrix
  std::vector<std::map<int, double> > erows_; // nodal rows of E matrix

  std::vector<std::map<int, double> > d2rows_; // nodal rows of master sided D matrix

  std::map<int, std::map<int, double> > derivt_; // directional derivative of nodal Tw-matrix values
  std::map<int, std::map<int, double> > derive_; // directional derivative of nodal E-matrix values

};
// class FriNodeWearDataContainer

/*!
 \brief A class for a frictional contact node derived from CONTACT::CoNode

 This class represents a finite element node capable of frictional contact.

 \author gitterle (gitterle@lnm.mw.tum.de)
 */
class FriNode: public CONTACT::CoNode
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
  FriNode(int id, const double* coords, const int owner, const int numdof,
      const std::vector<int>& dofs, const bool isslave, const bool initactive,
      const bool friplus);

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
  virtual ~FriNode()
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
    return FriNodeType::Instance().UniqueParObjectId();
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

  //! @name access methods

  /*!
   \brief Print this cnode
   */
  virtual void Print(std::ostream& os) const;

  /*!
   \brief Return of data container of this node

   This method returns the data container of this node where additional
   contact specific quantities/information are stored.

   */
  inline CONTACT::FriNodeDataContainer& FriData()
  {
    return *fridata_;
  }

  /*!
   \brief Return of additional data container of this node

   This method returns the additional data container of this node where additional
   wear specific quantities/information are stored.

   */
  inline CONTACT::FriNodeWearDataContainer& WearData()
  {
    return *weardata_;
  }

  /*!
   \brief calculate the apparent coefficient of friction

   in a TSI problem, this is dependent on the slave and master
   temperatures. Since the master temperature is displacement
   dependent (due to the involved projections), it requires
   a consistent linearization (see function below).
   */
  double FrCoeff(const double& frcoeff_in);

  /*!
   \brief calculate the derivative of apparent coefficient of friction

   in a TSI problem, this is dependent on the slave and master
   temperatures. Since the master temperature is displacement
   dependent (due to the involved projections), it requires
   a consistent linearization (see function below).
   */
  void derivFrCoeffTemp(
      const double& frcoeff_in,
      std::map<int,double>& derivT,
      std::map<int,double>& derivDisp);

  //@}

  //! @name Evaluation methods

  /*!
   \brief Add a value to the SNode set of this node

   */
  void AddSNode(int node);

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
  void AddMNode(int node);

  /*!
   \brief Add a value to the map of Jump derivatives of this node

   Note that derivjump_ here is a vector.

   \param row : local dof row id to add to (rowise)
   \param val : value to be added
   \param col : global dof column id of the value added

   */
  void AddDerivJumpValue(int& row, const int& col, double val);

  void AddJumpValue(double val, int k);

  void AddD2Value(int& row, int& col, double& val);

  /*!
   \brief Set deltawear to value

   */
  void AddDeltaWeightedWearValue(double& val);

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

  Teuchos::RCP<CONTACT::FriNodeDataContainer> fridata_;         // additional information of proc's friction nodes
  Teuchos::RCP<CONTACT::FriNodeWearDataContainer> weardata_;    // additional information to proc's data container

  bool wear_; // bool for wear
};
// class FriNode

}// namespace CONTACT

//// << operator
//std::ostream& operator << (std::ostream& os, const CONTACT::FriNode& cnode);

#endif  // #ifndef FRICTION_NODE_H
