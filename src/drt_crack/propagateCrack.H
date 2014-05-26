/*----------------------------------------------------------------------*/
/*!
\file PropagateCrack.H

\brief After each time step, check the structure field solution and
propagate crack in the structure if necessary conditions are met

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef PROPAGATECRACK_H_
#define PROPAGATECRACK_H_

#include "../drt_lib/drt_discret.H"

namespace DRT
{
namespace CRACK
{

class PropagateCrack
{
public:

  /*!
  \brief constructor
   */
  PropagateCrack( Teuchos::RCP<DRT::Discretization>& discret );

  ~PropagateCrack(){};

  /*!
  \brief Perform all the operations related to crack propagation
   */
  void propagateOperations( Teuchos::RCP<const Epetra_Vector>& displace,
                            Teuchos::RCP<std::vector<char> >& strdata );

  /*!
  \brief Return the map of old_node_id and the duplicated new_node_id
   */
  std::map<int,int> GetOldNewNodeIds(){ return oldnew_; }

	/*!
	\brief Return the crack tip nodes after crack has been propagated
	 */
  std::vector<int> GetCrackTipNodes(){ return tipnodes_; }

private:

  /*!
  \brief Initiate crack into the structure if not already present
   */
  void InitiateCrack( Teuchos::RCP<std::vector<char> >& strdata );

  /*!
  \brief Get attached nodes for each crack tip nodes, and construct neighboring node data structure which will be
  used to calculate stress intensity factors
   */
  void buildNeighborInfoTipNodes( DRT::Node * node1, DRT::Node * node2, DRT::Node * tipnode );

  /*!
  \brief Construct normal and tangential vectors at the crack tip
   */
  void findNormal( const DRT::Node * tipnode, const DRT::Node * surnode1, const DRT::Node * surnode2 );

  /*!
  \brief For each crack tip node, we get a neighboring node that falls on the crack surface
   */
  DRT::Node * findNeighboringCrackNode( const DRT::Node * tipnode, bool second, const DRT::Node * first_node );

  /*!
  \brief For each neighboring node, we find an attached node that satisfied certain criterion
   */
  DRT::Node * findAttachedNode( DRT::Node * neigh, DRT::Node * tipnode );

  /*!
  \brief Compute stress intensity factors at crack tip
   */
  void findStressIntensityFactor();

  /*!
  \brief Extract displacement at a particular node
   */
  std::vector<double> getDisplacementNode( const DRT::Node * node, Teuchos::RCP<Epetra_Vector>& disp );

  /*!
  \brief Returns true if the crack propagation criterion is satisfied
   */
  bool DoCrackPropagate();

  /*!
  \brief Calculate crack propagation angle from stress intensity factors
   */
  void decidePropagationAngle();

  /*!
  \brief Make all the modifications in the discretization related to crack propagation
   */
  void updateCrack( std::vector<int>& newTip );

  /*!
  \brief Get limiting angles for this crack configuration to decide to which elements to attach new nodes
   */
  std::vector<double> getLimitAngles( const std::vector<int>& newTip );

  /*!
  \brief Returns true if the criteron to replace the tipnode with a new duplicate node is satisfied
   */
  bool toReplaceNode( DRT::Element * ele, DRT::Node * tip, const std::vector<double>& lmtAngle );

  /*!
  \brief When crack propagated to the new location, we find the new crack tip nodes
   */
  std::vector<int> findNewCrackTip();
  std::vector<int> findNewCrackTip1();

  /*!
  \brief Get a particular surface that lies completely in the z-plane given by "coord"
   */
  Teuchos::RCP<DRT::Element> getSurfaceSameZplane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces,
                                                   const double * coord );

  /*!
  \brief Check whether the crack has completely split the structure into two, and do the operations necessary to stop the crack
  propagation after this point
   */
  void CheckCompleteSplit();

  /*!
  \brief When a new node is introduced into the discretization, it is mandatory to copy the conditions
  existing in old nodes into its duplicated new nodes
   */
  void AddConditions();

  void DeleteConditions();


  void printConditions(std::multimap<std::string,Teuchos::RCP<Condition> > allcondn);

  //! underlying structural discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! Holds the displacement in column vector
  Teuchos::RCP<Epetra_Vector> disp_col_;

  //! communicator
  const Epetra_Comm& comm_;

  //! all the nodes that defines the crack tip
  std::vector<int> tipnodes_;

  //! all the nodes that are on the boundary of the domain
  std::vector<int> boun_nodes_;

  //! unit normal vector at crack tip
  LINALG::Matrix<3,1> normal_;

  //! unit tangent vector at crack tip
  LINALG::Matrix<3,1> tangent_;

  //! all the nodes that are falling on the crack surface, including tip nodes
  std::set<int> cracknodes_;

  //! store tip nodes at the previous time step. This is used to construct normal vector
  std::set<int> oldTipnodes_;

  //! map of old crack tip node ids and new node ids
  std::map<int, int> oldnew_;

  //! stress intensity factors
  double K_I_, K_II_;

  //! critical stress intensity factors
  double critical_K_I_, critical_K_II_;

  //! Young's modulus
  double young_;

  //! Poisson's ratio
  double poisson_;

  double kappa_;

  double fac_;

  //! crack propagation angle in radians
  double propAngle_;

  // This is in RADIANS
  //! If the angle between an edge and the exact crack propagation angle is less than this value,
  //! then crack is assumed to propagate through this edge
  double min_angle_tol_;

  //! Whether the crack tip nodes should be moved in the current configuration
  bool moveNodes_;

  //! tip boundary condition displacements to be enforced
  std::map<int, std::vector<double> > tip_bc_disp_;

  //! Whether Dirichlet BC set at the crack tip nodes should be cleared
  bool clearCondns_;

  int myrank_;

  //! data structure that stores two neighboring crack surface points for each tip node
  std::map<int, int> tip_phi_;
  std::map<int, int> tip_mphi_;

  //! Whether crack has propagated completely to such a extent that the structure is completely split into two
  bool strIsSplit_;

  //! The Id of the first node that will added into discretization after crack propagation. This will be incremented
  //! each time we add a new element
  // This number is set from the input file. It should be more than the total no of nodes in all the discretization of the problem
  // included
  int startNewId_;

  //! value of PI
  double PI_;

  //! Is crack already initiated
  bool crackInitiated_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* PROPAGATECRACK_H_ */
