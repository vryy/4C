/*----------------------------------------------------------------------*/
/*!
\file PropagateCrack.H

\brief After each time step, check the structure field and propagate
this crack tip segment in the structure if necessary.

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef PROPAGATETIP_H_
#define PROPAGATETIP_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_crack.H"

namespace DRT {
namespace CRACK {

class PropagateTip {
public:

  /*!
  \brief constructor
   */
  PropagateTip( Teuchos::RCP<DRT::Discretization> discret,
                std::vector<int> tipnodes,
                std::set<int> cracknodes,
                int segment_id );

  /*!
  \brief destructor
   */
  ~PropagateTip(){};

  /*!
  \brief Set Cauchy stress values at Gauss points
   */
  void set_stress_values( Teuchos::RCP<std::vector<char> >& stressdata ){ gausspts_stress_ = stressdata; }

  /*!
  \brief Set Green-Lagrange strain values at Gauss points
   */
  void set_strain_values( Teuchos::RCP<std::vector<char> >& straindata ){ gausspts_strain_ = straindata; }

  /*!
  \brief Check whether crack tip propagates, if so find in the mesh which are the new crack tip nodes
   */
  void propagateThisTip( Teuchos::RCP<const Epetra_Vector>& displace );

  /*!
  \brief Get ALE displacement BC calculated at this crack segment at this step
   */
  std::map<int, std::vector<double> > getALEtipBC(){ return tip_bc_disp_; }

  /*!
  \brief Perform nodal releasing operations
  */
  void updateCrack( int& start_new_ele_id, int & start_new_node_id );

  /*!
  \brief Get old and new tip nodes after crack segment has propagated
   */
  std::map<int,int> getOldNewTipNodes(){ return oldnew_; }

  /*!
  \brief set all nodes that are falling on the crack surface
   */
  void setCrackNodes( std::set<int> cracknodes ){  cracknodes_ = cracknodes; }

  /*!
  \brief Do all operations related to completely splitting the body into two at this crack segment
   */
  void CheckCompleteSplit( int & start_new_node_id );

  /*!
  \brief Is the body completely split into two by this crack tip segment?
   */
  bool isSplit(){ return alreadySplit_; }

  /*!
  \brief Return true if the crack tip is propagated at this time step
   */
  bool isCrackPropagatedThisStep(){ return isProp_; }

  /*!
  \brief Add given node into the list of boundary nodes
   */
  void addBoundaryNode( int nodeid ){ boun_nodes_.insert( nodeid ); }

  /*!
  \brief Get all the nodes that are on the boundary (structure boundary + crack surfaces)
   */
  std::set<int> getAllBoundaryNodes(){ return boun_nodes_; }

  /*!
  \brief Replace the existing boundary nodes with the given nodes
   */
  void InitializeBoundaryNodes( std::set<int> bound ){ boun_nodes_ = bound; }

  /*!
  \brief Get segment tip nodes
   */
  std::vector<int> getSegmentTipNodes(){return tipnodes_;}

private:

  /*!
  \brief calculate stress intensity factors
   */
  void findStressIntensityFactor();

  /*!
  \brief Construct local coordinate system at crack tip
   */
  void findNormal( const DRT::Node * tipnode, const DRT::Node * surnode1, const DRT::Node * surnode2 );

  /*!
  \brief Build neighbor information to enable calculation of stress intensity factors using displacement correlation approach
   */
  void buildNeighborInfoTipNodes( DRT::Node * node1, DRT::Node * node2, DRT::Node * tipnode, int& attach1_id, int& attach2_id );

  /*!
  \brief Find a neighboring node that is on the crack surface
   */
  DRT::Node * findNeighboringCrackNode( DRT::Node * tipnode, bool second, const DRT::Node * first_node, int & attach_id );

  /*!
  \brief Find a node that is attached to this crack tip node
   */
  DRT::Node * findAttachedNode( DRT::Node * neigh, DRT::Node * tipnode );

  /*!
  \brief Get a surface which contains only one crack tip node, that is "tipid"
   */
  Teuchos::RCP<DRT::Element> getSurfaceThisPlane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces, int tipid );

  /*!
  \brief Returns true if mixed mode crack propagation criterion is met
   */
  bool DoCrackPropagate();

  /*!
  \brief Compute crack propagation angle
   */
  void decidePropagationAngle();

  /*!
  \brief Decide to which nodes in the mesh, the current crack tip should be propagated
   */
  std::vector<int> findNewCrackTip();

  /*!
  \brief Prepare data structure necessary to separate HEX into WEDGE elements
   */
  void SplitEleData( std::map<int,Teuchos::RCP<DRT::Element> >& zele, int add_id );

  /*!
  \brief Prepare data structure necessary to separate HEX into WEDGE elements
   */
  void SplitEleData( int ele_id, int node_id );

  /*!
  \brief Write location of crack tip
   */
  void WriteCrackTipLocation( const std::vector<int>& newTip );

  /*!
  \brief Write all points located on the crack surface
   */
  void WriteCrackSurfacePoints();

  /*!
  \brief Split appropriate HEX element into WEDGEs
   */
  void split_All_HEX_Elements( int &new_ele_id );

  /*!
  \brief Calculate two limiting angles those decides which elements get new node and which retains old tip nodes
   */
  std::vector<double> getLimitAngles( const std::vector<int>& newTip );

  /*!
  \brief Returns true if this element gets new node while performing nodal releasing approach
   */
  bool toReplaceNode( DRT::Element * ele, DRT::Node * tip, const std::vector<double>& lmtAngle );

  /*!
  \brief Returns true if the given node id exists on boundary
   */
  bool is_A_BoundaryNode( int & id );

  /*!
  \brief Data structure that stores details of how HEX element is split into WEDGEs
   */
  struct splitThisEle_
  {
    void print()
    {
      std::cout<<"split this element with id "<<element_id_<<" ";
      std::cout<<"Cut along diagonals with ids = "<<node_ids_ele_[0]<<" "<<node_ids_ele_[1]<<"\n";
    }

    int element_id_;
    std::vector<int> node_ids_ele_;
  };

  /*!
  \brief Compute vector J-integral around the crack tip
   */
  void Compute_J_integral();

  //! Collection of all HEX elements that must be split
  std::map<int,splitThisEle_> all_split_ele_;

  //! underlying structural discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! stress at Gauss points
  Teuchos::RCP<std::vector<char> > gausspts_stress_;

  //! Green-Lagrange strain at Gauss points
  Teuchos::RCP<std::vector<char> > gausspts_strain_;

  //! all the nodes that defines the crack tip
  std::vector<int> tipnodes_;

  //! store tip nodes at the previous time step. This is used to construct normal vector
  std::set<int> oldTipnodes_;

  //! all the nodes that are on the boundary of the domain
  std::set<int> boun_nodes_;

  //! all the nodes that are falling on the crack surface, including tip nodes
  std::set<int> cracknodes_;

  //! segment id of this crack tip node
  int segment_id_;

  //! Holds the displacement in column vector
  Teuchos::RCP<Epetra_Vector> disp_col_;

  //! Interface forces in column vector
  Teuchos::RCP<Epetra_Vector> iforce_col_;

  //! communicator
  const Epetra_Comm& comm_;

  //! processor rank
  int myrank_;

  //! critical stress intensity factors
  double critical_K_I_, critical_K_II_;

  //! critical J value -- used to check whether crack propagates or not
  double critical_J_;

  //! value of PI
  double PI_;

  //! value of zero angle -- set from crack_tolerance.H
  double min_angle_tol_;

  //! whether structure is completely split from this crack tip
  bool alreadySplit_;

  //! material id in the input file
  int material_id_;

  //! Young's modulus
  double young_;

  //! Poisson's ratio
  double poisson_;

  double kappa_;

  double fac_;

  //! unit normal vector at crack tip in current configuration
  LINALG::Matrix<3,1> normal_;

  //! unit tangent vector at crack tip in current configuration
  LINALG::Matrix<3,1> tangent_;

  //! unit normal vector at crack tip in reference configuration
  LINALG::Matrix<3,1> normal_ref_;

  //! unit tangent vector at crack tip in reference configuration
  LINALG::Matrix<3,1> tangent_ref_;

  //! data structure that stores two neighboring crack surface points for each tip node
  std::map<int, int> tip_phi_;
  std::map<int, int> tip_mphi_;

  //! stress intensity factors
  double K_I_, K_II_;

  //! crack propagation angle in radians
  double propAngle_;

  //! tip boundary condition displacements to be enforced
  std::map<int, std::vector<double> > tip_bc_disp_;

  //! file name to write crack tip location
  std::ostringstream filename_;

  //! file name to write J-integral
  std::ostringstream filename_jint_;

  //! All nodal ids that will be crack tip in next step
  std::vector<int> newTip_;

  //! map of old crack tip node ids and new node ids
  std::map<int, int> oldnew_;

  //! Is crack propagation conditions met?
  bool isProp_;

  //! crack propagation criterion
  INPAR::CRACK::propagationCriterion propcrit_;

  //! value of J-vector (stores only 2 components)
  std::vector<double> Jvector_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* PROPAGATETIP_H_ */
