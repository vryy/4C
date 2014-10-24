/*----------------------------------------------------------------------*/
/*!
\file initiateCrack.H

\brief Crack initiation procedure in a homogeneous flawless
structure.

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef INITIATECRACK_H_
#define INITIATECRACK_H_

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"

namespace DRT {
namespace CRACK {

class InitiateCrack {
public:

  /*!
  \brief Constructor
   */
  InitiateCrack( Teuchos::RCP<DRT::Discretization> dis,
                 Teuchos::RCP<std::vector<char> >& strdata,
                 std::string criterion )
  :discret_( dis ),
   gausspts_str_( strdata ),
   criterion_( criterion ),
   comm_( discret_->Comm() ),
   myrank_( comm_.MyPID() )
  {
  }

  /*!
  \brief destructor
   */
  ~InitiateCrack(){};

  /*!
  \brief Perform all operations related to crack initiation
   */
  bool initiateOperations( int & start_new_node_id );

  /*!
  \brief Set the value of stress or strain at which failure occur
   */
  void set_failure_value( double val ){ failure_value_ = val; }

  /*!
  \brief Set inner layer of nodes to determine split nodes
   */
  void set_inner_layer_nodes( std::vector<int>& nodes ){ inner_nodes_ = nodes; }

  /*!
  \brief Set outer layer of nodes to determine split nodes
   */
  void set_outer_layer_nodes( std::vector<int>& nodes ){ outer_nodes_ = nodes; }

  /*!
  \brief Get old and new crack nodes
   */
  std::map<int,int> getOldNewNodes(){ return oldnew_; }

  /*!
  \brief Get all nodes that are falling on crack surface
   */
  std::set<int> getCrackSurfaceNodes(){ return crackSurfNodes_; }

  /*!
  \brief Get crack tip segment nodes that are at positive x-distance from split nodes
   */
  std::vector<int> getPositiveCrackTipNodes(){ return posCrackTip_; }

  /*!
  \brief Get crack tip segment nodes that are at negative x-distance from split nodes
   */
  std::vector<int> getNegativeCrackTipNodes(){ return negCrackTip_; }

private:

  /*!
  \brief Determine the  nodal at which von Mises stress is maximum
   */
  void Determine_node_max_von_mises_stress( int & max_von_node, double & max_von_stress );

  /*!
  \brief Determine split nodes
   */
  void Determine_split_nodes( int& max_von_node );

  /*!
  \brief Determine all nodes that will be assigned as crack tip nodes
   */
  void Determine_crack_tip_nodes();

  /*!
  \brief For each split nodes, determine a connected node that is in max positive and max negative x-distance
   */
  void get_positive_negative_xnode( int & splnode, std::vector<int>& nodes, int& pos_node, int& neg_node );

  /*!
  \brief Get a surface in this element that has all the given nodes
   */
  Teuchos::RCP<DRT::Element> get_surface_all_in_nodes( DRT::Element * ele, std::vector<int>& nodes );

  /*!
  \brief Perform nodal releasing operation at split nodes
   */
  void nodal_release( int & start_new_node_id );

  /*!
  \brief Decide projection plane to which given 3D geometry is projected
   */
  int getProjectionPlane();

  /*!
  \brief Add appropriate nodes to crack surfaces
   */
  void build_tip_surface_nodes();

  /*!
  \brief Get limiting angles that will be used to maintain appropriate connectivity while performing nodal release
   */
  std::vector<double> FindLimitingAngles( const int &ind1, const int & ind2 );

  //! underlying structural discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! stress at Gauss points
  Teuchos::RCP<std::vector<char> >& gausspts_str_;

  //! stores stress or strain based criterion
  std::string criterion_;

  //! communicator
  const Epetra_Comm& comm_;

  //! processor id
  int myrank_;

  //! value of failure stress or strain
  double failure_value_;

  //! inner layer of nodes
  std::vector<int> inner_nodes_;

  //! outer layer of nodes
  std::vector<int> outer_nodes_;

  //! stores split nodes and its location whether it is in inner or outer layer nodes
  std::map<int,std::string> split_nodes_;

  //! stores split nodes, and corresponding positive and negative x-nodes
  std::map<int, std::pair<int,int> > split_pos_neg_;

  //! stores split nodes and their duplicated nodes
  std::map<int,int> oldnew_;

  //! crack tip segment nodes at positive distance
  std::vector<int> posCrackTip_;

  //! crack tip segment nodes at negative distance
  std::vector<int> negCrackTip_;

  //! all nodes that are on the crack surface after initiation
  std::set<int> crackSurfNodes_;

  //! When initiating crack, distance along this direction is considered
  int dist_index_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* INITIATECRACK_H_ */
