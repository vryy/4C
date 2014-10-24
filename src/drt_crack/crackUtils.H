/*----------------------------------------------------------------------*/
/*!
\file crackUtils.H

\brief Utility functions for crack propagation problem

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef CRACKUTILS_H_
#define CRACKUTILS_H_

#include "../drt_lib/drt_discret.H"

namespace DRT {
namespace CRACK {
namespace UTILS{


  /*!
  \brief Update the given EpetraVector for the new map of discretization
   */
  void UpdateThisEpetraVectorCrack( const Teuchos::RCP<DRT::Discretization>& discret, Teuchos::RCP<Epetra_Vector>& vec,
                                  const std::map<int,int>& oldnewIds );

  /*!
  \brief Return the angle of the vector that points from a to b
   */
  double FindAngle( const double a[3], const double b[3] );

  /*!
  \brief convert the angle calculated in the range [-pi,pi] to [0,2pi]
   */
  void convertAngleTo_02PI_range( double &ang );

  /*!
  \brief convert the angle calculated in the range [0,2pi] to [-pi,pi]
   */
  void convertAngleTo_PI_mPI_range( double & ang );

  /*!
  \brief Returns true when the given element has given nodeid
   */
  bool ElementHasThisNodeId( const DRT::Element* ele, int nodeid );

  /*!
  \brief Returns true when the given element has given nodeid
   */
  bool ElementHasThisNodeId( const Teuchos::RCP<DRT::Element>& ele, int nodeid );

  /*!
  \brief Delete all Dirichlet conditions that has only one associated node
   */
  void deleteConditions( Teuchos::RCP<DRT::Discretization> discret );

  /*!
  \brief Add the given nodes and corresponding displacement as PointDirichlet conditions in discretization
   */
  void AddConditions( Teuchos::RCP<DRT::Discretization> discret, const std::map<int, std::vector<double> >& ale_bc_nodes );

  /*!
  \brief Add the nodal numbers given in the map to the given condition
   */
  void addNodesToConditions( DRT::Condition * cond, std::map<int,int> oldnew );

  /*!
  \brief Add the nodal numbers given in the vector to the given condition
   */
  void addNodesToConditions( DRT::Condition * cond, std::vector<int> add );

  /*!
  \brief Add the nodal numbers given in the set to the given condition
   */
  void addNodesToConditions( DRT::Condition * cond, std::set<int> add );

  double ComputeDiracDelta( const double & dx, const double & dy, const double & h );

  double HatFunction( const double & r );

  double SineHatFunction( const double &r, double n );

  double QuadraticHatFunction( const double &r, double n );

  /*!
  \brief Modify element connectivity of given discretization
   */
  void ModifyElementConnectivity( Teuchos::RCP<DRT::Discretization>& discret, std::map<int,int> delEle, std::map<int,int> oldnew );

  /*!
  \brief When a new node is introduced into the discretization, it is mandatory to copy the conditions
  existing in old nodes into its duplicated new nodes
   */
  void AddConditions( Teuchos::RCP<DRT::Discretization>& discret, std::map<int,int> oldnew );

  /*!
  \brief For a the given node on the given surface element, get indices of two nodes that are directly connected with this node
   */
  void getConnectedNodeIdIndex( const Teuchos::RCP<DRT::Element> & surfele, int & tip_id, int & conn1_index, int & conn2_index );

  /*!
  \brief Get nodal stress or strain values from their Gauss point values
   */
  void get_nodal_values_from_gauss_point_val( Teuchos::RCP<DRT::Discretization>& discret,
                                              Teuchos::RCP<std::vector<char> >& gausspts_val,
                                              Teuchos::RCP<Epetra_MultiVector>& nodal_val );

  /*!
  \brief Get displacement at given node
   */
  std::vector<double> getDisplacementNode( Teuchos::RCP<DRT::Discretization>& discret, const int node_gid, Teuchos::RCP<Epetra_Vector>& disp );

  /*!
  \brief Get displacement at given node
   */
  std::vector<double> getDisplacementNode( Teuchos::RCP<DRT::Discretization>& discret, const DRT::Node * node, Teuchos::RCP<Epetra_Vector>& disp );

} /* namespace CRACK */
} /* namespace CRACK */
} /* namespace DRT */
#endif
