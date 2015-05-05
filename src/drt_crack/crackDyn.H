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

#ifndef CRACKDYN_H_
#define CRACKDYN_H_

#include "../drt_lib/drt_discret.H"
#include "propagateTip.H"

//#include "../drt_io/io.H"

namespace ALE
{
class Ale;
}

namespace IO
{
class DiscretizationWriter;
class DiscretizationReader;
}

namespace DRT
{
namespace CRACK
{
class CrackDyn
{
public:

  /*!
  \brief constructor
   */
  CrackDyn( Teuchos::RCP<DRT::Discretization>& discret );

  ~CrackDyn(){};

  /*!
  \brief Perform all the operations related to crack propagation
   */
  void propagateOperations( Teuchos::RCP<const Epetra_Vector>& displace,
                            Teuchos::RCP<std::vector<char> >& stressdata,
                            Teuchos::RCP<std::vector<char> >& straindata );

  /*!
  \brief Return the map of old_node_id and the duplicated new_node_id
   */
  std::map<int,int> GetOldNewNodeIds(){ return oldnew_; }

   /*!
   \brief Return the crack tip nodes after crack has been propagated
   */
  std::vector<int> GetCrackTipNodes(){ return tipnodes_; }

  /*!
  \brief Returns true if the Dirichlet conditions are cleared in this time step
   */
  bool DidIClearConditionsNow(){ return justClearedCondns_; }

  /*!
  \brief Write all data related to restart
   */
  void WriteRestartCrack( Teuchos::RCP<IO::DiscretizationWriter> outputWriter );

  /*!
  \brief Read restart data for crack
   */
  void ReadRestartCrack( IO::DiscretizationReader& reader );

private:

  /*!
  \brief Initiate crack into the structure if not already present
   */
  void InitiateCrack( Teuchos::RCP<std::vector<char> >& stressdata,
                      Teuchos::RCP<std::vector<char> >& straindata );

  /*!
  \brief Get a particular surface that lies completely in the z-plane given by "coord"
   */
  Teuchos::RCP<DRT::Element> getSurfaceSameZplane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces,
                                                   const double * coord );

  /*!
  \brief Perform ALE step of moving the crack tip to the desired location
   */
  void Perform_ALE_Step( Teuchos::RCP<const Epetra_Vector>& displace, std::map<int, std::vector<double> >& tipbc_disp );

  /*!
  \brief Analyze and clean discretization after adding new nodes
   */
  void analyzeAndCleanDiscretization();

  /*!
  \brief Delete all conditions that has only one node. These are the conditions added during last step of crack propagation
   */
  void DeleteConditions();

  /*!
  \brief Print all conditions of the discretization
   */
  void printConditions(std::multimap<std::string,Teuchos::RCP<Condition> > allcondn);

  /*!
  \brief Write nodes and elements in gmsh format for debugging
   */
  void WriteNodesDebug();

  /*!
  \brief Update the nodes that carries ALE boundary condition
   */
  void updateALE_BC_nodes();

  /*!
  \brief Update boundary nodes of discretization
   */
  void update_boundary_nodes();

  /*!
  \brief Update crack tip nodes
   */
  void update_crack_tip_nodes();

  //! underlying structural discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! communicator
  const Epetra_Comm& comm_;

  //! all the nodes that defines the crack tip
  std::vector<int> tipnodes_;

  //! all the nodes that are falling on the crack surface, including tip nodes
  std::set<int> cracknodes_;

  //! map of old crack tip node ids and new node ids
  std::map<int, int> oldnew_;

  //! new Point Dirichlet ALE boundary nodes and the corresponding displacement values
  std::set<int> new_ale_bc_nodes_;

  //! Whether Dirichlet BC set at the crack tip nodes should be cleared
  bool clearCondns_;

  //! rank of processor
  int myrank_;

  //! Whether crack has propagated completely to such a extent that the structure is completely split into two
  bool strIsSplit_;

  //! The Id of the first node that will added into discretization after crack propagation. This will be incremented
  //! each time we add a new element
  // This number is set from the input file. It should be more than the total no of nodes in all the discretization of the problem
  // included
  int startNewNodeId_;

  //! The Id of the first element that will added into discretization after crack propagation. This will be incremented
  //! each time we add a new element
  // This number is set from the input file. It should be more than the total no of nodes in all the discretization of the problem
  // included
  int startNewEleId_;

  //! Is crack already initiated
  bool crackInitiated_;

  //! cleared conditions in this time step
  bool justClearedCondns_;

  //! Nodes that carries ALE line Dirichlet conditions
  std::set<int> ALE_line_nodes_;

  //! Returns true for the first time crack propagation
  int step_;

  //! All crack tip segments
  std::vector<Teuchos::RCP<PropagateTip> > tip_segments_;

  //! inner layer of nodes used for crack initiation
  std::vector<int> inner_layer_nodes_;

  //! outer layer of nodes used for crack initiation
  std::vector<int> outer_layer_nodes_;

  //! whether crack just initiated at this step
  bool crackInitiatedThisStep_;

  //! bool to store the information whether it is a restart
  bool isRestart_;

};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* CRACKDYN_H_ */
