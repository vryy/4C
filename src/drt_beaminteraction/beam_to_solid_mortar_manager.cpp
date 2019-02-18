/*!
\file beam_to_solid_mortar_manager.cpp

\brief Manage the creation of additional DOFs for mortar couplings between beams and solids.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_mortar_manager.H"

#include "beam_contact_pair.H"
#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "beaminteraction_calc_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidMortarManager::BeamToSolidMortarManager(
    const INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions& mortar_shape_function)
    : is_setup_(false),
      is_local_maps_build_(false),
      index_base_(0),
      lambda_dof_rowmap_(Teuchos::null),
      node_gid_to_lambda_gid_(Teuchos::null),
      element_gid_to_lambda_gid_(Teuchos::null)
{
  // Get the number of Lagrange multiplier DOF on a beam node and on a beam element.
  switch (mortar_shape_function)
  {
    case INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::line2:
    {
      n_lambda_node_ = 1 * 3;
      n_lambda_element_ = 0 * 3;
      break;
    }
    case INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::line3:
    {
      n_lambda_node_ = 1 * 3;
      n_lambda_element_ = 1 * 3;
      break;
    }
    default:
      dserror("Mortar shape function not implemented!");
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::Setup(
    const Teuchos::RCP<DRT::Discretization> discret)
{
  // Get the global ids of all beam centerline nodes on this rank.
  std::vector<int> my_nodes_gid;
  for (int i_node = 0; i_node < discret->NodeRowMap()->NumMyElements(); i_node++)
  {
    DRT::Node const& node = *(discret->lRowNode(i_node));
    if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node)) my_nodes_gid.push_back(node.Id());
  }

  // Get the global ids of all beam elements on this rank.
  std::vector<int> my_elements_gid;
  for (int i_element = 0; i_element < discret->ElementRowMap()->NumMyElements(); i_element++)
  {
    DRT::Element const& element = *(discret->lRowElement(i_element));
    if (BEAMINTERACTION::UTILS::IsBeamElement(element)) my_elements_gid.push_back(element.Id());
  }

  // Calculate the local number of centerline nodes, beam elements and Lagrange multiplier DOF.
  const unsigned int n_nodes = my_nodes_gid.size();
  const unsigned int n_element = my_elements_gid.size();
  const unsigned int n_lambda_dof = n_nodes * n_lambda_node_ + n_element * n_lambda_element_;

  // Rowmap for the additional DOFs used by the mortar contact discretization.
  lambda_dof_rowmap_ = Teuchos::rcp(new Epetra_Map(-1, n_lambda_dof, index_base_, discret->Comm()));

  // We need to be able to get the global ids for a Lagrange multiplier DOF from the global id of a
  // node or element. To do so, we 'abuse' the Epetra_MultiVector as map between the global node /
  // element ids and the global Lagrange multiplier DOF ids.
  Teuchos::RCP<Epetra_Map> node_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_nodes, &my_nodes_gid[0], 0, discret->Comm()));
  Teuchos::RCP<Epetra_Map> element_gid_rowmap =
      Teuchos::rcp(new Epetra_Map(-1, n_element, &my_elements_gid[0], 0, discret->Comm()));

  // Map from global node / element ids to global lagrange multiplier ids. Only create the
  // multivector if it hase one or more columns.
  node_gid_to_lambda_gid_ = Teuchos::null;
  element_gid_to_lambda_gid_ = Teuchos::null;
  if (n_lambda_node_ > 0)
    node_gid_to_lambda_gid_ =
        Teuchos::rcp(new Epetra_MultiVector(*node_gid_rowmap, n_lambda_node_, true));
  if (n_lambda_element_ > 0)
    element_gid_to_lambda_gid_ =
        Teuchos::rcp(new Epetra_MultiVector(*element_gid_rowmap, n_lambda_element_, true));

  // Fill in the entries in the node / element global id to Lagrange multiplier global id vector.
  int error_code = 0;
  int lagrange_gid = -1;
  if (node_gid_to_lambda_gid_ != Teuchos::null)
  {
    for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
      for (unsigned int i_lambda = 0; i_lambda < n_lambda_node_; i_lambda++)
      {
        // Get the global Lagrange multiplier id for this node.
        lagrange_gid = lambda_dof_rowmap_->GID(i_node * n_lambda_node_ + i_lambda);
        if (lagrange_gid < index_base_) dserror("Local id not found on this processor!");

        // Set the global Lagrange multiplier id for this node.
        error_code = node_gid_to_lambda_gid_->ReplaceMyValue(i_node, i_lambda, lagrange_gid);
        if (error_code != 0) dserror("Got error code %d!", error_code);
      }
  }
  if (element_gid_to_lambda_gid_ != Teuchos::null)
  {
    for (unsigned int i_element = 0; i_element < n_element; i_element++)
      for (unsigned int i_lambda = 0; i_lambda < n_lambda_element_; i_lambda++)
      {
        // Get the global Lagrange multiplier id for this node.
        lagrange_gid = lambda_dof_rowmap_->GID(
            n_nodes * n_lambda_node_ + i_element * n_lambda_element_ + i_lambda);
        if (lagrange_gid < index_base_) dserror("Local id not found on this processor!");

        // Set the global Lagrange multiplier id for this node.
        error_code = element_gid_to_lambda_gid_->ReplaceMyValue(i_element, i_lambda, lagrange_gid);
        if (error_code != 0) dserror("Got error code %d!", error_code);
      }
  }

  // Set flag for successfull setup.
  is_setup_ = true;
  is_local_maps_build_ = false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::SetLocalMaps(
    const Teuchos::RCP<DRT::Discretization> discret,
    const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> contact_pairs)
{
  CheckSetup();

  // At this point the global multi vectors are filled up completely. To get the map for global
  // node element ids to the global lambda ids we need to be able to extract more than the local
  // values on this processor. Therefore we need a new map that contains all rows we want to
  // access in the global multi vector.
  std::vector<int> node_gid_needed;
  std::vector<int> element_gid_needed;

  // Loop over the pairs and get the global node and element indices needed on this rank.
  for (unsigned int i_pair = 0; i_pair < contact_pairs.size(); i_pair++)
  {
    const Teuchos::RCP<BEAMINTERACTION::BeamContactPair>& pair = contact_pairs[i_pair];

    // The first (beam) element should always be on the same processor as the pair.
    if (pair->Element1()->Owner() != discret->Comm().MyPID())
      dserror(
          "The current implementation need the first element of a beam contact pair to be on the "
          "same processor as the pair!");

    // Get the global id of the nodes / elements that the pairs on this rank need.
    if (n_lambda_node_ > 0)
      // There are nodal lambda DOFs, add the gid for the nodes in this element to the vector.
      // The first two nodes are the centerline nodes.
      for (unsigned int i_node = 0; i_node < 2; i_node++)
        node_gid_needed.push_back(pair->Element1()->Nodes()[i_node]->Id());

    if (n_lambda_element_ > 0)
      // There are element lambda DOFs, add the gid for this element to the vector.
      element_gid_needed.push_back(pair->Element1()->Id());
  }

  // Make the entries in the vectors unique.
  std::vector<int>::iterator it;
  std::sort(node_gid_needed.begin(), node_gid_needed.end());
  it = std::unique(node_gid_needed.begin(), node_gid_needed.end());
  node_gid_needed.resize(std::distance(node_gid_needed.begin(), it));
  std::sort(element_gid_needed.begin(), element_gid_needed.end());
  it = std::unique(element_gid_needed.begin(), element_gid_needed.end());
  element_gid_needed.resize(std::distance(element_gid_needed.begin(), it));

  // Create the maps for the extraction of the values.
  Teuchos::RCP<Epetra_Map> node_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, node_gid_needed.size(), &node_gid_needed[0], 0, discret->Comm()));
  Teuchos::RCP<Epetra_Map> element_gid_needed_rowmap = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, element_gid_needed.size(), &element_gid_needed[0], 0, discret->Comm()));

  // Create the Multivectors that will be filled with all values needed on this rank.
  Teuchos::RCP<Epetra_MultiVector> node_gid_to_lambda_gid_copy = Teuchos::rcp<Epetra_MultiVector>(
      new Epetra_MultiVector(*node_gid_needed_rowmap, n_lambda_node_, true));
  Teuchos::RCP<Epetra_MultiVector> element_gid_to_lambda_gid_copy =
      Teuchos::rcp<Epetra_MultiVector>(
          new Epetra_MultiVector(*element_gid_needed_rowmap, n_lambda_element_, true));

  // Export values from the global multi vector to the ones needed on this rank.
  if (node_gid_to_lambda_gid_ != Teuchos::null)
    LINALG::Export(*node_gid_to_lambda_gid_, *node_gid_to_lambda_gid_copy);
  if (element_gid_to_lambda_gid_ != Teuchos::null)
    LINALG::Export(*element_gid_to_lambda_gid_, *element_gid_to_lambda_gid_copy);

  // Fill in the local maps.
  node_gid_to_lambda_gid_map_.clear();
  element_gid_to_lambda_gid_map_.clear();
  std::vector<int> temp_node(n_lambda_node_);
  for (int i_node = 0; i_node < node_gid_needed_rowmap->NumMyElements(); i_node++)
  {
    for (unsigned int i_temp = 0; i_temp < n_lambda_node_; i_temp++)
      temp_node[i_temp] = (int)((*(*node_gid_to_lambda_gid_copy)(i_temp))[i_node]);
    node_gid_to_lambda_gid_map_[node_gid_needed_rowmap->GID(i_node)] = temp_node;
  }
  std::vector<int> temp_elements(n_lambda_element_);
  for (int i_element = 0; i_element < element_gid_needed_rowmap->NumMyElements(); i_element++)
  {
    for (unsigned int i_temp = 0; i_temp < n_lambda_element_; i_temp++)
      temp_elements[i_temp] = (int)((*(*element_gid_to_lambda_gid_copy)(i_temp))[i_element]);
    element_gid_to_lambda_gid_map_[element_gid_needed_rowmap->GID(i_element)] = temp_elements;
  }

  // Set flags for local maps.
  is_local_maps_build_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::LocationVector(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactPair> contact_pair,
    std::vector<int>& lambda_row) const
{
  CheckSetup();
  CheckLocalMaps();

  // Clear the output vectors.
  lambda_row.clear();

  // Get the global DOFs ids of the nodal Lagrange multipliers.
  if (n_lambda_node_ > 0)
  {
    for (int i_node = 0; i_node < contact_pair->Element1()->NumNode(); i_node++)
    {
      const DRT::Node& node = *(contact_pair->Element1()->Nodes()[i_node]);
      if (BEAMINTERACTION::UTILS::IsBeamCenterlineNode(node))
      {
        // Get the global id of the node.
        int node_id = node.Id();

        // Check if the id is in the map. If it is, add it to the output vector.
        auto search_key_in_map = node_gid_to_lambda_gid_map_.find(node_id);
        if (search_key_in_map == node_gid_to_lambda_gid_map_.end())
          dserror("Global node id %d not in map!", node_id);
        for (auto const& lambda_gid : search_key_in_map->second) lambda_row.push_back(lambda_gid);
      }
    }
  }

  // Get the global DOFs ids of the element Lagrange multipliers.
  if (n_lambda_element_ > 0)
  {
    if (BEAMINTERACTION::UTILS::IsBeamElement(*contact_pair->Element1()))
    {
      // Get the global id of the element.
      int element_id = contact_pair->Element1()->Id();

      // Check if the id is in the map. If it is, add it to the output vector.
      auto search_key_in_map = element_gid_to_lambda_gid_map_.find(element_id);
      if (search_key_in_map == element_gid_to_lambda_gid_map_.end())
        dserror("Global element id %d not in map!", element_id);
      for (auto const& lambda_gid : search_key_in_map->second) lambda_row.push_back(lambda_gid);
    }
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidMortarManager::EvaluateGlobalDM(
    const std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>& contact_pairs)
{
  // Local mortar matrices.
  LINALG::SerialDenseMatrix mortar_D;
  LINALG::SerialDenseMatrix mortar_M;

  // Flag if pair has a active contribution.
  bool mortar_is_active = false;

  for (auto& elepairptr : contact_pairs)
  {
    if (elepairptr->IsMortar())
    {
      // Check if this pair contributes to the mortar matrices.
      mortar_is_active = elepairptr->EvaluateDM(mortar_D, mortar_M);
    }
  }
}
