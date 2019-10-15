/*----------------------------------------------------------------------*/
/*! \file
\brief Performs ScaTra specifc functions not yet generalized for other fields.

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../linalg/linalg_utils.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"

#include "scatra_utils.H"


/*----------------------------------------------------------------------*
 | Calculate the Mean-average of gradient of a scalar       winter 04/17|
 *----------------------------------------------------------------------*/
template <const int dim>
Teuchos::RCP<Epetra_MultiVector> SCATRA::SCATRAUTILS::ComputeGradientAtNodesMeanAverage(
    Teuchos::RCP<DRT::Discretization> discret, const Teuchos::RCP<const Epetra_Vector> state,
    const int scatra_dofid)
{
  // number space dimensions
  const size_t nsd = dim;
  // const int scatra_dofid = 0; //<  this is the first DoFSet (i.e. the scalar one!!)

  if (nsd != 3) dserror("Only implemented for 3D elements. Should be simple enough to extend...");

  // DOF-COL-MAP
  const Teuchos::RCP<Epetra_Vector> phinp_col =
      Teuchos::rcp(new Epetra_Vector(*discret->DofColMap()));
  // export DofRowMap to DofColMap phinp
  LINALG::Export(*state, *phinp_col);

  // ---------------------------------------------------------------

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_MultiVector> gradphirow =
      Teuchos::rcp(new Epetra_MultiVector(*discret->DofRowMap(), nsd));
  gradphirow->PutScalar(0.0);

  // map of pointers to nodes which must be reconstructed by this processor <local id, node>
  std::map<int, const DRT::Node*> nodesToReconstruct;
  nodesToReconstruct.clear();

  //--------------------------------------------------------------------------------------------
  // PART I:  loop over all elements in column map to find all nodes which must be reconstructed
  // remark: intersected elements at processor boundary must be seen by all processors
  //--------------------------------------------------------------------------------------------
  for (int iele = 0; iele < discret->NumMyColElements(); iele++)
  {
    // get element from fluid discretization
    const DRT::Element* actele = discret->lColElement(iele);

    // get number of nodes of this element (number of vertices)
    const int numberOfNodes = actele->NumNode();
    // get vector of pointers of node (for this element)
    const DRT::Node* const* ele_vecOfPtsToNode = actele->Nodes();

    // loop nodes of this element
    for (int vec_it = 0; vec_it < numberOfNodes; vec_it++)
    {
      // get owner of the node to compare with my_rank
      int node_owner = (ele_vecOfPtsToNode[vec_it])->Owner();
      // check wheather this node is a row node, compare with actual processor id
      if (node_owner == discret->Comm().MyPID())
      {
        // insert in map (overwrite existing entry)
        int lid = ele_vecOfPtsToNode[vec_it]->LID();
        nodesToReconstruct[lid] = ele_vecOfPtsToNode[vec_it];
      }
      // remark: all non-row nodes are reconstructed by another processor
    }
  }  // end for loop over row elements

  //-------------------------------------------------------------------------------------------
  // PART II: reconstruct nodes
  // remark: row nodes adjacent to intersected elements must be reconstructed by this processor
  //-------------------------------------------------------------------------------------------

  typedef std::map<int, const DRT::Node*> Map_NodesToReconstruct;
  typedef Map_NodesToReconstruct::iterator Reconstruct_iterator;

  // loop over all nodes inserted into map 'nodesToReconstruct'
  for (Reconstruct_iterator it_node = nodesToReconstruct.begin();
       it_node != nodesToReconstruct.end(); it_node++)
  {
    static LINALG::Matrix<nsd, 1>
        node_gradphi_smoothed;  // set whole 3D vector also for 2D examples
    node_gradphi_smoothed.Clear();

    // get local processor id of current node and pointer to current node
    // int lid_node = it_node->first;
    const DRT::Node* ptToNode = it_node->second;
    const int nodegid = ptToNode->Id();

    // vector of elements located around this node
    std::vector<const DRT::Element*> elements;

    // get adjacent elements for this node
    const DRT::Element* const* adjelements = ptToNode->Elements();

    const DRT::Element::DiscretizationType DISTYPE =
        adjelements[0]->Shape();  // DRT::Element::hex8;

    for (int iele = 0; iele < ptToNode->NumElement(); iele++)
    {
      if (DISTYPE != adjelements[iele]->Shape())
        dserror("Discretization not with same elements!!!");

      elements.push_back(adjelements[iele]);
    }

    // -------------------------------------------------------------
    //--------------------------------------
    // add elements along perodic boundaries
    //--------------------------------------
    // boolean indicating whether this node is a pbc node
    bool pbcnode = false;
    std::set<int> coupnodegid;
    // loop all nodes with periodic boundary conditions (master nodes)
    Teuchos::RCP<std::map<int, std::vector<int>>> pbccolmap = discret->GetAllPBCCoupledColNodes();
    for (std::map<int, std::vector<int>>::const_iterator pbciter = (*pbccolmap).begin();
         pbciter != (*pbccolmap).end(); ++pbciter)
    {
      if (pbciter->first == nodegid)  // node is a pbc master node
      {
        pbcnode = true;
        // coupled node is the slave node; there can be more than one per master node
        for (unsigned int i = 0; i < pbciter->second.size(); i++)
          coupnodegid.insert(pbciter->second[i]);
      }
      else
      {
        // loop all slave nodes
        for (size_t islave = 0; islave < pbciter->second.size(); islave++)
        {
          if (pbciter->second[islave] == nodegid)  // node is a pbc slave node
          {
            pbcnode = true;
            // coupled node is the master node
            coupnodegid.insert(pbciter->first);

            // there can be multiple slaves -> add all other slaves
            for (size_t i = 0; i < pbciter->second.size(); ++i)
            {
              if (pbciter->second[islave] != pbciter->second[i])
                coupnodegid.insert(pbciter->second[i]);
            }
          }
        }
      }
    }

    // add elements located around the coupled pbc node
    if (pbcnode)
    {
      for (std::set<int>::const_iterator icnode = coupnodegid.begin(); icnode != coupnodegid.end();
           ++icnode)
      {
        // get coupled pbc node (master or slave)
        const DRT::Node* ptToCoupNode = discret->gNode(*icnode);
        // get adjacent elements of this node
        const DRT::Element* const* pbcelements = ptToCoupNode->Elements();
        // add elements to list
        for (int iele = 0; iele < ptToCoupNode->NumElement(); iele++)  // = ptToNode->Elements();
        {
          elements.push_back(pbcelements[iele]);
        }
      }
    }

    // -------------------------------------------------------------

    // FOR OTHER TYPES OF ELEMENT THAN HEX -> One needs to change DoMeanValueAveraging - code as
    // well.
    if (DISTYPE == DRT::Element::hex8)
      node_gradphi_smoothed = DoMeanValueAveragingOfElementGradientNode<nsd, DRT::Element::hex8>(
          discret, elements, phinp_col, nodegid, scatra_dofid);
    else if (DISTYPE == DRT::Element::hex27)
      node_gradphi_smoothed = DoMeanValueAveragingOfElementGradientNode<nsd, DRT::Element::hex27>(
          discret, elements, phinp_col, nodegid, scatra_dofid);
    else
      dserror("Element type not supported yet!");

    //----------------------------------------------------------------------------------------------------
    // set the global vector gradphirow holding the new reconstructed values of gradient of phi in
    // row map
    //----------------------------------------------------------------------------------------------------

    const std::vector<int> lm = discret->Dof(scatra_dofid, (it_node->second));
    if (lm.size() != 1) dserror("assume a unique level-set dof in ScaTra DoFset");

    int GID = lm[0];  // Global ID of DoF Map
    // get local processor id according to global node id
    const int lid = (*gradphirow).Map().LID(GID);
    if (lid < 0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector", (*gradphirow).Comm().MyPID(), GID);

    const int numcol = (*gradphirow).NumVectors();
    if (numcol != (int)nsd)
      dserror("number of columns in Epetra_MultiVector is not identically to nsd");

    // loop over dimensions (= number of columns in multivector)
    for (int col = 0; col < numcol; col++)
    {
      // get columns vector of multivector
      double* globalcolumn = (*gradphirow)[col];

      // set smoothed gradient entry of phi into column of global multivector
      globalcolumn[lid] = node_gradphi_smoothed(col, 0);
    }
  }  // end loop over nodes

  return gradphirow;
}

template Teuchos::RCP<Epetra_MultiVector> SCATRA::SCATRAUTILS::ComputeGradientAtNodesMeanAverage<3>(
    Teuchos::RCP<DRT::Discretization> discret, const Teuchos::RCP<const Epetra_Vector> state,
    const int scatra_dofid);



template <const int dim, DRT::Element::DiscretizationType DISTYPE>
LINALG::Matrix<dim, 1> SCATRA::SCATRAUTILS::DoMeanValueAveragingOfElementGradientNode(
    Teuchos::RCP<DRT::Discretization> discret, std::vector<const DRT::Element*> elements,
    Teuchos::RCP<Epetra_Vector> phinp_node, const int nodegid, const int scatra_dofid)
{
  const size_t nsd = dim;
  // number of nodes of this element for interpolation
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  LINALG::Matrix<dim, 1> node_gradphi_smoothed(true);

  // number of elements located around this node
  const int numberOfElements = elements.size();
  {
    //--------------------------------------------------------------------------
    // average (mean value) reconstruction for boundary nodes and as alternative
    //--------------------------------------------------------------------------
    // loop over elements located around this node
    for (int ele_current = 0; ele_current < numberOfElements; ele_current++)
    {
      // get current element
      const DRT::Element* ele_adj = elements[ele_current];

      const int* ptToNodeIds_adj = ele_adj->NodeIds();

      // get phi-values of current adjacent element ele_adj
      // create vector "ephinp" holding scalar phi values for this element
      Epetra_SerialDenseVector ephinp(numnode);  // local vector phi-values of adjacent element

      // which node in param space of element ele_adj has actnode
      int ID_param_space = -1;

      // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
      std::vector<int> nodeID_adj(numnode);
      std::vector<int> nodeDOFID_adj(numnode);
      for (size_t inode = 0; inode < numnode; inode++)
      {
        nodeID_adj[inode] = ptToNodeIds_adj[inode];

        const std::vector<int> lm = discret->Dof(scatra_dofid, (ele_adj->Nodes()[inode]));
        if (lm.size() != 1) dserror("assume a unique level-set dof in cutterdis-Dofset");
        nodeDOFID_adj[inode] = lm[0];

        // get local number of node actnode in ele_adj
        if (nodegid == ptToNodeIds_adj[inode]) ID_param_space = inode;
      }
      if (ID_param_space < 0) dserror("node not found in element");

      // extract the phi-values of adjacent element with local ids from global vector *phinp
      // get pointer to vector holding G-function values at the fluid nodes
      DRT::UTILS::ExtractMyValues(*phinp_node, ephinp, nodeDOFID_adj);
      LINALG::Matrix<numnode, 1> ephi_adj(ephinp);

      //-------------------------------------
      // compute gradient of phi at this node
      //-------------------------------------
      // get derivatives of shape functions evaluated at node in XYZ-coordinates
      static LINALG::Matrix<nsd, numnode> deriv3Dele_xyz;
      // get derivatives of shape functions evaluates at node in Xi-coordinates
      static LINALG::Matrix<nsd, numnode> deriv3Dele;

      // TODO: Implement for other elements than HEX
      // get Xi-coordinates of current node in current adjacent element
      static LINALG::Matrix<nsd, 1> node_Xicoordinates;
      node_Xicoordinates.Clear();
      for (size_t icomp = 0; icomp < nsd; ++icomp)
      {
        node_Xicoordinates(icomp) =
            DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[ID_param_space][icomp];
      }

      // get derivatives of shape functions at node
      switch (nsd)
      {
        case 3:
        {
          DRT::UTILS::shape_function_3D_deriv1(deriv3Dele, node_Xicoordinates(0),
              node_Xicoordinates(1), node_Xicoordinates(2), DISTYPE);
          break;
        }
        case 2:
        {
          DRT::UTILS::shape_function_2D_deriv1(
              deriv3Dele, node_Xicoordinates(0), node_Xicoordinates(1), DISTYPE);
          break;
        }
        case 1:
        {
          DRT::UTILS::shape_function_1D_deriv1(deriv3Dele, node_Xicoordinates(0), DISTYPE);
          break;
        }
        default:
          dserror("Only spacial dimension 1,2,3 are allowed!");
      }

      // reconstruct XYZ-gradient
      // get node coordinates of this element
      static LINALG::Matrix<nsd, numnode> xyze_adj;
      GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

      // get Jacobi-Matrix for transformation
      static LINALG::Matrix<nsd, nsd> xjm_ele_XiToXYZ;
      xjm_ele_XiToXYZ.MultiplyNT(deriv3Dele, xyze_adj);

      // inverse of jacobian
      static LINALG::Matrix<nsd, nsd> xji_ele_XiToXYZ;
      xji_ele_XiToXYZ.Invert(xjm_ele_XiToXYZ);

      // set XYZ-derivates of shapefunctions
      deriv3Dele_xyz.Multiply(xji_ele_XiToXYZ, deriv3Dele);

      //----------------------------------------------------
      // compute gradient of phi at node for current element
      //----------------------------------------------------
      static LINALG::Matrix<nsd, 1> nodal_grad_tmp;
      nodal_grad_tmp.Clear();

      // get xyz-gradient
      nodal_grad_tmp.Multiply(deriv3Dele_xyz, ephi_adj);

      //===============================================================================
      // add to vector with smoothed vector

      node_gradphi_smoothed.Update(1.0, nodal_grad_tmp, 1.0);

    }  // end loop over all adjacent elements

    // weight sum of nodal_grad_tmp 1/number_of_vectors to get an average value
    node_gradphi_smoothed.Scale(1.0 / numberOfElements);
  }
  return node_gradphi_smoothed;
}

// Templates for Mean value averaging -- For now only HEX-type elements allowed!
template LINALG::Matrix<3, 1>
SCATRA::SCATRAUTILS::DoMeanValueAveragingOfElementGradientNode<3, DRT::Element::hex8>(
    Teuchos::RCP<DRT::Discretization> discret, std::vector<const DRT::Element*> elements,
    Teuchos::RCP<Epetra_Vector> phinp_node, const int nodegid, const int scatra_dofid);

template LINALG::Matrix<3, 1>
SCATRA::SCATRAUTILS::DoMeanValueAveragingOfElementGradientNode<3, DRT::Element::hex27>(
    Teuchos::RCP<DRT::Discretization> discret, std::vector<const DRT::Element*> elements,
    Teuchos::RCP<Epetra_Vector> phinp_node, const int nodegid, const int scatra_dofid);
