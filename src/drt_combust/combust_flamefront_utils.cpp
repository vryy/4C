/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront_utils.cpp

 \brief

  detailed description in header file combust_interface.H

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "combust_flamefront.H"
#include "combust_defines.H"
#include "combust3_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_gmsh.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_MpiComm.h>
#ifdef DEBUG
#include "combust3.H"
#endif


/*------------------------------------------------------------------------------------------------*
 | call function for different types of computation of smoothed gradient field of G-function      |
 | (level set field) for use in fluid time                                                        |
 | integration scheme                                                                 henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CallSmoothGradPhi(const Teuchos::ParameterList& combustdyn)
{
  TEUCHOS_FUNC_TIME_MONITOR("SMOOTHING OF GRAD_PHI");
  if (gfuncdis_->Comm().MyPID()==0)
    IO::cout << "---  compute smoothed gradient of phi... " << IO::flush;

  // get type of reconstruction
  const INPAR::COMBUST::SmoothGradPhi SmoothGradPhi = DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI");


  // standard dimension is 3, real dimension for least squares reconstruction type is necessary

  size_t nsd_real = 3;
  if( SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx ||
      SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy ||
      SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz ) nsd_real = 2;

  // IO::cout reconstruction type
  if (fluiddis_->Comm().MyPID() == 0)
  {
    switch (SmoothGradPhi)
    {
    case INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx:
      IO::cout << "\n---  \t reconstruction with:\t LeastSquares_2Dx... " << IO::flush;
      break;
    case INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy:
      IO::cout << "\n---  \t reconstruction with:\t LeastSquares_2Dy... " << IO::flush;
      break;
    case INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz:
      IO::cout << "\n---  \t reconstruction with:\t LeastSquares_2Dz... " << IO::flush;
      break;
    case INPAR::COMBUST::smooth_grad_phi_leastsquares_3D:
      IO::cout << "\n---  \t reconstruction with:\t LeastSquares_3D... " << IO::flush;
      break;
    case INPAR::COMBUST::smooth_grad_phi_meanvalue:
      IO::cout << "\n---  \t reconstruction with:\t MeanValue... " << IO::flush;
      break;
    case INPAR::COMBUST::smooth_grad_phi_l2_projection:
      IO::cout << "\n---  \t reconstruction with:\t L2Projection... " << IO::flush;
      break;
    default:
      dserror("Smoothing strategy for gradient of phi expected!");
      break;
    }
  }

  if (SmoothGradPhi != INPAR::COMBUST::smooth_grad_phi_l2_projection)
  {
    // switch: real dimension of reconstruction, is real dimension 2D or 3D?
    switch (nsd_real)
    {
      case 3:
      ComputeSmoothGradPhi<3>(combustdyn);
      break;
      case 2:
      ComputeSmoothGradPhi<2>(combustdyn);
      break;
      default:
      dserror("wrong nsd_real");
      break;
    }
  }
  else
  {
    // get diffusion for additional smoothing
    const double eta_smooth = (combustdyn.sublist("COMBUSTION FLUID").get<double>("SMOOTHING_PARAMETER"));

    ComputeL2ProjectedGradPhi(eta_smooth);
    // currently unused to also smooth phi after reinitialization via signed distance function
//    const Teuchos::RCP<Epetra_Vector> phi_smoothed = ComputeL2ProjectedPhi(eta_smooth);
//    ComputeL2ProjectedGradPhi(eta_smooth, phi_smoothed);
    if (DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION FLUID"),"L2_PROJECTION_SECOND_DERIVATIVES"))
      //L2_Projection_Gradphi2();
      ComputeL2ProjectedGrad2Phi(eta_smooth);
    else
    {
      // set gradphi2 to zero
      gradphi2_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),9));
      gradphi2_->PutScalar(0.0);
    }
  }

  if (gfuncdis_->Comm().MyPID()==0)
    IO::cout << "done" << IO::endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute smoothed gradient field of G-function (level set field) for use in fluid time          |
 | integration scheme                                                                schott 06/10 |
 *------------------------------------------------------------------------------------------------*/
template <const size_t nsd_real>
void COMBUST::FlameFront::ComputeSmoothGradPhi(const Teuchos::ParameterList& combustdyn)
{
  // get type of reconstruction
  const INPAR::COMBUST::SmoothGradPhi SmoothGradPhi = DRT::INPUT::IntegralValue<INPAR::COMBUST::SmoothGradPhi>
  (combustdyn.sublist("COMBUSTION FLUID"),"SMOOTHGRADPHI");

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8;

  // number of nodes of this element for interpolation
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // gradphi_ gets a NodeColMap
  gradphi_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),nsd));
  gradphi_->PutScalar(0.0);
  // set also gradphi2
  gradphi2_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),3*nsd));
  gradphi2_->PutScalar(0.0);

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_MultiVector> gradphirow = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeRowMap(),nsd));
  gradphirow->PutScalar(0.0);

  // map of pointers to nodes which must be reconstructed by this processor <local id, node>
  std::map<int, const DRT::Node*> nodesToReconstruct;
  nodesToReconstruct.clear();

  //--------------------------------------------------------------------------------------------
  // PART I:  loop over all elements in column map to find all nodes which must be reconstructed
  // remark: intersected elements at processor boundary must be seen by all processors
  //--------------------------------------------------------------------------------------------
  for(int iele=0;iele<fluiddis_->NumMyColElements();iele++)
  {
    // get element from fluid discretization
    const DRT::Element *actele = fluiddis_->lColElement(iele);
    // get global Id of element
    //int actele_GID = actele->Id();
    // get vector of all integration cells of this element
    //GEO::DomainIntCells IntCells = myelementintcells_[actele_GID];
    // get vector of all boundary integration cells of this element -> also touched elements need reconstruction
    //GEO::BoundaryIntCells BoundIntCells = myboundaryintcells_[actele_GID];
    //    IO::cout << "actele->ID()" << actele->Id() << IO::endl;
    //    IO::cout << actele->Id() << "\t" << IntCells.size() << IO::endl;
    //    IO::cout << actele->Id() << "\t" <<BoundIntCells.size() << IO::endl;


    // TODO check if this can be put back in
    // find out whether this element is intersected or not by looking at number of integration cells
    // remark: this might include some uncut elements, if the GEO::Cut algorithm created domain
    //         integration cells for those elements
    //if (IntCells.size() > 1 || BoundIntCells.size() > 0) // element is intersected or touched
    //{
    //  IO::cout << "element with ID: " << actele->Id() << "is reconstructed "<< IO::endl;
    //-> assemble all adjacent nodes, these node values must be reconstructed

    // get number of nodes of this element (number of vertices)
    // remark: e.g. hex20 elements has GEO::DomainIntCells IntCells = myelementintcells_[actele_GID]; numnode=20 but numberOfNodes=8
    //         if actele is an element at processor boundary, actele will be less than 8
    const int numberOfNodes = actele->NumNode();
    // get vector of pointers of node (for this element)
    const DRT::Node* const* ele_vecOfPtsToNode = actele->Nodes();

    // loop nodes of this element
    for(int vec_it = 0; vec_it < numberOfNodes; vec_it++)
    {
      // get owner of the node to compare with my_rank
      int node_owner = (ele_vecOfPtsToNode[vec_it])->Owner();
      // check wheather this node is a row node, compare with actual processor id
      if(node_owner == fluiddis_->Comm().MyPID() )
      {
        // insert in map (overwrite existing entry)
        int lid = ele_vecOfPtsToNode[vec_it]->LID();
        nodesToReconstruct[lid] = ele_vecOfPtsToNode[vec_it];
      }
      // remark: all non-row nodes are reconstructed by another processor
    }

    //}
  }// end for loop over row elements


  //-------------------------------------------------------------------------------------------
  // PART II: reconstruct nodes
  // remark: row nodes adjacent to intersected elements must be reconstructed by this processor
  //-------------------------------------------------------------------------------------------

  typedef std::map<int, const DRT::Node*> Map_NodesToReconstruct;
  typedef Map_NodesToReconstruct::iterator Reconstruct_iterator;

  // loop over all nodes inserted into map 'nodesToReconstruct'
  for(Reconstruct_iterator it_node = nodesToReconstruct.begin(); it_node!= nodesToReconstruct.end(); it_node++)
  {
    // define vector for smoothed values (Phi, grad_Phi) of current node
    static LINALG::Matrix<nsd_real,1> PHI_SMOOTHED; // for real reconstruction 2D or 3D
    PHI_SMOOTHED.Clear();
    static LINALG::Matrix<nsd,1> PHI_SMOOTHED_3D; // set whole 3D vector also for 2D examples
    PHI_SMOOTHED_3D.Clear();

    // get local processor id of current node and pointer to current node
    //int lid_node = it_node->first;
    const DRT::Node* ptToNode = it_node->second;
    const int nodegid = ptToNode->Id();

    // vector of elements located around this node
    std::vector<const DRT::Element*> elements;

    // get adjacent elements for this node
    const DRT::Element*const* adjelements = ptToNode->Elements();
    for (int iele=0;iele<ptToNode->NumElement();iele++)
    {
      //const DRT::Element* ele = elements[iele];
      elements.push_back(adjelements[iele]);
    }

    //--------------------------------------
    // add elements along perodic boundaries
    //--------------------------------------
    // boolean indicating whether this node is a pbc node
    bool pbcnode = false;
    std::set<int> coupnodegid;
    // loop all nodes with periodic boundary conditions (master nodes)
    Teuchos::RCP<std::map<int,std::vector<int> > > pbccolmap = gfuncdis_->GetAllPBCCoupledColNodes();
    for (std::map<int, std::vector<int>  >::const_iterator pbciter= (*pbccolmap).begin(); pbciter != (*pbccolmap).end(); ++pbciter)
    {
      if (pbciter->first == nodegid) // node is a pbc master node
      {
        pbcnode = true;
        // coupled node is the slave node; there can be more than one per master node
        for (unsigned int i = 0; i < pbciter->second.size(); i++)
          coupnodegid.insert(pbciter->second[i]);
      }
      else
      {
        // loop all slave nodes
        for (size_t islave=0;islave<pbciter->second.size();islave++)
        {
          if (pbciter->second[islave] == nodegid) // node is a pbc slave node
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
      for (std::set<int>::const_iterator icnode = coupnodegid.begin(); icnode != coupnodegid.end(); ++icnode)
      {
        // get coupled pbc node (master or slave)
        const DRT::Node* ptToCoupNode = gfuncdis_->gNode(*icnode);
        // get adjacent elements of this node
        const DRT::Element*const* pbcelements = ptToCoupNode->Elements();
        // add elements to list
        for (int iele=0;iele<ptToCoupNode->NumElement();iele++)// = ptToNode->Elements();
        {
          elements.push_back(pbcelements[iele]);
        }
      }
    }

    // number of elements located around this node
    const int numberOfElements = elements.size();

    //--------------------------------------------
    // set type of reconstruction for current node
    //--------------------------------------------
    // initialize reconstruction type for current node with 'none'
    INPAR::COMBUST::SmoothGradPhi CurrTypeOfNodeReconst = INPAR::COMBUST::smooth_grad_phi_none;

    // set reconstruction node for current node with special cases (boundary elements)
    if (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_none)
      dserror("SmoothGradPhi() should not be called for reconstruction type Grad_phi_None");
    else if (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_meanvalue)
      CurrTypeOfNodeReconst = INPAR::COMBUST::smooth_grad_phi_meanvalue;
    else if (
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D  && numberOfElements <  int(nsd_real)) ||
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx && numberOfElements <= int(nsd_real)) ||
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy && numberOfElements <= int(nsd_real)) ||
        (SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz && numberOfElements <= int(nsd_real)) )
    {
      CurrTypeOfNodeReconst = INPAR::COMBUST::smooth_grad_phi_meanvalue;
      // remark: although least squares reconstruction is chosen in input file, boundary elements
      //         have to be reconstructed by the mean value (average) method because there are not
      //         enough adjacent elements to reconstruct the nodal phi gradient minimal number of
      //         adjacent elements for least_squares reconstruction is 3 for 3D and 2 for 2D,
      //         so nsd_real mean value method is the same for 2D and 3D
//      IO::cout << "/!\\" << "\t\tnode\t" << (it_node->second)->Id() << "\tis a boundary node with "
//           << numberOfElements << " elements and is reconstructed with average (mean value) method" << IO::endl;
    }
    else // type of reconstruction as chosen in input file
      CurrTypeOfNodeReconst = SmoothGradPhi;

    // remark:
    // special case for intersected elements which have boundary nodes
    // for boundary nodes there are not enough elements to reconstruct node gradient with least squares
    // => build average of gradients of adjacent elements instead of least squares reconstruction
    // at least nsd+1 adjacent elements necessary for least squares reconstruction, we take at least numnode!!!
    if(CurrTypeOfNodeReconst == INPAR::COMBUST::smooth_grad_phi_meanvalue)  //(int)numnode
    {
      //--------------------------------------------------------------------------
      // average (mean value) reconstruction for boundary nodes and as alternative
      //--------------------------------------------------------------------------
      // loop over elements located around this node
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get current element
        const DRT::Element* ele_adj = elements[ele_current];

        const int* ptToNodeIds_adj = ele_adj->NodeIds();

        // get phi-values of current adjacent element ele_adj
        // create vector "ephinp" holding scalar phi values for this element
        Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element

        // which node in param space of element ele_adj has actnode
        int ID_param_space = -1;

        // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
        std::vector<int> nodeID_adj(numnode);
        for (size_t inode=0; inode < numnode; inode++)
        {
          nodeID_adj[inode] = ptToNodeIds_adj[inode];
          // get local number of node actnode in ele_adj
          if(nodegid == ptToNodeIds_adj[inode]) ID_param_space = inode;
        }
        // node not found in this element -> must be a pbc element
        if (ID_param_space < 0)
        {
          // check if coupled pbc node is found instead
          for (size_t inode=0; inode < numnode; inode++)
          {
            nodeID_adj[inode] = ptToNodeIds_adj[inode];
            // get local number of node actnode in ele_adj
            if (coupnodegid.find(ptToNodeIds_adj[inode]) != coupnodegid.end())
              ID_param_space = inode;
          }
        }
        if (ID_param_space < 0) dserror("node not found in element");

        // extract the phi-values of adjacent element with local ids from global vector *phinp
        // get pointer to vector holding G-function values at the fluid nodes
        DRT::UTILS::ExtractMyValues(*phinp_, ephinp, nodeID_adj);
        LINALG::Matrix<numnode,1> ephi_adj(ephinp);

        //-------------------------------------
        // compute gradient of phi at this node
        //-------------------------------------
        // get derivatives of shape functions evaluated at node in XYZ-coordinates
        static LINALG::Matrix<nsd,numnode> deriv3Dele_xyz;
        // get derivatives of shape functions evaluates at node in Xi-coordinates
        static LINALG::Matrix<nsd,numnode> deriv3Dele;

        // get Xi-coordinates of current node in current adjacent element
        static LINALG::Matrix<nsd,1> node_Xicoordinates;
        node_Xicoordinates.Clear();
        for(size_t icomp = 0; icomp<nsd; ++icomp)
        {
          node_Xicoordinates(icomp) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[ID_param_space][icomp];
        }

        // get derivatives of shapefunctions at node
        DRT::UTILS::shape_function_3D_deriv1(deriv3Dele,node_Xicoordinates(0),node_Xicoordinates(1),node_Xicoordinates(2),DISTYPE);

        // reconstruct XYZ-gradient

        // get node coordinates of this element
        static LINALG::Matrix<nsd,numnode> xyze_adj;
        GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

        // get Jacobi-Matrix for transformation
        static LINALG::Matrix<nsd,nsd> xjm_ele_XiToXYZ;
        xjm_ele_XiToXYZ.MultiplyNT(deriv3Dele,xyze_adj);

        // inverse of jacobian
        static LINALG::Matrix<nsd,nsd> xji_ele_XiToXYZ;
        xji_ele_XiToXYZ.Invert(xjm_ele_XiToXYZ);

        // set XYZ-derivates of shapefunctions
        deriv3Dele_xyz.Multiply(xji_ele_XiToXYZ,deriv3Dele);

        //----------------------------------------------------
        // compute gradient of phi at node for current element
        //----------------------------------------------------
        static LINALG::Matrix<3,1> nodal_grad_tmp;
        nodal_grad_tmp.Clear();

        // get xyz-gradient
        nodal_grad_tmp.Multiply(deriv3Dele_xyz, ephi_adj);

        //===============================================================================
        // add to vector with smoothed vector
        // we dont need PHI_SMOOTHED(0,0)

        // full 3D reconstruction, also for 2D examples
        PHI_SMOOTHED_3D(0,0) += nodal_grad_tmp(0,0);
        PHI_SMOOTHED_3D(1,0) += nodal_grad_tmp(1,0);
        PHI_SMOOTHED_3D(2,0) += nodal_grad_tmp(2,0);

      }// end loop over all adjacent elements

#if(0)
      // special case for straight_bodyforce
      // this is a node of an intersected element
      // the element must be a boundary element!
      // assume a continued interface across the domain boundary
      if(numberOfElements == 2)
      {
        // set two vectors in y-direction
        PHI_SMOOTHED_3D(0,0) = 1.0;
        PHI_SMOOTHED_3D(1,0) = 1.0;
        PHI_SMOOTHED_3D(2,0) = 0.0;

        IO::cout << "\n\t !!! warning !!! (Rayleigh-Taylor modification) we modify the gradient of phi periodically at the domain boundary";
      }
#endif

#if(0)
      // special case for Rayleigh-Taylor
      // this is a node of an intersected element
      // the element must be a boundary element!
      // assume a continued interface across the domain boundary
      if(numberOfElements == 2)
      {
        // set two vectors in y-direction
        PHI_SMOOTHED_3D(0,0) = 0.0;
        PHI_SMOOTHED_3D(1,0) = 2.0;
        PHI_SMOOTHED_3D(2,0) = 0.0;

        IO::cout << "\n\t !!! warning !!! (Rayleigh-Taylor modification) we modify the gradient of phi periodically at the domain boundary";
      }
#endif

#ifdef COMBUST_2D
      //IO::cout << "/!\\ third component or normal vector set to 0 to keep 2D-character!" << IO::endl;
      PHI_SMOOTHED_3D(2,0) = 0.0;
#endif

      // weight sum of nodal_grad_tmp 1/number_of_vectors to get an average value
      PHI_SMOOTHED_3D.Scale(1.0/numberOfElements);

    } // end of average (mean value) reconstruction
    else // least squares reconstruction
    {
      //----------------------------------------------------------------------
      // standard reconstruction with least squares method, see Merchandise'07
      //----------------------------------------------------------------------

      // we need Epetra-Matrix (numberOfElements is not a valid template parameter)

      Epetra_SerialDenseMatrix RHS_LS(numberOfElements,1);
      Epetra_SerialDenseMatrix MAT_LS(numberOfElements, nsd_real);

      int xyz_2D_dim = -1;
      // set one dimension (xyz-dim) to zero for 2D LS-reconstruction
      if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dx) xyz_2D_dim = 0;
      if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dy) xyz_2D_dim = 1;
      if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_2Dz) xyz_2D_dim = 2;

      // loop over adjacent elements
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get current element
        const DRT::Element* ele_adj = elements[ele_current];

        const int* ptToNodeIds_adj = ele_adj->NodeIds();

        // get phi-values of current adjacent element ele_adj
        // create vector "ephinp" holding scalar phi values for this element
        Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element

        // which node in param space of element ele_adj has actnode
        int ID_param_space = -1;

        // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
        std::vector<int> nodeID_adj(numnode);
        for (size_t inode=0; inode < numnode; inode++){
          nodeID_adj[inode] = ptToNodeIds_adj[inode];
          // get local number of node actnode in ele_adj
          if(nodegid == ptToNodeIds_adj[inode]) ID_param_space = inode;
        }
        // node not found in this element -> must be a pbc element
        if (ID_param_space < 0)
        {
          // check if coupled pbc node is found instead
          for (size_t inode=0; inode < numnode; inode++)
          {
            nodeID_adj[inode] = ptToNodeIds_adj[inode];
            // get local number of node actnode in ele_adj
            if (coupnodegid.find(ptToNodeIds_adj[inode]) != coupnodegid.end())
              ID_param_space = inode;
          }
        }
        if (ID_param_space < 0) dserror ("node in current adjacent not found!!!");

        // extract the phi-values of adjacent element with local ids from global vector *phinp
        // get pointer to vector holding G-function values at the fluid nodes
        DRT::UTILS::ExtractMyValues(*phinp_, ephinp, nodeID_adj);
        LINALG::Matrix<numnode,1> ephi_adj(ephinp);

        // calculate center of gravity
        static LINALG::Matrix<numnode,1> funct;
        funct.Clear();
        static LINALG::Matrix<nsd,1> centerOfGravXi;
        centerOfGravXi.Clear();

        // xi-coord = 0.0.0 for hex8
        if (DISTYPE != DRT::Element::hex8) dserror("center of gravity implemented only for hex8 elements");
        centerOfGravXi(0) = 0.0;
        centerOfGravXi(1) = 0.0;
        centerOfGravXi(2) = 0.0;

        DRT::UTILS::shape_function_3D(funct,centerOfGravXi(0),centerOfGravXi(1),centerOfGravXi(2),DISTYPE);

        // get node coordinates of this element
        static LINALG::Matrix<nsd,numnode> xyze_adj;
        GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

        // interpolate center of gravity
        static LINALG::Matrix<nsd,1> centerOfGravXYZ;
        centerOfGravXYZ.Clear();

        centerOfGravXYZ.Multiply(xyze_adj,funct);

        // reconstruct xyz-gradient
        // calculate vector from current node to point of gravity (direction for taylor)
        static LINALG::Matrix<nsd,1> direction;
        direction.Clear();

        direction(0) = centerOfGravXYZ(0) - xyze_adj(0, ID_param_space);
        direction(1) = centerOfGravXYZ(1) - xyze_adj(1, ID_param_space);
        direction(2) = centerOfGravXYZ(2) - xyze_adj(2, ID_param_space);

        // cancel out the 2D direction
        if(xyz_2D_dim != -1) direction(xyz_2D_dim) = 0.0;

        // assemble direction into matrix A

        // calculate ephi at point of gravity via interpolation
        static LINALG::Matrix<1,1> phi_adj;
        phi_adj.Clear();
        phi_adj.MultiplyTN(ephi_adj,funct);

        // set RHS and MAT for least squares method

        if(SmoothGradPhi == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D){
          RHS_LS(ele_current,0) = phi_adj(0,0)-ephi_adj(ID_param_space,0);
          MAT_LS(ele_current,0) = direction(0);
          MAT_LS(ele_current,1) = direction(1);
          MAT_LS(ele_current,2) = direction(2);
        }
        else{ // assemble MAT_LS for 2D cases
          RHS_LS(ele_current,0) = phi_adj(0,0)-ephi_adj(ID_param_space,0);

          // 2D LEAST_SQUARES cases
          size_t MAT_isd = 0;
          for(size_t isd = 0; isd<nsd; isd++)
          {
            if(int(isd)!=xyz_2D_dim){// we have to not!!! assemble the xyz_2D_dim
              MAT_LS(ele_current,MAT_isd) = direction(isd);
              MAT_isd++;
            }
          } // end for
        } // end else

      } // end loop over all adjacent elements

      // the system MAT_LS * phi_smoothed = RHS_LS is only solvable in a least squares manner
      // MAT_LS is not square
      // -> solve MAT_LS^T * MAT_LS * phi_smoothed = MAT_LS^T * RHS_LS
      Epetra_SerialDenseMatrix MAT_tmp(nsd_real,nsd_real);
      Epetra_SerialDenseMatrix RHS_tmp(nsd_real,1);

      MAT_tmp.Multiply('T','N', 1.0, MAT_LS,MAT_LS, 0.0);
      RHS_tmp.Multiply('T','N', 1.0, MAT_LS,RHS_LS, 0.0);

      // set LINALG-Matrix with fixed size
      // initialize element matrix for A^T*A and RHS-vector to solve Normalengleichung

      LINALG::Matrix<nsd_real,nsd_real> MAT(MAT_tmp);
      LINALG::Matrix<nsd_real,1> RHS(RHS_tmp);

      //----------------------------------------------------
      // solve A^T*A* phi_smoothed = A^T*RHS (LEAST SQUARES)
      //----------------------------------------------------
      // solve the system for current node
      LINALG::FixedSizeSerialDenseSolver<nsd_real,nsd_real,1> solver; //1 is dimension of RHS
      solver.SetMatrix(MAT);
      solver.SetVectors(PHI_SMOOTHED,RHS);
      solver.Solve();

      //      IO::cout << " node:" << (it_node->second)->Id() << IO::endl;
      //      IO::cout << "PHI_SMOOTHED" << PHI_SMOOTHED << IO::endl;
      // set full 3D vector especially important for 2D leastsquares reconstruction
      if(CurrTypeOfNodeReconst == INPAR::COMBUST::smooth_grad_phi_leastsquares_3D)
      {
        PHI_SMOOTHED_3D(0,0) = PHI_SMOOTHED(0,0);
        PHI_SMOOTHED_3D(1,0) = PHI_SMOOTHED(1,0);
        PHI_SMOOTHED_3D(2,0) = PHI_SMOOTHED(2,0);
      }
      else // 2D-LeastSquares-cases
      {
        size_t PHI_2D_isd = 0;
        for(size_t isd = 0; isd<nsd; isd++)
        {
          if(int(isd)!=xyz_2D_dim){ // we have to not!!! assemble the xyz_2D_dim
            PHI_SMOOTHED_3D(isd,0) = PHI_SMOOTHED(PHI_2D_isd,0);
            PHI_2D_isd++;
          }
        }
      }

    } // end of standard (Least_squares-) reconstruction case


    //----------------------------------------------------------------------------------------------------
    // set the global vector gradphirow holding the new reconstructed values of gradient of phi in row map
    //----------------------------------------------------------------------------------------------------

    // get the global id for current node
    int GID = (it_node->second)->Id();

    // get local processor id according to global node id
    const int lid = (*gradphirow).Map().LID(GID);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphirow).Comm().MyPID(),GID);

    const int numcol = (*gradphirow).NumVectors();
    if( numcol != (int)nsd) dserror("number of columns in Epetra_MultiVector is not identically to nsd");

    // loop over dimensions (= number of columns in multivector)
    for(int col=0; col< numcol; col++)
    {
      // get columns vector of multivector
      double* globalcolumn = (*gradphirow)[col];

      // set smoothed gradient entry of phi into column of global multivector
      globalcolumn[lid] = PHI_SMOOTHED_3D(col,0);
    }

  }// end loop over nodes

  // export NodeRowMap to NodeColMap gradphi_
  LINALG::Export(*gradphirow,*gradphi_);

#if 0 // only for debug
  // additional gmsh output
  {
    // turn on/off screen output for writing process of Gmsh postprocessing file
    const bool screen_out = true;

    // create Gmsh postprocessing file
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("gradient_phi", 0, 10, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Grad Phi \" {" << std::endl;
      // draw scalar field 'Phinp' for every element
      IO::GMSH::VectorFieldNodeBasedToGmsh(fluiddis_,gradphirow,gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
  }
#endif

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute nodal curvature based on L2 projection                                 rasthofer 08/13 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ComputeL2ProjectedNodalCurvature(const Teuchos::ParameterList& combustdyn)
{

  const bool L2ProjSecDerivatives = DRT::INPUT::IntegralValue<int>
  (combustdyn.sublist("COMBUSTION FLUID"),"L2_PROJECTION_SECOND_DERIVATIVES");

  // get smoothing parameter
  const double eta_smooth = (combustdyn.sublist("COMBUSTION FLUID").get<double>("SMOOTHING_PARAMETER"));

  //get the right number of space dimensions
  const size_t nsd=3;

  // initialize curvature
  curvature_ = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap(),true));
  curvature_->PutScalar(0.0);

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_Vector> curvaturerow = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap(),true));
  curvaturerow->PutScalar(0.0);

  // safety check
  if (!gfuncdis_->Filled())   dserror("FillComplete() was not called");
  if (!gfuncdis_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  //get dof row map of scalar field
  const Epetra_Map* dofrowmap = gfuncdis_->DofRowMap();

  // safety check
  if (not fluiddis_->NodeRowMap()->PointSameAs(*(gfuncdis_->NodeRowMap())))
    dserror("Same node map expected");

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> matrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  // create right hand side
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);

  // element matrix and rhs
  Epetra_SerialDenseMatrix elemat;
  Epetra_SerialDenseVector elerhs;

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // get number of elements
  const int numele = gfuncdis_->NumMyColElements();

  //loop column elements
  for (int i = 0; i< numele;i++)
  {
      DRT::Element* actele = gfuncdis_->lColElement(i);

      // get element location vector, dirichlet flags and ownerships
      lm.clear();
      lmowner.clear();
      lmstride.clear();
      actele->LocationVector(*gfuncdis_,lm,lmowner,lmstride);

      // get shape
      DRT::Element::DiscretizationType distype = actele->Shape();
      if (distype != DRT::Element::hex8) dserror("Hex 8 expected!");

      //set element data: get the right number of nodes per element
      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

      //node coordinates in global dimensions
      LINALG::Matrix<3,numnode> xyze;

      GEO::fillInitialPositionArray<DRT::Element::hex8>(actele, xyze);

      // get Gauss rule and Gauss points for integration
      const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_hex_8point);

      // Gauss point vector
      Epetra_SerialDenseVector  gp(nsd);

      // get dimension of element matrices and vectors
      // reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      if (elemat.M()!=eledim or elemat.N()!=eledim)
        elemat.Shape(eledim,eledim);
      else
        memset(elemat.A(),0,eledim*eledim*sizeof(double));

      if (elerhs.Length()!=eledim)
        elerhs.Size(eledim);
      else
        memset(elerhs.Values(),0,eledim*sizeof(double));

      // loop integration points
      for (int iquad=0; iquad<intpoints.nquad; ++iquad)
      {
        // fill vector with Gauss points
        for (size_t rr=0;rr<nsd;++rr)
          gp(rr)=intpoints.qxg[iquad][rr];

        // shape function at Gauss point
        static LINALG::Matrix<numnode,1> funct_gp;
        DRT::UTILS::shape_function_3D(funct_gp,gp(0),gp(1),gp(2),distype);

        // derivatives evaluated at Gauss point
        static LINALG::Matrix<nsd,numnode> deriv_gp;
        DRT::UTILS::shape_function_3D_deriv1(deriv_gp,gp(0),gp(1),gp(2),distype);

        // Jacobian for integration over element domain
        static LINALG::Matrix<nsd,nsd> xjm_gp;
        xjm_gp.MultiplyNT(deriv_gp,xyze);

        // invert of Jacobian for integration over element domain
        static LINALG::Matrix<nsd,nsd> xji_gp;
        const double det = xji_gp.Invert(xjm_gp);

        // global derivatives
        static LINALG::Matrix<nsd,numnode> deriv_gp_xyz;
        deriv_gp_xyz.Multiply(xji_gp,deriv_gp);

        // check for degenerated elements
        if (det < 0.0)
          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

        // integration factor
        const double fac = intpoints.qwgt[iquad]*det;

        // smoothed gradient phi vector at this node
        LINALG::Matrix<3,numnode> egradephi(true);
        // extract local (element level) G-function values from global vector
        DRT::UTILS::ExtractMyNodeBasedValues(actele, egradephi, gradphi_,3);
        // get gradient of phi at Gauss point
        static LINALG::Matrix<nsd,1> grad_phi;
        grad_phi.Clear();
        for (std::size_t rr=0; rr<nsd; rr ++)
        {
          for(size_t inode = 0; inode < numnode; inode++)
            grad_phi(rr) += egradephi(rr,inode) *funct_gp(inode);
        }

        // smoothed second derivatives of phi at this node
        LINALG::Matrix<9,numnode> egrad2phi(true);
        if (L2ProjSecDerivatives)
          DRT::UTILS::ExtractMyNodeBasedValues(actele, egrad2phi, gradphi2_,9);
        // get second derivatives of phi at Gauss point
        static LINALG::Matrix<nsd*nsd,1> grad2_phi;
        grad2_phi.Clear();

        if (L2ProjSecDerivatives)
        {
          for(size_t inode = 0; inode < numnode; inode++)
          {
            grad2_phi(0) += egrad2phi(0,inode)*funct_gp(inode); // ,xx
            grad2_phi(1) += egrad2phi(1,inode)*funct_gp(inode);  // ,yy
            grad2_phi(2) += egrad2phi(2,inode)*funct_gp(inode);  // ,zz
            grad2_phi(3) += egrad2phi(3,inode)*funct_gp(inode);  // ,xy
            grad2_phi(4) += egrad2phi(4,inode)*funct_gp(inode);  // ,xz
            grad2_phi(5) += egrad2phi(5,inode)*funct_gp(inode);  // ,yz
            grad2_phi(6) += egrad2phi(6,inode)*funct_gp(inode);  // ,yx
            grad2_phi(7) += egrad2phi(7,inode)*funct_gp(inode);  // ,zx
            grad2_phi(8) += egrad2phi(8,inode)*funct_gp(inode);  // ,zy
          }
        }
        else
        {
          for(size_t inode = 0; inode < numnode; inode++)
          {
            grad2_phi(0) += deriv_gp_xyz(0,inode)*egradephi(0,inode); // ,xx
            grad2_phi(1) += deriv_gp_xyz(1,inode)*egradephi(1,inode); // ,yy
            grad2_phi(2) += deriv_gp_xyz(2,inode)*egradephi(2,inode); // ,zz
            grad2_phi(3) += deriv_gp_xyz(1,inode)*egradephi(0,inode); // ,xy
            grad2_phi(4) += deriv_gp_xyz(2,inode)*egradephi(0,inode); // ,xz
            grad2_phi(5) += deriv_gp_xyz(2,inode)*egradephi(1,inode); // ,yz
            grad2_phi(6) += deriv_gp_xyz(0,inode)*egradephi(1,inode); // ,yx
            grad2_phi(7) += deriv_gp_xyz(0,inode)*egradephi(2,inode); // ,zx
            grad2_phi(8) += deriv_gp_xyz(1,inode)*egradephi(2,inode); // ,zy
          }
        }

        // get curvature at Gauss point
        double curvature_gp = 0.0;
        // norm of gradphi
        double grad_phi_norm = grad_phi.Norm2();
        {
          // check norm of normal gradient
          if (fabs(grad_phi_norm < 1.0E-9))
          {
            std::cout << "grad phi is small -> set to 1.0E9" << grad_phi_norm << std::endl;
            // phi gradient too small -> there must be a local max or min in the level-set field
            // set curvature to a large value
            // note: we have 1.0E12 in CalcCurvature Function
            curvature_gp = 1.0E2;
          }
          else
          {
            double val = grad_phi_norm*grad_phi_norm*grad_phi_norm;
            double invval = 1.0 / val;
            curvature_gp = -invval * ( grad_phi(0)*grad_phi(0)*grad2_phi(0)
                                      + grad_phi(1)*grad_phi(1)*grad2_phi(1)
                                      + grad_phi(2)*grad_phi(2)*grad2_phi(2) )
                           -invval * ( grad_phi(0)*grad_phi(1)*( grad2_phi(3) + grad2_phi(6) )
                                      + grad_phi(0)*grad_phi(2)*( grad2_phi(4) + grad2_phi(7) )
                                      + grad_phi(1)*grad_phi(2)*( grad2_phi(5) + grad2_phi(8)) )
                           +1.0/grad_phi_norm * ( grad2_phi(0) + grad2_phi(1) + grad2_phi(2) );

            // *-1.0 because of direction of gradient phi vector, which is minus the normal on the interface pointing from + to -
            curvature_gp *= -1.0;
          }
        }

       // check for norm of grad phi almost zero
       if (fabs(grad_phi_norm < 1.0E-9))
       {
         std::cout << "grad phi is small -> set to 1.0E9" << grad_phi_norm << std::endl;
         grad_phi_norm = 1.0;
       }

       // element matrix and rhs
       for (size_t vi=0; vi<numnode; ++vi) // loop rows  (test functions)
       {
         for (size_t ui=0; ui<numnode; ++ui) // loop columns  (test functions)
         {
           const double tmp=fac* (funct_gp(ui)*funct_gp(vi)
                                 + eta_smooth * (deriv_gp_xyz(0,vi) * deriv_gp_xyz(0,ui)
                                               + deriv_gp_xyz(1,vi) * deriv_gp_xyz(1,ui)
                                               + deriv_gp_xyz(2,vi) * deriv_gp_xyz(2,ui)));
           elemat(vi,ui)+=tmp;
         }

           elerhs(vi) += fac*funct_gp(vi)*curvature_gp;
       }
     }//end Gauss loop

     //get global ID of current element
     int eid = actele->Id();
     // assembling
     matrix->Assemble(eid,elemat,lm,lmowner);
     LINALG::Assemble(*rhs,elerhs,lm,lmowner);
  }//end element loop

  // finalize the system matrix
  matrix->Complete();

  // create a solver
  // remark: we take a new here, which is assumed to have number 3
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(3);
  Teuchos::RCP<LINALG::Solver>  solver =
//      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->UMFPACKSolverParams(),
//                                      gfuncdis_->Comm(),
//                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  Teuchos::rcp(new LINALG::Solver(solverparams,
                                  gfuncdis_->Comm(),
                                  DRT::Problem::Instance()->ErrorFile()->Handle()));

  // use only UMFPACKSolver oder ILU PreCond!!!
  // ML requires computation of nullspace
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams,"SOLVER");
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
  if (solvertype != INPAR::SOLVER::umfpack and solvertype != INPAR::SOLVER::aztec_msr)
    dserror("You have to choose UMFPACK or AZPREC_MSR for gradient projection.");
  if (solvertype != INPAR::SOLVER::umfpack and prectyp != INPAR::SOLVER::azprec_ILU)
    dserror("You have to choose AZPREC_MSR with ILU preconditioner for gradient projection.");

  // always refactor and reset the matrix before a single new solver call
  bool refactor=true;
  bool reset   =true;

  // solution vector on dof row map of gfunc discretization
  Teuchos::RCP<Epetra_Vector> solcurve = Teuchos::rcp(new Epetra_Vector(*gfuncdis_->DofRowMap(),true));
  // solver call
  solver->Solve(matrix->EpetraOperator(),
                 solcurve,
                 rhs,
                 refactor,
                 reset);
  solver->Reset();

  // store result in vector
  // loop all nodes on the processor
  for(int lfluidnodeid=0;lfluidnodeid<fluiddis_->NumMyRowNodes();lfluidnodeid++)
  {
    // get the processor's scatra node
    // remark: we rely on identical parallel distributions of fluid and G-function discretizations, here
    DRT::Node* gfuncnode = gfuncdis_->lRowNode(lfluidnodeid);

    const int dofgid = gfuncdis_->Dof(0,gfuncnode,0);

    const int lid = solcurve->Map().LID(dofgid);

    const int nlid = curvaturerow->Map().LID(gfuncnode->Id());

    if (lid<0) dserror("Proc %d: Cannot find gid=%d (lid = %d) in Epetra_Vector",(*curvaturerow).Comm().MyPID(),gfuncnode->Id(),lid);

    double err = curvaturerow->ReplaceMyValues(1,&(*solcurve)[lid],&nlid);
    if (err) dserror("Could not set gradient of phi!");
  }

  // export gradient to column map
  LINALG::Export(*curvaturerow,*curvature_);

#if 0 // only for debug
  // additional gmsh output
  {
    // turn on/off screen output for writing process of Gmsh postprocessing file
    const bool screen_out = true;

    // create Gmsh postprocessing file
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("curvature", 0, 10, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Curvature \" {" << std::endl;
      // draw scalar field 'Phinp' for every element
      IO::GMSH::ScalarFieldNodeBasedToGmsh(fluiddis_,curvaturerow,gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
  }
#endif

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute smoothed gradient of phi based on L2 projection                        rasthofer 08/13 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ComputeL2ProjectedGradPhi(const double eta_smooth,
                                                    const Teuchos::RCP<Epetra_Vector> phi_smoothed)
{
  //get the right number of space dimensions
  const size_t nsd=3;

  // gradphi_ gets a NodeColMap
  gradphi_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),nsd));
  gradphi_->PutScalar(0.0);

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_MultiVector> gradphirow = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeRowMap(),3));
  gradphirow->PutScalar(0.0);

  //get dof row map of scalar field
  const Epetra_Map* dofrowmap = gfuncdis_->DofRowMap();

  // safety check
  if (not fluiddis_->NodeRowMap()->SameAs(*(gfuncdis_->NodeRowMap())))
    dserror("Same node map expected");

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> matrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  // create right hand side
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);

  // element matrix and rhs
  Epetra_SerialDenseMatrix elemat;
  Epetra_SerialDenseVector elerhs;

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // get number of elements
  const int numele = gfuncdis_->NumMyColElements();

  //loop all three dimensions
  // we do all spacial directions separately since we work on the dof row map of the gfunc discretization,
  // which has one dof per node
  // the three directions are decoupled anyway
  for (size_t idim=0;idim<nsd;idim++)
  {
    rhs->PutScalar(0.0);
    matrix->PutScalar(0.0);

    if (!gfuncdis_->Filled())   dserror("FillComplete() was not called");
    if (!gfuncdis_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    //loop column elements
    for (int iele = 0; iele< numele;iele++)
    {
        DRT::Element* actele = gfuncdis_->lColElement(iele);

        // get element location vector, dirichlet flags and ownerships
        lm.clear();
        lmowner.clear();
        lmstride.clear();
        actele->LocationVector(*gfuncdis_,lm,lmowner,lmstride);

        // get shape
        DRT::Element::DiscretizationType distype = actele->Shape();
        if (distype != DRT::Element::hex8) dserror("Hex 8 expected!");

        //set element data: get the right number of nodes per element
        const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

        //node coordinates in global dimensions
        LINALG::Matrix<3,numnode> xyze;

        GEO::fillInitialPositionArray<DRT::Element::hex8>(actele, xyze);

        // get Gauss rule and Gauss points for integration
        const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_hex_8point);

        // Gauss point vector
        Epetra_SerialDenseVector  gp(nsd);

        // get dimension of element matrices and vectors
        // reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();

        if (elemat.M()!=eledim or elemat.N()!=eledim)
          elemat.Shape(eledim,eledim);
        else
          memset(elemat.A(),0,eledim*eledim*sizeof(double));

        if (elerhs.Length()!=eledim)
          elerhs.Size(eledim);
        else
          memset(elerhs.Values(),0,eledim*sizeof(double));

        // loop integration points
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
          // fill vector with Gauss points
          for (size_t rr=0;rr<nsd;++rr)
            gp(rr)=intpoints.qxg[iquad][rr];

          // shape function at Gauss point
          static LINALG::Matrix<numnode,1> funct_gp;
          funct_gp.Clear();
          DRT::UTILS::shape_function_3D(funct_gp,gp(0),gp(1),gp(2),distype);

          // derivatives evaluated at Gauss point
          static LINALG::Matrix<nsd,numnode> deriv_gp;
          deriv_gp.Clear();
          DRT::UTILS::shape_function_3D_deriv1(deriv_gp,gp(0),gp(1),gp(2),distype);

          // Jacobian for integration over element domain
          static LINALG::Matrix<nsd,nsd> xjm_gp;
          xjm_gp.Clear();
          xjm_gp.MultiplyNT(deriv_gp,xyze);

          // invert of Jacobian for integration over element domain
          static LINALG::Matrix<nsd,nsd> xji_gp;
          xji_gp.Clear();
          xji_gp.Invert(xjm_gp);

          // global derivatives
          static LINALG::Matrix<nsd,numnode> deriv_gp_xyz;
          deriv_gp_xyz.Clear();
          deriv_gp_xyz.Multiply(xji_gp,deriv_gp);

          // Jacobian determinant
          const double det = xjm_gp.Determinant();

          // check for degenerated elements
          if (det < 0.0)
            dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

          // integration factor
          const double fac = intpoints.qwgt[iquad]*det;

         // get phi at integration point
         double grad_phi=0.0;
         // element phi vector
         std::vector<double> myephi(numnode);
         if (phi_smoothed == Teuchos::null)
           DRT::UTILS::ExtractMyNodeBasedValues(actele,myephi,*phinp_);
         else
           DRT::UTILS::ExtractMyNodeBasedValues(actele,myephi,*phi_smoothed);
         // get gradient of phi at Gauss point
         for (size_t inode = 0; inode< numnode; inode++)
           grad_phi += myephi[inode] *deriv_gp_xyz(idim,inode);

         // element matrix and rhs
         for (size_t vi=0; vi<numnode; ++vi) // loop rows  (test functions)
         {
           for (size_t ui=0; ui<numnode; ++ui) // loop columns  (test functions)
           {
             const double tmp=fac* (funct_gp(ui)*funct_gp(vi)
                                   + eta_smooth * (deriv_gp_xyz(0,vi) * deriv_gp_xyz(0,ui)
                                                 + deriv_gp_xyz(1,vi) * deriv_gp_xyz(1,ui)
                                                 + deriv_gp_xyz(2,vi) * deriv_gp_xyz(2,ui)));
             elemat(vi,ui)+=tmp;
           }

           elerhs(vi)+=fac*funct_gp(vi)*grad_phi;
         }
       }//end Gauss loop

       //get global ID of current element
       int eid = actele->Id();
       // assembling
       matrix->Assemble(eid,elemat,lm,lmowner);
       LINALG::Assemble(*rhs,elerhs,lm,lmowner);

    }//end element loop

    // finalize the system matrix
    matrix->Complete();

    // create a solver
    // remark: we take a new here, which is assumed to have number 3
    const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(3);
    Teuchos::RCP<LINALG::Solver>  solver_ =
//      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->UMFPACKSolverParams(),
//                                      gfuncdis_->Comm(),
//                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
    Teuchos::rcp(new LINALG::Solver(solverparams,
                                    gfuncdis_->Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle()));
    // use only UMFPACKSolver oder ILU PreCond!!!
    // ML requires computation of nullspace
    const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams,"SOLVER");
    const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
    if (solvertype != INPAR::SOLVER::umfpack and solvertype != INPAR::SOLVER::aztec_msr)
      dserror("You have to choose UMFPACK or AZPREC_MSR for gradient projection.");
    if (solvertype != INPAR::SOLVER::umfpack and prectyp != INPAR::SOLVER::azprec_ILU)
      dserror("You have to choose AZPREC_MSR with ILU preconditioner for gradient projection.");

    // always refactor and reset the matrix before a single new solver call
    bool refactor=true;
    bool reset   =true;

    // solution vector on dof row map of gfunc discretization
    Teuchos::RCP<Epetra_Vector> solgradphi = Teuchos::rcp(new Epetra_Vector(*gfuncdis_->DofRowMap(),true));
    // solver call
    solver_->Solve(matrix->EpetraOperator(),
                   solgradphi,
                   rhs,
                   refactor,
                   reset);
    solver_->Reset();

    // store result in MultiVector
    // loop all nodes on the processor
    for(int lfluidnodeid=0;lfluidnodeid<fluiddis_->NumMyRowNodes();lfluidnodeid++)
    {
      // get the processor's scatra node
      // remark: we rely on identical parallel distributions of fluid and G-function discretizations, here
      DRT::Node* gfuncnode = gfuncdis_->lRowNode(lfluidnodeid);

      const int dofgid = gfuncdis_->Dof(0,gfuncnode,0);

      const int lid = solgradphi->Map().LID(dofgid);

      const int nlid = gradphirow->Map().LID(gfuncnode->Id());

      if (lid<0) dserror("Proc %d: Cannot find gid=%d (lid = %d) in Epetra_Vector",(*gradphirow).Comm().MyPID(),gfuncnode->Id(),lid);

      const int numcol = (*gradphirow).NumVectors();
      if( numcol != (int)nsd) dserror("Number of columns in Epetra_MultiVector is not identically to nsd");

      double err = gradphirow->ReplaceMyValue(nlid,idim,(*solgradphi)[lid]);
      if (err) dserror("Could not set gradient of phi!");
    }

  } // loop dimensions

  // export gradient to column map
  LINALG::Export(*gradphirow,*gradphi_);

#if 0 // only for debug
  // additional gmsh output
  {
    // turn on/off screen output for writing process of Gmsh postprocessing file
    const bool screen_out = true;

    // create Gmsh postprocessing file
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("gradient_phi", 0, 10, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Grad Phi \" {" << std::endl;
      // draw scalar field 'Phinp' for every element
      IO::GMSH::VectorFieldNodeBasedToGmsh(fluiddis_,gradphirow,gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
  }
#endif

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute smoothed second derivatives of phi based on L2 projection              rasthofer 08/13 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ComputeL2ProjectedGrad2Phi(const double eta_smooth)
{
  //get the right number of space dimensions
  const size_t nsd=3;

  // gradphi2_ gets a NodeColMap and with following row vectors(derivatives):xx,yy,zz,xy,xz,yz,yx,zx,zy
  gradphi2_ = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeColMap(),3*nsd));
  gradphi2_->PutScalar(0.0);

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_MultiVector> gradphirow2 = Teuchos::rcp(new Epetra_MultiVector(*fluiddis_->NodeRowMap(),3*nsd));
  gradphirow2->PutScalar(0.0);

  //get dof row map of scalar field
  const Epetra_Map* dofrowmap = gfuncdis_->DofRowMap();

  // safety check
  if (not fluiddis_->NodeRowMap()->SameAs(*(gfuncdis_->NodeRowMap())))
    dserror("Same node map expected");

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> matrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  // create right hand side
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);

  // element matrix and rhs
  Epetra_SerialDenseMatrix elemat;
  Epetra_SerialDenseVector elerhs;

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // get number of elements
  const int numele = gfuncdis_->NumMyColElements();

  //loop all three dimensions
  // we do all special directions separately since we work on the dof row map of the gfunc discretization,
  // which has one dof per node
  // the three directions are decoupled anyway
  for (size_t ijdim=0;ijdim<(3*nsd);ijdim++)
  {
    rhs->PutScalar(0.0);
    matrix->PutScalar(0.0);

    if (!gfuncdis_->Filled())   dserror("FillComplete() was not called");
    if (!gfuncdis_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    //loop column elements
    for (int iele = 0; iele < numele; iele++)
    {
        DRT::Element* actele = gfuncdis_->lColElement(iele);

        // get element location vector, dirichlet flags and ownerships
        lm.clear();
        lmowner.clear();
        lmstride.clear();
        actele->LocationVector(*gfuncdis_,lm,lmowner,lmstride);

        // get shape
        DRT::Element::DiscretizationType distype = actele->Shape();
        if (distype != DRT::Element::hex8) dserror("Hex 8 expected!");

        //set element data: get the right number of nodes per element
        const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

        //node coordinates in global dimensions
        LINALG::Matrix<3,numnode> xyze;

        GEO::fillInitialPositionArray<DRT::Element::hex8>(actele, xyze);

        // get Gauss rule and Gauss points for integration
        const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_hex_8point);

        // Gauss point vector
        Epetra_SerialDenseVector  gp(nsd);

        // get dimension of element matrices and vectors
        // reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();

        if (elemat.M()!=eledim or elemat.N()!=eledim)
          elemat.Shape(eledim,eledim);
        else
          memset(elemat.A(),0,eledim*eledim*sizeof(double));

        if (elerhs.Length()!=eledim)
          elerhs.Size(eledim);
        else
          memset(elerhs.Values(),0,eledim*sizeof(double));

        // loop integration points
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
          // fill vector with Gauss points
          for (size_t rr=0;rr<nsd;++rr)
            gp(rr)=intpoints.qxg[iquad][rr];

          // shape function at Gauss point
          static LINALG::Matrix<numnode,1> funct_gp;
          funct_gp.Clear();
          DRT::UTILS::shape_function_3D(funct_gp,gp(0),gp(1),gp(2),distype);

          // derivatives evaluated at Gauss point
          static LINALG::Matrix<nsd,numnode> deriv_gp;
          deriv_gp.Clear();
          DRT::UTILS::shape_function_3D_deriv1(deriv_gp,gp(0),gp(1),gp(2),distype);

          // Jacobian for integration over element domain
          static LINALG::Matrix<nsd,nsd> xjm_gp;
          xjm_gp.Clear();
          xjm_gp.MultiplyNT(deriv_gp,xyze);

          // invert of Jacobian for integration over element domain
          static LINALG::Matrix<nsd,nsd> xji_gp;
          xji_gp.Clear();
          xji_gp.Invert(xjm_gp);

          // global derivatives
          static LINALG::Matrix<nsd,numnode> deriv_gp_xyz;
          deriv_gp_xyz.Clear();
          deriv_gp_xyz.Multiply(xji_gp,deriv_gp);

          // Jacobian determinant
          const double det = xjm_gp.Determinant();

          // check for degenerated elements
          if (det < 0.0)
            dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

          // integration factor
          const double fac = intpoints.qwgt[iquad]*det;

         double grad2_phi=0.0;
         // projected nodal gradient at this node
         LINALG::Matrix<nsd,numnode> egradphi(true);
         DRT::UTILS::ExtractMyNodeBasedValues(actele, egradphi, gradphi_,3);
         static LINALG::Matrix<nsd*nsd,1> grad_phi2;
         grad_phi2.Clear();
         // second derivatives at Gauss point
         for(size_t inode = 0; inode < numnode; inode++)
         {
             grad_phi2(0) += egradphi(0,inode) *deriv_gp_xyz(0,inode); //xx
             grad_phi2(1) += egradphi(1,inode) *deriv_gp_xyz(1,inode); //yy
             grad_phi2(2) += egradphi(2,inode) *deriv_gp_xyz(2,inode); //zz
             grad_phi2(3) += egradphi(0,inode) *deriv_gp_xyz(1,inode); //xy
             grad_phi2(4) += egradphi(0,inode) *deriv_gp_xyz(2,inode); //xz
             grad_phi2(5) += egradphi(1,inode) *deriv_gp_xyz(2,inode); //yz
             grad_phi2(6) += egradphi(1,inode) *deriv_gp_xyz(0,inode); //yx
             grad_phi2(7) += egradphi(2,inode) *deriv_gp_xyz(0,inode); //zx
             grad_phi2(8) += egradphi(2,inode) *deriv_gp_xyz(1,inode); //zy
         }

         grad2_phi = grad_phi2(ijdim);

         // element matrix and rhs
         for (size_t vi=0; vi<numnode; ++vi) // loop rows  (test functions)
         {
           for (size_t ui=0; ui<numnode; ++ui) // loop columns  (test functions)
           {
             const double tmp=fac* (funct_gp(ui)*funct_gp(vi)
                                   + eta_smooth * (deriv_gp_xyz(0,vi) * deriv_gp_xyz(0,ui)
                                                 + deriv_gp_xyz(1,vi) * deriv_gp_xyz(1,ui)
                                                 + deriv_gp_xyz(2,vi) * deriv_gp_xyz(2,ui)));
             elemat(vi,ui)+=tmp;
           }

           elerhs(vi)+=fac*funct_gp(vi)*grad2_phi;
         }
       }//end Gauss loop

       //get global ID of current element
       int eid = actele->Id();
       // assembling
       matrix->Assemble(eid,elemat,lm,lmowner);
       LINALG::Assemble(*rhs,elerhs,lm,lmowner);
    }//end element loop

    // finalize the system matrix
    matrix->Complete();

    // create a solver
    // remark: we take a new here, which is assumed to have number 3
    const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(3);
    Teuchos::RCP<LINALG::Solver>  solver_ =
//      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->UMFPACKSolverParams(),
//                                      gfuncdis_->Comm(),
//                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
    Teuchos::rcp(new LINALG::Solver(solverparams,
                                    gfuncdis_->Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle()));
    // use only UMFPACKSolver oder ILU PreCond!!!
    // ML requires computation of nullspace
    const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams,"SOLVER");
    const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
    if (solvertype != INPAR::SOLVER::umfpack and solvertype != INPAR::SOLVER::aztec_msr)
      dserror("You have to choose UMFPACK or AZPREC_MSR for gradient projection.");
    if (solvertype != INPAR::SOLVER::umfpack and prectyp != INPAR::SOLVER::azprec_ILU)
      dserror("You have to choose AZPREC_MSR with ILU preconditioner for gradient projection.");

    // always refactor and reset the matrix before a single new solver call
    bool refactor=true;
    bool reset   =true;

    // solution vector on dof row map of gfunc discretization
    Teuchos::RCP<Epetra_Vector> solgrad2phi = Teuchos::rcp(new Epetra_Vector(*gfuncdis_->DofRowMap(),true));
    // solver call
    solver_->Solve(matrix->EpetraOperator(),
                   solgrad2phi,
                   rhs,
                   refactor,
                   reset);
    solver_->Reset();

    // store result in MultiVector
    // loop all nodes on the processor
    for(int lfluidnodeid=0;lfluidnodeid<fluiddis_->NumMyRowNodes();lfluidnodeid++)
    {
      // get the processor's scatra node
      // remark: we rely on identical parallel distributions of fluid and G-function discretizations, here
      DRT::Node* gfuncnode = gfuncdis_->lRowNode(lfluidnodeid);

      const int dofgid = gfuncdis_->Dof(0,gfuncnode,0);

      const int lid = solgrad2phi->Map().LID(dofgid);

      const int nlid = gradphirow2->Map().LID(gfuncnode->Id());

      if (lid<0) dserror("Proc %d: Cannot find gid=%d (lid = %d) in Epetra_Vector",(*gradphirow2).Comm().MyPID(),gfuncnode->Id(),lid);

      const int numcol = (*gradphirow2).NumVectors();
      if( numcol != (3*(int)nsd)) dserror("Number of columns in Epetra_MultiVector is not identically to 2*nsd");

      // commented lines allow for enforcing symmetry of Hesse matrix
//      if (ijdim < (2*nsd))
//      {
        double err = gradphirow2->ReplaceMyValue(nlid,ijdim,(*solgrad2phi)[lid]);
        if (err) dserror("Could not set gradient of phi!");
//      }
//      if (ijdim > (nsd-1) and ijdim < (2*nsd))
//      {
//        int pos = ijdim + 3;
//        double err = gradphirow2->ReplaceMyValue(lid,pos,(*solgrad2phi)[lid]);
//        if (err) dserror("Could not set gradient of phi!");
//      }
    }

  } // loop dimensions

  // export gradient to column map
  LINALG::Export(*gradphirow2,*gradphi2_);

#if 0 // only for debug
  // additional gmsh output
  {
    // turn on/off screen output for writing process of Gmsh postprocessing file
    const bool screen_out = true;

    // create Gmsh postprocessing file
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("gradient2_phi", 0, 10, screen_out, fluiddis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Hesse Phi \" {" << std::endl;
      // draw scalar field 'Phinp' for every element
      IO::GMSH::VectorFieldNodeBasedToGmsh(fluiddis_,gradphirow2,gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
  }
#endif

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute nodal curvature based on L2 projection                                 rasthofer 08/13 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector>  COMBUST::FlameFront::ComputeL2ProjectedPhi(const double eta_smooth)
{

  //get the right number of space dimensions
  const size_t nsd=3;

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_Vector> phiprojectedrow = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap(),true));
  phiprojectedrow->PutScalar(0.0);

  //get dof row map of scalar field
  const Epetra_Map* dofrowmap = gfuncdis_->DofRowMap();

  // safety check
  if (not fluiddis_->NodeRowMap()->SameAs(*(gfuncdis_->NodeRowMap())))
    dserror("Same node map expected");

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> matrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  // create right hand side
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);

  // element matrix and rhs
  Epetra_SerialDenseMatrix elemat;
  Epetra_SerialDenseVector elerhs;

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // get number of elements
  const int numele = gfuncdis_->NumMyColElements();

  rhs->PutScalar(0.0);
  matrix->PutScalar(0.0);

  if (!gfuncdis_->Filled())   dserror("FillComplete() was not called");
  if (!gfuncdis_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  //loop column elements
  for (int i = 0; i< numele;i++)
  {
      DRT::Element* actele = gfuncdis_->lColElement(i);

      // get element location vector, dirichlet flags and ownerships
      lm.clear();
      lmowner.clear();
      lmstride.clear();
      actele->LocationVector(*gfuncdis_,lm,lmowner,lmstride);

      // get shape
      DRT::Element::DiscretizationType distype = actele->Shape();
      if (distype != DRT::Element::hex8) dserror("Hex 8 expected!");

      //set element data: get the right number of nodes per element
      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

      //node coordinates in global dimensions
      LINALG::Matrix<3,numnode> xyze;

      GEO::fillInitialPositionArray<DRT::Element::hex8>(actele, xyze);

      // get Gauss rule and Gauss points for integration
      const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_hex_8point);

      // Gauss point vector
      Epetra_SerialDenseVector  gp(nsd);

      // get dimension of element matrices and vectors
      // reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      if (elemat.M()!=eledim or elemat.N()!=eledim)
        elemat.Shape(eledim,eledim);
      else
        memset(elemat.A(),0,eledim*eledim*sizeof(double));

      if (elerhs.Length()!=eledim)
        elerhs.Size(eledim);
      else
        memset(elerhs.Values(),0,eledim*sizeof(double));

      // loop integration points
      for (int iquad=0; iquad<intpoints.nquad; ++iquad)
      {
        // fill vector with Gauss points
        for (size_t rr=0;rr<nsd;++rr)
          gp(rr)=intpoints.qxg[iquad][rr];

        // shape function at Gauss point
        static LINALG::Matrix<numnode,1> funct_gp;
        funct_gp.Clear();
        DRT::UTILS::shape_function_3D(funct_gp,gp(0),gp(1),gp(2),distype);

        // derivatives evaluated at Gauss point
        static LINALG::Matrix<nsd,numnode> deriv_gp;
        deriv_gp.Clear();
        DRT::UTILS::shape_function_3D_deriv1(deriv_gp,gp(0),gp(1),gp(2),distype);

        // Jacobian for integration over element domain
        static LINALG::Matrix<nsd,nsd> xjm_gp;
        xjm_gp.Clear();
        xjm_gp.MultiplyNT(deriv_gp,xyze);

        // invert of Jacobian for integration over element domain
        static LINALG::Matrix<nsd,nsd> xji_gp;
        xji_gp.Clear();
        xji_gp.Invert(xjm_gp);

        // global derivatives
        static LINALG::Matrix<nsd,numnode> deriv_gp_xyz;
        deriv_gp_xyz.Clear();
        deriv_gp_xyz.Multiply(xji_gp,deriv_gp);

        // Jacobian determinant
        const double det = xjm_gp.Determinant();

        // check for degenerated elements
        if (det < 0.0)
          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

        // integration factor
        const double fac = intpoints.qwgt[iquad]*det;

        // get phi at Gauss point
        double phi_gp = 0.0;
        // element phi vector
        std::vector<double> myephi(numnode);
        DRT::UTILS::ExtractMyNodeBasedValues(actele,myephi,*phinp_);
        for (size_t inode = 0; inode< numnode; inode++)
                  phi_gp += myephi[inode] *funct_gp(inode);

       // element matrix and rhs
       for (size_t vi=0; vi<numnode; ++vi) // loop rows  (test functions)
       {
         for (size_t ui=0; ui<numnode; ++ui) // loop columns  (test functions)
         {
           const double tmp=fac* (funct_gp(ui)*funct_gp(vi)
                                 + eta_smooth * (deriv_gp_xyz(0,vi) * deriv_gp_xyz(0,ui)
                                               + deriv_gp_xyz(1,vi) * deriv_gp_xyz(1,ui)
                                               + deriv_gp_xyz(2,vi) * deriv_gp_xyz(2,ui)));
           elemat(vi,ui)+=tmp;
         }

         elerhs(vi) += fac*funct_gp(vi)*phi_gp;
       }
     }//end Gauss loop

     //get global ID of current element
     int eid = actele->Id();
     // assembling
     matrix->Assemble(eid,elemat,lm,lmowner);
     LINALG::Assemble(*rhs,elerhs,lm,lmowner);
  }//end element loop

  // finalize the system matrix
  matrix->Complete();

  // create a solver
  // remark: we take a new here, which is assumed to have number 3
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(3);
  Teuchos::RCP<LINALG::Solver>  solver_ =
//      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->UMFPACKSolverParams(),
//                                      gfuncdis_->Comm(),
//                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  Teuchos::rcp(new LINALG::Solver(solverparams,
                                  gfuncdis_->Comm(),
                                  DRT::Problem::Instance()->ErrorFile()->Handle()));
  // use only UMFPACKSolver oder ILU PreCond!!!
  // ML requires computation of nullspace
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams,"SOLVER");
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
  if (solvertype != INPAR::SOLVER::umfpack and solvertype != INPAR::SOLVER::aztec_msr)
    dserror("You have to choose UMFPACK or AZPREC_MSR for gradient projection.");
  if (solvertype != INPAR::SOLVER::umfpack and prectyp != INPAR::SOLVER::azprec_ILU)
    dserror("You have to choose AZPREC_MSR with ILU preconditioner for gradient projection.");

  // always refactor and reset the matrix before a single new solver call
  bool refactor=true;
  bool reset   =true;

  // solution vector on dof row map of gfunc discretization
  Teuchos::RCP<Epetra_Vector> solphi = Teuchos::rcp(new Epetra_Vector(*gfuncdis_->DofRowMap(),true));
  // solver call
  solver_->Solve(matrix->EpetraOperator(),
                 solphi,
                 rhs,
                 refactor,
                 reset);
  solver_->Reset();

  // loop all nodes on the processor
  for(int lfluidnodeid=0;lfluidnodeid<fluiddis_->NumMyRowNodes();lfluidnodeid++)
  {
    // get the processor's scatra node
    // remark: we rely on identical parallel distributions of fluid and G-function discretizations, here
    DRT::Node* gfuncnode = gfuncdis_->lRowNode(lfluidnodeid);

    // find out the local dof id of the dof at the scatra node
    const int gfuncdofgid = gfuncdis_->Dof(0,gfuncnode,0);

    const int lgfuncdofidn = solphi->Map().LID(gfuncdofgid);

    if (lgfuncdofidn < 0) dserror("local dof id not found in map for given global dof id");

    // now copy the values
    const double val = (*solphi)[lgfuncdofidn];
    const int err = phiprojectedrow->ReplaceMyValue(lfluidnodeid,0,val);
    if (err) dserror("error while inserting value into phinrow");
  }

  const Teuchos::RCP<Epetra_Vector> phiprojectedcol = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap(),true));
  // export to col map
  LINALG::Export(*phiprojectedrow,*phiprojectedcol);

  return phiprojectedcol;
}


/*------------------------------------------------------------------------------------------------*
 | compute curvature based on G-function                                              henke 05/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ComputeCurvatureForCombustion(const Teuchos::ParameterList& combustdyn)
{
  // gradphi_ lives on fluid NodeColMap
  curvature_ = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap(),true));

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  const Teuchos::RCP<Epetra_Vector> rowcurv = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap(),true));

  // loop over nodes on this processor
  for(int lnodeid=0; lnodeid<fluiddis_->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    DRT::Node*  lnode = fluiddis_->lRowNode(lnodeid);
    const size_t numele = lnode->NumElement();
    // get list of adjacent elements of this node
    DRT::Element** adjeles = lnode->Elements();

    double avcurv = 0.0;

    for (size_t iele=0; iele<numele;iele++)
    {
#ifdef ORACLES
      // skip the part of the domain left of the expansion
      if (lnode->X()[0] < 0.0)
        continue;
#endif
      DRT::Element* adjele = adjeles[iele];

      // number of nodes of this element
      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

      LINALG::Matrix<3,1> posXiDomain(true);
      {
      bool nodefound = false;
      // find out which node in the element is my local node lnode
      for (size_t inode=0; inode<numnode; ++inode)
      {
        if (adjele->NodeIds()[inode] == lnode->Id())
        {
          // get local (element) coordinates of this node
          posXiDomain = DRT::UTILS::getNodeCoordinates(inode,DRT::Element::hex8);
          nodefound = true;
        }
      }
      if (nodefound==false)
        dserror("node was not found in list of elements");
      }

      // smoothed normal vector at this node
      LINALG::Matrix<3,numnode> mygradphi(true);

      // extract local (element level) G-function values from global vector
      DRT::UTILS::ExtractMyNodeBasedValues(adjele, mygradphi, gradphi_,3);

      // get node coordinates of the current element
      LINALG::Matrix<3,numnode> xyze;
      GEO::fillInitialPositionArray<DRT::Element::hex8>(adjele, xyze);

      double curvature=0.0;
      LINALG::Matrix<6,numnode> mygradphi2(true);
      COMBUST::CalcCurvature<DRT::Element::hex8>(curvature,posXiDomain,xyze,mygradphi,mygradphi2);
      //----------------------------------------------------------------
      // cut off curvature at a maximum value to prevent singular values
      //----------------------------------------------------------------
      // calculate largest element diameter
      const double elesize = COMBUST::getEleDiameter<DRT::Element::hex8>(xyze);
      // use 1/h as the maximum admissible curvature
      if (fabs(curvature) > (5.0/elesize) )
      {
        if (curvature > 0.0)
          curvature = -5.0/elesize;
        else
          curvature = 5.0/elesize;

        IO::cout << "curvature cut off at value " << 5.0/elesize << " in element " << adjele->Id() << IO::endl;
      }

      avcurv += curvature;
    }
    avcurv /= numele;
    avcurv *= -1.0;
    const int gid = lnode->Id();
    int nodelid = fluiddis_->NodeRowMap()->LID(gid);
    // insert velocity value into node-based vector
    const int err = rowcurv->ReplaceMyValues(1, &avcurv, &nodelid);
    if (err) dserror("could not insert values for curvature");

  }

  // export NodeRowMap to NodeColMap gradphi_
  LINALG::Export(*rowcurv,*curvature_);

}


/*------------------------------------------------------------------------------------------------*
 | compute curvature based on G-function                                          rasthofer 08/13 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ComputeCurvatureForSurfaceTension(const Teuchos::ParameterList& combustdyn)
{
  // this is only needed in case we use the node based curvature
  if (DRT::INPUT::IntegralValue<INPAR::COMBUST::SurfaceTensionApprox>(combustdyn.sublist("COMBUSTION FLUID"),"SURFTENSAPPROX") == INPAR::COMBUST::surface_tension_approx_nodal_curvature)
  {
    if (DRT::INPUT::IntegralValue<INPAR::COMBUST::NodalCurvatureCalc>(combustdyn.sublist("COMBUSTION FLUID"),"NODAL_CURVATURE") == INPAR::COMBUST::l2_projected)
      ComputeL2ProjectedNodalCurvature(combustdyn);
    else
      ComputeAveragedNodalCurvature(combustdyn);
  }
  else
   curvature_ = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap(),true));


  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute nodal curvature based averages of adjacent elements                     wichmann 02/12 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ComputeAveragedNodalCurvature(const Teuchos::ParameterList& combustdyn)
{
  // gradphi_ lives on fluid NodeColMap
  curvature_ = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap(),true));

  // before we can export to NodeColMap we need reconstruction with a NodeRowMap
  // this is only used in this method, so there is no need to make it an Teuchos::rcp
  Epetra_Vector rowcurv(*fluiddis_->NodeRowMap(),true);

    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid<fluiddis_->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode = fluiddis_->lRowNode(lnodeid);
      const size_t numele = lnode->NumElement();
      // get list of adjacent elements of this node
      DRT::Element** adjeles = lnode->Elements();
      const int gid = lnode->Id();

      bool isPBCMaster = false;
      bool isPBCSlave  = false;

      /********************************************************************
       * Determine if the current node is a master, slave or regular node *
       * regular: average over adjacent elements                          *
       * master:  also consider slave nodes' elements                     *
       * slave:   skip                                                    *
       ********************************************************************/

      Teuchos::RCP<std::map<int, std::vector<int> > > pbccolmap = gfuncdis_->GetAllPBCCoupledColNodes();
      for (std::map<int,std::vector<int> >::const_iterator ipbcmap = pbccolmap->begin(); ipbcmap != pbccolmap->end(); ++ipbcmap)
      {
        if (ipbcmap->first == gid)
        {
          isPBCMaster = true;
          break;
        }
        else
        {
          for (std::vector<int>::const_iterator islave = ipbcmap->second.begin(); islave != ipbcmap->second.end(); ++islave)
          {
            if (*islave == gid)
            {
              isPBCSlave = true;
              break;
            }
          }
        }
      }

      // slave nodes will get their curvature from the master, so we will skip them
      if (isPBCSlave)
        continue;

      /**********************************************************************
       * Determine if any of the adjacent elements (incl. slave's elements) *
       * is cut. If none of them are cut we do not need the curvature       *
       **********************************************************************/

      // the curvature is only needed for elements with boundary cells
      bool iscut = false;
      for (size_t i=0; i < numele; ++i)
      {
        if (InterfaceHandle()->ElementIntersected(adjeles[i]->Id()))
        {
          iscut = true;
          break;
        }
      }
      // now the PBC slave nodes
      if (!iscut and isPBCMaster)
      {
        std::vector<int> slaveids = pbccolmap->find(gid)->second;
        for (std::vector<int>::const_iterator islave = slaveids.begin(); islave != slaveids.end(); ++islave)
        {
          size_t slavenumele = fluiddis_->gNode(*islave)->NumElement();
          DRT::Element** slaveadjeles = fluiddis_->gNode(*islave)->Elements();
          for (size_t i=0; i < slavenumele; ++i)
          {
            if (InterfaceHandle()->ElementIntersected((slaveadjeles[i]->Id())))
            {
              iscut = true;
              break;
            }
          }
          if (iscut)
            break;
        }
      }
      if (not iscut)
        continue;

      /*********************************************************************
       * Calculate the average curvature over all adjacent elements (incl. *
       * slave node's elements)                                            *
       *********************************************************************/

      // do the actual curvature evaluation
      double sumcurv = 0.0;
      int    sumele  = 0;

      for (size_t iele=0; iele<numele;iele++)
      {
        DRT::Element* adjele = adjeles[iele];

        // number of nodes of this element
        const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

        LINALG::Matrix<3,1> posXiDomain(true);
        {
          bool nodefound = false;
          // find out which node in the element is my local node lnode
          for (size_t inode=0; inode<numnode; ++inode)
          {
            if (adjele->NodeIds()[inode] == gid)
            {
              // get local (element) coordinates of this node
              posXiDomain = DRT::UTILS::getNodeCoordinates(inode,DRT::Element::hex8);
              nodefound = true;
              break;
            }
          }
          if (nodefound==false)
            dserror("node was not found in list of elements");
        }

        // smoothed normal vector at this node
        LINALG::Matrix<3,numnode> mygradphi(true);

        // extract local (element level) G-function values from global vector
        DRT::UTILS::ExtractMyNodeBasedValues(adjele, mygradphi, gradphi_,3);

        // smoothed second derivatives at this node
        LINALG::Matrix<9,numnode> mygradphi2(true);
        bool include_2deriv = false;
        if (DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION FLUID"),"L2_PROJECTION_SECOND_DERIVATIVES"))
        {
          // extract local (element level) G-function values from global vector
          DRT::UTILS::ExtractMyNodeBasedValues(adjele, mygradphi2, gradphi2_,9);
          include_2deriv = true;
        }

        // get node coordinates of the current element
        LINALG::Matrix<3,numnode> xyze;
        GEO::fillInitialPositionArray<DRT::Element::hex8>(adjele, xyze);

        double curvature=0.0;
        COMBUST::CalcCurvature<DRT::Element::hex8>(curvature,posXiDomain,xyze,mygradphi,mygradphi2,include_2deriv);

        sumcurv += curvature;
        sumele++;
      }
      // now the PBC slave nodes
      if (isPBCMaster)
      {
        std::vector<int> slaveids = pbccolmap->find(gid)->second;
        for (std::vector<int>::const_iterator islave = slaveids.begin(); islave != slaveids.end(); ++islave)
        {
          size_t slavenumele = fluiddis_->gNode(*islave)->NumElement();
          DRT::Element** slaveadjeles = fluiddis_->gNode(*islave)->Elements();
          for (size_t iele=0; iele < slavenumele; ++iele)
          {
            DRT::Element* adjele = slaveadjeles[iele];

            // number of nodes of this element
            const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

            LINALG::Matrix<3,1> posXiDomain(true);
            {
            bool nodefound = false;
            // find out which node in the element is my local node lnode
            for (size_t inode=0; inode<numnode; ++inode)
            {
              if (adjele->NodeIds()[inode] == *islave)
              {
                // get local (element) coordinates of this node
                posXiDomain = DRT::UTILS::getNodeCoordinates(inode,DRT::Element::hex8);
                nodefound = true;
                break;
              }
            }
            if (nodefound==false)
              dserror("PBC slave node was not found in list of its elements");
            }

            // smoothed normal vector at this node
            LINALG::Matrix<3,numnode> mygradphi(true);

            // extract local (element level) G-function values from global vector
            DRT::UTILS::ExtractMyNodeBasedValues(adjele, mygradphi, gradphi_,3);

            // smoothed second derivatives at this node
            LINALG::Matrix<9,numnode> mygradphi2(true);
            bool include_2deriv = false;
            if (DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION FLUID"),"L2_PROJECTION_SECOND_DERIVATIVES"))
            {
              // extract local (element level) G-function values from global vector
              DRT::UTILS::ExtractMyNodeBasedValues(adjele, mygradphi2, gradphi2_,9);
              include_2deriv = true;
            }

            // get node coordinates of the current element
            LINALG::Matrix<3,numnode> xyze;
            GEO::fillInitialPositionArray<DRT::Element::hex8>(adjele, xyze);

            double curvature=0.0;
            COMBUST::CalcCurvature<DRT::Element::hex8>(curvature,posXiDomain,xyze,mygradphi,mygradphi2,include_2deriv);

            sumcurv += curvature;
            sumele++;
          }
        }
      }

      /*********************************************************************
       * Write average curvature to the epetra vector. If this is a master *
       * node we also write the curvature to all its slave nodes (, which  *
       * we previously skipped)                                            *
       ********************************************************************/

      double avcurv = -sumcurv / sumele;

      if (fabs(avcurv)<1.0E-9) // spurious velocities for almost planar surfaces observed
      {                        // -> set curvature to zero
        avcurv=0.0;
//        IO::cout << "small curvature value < 1.0e-9 set to zero" << IO::endl;
      }

      const int nodelid = rowcurv.Map().LID(gid);
      if (nodelid < 0)
        dserror("could not insert values for curvature");
      // insert velocity value into node-based vector
      rowcurv[nodelid] = avcurv;

      // now the PBC slave nodes
      if (isPBCMaster)
      {
        std::vector<int> slaveids = pbccolmap->find(gid)->second;
        for (std::vector<int>::const_iterator islave = slaveids.begin(); islave != slaveids.end(); ++islave)
        {
          const int slid = rowcurv.Map().LID(*islave);
          if (slid < 0)
            dserror("PBC slave nodes could not be found on the master proc.");
          rowcurv[slid] = avcurv;
        }
      }

    } // loop nodes on this proc

  // export NodeRowMap to NodeColMap gradphi_
  LINALG::Export(rowcurv,*curvature_);

}
