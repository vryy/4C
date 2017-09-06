/*----------------------------------------------------------------------*/
/*!
\file xfem_mesh_projector.cpp

\brief Projection of state vectors between overlapping meshes

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>

*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "xfem_mesh_projector.H"
#include "xfem_discretization_utils.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"

#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_position.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

XFEM::MeshProjector::MeshProjector(
    Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis,
    const Teuchos::ParameterList &          params,
    Teuchos::RCP<const Epetra_Vector>       sourcedisp ) :
    sourcedis_(sourcedis),
    targetdis_(targetdis),
    searchradius_fac_(params.sublist("XFLUID DYNAMIC/GENERAL").get<double>("XFLUIDFLUID_SEARCHRADIUS"))
{
  SetSourcePositionVector(sourcedisp);
  // in case the source discretization is empty on this proc
  if (!sourcedis_->NumMyRowElements())
  {
    searchradius_ = 0.0;
    return;
  }

  // determine the radius of the search tree - grab an arbitrary element to
  // find a characteristic size -dependent length scale of the fluid mesh
  // (not the best choice)
  switch (sourcedis_->lRowElement(0)->Shape())
  {
  case DRT::Element::hex8:
    FindSearchRadius<DRT::Element::hex8>();
    break;
  case DRT::Element::hex20:
    FindSearchRadius<DRT::Element::hex20>();
    break;
  case DRT::Element::hex27:
    FindSearchRadius<DRT::Element::hex27>();
    break;
  default:
    searchradius_ = searchradius_fac_; // avoid a
    break;
  }
}

void XFEM::MeshProjector::SetSourcePositionVector(
  Teuchos::RCP<const Epetra_Vector> sourcedisp
)
{
  src_nodepositions_n_.clear();
  // set position of source nodes
  // we run over the col nodes, as we need the full src_nodepositions
  // for all nodes of an element on each proc
  for (int lid = 0; lid < sourcedis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = sourcedis_->lColNode(lid);
    std::vector<int> src_dofs(4);
    std::vector<double> mydisp(3, 0.0);

    if (sourcedisp != Teuchos::null)
    {
      // get the current displacement
      sourcedis_->Dof(node,0,src_dofs);
      DRT::UTILS::ExtractMyValues(*sourcedisp,mydisp,src_dofs);
    }

    for (int d=0; d < 3; ++d)
      src_nodepositions_n_[node->Id()](d) = node->X()[d]+mydisp.at(d);
  }
}

template<DRT::Element::DiscretizationType distype>
void XFEM::MeshProjector::FindSearchRadius()
{
  DRT::Element* actele = sourcedis_->lRowElement(0);
  const DRT::Node* const* nodes = actele->Nodes();

  // problem dimension
  const unsigned int dim = DRT::UTILS::DisTypeToDim<distype>::dim;

  // we are looking for the maximum diameter of the source element
  // as an estimate for the search radius
  //REMARK: the selection of the embedded element for this estimate is still
  //arbitrary --> choose a sufficiently large safety factor in the input file
  double max_diameter = 0.0;

  // build connectivity matrix for every surface of the embedded element
  std::vector< std::vector<int> > connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(distype);

  //-----------------------------------------------------------------------------------
  // We have hex elements & the faces are quads:
  // the first 4 nodes in the element node numbering vector for a given surface are the
  // corner nodes (equally numbered for hex8/20/hex27), followed by the central nodes
  // in case of hex20/27; in approx. diameter estimation, mid nodes are neglected
  //-----------------------------------------------------------------------------------

  // loop over element surfaces
  for (std::vector< std::vector<int> >::const_iterator ic = connectivity.begin();
       ic != connectivity.end(); ++ ic)
  {
    // get the set of nodes (connected in sequence) for the current surface
    const std::vector<int> & surf_nodeset = *ic;

    // compute the connections 0th->2nd, 1st->3rd corner node
    for (unsigned int icn = 0; icn < 2; ++ icn)
    {
      // next but one node position in vector
      const unsigned icnn = icn+2;

      // compute the distance
      double dist_square = 0.0;
      for (unsigned int isd=0; isd<dim; isd++)
      {
        double dx = nodes[surf_nodeset[icnn]]->X()[isd] - nodes[icn]->X()[isd];
        dist_square += dx*dx;
      }

      double dist = sqrt(dist_square);

      // new maximum?
      if (dist > max_diameter)
        max_diameter = dist;
    }
  } // done with the surface elements

  // the spatial diagonals

  const unsigned ncn_face = 4;
  for (unsigned icn = 0; icn < 1; ++ icn)
  {
    // diagonally opposite (0-6, 1-7)
    {
      const unsigned icn_opp = icn + 2 + ncn_face;
      double dist_square = 0.0;
      for (unsigned int isd=0; isd<dim; isd++)
      {
        double dx = nodes[icn_opp]->X()[isd] - nodes[icn]->X()[isd];
        dist_square += dx*dx;
      }
      double dist = sqrt(dist_square);
      if (dist > max_diameter)
        max_diameter = dist;
    }

    // diagonally opposite (2-4, 3-5)
    {
      const unsigned icn_opp = icn + ncn_face;
      double dist_square = 0.0;
      for (unsigned int isd=0; isd<dim; isd++)
      {
        double dx = nodes[icn_opp]->X()[isd] - nodes[icn+2]->X()[isd];
        dist_square += dx*dx;
      }
      double dist = sqrt(dist_square);
      if (dist > max_diameter)
        max_diameter = dist;
    }
  }

  // TODO: tets are not yet supported by this framework!
  searchradius_ =  searchradius_fac_*max_diameter;
}

void XFEM::MeshProjector::SetupSearchTree()
{
  // init of 3D search tree
  searchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // find the bounding box of all elements of source discretization
  const LINALG::Matrix<3,2> sourceEleBox  = GEO::getXAABBofPositions(src_nodepositions_n_);
  searchTree_->initializeTree(sourceEleBox,*sourcedis_,GEO::TreeType(GEO::OCTTREE));

  // TODO: find the bounding box of the nodes from the target discretization, that demand projection,
  // intersect the bounding boxes to obtain a smaller one
}

void XFEM::MeshProjector::Project(
  std::map<int, std::set<int > >&           projection_nodeToDof,
  std::vector<Teuchos::RCP<Epetra_Vector> > target_statevecs,
  Teuchos::RCP<const Epetra_Vector>         targetdisp
)
{
  //TEUCHOS_FUNC_TIME_MONITOR( "XFEM::MeshProjector::Project" );

  const unsigned num_projection_nodes = projection_nodeToDof.size();
  // size of a fluid dofset
  const unsigned numdofperset = 4;

  targetnodeToParent_.clear();

  // vector of node ids to be projected
  std::vector<int> projection_targetnodes;
  projection_targetnodes.reserve(num_projection_nodes);

  // target node positions (in sequence of projection_targetnodes)
  std::vector<LINALG::Matrix<3,1> > tar_nodepositions_n;
  tar_nodepositions_n.reserve(num_projection_nodes);

  // state vectors veln and accn (in sequence of projection_targetnodes)
  std::vector<LINALG::Matrix<8,1> > interpolated_vecs;
  interpolated_vecs.reserve(num_projection_nodes);

  // set position of nodes in target cloud
  for (std::map<int,std::set<int> >::const_iterator i = projection_nodeToDof.begin();
       i != projection_nodeToDof.end(); ++i)
  {
    const DRT::Node* node = targetdis_->gNode(i->first);

    std::vector<int> tar_dofs(4);
    std::vector<double> mydisp(4, 0.0);

    if (targetdisp != Teuchos::null)
    {
      // get the current displacement
      targetdis_->Dof(node,0,tar_dofs);
      DRT::UTILS::ExtractMyValues(*targetdisp,mydisp,tar_dofs);
    }

    LINALG::Matrix<3,1> pos;
    for (int d=0; d < 3; ++d)
    {
      pos(d) = node->X()[d]+mydisp.at(d);
    }

    tar_dofs.clear();
    mydisp.clear();

    tar_nodepositions_n.push_back(pos);
    projection_targetnodes.push_back(i->first);
    interpolated_vecs.push_back(LINALG::Matrix<8,1>(true));
  }

  SetupSearchTree();

  // vector which identifies if a target node has already interpolated values (initialize to false)
  std::vector<int> have_values(projection_targetnodes.size(),0);
  if (sourcedis_->Comm().NumProc() > 1)
    CommunicateNodes(
        tar_nodepositions_n,
        interpolated_vecs,
        projection_targetnodes,
        have_values);
  else
  {
    FindCoveringElementsAndInterpolateValues(
        tar_nodepositions_n,
        interpolated_vecs,
        projection_targetnodes,
        have_values);
  }

  for (unsigned ni = 0; ni < projection_targetnodes.size(); ++ ni)
  {
    const int node_id = projection_targetnodes[ni];
    const DRT::Node* node = targetdis_->gNode(node_id);

    if (!have_values.at(ni))
    {
#ifdef DEBUG
      if (targetdis_->Comm().MyPID() == 0)
        IO::cout << "WARNING: Found no parent for node: " << node_id << IO::endl;
#endif
      continue;
    }

    const std::set<int>& dofsets = projection_nodeToDof.at(node_id);
    int offset = 0;
    for (size_t iv = 0; iv < target_statevecs.size(); ++iv)
    {
      if (target_statevecs[iv] == Teuchos::null)
        continue;

      std::vector<int> dofs;
      dofs.reserve(dofsets.size()*numdofperset);

      for (std::set<int>::const_iterator iset=dofsets.begin(); iset != dofsets.end(); ++iset)
      {
        targetdis_->Dof(dofs,node,0,*iset);

        for (unsigned isd = 0; isd < numdofperset; ++isd)
        {
          (*target_statevecs[iv])[target_statevecs[iv]->Map().LID(dofs[isd])] = interpolated_vecs[ni](isd+offset);
        }
        dofs.clear();
      }
      offset += numdofperset;
    }
    // if projection was successful, remove the node from the projection map
    projection_nodeToDof.erase(node_id);
  }
}

void XFEM::MeshProjector::ProjectInFullTargetDiscretization(
  std::vector<Teuchos::RCP<Epetra_Vector> > target_statevecs,
  Teuchos::RCP<const Epetra_Vector>         targetdisp
)
{
  //this routine supports only non-XFEM discretizations!
  Teuchos::RCP<const DRT::DiscretizationXFEM> xdiscret = Teuchos::rcp_dynamic_cast<const DRT::DiscretizationXFEM>(targetdis_);
  if (xdiscret != Teuchos::null)
    dserror("Value projection for between different mesh deformation states does not support DiscretizationXFEM.");
  std::map<int,std::set<int> > projection_nodeToDof;
  for (int ni = 0; ni < targetdis_->NumMyRowNodes(); ++ ni)
  {
    const DRT::Node* node = targetdis_->lRowNode(ni);
    // set of dofset indices
    std::set<int> dofsets;
    dofsets.insert(0);

    projection_nodeToDof[node->Id()] = dofsets;
  }

  Project(
      projection_nodeToDof,
      target_statevecs,
      targetdisp);
}

template<DRT::Element::DiscretizationType distype>
bool XFEM::MeshProjector::CheckPositionAndProject(
  const DRT::Element *        src_ele,
  const LINALG::Matrix<3,1> & node_xyz,
  LINALG::Matrix<8,1> &       interpolatedvec
)
{
  // number of element's nodes
  const unsigned int src_numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  // nodal coordinates
  LINALG::Matrix<3,src_numnodes> src_xyze(true);

  for (int in=0; in < src_ele->NumNode(); ++ in)
  {
    const unsigned nid = src_ele->NodeIds()[in];

    for (int d=0; d<3; ++d)
    {
      src_xyze(d,in) = src_nodepositions_n_.at(nid)(d);
    }
  }

  // compute node position w.r.t. embedded element
  Teuchos::RCP<GEO::CUT::Position> pos =
        GEO::CUT::PositionFactory::BuildPosition<3,distype>(src_xyze,node_xyz);
  bool inside = pos->Compute();

  if (inside)
  {
    // node position in covering element's local coordinates
    LINALG::Matrix<3,1> xsi;
    pos->LocalCoordinates(xsi);

    // Evaluate elements shape function at this point and fill values
    LINALG::SerialDenseVector shp(src_numnodes);
    DRT::UTILS::shape_function_3D( shp, xsi(0,0), xsi(1,0), xsi(2,0), distype );

    // extract state values and interpolate
    for (int in=0; in < src_ele->NumNode(); ++ in)
    {
      const DRT::Node* node = src_ele->Nodes()[in];
      const unsigned numdofpernode = src_ele->NumDofPerNode(*node);

      std::vector<double> myval(numdofpernode);
      std::vector<int> src_dofs(numdofpernode);

      sourcedis_->Dof(node,0,src_dofs);
      unsigned offset = 0;
      for (size_t iv = 0; iv < source_statevecs_.size(); ++ iv)
      {
        if (source_statevecs_[iv] == Teuchos::null)
          continue;

        DRT::UTILS::ExtractMyValues(*source_statevecs_[iv],myval,src_dofs);
        for (unsigned isd = 0; isd < numdofpernode; ++isd)
        {
          interpolatedvec(isd+offset) += myval[isd]*shp(in);
        }

        offset += myval.size();

        myval.clear();
      }

      src_dofs.clear();
    }
  }

  return inside;
}

void XFEM::MeshProjector::FindCoveringElementsAndInterpolateValues(
  std::vector<LINALG::Matrix<3,1> >  & tar_nodepositions,
  std::vector<LINALG::Matrix<8,1> >  & interpolated_vecs,
  std::vector<int>                   & projection_targetnodes,
  std::vector<int>                   & have_values
)
{
  // loop over the nodes (coordinates)
  for (unsigned int ni = 0; ni < projection_targetnodes.size(); ++ ni)
  {
    bool insideelement = false;

    // node coordinate
    const LINALG::Matrix<3,1> & node_xyz = tar_nodepositions.at(ni);
    // interpolated vector which is zero at the beginning
    LINALG::Matrix<8,1> interpolatedvec(true);

    //search for near elements
    std::map<int,std::set<int> > closeeles =
        searchTree_->searchElementsInRadius(*sourcedis_,src_nodepositions_n_,node_xyz,searchradius_,0);

    if (closeeles.empty())
    {
      continue;
    }

    // loop over the map of target node-IDs and source elements within the search radius
    for (std::map<int, std::set<int> >::const_iterator closele = closeeles.begin(); closele != closeeles.end(); closele++)
    {
      // loop over the set of source elements within the search radius
      for (std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
      {
        DRT::Element* pele = sourcedis_->gElement(*eleIter);
        // determine values for target fluid node
        switch (pele->Shape())
        {
        case DRT::Element::hex8:
          insideelement = CheckPositionAndProject<DRT::Element::hex8>(
             pele,node_xyz,interpolatedvec);
          break;
        case DRT::Element::hex20:
          insideelement = CheckPositionAndProject<DRT::Element::hex20>(
              pele,node_xyz,interpolatedvec);
          break;
        case DRT::Element::hex27:
          insideelement = CheckPositionAndProject<DRT::Element::hex27>(
              pele,node_xyz,interpolatedvec);
          break;
        default:
          dserror("Unsupported element shape %s!", DRT::DistypeToString(pele->Shape()).c_str()); break;
        }

        if (insideelement)
        {
          targetnodeToParent_[projection_targetnodes[ni]] = pele->Id();
          break;
        }
      }
      if (insideelement)
      {
        if (have_values.at(ni)==0)
        {
          have_values[ni]=1;
          interpolated_vecs.at(ni) = interpolatedvec;
        }
        break;
      }
    }
  }

  return;
}

void XFEM::MeshProjector::CommunicateNodes(
  std::vector<LINALG::Matrix<3,1> >  & tar_nodepositions,
  std::vector<LINALG::Matrix<8,1> >  & interpolated_vecs,
  std::vector<int>                   & projection_targetnodes,
  std::vector<int>                   & have_values
)
{
  // get number of processors and the current processors id
  const int numproc=sourcedis_->Comm().NumProc();

  //information how many processors work at all
  std::vector<int> allproc(numproc);

  // create an exporter for point to point comunication
  DRT::Exporter exporter(sourcedis_->Comm());

  // necessary variables
  MPI_Request request;

  // define send and receive blocks
  std::vector<char> sblock;
  std::vector<char> rblock;

  //----------------------------------------------------------------------
  // communication is done in a round robin loop
  //----------------------------------------------------------------------
  for (int np=0; np<numproc+1; ++np)
  {
    // in the first step, we cannot receive anything
    if (np > 0)
    {
      ReceiveBlock(rblock,exporter,request);

      std::vector<char>::size_type position = 0;
      DRT::ParObject::ExtractfromPack(position,rblock,tar_nodepositions);
      DRT::ParObject::ExtractfromPack(position,rblock,interpolated_vecs);
      DRT::ParObject::ExtractfromPack(position,rblock,projection_targetnodes);
      DRT::ParObject::ExtractfromPack(position,rblock,have_values);
    }

    // in the last step, we keep everything on this proc
    if (np < numproc)
    {
      // -----------------------
      // do what we wanted to do
      FindCoveringElementsAndInterpolateValues(
          tar_nodepositions,
          interpolated_vecs,
          projection_targetnodes,
          have_values
      );

      // Pack info into block to send it
      PackValues(tar_nodepositions,interpolated_vecs,projection_targetnodes,have_values,sblock);

      // add size to sendblock
      SendBlock(sblock,exporter,request);
    }
  } // end of loop over processors
}

void XFEM::MeshProjector::ReceiveBlock( std::vector<char> &   rblock,
                                        DRT::Exporter  &      exporter,
                                        MPI_Request    &      request)
{
  // get number of processors and the current processors id
  int numproc=sourcedis_->Comm().NumProc();
  int myrank =sourcedis_->Comm().MyPID();

  // necessary variables
  int length =-1;
  int frompid=(myrank+numproc-1)%numproc;
  int tag    =frompid;

  // receive from predecessor
  exporter.ReceiveAny(frompid,tag,rblock,length);

#ifdef DEBUG
  //IO::cout << "----receiving " << rblock.size() <<  " bytes: to proc " << myrank << " from proc " << frompid << IO::endl;
#endif

  if (tag!=(myrank+numproc-1)%numproc)
  {
    dserror("received wrong message (ReceiveAny)");
  }

  exporter.Wait(request);

  // for safety
  exporter.Comm().Barrier();

  return;
}

void XFEM::MeshProjector::SendBlock( std::vector<char>  & sblock  ,
                                     DRT::Exporter & exporter,
                                     MPI_Request   & request )
{
  // get number of processors and the current processors id
  int numproc=sourcedis_->Comm().NumProc();
  int myrank =sourcedis_->Comm().MyPID();

  // Send block to next proc.
  int tag    =myrank;
  int frompid=myrank;
  int topid  =(myrank+1)%numproc;

#ifdef DEBUG
   //IO::cout << "----sending " << sblock.size() <<  " bytes: from proc " << myrank << " to proc " << topid << IO::endl;
#endif

  exporter.ISend(frompid,topid,
                 &(sblock[0]),sblock.size(),
                 tag,request);

  // for safety
  exporter.Comm().Barrier();

  return;
}

void XFEM::MeshProjector::PackValues(
  std::vector<LINALG::Matrix<3,1> >  & tar_nodepositions,
  std::vector<LINALG::Matrix<8,1> >  & interpolated_vecs,
  std::vector<int>                   & projection_targetnodes,
  std::vector<int>                   & have_values,
  std::vector<char>                  & sblock)
{
  // Pack info into block to send
  DRT::PackBuffer data;
  DRT::ParObject::AddtoPack(data,tar_nodepositions);
  DRT::ParObject::AddtoPack(data,interpolated_vecs);
  DRT::ParObject::AddtoPack(data,projection_targetnodes);
  DRT::ParObject::AddtoPack(data,have_values);
  data.StartPacking();

  DRT::ParObject::AddtoPack(data,tar_nodepositions);
  DRT::ParObject::AddtoPack(data,interpolated_vecs);
  DRT::ParObject::AddtoPack(data,projection_targetnodes);
  DRT::ParObject::AddtoPack(data,have_values);
  swap( sblock, data() );
}

void XFEM::MeshProjector::GmshOutput(
  int step,
  Teuchos::RCP<const Epetra_Vector> targetdisp
)
{
  // output of source discretization with element numbers and target nodes together with element id
  // of source element for value projection
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("tarnode_to_src_ele", step, 30, 0, targetdis_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    XFEM::UTILS::PrintDiscretizationToStream(Teuchos::rcp_const_cast<DRT::Discretization>(sourcedis_),sourcedis_->Name(),true,false,false,false,false,false,gmshfilecontent,&src_nodepositions_n_);

    gmshfilecontent << "View \" " << "nodeToEle n\" {\n";

    std::vector<int> tar_dofs(3);
    std::vector<double> mydisp(3, 0.0);
    for (int i=0; i<targetdis_->NumMyColNodes(); ++i)
    {
      const DRT::Node* actnode = targetdis_->lColNode(i);
      LINALG::Matrix<3,1> pos(actnode->X(),false);
      if (targetdisp != Teuchos::null)
      {
        // get the current displacement
        targetdis_->Dof(actnode,0,tar_dofs);
        DRT::UTILS::ExtractMyValues(*targetdisp,mydisp,tar_dofs);
        for (unsigned isd = 0; isd < 3; ++ isd)
        {
          pos(isd,0) += mydisp[isd];
        }
        mydisp.clear();
      }
      tar_dofs.clear();

      std::map<int,int>::const_iterator iter = targetnodeToParent_.find(actnode->Id());

      if (iter != targetnodeToParent_.end())
      {
        int id = iter->second;
        IO::GMSH::ScalarToStream(pos,id,gmshfilecontent);
      }
    }
    gmshfilecontent << "};\n";
  }
}
