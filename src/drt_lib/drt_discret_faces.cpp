/*---------------------------------------------------------------------*/
/*! \file

\brief a class to manage an enhanced discretization including the faces between elements

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>


#include "drt_discret_faces.H"
#include "../linalg/linalg_mapextractor.H"

#include "drt_exporter.H"

#include "drt_utils.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_intfaces_calc.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_discsh3/discsh3.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/material.H"

#include "../drt_inpar/inpar_xfem.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                    schott 03/12|
 *----------------------------------------------------------------------*/
DRT::DiscretizationFaces::DiscretizationFaces(
    const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : Discretization(name, comm),  // use base class constructor
      extension_filled_(false),
      doboundaryfaces_(false){};

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                          schott 03/12|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationFaces::FillCompleteFaces(bool assigndegreesoffreedom, bool initelements,
    bool doboundaryconditions, bool createinternalfaces)

{
  // call standard FillComlete of base class
  DRT::Discretization::FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions);

  if (createinternalfaces)
  {
    CreateInternalFacesExtension();
  }

  return 0;
}



/*----------------------------------------------------------------------*
 |  Build internal faces extension (public)                 schott 03/12|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationFaces::CreateInternalFacesExtension(const bool verbose)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::DiscretizationFaces::CreateInternalFaces");

  // create internal faces for stabilization along edges
  BuildFaces(verbose);

  // (re)build map of internal faces
  BuildFaceRowMap();
  BuildFaceColMap();

  extension_filled_ = true;

  if (verbose)
  {
    int summyfaces = facerowptr_.size();
    int summall = 0;
    comm_->SumAll(&summyfaces, &summall, 1);

    if (comm_->MyPID() == 0)
      std::cout << "number of created faces:   " << summall << "\n" << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate edge-based integrals (public)               rasthofer 12/12|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationFaces::EvaluateEdgeBased(Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::ParameterList edgebasedparams)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::DiscretizationFaces::EvaluateEdgeBased");


  Teuchos::RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*(this->DofColMap()), true);

  const Epetra_Map* rmap = NULL;
  //  const Epetra_Map* dmap = NULL;

  Teuchos::RCP<Epetra_FECrsMatrix> sysmat_FE;
  if (systemmatrix1 != Teuchos::null)
  {
    rmap = &(systemmatrix1->OperatorRangeMap());
    //    dmap = rmap;
    sysmat_FE = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, *rmap, 256, false));
  }
  else
    dserror("sysmat is NULL!");

  Teuchos::RCP<LINALG::SparseMatrix> sysmat_linalg =
      Teuchos::rcp(new LINALG::SparseMatrix(Teuchos::rcp_static_cast<Epetra_CrsMatrix>(sysmat_FE),
          LINALG::View, true, false, LINALG::SparseMatrix::FE_MATRIX));

  const int numrowintfaces = NumMyRowFaces();

  for (int i = 0; i < numrowintfaces; ++i)
  {
    DRT::Element* actface = lRowFace(i);

    if (actface->ElementType() ==
        DRT::ELEMENTS::DiscSh3LineType::Instance())  // Discrete Structural Shell
    {
      DRT::ELEMENTS::DiscSh3Line* ele = dynamic_cast<DRT::ELEMENTS::DiscSh3Line*>(actface);
      if (ele == NULL) dserror("expect DiscSh3Line element");


      // get the parent Shell elements
      DRT::Element* p_master = ele->ParentMasterElement();
      DRT::Element* p_slave = ele->ParentSlaveElement();


      size_t p_master_numnode = p_master->NumNode();
      size_t p_slave_numnode = p_slave->NumNode();

      std::vector<int> nds_master;
      nds_master.reserve(p_master_numnode);

      std::vector<int> nds_slave;
      nds_slave.reserve(p_slave_numnode);

      for (size_t i = 0; i < p_master_numnode; i++) nds_master.push_back(0);

      for (size_t i = 0; i < p_slave_numnode; i++) nds_slave.push_back(0);


      // Set master ele to the Material for evaluation.
      Teuchos::RCP<MAT::Material> material = p_master->Material();

      // input parameters for structural dynamics
      const Teuchos::ParameterList& params = DRT::Problem::Instance()->StructuralDynamicParams();

      // call the egde-based assemble and evaluate routine
      ele->AssembleInternalFacesUsingNeighborData(
          params, ele, material, nds_master, nds_slave, *this, sysmat_linalg, residual_col);
    }
    else  // Fluid
    {
      DRT::ELEMENTS::FluidIntFace* ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace*>(actface);
      if (ele == NULL) dserror("expect FluidIntFace element");

      // get the parent fluid elements
      DRT::ELEMENTS::Fluid* p_master = ele->ParentMasterElement();
      DRT::ELEMENTS::Fluid* p_slave = ele->ParentSlaveElement();

      size_t p_master_numnode = p_master->NumNode();
      size_t p_slave_numnode = p_slave->NumNode();


      std::vector<int> nds_master;
      nds_master.reserve(p_master_numnode);

      std::vector<int> nds_slave;
      nds_slave.reserve(p_slave_numnode);

      {
        TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");

        for (size_t i = 0; i < p_master_numnode; i++) nds_master.push_back(0);

        for (size_t i = 0; i < p_slave_numnode; i++) nds_slave.push_back(0);
      }

      // call the internal faces stabilization routine for the current side/surface
      TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: AssembleEdgeStabGhostPenalty");

      // set action for elements
      edgebasedparams.set<int>("action", FLD::EOS_and_GhostPenalty_stabilization);

      // Set master ele to the Material for evaluation.
      Teuchos::RCP<MAT::Material> material = p_master->Material();

#ifdef DEBUG
      // Set master ele to the Material for slave.
      Teuchos::RCP<MAT::Material> material_s = p_slave->Material();

      // Test whether the materials for the parent and slave element are the same.
      if (material->MaterialType() != material_s->MaterialType())
        dserror(" not the same material for master and slave parent element");
#endif

      // call the egde-based assemble and evaluate routine
      DRT::ELEMENTS::FluidIntFaceImplInterface::Impl(ele)->AssembleInternalFacesUsingNeighborData(
          ele, material, nds_master, nds_slave, INPAR::XFEM::face_type_std, edgebasedparams, *this,
          sysmat_linalg, residual_col);
    }
  }

  sysmat_linalg->Complete();

  // if the fluid system matrix is of type BlockSparseMatrix, we cannot add
  // and have to split sysmat_linalg - therefore, we try to cast the fluid system matrix!
  // we need RTTI here - the type-IDs are compared and the dynamic cast is only performed,
  // if we really have an underlying BlockSparseMatrix; hopefully that saves some
  // runtime.. (kruse, 09/14)
  if (typeid(*systemmatrix1) == typeid(*sysmat_linalg))
  {
    (systemmatrix1)->Add(*sysmat_linalg, false, 1.0, 1.0);
  }
  else
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> block_sysmat =
        Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix1, false);
    if (block_sysmat == Teuchos::null)
      dserror("Expected fluid system matrix as BlockSparseMatrix. Failed to cast to it.");
    Teuchos::RCP<LINALG::SparseMatrix> f00, f01, f10, f11;
    Teuchos::RCP<Epetra_Map> domainmap_00 =
        Teuchos::rcp(new Epetra_Map(block_sysmat->DomainMap(0)));
    Teuchos::RCP<Epetra_Map> domainmap_11 =
        Teuchos::rcp(new Epetra_Map(block_sysmat->DomainMap(1)));

    // Split sparse system matrix into blocks according to the given maps
    LINALG::SplitMatrix2x2(
        sysmat_linalg, domainmap_00, domainmap_11, domainmap_00, domainmap_11, f00, f01, f10, f11);
    // add the blocks subsequently
    block_sysmat->Matrix(0, 0).Add(*f00, false, 1.0, 1.0);
    block_sysmat->Matrix(0, 1).Add(*f01, false, 1.0, 1.0);
    block_sysmat->Matrix(1, 0).Add(*f10, false, 1.0, 1.0);
    block_sysmat->Matrix(1, 1).Add(*f11, false, 1.0, 1.0);
  }

  //------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector res_tmp(systemvector1->Map(), false);
  Epetra_Export exporter(residual_col->Map(), res_tmp.Map());
  int err2 = res_tmp.Export(*residual_col, exporter, Add);
  if (err2) dserror("Export using exporter returned err=%d", err2);
  systemvector1->Update(1.0, res_tmp, 1.0);

  return;
}


/*----------------------------------------------------------------------*
 |  Build internal faces geometry (public)                  schott 03/12|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationFaces::BuildFaces(const bool verbose)
{
  faces_.clear();

  if (verbose and comm_->MyPID() == 0)
  {
    std::cout << "Create internal faces ..." << std::endl;
  }

  //----------------------------------------------------------------------
  /* First: Create the surface objects between to elements . */

  // map of surfaces in this cloud: (sorted node_ids) -> (surface)
  std::map<std::vector<int>, InternalFacesData> surfmapdata;

  // loop col elements and find all surfaces attached to them
  //
  // REMARK: in a first step: find all surfaces and adjacent elements and fill InternalFacesData
  //         without creating the internal faces elements

  std::vector<DRT::Element*>::iterator fool;

  for (fool = elecolptr_.begin(); fool != elecolptr_.end(); ++fool)
  {
    DRT::Element* ele = *fool;

    //-------------------------------------------
    // create

    DRT::UTILS::BoundaryBuildType buildtype = DRT::UTILS::buildNothing;

    // 3D elements
    if (ele->NumSurface() > 1)  // 2D boundary element and 3D parent element
    {
      buildtype = DRT::UTILS::buildSurfaces;
    }
    else if (ele->NumSurface() == 1)  // 1D boundary element and 2D parent element
    {
      buildtype = DRT::UTILS::buildLines;
    }
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented yet");


    // get node connectivity for specific distype of parent element
    unsigned int nele = 0;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    std::vector<std::vector<int>> connectivity;
    switch (buildtype)
    {
      case DRT::UTILS::buildSurfaces:
      {
        nele = ele->NumSurface();
        connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
        break;
      }
      case DRT::UTILS::buildLines:
      {
        nele = ele->NumLine();
        connectivity = DRT::UTILS::getEleNodeNumberingLines(distype);
        break;
      }
      default:
        dserror("DRT::UTILS::build... not supported");
        break;
    }


    // does DRT::UTILS convention match your implementation of NumSurface() or NumLine()?
    if (nele != connectivity.size()) dserror("number of surfaces or lines does not match!");

    // now, get the nodal information for the new surface/line faces
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = connectivity[iele].size();  // this number changes for pyramids or wedges
      std::vector<int> nodeids(nnode);
      std::vector<DRT::Node*> nodes(nnode);

      // get connectivity info
      for (unsigned int inode = 0; inode < nnode; inode++)
      {
        nodeids[inode] = ele->NodeIds()[connectivity[iele][inode]];
        nodes[inode] = ele->Nodes()[connectivity[iele][inode]];
      }

      // sort the nodes. Used to identify surfaces that are created multiple
      std::sort(nodeids.begin(), nodeids.end());

      // find existing InternalFacesData
      std::map<std::vector<int>, InternalFacesData>::iterator surf_it = surfmapdata.find(nodeids);
      if (surf_it == surfmapdata.end())
      {
        // not found -> generate new Data
        // add the faces information to the map (key is the sorted vector of nodeids)
        surfmapdata.insert(std::pair<std::vector<int>, InternalFacesData>(
            nodeids, InternalFacesData(ele->Id(), nodes, iele)));
      }
      else
      {
        if (surf_it->second.GetSlavePeid() != -1) dserror("slave peid should not be set!!!");
        // if found -> add second neighbor data to existing data
        surf_it->second.SetSlavePeid(ele->Id());
        surf_it->second.SetLSurfaceSlave(iele);

        std::vector<int> localtrafomap;

        // get the face's nodes sorted w.r.t local coordinate system of the parent's face element
        const std::vector<DRT::Node*> nodes_face_master = surf_it->second.GetNodes();
        if (nodes_face_master.size() != nnode)
          dserror(
              "the number of the face w.r.t parent element and slave element are not the same. "
              "That is wrong!");

        // find the nodes given with the master element node numbering also for the slave element
        // to define a connectivity map between the local face's coordinate systems
        for (unsigned int inode = 0; inode < nnode; inode++)  // master face nodes
        {
          int position = -1;
          for (std::size_t knode = 0; knode < nodes.size(); knode++)
          {
            if (nodes[knode] == nodes_face_master[inode]) position = knode;
          }

          if (position >= 0)
            localtrafomap.push_back(position);
          else
            dserror("face's node from master's face element not found in slave's face element!");
        }

        surf_it->second.SetLocalNumberingMap(localtrafomap);
      }
    }  // loop iele

  }  // loop elecolptr_

  //----------------------------------------------------------------------
  // in a second step: create the internal faces elements ( sorted nids -> surface element)
  // REMARK: internal faces are created and distributed on procs as following:
  // * faces are created whenever two adjacent elements are available on this proc (sometimes faces
  // are created multiply on more procs)
  // * each face is created at least once (at least one node of the surface is on a proc a row node
  // and a 1-ring of elements around
  //   this node is available as col elements)
  // * how to set the owner for this face on all procs equally?
  //    -> if one set has been created on a proc, there are both parent elements available as row or
  //    col elements
  //    -> therefore for each node of this surface both parent elements are available
  //    -> choose the node with smallest global id
  //    -> the owner of this node will be the owner for the face
  //       (this criterion is working in the same way on all procs holding this face)

  std::map<std::vector<int>, Teuchos::RCP<DRT::Element>> faces;

  // get pbcs
  Teuchos::RCP<std::map<int, std::vector<int>>> col_pbcmapmastertoslave =
      GetAllPBCCoupledColNodes();

  std::map<std::vector<int>, InternalFacesData>::iterator face_it;
  for (face_it = surfmapdata.begin(); face_it != surfmapdata.end(); ++face_it)
  {
    int master_peid = face_it->second.GetMasterPeid();
    int slave_peid = face_it->second.GetSlavePeid();
    if (master_peid == -1) dserror("Face master expected!");

    dsassert(master_peid == gElement(master_peid)->Id(), "Internal error");
    dsassert(slave_peid == -1 || slave_peid == gElement(slave_peid)->Id(), "Internal error");

    // check for potential periodic boundary conditions and connect respective faces/elements
    std::map<int, std::vector<int>>::iterator masternodes = col_pbcmapmastertoslave->begin();
    if (masternodes != col_pbcmapmastertoslave->end())
    {
      // unconnected face is potential pbc face
      if (slave_peid == -1)
      {
        // get node ids of current face
        std::vector<int> mynodeids = face_it->first;

        // get periodic surface boundary conditions
        // number of pairs of periodic boundary conditions
        int numpbcpairs;
        // vector of periodic surface boundary conditions
        std::vector<DRT::Condition*> mypbcs;
        GetCondition("SurfacePeriodic", mypbcs);
        if (mypbcs.empty())
        {
          GetCondition("LinePeriodic", mypbcs);
        }
        // set number of pairs of periodic boundary conditions
        numpbcpairs = mypbcs.size() / 2;

        // sets of pbc id and related node ids
        // for master and slave
        std::map<int, std::set<int>> mastertopbcset;
        std::map<int, std::set<int>> slavetopbcset;
        for (std::size_t numcond = 0; numcond < mypbcs.size(); numcond++)
        {
          const std::vector<int>* myid =
              mypbcs[numcond]->Get<std::vector<int>>("Id of periodic boundary condition");

          const std::string* mymasterslavetoggle =
              mypbcs[numcond]->Get<std::string>("Is slave periodic boundary condition");

          if (*mymasterslavetoggle == "Master")
          {
            // get global master node ids
            const std::vector<int>* masteridstoadd = mypbcs[numcond]->Nodes();

            // store them in list depending on the pbc id
            for (std::vector<int>::const_iterator idtoadd = (*masteridstoadd).begin();
                 idtoadd != (*masteridstoadd).end(); ++idtoadd)
            {
              (mastertopbcset[myid[0][0]]).insert(*idtoadd);
            }
          }
          else if (*mymasterslavetoggle == "Slave")
          {
            // get global slave node ids
            const std::vector<int>* slaveidstoadd = mypbcs[numcond]->Nodes();

            // store them in list depending on the pbc id
            for (std::vector<int>::const_iterator idtoadd = (*slaveidstoadd).begin();
                 idtoadd != (*slaveidstoadd).end(); ++idtoadd)
            {
              (slavetopbcset[myid[0][0]]).insert(*idtoadd);
            }
          }
          else
            dserror("Unknown type for pbc!");
        }

        // provide vectors for master and slave node ids
        std::vector<int> mymasternodeids;
        std::vector<int> myslavenodeids;
        // provide vector for undefined nodes
        // i.e., special nodes on edges or in corners
        // for multiple pbc sets master nodes of boundary condition
        // become slave nodes
        // e.g. for two sets two master nodes at the corners become slave nodes
        //
        //                PBC M surface
        //           M------------------------S
        //           |                        |
        //  PBC M    |                        | PBC S
        //  surface  |                        | surface
        //           |                        |
        //           S------------------------S
        //                PBC S surface
        // these nodes are not contained in the list col_pbcmapmastertoslave as master nodes
        // but result in more than one slave node for the corner or edge master
        std::vector<int> myfurthermasternodeids;

        // local (or face) master to slave coupling
        std::map<int, int> local_pbcmapmastertoslave;

        // bool to indicate if slave element has been found and should be added to the patch
        bool add_salve_ele_to_face = true;

        // loop node ids of current face and check if they are contained in the list
        // of all master node ids
        for (std::size_t inode = 0; inode < mynodeids.size(); inode++)
        {
          if (col_pbcmapmastertoslave->find(mynodeids[inode]) != col_pbcmapmastertoslave->end())
          {
            // add node id to list of current master nodes
            mymasternodeids.push_back(mynodeids[inode]);
          }
          else
          {
            // if node is not in (master) list col_pbcmapmastertoslave, it may be special
            // node as explained above
            // check whether node is master and slave due to several pbcs
            bool found = false;
            // loop all master sets
            for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
            {
              if ((mastertopbcset[ipbc]).find(mynodeids[inode]) != (mastertopbcset[ipbc]).end())
                found = true;
            }

            // yes, we have a master node here
            // add to list of further master nodes which require special care
            if (found) myfurthermasternodeids.push_back(mynodeids[inode]);
          }
        }

        if (mymasternodeids.size() > 0)
        {
          //          std::cout << "current sets" << std::endl;
          //          std::cout << "master nodes" << std::endl;
          //          for (std::size_t rr=0; rr < mymasternodeids.size(); rr++)
          //            std::cout << mymasternodeids[rr] << std::endl;
          //          std::cout << "further master nodes" << std::endl;
          //          for (std::size_t rr=0; rr < myfurthermasternodeids.size(); rr++)
          //            std::cout << myfurthermasternodeids[rr] << std::endl;

          // check if all nodes of the face are masters of pbcs
          // -> this is a master face
          if ((mymasternodeids.size() + myfurthermasternodeids.size()) == mynodeids.size())
          {
            // get corresponding slave ids
            // do the standard master nodes of col_pbcmapmastertoslave first
            for (std::size_t rr = 0; rr < mymasternodeids.size(); rr++)
            {
              // this master node has one slave node
              if (((*col_pbcmapmastertoslave)[mymasternodeids[rr]]).size() == 1)
              {
                myslavenodeids.push_back(((*col_pbcmapmastertoslave)[mymasternodeids[rr]])[0]);
                local_pbcmapmastertoslave[mymasternodeids[rr]] =
                    ((*col_pbcmapmastertoslave)[mymasternodeids[rr]])[0];
              }
              // this master node has several slave nodes
              // it is a corner or edge node of two or three pbc sets
              else
              {
                // this is only possible for multiple pbcs
                if (numpbcpairs == 1) dserror("Two or three pbs sets expected");

                // identify the pbc condition (i.e., pbc id) to which the current face belongs

                std::map<int, int> pbcspermaster;
                // initialize with zeros
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++) pbcspermaster[ipbc] = 0;

                // identify pbc set to which master nodes belong
                for (std::size_t imnode = 0; imnode < mymasternodeids.size(); imnode++)
                {
                  for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                  {
                    std::set<int>::iterator iter =
                        (mastertopbcset[ipbc]).find(mymasternodeids[imnode]);
                    if (iter != (mastertopbcset[ipbc]).end()) pbcspermaster[ipbc] += 1;
                  }
                }

                // all master nodes of current surface share the same pbc id
                int masterpbcid = -1;
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                {
                  if (pbcspermaster[ipbc] == (int)mymasternodeids.size())
                  {
                    masterpbcid = ipbc;
                    break;
                  }
                }

                // find the corresponding slave of the current master node

                // the corresponding slave node is
                // for 2 pbc sets
                // (i) slave node with respect to the pbc id of the master face
                // (ii) master with respect to the remaining pbc sets
                // for 3 pbc sets
                // here to cases may occur
                // master has 7 slaves -> corner node
                // this results as for 2 sets in
                // (i) slave node with respect to the pbc id of the master face
                // (ii) master with respect to the remaining pbc sets
                // master has 3 slaves -> edge node
                // this results in
                // (i) slave node with respect to the pbc id of the master face
                // (ii) master with respect to one of the two remaining pbc sets
                //  this special case is marked by flag
                bool three_sets_edge_node = false;
                if (numpbcpairs == 3 and
                    ((*col_pbcmapmastertoslave)[mymasternodeids[rr]]).size() == 3)
                  three_sets_edge_node = true;

                // pbc id of master face also for the slave
                int slavepbcid = masterpbcid;
                // identify the remaining pbc sets via their id
                std::vector<int> remainingmasterpbcids;
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                {
                  if (ipbc != slavepbcid) remainingmasterpbcids.push_back(ipbc);
                }

                // loop all slave nodes of the current master and
                // check which node fulfills above conditions
                int actslaveid = -999;
                for (std::size_t islave = 0;
                     islave < ((*col_pbcmapmastertoslave)[mymasternodeids[rr]]).size(); islave++)
                {
                  // get id
                  actslaveid = ((*col_pbcmapmastertoslave)[mymasternodeids[rr]])[islave];

                  // check first criterion -> (i)
                  if ((slavetopbcset[slavepbcid]).find(actslaveid) !=
                      (slavetopbcset[slavepbcid]).end())
                  {
                    std::size_t found = 0;
                    // if satisfied
                    // check second criterion -> (ii)
                    for (std::size_t k = 0; k < remainingmasterpbcids.size(); k++)
                    {
                      if ((mastertopbcset[remainingmasterpbcids[k]]).find(actslaveid) !=
                          (mastertopbcset[remainingmasterpbcids[k]]).end())
                        found++;
                    }

                    if ((not three_sets_edge_node) and found == remainingmasterpbcids.size())
                      break;
                    else if (three_sets_edge_node and found == 1)
                      break;
                  }
                }

                // store in list
                myslavenodeids.push_back(actslaveid);
                local_pbcmapmastertoslave[mymasternodeids[rr]] = actslaveid;
              }
            }

            // next go to the special masters which occur as slaves in the list
            // col_pbcmapmastertoslave and are indeed edge or corner nodes of master surfaces
            if (myfurthermasternodeids.size() > 0)
            {
              // identify the pbc condition (i.e., id) to which the current face belongs
              // perform as explained above of the special master nodes with several slaves

              std::map<int, int> pbcspermaster;
              for (int ipbc = 0; ipbc < numpbcpairs; ipbc++) pbcspermaster[ipbc] = 0;

              // identify pbc set to which master nodes belong
              for (std::size_t imnode = 0; imnode < mymasternodeids.size(); imnode++)
              {
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                {
                  std::set<int>::iterator iter =
                      (mastertopbcset[ipbc]).find(mymasternodeids[imnode]);
                  if (iter != (mastertopbcset[ipbc]).end()) pbcspermaster[ipbc] += 1;
                }
              }
              // identify pbc set to which additional master nodes belong
              for (std::size_t imnode = 0; imnode < myfurthermasternodeids.size(); imnode++)
              {
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                {
                  std::set<int>::iterator iter =
                      (mastertopbcset[ipbc]).find(myfurthermasternodeids[imnode]);
                  if (iter != (mastertopbcset[ipbc]).end()) pbcspermaster[ipbc] += 1;
                }
              }

              // all master nodes of current surface share the same pbc id
              int masterpbcid = -1;
              for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
              {
                if (pbcspermaster[ipbc] ==
                    (int)(mymasternodeids.size() + myfurthermasternodeids.size()))
                {
                  masterpbcid = ipbc;
                  break;
                }
              }

              // find the corresponding slaves of the additional master nodes

              for (std::size_t ifnode = 0; ifnode < myfurthermasternodeids.size(); ifnode++)
              {
                // get node id of master
                int actnodeid = myfurthermasternodeids[ifnode];

                // get list of all potential slave nodes
                // i.e., further edge or corner nodes
                std::vector<int> mypotslaveids;

                // first, look for masters with more than one slave
                // then, check whether node id (i.e. actnodeid) is contained in slave list
                std::map<int, std::vector<int>>::iterator master_it;
                for (master_it = col_pbcmapmastertoslave->begin();
                     master_it != col_pbcmapmastertoslave->end(); master_it++)
                {
                  if ((master_it->second).size() > 1)
                  {
                    bool found = false;

                    for (std::size_t k = 0; k < (master_it->second).size(); k++)
                    {
                      if ((master_it->second)[k] == actnodeid)
                      {
                        found = true;
                      }
                    }

                    if (found)
                    {
                      for (std::size_t k = 0; k < (master_it->second).size(); k++)
                        mypotslaveids.push_back((master_it->second)[k]);
                    }

                    if (found) break;
                  }
                }

                if (mypotslaveids.size() == 0) dserror("Expected to find node!");

                // find the corresponding slave of the current master node

                // the corresponding slave node is
                // for 2 pbc sets
                // (i) slave node with respect to the pbc id of the master face
                // (ii) slave with respect to the remaining pbc sets
                // for 3 pbc sets
                // master has 3 slaves -> edge node (special case 1)
                // this results in
                // (i) slave node with respect to the pbc id of the master face
                // (ii) slave with respect to one of the two remaining pbc sets
                // master has 7 slaves -> corner node (special case 2)
                // (i) slave node with respect to the pbc id of the master face
                // (ii) slave or master with respect to the remaining pbc sets
                //      depending on the status of the current node with respect to those pbcs
                // special case 1 marked by flag
                bool three_sets_edge_node = false;
                if (numpbcpairs == 3 and mypotslaveids.size() == 3) three_sets_edge_node = true;

                int slavepbcid = masterpbcid;

                std::vector<int> remainingslavepbcids;
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                {
                  if (ipbc != slavepbcid) remainingslavepbcids.push_back(ipbc);
                }

                std::vector<int> remainingmasterpbcids;
                for (int ipbc = 0; ipbc < numpbcpairs; ipbc++)
                {
                  if (ipbc != slavepbcid) remainingmasterpbcids.push_back(ipbc);
                }

                // special case 2 marked by flag
                bool corner_node = false;
                int furthermastercond = -1;
                int slavecond_1 = -1;

                // get status of corresponding slave with respect to the
                // remaining to pbc sets
                // may be master and slave or slave and slave
                if (numpbcpairs == 3 and mypotslaveids.size() == 7)
                {
                  corner_node = true;

                  // set the sets to check master or slave
                  for (std::size_t k = 0; k < remainingmasterpbcids.size(); k++)
                  {
                    if ((mastertopbcset[remainingmasterpbcids[k]]).find(actnodeid) !=
                        (mastertopbcset[remainingmasterpbcids[k]]).end())
                      furthermastercond = remainingmasterpbcids[k];
                  }

                  // corresponding slave is pure slave
                  // -> we just have to check the slave lists
                  // -> similar to the 2 pbc-sets case
                  if (furthermastercond == -1)
                  {
                    corner_node = false;
                  }
                  // corresponding slave is master and slave
                  // -> we set the slave list to check
                  else
                  {
                    if (furthermastercond == 0)
                    {
                      if (slavepbcid == 1)
                        slavecond_1 = 2;
                      else if (slavepbcid == 2)
                        slavecond_1 = 1;
                      else
                        dserror("Same pbc ids?");
                    }
                    else if (furthermastercond == 1)
                    {
                      if (slavepbcid == 0)
                        slavecond_1 = 2;
                      else if (slavepbcid == 2)
                        slavecond_1 = 0;
                      else
                        dserror("Same pbc ids?");
                    }
                    else if (furthermastercond == 2)
                    {
                      if (slavepbcid == 0)
                        slavecond_1 = 1;
                      else if (slavepbcid == 1)
                        slavecond_1 = 0;
                      else
                        dserror("Same pbc ids?");
                    }
                    else
                      dserror("Unknown pbc id!");
                  }
                }

                //                std::cout << "furthermastercond " <<furthermastercond<< std::endl;
                //                std::cout << "slavecond_1 " <<slavecond_1<< std::endl;

                //                std::cout<< "master 1 \n" << masterpbcid << std::endl;
                //                std::set<int>::iterator myiter;
                //                for (myiter=mastertopbcset[remainingslavepbcids[0]].begin();
                //                myiter!=mastertopbcset[remainingslavepbcids[0]].end(); myiter++)
                //                {
                //                  std::cout<< *myiter << std::endl;
                //                }
                //                std::cout<< "master 2 \n" << masterpbcid << std::endl;
                //                for (myiter=mastertopbcset[remainingslavepbcids[1]].begin();
                //                myiter!=mastertopbcset[remainingslavepbcids[1]].end(); myiter++)
                //                  std::cout<< *myiter << std::endl;

                // loop all slave nodes of the current master and
                // check which node fulfills above conditions
                int actslaveid = -999;
                for (std::size_t islave = 0; islave < mypotslaveids.size(); islave++)
                {
                  // get slave node id
                  actslaveid = mypotslaveids[islave];

                  // check first criterion
                  if ((slavetopbcset[slavepbcid]).find(actslaveid) !=
                      (slavetopbcset[slavepbcid]).end())
                  {
                    std::size_t found = 0;
                    // if satisfied
                    // check second criterion
                    if (not corner_node)
                    {
                      // check to slave conditions
                      for (std::size_t k = 0; k < remainingslavepbcids.size(); k++)
                      {
                        if ((slavetopbcset[remainingslavepbcids[k]]).find(actslaveid) !=
                            (slavetopbcset[remainingslavepbcids[k]]).end())
                          found++;
                      }
                    }
                    else
                    {
                      // check a master and a slave condition
                      if ((mastertopbcset[furthermastercond]).find(actslaveid) !=
                          (mastertopbcset[furthermastercond]).end())
                        found++;
                      if ((slavetopbcset[slavecond_1]).find(actslaveid) !=
                          (slavetopbcset[slavecond_1]).end())
                        found++;
                    }

                    if ((not three_sets_edge_node) and found == remainingslavepbcids.size())
                      break;
                    else if (three_sets_edge_node and found == 1)
                      break;
                  }
                }

                // store in list
                myslavenodeids.push_back(actslaveid);
                local_pbcmapmastertoslave[myfurthermasternodeids[ifnode]] = actslaveid;
              }
            }

            // sort the slave node ids
            std::sort(myslavenodeids.begin(), myslavenodeids.end());
          }
          else
          {
            add_salve_ele_to_face = false;
          }

          // this criterion ensures that slave set is available on this proc
          int counter = 0;
          for (std::size_t kk = 0; kk < mymasternodeids.size(); kk++)
          {
            std::vector<DRT::Node*>::iterator nofool;
            for (nofool = noderowptr_.begin(); nofool != noderowptr_.end(); ++nofool)
            {
              if ((*nofool)->Id() == mymasternodeids[kk]) counter++;
            }
          }
          for (std::size_t kk = 0; kk < myfurthermasternodeids.size(); kk++)
          {
            std::vector<DRT::Node*>::iterator nofool;
            for (nofool = noderowptr_.begin(); nofool != noderowptr_.end(); ++nofool)
            {
              if ((*nofool)->Id() == myfurthermasternodeids[kk]) counter++;
            }
          }
          if (counter == 0)
          {
            add_salve_ele_to_face = false;
          }

          // add slave element to the patch
          if (add_salve_ele_to_face)
          {
            // get master element
            DRT::Element* master_ele = elecolptr_[0];
            for (fool = elecolptr_.begin(); fool != elecolptr_.end(); ++fool)
            {
              if ((*fool)->Id() == master_peid) master_ele = *fool;
            }

            // look for the corresponding slave face in the list of all faces
            std::map<std::vector<int>, InternalFacesData>::iterator pbc_surf_it =
                surfmapdata.find(myslavenodeids);
            if (pbc_surf_it == surfmapdata.end())
            {
              // print some helpful information first
              master_ele->Print(std::cout);

              std::cout << "\n slave " << std::endl;
              for (std::size_t kk = 0; kk < myslavenodeids.size(); kk++)
                std::cout << myslavenodeids[kk] << std::endl;

              dserror("Expected to find slave face!");
            }

            // add slave data to master data
            face_it->second.SetSlavePeid(pbc_surf_it->second.GetMasterPeid());
            slave_peid = face_it->second.GetSlavePeid();
            face_it->second.SetLSurfaceSlave(pbc_surf_it->second.GetLSurfaceMaster());

            // add connection of coordinate systems for master and slave
            std::vector<int> localtrafomap;

            // get the face's nodes sorted w.r.t local coordinate system of the parent's face
            // element
            const std::vector<DRT::Node*> nodes_face_master = face_it->second.GetNodes();
            // get number of nodes
            unsigned int nnode = nodes_face_master.size();

            // get slave nodes
            std::vector<DRT::Node*> slave_nodes = pbc_surf_it->second.GetNodes();

            // find the nodes given with the master element node numbering also for the slave
            // element to define a connectivity map between the local face's coordinate systems
            for (unsigned int inode = 0; inode < nnode; inode++)  // master face nodes
            {
              int position = -1;

              for (std::size_t knode = 0; knode < slave_nodes.size(); knode++)
              {
                if (slave_nodes[knode]->Id() ==
                    local_pbcmapmastertoslave[nodes_face_master[inode]->Id()])
                  position = knode;
              }

              if (position >= 0)
                localtrafomap.push_back(position);
              else
                dserror(
                    "face's node from master's face element not found in slave's face element!");
            }
            // set in face
            face_it->second.SetLocalNumberingMap(localtrafomap);

            //            if (comm_->MyPID()==1)
            //            {
            //            std::cout << "\n added pbc face "  << std::endl;
            //
            //            std::cout << "master nodes" << std::endl;
            //            for (std::size_t rr=0; rr < mymasternodeids.size(); rr++)
            //              std::cout << mymasternodeids[rr] << std::endl;
            //
            //            std::cout << "further master nodes" << std::endl;
            //            for (std::size_t rr=0; rr < myfurthermasternodeids.size(); rr++)
            //              std::cout << myfurthermasternodeids[rr] << std::endl;
            //
            //              std::cout << "slave node ids "  << std::endl;
            //            for (std::size_t kk=0; kk<myslavenodeids.size(); kk++)
            //               std::cout << myslavenodeids[kk] << std::endl;
            //
            ////            std::cout << "local trafo map  " << std::endl;
            ////            for (std::size_t kk=0; kk<localtrafomap.size(); kk++)
            ////            {
            ////              std::cout << "master node id " << nodes_face_master[kk]->Id() <<
            /// std::endl; /              std::cout << "slave node id " << slave_nodes[kk]->Id() <<
            /// std::endl; /              std::cout << "slave node position " << localtrafomap[kk]
            /// << std::endl; /            }
            //            }
          }
        }
      }
    }

    // create faces
    if (doboundaryfaces_ || (master_peid != -1 && slave_peid != -1))
    {
      dsassert(master_peid != -1, "At least the master element should be present");
      DRT::Element* parent_master = gElement(master_peid);
      DRT::Element* parent_slave = slave_peid != -1 ? gElement(slave_peid) : NULL;

      dsassert(master_peid == parent_master->Id(), "Internal error");
      dsassert(slave_peid == -1 || slave_peid == parent_slave->Id(), "Internal error");

      // get the unsorted nodes
      std::vector<DRT::Node*> nodes = face_it->second.GetNodes();

      // get corresponding nodeids
      std::vector<int> nodeids(nodes.size());
      std::transform(nodes.begin(), nodes.end(), nodeids.begin(), std::mem_fun(&DRT::Node::Id));

      // create the internal face element
      Teuchos::RCP<DRT::FaceElement> surf = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(
          parent_master->CreateFaceElement(parent_slave, nodeids.size(), &nodeids[0], &nodes[0],
              face_it->second.GetLSurfaceMaster(), face_it->second.GetLSurfaceSlave(),
              face_it->second.GetLocalNumberingMap()),
          true);
      dsassert(surf != Teuchos::null,
          "Creating a face element failed. Check overloading of CreateFaceElement");

      // create a clone (the internally created element does not exist anymore when all
      // Teuchos::RCP's finished)
      Teuchos::RCP<DRT::FaceElement> surf_clone =
          Teuchos::rcp(dynamic_cast<DRT::FaceElement*>(surf->Clone()));
      if (surf_clone.get() == NULL) dserror("Invalid element detected. Expected face element");

      // Set owning process of surface to node with smallest gid
      // REMARK: see below
      sort(nodeids.begin(), nodeids.end());
      int owner = gNode(nodeids[0])->Owner();

      // set the owner
      surf_clone->SetOwner(owner);

      // insert the newly created element
      faces.insert(
          std::pair<std::vector<int>, Teuchos::RCP<DRT::Element>>(face_it->first, surf_clone));

      // set face to elements
      parent_master->SetFace(face_it->second.GetLSurfaceMaster(), surf_clone.get());
      if (slave_peid != -1)
        parent_slave->SetFace(face_it->second.GetLSurfaceSlave(), surf_clone.get());
    }
  }

  // Surfaces be added to the faces_-map: (line_id) -> (surface).
  // this clear is important to have here
  // if the discretization has been redistributed (combustion module), we have to
  // rebuild the faces and therefore we have to be sure that the map faces_ is clear
  // therefore, the old faces are deleted and replaced by new ones
  std::map<int, Teuchos::RCP<DRT::Element>> finalFaces;
  AssignGlobalIDs(Comm(), faces, finalFaces);
  for (std::map<int, Teuchos::RCP<DRT::Element>>::iterator faceit = finalFaces.begin();
       faceit != finalFaces.end(); ++faceit)
    faces_[faceit->first] = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(faceit->second, true);

  if (verbose and comm_->MyPID() == 0)
  {
    std::cout << "... done!" << std::endl;
  }

  return;
}  // DRT::DiscretizationFaces::BuildInternalFaces



/*----------------------------------------------------------------------*
 |  Build intfacerowmap_ (private)                          schott 03/12|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationFaces::BuildFaceRowMap()
{
  const int myrank = Comm().MyPID();
  int nummyeles = 0;
  std::map<int, Teuchos::RCP<DRT::FaceElement>>::iterator curr;
  for (curr = faces_.begin(); curr != faces_.end(); ++curr)
    if (curr->second->Owner() == myrank) nummyeles++;
  std::vector<int> eleids(nummyeles);
  facerowptr_.resize(nummyeles);
  int count = 0;
  for (curr = faces_.begin(); curr != faces_.end(); ++curr)
    if (curr->second->Owner() == myrank)
    {
      eleids[count] = curr->second->Id();
      facerowptr_[count] = curr->second.get();
      ++count;
    }
  if (count != nummyeles) dserror("Mismatch in no. of internal faces");
  facerowmap_ = Teuchos::rcp(new Epetra_Map(-1, nummyeles, &eleids[0], 0, Comm()));
  return;
}


/*----------------------------------------------------------------------*
 |  Build intfacecolmap_ (private)                          schott 03/12|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationFaces::BuildFaceColMap()
{
  int nummyeles = (int)faces_.size();
  std::vector<int> eleids(nummyeles);
  facecolptr_.resize(nummyeles);
  std::map<int, Teuchos::RCP<DRT::FaceElement>>::iterator curr;
  int count = 0;
  for (curr = faces_.begin(); curr != faces_.end(); ++curr)
  {
    eleids[count] = curr->second->Id();
    facecolptr_[count] = curr->second.get();
    curr->second->SetLID(count);
    ++count;
  }
  if (count != nummyeles) dserror("Mismatch in no. of elements");
  facecolmap_ = Teuchos::rcp(new Epetra_Map(-1, nummyeles, &eleids[0], 0, Comm()));
  return;
}


/*----------------------------------------------------------------------*
 |  get internal faces row map (public)                     schott 03/12|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DiscretizationFaces::FaceRowMap() const
{
  dsassert(Filled(), "FillComplete() must be called before call to FaceRowMap()");
  return facerowmap_.get();
}


/*----------------------------------------------------------------------*
 |  get internal faces col map (public)                     schott 03/12|
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DiscretizationFaces::FaceColMap() const
{
  dsassert(Filled(), "FillComplete() must be called before call to FaceColMap()");
  return facecolmap_.get();
}


/*----------------------------------------------------------------------*
 |  get global no of internal faces (public)                schott 03/12|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationFaces::NumGlobalFaces() const
{
  dsassert(Filled(), "FillComplete() must be called before call to NumGlobalFaces()");
  return FaceRowMap()->NumGlobalElements();
}


/*----------------------------------------------------------------------*
 |  get no of my row internal faces (public)                schott 03/12|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationFaces::NumMyRowFaces() const
{
  dsassert(Filled(), "FillComplete() must be called before call to NumMyRowFaces()");
  return FaceRowMap()->NumMyElements();
}


/*----------------------------------------------------------------------*
 |  get no of my column internal faces (public)             schott 03/12|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationFaces::NumMyColFaces() const
{
  if (Filled())
    return FaceColMap()->NumMyElements();
  else
    return (int)faces_.size();
}

/*----------------------------------------------------------------------*
 |  query existance of element (public)                kronbichler 05/13|
 *----------------------------------------------------------------------*/
bool DRT::DiscretizationFaces::HaveGlobalFace(int gid) const
{
  std::map<int, Teuchos::RCP<DRT::FaceElement>>::const_iterator curr = faces_.find(gid);
  if (curr == faces_.end())
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------*
 |  get element with global id (public)                kronbichler 05/13|
 *----------------------------------------------------------------------*/
DRT::FaceElement* DRT::DiscretizationFaces::gFace(int gid) const
{
  std::map<int, Teuchos::RCP<DRT::FaceElement>>::const_iterator curr = faces_.find(gid);
#ifdef DEBUG
  if (curr == faces_.end()) dserror("Face with gobal id gid=%d not stored on this proc", gid);
#endif
  return curr->second.get();
}


/*----------------------------------------------------------------------*
 |  << operator                                             schott 03/12|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::DiscretizationFaces& dis)
{
  // print standard discretization info
  dis.Print(os);
  // print additional info about internal faces
  dis.PrintFaces(os);

  return os;
}


/*----------------------------------------------------------------------*
 |  Print internal faces discretization (public)            schott 03/12|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationFaces::PrintFaces(std::ostream& os) const
{
  int numglobalfaces = 0;
  if (Filled())
  {
    numglobalfaces = NumGlobalFaces();
  }
  else
  {
    int nummyfaces = 0;
    std::map<int, Teuchos::RCP<DRT::FaceElement>>::const_iterator ecurr;
    for (ecurr = faces_.begin(); ecurr != faces_.end(); ++ecurr)
      if (ecurr->second->Owner() == Comm().MyPID()) nummyfaces++;

    Comm().SumAll(&nummyfaces, &numglobalfaces, 1);
  }

  // print head
  if (Comm().MyPID() == 0)
  {
    os << "--------------------------------------------------\n";
    os << "Discretization: " << Name() << std::endl;
    os << "--------------------------------------------------\n";
    os << numglobalfaces << " Faces (global)\n";
    os << "--------------------------------------------------\n";
    if (Filled())
      os << "Filled() = true\n";
    else
      os << "Filled() = false\n";
    os << "--------------------------------------------------\n";
  }
  // print elements
  for (int proc = 0; proc < Comm().NumProc(); ++proc)
  {
    if (proc == Comm().MyPID())
    {
      if ((int)faces_.size()) os << "-------------------------- Proc " << proc << " :\n";
      std::map<int, Teuchos::RCP<DRT::FaceElement>>::const_iterator curr;
      for (curr = faces_.begin(); curr != faces_.end(); ++curr)
      {
        os << *(curr->second);
        os << std::endl;
      }
      os << std::endl;
    }
    Comm().Barrier();
  }

  return;
}
