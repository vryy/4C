/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_serv.cpp

\brief Helper class for coupled simulations (INCA - BACI)

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "aero_tfsi_serv.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_thermo/thermo_element.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_io/io_pstream.H"

#include <Epetra_FEVector.h>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::UTILS::AeroCouplingUtils::AeroCouplingUtils
(
  Teuchos::RCP<DRT::Discretization> structdis,
  Teuchos::RCP<DRT::Discretization> thermodis
) :
istructdis_(0),
ithermodis_(0),
istructdofredumap_(0),
ithermodofredumap_(0),
structreduelements_(0),
structrowmapext_(0),
thermorowmapext_(0),
DinvM_(0),
mastermaxdof_(0),
mastermaxnodeid_(0),
mastermaxeleid_(0)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if(ndim != 3)
    dserror("Coupling with Aero code is only implemented for 3D problems!!");

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  lengthscaling_ = Teuchos::getIntegralValue<double>(tsidyn,"TFSI_length_unit");
  timescaling_ = Teuchos::getIntegralValue<double>(tsidyn,"TFSI_time_unit");

  // declare struct objects in interface
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > > structelements;
  std::map<int, std::map<int, DRT::Node*> > dummy2;  // dummy map
  std::map<int, std::map<int, DRT::Node*> > structgnodes; // col map of structure nodes

  //initialize struct objects in interface
  DRT::UTILS::FindConditionObjects(*structdis, dummy2, structgnodes, structelements,"FSICoupling");

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit;
  std::map<int, Teuchos::RCP<DRT::Element> >::iterator eit;

  int max_id = 0;
  //determine number of coupling surfaces
  for ( meit=structelements.begin(); meit != structelements.end(); meit++ )
  {
    if (meit->first > max_id)
      max_id = meit->first;
  }
  structdis->Comm().MaxAll(&max_id, &maxid_, 1);

  //create a new discretization of the structural/thermal interface
  if (!structdis->Filled())
  {
    structdis->FillComplete();
  }
  if (!thermodis->Filled())
  {
    thermodis->FillComplete();
  }

  // initialize new discretizations
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(structdis->Comm().Clone());
  // structure
  const std::string discret_name1 = "istructdis";
  // thermo
  const std::string discret_name2 = "ithermodis";
  const int myrank = com->MyPID();
  
  for(int interf=0; interf <= maxid_; interf++)
  {
    // initialize new interface discretizations
    Teuchos::RCP<DRT::Discretization> istructnewdis = Teuchos::rcp(new DRT::Discretization(discret_name1,com));
    Teuchos::RCP<DRT::Discretization> ithermonewdis = Teuchos::rcp(new DRT::Discretization(discret_name2,com));

    std::map<int, std::map<int, DRT::Node*> >::iterator mnit;
    std::map<int, DRT::Node* >::iterator nit;

    std::vector<int> nodeids;

    // nodes filled in new interface discretizations
    {
      std::map<int, DRT::Node*> structgnodesinterf = structgnodes[interf];
      for ( nit=structgnodesinterf.begin(); nit != structgnodesinterf.end(); nit++ )
      {
        DRT::Node* currnode = (*nit).second;
        //ids for row map are gathered
        if (currnode->Owner() == myrank)
          nodeids.push_back(currnode->Id());

        //discretizations are feeded
        istructnewdis->AddNode(Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
        ithermonewdis->AddNode(Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
      }
    }
    //row node map of struct interface
    Teuchos::RCP<Epetra_Map> istructnoderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,istructnewdis->Comm()));
    //fully overlapping node map
    Teuchos::RCP<Epetra_Map> istructrednodecolmap = LINALG::AllreduceEMap(*istructnoderowmap);

    //row node map of thermo interface
    Teuchos::RCP<Epetra_Map> ithermonoderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,istructnewdis->Comm()));
    //fully overlapping node map
    Teuchos::RCP<Epetra_Map> ithermorednodecolmap = LINALG::AllreduceEMap(*ithermonoderowmap);


    //care about eles in the interface
    std::vector<int> eleids;
    //eles filled in new interface discretizations
    {
      std::map<int, Teuchos::RCP<DRT::Element> > structelementsinterf = structelements[interf];
      for ( eit=structelementsinterf.begin(); eit != structelementsinterf.end(); eit++ )
      {
        Teuchos::RCP<DRT::Element> currele = eit->second;
        if (currele->Owner() == myrank)
        {
          // each interface stores its row elements in the interface
          eleids.push_back(currele->Id() );

          // structural surface elements cannot be distributed --> Bele3 element is used
          Teuchos::RCP<DRT::Element> istructele = DRT::UTILS::Factory("BELE3","Polynomial", currele->Id(), currele->Owner());
          istructele->SetNodeIds(currele->NumNode(), currele->NodeIds());
          istructnewdis->AddElement( istructele );
          // thermo interface elements
          Teuchos::RCP<DRT::Element> ithermoele = DRT::UTILS::Factory("THERMO","Polynomial", currele->Id(), currele->Owner());
          ithermoele->SetNodeIds(currele->NumNode(), currele->NodeIds());
          Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Thermo>(ithermoele,true)->SetDisType( istructele->Shape() );
          ithermonewdis->AddElement( ithermoele );
        }
      }
    }

    //row ele map of this structural interface
    Teuchos::RCP<Epetra_Map> istructelerowmap = Teuchos::rcp(new Epetra_Map(-1,eleids.size(),&eleids[0],0,istructnewdis->Comm()));
    //fully overlapping ele map
    Teuchos::RCP<Epetra_Map> istructredelecolmap = LINALG::AllreduceEMap(*istructelerowmap);

    //row ele map of this thermo interface
    Teuchos::RCP<Epetra_Map> ithermoelerowmap = Teuchos::rcp(new Epetra_Map(-1,eleids.size(),&eleids[0],0,istructnewdis->Comm()));
    //fully overlapping ele map
    Teuchos::RCP<Epetra_Map> ithermoredelecolmap = LINALG::AllreduceEMap(*ithermoelerowmap);

    //do the fully overlapping ghosting of the structural TFSI interface to have everything redundant
    istructnewdis->ExportColumnNodes(*istructrednodecolmap);
    istructnewdis->ExportColumnElements(*istructredelecolmap);

    //do the fully overlapping ghosting of the thermal TFSI interface to have everything redundant
    ithermonewdis->ExportColumnNodes(*ithermorednodecolmap);
    ithermonewdis->ExportColumnElements(*ithermoredelecolmap);

    //find out if we are parallel; needed for TransparentDofSet
    bool parallel = false;
    if(istructnewdis->Comm().NumProc() != 1)
      parallel = true;

    //dofs of the original discretizations are used to set same dofs for the new interface discretization
    Teuchos::RCP<DRT::DofSet> newdofset1=Teuchos::rcp(new DRT::TransparentDofSet(structdis,parallel));
    istructnewdis->ReplaceDofSet(newdofset1);
    newdofset1=Teuchos::null;

    //final fill complete to reorganize everything in the discretization
    istructnewdis->FillComplete(true, false, false);


    //dofs of the original discretizations are used to set same dofs for the new interface discretization
    Teuchos::RCP<DRT::DofSet> newdofset2=Teuchos::rcp(new DRT::TransparentDofSet(thermodis,parallel));
    ithermonewdis->ReplaceDofSet(newdofset2);
    newdofset2=Teuchos::null;

    //final fill complete to reorganize everything in the discretization
    ithermonewdis->FillComplete(true, false, false);

    //storage of the newly created interface discretization
    istructdis_.push_back(istructnewdis);
    ithermodis_.push_back(ithermonewdis);
    //storage of the full redundant structural elements in each interface
    std::map<int, Teuchos::RCP<DRT::Element> > streles;
    for (int lid=0; lid<istructnewdis->NumMyColElements(); lid++)
    {
      DRT::Element* strele = istructnewdis->lColElement(lid);
      streles[strele->Id()] = Teuchos::rcp(strele, false);
    }
    structreduelements_.push_back(streles);

    //storage of useful dof maps
    istructdofredumap_.push_back(Teuchos::rcp(istructnewdis->DofColMap(), false));
    ithermodofredumap_.push_back(Teuchos::rcp(ithermonewdis->DofColMap(), false));

    // map extractor for the communication between the full thermo field and the newly created thermo surf discret
    structrowmapext_.push_back(Teuchos::rcp(new LINALG::MapExtractor(*(structdis->DofRowMap()),Teuchos::rcp(istructnewdis->DofRowMap(), false))));

    // map extractor for the communication between the full thermo field and the newly created thermo surf discret
    thermorowmapext_.push_back(Teuchos::rcp(new LINALG::MapExtractor(*(thermodis->DofRowMap()),Teuchos::rcp(ithermonewdis->DofRowMap(), false))));

    // importer in order to contract result values on proc0 for communication with INCA (proc0 gets everything, other procs empty)
    Teuchos::RCP<Epetra_Map> proc0datamap = LINALG::AllreduceEMap(*istructnewdis->DofRowMap(),0);
    istructproc0importer_.push_back(Teuchos::rcp(new Epetra_Import (*proc0datamap,*istructdis_[interf]->DofRowMap())));
    proc0datamap = LINALG::AllreduceEMap(*ithermonewdis->DofRowMap(),0);
    ithermoproc0importer_.push_back(Teuchos::rcp(new Epetra_Import (*proc0datamap,*ithermodis_[interf]->DofRowMap())));

    // find largest dof/node/element id in TSI problem for offset in mortar
    mastermaxdof_ = std::max(mastermaxdof_, thermodis->DofRowMap()->MaxAllGID() + 1);
    mastermaxnodeid_ = std::max(mastermaxnodeid_, thermodis->NodeRowMap()->MaxAllGID() + 1);
    mastermaxeleid_ = std::max(mastermaxeleid_, istructdis_[interf]->ElementRowMap()->MaxAllGID() + 1);

    // build mortar interface and fill it initially with fully redundant master side (=structural side)
    const Epetra_Comm& lcomm = istructdis_[interf]->Comm();

    // transfer of four degrees of freedom per node
    const int dim = 3;
    int dofpernode = 4;
    Teuchos::ParameterList input(DRT::Problem::Instance()->MortarCouplingParams());
    // set options for mortar coupling
    input.set<std::string>("SEARCH_ALGORITHM","Binarytree");
    input.set<double>("SEARCH_PARAM", 0.1);
    input.set<std::string>("SEARCH_USE_AUX_POS", "no");
    input.set<std::string>("PARALLEL_REDIST","no");
    input.set<std::string>("SHAPEFCN","dual");
    // master side is given to the mortar interface once in the beginning with full overlap
    input.set<std::string>("REDUNDANT_STORAGE","Master");
    // for triangles which will be almost all lying in a quad element on the master side, this should be a good set of params
    input.set<std::string>("INTTYPE","Elements");
    input.set<int>("NUMGP_PER_DIM",3);
    input.set<bool>("NURBS",false);

    // create mortar interface for each coupling interface with INCA
    Teuchos::RCP<MORTAR::MortarInterface> interface = Teuchos::rcp(
        new MORTAR::MortarInterface(interf, lcomm, dim, input, INPAR::MORTAR::redundant_master));

    // feeding master nodes to the interface including ghosted nodes
    for(int lid=0; lid<istructdis_[interf]->NumMyColNodes(); lid++)
    {
      DRT::Node* node = istructdis_[interf]->lColNode(lid);
      std::vector<int> dofids(dofpernode);
      // get the dofs of the structural node
      for (int k=0; k<dim; ++k)
       dofids[k] = istructdis_[interf]->Dof(node)[k];
      // fill in remaining thermal dof
      dofids[dim] = ithermodis_[interf]->Dof(ithermodis_[interf]->lColNode(lid))[0];

      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(new MORTAR::MortarNode(node->Id(),
         node->X(), node->Owner(), dofpernode, dofids, false));

      interface->AddMortarNode(mrtrnode);
    }

    // feeding master elements to the interface
    for(int lid=0; lid<istructdis_[interf]->NumMyColElements(); ++lid)
    {
     DRT::Element* currele = istructdis_[interf]->lColElement(lid);
     Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(new MORTAR::MortarElement(currele->Id(),
         currele->Owner(), currele->Shape(), currele->NumNode(), currele->NodeIds(), false));

     interface->AddMortarElement(mrtrele);
    }

    mortarinterface_.push_back(interface);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::BuildMortarCoupling
(
  int interf,
  Teuchos::RCP<Epetra_Vector> idispn,
  std::vector<LINALG::Matrix<3,1> >& aerocoords,
  const int nodeoffset
)
{
  // new mortar coupling matrices are needed in every time step
  if(interf == 0)
    DinvM_.clear();

  // remove all slave nodes and elements for new coupling
  mortarinterface_[interf]->RemoveSingleInterfaceSide(true);

  int myrank = mortarinterface_[interf]->Comm().MyPID();

  // transfer of four degrees of freedom per node (3 str + 1 thr) in a 3D problem
  const int dim = 3;
  const int dofpernode = 4;

  // feeding slave nodes to the interface
  // Note: offset is needed for slave dofs and slave node ids
  for (int nodeiter = 0; nodeiter<(int)aerocoords.size(); ++nodeiter)
  {
    double nodalpos[dim];
    for (int k=0; k<dim; ++k)
      nodalpos[k] = aerocoords[nodeiter](k);
    std::vector<int> dofids(dofpernode);
    for (int k=0; k<dofpernode; ++k)
      dofids[k] = mastermaxdof_ + (nodeiter+nodeoffset)*dofpernode + k;
    Teuchos::RCP<MORTAR::MortarNode> mrtrnode =
        Teuchos::rcp(new MORTAR::MortarNode(nodeiter+nodeoffset+mastermaxnodeid_, nodalpos, myrank, dofpernode, dofids, true));

    mortarinterface_[interf]->AddMortarNode(mrtrnode);
  }

  // feeding slave elements (tris) to the interface
  int numtris = aerocoords.size()/3;
  for (int elemiter = 0; elemiter<numtris; ++elemiter)
  {
    std::vector<int> nodeids(3);
    for(int k=0; k<3; ++k)
    {
      int nodeiter = elemiter*3 + k;
      nodeids[k] = nodeiter + nodeoffset + mastermaxnodeid_;
    }

    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(new MORTAR::MortarElement(elemiter + nodeoffset/3 + mastermaxeleid_,
        myrank, DRT::Element::tri3, (int)nodeids.size(), &nodeids[0], true));

    mortarinterface_[interf]->AddMortarElement(mrtrele);
  }

  // finalize the mortar interface without new ghosting
  mortarinterface_[interf]->FillComplete(0, false);

  // print parallel distribution
//  mortarinterface_[interf]->PrintParallelDistribution(interf);

  // we are finally interested in D and M
  Teuchos::RCP<Epetra_Map> slavedofrowmap  = Teuchos::rcp(new Epetra_Map(*mortarinterface_[interf]->SlaveRowDofs()));
  Teuchos::RCP<Epetra_Map> masterdofrowmap = Teuchos::rcp(new Epetra_Map(*mortarinterface_[interf]->MasterRowDofs()));

  // interface displacement has to be merged from slave and master side
  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(masterdofrowmap,slavedofrowmap, false);
  Teuchos::RCP<Epetra_Vector> mergeddispn = LINALG::CreateVector(*dofrowmap, true);

  // insert current displacement status of master side; slave side is placed correctly
  for(int i=0; i<idispn->MyLength(); ++i)
  {
    int gid = idispn->Map().GID(i);
    int lid = mergeddispn->Map().LID(gid);
    (*mergeddispn)[lid] = (*idispn)[i];
  }

  // set displacement state in mortar interface
  mortarinterface_[interf]->SetState("displacement", mergeddispn);

  // create binary search tree
  mortarinterface_[interf]->CreateSearchTree();

  // in the following two steps MORTAR does all the work
  mortarinterface_[interf]->Initialize();
  mortarinterface_[interf]->Evaluate();

  // preparation for AssembleDM
  Teuchos::RCP<LINALG::SparseMatrix> D = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap, 10));
  Teuchos::RCP<LINALG::SparseMatrix> M = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap, 100));
  mortarinterface_[interf]->AssembleDM(*D, *M);

  // Complete() global Mortar matrices
  D->Complete();
  M->Complete(*masterdofrowmap, *slavedofrowmap);

  // extract diagonal of D into diag
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*slavedofrowmap,true);
  D->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
  if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  diag->Reciprocal(*diag);

  // left scale instead of matrix-matrix multiplication Dinv * M
  int err = M->LeftScale(*diag);
  if(err)
    dserror("LeftScale of M matrix failed");

  // store DinvM
  DinvM_.push_back(M);

  // visualize interface with gmsh
//  if(myrank==0)
//    cout << "mortar interface is printed in gmsh file format" << endl;
//  mortarinterface_[interf]->VisualizeGmsh(0,0);

  return;
}


/*----------------------------------------------------------------------*
 | transfer load vector from fluid to structure             ghamm 12/13 |
 *----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::TransferFluidLoadsToStruct
(
  int interf,
  std::vector<LINALG::Matrix<4,1> >& aeroforces,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  Teuchos::RCP<Epetra_Vector> iload_fl = LINALG::CreateVector(DinvM_[interf]->RowMap(),true);

  // NOTE: There is no dofset available for artificially created fluid nodes
  // Hence: Order is important and values from aeroforces can directly be inserted into interface load vector
  for(size_t nodeiter=0; nodeiter<aeroforces.size(); ++nodeiter)
  {
    LINALG::Matrix<4,1> currload = aeroforces[nodeiter];
    for(size_t k=0; k<4; ++k)
    {
      (*iload_fl)[nodeiter*4 + k] = currload(k);
    }
  }

  // compute corresponding structural values using mortar
  // return type is Epetra_MultiVector because Epetra_FEVector is used in SlaveToMaster
  Teuchos::RCP<Epetra_MultiVector> iload_solid = SlaveToMaster(interf, iload_fl);

  // split data in iload_solid into structural and thermal load (4 values per node)
  int numnodes = iload_solid->MyLength()*0.25;
  for(int nodeiter=0; nodeiter<numnodes; ++nodeiter)
  {
    // structural load
    for(int k=0; k<3; ++k)
    {
      int lid = nodeiter*4 + k;
      double val = (*(*iload_solid)(0))[lid];
      int gid = iload_solid->Map().GID(lid);
      (*iforce)[ iforce->Map().LID(gid) ] = val;
    }
    // thermal load
    int lid = nodeiter*4 + 3;
    double val = (*(*iload_solid)(0))[lid];
    int gid = iload_solid->Map().GID(lid);
    (*ithermoload)[ ithermoload->Map().LID(gid) ] = val;
  }

  return;
}


/*----------------------------------------------------------------------*
 | slave to master transfer helper                          ghamm 12/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> FS3I::UTILS::AeroCouplingUtils::SlaveToMaster
(
  int interf,
  Teuchos::RCP<Epetra_Vector> sv
) const
{
  if ( not (DinvM_[interf]->RowMap().SameAs( sv->Map())) )
    dserror("Vector with slave dof map expected");

  // Matrix vector product -(D^{-1}*M)^{T}*sv needs to be computed manually because col and domain map of DinvM have
  // equal entries but may be shuffled.
  // Internal importers and exporters are not built when calling DinvM_->Multiply(..) because col and domain map
  // are tested via SameAs which gives the correct answer regarding the entries but not the order of the entries.
  // Background: col and domain are equal because slave elements are handled in a DG like fashion meaning that there
  // is no node that is hold by several elements. We have distinct nodes at the same position for neighboring elements.
  // Hence there is no ghosting needed and col and domain map have equal entries

#ifdef DEBUG
  if ( not (DinvM_[interf]->DomainMap().PointSameAs( *mortarinterface_[interf]->MasterRowDofs())) )
    dserror("Vector with master dofs does not match domain map of coupling matrix");
#endif

  // FE_Vector for assembling -(D^{-1}*M)^{T}*sv
  Teuchos::RCP<Epetra_FEVector> mv = Teuchos::rcp(new Epetra_FEVector(*mortarinterface_[interf]->MasterRowDofs()));

  int DinvM_length = DinvM_[interf]->MaxNumEntries();
  int DinvM_numentries = 0;
  double DinvM_values[DinvM_length];
  int DinvM_indices[DinvM_length];
  double results[DinvM_length];
  for(int irow=0; irow<DinvM_[interf]->EpetraMatrix()->NumMyRows(); ++irow)
  {
    int rowgid = DinvM_[interf]->RowMap().GID(irow);
    DinvM_[interf]->EpetraMatrix()->ExtractGlobalRowCopy(rowgid, DinvM_length, DinvM_numentries, &DinvM_values[0], &DinvM_indices[0]);

    // transposed of DinvM_ is used
    for(int icol=0; icol<DinvM_numentries; ++icol)
    {
      results[icol] = - DinvM_values[icol] * (*sv)[irow];
    }

    int err = mv->SumIntoGlobalValues(DinvM_numentries, DinvM_indices, results);
    if (err<0)
      dserror("summing into Epetra_FEVector failed");
  }

  // call global assemble
  int err = mv->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble of mv failed");

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::TransferStructValuesToFluid
(
  int interf,
  Teuchos::RCP<const Epetra_Vector> idisp,
  Teuchos::RCP<const Epetra_Vector> ivel,
  Teuchos::RCP<const Epetra_Vector> itemp,
  std::vector<double>& packeddata,
  bool writedata
)
{
  // make all interface data available on proc 0
  // idisp and ivel
  Teuchos::RCP<Epetra_Vector> proc0disp = Teuchos::rcp(new Epetra_Vector(istructproc0importer_[interf]->TargetMap()));
  int err = proc0disp->Import(*idisp,*istructproc0importer_[interf],Insert);
  if (err>0)
    dserror("Importing everything to proc 0 went wrong. Import returns %d",err);
  Teuchos::RCP<Epetra_Vector> proc0vel = Teuchos::rcp(new Epetra_Vector(istructproc0importer_[interf]->TargetMap()));
  err = proc0vel->Import(*ivel,*istructproc0importer_[interf],Insert);
  if (err>0)
    dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  // itemp
  Teuchos::RCP<Epetra_Vector> proc0temp = Teuchos::rcp(new Epetra_Vector(ithermoproc0importer_[interf]->TargetMap()));
  err = proc0temp->Import(*itemp,*ithermoproc0importer_[interf],Insert);
  if (err>0)
    dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  // gather data on TSI side on proc 0
  if(istructdis_[interf]->Comm().MyPID() == 0)
  {
    if(istructdis_[interf]->lColElement(0)->NumNode() != 4)
      dserror("current implementation only for quad4 elements in the coupling interface");

    int numcolele = istructdis_[interf]->NumMyColElements();
    packeddata.resize(numcolele*35);
    for(int lid=0; lid<istructdis_[interf]->NumMyColElements(); ++lid)
    {
      DRT::Element* struele = istructdis_[interf]->lColElement(lid);
      // get element location vector of structural problem
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      struele->LocationVector(*istructdis_[interf],lm,lmowner,lmstride);

      // get element data
      // disp
      std::vector<double> eledisp(lm.size());
      DRT::UTILS::ExtractMyValues(*proc0disp, eledisp ,lm);
      // vel
      std::vector<double> elevel(lm.size());
      DRT::UTILS::ExtractMyValues(*proc0vel, elevel ,lm);

      // get element location vector of thermal problem
      DRT::Element* threle = ithermodis_[interf]->lColElement(lid);
      threle->LocationVector(*ithermodis_[interf],lm,lmowner,lmstride);

      // get element data
      // temp
      std::vector<double> eletemp(lm.size());
      DRT::UTILS::ExtractMyValues(*proc0temp, eletemp ,lm);

      // fill in position, velocity and temperature for quad4 (fifth point is center point)
      // x1y1z1 x2y2z2 x3y3z3 x4y4z4 x5y5z5  u1v1w1 u2v2w2 u3v3w3 u4v4w4 u5v5w5  T1T2T3T4T5
      int eleoffset = lid*35;
      std::vector<double> xcenter(3, 0.0);
      std::vector<double> velcenter(3, 0.0);
      double Tcenter = 0.0;
      // loop over nodes of element
      for(int k=0; k<4; ++k)
      {
        DRT::Node** nodes = struele->Nodes();

        for(int d=0; d<3; ++d)
        {
          // pos
          double pos = (nodes[k]->X()[d] + eledisp[k*3 + d]) / LengthScaling();
          packeddata[eleoffset + k*3 + d] = pos;
          xcenter[d] += pos/4.0;
          // vel
          packeddata[eleoffset + (k+5)*3 + d] = (elevel[k*3 + d]) / LengthScaling() * TimeScaling();
          velcenter[d] += elevel[k*3 + d]/4.0;
        }
        // temp
        packeddata[eleoffset + 10*3 + k] = eletemp[k];
        Tcenter += eletemp[k]/4.0;
      }
      // center point is added
      for(int d=0; d<3; ++d)
      {
        // pos
        packeddata[eleoffset + 4*3 + d] = xcenter[d];
        // vel
        packeddata[eleoffset + 9*3 + d] = velcenter[d];
      }
      // temp
      packeddata[eleoffset + 10*3 + 4] = Tcenter;
    }

    // output
//    for(int i=0; i<(int)packeddata.size(); i++)
//    {
//      if(i%35 == 0)
//        std::cout << "\n element " << i/35 << " is packed: ";
//      std::cout << packeddata[i] << "  " ;
//    }

//    // node 305 is observed
//    DRT::Node* node305 = istructdis_[interf]->gNode(305);
//    std::vector<int> lm = istructdis_[interf]->Dof(node305);
//    std::vector<double> nodedisp(lm.size());
//    DRT::UTILS::ExtractMyValues(*proc0disp, nodedisp ,lm);
//
//    if(writedata == true)
//    {
//      FILE *outFile;
//      outFile = fopen("interfaceDisp.txt", "a");
//      fprintf(outFile, "%.8e  %.8e  %.8e\n", nodedisp[0]/LengthScaling(), nodedisp[1]/LengthScaling(), nodedisp[2]/LengthScaling());
//      fclose(outFile);
//    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::TransferStructValuesToFluidConforming
(
  int interf,
  Teuchos::RCP<const Epetra_Vector> idisp,
  Teuchos::RCP<const Epetra_Vector> itemp,
  std::vector<double>& packeddata,
  bool writedata
)
{
  // make all interface data available on proc 0
  // idisp
  Teuchos::RCP<Epetra_Vector> proc0disp = Teuchos::rcp(new Epetra_Vector(istructproc0importer_[interf]->TargetMap()));
  int err = proc0disp->Import(*idisp,*istructproc0importer_[interf],Insert);
  if (err>0)
    dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  // itemp
  Teuchos::RCP<Epetra_Vector> proc0temp = Teuchos::rcp(new Epetra_Vector(ithermoproc0importer_[interf]->TargetMap()));
  err = proc0temp->Import(*itemp,*ithermoproc0importer_[interf],Insert);
  if (err>0)
    dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  // gather data on TSI side on proc 0
  if(istructdis_[interf]->Comm().MyPID() == 0)
  {
    if(istructdis_[interf]->lColElement(0)->NumNode() != 4)
      dserror("current implementation only for quad4 elements in the coupling interface");

    int numcolele = istructdis_[interf]->NumMyColElements();
    packeddata.resize(numcolele*4);
    for(int lid=0; lid<istructdis_[interf]->NumMyColElements(); ++lid)
    {
      DRT::Element* struele = istructdis_[interf]->lColElement(lid);
      // get element location vector of structural problem
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      struele->LocationVector(*istructdis_[interf],lm,lmowner,lmstride);

      // get element data
      // disp
      std::vector<double> eledisp(lm.size());
      DRT::UTILS::ExtractMyValues(*proc0disp, eledisp ,lm);

      // get element location vector of thermal problem
      DRT::Element* threle = ithermodis_[interf]->lColElement(lid);
      threle->LocationVector(*ithermodis_[interf],lm,lmowner,lmstride);

      // get element data
      // temp
      std::vector<double> eletemp(lm.size());
      DRT::UTILS::ExtractMyValues(*proc0temp, eletemp ,lm);

      // fill in center position and temperature of quad4
      // x1y1z1 T1
      int eleoffset = lid*4;
      std::vector<double> xcenter(3, 0.0);
      double Tcenter = 0.0;
      // loop over nodes of element
      for(int k=0; k<4; ++k)
      {
        DRT::Node** nodes = struele->Nodes();

        for(int d=0; d<3; ++d)
        {
          // pos
          double pos = (nodes[k]->X()[d] + eledisp[k*3 + d]) / LengthScaling();
          xcenter[d] += pos/4.0;
        }
        // temp
        Tcenter += eletemp[k]/4.0;
      }
      // center point is written
      for(int d=0; d<3; ++d)
      {
        // pos
        packeddata[eleoffset + d] = xcenter[d];
      }
      // temp
      packeddata[eleoffset + 3] = Tcenter;
    }

    // output
//    for(int i=0; i<(int)packeddata.size(); i++)
//    {
//      if(i%4 == 0)
//        std::cout << "\n element " << i/4 << " is packed: ";
//      std::cout << packeddata[i] << "  " ;
//    }
  }

  return;
}


/// extract structural interface values for a given full field
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::StrExtractInterfaceVal
(
  int interf,
  Teuchos::RCP<const Epetra_Vector> fullvector
)
{
  return structrowmapext_[interf]->ExtractCondVector(fullvector);
}


/// extract thermal interface values for a given full field
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::ThrExtractInterfaceVal
(
  int interf,
  Teuchos::RCP<Epetra_Vector> fullvector
)
{
  return thermorowmapext_[interf]->ExtractCondVector(fullvector);
}


/// apply interface tractions from the interface to a full field
void FS3I::UTILS::AeroCouplingUtils::StrInsertInterfaceVal
(
  int interf,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> force
)
{
  structrowmapext_[interf]->InsertCondVector(iforce, force);
  return;
}


/// apply interface heat flux from the interface to a full field
void FS3I::UTILS::AeroCouplingUtils::ThrInsertInterfaceVal
(
  int interf,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> force
)
{
  thermorowmapext_[interf]->InsertCondVector(iforce, force);
  return;
}
