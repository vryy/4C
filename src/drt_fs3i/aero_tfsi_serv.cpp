/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_serv.cpp

\brief Helper class for coupled simulations (INCA - BACI)

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "aero_tfsi_serv.H"
#include "aero_tfsi_delaunay.H"
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
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_io/io_pstream.H"

#include <Epetra_MpiComm.h>


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
serialslrownoderestr_(0),
D_(0),
M_(0),
shapefcn_(INPAR::MORTAR::shape_undefined)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if(ndim != 3)
    dserror("Coupling with Aero code is only implemented for 3D problems!!");

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
    istructdofredumap_.push_back(LINALG::AllreduceEMap(*(istructnewdis->DofRowMap())));
    ithermodofredumap_.push_back(LINALG::AllreduceEMap(*(ithermonewdis->DofRowMap())));

    // map extractor for the communication between the full thermo field and the newly created thermo surf discret
    structrowmapext_.push_back(Teuchos::rcp(new LINALG::MapExtractor(*(structdis->DofRowMap()),Teuchos::rcp(istructnewdis->DofRowMap(), false))));

    // map extractor for the communication between the full thermo field and the newly created thermo surf discret
    thermorowmapext_.push_back(Teuchos::rcp(new LINALG::MapExtractor(*(thermodis->DofRowMap()),Teuchos::rcp(ithermonewdis->DofRowMap(), false))));

  }

  // initial search radius in search tree
  maxmindist_ = 1.0e-2;

  // find largest dof id and node id in TSI problem for offset in mortar
  int localmastermaxdof = thermodis->DofRowMap()->MaxAllGID();
  localmastermaxdof++;
  com->MaxAll(&localmastermaxdof, &mastermaxdof_, 1);
  int localmastermaxnodeid = thermodis->NodeRowMap()->MaxAllGID();
  localmastermaxnodeid++;
  // to get a nice number, starting from at least 10^5
  localmastermaxnodeid = ((localmastermaxnodeid/100000)+1)*100000;
  com->MaxAll(&localmastermaxnodeid, &mastermaxnodeid_, 1);
  
//  // this can be interesting when interface is allowed to move
//  // epetra_vector with material configuration of row nodes is created
//  matconfig_ = LINALG::CreateVector(*(istructdis_->DofRowMap()), true);
//  for(int myl=0; myl<istructnoderowmap->NumMyElements(); myl++)
//  {
//    int gid = istructnoderowmap->GID(myl);
//    const DRT::Node* node = istructdis_->gNode(gid);
//
//    std::vector<int> lm;
//    lm.reserve(3);
//    // extract global dof ids
//    istructdis_->Dof(node, lm);
//
//    for (int a=0; a<3; a++)
//    {
//      int gid = lm[a];
//      double val = node->X()[a];
//      matconfig_->ReplaceGlobalValues(1, &val, &gid);
//    }
//  }



}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int,LINALG::Matrix<3,1> > FS3I::UTILS::AeroCouplingUtils::CurrentInterfacePos
(
  int interf,
  Teuchos::RCP<Epetra_Vector> redudisp
)
{
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleiter;

  //loop over all elements in the interface due to fully redundant storage
  for(eleiter=structreduelements_[interf].begin(); eleiter!=structreduelements_[interf].end(); eleiter++)
  {
    Teuchos::RCP<DRT::Element> tmpele = eleiter->second;

    const int* n = tmpele->NodeIds();

    // fill currentpositions
    for (int j=0; j < tmpele->NumNode(); j++)
    {
      const int gid = n[j];
      const DRT::Node* node = istructdis_[interf]->gNode(gid);
      std::vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      istructdis_[interf]->Dof(node, lm);
      std::vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;

      DRT::UTILS::ExtractMyValues(*redudisp,mydisp,lm);

      for (int a=0; a<3; a++)
      {
        currpos(a,0) = node->X()[a] + mydisp[a];
      }
      currentpositions[node->Id()] = currpos;
    }
  }

    return currentpositions;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::ProjectForceOnStruct
(
  int interf,
  Teuchos::RCP<Epetra_Vector> idispn,
  std::map<int, LINALG::Matrix<3,1> >& aerocoords,
  std::map<int, LINALG::Matrix<4,1> >& aeroforces,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  // Redistribute displacement of structural nodes on the interface to all processors --> redundant disp needed
  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_[interf],*(istructdis_[interf]->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_[interf],true);
  redudisp -> Import(*idispn,*interimpo,Add);

  //current position of structural nodes for the search tree (always 3 coordinates)
  std::map<int,LINALG::Matrix<3,1> > currentpositions = CurrentInterfacePos(interf, redudisp);

  //projection of the aero forces onto the closest element of each structural interface
  std::map<int, LINALG::Matrix<4,1> >::iterator forcesiter;
  for(forcesiter = aeroforces.begin(); forcesiter != aeroforces.end(); ++forcesiter)
  {
    LINALG::Matrix<3,1> forcecoords = aerocoords[forcesiter->first];
    //init of 3D search tree
    Teuchos::RCP<GEO::SearchTree> searchTree = Teuchos::rcp(new GEO::SearchTree(8));
    const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofEles(structreduelements_[interf], currentpositions);
    searchTree->initializeTreeSlideALE(rootBox, structreduelements_[interf], GEO::TreeType(GEO::OCTTREE));

    //search for near elements to the aero coords
    std::map<int,std::set<int> >  closeeles =
        searchTree->searchElementsInRadius(*istructdis_[interf],currentpositions,forcecoords,maxmindist_,0);

    //if no close elements could be found, try with a much larger radius and print a warning
    if (closeeles.empty())
    {
      cout<<"WARNING: no elements found in radius r="<<maxmindist_<<". Will try once with a bigger radius!"<<endl;
      closeeles = searchTree->searchElementsInRadius(*istructdis_[interf],currentpositions,forcecoords,100.0*maxmindist_,0);
      maxmindist_ *= 10.0;

      // if still no element is found, complain about it!
      if (closeeles.empty())
        dserror("No elements in a large radius! Should not happen!");
    }

    //search for the closest object, more exactly it's coordinates and the corresponding surface id
    LINALG::Matrix<3,1> minDistCoords;
    int eleid = GEO::nearest3DObjectInNode(istructdis_[interf], structreduelements_[interf], currentpositions,
        closeeles, forcecoords, minDistCoords);
    if(eleid == -1)
      dserror("Surface id couldn't be found. Weird! Hence, a distribution of the aero force is not performed!");

    //corresponding forces (from the aero code) to the coords that are about to be distributed
    LINALG::Matrix<4,1> force = forcesiter->second;

    //distribution of forces takes place in here
    DistributeForceToNodes(eleid, interf, currentpositions, minDistCoords, force, iforce, ithermoload);

  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::DistributeForceToNodes
(
  int eleid,
  int interf,
  std::map<int,LINALG::Matrix<3,1> >& currentpositions,
  LINALG::Matrix<3,1>& projpoint,
  LINALG::Matrix<4,1>& aeroforce,
  Teuchos::RCP<Epetra_Vector> distributedforces,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  Teuchos::RCP<DRT::Element> element = structreduelements_[interf][eleid];

  LINALG::Matrix<2,1> elecoord(true);
  const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element, currentpositions));

  //get coordinates of the projection point in parameter space of the element
  GEO::CurrentToSurfaceElementCoordinates(element->Shape(), xyze, projpoint, elecoord);

  const int numnode = element->NumNode();
  Epetra_SerialDenseVector   funct(numnode);

  // get shape functions of the element; evaluated at the projection point --> distribution
  DRT::UTILS::shape_function_2D(funct,elecoord(0,0),elecoord(1,0),element->Shape());

  for(int iter=0; iter<numnode; iter++)
  {
    // firstly, care about forces
    DRT::Node* snode = element->Nodes()[iter];

    // prepare assembly
    int num = 3;
    Epetra_SerialDenseVector val(num);
    std::vector<int> slm;
    slm.reserve(num);
    std::vector<int> lmowner(num);
    // extract global dof ids
    istructdis_[interf]->Dof(snode, slm);
    // aeroforce[0...2] are forces
    for(int k=0; k<num; k++)
    {
      val[k] = funct[iter] * aeroforce(k);
      lmowner[k] = snode->Owner();
    }

    // do assembly of structural forces
    LINALG::Assemble(*distributedforces,val,slm,lmowner);

    // secondly, care about heat fluxes
    DRT::Node* tnode = ithermodis_[interf]->gNode(snode->Id());

    // prepare assembly
    num = 1;
    Epetra_SerialDenseVector tval(num);
    std::vector<int> tlm;
    tlm.reserve(num);
    std::vector<int> tlmowner(num);
    // extract global dof ids
    ithermodis_[interf]->Dof(tnode, tlm);
    // aeroforce[3] is a thermal flux
    tval[0] = funct[iter] * aeroforce(3);
    tlmowner[0] = tnode->Owner();

    // do assembly of thermal fluxes
    LINALG::Assemble(*ithermoload,tval,tlm,tlmowner);

  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::BuildMortarCoupling
(
  int interf,
  Teuchos::RCP<Epetra_Vector> idispn,
  std::map<int, LINALG::Matrix<3,1> >& aerocoords
)
{
  // new mortar coupling matrices are needed in every time step
//  if(interf == 0)
//    DinvM_.clear();
  // CAREFUL HERE!
  // As long as the interface is fixed, it is enough to build the coupling matrices just once
  if(DinvM_.size() != 0)
    return;

  // Redistribute displacement of structural nodes on the interface to all processors --> redundant disp needed
  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_[interf],*(istructdis_[interf]->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_[interf],true);
  redudisp -> Import(*idispn,*interimpo,Add);

  //current position of fully redundant structural interface
  std::map<int,LINALG::Matrix<3,1> > currentpositions = CurrentInterfacePos(interf, redudisp);

  //init of 3D search tree with fluid cloud -> zeroth interface at the moment
  Teuchos::RCP<GEO::SearchTree> searchTree = Teuchos::rcp(new GEO::SearchTree(8));
  const LINALG::Matrix<3,2> nodeBox = GEO::getXAABBofPositions(aerocoords);
  searchTree->initializePointTree(nodeBox, aerocoords,GEO::TreeType(GEO::OCTTREE));

  // input data for Mortar interface which is build on each proc
  int interfID = 0;
  Teuchos::RCP<Epetra_MpiComm> serialcomm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF));
  // transfer of four degrees of freedom per node
  int dim = 3;
  int dofpernode = 1; // = 4; (for mechanical and thermal coupling)
  Teuchos::ParameterList input(DRT::Problem::Instance()->MortarCouplingParams());
  // brute force search is the only useful option because a single master element with its neighboring slave elements is in the interface
  input.set<string>("SEARCH_ALGORITHM","bruteforce");
  input.set<string>("PARALLEL_REDIST","no");
  // overwrite shape function type
  INPAR::TSI::BaciIncaCoupling tfsi_coupling =
      DRT::INPUT::IntegralValue<INPAR::TSI::BaciIncaCoupling>(DRT::Problem::Instance()->TSIDynamicParams(),"TFSI_COUPALGO");
  switch(tfsi_coupling)
  {
  case INPAR::TSI::mortar_mortar_std :
  case INPAR::TSI::proj_mortar_std :
  {
    input.set<string>("SHAPEFCN","std");
    shapefcn_ = INPAR::MORTAR::shape_standard;
    break;
  }
  case INPAR::TSI::mortar_mortar_dual :
  {
    input.set<string>("SHAPEFCN","dual");
    shapefcn_ = INPAR::MORTAR::shape_dual;
    break;
  }
  default:
    dserror("ERROR: Interface must either have dual or standard shape functions.");
  }

  // generate serialmasterdofrowmap_, no overlap is allowed for merged dof map (not needed for pure thermal coupling)
//  Teuchos::RCP<Epetra_Map> serialmasterdofrowmap = LINALG::MergeMap(istructdofredumap_, ithermodofredumap_,false);
  Teuchos::RCP<Epetra_Map> serialmasterdofrowmap = Teuchos::rcp(new Epetra_Map(*ithermodofredumap_[interf]));
  // gids and node ids on master and slave side must not be overlapping
  Teuchos::RCP<Epetra_Map> serialslavedofrowmap = Teuchos::rcp(new Epetra_Map((int)aerocoords.size()*dofpernode, mastermaxdof_+1, *serialcomm));

  // preparation for global AssembleDM
  std::set<int> slaverowdofsrestricted;
  std::set<int> slaverownodesrestricted;
  Teuchos::RCP<LINALG::SparseMatrix> D = Teuchos::rcp(new LINALG::SparseMatrix(*serialslavedofrowmap, 10));
  Teuchos::RCP<LINALG::SparseMatrix> M = Teuchos::rcp(new LINALG::SparseMatrix(*serialslavedofrowmap, 100));

  for(int lid=0; lid<istructdis_[interf]->NumMyColElements(); lid++)
  {
    DRT::Element* currele = istructdis_[interf]->lColElement(lid);
    // centroid of currele which is used for the search tree
    LINALG::Matrix<3,1> querypoint(true);

    // create a mortar interface for each structural interface element and the corresponding nodal fluid cloud
    MORTAR::MortarInterface interface(interfID, *serialcomm, dim, input, INPAR::MORTAR::redundant_all);

    // feeding master nodes to the interface including ghosted nodes
    // only consider the first 'dim' dofs
    DRT::Node** nodes = currele->Nodes();
    int numnode = currele->NumNode();
    for (int inode = 0; inode < numnode; ++inode)
    {
      DRT::Node* node = nodes[inode];
      // structural dofs are inserted first, followed by thermal dof
      std::vector<int> dofids(dofpernode);
      /* so far only one dof is needed (for thermo)
      for (int k=0;k<3;++k)
        dofids[k] = istructdis_->Dof(node)[k];
      dofids[3] = ithermodis_->Dof(ithermodis_->gNode(node->Id()))[0];
      */
      dofids[0] = ithermodis_[interf]->Dof(ithermodis_[interf]->gNode(node->Id()))[0];
      double nodalpos[3];
      for (int k=0;k<3;++k)
        nodalpos[k] = currentpositions[node->Id()](k);
      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(new MORTAR::MortarNode(node->Id(),
          nodalpos, /*node->Owner()*/ serialcomm->MyPID(), dofpernode, dofids, false));

      interface.AddMortarNode(mrtrnode);

      // calculate centroid of element
      for (int k=0;k<3;++k)
        querypoint(k) += nodalpos[k]/numnode;
    }

    // get an idea of the size of the element to have a good guess for the searchtree
    const int* nodeids = currele->NodeIds();
    LINALG::Matrix<3,1> edgelength = currentpositions[nodeids[0]];
    edgelength -= currentpositions[nodeids[2]];

    // TODO: FIND GOOD VALUE HERE: CURRENT 1.5
    double radius = 1.5 * sqrt( edgelength(0)*edgelength(0) + edgelength(1)*edgelength(1) + edgelength(2)*edgelength(2) );
    // search close fluid nodes and triangulate them
    std::vector<int> closefluidnodes = searchTree->searchPointsInRadius(aerocoords, querypoint, radius);
    if(closefluidnodes.size() < 9)
      dserror("Radius to find fluid nodes is probably to small. There are just %i nodes to triangulate!", closefluidnodes.size());

    // input data for delaunay triangulation
    std::map<int, LINALG::Matrix<3,1> > closeaerocoords;
    for(size_t i=0;i<closefluidnodes.size(); i++)
      closeaerocoords[closefluidnodes[i]] = aerocoords[closefluidnodes[i]];
    std::map<int, LINALG::Matrix<3,1> > structure;
    std::map<int, std::vector<int> > tris;
    // arbitrarily, the first three nodes of the structural element are chosen, needed for projection before triangulation
    for(int k=0; k<3; k++)
      structure[k] = currentpositions[nodeids[k]];

    // do Delaunay triangulation of close fluid nodes
    FS3I::UTILS::DelaunayTriangulation(closeaerocoords, structure, tris);

    // feeding slave nodes to the interface including ghosted nodes
    // Note: offset is needed for slave dofs and slave node ids
    for (size_t nodeiter = 0; nodeiter<closefluidnodes.size(); ++nodeiter)
    {
      double pos[3];
      int nodeid = closefluidnodes[nodeiter];
      for (int k=0;k<3;++k) pos[k] = aerocoords[nodeid](k);
      std::vector<int> dofids(dofpernode);
      for (int k=0;k<dofpernode;++k) dofids[k] = (mastermaxdof_+1) + nodeid*dofpernode+k;
      Teuchos::RCP<MORTAR::MortarNode> mrtrnode =
          Teuchos::rcp(new MORTAR::MortarNode(nodeid+mastermaxnodeid_, pos, serialcomm->MyPID(), dofpernode, dofids, true));

      interface.AddMortarNode(mrtrnode);
    }

    // max master element ID needed for unique eleIDs in interface discretization
    // will be used as offset for slave elements
    int EleOffset = istructdis_[interf]->ElementRowMap()->MaxAllGID()+1;

    // feed single master element to the interface
    {
      Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                  new MORTAR::MortarElement(currele->Id(), /*currele->Owner()*/ serialcomm->MyPID(), currele->Shape(),
                      currele->NumNode(), currele->NodeIds(), false));

      interface.AddMortarElement(mrtrele);
    }

    // feeding slave elements to the interface; take care about offset (just within mortar)
    for (size_t elemiter = 0; elemiter < tris.size(); ++elemiter)
    {
      std::vector<int> nodeids = tris[elemiter];
      std::vector<int> nodeidswithoffset(nodeids.size());
      for(size_t n=0; n<nodeids.size(); n++)
        nodeidswithoffset[n] = nodeids[n] + mastermaxnodeid_;
      Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(new MORTAR::MortarElement(elemiter + EleOffset,
          serialcomm->MyPID(), DRT::Element::tri3, (int)nodeidswithoffset.size(), &nodeidswithoffset[0], true));

      interface.AddMortarElement(mrtrele);
    }

    // finalize the mortar interface
    interface.FillComplete();

    // store maps
    Teuchos::RCP<Epetra_Map> partialslavedofrowmap_  = interface.SlaveRowDofs();
    Teuchos::RCP<Epetra_Map> partialmasterdofrowmap_ = interface.MasterRowDofs();

    // all the following stuff has to be done once in setup
    // in order to get initial D and M

    // interface displacement (=0) has to be merged from slave and master discretization
    Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*partialmasterdofrowmap_,*partialslavedofrowmap_, false);
    Teuchos::RCP<Epetra_Vector> dispn = LINALG::CreateVector(*dofrowmap, true);

    // set displacement state in mortar interface
    interface.SetState("displacement", dispn);

    //in the following two steps MORTAR does all the work
    interface.Initialize();
    interface.Evaluate();

    // restrict mortar treatment to actual meshtying zone (drop untied slave nodes in the additional boundary layer)
    {
      // Step 1: detect tied slave nodes on all interfaces
      int globalfounduntied = 0; // Note: MPI_COMM_SELF
      interface.DetectTiedSlaveNodes(globalfounduntied);

      // get out of here if the whole slave surface is tied
      if (globalfounduntied != 0 and shapefcn_ == INPAR::MORTAR::shape_standard)
      {
        interface.RestrictSlaveSets();
      }
      // extract tied dofs
      Teuchos::RCP<Epetra_Map> slrowdofs = interface.SlaveRowDofs();
      for(int i=0; i<slrowdofs->NumMyElements(); i++)
        slaverowdofsrestricted.insert(slrowdofs->GID(i));
      // extract tied nodes
      Teuchos::RCP<Epetra_Map> slrownodes = interface.SlaveRowNodes();
      for(int i=0; i<slrownodes->NumMyElements(); i++)
        slaverownodesrestricted.insert(slrownodes->GID(i)-mastermaxnodeid_);
    }

    // preparation for local AssembleDM
    Teuchos::RCP<LINALG::SparseMatrix> dmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*interface.SlaveRowDofs(), 10));
    Teuchos::RCP<LINALG::SparseMatrix> mmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*interface.SlaveRowDofs(), 100));
    interface.AssembleDM(*dmatrix, *mmatrix);
    dmatrix->Complete();
    mmatrix->Complete(*partialmasterdofrowmap_,*interface.SlaveRowDofs());

    // add local D and M matrix entries to the global ones
    D->Add(*dmatrix,false,1.0,1.0);
    M->Add(*mmatrix,false,1.0,1.0);

  } // END LOOP OVER ALL ELEMENTS

  // build DINVM_ once at the end

  // Complete() global Mortar matrices
  D->Complete();
  M->Complete(*serialmasterdofrowmap, *serialslavedofrowmap);

  // restricted slave node map to extract correct data from aeroforces later
  std::vector<int> slrownodes(slaverownodesrestricted.begin(),slaverownodesrestricted.end());
  Teuchos::RCP<Epetra_Map> serialslrownoderestr = Teuchos::rcp(new Epetra_Map(-1,(int)slrownodes.size(),&slrownodes[0],0,*serialcomm));
  serialslrownoderestr_.push_back(serialslrownoderestr);

  if(shapefcn_ == INPAR::MORTAR::shape_dual)
  {
    // Build Dinv
    Teuchos::RCP<LINALG::SparseMatrix> Dinv = Teuchos::rcp(new LINALG::SparseMatrix(*D));

    // extract diagonal of invd into diag
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*serialslavedofrowmap,true);
    Dinv->ExtractDiagonalCopy(*diag);

    // set zero diagonal values to dummy 1.0
    for (int i=0;i<diag->MyLength();++i)
      if ((*diag)[i]==0.0) (*diag)[i]=1.0;

    // scalar inversion of diagonal values
    diag->Reciprocal(*diag);
    Dinv->ReplaceDiagonalValues(*diag);

    Dinv->Complete( D->RangeMap(), D->DomainMap() );

    // MLMultiply is not prepared for nested parallelism  --> Multiply must be enough
    Teuchos::RCP<LINALG::SparseMatrix> DinvM = Multiply(*Dinv,false,*M,false,true);

    DinvM_.push_back(DinvM);
  }
  else
  {
    // restricted slave dof map must be used
    std::vector<int> slrowdofs(slaverowdofsrestricted.begin(),slaverowdofsrestricted.end());
    Teuchos::RCP<Epetra_Map> serialslrowdofrestr = Teuchos::rcp(new Epetra_Map(-1,(int)slrowdofs.size(),&slrowdofs[0],0,*serialcomm));

    // copy matrix to correct row map
    Teuchos::RCP<LINALG::SparseMatrix> D_restr = Teuchos::rcp(new LINALG::SparseMatrix(*serialslrowdofrestr, 10));
    Teuchos::RCP<LINALG::SparseMatrix> M_restr = Teuchos::rcp(new LINALG::SparseMatrix(*serialslrowdofrestr, 100));

    int D_length = D->MaxNumEntries();
    int D_numentries = 0;
    double D_values[D_length];
    int D_indices[D_length];
    int M_length = M->MaxNumEntries();
    int M_numentries = 0;
    double M_values[M_length];
    int M_indices[M_length];
    for(int irow=0; irow<serialslrowdofrestr->NumMyElements(); irow++)
    {
      int rowgid = serialslrowdofrestr->GID(irow);
      (D->EpetraMatrix())->ExtractGlobalRowCopy(rowgid, D_length, D_numentries, &D_values[0], &D_indices[0]);
      (D_restr->EpetraMatrix())->InsertGlobalValues(rowgid, D_numentries, D_values, D_indices);
      (M->EpetraMatrix())->ExtractGlobalRowCopy(rowgid, M_length, M_numentries, &M_values[0], &M_indices[0]);
      (M_restr->EpetraMatrix())->InsertGlobalValues(rowgid, M_numentries, M_values, M_indices);
    }
    // now we have correctly restricted matrices
    D_restr->Complete(*serialslrowdofrestr,*serialslrowdofrestr);
    M_restr->Complete(*serialmasterdofrowmap, *serialslrowdofrestr);

    // just store everything
    D_.push_back(D_restr);
    M_.push_back(M_restr);
  }

  return;
}

// for dual shape functions in mortar
void FS3I::UTILS::AeroCouplingUtils::TransferFluidLoadsToStructDual
(
  int interf,
  std::map<int, LINALG::Matrix<4,1> >& aeroforces,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{

  Teuchos::RCP<Epetra_Vector> slavevalues = LINALG::CreateVector(DinvM_[interf]->RowMap(),true);
  // version with fluxes for the additional boundary layer (which is not conservative)
  /*
  if(slavevalues->MyLength() > (int)aeroforces.size())
    dserror("size of slavevalues (%i) is larger than incoming data (%i)", slavevalues->MyLength(), (aeroforces.size()));

  // so far only thermo fluxes need to be transferred (4th entry in aeroforces)
  for(int i=0; i<slavevalues->MyLength(); i++)
    (*slavevalues)[i] = aeroforces[serialslrownoderestr_[interf]->GID(i)](3);
  */

  // so far only thermo fluxes need to be transferred (4th entry in aeroforces)
  // NOTE: There is no dofset available for artificially created fluid nodes
  // Assumption: Meshtying restriction is only needed for additional boundary layer
  // Hence: Here in serial (LID=GID), values from aeroforces can directly be inserted into slavevalues vector
  std::map<int, LINALG::Matrix<4,1> >::const_iterator iter;
  for(iter=aeroforces.begin(); iter!=aeroforces.end(); ++iter)
  {
    (*slavevalues)[iter->first] = iter->second(3);
  }

  // currently the interface only lives within MPI_COMM_SELF --> not parallel
  Teuchos::RCP<Epetra_Vector> tmpmastervalues = FluidToStruct(interf, slavevalues);

  // hard copy of each value which actually lives on this proc; ithermoload is of type dofrowmap
  // lid (serial vector) --> gid (identical)--> lid (distributed vector)
  for(int lid=0; lid<ithermoload->MyLength(); lid++)
    (*ithermoload)[lid] = (*tmpmastervalues)[tmpmastervalues->Map().LID( ithermoload->Map().GID(lid) )];

  // forces are zero so far
  iforce->PutScalar(0.0);

  slavevalues = Teuchos::null;
  tmpmastervalues = Teuchos::null;
  return;
}


// for standard shape functions in mortar
void FS3I::UTILS::AeroCouplingUtils::TransferFluidLoadsToStructStd
(
  int interf,
  std::map<int, LINALG::Matrix<4,1> >& aeroforces,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  Teuchos::RCP<Epetra_Vector> slavevalues = LINALG::CreateVector(D_[interf]->RowMap(),true);
  /*
  if(slavevalues->MyLength() > (int)aeroforces.size())
    dserror("size of slavevalues (%i) is larger than incoming data (%i)", slavevalues->MyLength(), (aeroforces.size()));

  // so far only thermo fluxes need to be transferred (4th entry in aeroforces)
  for(int i=0; i<slavevalues->MyLength(); i++)
    (*slavevalues)[i] = aeroforces[serialslrownoderestr_[interf]->GID(i)](3);
  */

  // so far only thermo fluxes need to be transferred (4th entry in aeroforces)
  // NOTE: There is no dofset available for artificially created fluid nodes
  // Assumption: Meshtying restriction is only needed for additional boundary layer
  // Hence: Here in serial (LID=GID), values from aeroforces can directly be inserted into slavevalues vector
  std::map<int, LINALG::Matrix<4,1> >::const_iterator iter;
  for(iter=aeroforces.begin(); iter!=aeroforces.end(); ++iter)
  {
    (*slavevalues)[iter->first] = iter->second(3);
  }

  const int linsolvernumber = 9;
  // get solver parameter list of linear solver
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  Teuchos::RCP<Epetra_MpiComm> serialcomm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF));

  // create solver
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(solverparams, *serialcomm, NULL));

  // solve D^T lambda = forces_fl for lambda
  Teuchos::RCP<Epetra_Operator> Dtransp = LINALG::Transpose(*D_[interf]->EpetraMatrix());
  Teuchos::RCP<Epetra_Vector> lambda = LINALG::CreateVector(Dtransp->OperatorDomainMap(), true);
  solver->Solve(Dtransp, lambda, slavevalues, true, true);

  // M^T lambda = forces_str
  // Note: DomainMap of M is equal to RowMap of M^T
  Teuchos::RCP<Epetra_Vector> tmpmastervalues = Teuchos::rcp(new Epetra_Vector(M_[interf]->DomainMap(), true));
  M_[interf]->Multiply(true, *lambda, *tmpmastervalues);

  // hard copy of each value which actually lives on this proc; ithermoload is of type dofrowmap
  // lid (serial vector) --> gid (identical)--> lid (distributed vector)
  for(int lid=0; lid<ithermoload->MyLength(); lid++)
    (*ithermoload)[lid] = (*tmpmastervalues)[tmpmastervalues->Map().LID( ithermoload->Map().GID(lid) )];

  slavevalues = Teuchos::null;
  tmpmastervalues = Teuchos::null;
  return;
}


/*----------------------------------------------------------------------*
 | master to slave transfer helper                          ghamm 11/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::StructToFluid
(
  int interf,
  Teuchos::RCP<Epetra_Vector> mv
) const
{
  if ( not (DinvM_[interf]->DomainMap().SameAs( mv->Map())) )
    dserror("Vector with master dof map expected");

  Teuchos::RCP<Epetra_Vector> sv = Teuchos::rcp(new Epetra_Vector(DinvM_[interf]->RowMap()));

  if ( DinvM_[interf]->Multiply( false, *mv, *sv ) )
    dserror( "D^{-1}*M*mv multiplication failed" );

  return sv;
}


/*----------------------------------------------------------------------*
 | slave to master transfer helper                          ghamm 11/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::FluidToStruct
(
  int interf,
  Teuchos::RCP<Epetra_Vector> sv
) const
{
  if ( not (DinvM_[interf]->RowMap().SameAs( sv->Map())) )
    dserror("Vector with slave dof map expected");

  Teuchos::RCP<Epetra_Vector> mv = Teuchos::rcp(new Epetra_Vector(DinvM_[interf]->DomainMap()));
  if (DinvM_[interf]->Multiply(true, *sv, *mv))
    dserror( "-(D^{-1}*M)^{T}*sv multiplication failed" );

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::PackData
(
  int interf,
  Teuchos::RCP<Epetra_Vector> idispnp,
  Teuchos::RCP<Epetra_Vector> itempnp,
  std::vector<double>& packeddata,
  bool writedata
)
{
  // make all interface data full redundant
  // idispn
  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_[interf],*(istructdis_[interf]->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_[interf],true);
  redudisp -> Import(*idispnp,*interimpo,Add);

  // itemp
  Teuchos::RCP<Epetra_Import> interimpothermo = Teuchos::rcp(new Epetra_Import(*ithermodofredumap_[interf],*(ithermodis_[interf]->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redutemp = LINALG::CreateVector(*ithermodofredumap_[interf],true);
  redutemp -> Import(*itempnp,*interimpothermo,Add);

  //current position of structural nodes
  std::map<int,LINALG::Matrix<3,1> > currentpositions = CurrentInterfacePos(interf, redudisp);

  int dim = 3;

  std::vector<double> INCAdata;

  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleiter;
  for(eleiter=structreduelements_[interf].begin(); eleiter!=structreduelements_[interf].end(); eleiter++)
  {
    Teuchos::RCP<DRT::Element> tmpele = eleiter->second;
    const int* nodeids = tmpele->NodeIds();

    std::vector<double> midposition(true);
    midposition.reserve(3);
    for (int k=0;k<3;++k) midposition[k] = 0.0;

    std::vector<double> midtemp(true);
    midtemp.reserve(1);
    midtemp[0]=0.0;

    // insert nodal positions and temperature for each quad element
    for(int j=0; j<4; j++)
    {
      int nodalgid = nodeids[j];

      // insert nodal positions
      for(int i=0; i<dim; i++)
      {
        midposition[i] += (0.25 * currentpositions[nodalgid](i));
      }

      // insert nodal temperature
      DRT::Node* tnode = ithermodis_[interf]->gNode(nodalgid);
      std::vector<int> lmt;
      lmt.reserve(1);
      ithermodis_[interf]->Dof(tnode, lmt);

      int lid = ithermodofredumap_[interf]->LID(lmt[0]);

      midtemp[0] += (0.25 * (redutemp->operator[](lid)));

    } // end of loop over four nodes

    // add the gathered data
    for(int i=0; i<dim; i++)
    {
      INCAdata.push_back(midposition[i]);
    }
    INCAdata.push_back(midtemp[0]);

  } // end elements loop

  // do some reordering for INCA (xyzTxyzTxyzTxyzT --> xxxxyyyyzzzzTTTT)
  size_t size = INCAdata.size();
  packeddata.resize(size);
  size_t sizequarter = size/4;
  for(size_t out=0; out<sizequarter; out++)
  {
    for(size_t in=0; in<4; in++)
    {
      packeddata[in*sizequarter + out] = INCAdata[out*4 + in];
    }
  }

  if(istructdis_[0]->Comm().MyPID() == 0 && interf == 0 && writedata)
  {
    FILE *outFile;
    outFile = fopen("interfaceTemp.txt", "a");
    fprintf(outFile, "%.8e\n", packeddata[699]);
    fclose(outFile);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::TransferStructValuesToFluidDual
(
  int interf,
  int sizeofphysicaldomain,
  std::map<int, LINALG::Matrix<3,1> >& aerocoords,
  Teuchos::RCP<Epetra_Vector> itempnp,
  std::vector<double>& packeddata,
  bool writedata /* = true*/
)
{
  // make all interface data full redundant (after mapping with StructToFluid!!!!!!) here: fixed interface
  // idispn
//  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_[interf],*(istructdis_[interf]->DofRowMap())));
//  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_[interf],true);
//  redudisp -> Import(*idispnp,*interimpo,Add);

  // itemp
  Teuchos::RCP<Epetra_Import> interimpothermo = Teuchos::rcp(new Epetra_Import(*ithermodofredumap_[interf],*(ithermodis_[interf]->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redutemp = LINALG::CreateVector(*ithermodofredumap_[interf],true);
  redutemp -> Import(*itempnp,*interimpothermo,Add);

  // currently the interface only lives within MPI_COMM_SELF --> not parallel
  Teuchos::RCP<Epetra_Vector> slavevalues = StructToFluid(interf, redutemp);

  if(slavevalues->MyLength() != (int)aerocoords.size())
    dserror("this is not possible because aerocoords were used for building mortar matrices");

  if(sizeofphysicaldomain > (int)aerocoords.size())
    dserror("Physical domain is larger (%i) than size of aerocoords (%i); is the additional boundary layer missing?", sizeofphysicaldomain, aerocoords.size());

  std::vector<double> INCAdata;
  for(int i=0; i<sizeofphysicaldomain; i++)
  {
    // add the gathered data
    for(int dim=0; dim<3; dim++)
    {
      INCAdata.push_back(aerocoords[i](dim));
    }
    INCAdata.push_back((*slavevalues)[i]);
  }

  // do some reordering for INCA (xyzTxyzTxyzTxyzT --> xxxxyyyyzzzzTTTT)
  size_t size = INCAdata.size();
  packeddata.resize(size);
  size_t sizequarter = size/4;
  for(size_t out=0; out<sizequarter; out++)
  {
    for(size_t in=0; in<4; in++)
    {
      packeddata[in*sizequarter + out] = INCAdata[out*4 + in];
    }
  }

  if(istructdis_[0]->Comm().MyPID() == 0 && interf == 0 && writedata)
  {
    FILE *outFile;
    outFile = fopen("interfaceTemp.txt", "a");
    fprintf(outFile, "%.8e\n", packeddata[699]);
    fclose(outFile);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::TransferStructValuesToFluidStd
(
  int interf,
  int sizeofphysicaldomain,
  std::map<int, LINALG::Matrix<3,1> >& aerocoords,
  Teuchos::RCP<Epetra_Vector> itempnp,
  std::vector<double>& packeddata,
  bool writedata /* = true*/
)
{
  // make all interface data full redundant (after mapping with StructToFluid!!!!!!) here: fixed interface
  // idispn
//  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_[interf],*(istructdis_[interf]->DofRowMap())));
//  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_[interf],true);
//  redudisp -> Import(*idispnp,*interimpo,Add);

  // itemp
  Teuchos::RCP<Epetra_Import> interimpothermo = Teuchos::rcp(new Epetra_Import(*ithermodofredumap_[interf],*(ithermodis_[interf]->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redutemp = LINALG::CreateVector(*ithermodofredumap_[interf],true);
  redutemp -> Import(*itempnp,*interimpothermo,Add);

  const int linsolvernumber = 9;
  // get solver parameter list of linear solver
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  Teuchos::RCP<Epetra_MpiComm> serialcomm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF));

  // create solver
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(solverparams, *serialcomm, NULL));

  // M d_str = rhs
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(M_[interf]->RangeMap(), true));
  M_[interf]->Multiply(false, *redutemp, *rhs);

  // solve D temp_fl = rhs for temp_fl
  Teuchos::RCP<Epetra_Vector> slavevalues = LINALG::CreateVector(D_[interf]->DomainMap(), true);
  solver->Solve(D_[interf]->EpetraOperator(), slavevalues, rhs, true, true);

  if(sizeofphysicaldomain > (int)aerocoords.size())
    dserror("Physical domain is larger (%i) than size of aerocoords (%i); is the additional boundary layer missing?", sizeofphysicaldomain, aerocoords.size());

  // size of physical domain is always within restricted zone --> restriction does not have to be considered here
  std::vector<double> INCAdata;
  for(int i=0; i<sizeofphysicaldomain; i++)
  {
    // add the gathered data
    for(int dim=0; dim<3; dim++)
    {
      INCAdata.push_back(aerocoords[i](dim));
    }
    INCAdata.push_back((*slavevalues)[i]);
  }

  // do some reordering for INCA (xyzTxyzTxyzTxyzT --> xxxxyyyyzzzzTTTT)
  size_t size = INCAdata.size();
  packeddata.resize(size);
  size_t sizequarter = size/4;
  for(size_t out=0; out<sizequarter; out++)
  {
    for(size_t in=0; in<4; in++)
    {
      packeddata[in*sizequarter + out] = INCAdata[out*4 + in];
    }
  }

  if(istructdis_[0]->Comm().MyPID() == 0 && interf == 0 && writedata)
  {
    FILE *outFile;
    outFile = fopen("interfaceTemp.txt", "a");
    fprintf(outFile, "%.8e\n", packeddata[699]);
    fclose(outFile);
  }

  return;
}


/// extract structural interface values for a given full field
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::StrExtractInterfaceVal
(
  int interf,
  Teuchos::RCP<Epetra_Vector> fullvector
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
void FS3I::UTILS::AeroCouplingUtils::StrApplyInterfaceVal
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
void FS3I::UTILS::AeroCouplingUtils::ThrApplyInterfaceVal
(
  int interf,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> force
)
{
  thermorowmapext_[interf]->InsertCondVector(iforce, force);
  return;
}
