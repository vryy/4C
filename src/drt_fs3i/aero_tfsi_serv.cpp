/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_serv.cpp

\brief Helper class for coupled simulations (INCA - BACI)

\level 3

\maintainer Georg Hammerl
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "aero_tfsi_serv.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element_geo_decoupl.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"

#include <Epetra_FEVector.h>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 | implementation with 4 dofs per mortar node = TFSI        ghamm 12/11 |
 *----------------------------------------------------------------------*/
FS3I::UTILS::AeroCouplingUtils::AeroCouplingUtils
(
  Teuchos::RCP<DRT::Discretization> structdis,
  Teuchos::RCP<DRT::Discretization> thermodis
) :
dim_(3),
myrank_(structdis->Comm().MyPID()),
istructdis_(0),
ithermodis_(0),
istructdofredumap_(0),
ithermodofredumap_(0),
structrowmapext_(0),
thermorowmapext_(0),
DinvM_(0),
dofpernode_(4),
mechanicalcoupling_(true),
thermalcoupling_(true),
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
        if (currnode->Owner() == myrank_)
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
        if (currele->Owner() == myrank_)
        {
          // each interface stores its row elements in the interface
          eleids.push_back(currele->Id() );

          // structural surface elements cannot be distributed --> BELE3_3 element is used
          Teuchos::RCP<DRT::Element> istructele = DRT::UTILS::Factory("BELE3_3","Polynomial", currele->Id(), currele->Owner());
          istructele->SetNodeIds(currele->NumNode(), currele->NodeIds());
          istructnewdis->AddElement( istructele );
          // thermo interface elements --> represented by BELE3_1 element
          Teuchos::RCP<DRT::Element> ithermoele = DRT::UTILS::Factory("BELE3_1","Polynomial", currele->Id(), currele->Owner());
          ithermoele->SetNodeIds(currele->NumNode(), currele->NodeIds());
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

    // set useful parameters for mortar coupling
    Teuchos::ParameterList input;
    const Teuchos::ParameterList& mortar = DRT::Problem::Instance()->MortarCouplingParams();
    const Teuchos::ParameterList& cmortar = DRT::Problem::Instance()->ContactDynamicParams();
    input.setParameters(cmortar);
    input.setParameters(mortar);
    MortarParams(input);

    // create mortar interface for each coupling interface with INCA
    Teuchos::RCP<MORTAR::MortarInterface> interface =
        MORTAR::MortarInterface::Create(interf, lcomm, ndim, input,
            INPAR::MORTAR::redundant_master);

    // feeding master nodes to the interface including ghosted nodes
    // transfer of four degrees of freedom per node
    std::vector<int> dofids(dofpernode_);
    for(int lid=0; lid<istructdis_[interf]->NumMyColNodes(); lid++)
    {
      DRT::Node* node = istructdis_[interf]->lColNode(lid);
      // get the dofs of the structural node
      for (int k=0; k<ndim; ++k)
       dofids[k] = istructdis_[interf]->Dof(node)[k];
      // fill in remaining thermal dof
      dofids[ndim] = ithermodis_[interf]->Dof(ithermodis_[interf]->lColNode(lid))[0];

      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(new MORTAR::MortarNode(node->Id(),
         node->X(), node->Owner(), dofpernode_, dofids, false));

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


/*----------------------------------------------------------------------*
 | implementation with 3 or 1 dof(s) per mortar node =      ghamm 06/14 |
 |   FSI or Conjugate heat transfer problem                             |
 *----------------------------------------------------------------------*/
FS3I::UTILS::AeroCouplingUtils::AeroCouplingUtils
(
  Teuchos::RCP<DRT::Discretization> dis,
  bool mechanical
) :
dim_(3),
myrank_(dis->Comm().MyPID()),
istructdis_(0),
ithermodis_(0),
istructdofredumap_(0),
ithermodofredumap_(0),
structrowmapext_(0),
thermorowmapext_(0),
DinvM_(0),
dofpernode_(mechanical ? 3 : 1),
mechanicalcoupling_(mechanical ? true : false),
thermalcoupling_(mechanical ? false : true),
mastermaxdof_(0),
mastermaxnodeid_(0),
mastermaxeleid_(0)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if(ndim != 3)
    dserror("Coupling with Aero code is only implemented for 3D problems!!");

  // build references to the discret, map and map extractors for mechanical coupling ...
  std::vector<Teuchos::RCP<DRT::Discretization> >* idis  = &istructdis_;
  std::vector<Teuchos::RCP<const Epetra_Map> >* idofredumap = &istructdofredumap_;
  std::vector<Teuchos::RCP<const LINALG::MapExtractor> >* rowmapext = &structrowmapext_;
  std::vector<Teuchos::RCP<const Epetra_Import> >* iproc0importer = &istructproc0importer_;
  // ... or for thermal coupling
  if(mechanical == false)
  {
    idis = &ithermodis_;
    idofredumap = &ithermodofredumap_;
    rowmapext = &thermorowmapext_;
    iproc0importer = &ithermoproc0importer_;
  }

  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  lengthscaling_ = Teuchos::getIntegralValue<double>(tsidyn,"TFSI_length_unit");
  timescaling_ = Teuchos::getIntegralValue<double>(tsidyn,"TFSI_time_unit");

  // declare struct objects in interface
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > > structelements;
  std::map<int, std::map<int, DRT::Node*> > dummy2;  // dummy map
  std::map<int, std::map<int, DRT::Node*> > structgnodes; // col map of structure nodes

  //initialize struct objects in interface
  DRT::UTILS::FindConditionObjects(*dis, dummy2, structgnodes, structelements,"FSICoupling");

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit;
  std::map<int, Teuchos::RCP<DRT::Element> >::iterator eit;

  int max_id = 0;
  //determine number of coupling surfaces
  for ( meit=structelements.begin(); meit != structelements.end(); meit++ )
  {
    if (meit->first > max_id)
      max_id = meit->first;
  }
  dis->Comm().MaxAll(&max_id, &maxid_, 1);

  //create a new discretization of the structural/thermal interface
  if (!dis->Filled())
  {
    dis->FillComplete();
  }

  // initialize the new discretization
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(dis->Comm().Clone());
  // structure
  const std::string discret_name = "idis";

  for(int interf=0; interf <= maxid_; interf++)
  {
    // initialize new interface discretization
    Teuchos::RCP<DRT::Discretization> inewdis = Teuchos::rcp(new DRT::Discretization(discret_name,com));

    std::map<int, std::map<int, DRT::Node*> >::iterator mnit;
    std::map<int, DRT::Node* >::iterator nit;

    std::vector<int> nodeids;

    // nodes filled in new interface discretization
    {
      std::map<int, DRT::Node*> structgnodesinterf = structgnodes[interf];
      for ( nit=structgnodesinterf.begin(); nit != structgnodesinterf.end(); nit++ )
      {
        DRT::Node* currnode = (*nit).second;
        //ids for row map are gathered
        if (currnode->Owner() == myrank_)
          nodeids.push_back(currnode->Id());

        //discretization is fed
        inewdis->AddNode(Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
      }
    }
    //row node map of interface
    Teuchos::RCP<Epetra_Map> inoderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,inewdis->Comm()));
    //fully overlapping node map
    Teuchos::RCP<Epetra_Map> irednodecolmap = LINALG::AllreduceEMap(*inoderowmap);

    //care about eles in the interface
    std::stringstream elename;
    elename << "BELE3_" << dofpernode_;

    std::vector<int> eleids;
    //eles filled in new interface discretization
    {
      std::map<int, Teuchos::RCP<DRT::Element> > structelementsinterf = structelements[interf];
      for ( eit=structelementsinterf.begin(); eit != structelementsinterf.end(); eit++ )
      {
        Teuchos::RCP<DRT::Element> currele = eit->second;
        if (currele->Owner() == myrank_)
        {
          // each interface stores its row elements in the interface
          eleids.push_back(currele->Id() );

          // BELE3_3 or BELE3_1 element is used (mechanical or thermal coupling)
          Teuchos::RCP<DRT::Element> iele = DRT::UTILS::Factory(elename.str(),"Polynomial", currele->Id(), currele->Owner());
          iele->SetNodeIds(currele->NumNode(), currele->NodeIds());
          inewdis->AddElement( iele );
        }
      }
    }

    //row ele map of this interface
    Teuchos::RCP<Epetra_Map> ielerowmap = Teuchos::rcp(new Epetra_Map(-1,eleids.size(),&eleids[0],0,inewdis->Comm()));
    //fully overlapping ele map
    Teuchos::RCP<Epetra_Map> iredelecolmap = LINALG::AllreduceEMap(*ielerowmap);

    //do the fully overlapping ghosting of the master side of the FSI interface to have everything redundant
    inewdis->ExportColumnNodes(*irednodecolmap);
    inewdis->ExportColumnElements(*iredelecolmap);

    //find out if we are parallel; needed for TransparentDofSet
    bool parallel = false;
    if(inewdis->Comm().NumProc() != 1)
      parallel = true;

    //dofs of the original discretizations are used to set same dofs for the new interface discretization
    Teuchos::RCP<DRT::DofSet> newdofset1=Teuchos::rcp(new DRT::TransparentDofSet(dis,parallel));
    inewdis->ReplaceDofSet(newdofset1);
    newdofset1=Teuchos::null;

    //final fill complete to reorganize everything in the discretization
    inewdis->FillComplete(true, false, false);

    //storage of the newly created interface discretization
    idis->push_back(inewdis);

    //storage of useful dof maps
    idofredumap->push_back(Teuchos::rcp(inewdis->DofColMap(), false));

    // map extractor for the communication between the full structural/thermo field and the newly created structural/thermo surf discret
    rowmapext->push_back(Teuchos::rcp(new LINALG::MapExtractor(*(dis->DofRowMap()),Teuchos::rcp(inewdis->DofRowMap(), false))));

    // importer in order to contract result values on proc0 for communication with INCA (proc0 gets everything, other procs empty)
    Teuchos::RCP<Epetra_Map> proc0datamap = LINALG::AllreduceEMap(*inewdis->DofRowMap(),0);
    iproc0importer->push_back(Teuchos::rcp(new Epetra_Import (*proc0datamap,*(*idis)[interf]->DofRowMap())));

    // find largest dof/node/element id in TSI problem for offset in mortar
    mastermaxdof_ = std::max(mastermaxdof_, dis->DofRowMap()->MaxAllGID() + 1);
    mastermaxnodeid_ = std::max(mastermaxnodeid_, dis->NodeRowMap()->MaxAllGID() + 1);
    mastermaxeleid_ = std::max(mastermaxeleid_, (*idis)[interf]->ElementRowMap()->MaxAllGID() + 1);

    // build mortar interface and fill it initially with fully redundant master side (=structural side)
    const Epetra_Comm& lcomm = (*idis)[interf]->Comm();

    // get mortar coupling parameters
    Teuchos::ParameterList input;
    const Teuchos::ParameterList& mortar = DRT::Problem::Instance()->MortarCouplingParams();
    const Teuchos::ParameterList& cmortar = DRT::Problem::Instance()->ContactDynamicParams();
    input.setParameters(cmortar);
    input.setParameters(mortar);
    MortarParams(input);

    // create mortar interface for each coupling interface with INCA
    Teuchos::RCP<MORTAR::MortarInterface> interface = MORTAR::MortarInterface::Create(
        interf, lcomm, ndim, input, INPAR::MORTAR::redundant_master);

    // feeding master nodes to the interface including ghosted nodes
    // transfer of three degrees of freedom per node
    std::vector<int> dofids(dofpernode_);
    for(int lid=0; lid<(*idis)[interf]->NumMyColNodes(); lid++)
    {
      DRT::Node* node = (*idis)[interf]->lColNode(lid);
      // get the dofs of the structural node
      for (int k=0; k<dofpernode_; ++k)
       dofids[k] = (*idis)[interf]->Dof(node)[k];

      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(new MORTAR::MortarNode(node->Id(),
         node->X(), node->Owner(), dofpernode_, dofids, false));

      interface->AddMortarNode(mrtrnode);
    }

    // feeding master elements to the interface
    for(int lid=0; lid<(*idis)[interf]->NumMyColElements(); ++lid)
    {
     DRT::Element* currele = (*idis)[interf]->lColElement(lid);
     Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(new MORTAR::MortarElement(currele->Id(),
         currele->Owner(), currele->Shape(), currele->NumNode(), currele->NodeIds(), false));

     interface->AddMortarElement(mrtrele);
    }

    mortarinterface_.push_back(interface);
  }

  return;
}


/*----------------------------------------------------------------------*
 | useful parameters for mortar coupling                    ghamm 06/14 |
 *----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::MortarParams
(
  Teuchos::ParameterList& mortarparams
)
{
  // set options for mortar coupling
  mortarparams.set<std::string>("SEARCH_ALGORITHM","Binarytree");
  mortarparams.set<double>("SEARCH_PARAM", 0.1);
  mortarparams.set<std::string>("SEARCH_USE_AUX_POS", "no");
  mortarparams.set<std::string>("PARALLEL_REDIST","no");
  mortarparams.set<std::string>("LM_SHAPEFCN","dual");
  // master side is given to the mortar interface once in the beginning with full overlap
  mortarparams.set<std::string>("REDUNDANT_STORAGE","Master");
  // for triangles which will all be lying in a quad element on the master side, this should be a good set of params
  mortarparams.set<std::string>("INTTYPE","Elements");
  // constant load per fluid element assumed --> p0 element for mortar coupling only needs 3 gp per tri
  mortarparams.set<int>("NUMGP_PER_DIM",1);
  mortarparams.set<bool>("NURBS",false);
  // elements with decoupled geometry and node information are used
  mortarparams.set<bool>("GEO_DECOUPLED", true);

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

  // feeding slave nodes and elements (tris) to the interface
  // Note: offset is needed for slave dofs and slave node ids
  double nodalpos[dim_];
  std::vector<int> dofids(dofpernode_);
  // 3 points and 1 node per tri need unique ids
  int nodepointids[4];
  const int nodepointoffset = 4*nodeoffset/3;  // point/node offset from proc to proc
  const int nodeoffsetthird = nodeoffset/3;  // ele offset from proc to proc
  // three points and one node is necessary for each slave element
  const int numtris = aerocoords.size()/3;
  for (int eleiter = 0; eleiter<numtris; ++eleiter)
  {
    // add single node for the dof information of the mortar element

    for (int k=0; k<dim_; ++k)
      nodalpos[k] = 0.0; // dummy values here
    for (int k=0; k<dofpernode_; ++k)
      dofids[k] = mastermaxdof_ + (eleiter*4+0+nodepointoffset)*dofpernode_ + k;
    nodepointids[0] = eleiter*4+0+nodepointoffset+mastermaxnodeid_;

    Teuchos::RCP<MORTAR::MortarNode> mrtrnode =
        Teuchos::rcp(new MORTAR::MortarNode(nodepointids[0], nodalpos, myrank_, dofpernode_, dofids, true));

    mortarinterface_[interf]->AddMortarNode(mrtrnode);

    // insert three points for the geometry of the mortar element

    for (int nodeiter=1; nodeiter<4; ++nodeiter)
    {
      for (int k=0; k<dim_; ++k)
        nodalpos[k] = aerocoords[eleiter*3+(nodeiter-1)](k);
      for (int k=0; k<dofpernode_; ++k)
        dofids[k] = mastermaxdof_ + (eleiter*4+nodeiter+nodepointoffset)*dofpernode_ + k;
      nodepointids[nodeiter] = eleiter*4+nodeiter+nodepointoffset+mastermaxnodeid_;

      Teuchos::RCP<MORTAR::MortarNode> mrtrnode =
          Teuchos::rcp(new MORTAR::MortarNode(nodepointids[nodeiter], nodalpos, myrank_, dofpernode_, dofids, true));

      mortarinterface_[interf]->AddMortarPoint(mrtrnode);
    }

    Teuchos::RCP<MORTAR::MortarElementGeoDecoupl> mrtrele =
        Teuchos::rcp(new MORTAR::MortarElementGeoDecoupl(eleiter + nodeoffsetthird + mastermaxeleid_,
        myrank_, DRT::Element::tri3, 1, &nodepointids[0], 3, &nodepointids[1], true));

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

  if(mechanicalcoupling_ == true)
  {
    // insert current displacement status of master side; slave side is placed correctly
    for(int i=0; i<idispn->MyLength(); ++i)
    {
      int gid = idispn->Map().GID(i);
      int lid = mergeddispn->Map().LID(gid);
      (*mergeddispn)[lid] = (*idispn)[i];
    }
  }

  // set displacement state in mortar interface
  mortarinterface_[interf]->SetState(MORTAR::state_new_displacement, *mergeddispn);

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
//  if(myrank_==0)
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
  int offset = 0;
  if (mechanicalcoupling_ == false && thermalcoupling_ == true)
    offset = dim_;

  Teuchos::RCP<Epetra_Vector> iload_fl = LINALG::CreateVector(DinvM_[interf]->RowMap(),true);

  // NOTE: There is no dofset available for artificially created fluid nodes
  // Hence: Order is important and values from aeroforces can directly be inserted into interface load vector
  for(size_t eleiter=0; eleiter<aeroforces.size(); ++eleiter)
  {
    // currload is always filled with three mechanical dofs and one thermal dof
    // which is assumed constant within one mortar coupling element
    LINALG::Matrix<4,1>& currload = aeroforces[eleiter];
    for(int k=0; k<dofpernode_; ++k)
    {
      (*iload_fl)[eleiter*dofpernode_ + k] = currload(k+offset);
    }
  }

  // compute corresponding structural values using mortar
  // return type is Epetra_MultiVector because Epetra_FEVector is used in SlaveToMaster
  Teuchos::RCP<Epetra_MultiVector> iload_solid = SlaveToMaster(interf, iload_fl);

  const int numnodes = iload_solid->MyLength()/dofpernode_;
  // split data in iload_solid into structural and thermal load (4 values per node)
  if (mechanicalcoupling_ == true && thermalcoupling_ == true)
  {
    for(int nodeiter=0; nodeiter<numnodes; ++nodeiter)
    {
      // structural load
      for(int k=0; k<dim_; ++k)
      {
        const int lid = nodeiter*dofpernode_ + k;
        const double val = (*(*iload_solid)(0))[lid];
        const int gid = iload_solid->Map().GID(lid);
        (*iforce)[ iforce->Map().LID(gid) ] = val;
      }
      // thermal load
      const int lid = nodeiter*dofpernode_ + dim_;
      const double val = (*(*iload_solid)(0))[lid];
      const int gid = iload_solid->Map().GID(lid);
      (*ithermoload)[ ithermoload->Map().LID(gid) ] = val;
    }
  }
  else if(mechanicalcoupling_ == true && thermalcoupling_ == false)
  {
    for(int nodeiter=0; nodeiter<numnodes; ++nodeiter)
    {
      // structural load
      for(int k=0; k<dim_; ++k)
      {
        const int lid = nodeiter*dofpernode_ + k;
        const double val = (*(*iload_solid)(0))[lid];
        const int gid = iload_solid->Map().GID(lid);
        (*iforce)[ iforce->Map().LID(gid) ] = val;
      }
    }
  }
  else if(mechanicalcoupling_ == false && thermalcoupling_ == true)
  {
    for(int nodeiter=0; nodeiter<numnodes; ++nodeiter)
    {
      // thermal load
      const int lid = nodeiter;
      const double val = (*(*iload_solid)(0))[lid];
      const int gid = iload_solid->Map().GID(lid);
      (*ithermoload)[ ithermoload->Map().LID(gid) ] = val;
    }
  }
  else
    dserror("coupling type not available");

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


/*----------------------------------------------------------------------*
 | pack data for quad4 elements on Baci side                ghamm 12/11 |
 *----------------------------------------------------------------------*/
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
  const int quadnodes = 4;
  const double invlengthscaling = 1.0 / LengthScaling();
  Teuchos::RCP<Epetra_Vector> proc0disp = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> proc0vel = Teuchos::null;
  if(mechanicalcoupling_ == true)
  {
    // make all interface data available on proc 0
    // idisp and ivel
    proc0disp = Teuchos::rcp(new Epetra_Vector(istructproc0importer_[interf]->TargetMap()));
    int err = proc0disp->Import(*idisp,*istructproc0importer_[interf],Insert);
    if (err>0)
      dserror("Importing everything to proc 0 went wrong. Import returns %d",err);
    proc0vel = Teuchos::rcp(new Epetra_Vector(istructproc0importer_[interf]->TargetMap()));
    err = proc0vel->Import(*ivel,*istructproc0importer_[interf],Insert);
    if (err>0)
      dserror("Importing everything to proc 0 went wrong. Import returns %d",err);
  }

  Teuchos::RCP<Epetra_Vector> proc0temp = Teuchos::null;
  if(thermalcoupling_ == true)
  {
    // itemp
    proc0temp = Teuchos::rcp(new Epetra_Vector(ithermoproc0importer_[interf]->TargetMap()));
    int err = proc0temp->Import(*itemp,*ithermoproc0importer_[interf],Insert);
    if (err>0)
      dserror("Importing everything to proc 0 went wrong. Import returns %d",err);
  }

  // gather data on TSI side on proc 0
  if(myrank_ == 0)
  {
    int numcolele = -1;
    if(mechanicalcoupling_ == true)
    {
      if(istructdis_[interf]->lColElement(0)->NumNode() != quadnodes)
        dserror("current implementation only for quad4 elements in the coupling interface");

      numcolele = istructdis_[interf]->NumMyColElements();
    }
    else
    {
      if(ithermodis_[interf]->lColElement(0)->NumNode() != quadnodes)
        dserror("current implementation only for quad4 elements in the coupling interface");

      numcolele = ithermodis_[interf]->NumMyColElements();
    }

    packeddata.resize(numcolele*35);
    for(int lid=0; lid<numcolele; ++lid)
    {
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      std::vector<double> eledisp(quadnodes*dim_, 0.0);
      std::vector<double> elevel(quadnodes*dim_, 0.0);
      DRT::Element* ele = NULL;
      if(mechanicalcoupling_ == true)
      {
        ele = istructdis_[interf]->lColElement(lid);
        // get element location vector of structural problem
        ele->LocationVector(*istructdis_[interf],lm,lmowner,lmstride);

        // get element data
        // disp
        eledisp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*proc0disp, eledisp ,lm);
        // vel
        elevel.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*proc0vel, elevel ,lm);
      }

      std::vector<double> eletemp(quadnodes, 0.0);
      if(thermalcoupling_ == true)
      {
        // get element location vector of thermal problem
        ele = ithermodis_[interf]->lColElement(lid);
        ele->LocationVector(*ithermodis_[interf],lm,lmowner,lmstride);

        // get element data
        // temp
        eletemp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*proc0temp, eletemp ,lm);
      }

      // fill in position, velocity and temperature for quad4 (fifth point is center point)
      // x1y1z1 x2y2z2 x3y3z3 x4y4z4 x5y5z5  u1v1w1 u2v2w2 u3v3w3 u4v4w4 u5v5w5  T1T2T3T4T5
      const int eleoffset = lid*35;
      std::vector<double> xcenter(dim_, 0.0);
      std::vector<double> velcenter(dim_, 0.0);
      double Tcenter = 0.0;
      // loop over nodes of element
      for(int k=0; k<quadnodes; ++k)
      {
        DRT::Node** nodes = ele->Nodes();

        for(int d=0; d<dim_; ++d)
        {
          // pos
          double pos = (nodes[k]->X()[d] + eledisp[k*3 + d]) * invlengthscaling;
          packeddata[eleoffset + k*dim_ + d] = pos;
          xcenter[d] += pos*0.25;
          // vel
          packeddata[eleoffset + (k+5)*3 + d] = (elevel[k*3 + d]) * invlengthscaling * TimeScaling();
          velcenter[d] += elevel[k*dim_ + d]*0.25;
        }
        // temp
        packeddata[eleoffset + 10*dim_ + k] = eletemp[k];
        Tcenter += eletemp[k]*0.25;
      }
      // center point is added
      for(int d=0; d<dim_; ++d)
      {
        // pos
        packeddata[eleoffset + 4*dim_ + d] = xcenter[d];
        // vel
        packeddata[eleoffset + 9*dim_ + d] = velcenter[d];
      }
      // temp
      packeddata[eleoffset + 10*dim_ + 4] = Tcenter;
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
  Teuchos::RCP<const Epetra_Vector> fullvector
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
  if(iforce != Teuchos::null and force != Teuchos::null)
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
  if(iforce != Teuchos::null and force != Teuchos::null)
    thermorowmapext_[interf]->InsertCondVector(iforce, force);
  return;
}
