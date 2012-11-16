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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::UTILS::AeroCouplingUtils::AeroCouplingUtils
(
  Teuchos::RCP<DRT::Discretization> structdis,
  Teuchos::RCP<DRT::Discretization> thermodis
)
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
  if(maxid_ > 0)
    dserror("currently only 1 coupling surface with the Aero-Code is allowed!");

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
  RCP<Epetra_Comm> com = Teuchos::rcp(structdis->Comm().Clone());
  // structure
  const string discret_name1 = "istructdis";
  Teuchos::RCP<DRT::Discretization> istructnewdis = Teuchos::rcp(new DRT::Discretization(discret_name1,com));
  const int myrank = istructnewdis->Comm().MyPID();
  // thermo
  const string discret_name2 = "ithermodis";
  Teuchos::RCP<DRT::Discretization> ithermonewdis = Teuchos::rcp(new DRT::Discretization(discret_name2,com));

  
  std::map<int, std::map<int, DRT::Node*> >::iterator mnit;
  std::map<int, DRT::Node* >::iterator nit;

  std::vector<int> nodeids;
  //nodes filled in new interface discretizations
  for ( mnit=structgnodes.begin(); mnit != structgnodes.end(); mnit++ )
  {
    for ( nit=mnit->second.begin(); nit != mnit->second.end(); nit++ )
    {
      //ids for row map are gathered
      if ((*nit).second->Owner() == myrank)
        nodeids.push_back( (*nit).second->Id() );

      //discretizations are feeded
      istructnewdis->AddNode(Teuchos::rcp(new DRT::Node((*nit).second->Id(), (*nit).second->X(), (*nit).second->Owner())));
      ithermonewdis->AddNode(Teuchos::rcp(new DRT::Node((*nit).second->Id(), (*nit).second->X(), (*nit).second->Owner())));
    }
  }
  //row node map of struct interface
  Teuchos::RCP<Epetra_Map> istructnoderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,istructnewdis->Comm()));
  //fully overlapping node map
  Teuchos::RCP<Epetra_Map> istructrednodecolmap =  LINALG::AllreduceEMap(*istructnoderowmap);

  //row node map of thermo interface
  Teuchos::RCP<Epetra_Map> ithermonoderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,istructnewdis->Comm()));
  //fully overlapping node map
  Teuchos::RCP<Epetra_Map> ithermorednodecolmap =  LINALG::AllreduceEMap(*ithermonoderowmap);


  //care about eles in the interface
  std::vector<int> eleids;
  std::map<int, std::vector<int> > interfaceeleids;
  //eles filled in new interface discretizations
  for ( meit=structelements.begin(); meit != structelements.end(); meit++ )
  {
    for ( eit=meit->second.begin(); eit != meit->second.end(); eit++ )
    {
      if ((*eit).second->Owner() == myrank)
      {
        // each interface stores its row elements in the interface
        interfaceeleids[((*meit).first)].push_back( (*eit).second->Id() );
        // all ele ids are gathered for the discretisation
        eleids.push_back( (*eit).second->Id() );

        // structural surface elements cannot be distributed --> Bele3 element is used
        Teuchos::RCP<DRT::Element> istructele = DRT::UTILS::Factory("BELE3","Polynomial", (*eit).second->Id(), (*eit).second->Owner());
        istructele->SetNodeIds( ((*eit).second->NumNode()), ((*eit).second->NodeIds()) );
        istructnewdis->AddElement( istructele );
        // thermo interface elements
#ifdef D_THERMO
        Teuchos::RCP<DRT::Element> ithermoele = DRT::UTILS::Factory("THERMO","Polynomial", (*eit).second->Id(), (*eit).second->Owner());
        ithermoele->SetNodeIds( ((*eit).second->NumNode()), ((*eit).second->NodeIds()) );
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Thermo>(ithermoele,true)->SetDisType( istructele->Shape() );
        ithermonewdis->AddElement( ithermoele );
#else
        dserror("D_THERMO flag is switched off! You need thermo for Aero_TFSI.");
#endif
      }
    }
  }

  //row ele map of structural interface
  Teuchos::RCP<Epetra_Map> istructelerowmap = Teuchos::rcp(new Epetra_Map(-1,eleids.size(),&eleids[0],0,istructnewdis->Comm()));
  //fully overlapping ele map
  Teuchos::RCP<Epetra_Map> istructredelecolmap =  LINALG::AllreduceEMap(*istructelerowmap);

  //row ele map of thermo interface
  Teuchos::RCP<Epetra_Map> ithermoelerowmap = Teuchos::rcp(new Epetra_Map(-1,eleids.size(),&eleids[0],0,istructnewdis->Comm()));
  //fully overlapping ele map
  Teuchos::RCP<Epetra_Map> ithermoredelecolmap =  LINALG::AllreduceEMap(*ithermoelerowmap);

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
  istructdis_ = istructnewdis;
  ithermodis_ = ithermonewdis;
  //storage of useful dof maps
  istructdofredumap_ =  LINALG::AllreduceEMap(*(istructdis_->DofRowMap()));
  ithermodofredumap_ =  LINALG::AllreduceEMap(*(ithermodis_->DofRowMap()));

  //storage of the redundant structural interface elements for each single interface
  std::map<int, Teuchos::RCP<Epetra_Map> > singleinterfstructelerowmap;
  std::map<int, Teuchos::RCP<Epetra_Map> > singleinterfstructredelecolmap;
  for(int interf=0; interf<(maxid_+1); interf++)
  {
    //row ele map of each single structural interface
    singleinterfstructelerowmap[interf] = Teuchos::rcp(new Epetra_Map(-1,interfaceeleids[interf].size(),
                          &interfaceeleids[interf][0],0,istructnewdis->Comm()));
    //fully overlapping ele map of each single structural interface
    singleinterfstructredelecolmap[interf] =  LINALG::AllreduceEMap(*(singleinterfstructelerowmap[interf]));

    int nummycoleles = singleinterfstructredelecolmap[interf]->NumMyElements();
    for (int it=0; it<nummycoleles; it++)
    {
      // first, find elegid in map of each single interface, then find lid in istructdis_
      int elegid = singleinterfstructredelecolmap[interf]->GID(it);
      int elelid = istructredelecolmap->LID(elegid);
      structreduelements_[interf][it] = Teuchos::rcp(istructdis_->lColElement(elelid), false);
    }
  }

  // initial search radius in search tree
  maxmindist_ = 1.0e-1;

  // epetra_vector with material configuration of row nodes is created
  matconfig_ = LINALG::CreateVector(*(istructdis_->DofRowMap()), true);
  for(int myl=0; myl<istructnoderowmap->NumMyElements(); myl++)
  {
    int gid = istructnoderowmap->GID(myl);
    const DRT::Node* node = istructdis_->gNode(gid);

    std::vector<int> lm;
    lm.reserve(3);
    // extract global dof ids
    istructdis_->Dof(node, lm);

    for (int a=0; a<3; a++)
    {
      int gid = lm[a];
      double val = node->X()[a];
      matconfig_->ReplaceGlobalValues(1, &val, &gid);
    }
  }
  
  // map extractor for the communication between the full thermo field and the newly created thermo surf discret
  thermorowmapext_ = Teuchos::rcp(new LINALG::MapExtractor(*(thermodis->DofRowMap()),Teuchos::rcp(ithermodis_->DofRowMap(), false)));

  // set-up FSI interface
  interface_ = Teuchos::rcp(new STR::AUX::MapExtractor);
  interface_->Setup(*structdis, *structdis->DofRowMap());

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int,LINALG::Matrix<3,1> > FS3I::UTILS::AeroCouplingUtils::CurrentInterfacePos
(
  Teuchos::RCP<Epetra_Vector> redudisp
)
{
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleiter;

  // multiple interfaces
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::const_iterator interfaceiter;

  //loop over all elements in the interface due to fully redundant storage
  for (interfaceiter = structreduelements_.begin(); interfaceiter != structreduelements_.end(); interfaceiter++)
  {
    for(eleiter=interfaceiter->second.begin(); eleiter!=interfaceiter->second.end(); eleiter++)
    {
      Teuchos::RCP<DRT::Element> tmpele = eleiter->second;

      const int* n = tmpele->NodeIds();

      // fill currentpositions
      for (int j=0; j < tmpele->NumNode(); j++)
      {
        const int gid = n[j];
        const DRT::Node* node = istructdis_->gNode(gid);
        std::vector<int> lm;
        lm.reserve(3);
        // extract global dof ids
        istructdis_->Dof(node, lm);
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
  }

    return currentpositions;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::ProjectForceOnStruct
(
  Teuchos::RCP<Epetra_Vector> idispn,
  std::map<int, std::map<int, LINALG::Matrix<3,1> > >& aerocoords,
  std::map<int, std::map<int, LINALG::Matrix<4,1> > >& aeroforces,
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  // Redistribute displacement of structural nodes on the interface to all processors --> redundant disp needed
  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_,*(istructdis_->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_,true);
  redudisp -> Import(*idispn,*interimpo,Add);

  //current position of structural nodes for the search tree (always 3 coordinates)
  std::map<int,LINALG::Matrix<3,1> > currentpositions = CurrentInterfacePos(redudisp);

  //projection of the aero forces onto the closest element of each structural interface
  for(int interf=0; interf<(maxid_+1); interf++)
  {
    std::map<int, LINALG::Matrix<3,1> >::iterator forcesiter;
    for(forcesiter = aerocoords[interf].begin(); forcesiter != aerocoords[interf].end(); ++forcesiter)
    {
      //init of 3D search tree
      Teuchos::RCP<GEO::SearchTree> searchTree = Teuchos::rcp(new GEO::SearchTree(0));
      const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofEles(structreduelements_[interf], currentpositions);
      searchTree->initializeTreeSlideALE(rootBox, structreduelements_[interf], GEO::TreeType(GEO::OCTTREE));

      //search for near elements to the aero coords
      std::map<int,std::set<int> >  closeeles =
          searchTree->searchElementsInRadius(*istructdis_,currentpositions,forcesiter->second,maxmindist_,0);
      //if no close elements could be found, try with a much larger radius and print a warning
      if (closeeles.empty())
      {
        cout<<"WARNING: no elements found in radius r="<<maxmindist_<<". Will try once with a bigger radius!"<<endl;
        closeeles = searchTree->searchElementsInRadius(*istructdis_,currentpositions,forcesiter->second,100.0*maxmindist_,0);
        maxmindist_ *= 10.0;

        // if still no element is found, complain about it!
        if (closeeles.empty())
          dserror("No elements in a large radius! Should not happen!");
      }

      //search for the closest object, more exactly it's coordinates and the corresponding surface id
      LINALG::Matrix<3,1> minDistCoords;
      int eleid = GEO::nearest3DObjectInNode(istructdis_, structreduelements_[interf], currentpositions,
          closeeles, forcesiter->second, minDistCoords);
      if(eleid == -1)
        dserror("Surface id couldn't be found. Weird! Hence, a distribution of the aero force is not performed!");

      //corresponding forces (from the aero code) to the coords that are about to be distributed
      LINALG::Matrix<4,1> force = aeroforces[interf][forcesiter->first];

      //distribution of forces takes place in here
      DistributeForceToNodes(eleid, interf, currentpositions, minDistCoords, force, iforce, ithermoload);
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::UTILS::AeroCouplingUtils::DistributeForceToNodes
(
  int eleid,
  int interface,
  std::map<int,LINALG::Matrix<3,1> >& currentpositions,
  LINALG::Matrix<3,1>& projpoint,
  LINALG::Matrix<4,1>& aeroforce,
  Teuchos::RCP<Epetra_Vector> distributedforces,
  Teuchos::RCP<Epetra_Vector> ithermoload
)
{
  Teuchos::RCP<DRT::Element> element = structreduelements_[interface][eleid];

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
    istructdis_->Dof(snode, slm);
    // aeroforce[0...2] are forces
    for(int k=0; k<num; k++)
    {
      val[k] = funct[iter] * aeroforce(k);
      lmowner[k] = snode->Owner();
    }

    // do assembly of structural forces
    LINALG::Assemble(*distributedforces,val,slm,lmowner);

    // secondly, care about heat fluxes
    DRT::Node* tnode = ithermodis_->gNode(snode->Id());

    // prepare assembly
    num = 1;
    Epetra_SerialDenseVector tval(num);
    std::vector<int> tlm;
    tlm.reserve(num);
    std::vector<int> tlmowner(num);
    // extract global dof ids
    ithermodis_->Dof(tnode, tlm);
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
void FS3I::UTILS::AeroCouplingUtils::PackData
(
  Teuchos::RCP<Epetra_Vector> idispnp,
  Teuchos::RCP<Epetra_Vector> itempnp,
  std::vector<double>& packeddata
)
{
  // make all interface data full redundant
  // idispn
  Teuchos::RCP<Epetra_Import> interimpo = Teuchos::rcp(new Epetra_Import(*istructdofredumap_,*(istructdis_->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redudisp = LINALG::CreateVector(*istructdofredumap_,true);
  redudisp -> Import(*idispnp,*interimpo,Add);

  // itemp
  Teuchos::RCP<Epetra_Import> interimpothermo = Teuchos::rcp(new Epetra_Import(*ithermodofredumap_,*(ithermodis_->DofRowMap())));
  Teuchos::RCP<Epetra_Vector> redutemp = LINALG::CreateVector(*ithermodofredumap_,true);
  redutemp -> Import(*itempnp,*interimpothermo,Add);

  //current position of structural nodes
  std::map<int,LINALG::Matrix<3,1> > currentpositions = CurrentInterfacePos(redudisp);

  int dim = 3;

  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleiter;
  for(eleiter=structreduelements_[0].begin(); eleiter!=structreduelements_[0].end(); eleiter++)
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
      DRT::Node* tnode = ithermodis_->gNode(nodalgid);
      std::vector<int> lmt;
      lmt.reserve(1);
      ithermodis_->Dof(tnode, lmt);

      int lid = ithermodofredumap_->LID(lmt[0]);

      midtemp[0] += (0.25 * (redutemp->operator[](lid)));

    } // end of loop over four nodes

    // add the gathered data
    for(int i=0; i<dim; i++)
    {
      packeddata.push_back(midposition[i]);
    }
    packeddata.push_back(midtemp[0]);

  } // end elements loop

  return;
}


/// extract structural interface values for a given full field
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::StrExtractInterfaceVal
(
  Teuchos::RCP<Epetra_Vector> fullvector)
{
  return interface_->ExtractFSICondVector(fullvector);
}


/// extract thermal interface values for a given full field
Teuchos::RCP<Epetra_Vector> FS3I::UTILS::AeroCouplingUtils::ThrExtractInterfaceVal
(
  Teuchos::RCP<Epetra_Vector> fullvector)
{
  return thermorowmapext_->ExtractCondVector(fullvector);
}


/// apply interface tractions from the interface to a full field
void FS3I::UTILS::AeroCouplingUtils::StrApplyInterfaceVal
(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> force
)
{
  interface_->AddFSICondVector(iforce, force);
  return;
}


/// apply interface heat flux from the interface to a full field
void FS3I::UTILS::AeroCouplingUtils::ThrApplyInterfaceVal
(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> force
)
{
  thermorowmapext_->InsertCondVector(iforce, force);
  return;
}
