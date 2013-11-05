/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsi_crack.cpp

\brief Adapter Layer for FSI with cracking structure

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15267
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
//#include "../drt_structure/strtimint_create.H"
#include "ad_str_fsi_crack.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_factory.H"
//#include "../linalg/linalg_utils.H"
//#include "../drt_constraint/constraint_manager.H"
//#include "../drt_structure/stru_aux.H"

#include "../drt_crack/InsertCohesiveElements.H"
#include "../drt_crack/dcohesive.H"

/*======================================================================*/
/* constructor */
ADAPTER::FSICrackingStructure::FSICrackingStructure
(
  Teuchos::RCP<Structure> stru
)
:FSIStructureWrapper(stru)
{

  allDone_ = false;

  // make sure
  if (structure_ == Teuchos::null)
    dserror("Failed to create the underlying structural adapter");

  // TODO: At the moment, we insert the cohesive elements before building the
  // structural time integration. This is because once we build timint, dofrowmap
  // of structural discretization is lost, and there is no possibility of building it again.
  // After inserting all cohesive elements, we loop over each element in discret to get
  // cohesive elements and build all other details here.
  // This is not possible in case, if we want to add cohesive elemets during the simulation.
  // TO-DO is to find a suitable way of including nodes and elements into structural discretization
  // after time integration is built

  {
    structdis_ = DRT::Problem::Instance()->GetDis("structure");

    /*******************************************************/
    const int totEle = structdis_->NumMyRowElements();

    for (int i=0; i<totEle; i++)
    {
      DRT::Element* Ele = structdis_->lRowElement(i);

      switch( Ele->Shape() )
      {
      case DRT::Element::line2:
      {
        const int* nodes = Ele->NodeIds();
        coheEleMasSla_[Ele->Id()] = std::make_pair( nodes[0], nodes[1] );

        break;
      }
      default:
        break;
      }
    }
    /*******************************************************/
  }

  // construct master and slave crack nodes discretization
  std::vector<std::string> conditions_to_copy;
  std::string element_name = "BELE3";
  masterCrackDis_ = DRT::UTILS::CreateDiscretizationFromCondition(structdis_, "masterCrackSurface", "boundary",
                                                                                                  element_name, conditions_to_copy);
  slaveCrackDis_ = DRT::UTILS::CreateDiscretizationFromCondition(structdis_, "slaveCrackSurface", "boundary",
                                                                                                  element_name, conditions_to_copy);
  if (masterCrackDis_->NumGlobalNodes() == 0)
  {
    dserror("Empty master crack discretization detected\n");
  }

  if (slaveCrackDis_->NumGlobalNodes() == 0)
  {
    dserror("Empty slave crack discretization detected\n");
  }

  masterCrackDis_->FillComplete();
  slaveCrackDis_->FillComplete();

  // compute all the possible master-slave combination of crack surfaces
  // build the structure "possSurfaces_" which is the basic framework of
  // finding the crack surfaces that are to be added
  const int numMasSur = masterCrackDis_->NumMyRowElements();
  const int numSlaSur = slaveCrackDis_->NumMyRowElements();

  //structdis_->Print(std::cout);//blockkk
  //std::cout<<"--------printing Master crack discretization---------\n";//blockkk
  //masterCrackDis_->Print(std::cout);//blockkk
  //std::cout<<"--------printing Slave crack discretization---------\n";//blockkk
  //slaveCrackDis_->Print(std::cout);//blockkk

  //dserror("done");//blockkk

  for ( int i=0; i<numMasSur; i++ )
  {
    DRT::Element* mSurEle = masterCrackDis_->lRowElement(i);

    // ------------------------------------------------------------------------------------
    // STEP 1 : find all dcohesive elements connected to this surface
    // ------------------------------------------------------------------------------------
    const int* mEleNodes = mSurEle->NodeIds();

    //std::cout<<"masster nodes = "<<mEleNodes[0]<<"\t"<<mEleNodes[1]<<"\t"<<mEleNodes[2]<<"\t"<<mEleNodes[3]<<"\n";;//blockkk
    //dserror("done");//blockkk

    // will store all cohesive elements associated with this master crack surface element
    std::vector<int> dcohMasEle;

    // for this master surface element, the corresponding slave element should have these nodes
    std::vector<int> reqdSlaNod;
    std::vector<int> tempMas;

    for( std::map<int, std::pair<int,int> >::iterator cohit = coheEleMasSla_.begin();
                                                      cohit != coheEleMasSla_.end(); cohit ++ )
    {
      std::pair<int,int> pai = cohit->second;
      int mas_coh = pai.first;

      for( int mnod=0; mnod<mSurEle->NumNode(); mnod++ )
      {
        if( mas_coh == mEleNodes[mnod] )
        {
          tempMas.push_back( mas_coh );
          dcohMasEle.push_back( cohit->first );
          reqdSlaNod.push_back( pai.second );

          break;
        }
      }
    }

    /*--------------------------------------------------------------*/ //blockkk
    /*std::cout<<"found nodes\n";
    for( int ik=0;ik<tempMas.size();ik++ )
      std::cout<<tempMas[ik]<<"\t";
    std::cout<<"\n";*/
    /*--------------------------------------------------------------*/

    // For elements that share Dirichlet boundary conditions,
    // all nodes of master (or slave) crack has a dcohesive elements
    // attached to it
    if( mSurEle->NumNode() != (int)reqdSlaNod.size() )
    {
      for( int mnod=0; mnod<mSurEle->NumNode(); mnod++ )
      {
        if( std::find( tempMas.begin(), tempMas.end(), mEleNodes[mnod] ) == tempMas.end() )
          reqdSlaNod.push_back( mEleNodes[mnod] );
      }
    }

    //std::cout<<"reqdSlaNod = "<<reqdSlaNod[0]<<"\t"<<reqdSlaNod[1]<<"\t"<<reqdSlaNod[2]<<"\t"<<reqdSlaNod[3]<<"\n";//blockkk

    if( dcohMasEle.size() == 0 )
      dserror("not a single dcohesive element attached with this element\n");

    crackSurface_ crs;
    crs.setMaster( mSurEle );
    crs.setProcInfo( false );
    crs.setAttachedCohesiveElem( dcohMasEle );

    // ------------------------------------------------------------------------------------
    // STEP 2 : find the corresponding slave side that is connected to this master side
    //          via the dcohesive spring elements. If two sides share common dcohesive
    //          elements, then they form a master-slave combination
    // ------------------------------------------------------------------------------------

    // sort so that comparison is easier
    std::sort( reqdSlaNod.begin(), reqdSlaNod.end() );

    /******************************************/ //blockkk
    /*std::cout<<"required slave nodes = ";
    for( unsigned imm=0;imm<reqdSlaNod.size();imm++ )
      std::cout<<reqdSlaNod[imm]<<"\t";
    std::cout<<"\n";*/
    /******************************************/ //blockkk

    bool foundSlave = false;
    for( int j=0; j<numSlaSur; j++ )
    {
      DRT::Element* sSurEle = slaveCrackDis_->lRowElement(j);
      const int* sEleNodes = sSurEle->NodeIds();

      // only tri-tri or quad-quad combination can be master-slave
      if( mSurEle->NumNode() != sSurEle->NumNode() )
        continue;

      std::vector<int> sNodes;
      for( int snod=0; snod<sSurEle->NumNode(); snod++ )
        sNodes.push_back( sEleNodes[snod] );

      // it is sorted so that we can compare whether any slave side has same nodes
      std::sort( sNodes.begin(), sNodes.end() );

      /******************************************/ //blockkk
      /*std::cout<<"found slave nodes = ";
      for( unsigned imm=0;imm<sNodes.size();imm++ )
        std::cout<<sNodes[imm]<<"\t";
      std::cout<<"\n";*/
      /******************************************/ //blockkk

      if( reqdSlaNod == sNodes )
      {
        foundSlave = true;
      }

      /*//
      // in some cases for example at the end of structure boundary
      // two nodes can be joined for application of Dirichlet boundary condition
      //          \     /
      //           \   /
      //            \ /
      //             + ----> Dirichlet boundary node
      //
      if( reqdSlaNod.size() < sNodes.size() )
      {
        unsigned nofound = 0;
        for( unsigned reqid=0; reqid < reqdSlaNod.size(); reqid++ )
        {
          if( std::find( sNodes.begin(), sNodes.end(), reqdSlaNod[reqid] ) != sNodes.end() )
          {
            nofound++;
          }
        }

        if( nofound == reqdSlaNod.size() )
          foundSlave = true;
      }*/

      if( foundSlave )
      {
        crs.setSlave( sSurEle );
        break;
      }
    }

    if( not foundSlave )
      dserror( "slave element not found for this master element\n" );

    possSurfaces_[mSurEle] = crs;
  //  dserror("atleast one done");
  }

  /******************************************************************///blockkk
  /*std::cout<<"number of poss surf = "<<possSurfaces_.size()<<"\n";//blockkk
  for( std::map<DRT::Element*,crackSurface_>::iterator ip=possSurfaces_.begin();ip!=possSurfaces_.end();ip++ )
  {
    crackSurface_ crs1 = ip->second;
    std::cout<<"number of cohesive elements = "<<crs1.dcoh_.size()<<"\n";
    std::cout<<"they are = "<<crs1.dcoh_[0]<<"\t"<<crs1.dcoh_[1]<<"\n";
    std::cout<<"master element id = "<<crs1.mas_->Id()<<"\n";
    std::cout<<"number of nodes = "<<crs1.mas_->NumNode()<<"\n";
    crs1.mas_->Print(std::cout);
    crs1.sla_->Print(std::cout);
  }
  dserror("done printing\n");//blockkk*/
  /******************************************************************///blockkk

  /********************/ //blockkk
  /*for( std::map<int, std::pair<int,int> >::iterator cohit = coheEleMasSla_.begin();
                                                        cohit != coheEleMasSla_.end(); cohit ++ )
  {
    std::cout<<"ele id = "<<cohit->first<<" mas sla = "<<cohit->second.first<<" "<<cohit->second.second<<"\n";
  }
  dserror("done");*/
  /********************/ //blockkk
}

/*-------------------------------------------------------------------------------------------*
 *  check whether all possible crack surfaces are already added to the cutsides       sudhakar 08/13
 *-------------------------------------------------------------------------------------------*/
void ADAPTER::FSICrackingStructure::checkAllDone()
{
#ifdef DEBUG
  if( allDone_ )
    dserror("All dcohesive elements are added in the last step itself. How is this called again?\n");
#endif

  for( std::map<DRT::Element*,crackSurface_>::iterator it=possSurfaces_.begin();
                                                       it!=possSurfaces_.end();it++ )
  {
    bool ispro = it->second.isDone_;
    if( not ispro )
      return;
  }

  allDone_ = true;


}

/*void ADAPTER::FSICrackingStructure::addCrackSurfacesToCutSides( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                                std::map<int, Teuchos::RCP<DRT::Element> >& tipele )
{
  if( allDone_ )
  {
    return;
  }

  std::vector<DRT::Element*> crk_sur;           // crack surface elements to be added to cut sides
  int tipeleno = 0;                             // counter to see how many tip elements are added

  std::cout<<"----------------------I am inserting new crack surface elemts--------------------\n";//blockkk

  std::cout<<"number of crack surfaces = "<<possSurfaces_.size()<<"\n";//blockkk

  for( std::map<DRT::Element*,crackSurface_>::iterator it = possSurfaces_.begin();
                                                       it != possSurfaces_.end(); it++ )
  {
    crackSurface_& crs1 = it->second;

    // this master-slave combination is already processed
    if( crs1.isDone_ )
      continue;

    // element ids of all cohesive elements attached with this surface
    std::vector<int> cohele = crs1.dcoh_;

    int failno = 0;             // no of failed cohesive elements
    bool allfailed = true;      // if all cohesive elements failed
    std::vector<int> workele;   // Ids of cohesive elements that are still working

    //structdis_->Print(std::cout);//blockkk

    for( unsigned cohno=0; cohno<cohele.size();cohno++ )
    {
      int cohid = cohele[cohno];
      DRT::ELEMENTS::Dcohesive* ele = dynamic_cast<DRT::ELEMENTS::Dcohesive*>(structdis_->gElement( cohid ));

      if( ele->isFailed() )
      {
        failno++;
        std::cout<<"fail no = "<<failno<<"\n";//blockkk
      }
      else
      {
        allfailed = false;
        workele.push_back( cohid );
      }
    }

    std::cout<<"workele = "<<workele.size()<<"\tnumnode = "<<crs1.mas_->NumNode()<<"\n";//blockkk
    crs1.mas_->Print(std::cout);//blockkk

    // no cohesive elements attached with this crack surface fails
    //if( workele.size() == crs1.mas_->NumNode() )
    //  continue;

    std::cout<<"number faile = "<<failno<<"\t working = "<<workele.size()<<"\n";//blockkk

    std::cout<<"failed ele = "<<failno<<"\n";//blockkk

    //TODO: Check whether this strategy works for triangular elements
    if( allfailed or failno > 1 )
    {
      crk_sur.push_back( crs1.mas_ );
      crk_sur.push_back( crs1.sla_ );
      crs1.setProcInfo( true );
    }

    // TODO: this works only for quad elements
    // this means that the virtual crack closing elements should be added with these spring elements
    if( workele.size() == 2 )
    {
      int nodeids[4];

      const int * ids1 = structdis_->gElement( workele[0] )->NodeIds();
      const int * ids2 = structdis_->gElement( workele[1] )->NodeIds();

      // always [0] denotes master and [1] denotes slave nodes of spring element
      // we want to make sure here that the nodes are ordered to form non-intersecting Quad
      // This is possible only if we go master1(m1)-slave1(s1)-slave2(s2)-master2(m2) or in reverse
      //
      //                 !             !
      //                 !             !
      //                m2---/\/\/\---s2
      //                /             /
      //             ! /           ! /
      //             !/            !/
      //             m1---/\/\/\---s1
      //

      nodeids[0] = ids1[0]; nodeids[1] = ids1[1];
      nodeids[2] = ids2[1]; nodeids[3] = ids2[0];

      int neweleid = boundary_dis->NumGlobalElements();
      Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3","quad4", neweleid, structdis_->gElement( workele[0] )->Owner() );
      spr->SetNodeIds( 4, nodeids );

      tipeleno++;
    }
  }

  //if( tipeleno > 0 )    //unblockkk
  //{
    deleteOldTipElements( boundary_dis, tipele );
    boundary_dis->FillComplete();
 // }

  boundary_dis->Print(std::cout);//blockkk
  std::cout<<"number of crack surfaces to be = "<<crk_sur.size()<<"\n";//blockkk

  //for( unsigned ma=0; ma<mascrk.size();ma++ )
  //{
   // DRT::Element* masele = masterCrackDis_->gElement( mascrk[ma] );
   // this->addThisElementBoundary( boundary_dis, masele );
  //}

  //-------------
  std::map<std::string,Teuchos::RCP<DRT::Condition> > fool;
  const DRT::Condition* co1 = boundary_dis->GetCondition( "FSICoupling" );
  Teuchos::RCP<DRT::Condition> condi1 = Teuchos::rcp(new DRT::Condition(*co1));

  const DRT::Condition* co2 = boundary_dis->GetCondition( "XFEMCoupling" );
  Teuchos::RCP<DRT::Condition> condi2 = Teuchos::rcp(new DRT::Condition(*co2));

  fool["FSICoupling"] = condi1;
  fool["XFEMCoupling"] = condi2;
  //-------------

  for( std::vector<DRT::Element*>::iterator itcrk = crk_sur.begin(); itcrk!=crk_sur.end(); itcrk++ )
  {
    DRT::Element* cr = *itcrk;

    cr->Print(std::cout);//blockkk
    std::cout<<"crack element added\n";//blockkk

    int neweleid = boundary_dis->NumGlobalElements();
    Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3","quad4", neweleid, cr->Owner() );

    spr->SetNodeIds( cr->NumNode(), cr->NodeIds() );

    spr->Print(std::cout);//blockkk
    //dserror("got this new elt");//blockkk

    //const int* nodesele = cr->NodeIds();
    DRT::Node** nodesele = cr->Nodes();
    for( int noid = 0; noid<cr->NumNode(); noid++ )
    {

      const DRT::Node* nod = nodesele[noid];
      int gloid = nod->Id();
      if( not boundary_dis->HaveGlobalNode( gloid ) )
      {
        std::cout<<"i got old node\n";//blockkk
        Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp( new DRT::Node( gloid, nod->X(), nod->Owner() ) );

        std::cout<<"new node created\n";//blockkk
        boundary_dis->AddNode( newnode );

        newnode->SetCondition( fool.begin()->first, fool.begin()->second );
        newnode->SetCondition( (fool.begin()++)->first, (fool.begin()++)->second );
      }
    }
    std::cout<<"added the required nodes\n";//blockkk

    boundary_dis->AddElement( spr );

    boundary_dis->FillComplete();

    std::cout<<"element added once\n";//blockkk
  }

  std::cout<<"have I added all elements?\n";//blockkk

  // after adding the elements, check whether all possible crack surfaces are added
  this->checkAllDone();

  // if to_add > 0 and tipeleno > 0 -- meaning we are adding atleast one crack surface into cut surface
  // unless all crack surfaces are processed, atleast one virtual tip closing element should be added
  //if( (not allDone_) and tipeleno == 0 and crk_sur.size() > 0 )
  //{
   // dserror( "either all crack surfaces should have been processed or at least one virtual tip closing element should"
    //    "be added at each time step" );
  //}

  //boundary_dis->Print(std::cout);//blockkk
  //dserror("done");//blockkk

  boundary_dis->Print(std::cout);//blockkk
  //dserror("done");//blockkk
}*/

void ADAPTER::FSICrackingStructure::addCrackSurfacesToCutSides( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                                std::map<int, LINALG::Matrix<3,1> >& tip_nodes )
{
  if( allDone_ )
  {
    return;
  }

  std::vector<DRT::Element*> crk_sur;           // crack surface elements to be added to cut sides

  std::cout<<"----------------------I am inserting new crack surface elemts--------------------\n";//blockkk

  std::cout<<"number of crack surfaces = "<<possSurfaces_.size()<<"\n";//blockkk

  for( std::map<DRT::Element*,crackSurface_>::iterator it = possSurfaces_.begin();
                                                       it != possSurfaces_.end(); it++ )
  {
    crackSurface_& crs1 = it->second;

    // this master-slave combination is already processed
    if( crs1.isDone_ )
      continue;

    // element ids of all cohesive elements attached with this surface
    std::vector<int> cohele = crs1.dcoh_;
                     // no of failed cohesive elements
    bool allfailed = true;              // if all cohesive elements failed
    std::vector<int> workele, failele;           // Ids of cohesive elements that are still working

    //structdis_->Print(std::cout);//blockkk

    for( unsigned cohno=0; cohno<cohele.size();cohno++ )
    {
      int cohid = cohele[cohno];
      DRT::ELEMENTS::Dcohesive* ele = dynamic_cast<DRT::ELEMENTS::Dcohesive*>(structdis_->gElement( cohid ));

      if( ele->isFailed() )
      {
        failele.push_back( cohid );
        //std::cout<<"fail no = "<<failele.size()<<"\n";//blockkk
      }
      else
      {
        allfailed = false;
        workele.push_back( cohid );
      }
    }

    //std::cout<<"workele = "<<workele.size()<<"\tnumnode = "<<crs1.mas_->NumNode()<<"\n";//blockkk
    //crs1.mas_->Print(std::cout);//blockkk

    // no cohesive elements attached with this crack surface fails
    //if( workele.size() == crs1.mas_->NumNode() )
    //  continue;

    //std::cout<<"number faile = "<<failele.size()<<"\t working = "<<workele.size()<<"\n";//blockkk

    //std::cout<<"failed ele = "<<failele.size()<<"\n";//blockkk

    //TODO: Check whether this strategy works for triangular elements
    if( allfailed or failele.size() > 1 )
    {
      crk_sur.push_back( crs1.mas_ );
      crk_sur.push_back( crs1.sla_ );
      crs1.setProcInfo( true );
    }

    // If a cohesive element has failed, the nodes of these elements cannot be
    // defining crack tip
    for( unsigned failno = 0; failno < failele.size(); failno++ )
    {
      int coheleno = failele[failno];
      std::pair<int,int> pai = coheEleMasSla_[coheleno];

      std::map<int, LINALG::Matrix<3,1> >::iterator it1 = tip_nodes.find( pai.first );
      std::map<int, LINALG::Matrix<3,1> >::iterator it2 = tip_nodes.find( pai.second );

      if( it1 == tip_nodes.end() and it2 == tip_nodes.end() )
        continue;

#ifdef DEBUG
      if( ( it1 != tip_nodes.end() and it2 == tip_nodes.end() ) or
          ( it1 == tip_nodes.end() and it2 != tip_nodes.end() ) )
      {
        dserror( "if one node of a cohesive element is on crack tip, another node should also be so" );
      }
#endif

      tip_nodes.erase( it1 );
      tip_nodes.erase( it2 );
    }
  }

  //std::cout<<"number of crack surfaces to be = "<<crk_sur.size()<<"\n";//blockkk
  //if( crk_sur.size() > 0 )
  //  dserror("done");

  /**************************************///blockkk
  /*for( std::vector<DRT::Element*>::iterator itcrk = crk_sur.begin(); itcrk!=crk_sur.end(); itcrk++ )
  {
    DRT::Element* cr = *itcrk;

    cr->Print(std::cout);//blockkk
  }
  dserror("okay");*/
  /**************************************///blockkk

  //for( unsigned ma=0; ma<mascrk.size();ma++ )
  //{
   // DRT::Element* masele = masterCrackDis_->gElement( mascrk[ma] );
   // this->addThisElementBoundary( boundary_dis, masele );
  //}

  //-------------
  std::map<std::string,Teuchos::RCP<DRT::Condition> > fool;
  const DRT::Condition* co1 = boundary_dis->GetCondition( "FSICoupling" );
  Teuchos::RCP<DRT::Condition> condi1 = Teuchos::rcp(new DRT::Condition(*co1));

  const DRT::Condition* co2 = boundary_dis->GetCondition( "XFEMCoupling" );
  Teuchos::RCP<DRT::Condition> condi2 = Teuchos::rcp(new DRT::Condition(*co2));

  fool["FSICoupling"] = condi1;
  fool["XFEMCoupling"] = condi2;
  //-------------

  //TODO: Check whether the newly added elements preserves correct normal direction
  // if not decide and set this normal appropriately
  for( std::vector<DRT::Element*>::iterator itcrk = crk_sur.begin(); itcrk!=crk_sur.end(); itcrk++ )
  {
    DRT::Element* cr = *itcrk;

    //cr->Print(std::cout);//blockkk
    //std::cout<<"crack element added\n";//blockkk

    int neweleid = boundary_dis->NumGlobalElements();
    Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3","quad4", neweleid, cr->Owner() );

    spr->SetNodeIds( cr->NumNode(), cr->NodeIds() );

    //spr->Print(std::cout);//blockkk
    //dserror("got this new elt");//blockkk

    //const int* nodesele = cr->NodeIds();
    DRT::Node** nodesele = cr->Nodes();
    for( int noid = 0; noid<cr->NumNode(); noid++ )
    {

      const DRT::Node* nod = nodesele[noid];
      int gloid = nod->Id();
      if( not boundary_dis->HaveGlobalNode( gloid ) )
      {
        Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp( new DRT::Node( gloid, nod->X(), nod->Owner() ) );

        boundary_dis->AddNode( newnode );

        newnode->SetCondition( fool.begin()->first, fool.begin()->second );
        newnode->SetCondition( (fool.begin()++)->first, (fool.begin()++)->second );
      }
    }

    boundary_dis->AddElement( spr );

    boundary_dis->FillComplete();
  }

  // after adding the elements, check whether all possible crack surfaces are added
  this->checkAllDone();

  // if to_add > 0 and tipeleno > 0 -- meaning we are adding atleast one crack surface into cut surface
  // unless all crack surfaces are processed, atleast one virtual tip closing element should be added
  //if( (not allDone_) and tipeleno == 0 and crk_sur.size() > 0 )
  //{
   // dserror( "either all crack surfaces should have been processed or at least one virtual tip closing element should"
    //    "be added at each time step" );
  //}
}

void ADAPTER::FSICrackingStructure::deleteOldTipElements( Teuchos::RCP<DRT::Discretization>& boundary,
                                                          std::map<int, Teuchos::RCP<DRT::Element> >& tipele )
{

  for( std::map<int, Teuchos::RCP<DRT::Element> >::iterator it = tipele.begin(); it != tipele.end(); it++ )
  {
    Teuchos::RCP<DRT::Element> ti = it->second;

#ifdef DEBUG
    if( not boundary->HaveGlobalElement( ti->Id() ) )
      dserror("tip element not found in boundary discretization\n");
#endif

    boundary->DeleteElement( ti );
  }

  tipele.clear();
}

void ADAPTER::FSICrackingStructure::addThisElementBoundary( Teuchos::RCP<DRT::Discretization>& boundary_dis,
                                                            DRT::Element* ele )
{
  int neweleid = boundary_dis->NumGlobalElements();
  Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("BELE3","quad4", neweleid, ele->Owner() );

  spr->SetNodeIds( ele->NumNode(), ele->NodeIds() );



  DRT::Node** nodes = spr->Nodes();
  for( int noid = 0; noid<spr->NumNode(); noid++ )
  {
    DRT::Node* node = nodes[noid];
    if( not boundary_dis->HaveGlobalNode( node->Id() ) )
    {
      //Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp( new DRT::Node( old->Id(), old->X(), old->Owner() ) );
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp( node->Clone() );
      boundary_dis->AddNode( newnode );
    }
  }

  //-----------------
  std::map<std::string,Teuchos::RCP<DRT::Condition> > cond;
  DRT::Element* el = boundary_dis->lRowElement(0);
  DRT::Condition* co1 = el->GetCondition( "FSICoupling" );
  Teuchos::RCP<DRT::Condition> condi = Teuchos::rcp(new DRT::Condition(*co1));
  spr->SetCondition( "FSICoupling", condi );
  //-----------------

  boundary_dis->AddElement( spr );


  boundary_dis->FillComplete();
}
