/*----------------------------------------------------------------------*/
/*!
\file InsertCohesiveElements.cpp

\brief Insertion of cohesive elements in to the discretization

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
#include "InsertCohesiveElements.H"


#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_cut/cut_boundarycell.H"


/*------------------------------------------------------------------------------------*
 * Insert discrete cohesive spring elements into existing structural discretization
 * Check from input parameters whether to insert springs in between the given         sudhakar 04/13
 * surfaces or throughout the discretization, and call appropriate functions
 *------------------------------------------------------------------------------------*/
DRT::CRACK::InsertCohesiveElements::InsertCohesiveElements( Teuchos::RCP<DRT::Discretization> dis )
          : discret_ ( dis )
{
  std::cout<<"I am inserting cohesive elements \n";

  // see if crack surfaces are predefined
  const Teuchos::ParameterList& params = DRT::Problem::Instance()->CrackParams();
  predefined_ = DRT::INPUT::IntegralValue<int>(params,"IS_CRACK_PREDEFINED")==1;

  // insert spring elements connecting predefined surfaces
  if( predefined_ )
    specifiedSurfaces();

  // insert spring elements at all nodes in the discretization
  else
    WholeDiscret();

  //dserror("inserted cohesive elements\n");//blockkk
}

/*------------------------------------------------------------------------------------*
 * Insert spring (discrete cohesive) elements over the specified surfaces
 * given in the input files                                                   sudhakar 04/13
 * Just inserting elements is sufficient because in this case the crack
 * already has two separate surfaces --> no need to modify existing elements
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::InsertCohesiveElements::specifiedSurfaces()
{
  // get master and slave nodes
  DRT::Condition* mas = discret_->GetCondition( "masterCrackSurface" );
  DRT::Condition* sla = discret_->GetCondition( "slaveCrackSurface" );

  master_ = const_cast<std::vector<int>* >(mas->Nodes());
  slave_ = const_cast<std::vector<int>* >(sla->Nodes());

  /****************************************************/
  /*std::cout<<"master nodes are = ";
  for( std::vector<int>::const_iterator i=master_->begin(); i!=master_->end();i++ )
  {
    const int m = *i;
    std::cout<<m<<"\n";
  }*/
  /****************************************************/

  /*---------------------------------------------------------------------------------*/
  /*std::cout<<"master nodes = "<<master->size()<<"\n";
  std::cout<<"slave nodes = "<<slave->size()<<"\n";

  std::cout<<"master nodes\n";
  for(std::vector<int>::const_iterator i=master->begin(); i!=master->end();i++ )
  {
    const int k = *i;
    std::cout<<k<<"\n";
  }

  std::cout<<"slave nodes\n";
  for(std::vector<int>::const_iterator i=slave->begin(); i!=slave->end();i++ )
  {
    const int k = *i;
    std::cout<<k<<"\n";
  }*/
  /*---------------------------------------------------------------------------------*/ //blockkk

  if( master_->size() != slave_->size() )
    dserror("both master and slave crack surfaces should have same number of nodes\n");

  // for each master node there is a corresponding slave node at the same position
  // connect these master-slave pair with a spring element
  bool found = false;
  for(std::vector<int>::const_iterator i=master_->begin(); i!=master_->end();i++ )
  {
    found = false;
    const int m = *i;
    // get master node
    const DRT::Node* nMas = discret_->gNode( m );

    // search whether any of the slave nodes is in same position as master node
    for(std::vector<int>::const_iterator j=slave_->begin(); j!=slave_->end();j++ )
    {
      const int s = *j;
      const DRT::Node* nSla = discret_->gNode( s );

      found = AreAtSameLocation( nMas, nSla );
      // if master-slave pair found, add spring element
      if( found )
      {
        AddSpringWithThisNodes( m, s );
        break;
      }
    }

    if( not found )
      dserror("All master crack nodes should have a corresponding slave node\n");
  }

  // only after fill complete, discretization know the new elements
  discret_->FillComplete( true, true, true );
}

/*----------------------------------------------------------------------------------------------*
 * Insert spring (discrete cohesive) elements at every surface in discret_
 * We need to insert duplicate nodes at each existing nodes
 * insert spring elements and modify the connectivity of existing elements
 *----------------------------------------------------------------------------------------------*/
void DRT::CRACK::InsertCohesiveElements::WholeDiscret()
{
  // we go through all the nodes in the discretization
  /*static int numnode = discret_->NumGlobalNodes();

  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    // get row node via local id
    const DRT::Node* nod = discret_->lRowNode(i);

    int numele = nod->NumElement();
    // node on the corners, so nothing to do
    if( numele == 1 )
      continue;

    DRT::Element** adjeles = nod->Elements();
    for(int iele=0; iele<numele; iele++)
    {
      DRT::Element * ele1 = adjeles[iele];
      //if(adjeles[iele]->Owner() == myrank_)
      //  fluideles[binId].insert(adjeles[iele]->Id());
    }
  }*/
  // In this case while dealing with each node, it has to be checked whether it holds any
  // Dirichlet or Neumann boundary conditions -- and all duplicate nodes generated from this node
  // should be added to the boundary condition
  // for Dirichlet bc -- leave this node. this should not carry any springs
  // for Neumann boundary condition -- dont know what to do????
  dserror("at the moment you can simulate only cracks with predefined crack surfaces\n");
}

/*----------------------------------------------------------------------------------------------*
 * Check whether two nodes given as input have the same position coordinates          sudhakar 04/13
 *----------------------------------------------------------------------------------------------*/
bool DRT::CRACK::InsertCohesiveElements::AreAtSameLocation( const DRT::Node* n1, const DRT::Node* n2 )
{
  const double* masco = n1->X();
  const double* slaco = n2->X();

  if( fabs(masco[0]-slaco[0]) > 1e-14 or
      fabs(masco[1]-slaco[1]) > 1e-14 or
      fabs(masco[2]-slaco[2]) > 1e-14 )
  {
    return false;
  }
  return true;
}

/*----------------------------------------------------------------------------------------------*
 * Create spring (discrete cohesive) element connecting the given two nodes             sudhakar 04/13
 * add this element to the discretization
 *----------------------------------------------------------------------------------------------*/
void DRT::CRACK::InsertCohesiveElements::AddSpringWithThisNodes( const int id1,
                                                                 const int id2 )
{
  // defined as static because once we add an element NumGlobalElements()
  // is no more available
  static int idele = discret_->NumGlobalElements();
  Teuchos::RCP<DRT::Element> spr = DRT::UTILS::Factory("DCOHESIVE","line2", idele,
                                                        discret_->gNode(id1)->Owner() );
  idele++;

  int ids[2];

  ids[0] = id1;
  ids[1] = id2;

  // set nodal ids of this spring element and add it to discretization
  spr->SetNodeIds( 2, ids );
  discret_->AddElement( spr );

  // set normal vector for the element
  LINALG::Matrix<3,1> normal;
  normal(0,0) = 1.0;
  Teuchos::RCP<DRT::ELEMENTS::Dcohesive> coh = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Dcohesive>(spr);
  coh->setRefNormal( normal );


  // compute and set the reference area of this element
  double area = determineReferenceArea( id1 );
  coh->setArea( area );
  if( area < 0.0 or fabs(area) < 1e-14 )
  {
    std::cout<<"area of spring element = "<<area<<"\n";
    dserror("zero or negative area found for a spring element\n");
  }
  //std::cout<<area<<"\n";//blockkk

}

/*----------------------------------------------------------------------------------------------*
 * Calculate the area over which this element acts and initialize it                       sudhakar 04/13
 *----------------------------------------------------------------------------------------------*/
double DRT::CRACK::InsertCohesiveElements::determineReferenceArea( const int node_id )
{
  double area = 0.0;

  // when crack propagation is specified, at each node there is only
  // one cohesive element
  // WARNING :: This is not true if crack branching exist in which case our method is wrong
  if( predefined_ )
  {
    //-----------------------------------------------------------------------
    // STEP 1 -- get all elements connected by this node
    //-----------------------------------------------------------------------
    DRT::Node* nod = discret_->gNode( node_id );
    DRT::Element** elems = nod->Elements();
    int numelems = nod->NumElement();

    for( int ne=0; ne < numelems; ne++  )
    {
      //-----------------------------------------------------------------------
      // STEP 2 -- get all surfaces of elements to which this node is part of
      //-----------------------------------------------------------------------
      DRT::Element* ele = elems[ne];
      std::vector< Teuchos::RCP< DRT::Element > > surf = ele->Surfaces();

      for (unsigned int surfele = 0; surfele < surf.size(); surfele++)
      {
        bool found = true;
        //-------------------------------------------------------------------------
        // STEP 3 -- Spring can act only one surface of an element, and this surface
        //---------- must springs at all nodes --> check this
        //-------------------------------------------------------------------------
        const int* searchnodes = surf[surfele]->NodeIds();
        const int no_nodes = surf[surfele]->NumNode();
#ifdef DEBUG
        if (searchnodes==NULL) dserror("No Nodes in Surface Element");
#endif
        /************************************************************/
        /*std::cout<<"search nodes are = ";
        for (int num = 0; num < no_nodes; num++)
        {
          int thisnode = searchnodes[num];
          std::cout<<thisnode<<"\t";
        }
        std::cout<<"\n";*/
        /************************************************************/
        for (int num = 0; num < no_nodes; num++)
        {

          int thisnode = searchnodes[num];
          if( thisnode == node_id )
            continue;
          if( not (std::find(master_->begin(),master_->end(), thisnode) != master_->end()) )
          {
            found = false;
            break;
          }
        }

        //--------------------------------------------------------------------------------------
        // STEP 4 -- calculate the area contribution for this surface = surface area/no.of nodes
        //--------------------------------------------------------------------------------------
        if( found )
        {
          area += areaSurface( searchnodes, no_nodes )/no_nodes;
          // only one surface of an element can contribute to area of spring elt
          // we found it and now we go to the next element
          break;
        }
      }
    }
  }
  else
    dserror( "you can simulate only known crack paths at the moment\n" );

  //area = 0.0001;//blockkk
  return area;
}

/*------------------------------------------------------------------------------------------------*
 * calculate the area of surface given by these nodes                                 sudhakar 04/13
 *------------------------------------------------------------------------------------------------*/
double DRT::CRACK::InsertCohesiveElements::areaSurface( const int* surfnodes, int num )
{
  double area = 0.0;
  Epetra_SerialDenseMatrix xyz( 3, num );

  // get coordinates of all nodes of a surface
  for (int i = 0; i < num; i++)
  {
    int nodeid = surfnodes[i];
    const double * coo = discret_->gNode( nodeid )->X();

    for( unsigned j=0;j<3;j++ )
      xyz( j,i ) = coo[j];
  }

  // calculate area of surface
  if( num == 3 )
  {
    GEO::CUT::Tri3BoundaryCell * bc = new GEO::CUT::Tri3BoundaryCell( xyz, NULL, std::vector<GEO::CUT::Point*>() );
    area = bc->Area();
  }
  else if( num == 4 )
  {
    GEO::CUT::Quad4BoundaryCell * bc = new GEO::CUT::Quad4BoundaryCell( xyz, NULL, std::vector<GEO::CUT::Point*>() );
    area = bc->Area();
  }
  else
    dserror( "surface of an FEM element should be a Tri or Quad\n" );

  return area;
}
