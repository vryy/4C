/*----------------------------------------------------------------------*/
/*!
\file SplitHexIntoTwoWedges.H

\brief Split one Hex element in a discretization into two Wedge elements
and rebuild the discretization

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "SplitHexIntoTwoWedges.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_lib/drt_utils_factory.H"

#include <iostream>
#include <sstream>

/*---------------------------------------------------------------------------------------------*
 * Do all operations on splitting HEX into WEDGE two WEDGE elements
 * First two Quad elements that are on z-plane are identified and split into Tri           sudhakar 08/14
 * Then these Tri are joined approp. to generate WEDGE
 * REMEMBER: It has already been checked both tipnodes and splitnodes have only 2 values
 *---------------------------------------------------------------------------------------------*/
void DRT::CRACK::SplitHexIntoTwoWedges::DoAllSplittingOperations( DRT::Element * ele,
                                                                  std::vector<int> tipnodes,
                                                                  std::vector<int> splitnodes,
                                                                  int neweleid1,
                                                                  int neweleid2 )
{
  const int * elenodes = ele->NodeIds();

  std::vector<int> face1, face2;

#if 0
  // We assume that first 4 nodes contains a face in z-plane
  // next 4 nodes has another face
  // not working for many cases
  for( int i=0;i<4;i++ )
    face1.push_back( elenodes[i] );
  for( int i=4;i<8;i++ )
    face2.push_back( elenodes[i] );
#else
  GetTwoCorrectFaces( elenodes, tipnodes, splitnodes, face1, face2 );
#endif

  std::vector<int> tri1, tri2, tri3, tri4;
  SplitThisQuad( face1, tipnodes, splitnodes, tri1, tri2 );
  SplitThisQuad( face2, tipnodes, splitnodes, tri3, tri4 );

  std::vector<int> wedge1_nodes, wedge2_nodes;
  for( unsigned i=0; i<tri1.size(); i++ )
    wedge1_nodes.push_back( tri1[i] );
  for( unsigned i=0; i<tri3.size(); i++ )
    wedge1_nodes.push_back( tri3[i] );

  for( unsigned i=0; i<tri2.size(); i++ )
    wedge2_nodes.push_back( tri2[i] );
  for( unsigned i=0; i<tri4.size(); i++ )
    wedge2_nodes.push_back( tri4[i] );


  int owner = ele->Owner();

  AddThisWedge( owner, neweleid1, wedge1_nodes );
  AddThisWedge( owner, neweleid2, wedge2_nodes );

#if 0
  // We should not get surfaces from element and proceed with it
  // Because when we get surfaces, the ordering of nodes are changed
  DRT::Node * splnode = discret_->gNode( spl_id );
  std::vector< Teuchos::RCP< DRT::Element > > elesurfaces = ele->Surfaces();
  Teuchos::RCP<DRT::Element> surele = getSurfaceSameZplane( elesurfaces, splnode->X() );

  const int* surnodes = surele->NodeIds();
  int numsurnodes = surele->NumNode();

  int tip_index = 0;
  bool found_tip_index = false;
  for( int i=0; i<numsurnodes; i++ )
  {
    if( std::find( tipnodes_.begin(), tipnodes_.end(), surnodes[i] ) != tipnodes_.end() )
    {
      tip_index = i;
      found_tip_index = true;
    }
  }

  if( not found_tip_index )
    dserror("Each surface should have a tip node\n");

  std::vector<int> tri1, tri2;
  tri1.push_back( surnodes[tip_index] );
  tri1.push_back( surnodes[(tip_index+1)%numsurnodes] );
  tri1.push_back( surnodes[(tip_index+2)%numsurnodes] );

  tri2.push_back( surnodes[(tip_index+2)%numsurnodes] );
  tri2.push_back( surnodes[(tip_index+3)%numsurnodes] );
  tri2.push_back( surnodes[tip_index] );
#endif
}

/*-----------------------------------------------------------------------------------------------*
 * From given HEX element nodes, get two Quad faces that will be split                   sudhakar 09/14
 * into Tri while forming WEDGE element
 *-----------------------------------------------------------------------------------------------*/
void DRT::CRACK::SplitHexIntoTwoWedges::GetTwoCorrectFaces( const int * elenodes,
                                                            std::vector<int>& tipnodes,
                                                            std::vector<int>& splitnodes,
                                                            std::vector<int>& face1,
                                                            std::vector<int>& face2 )
{
  // The faces are HEX8 are ordered in the following fashion in such a way
  // that this nodal arrangement always give outward pointing normal
  // Face0 -- 0 3 2 1 ----- (a)
  // Face1 -- 0 1 5 4 ----- (b)
  // Face2 -- 1 2 6 5 ----- (c)
  // Face3 -- 2 3 7 6 ----- (b)
  // Face4 -- 0 4 7 3 ----- (c)
  // Face5 -- 4 5 6 7 ----- (a)
  // (a), (b) and (c) mark paried surfaces
  // One can see that face1 must be from first 3 sides, and face2 must be from last 3 sides
  // Refer to BACI report -- single_report_conventions.pdf -- for more details

  int number = 0; // not even one surface found
  bool found_face2 = false;
  for( unsigned numsurf = 3; numsurf < 6; numsurf++ )
  {
    std::vector<int> thisfacenodes;
    if( numsurf == 3 )
    {
      thisfacenodes.push_back(elenodes[2]);
      thisfacenodes.push_back(elenodes[3]);
      thisfacenodes.push_back(elenodes[7]);
      thisfacenodes.push_back(elenodes[6]);
    }
    else if( numsurf == 4 )
    {
      thisfacenodes.push_back(elenodes[0]);
      thisfacenodes.push_back(elenodes[4]);
      thisfacenodes.push_back(elenodes[7]);
      thisfacenodes.push_back(elenodes[3]);
    }
    else if( numsurf == 5 )
    {
      thisfacenodes.push_back(elenodes[4]);
      thisfacenodes.push_back(elenodes[5]);
      thisfacenodes.push_back(elenodes[6]);
      thisfacenodes.push_back(elenodes[7]);
    }
    else
      dserror( "HEX8 should have only 6 surfaces\n" );

    if( std::find( thisfacenodes.begin(), thisfacenodes.end(), tipnodes[1] ) != thisfacenodes.end() and
        std::find( thisfacenodes.begin(), thisfacenodes.end(), tipnodes[0] ) == thisfacenodes.end() )
    {
      if( std::find( thisfacenodes.begin(), thisfacenodes.end(), splitnodes[0] ) != thisfacenodes.end() or
          std::find( thisfacenodes.begin(), thisfacenodes.end(), splitnodes[1] ) != thisfacenodes.end() )
      {
        number = numsurf;
        face2 = thisfacenodes;
        found_face2 = true;
      }
    }
  }

  if( not found_face2 )
    dserror( "face2 is not found\n" );

  if( number == 3 )
  {
    face1.push_back( elenodes[1] );
    face1.push_back( elenodes[0] );
    face1.push_back( elenodes[4] );
    face1.push_back( elenodes[5] );
  }
  else if( number == 4 )
  {
    face1.push_back( elenodes[1] );
    face1.push_back( elenodes[5] );
    face1.push_back( elenodes[6] );
    face1.push_back( elenodes[2] );
  }
  else if( number == 5 )
  {
    face1.push_back( elenodes[0] );
    face1.push_back( elenodes[1] );
    face1.push_back( elenodes[2] );
    face1.push_back( elenodes[3] );
  }
  else
    dserror( "number should be between this given values\n" );

  if( face1.size() != 4 )
  {
    std::cout<<"face1 nodes are = \n";
    for( std::vector<int>::iterator it = face1.begin(); it != face1.end(); it++ )
      std::cout<<*it<<" ";
    std::cout<<"\n";
    dserror( "face1 does not have 4 nodes\n" );
  }
  if( face2.size() != 4 )
  {
    std::cout<<"face2 nodes are = \n";
    for( std::vector<int>::iterator it = face2.begin(); it != face2.end(); it++ )
      std::cout<<*it<<" ";
    std::cout<<"\n";
    dserror( "face2 does not have 4 nodes\n" );
  }
}

/*-------------------------------------------------------------------------------------------*
 * Give Quad element is split into two Tri elements with a diagonal
 * which has one tip node and one split node at its ends
 *
 *
 *      -------*                 * splitnode                                          sudhakar 08/14
 *      |     /|                 o tip node
 *      |    / |
 *      |   /  |
 *      |  /   |
 *      | /    |
 *      |/     |
 *      o------
 *-------------------------------------------------------------------------------------------*/
void DRT::CRACK::SplitHexIntoTwoWedges::SplitThisQuad( std::vector<int> & face,
                                                       std::vector<int> & tipnodes,
                                                       std::vector<int> & splitnodes,
                                                       std::vector<int> & tri1,
                                                       std::vector<int> & tri2 )

{
  int tip_index = 0;
  bool found_tip_index = false;
  for( unsigned i=0; i<face.size(); i++ )
  {
    if( std::find( tipnodes.begin(), tipnodes.end(), face[i] ) != tipnodes.end() )
    {
      tip_index = i;
      found_tip_index = true;
    }
  }

  if( not found_tip_index )
    dserror("Each face must have a tip node\n");

  tri1.push_back( face[tip_index] );
  tri1.push_back( face[(tip_index+1)%face.size()] );
  tri1.push_back( face[(tip_index+2)%face.size()] );

  tri2.push_back( face[(tip_index+2)%face.size()] );
  tri2.push_back( face[(tip_index+3)%face.size()] );
  tri2.push_back( face[tip_index] );
}

/*--------------------------------------------------------------------------------------------*
 * Add the generated wedge element into discretization                                sudhakar 08/14
 * It is mandatory to set the material here as this element is evaluated as well
 *--------------------------------------------------------------------------------------------*/
void DRT::CRACK::SplitHexIntoTwoWedges::AddThisWedge( int owner,
                                                      int eleid,
                                                      std::vector<int> elenodes )
{
  Teuchos::RCP<DRT::Element> wedge = DRT::UTILS::Factory("SOLIDW6","Polynomial", eleid, owner );
  wedge->SetNodeIds( 6, &elenodes[0] );

  wedge->SetMaterial( material_id_ );

  discret_->AddElement( wedge );
}
