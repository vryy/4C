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
#include "../drt_lib/drt_linedefinition.H"

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

  // We assume that first 4 nodes contains a face in z-plane
  // next 4 nodes has another face
  // TODO: Check whether this is always true with HEX
  std::vector<int> face1, face2;

  for( int i=0;i<4;i++ )
    face1.push_back( elenodes[i] );
  for( int i=4;i<8;i++ )
    face2.push_back( elenodes[i] );

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

/*-------------------------------------------------------------------------------------------*
 * Give Quad element is split into two Tri elements with a diagonal                sudhakar 08/14
 * which has one tip node and one split node at its ends
 *
 *
 *      -------*                 * splitnode
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
