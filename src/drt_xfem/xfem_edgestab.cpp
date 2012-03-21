/*!
\file xfem_edgestab.cpp

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented, continuous interior penalty) scheme

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/


#include "xfem_edgestab.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_cut/cut_elementhandle.H"
//#include "../drt_cut/cut_side.H"
#include "../drt_cut/cut_volumecell.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_boundary_calc.H"
#include "xfem_fluidwizard.H"


XFEM::XFEM_EdgeStab::XFEM_EdgeStab(
  Teuchos::RCP<XFEM::FluidWizard>              wizard
  ) :
  wizard_(wizard)
  {

    return;

  } // end constructor



void XFEM::XFEM_EdgeStab::findSurface(std::vector<int>            nodeids,
                                           DRT::ELEMENTS::Fluid3 *     actele,
                                           vector<RCP<DRT::Element> >  surfaces,
                                           int & surfid
                                           )
{


  // loop over surfaces
  for (unsigned int surf=0; surf<surfaces.size(); ++surf)
  {
    RCP< DRT::Element > act_surface = surfaces[surf];

    DRT::Node** surf_NodesPtr = act_surface->Nodes();
    int numsurfnodes = act_surface->NumNode();

    // assume the current surface is the right one
    bool surf_found = true;

    for(int n=0; n< (int)nodeids.size(); n++)
    {
      bool matching_node = false;
      for(int surf_n=0; surf_n<numsurfnodes; surf_n++)
      {
        if(nodeids[n] == (surf_NodesPtr[surf_n])->Id())
        {
          matching_node = true;
          break; // break the inner for loop
        }
      }
      if(!matching_node)
      {
        surf_found = false;
      }
    }

    if(surf_found)
    {
      surfid = surf;
      break; // break the surface loop
    }

  }



  return;
}



void XFEM::XFEM_EdgeStab::findSide(   std::vector<int>                nodeids, // node ids of surface
                                           DRT::ELEMENTS::Fluid3 *         actele,
                                           const vector<GEO::CUT::Side*> & sides,
                                           int & side_index
                                           )
{


  // loop over surfaces
  for (unsigned int side=0; side<sides.size(); ++side)
  {
    const GEO::CUT::Side* actside = sides[side];

    std::vector<GEO::CUT::Node* > side_nodes = actside->Nodes();

    // assume the current surface is the right one
    bool side_found = true;

    for(int n=0; n< (int)nodeids.size(); n++)
    {
      bool matching_node = false;
      for(size_t side_n = 0; side_n<side_nodes.size(); side_n++)
      {
        if(nodeids[n] == (side_nodes[side_n])->Id())
        {
          matching_node = true;
          break; // break the inner for loop
        }
      }
      if(!matching_node)
      {
        side_found = false;
      }
    }

    if(side_found)
    {
      side_index = side;
      break; // break the surface loop
    }

  }



  return;
}


void XFEM::XFEM_EdgeStab::findNeighborElement(DRT::ELEMENTS::Fluid3 * actele,
                                                   RCP<DRT::Element> surface,
                                                   bool & neighbor_found,
                                                   int & neighbor_id)
{
  // find for each surface its unique neighbor if there is a neighbor
  DRT::Node** NodesPtr = surface->Nodes();
//        const int* nodeids = surface->NodeIds();
  int numsurfnodes = surface->NumNode();

  for (int surfnode = 0; surfnode < surface->NumNode(); ++surfnode)
  {
    DRT::Node* node = NodesPtr[surfnode];

    // get adjacent element to this surface node
    DRT::Element** adjeles = node->Elements();
    int numadjele = node->NumElement();

    // search for the right neighbor element
    for(int adjele = 0; adjele < numadjele; adjele++)
    {
      DRT::Element* actadjele = adjeles[adjele];

      DRT::ELEMENTS::Fluid3 * adjele = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( actadjele );


      // decide if adjele is the right neighbor
      if(adjele->Id() == actele->Id());// cout << "found the element itself" << endl ;// element found itself
      else
      {
        // do adjele and ele share a common surface ?
        vector<RCP<DRT::Element> > adjele_surfaces = adjele->Surfaces();

        // loop over surfaces of adjacent element
        for (unsigned int adj_surf=0; adj_surf<adjele_surfaces.size(); ++adj_surf)
        {
          RCP< DRT::Element > adj_surface = adjele_surfaces[adj_surf];
          int numadjsurfnodes = adj_surface->NumNode();
          DRT::Node** AdjNodesPtr = adj_surface->Nodes();

          if(numsurfnodes != numadjsurfnodes)
          {
            neighbor_found = false;
          }
          else
          {

            // if adj_surface and surface share the same nodes, then adjele is the neighbor element
            int matching_node = 0;
            for(int i=0; i< numsurfnodes; ++i)
            {
              for(int j=0; j< numadjsurfnodes; ++j)
              {
                if(NodesPtr[i]->Id() == AdjNodesPtr[j]->Id()) matching_node++;
              }
            }
            if(matching_node == (int)numsurfnodes)
            {
              neighbor_found =  true;
              neighbor_id    =  adjele->Id();
            }
          }

        }


      }

    }
  }

}



void XFEM::XFEM_EdgeStab::stabilizeStandardElement( DRT::AssembleStrategy&   strategy,
                                                             DRT::ELEMENTS::Fluid3 *  actele,
                                                             DRT::Discretization &    discret)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::EdgeBasedStabilization" );

  // all surfaces of Fluid3 elements are Fluid3BoundaryElements
  // better: create Fluid3InternalSurfaces elements
  vector<RCP<DRT::Element> > surfaces = actele->Surfaces();

  // loop over surfaces
  for (unsigned int surf=0; surf<surfaces.size(); ++surf)
  {
//    cout << "surface number " << surf << endl;
    RCP< DRT::Element > surface = surfaces[surf];

    bool neighbor_found = false;
    int  neighbor_id = -1;


    findNeighborElement(actele, surface, neighbor_found, neighbor_id);


    // if there exists a neighbor
    if(neighbor_id > -1)
    {

      // get the neighboring element
      DRT::Element* neighbor = discret.gElement(neighbor_id);

      // is there an element handle for the neigboring element
      GEO::CUT::ElementHandle * e_neighbor = wizard_->GetElement( neighbor );

      bool edge_based_stab = false;
      bool ghost_penalty   = false;

      // neighbor element does not have an element handle
      if(e_neighbor == NULL)
      {
        // two uncut elements / standard fluid case
        edge_based_stab = true;
        ghost_penalty   = false;

        // the current element must be the parent element

        // get the parent element
        int ele1_id = actele->Id();
        int ele2_id = neighbor_id;

        std::vector<int> nds_1;
        std::vector<int> nds_2;

        nds_1.clear();
        for(int i=0; i< discret.gElement(ele1_id)->NumNode(); i++)
        {
          nds_1.push_back(0);
        }

        nds_2.clear();
        for(int i=0; i< discret.gElement(neighbor_id)->NumNode(); i++)
        {
          nds_2.push_back(0);
        }


        DRT::Element* ele1 = discret.gElement(ele1_id);
        DRT::Element* ele2 = discret.gElement(ele2_id);


        // call evaluate routine
        callEdgeStabandGhostPenalty( edge_based_stab,
                                     ghost_penalty,
                                     surface,
                                     ele1,
                                     ele2,
                                     nds_1,
                                     nds_2,
                                     discret,
                                     strategy );



      }
      else // neighbor element has an element handle
      {
        // find the side equal to the surface, get its unique! facet and volumecell -> nds

        // find the right side
        GEO::CUT::plain_element_set subelements;
        e_neighbor->CollectElements(subelements);

        if(e_neighbor->Shape() == DRT::Element::hex8)
        {
          if(subelements.size() != 1) dserror("wrong number of subelements in hex8 element");

          GEO::CUT::Element * e = subelements[0];

          const std::vector<GEO::CUT::Side* > sides = e->Sides();

          int side_index = -1;

          std::vector<int> nodeids;

          // find for each surface its unique neighbor if there is a neighbor
          DRT::Node** NodesPtr = surface->Nodes();

          int numsurfnodes = surface->NumNode();

          for(int i=0; i<numsurfnodes; i++)
          {
            nodeids.push_back(NodesPtr[i]->Id());
          }

          findSide(nodeids, actele, sides, side_index);


          const GEO::CUT::Side* side = sides[side_index];

          const std::vector<GEO::CUT::Facet*> facets = side->Facets();

          if(facets.size() != 1) dserror("more than one facet on side between standard element and element with handle");

          // get unique facet
          GEO::CUT::Facet* f = facets[0];

          if(f->Position() == GEO::CUT::Point::outside)
          {
            const GEO::CUT::plain_volumecell_set vcs = f->Cells();

            if(vcs.size() != 1) dserror("more than one vc on side between standard element and element with handle");

            GEO::CUT::VolumeCell* vc = vcs[0];
            if(vc->ParentElement()->Id() != neighbor_id) dserror("vc is not a cell of the element with handle!");

            // get the parent element
            int ele1_id = actele->Id();
            int ele2_id = neighbor_id;


            std::vector<int> nds_1;
            std::vector<int> nds_2;

            nds_1.clear();
            nds_2.clear();

            nds_1.clear();
            // parent element has no element handle
            for(int i=0; i< discret.gElement(ele1_id)->NumNode(); i++)
            {
              nds_1.push_back(0);
            }

            nds_2 = vc->NodalDofSet();

//            // TODO:
//            // remove this later
//            nds_2.clear();
//            for(int i=0; i< discret.gElement(ele2_id)->NumNode(); i++)
//            {
//              nds_2.push_back(0);
//            }


            edge_based_stab = true;

            // perform ghost penalty if the element with element handle is cut
            if(e->IsCut()) ghost_penalty = true;


            DRT::Element* ele1 = discret.gElement(ele1_id);
            DRT::Element* ele2 = discret.gElement(ele2_id);



            // call evaluate routine
            callEdgeStabandGhostPenalty( edge_based_stab,
                ghost_penalty,
                surface,
                ele1,
                ele2,
                nds_1,
                nds_2,
                discret,
                strategy );
          }
          else if(f->Position() == GEO::CUT::Point::undecided)
          {
            dserror("the position of this facet is undecided, how to stabilize???");
          }
          else if(f->Position() == GEO::CUT::Point::oncutsurface)
          {
            dserror("the position of this facet is oncutsurface, how to stabilize???");
          }
          else
          {
            // facet is inside!
          }





        }
        else
        {
          dserror("not implemented for this element type");
        }




      }
    } // if(neighbor)

  } // surfaces
}





void XFEM::XFEM_EdgeStab::EvaluateEdgeStabandGhostPenalty( DRT::Discretization &   discret,
                                                                DRT::AssembleStrategy&  strategy,
                                                                DRT::ELEMENTS::Fluid3 * ele)
{
#ifdef D_FLUID3

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::EdgeStabandGhostPenalty" );



  //====================================================================================================
  // implementation of edge-based fluid stabilization and ghost penalty
  //====================================================================================================



  // EDGE-based stabilization
  // REMARK: the current implementation of edge-based stabilization uses the volume-element routine
  // this means that the sides/surfaces of the current element has to be checked whether stabilization is required
  // -> better: there should be an own side-evaluate routine that knows its neighboring elements

  // we distinguish to different stabilization cases
  // the general non-XFEM case: the current element and its neighboring element (connected via current side)
  //                            do not have an elementhandle -> then the side/surface has to be stabilized
  //                            later: the side has to know its elements/elementhandles
  // the XFEM-case: at least one element adjacent to the side (either the current element itself or the neighbor element)
  //                have an elementhandle
  //                2 ehandles: stabilize for each facet and its adjacent volumecells (connected set of volumecells)
  //                1 ehandle: (is there a neighboring element (without an elementhandle?)
  //                           for the unique facet of this side stabilize the edge (nds for element without handle, nds = 000...)
  // additional stabilization for XFEM case: 2 ehandles and at least one element is Cut -> stabilize interface
  //                                         1 ehandle and its element is Cut -> stabilize interface


  // has the current element an element handle
  GEO::CUT::ElementHandle * e_parent   = wizard_->GetElement( ele );


  if(e_parent != NULL)
  {
    // get all its elements (for hex20 element we get 8 hex8 subelements)
    GEO::CUT::plain_element_set curr_elements;
    e_parent->CollectElements(curr_elements);

    // hex8 element
    if(e_parent->Shape() == DRT::Element::hex8)
    {
      if (curr_elements.size() != 1) dserror("there should be only one GEO::CUT::(sub)element for hex8 elements");

      GEO::CUT::Element* e = curr_elements[0];

      // get the sides
      const std::vector<GEO::CUT::Side*> sides = e->Sides();

      // TODO: later: here will the side-routine start

      for(std::vector<GEO::CUT::Side*>::const_iterator i=sides.begin(); i!=sides.end(); i++)
      {
        // find the corresponding Fluid3Boundary element ------------------------

        const std::vector<GEO::CUT::Node* > side_nodes = (*i)->Nodes();
        std::vector<int> side_node_ids;

        for(std::vector<GEO::CUT::Node*>::const_iterator it = side_nodes.begin(); it!=side_nodes.end(); it++)
        {
          side_node_ids.push_back((*it)->Id());
        }

        // find all sides/surfaces
        vector<RCP<DRT::Element> > surfaces = ele->Surfaces();

        int surfid = -1;

        // TODO: this part can be removed when the side knows its parent elements!
        findSurface(side_node_ids, ele, surfaces, surfid);

        if(surfid == -1) dserror("surface for side not found!!!");

        RCP<DRT::Element> surface = surfaces[surfid];

        //------------------------------------------------------------------------



        // find the neighboring elements
        // TODO: later: only the nds vector have to be found for the adjacent elements (these will be known by the side)
        DRT::Element* ele1 = NULL;
        DRT::Element* ele2 = NULL;


        // facet of current side
        const std::vector<GEO::CUT::Facet*> facets = (*i)->Facets();

        // each facet should have 2 volumecells
        for(std::vector<GEO::CUT::Facet*>::const_iterator f=facets.begin(); f!=facets.end(); f++)
        {
          if((*f)->Position() == GEO::CUT::Point::outside /*or (*f)->Position() == GEO::CUT::Point::oncutsurface*/)
          {
//            if((*f)->Position() == GEO::CUT::Point::oncutsurface) cout << "parent ele: " << ele->Id() << endl;

            bool edge_based_stab = false;
            bool ghost_penalty   = false;

            std::vector<int> nds_1;
            std::vector<int> nds_2;

            GEO::CUT::plain_volumecell_set vcs = (*f)->Cells();

            // how many volumecells found?
            if(vcs.size() == 2) // standard XFEM case (facet between two vcs of two neighbouring cut elements
            {
              GEO::CUT::VolumeCell* vc1 = vcs[0];
              GEO::CUT::VolumeCell* vc2 = vcs[1];


              // get the parent element
              int ele1_id = vc1->ParentElement()->Id();
              int ele2_id = vc2->ParentElement()->Id();

              // which element is the parent element
              if(ele1_id == ele->Id())
              {
                ele1 = discret.gElement(ele1_id);
                ele2 = discret.gElement(ele2_id);

                nds_1 = vc1->NodalDofSet();
                nds_2 = vc2->NodalDofSet();
              }
              else if(ele2_id == ele->Id())
              { // switch ele 1 <-> ele 2
                ele1 = discret.gElement(ele2_id);
                ele2 = discret.gElement(ele1_id);

                nds_1 = vc2->NodalDofSet();
                nds_2 = vc1->NodalDofSet();
              }
              else dserror("no element (ele1 and ele2) is the parent element!!! WHY?");

              //
              //            // check if both volume cells are outside
              //            if(vc1->Position() == GEO::CUT::Point::inside or
              //                vc2->Position() == GEO::CUT::Point::inside)
              //            {
              //              cout << "ele1 " << ele1_id << " ele2 " << ele2_id << endl;
              //              cout << "vc1_pos " << vc1->Position() << "vc2_pos " << vc2->Position() << endl;
              //              //              dserror("at least one vc of two is inside !");
              //            }
              //            // check if both volume cells are outside
              //            if(vc1->Position() == GEO::CUT::Point::inside and
              //                vc2->Position() == GEO::CUT::Point::inside)
              //            {
              //              cout << "ele1 " << ele1_id << " ele2 " << ele2_id << endl;
              //              cout << "vc1_pos " << vc1->Position() << "vc2_pos " << vc2->Position() << endl;
              //              dserror("both vcs are inside!");
              //            }



//              nds_1.clear();
//              nds_2.clear();
//              // remove this for more than one dofset!!!
//              for (int t=0;t<8; ++t)
//              {
//                nds_1.push_back(0);
//                nds_2.push_back(0);
//              }
//              ele1 = discret.gElement(ele1_id);
//              ele2 = discret.gElement(ele2_id);

              edge_based_stab = true;

              GEO::CUT::ElementHandle * eh_1   = wizard_->GetElement( ele1 );
              GEO::CUT::ElementHandle * eh_2   = wizard_->GetElement( ele2 );
              // at least one element has to be cut
              if(eh_1->IsCut() or eh_2->IsCut()) ghost_penalty   = true;


              // call evaluate routine
              callEdgeStabandGhostPenalty( edge_based_stab,
                  ghost_penalty,
                  surface,
                  ele1,
                  ele2,
                  nds_1,
                  nds_2,
                  discret,
                  strategy );

            }
            else if(vcs.size() == 1)
            {
              // only the volumecell of one element found
              GEO::CUT::VolumeCell* vc = vcs[0];

              // find out, is there another element adjacent to this facet, or is this a facet of a boundary element, or a element at the processors limit


              bool neighbor_found = false;
              int  neighbor_id = -1;

              // find out, is there another neighboring element
              findNeighborElement(ele, surface, neighbor_found, neighbor_id);

              if(neighbor_id == -1)
              { // no neighbor element found
                edge_based_stab = false;
                ghost_penalty   = false;
              }
              else
              {
                // side is cut (because of the facet)
                // side has to be stabilized
                edge_based_stab = true;
                if( e_parent->IsCut() ) ghost_penalty   = true;

                // the current element must be the parent element
                if(ele->Id() != vc->ParentElement()->Id()) dserror("the current element with ehandle must be the vc's element");

                // get the parent element
                int ele1_id = vc->ParentElement()->Id();
                int ele2_id = neighbor_id;

                nds_1 = vc->NodalDofSet();

//                nds_1.clear();
//                for(int i=0; i< discret.gElement(ele1_id)->NumNode(); i++)
//                {
//                  nds_1.push_back(0);
//                }

                nds_2.clear();
                for(int i=0; i< discret.gElement(ele2_id)->NumNode(); i++)
                {
                  nds_2.push_back(0);
                }


                ele1 = discret.gElement(ele1_id);
                ele2 = discret.gElement(ele2_id);


                // call evaluate routine
                callEdgeStabandGhostPenalty( edge_based_stab,
                    ghost_penalty,
                    surface,
                    ele1,
                    ele2,
                    nds_1,
                    nds_2,
                    discret,
                    strategy );

              }

              //TODO:
            }
            else dserror("there should be either 2 or 1 volumecell for a facet");

          }
          else if((*f)->Position() == GEO::CUT::Point::undecided)
          {
            dserror("the position of this facet is undecided, how to stabilize???");
          }
          else if((*f)->Position() == GEO::CUT::Point::oncutsurface)
          {
            cout << "the position of this facet of element " << ele->Id() << " is oncutsurface, how to stabilize??? surfid: " << surfid << endl;;
          }
          else
          {
            // facet is inside!
          }

        }



        } // loop over sides

    }
    else if(e_parent->Shape() == DRT::Element::hex20)
    {
      dserror("not yet implemented for hex20 elements");
    }
    else if(e_parent->Shape() == DRT::Element::hex27)
    {
      dserror("not yet implemented for hex27 elements");
    }
    else dserror("not yet implemented for this element type");

    // stabilize the side for all its facets
  }
  else // no element handle for current element
  {
    // find neigbor element if available
    stabilizeStandardElement( strategy,
                              ele,
                              discret);
  }


#endif
  return;
}



void XFEM::XFEM_EdgeStab::callEdgeStabandGhostPenalty( bool & edge_based_stab,
                                                            bool & ghost_penalty,
                                                            RCP<DRT::Element>       surface,
                                                            DRT::Element* ele_1,
                                                            DRT::Element* ele_2,
                                                            std::vector<int> & nds_1,
                                                            std::vector<int> & nds_2,
                                                            DRT::Discretization &   discret,
                                                            DRT::AssembleStrategy&  strategy)
{
#ifdef D_FLUID3
  //======================================================================================
  // call the internal faces stabilization routine for the current side/surface

  // call edge-based stabilization and ghost penalty
  ParameterList edgebasedparams;

  // set action for elements
  edgebasedparams.set("action","edge_based_stabilization");
  edgebasedparams.set("edge_based_stab", edge_based_stab);
  edgebasedparams.set("ghost_penalty",   ghost_penalty);


  DRT::ELEMENTS::Fluid3Boundary * side_ele = dynamic_cast<DRT::ELEMENTS::Fluid3Boundary *>( &*surface );

  DRT::ELEMENTS::Fluid3 * fele_1 = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( ele_1 );
  DRT::ELEMENTS::Fluid3 * fele_2 = dynamic_cast<DRT::ELEMENTS::Fluid3 *>( ele_2 );

  // call the egde-based routine
  DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(&*surface)->EvaluateInternalFacesUsingNeighborData( &*side_ele,
                                                                                                       fele_1,
                                                                                                       fele_2,
                                                                                                       nds_1,
                                                                                                       nds_2,
                                                                                                       edgebasedparams,
                                                                                                       discret,
                                                                                                       strategy.Systemmatrix1(),
                                                                                                       strategy.Systemmatrix2(),
                                                                                                       strategy.Systemvector1(),
                                                                                                       strategy.Systemvector2(),
                                                                                                       strategy.Systemvector3() );

#endif

  return;
}
