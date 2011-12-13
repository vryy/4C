/*!
\file xfluidfluid_timeInt.classes

\brief provides the xfluidfluid timeIntegration

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
#ifdef CCADISCRET

#include "xfluidfluid_timeInt.H"
#include "xfem_fluidwizard.H"

#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_element.H"

#include <iostream>

XFEM::XFluidFluidTimeIntegration::XFluidFluidTimeIntegration(
  const RCP<DRT::Discretization> bgdis,
  const RCP<DRT::Discretization> embdis,
  XFEM::FluidWizard              wizard,
  int                            step
  ) :
  embdis_(embdis),
  step_(step)
  {

    CreateBgNodeMaps(bgdis,wizard);
    currentbgdofmap_ =  bgdis->DofRowMap();
    return;

  } // end constructor

// -------------------------------------------------------------------
// map of standard node ids and their dof-gids in for this time step
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::CreateBgNodeMaps(const RCP<DRT::Discretization> bgdis,
                                                        XFEM::FluidWizard              wizard)
{
  const Epetra_Map* noderowmap = bgdis->NodeRowMap();
  // map of standard nodes and their dof-ids

  for (int lid=0; lid<noderowmap->NumGlobalPoints(); lid++)
  {
    int gid;
    // get global id of a node
    gid = noderowmap->GID(lid);
    // get the node
    DRT::Node * node = bgdis->gNode(gid);
    GEO::CUT::Node * n = wizard.GetNode(node->Id());
    if (n!=NULL) // xfem nodes
    {
      GEO::CUT::Point * p = n->point();
      GEO::CUT::Point::PointPosition pos = p->Position();
      if (pos==GEO::CUT::Point::outside and bgdis->NumDof(node) != 0) //std
      {
        //cout << " outside " << pos <<  " "<< node->Id() << endl;
        vector<int> gdofs = bgdis->Dof(node);
        stdnodenp_[gid] = gdofs;
      }
      else if (pos==GEO::CUT::Point::inside and  bgdis->NumDof(node) == 0) //void
      {
        //cout << " inside " <<  pos << " " << node->Id() << endl;
      }
      else if (pos==GEO::CUT::Point::inside and  bgdis->NumDof(node) != 0) //enriched
      {
        //cout << " inside enriched" <<  pos << " " << node->Id() << endl;
        vector<int> gdofs = bgdis->Dof(node);
        enrichednodenp_[gid] = gdofs;
      }
      else if (pos==GEO::CUT::Point::oncutsurface and bgdis->NumDof(node) == 0)
      {
        cout << " oncutsurface " << node->Id() << endl;
      }
      else
      {
        cout << "  hier ?! " <<  pos << " " <<  node->Id() <<  "numdof: " << bgdis->NumDof(node) << endl;
      }
    }
    else if( bgdis->NumDof(node) != 0) // no xfem node
    {
      vector<int> gdofs = bgdis->Dof(node);
      stdnodenp_[gid] = gdofs;
    }
    else
      cout << " why here? " << "node "<<  node->Id()<< " " <<  bgdis->NumDof(node) <<  endl;
  }

//     //debug output
//     for (int i=0; i<bgdis->NumMyColNodes(); ++i)
//     {
//       const DRT::Node* actnode = bgdis->lColNode(i);
//       map<int, vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
//       map<int, vector<int> >::const_iterator iter2 = enrichednoden_.find(actnode->Id());
//       map<int, vector<int> >::const_iterator iter3 = stdnodenp_.find(actnode->Id());
//       map<int, vector<int> >::const_iterator iter4 = enrichednodenp_.find(actnode->Id());
//       if (iter2 != enrichednoden_.end()) cout  << " enrichned n : " << actnode->Id() << " "  ;
//       if (iter2 == enrichednoden_.end() and iter == stdnoden_.end()) cout  << " void n :" <<  actnode->Id() << " "  ;
//       if (iter4 != enrichednodenp_.end()) cout  << " enrichned np : " << actnode->Id() << " "  ;
//       if (iter4 == enrichednodenp_.end() and iter3 == stdnodenp_.end()) cout  << " void np :" <<  actnode->Id() << " "  ;
//     }
}

// -------------------------------------------------------------------
// save the old maps
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SaveBgNodeMaps()
{
  // save the old maps and clear the maps for the new cut
  // (all maps are related to the background fluid)
  stdnoden_ = stdnodenp_;
  enrichednoden_ = enrichednodenp_;
  stdnodenp_.clear();
  enrichednodenp_.clear();
}

// -------------------------------------------------------------------
// - Save the old maps of bg nodes
// - Create new map of bg nodes
// - Gmsh-output
// -------------------------------------------------------------------
int XFEM::XFluidFluidTimeIntegration::SaveAndCreateNewBgNodeMaps(RCP<DRT::Discretization> bgdis,
                                                                  XFEM::FluidWizard        wizard)
{
  // save the old maps and clear the maps for the new cut
  // (all maps are related to the background fluid)
  SaveBgNodeMaps();

  // Create new maps
  CreateBgNodeMaps(bgdis,wizard);

  oldbgdofmap_ = currentbgdofmap_;
  currentbgdofmap_ =  bgdis->DofRowMap();
  if ((oldbgdofmap_->SameAs(*currentbgdofmap_) and (stdnoden_ == stdnodenp_) and
       (enrichednoden_ == enrichednodenp_)))
    samemaps_ = true;
  else samemaps_ = false;

  GmshOutput(bgdis);

  return samemaps_;
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorAndProjectEmbToBg(const RCP<DRT::Discretization>        bgdis,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn,
                                                                            Teuchos::RCP<Epetra_Vector>           aledispn)
{
  for (int lnid=0; lnid<bgdis->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis->lRowNode(lnid);
    map<int, vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != stdnoden_.end() and iterstnp != stdnodenp_.end()) or
        (iterstn != stdnoden_.end() and iterenp != enrichednodenp_.end()))
    {
      vector<int> gdofsn = iterstn->second;
      (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[0])] =
        (*bgstatevn)[bgstatevn->Map().LID(gdofsn[0])];
      (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[1])] =
        (*bgstatevn)[bgstatevn->Map().LID(gdofsn[1])];
      (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[2])] =
        (*bgstatevn)[bgstatevn->Map().LID(gdofsn[2])];
      (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[3])] =
        (*bgstatevn)[bgstatevn->Map().LID(gdofsn[3])];
    }
    // Project dofs from embdis to bgdis:
    // n:void -> n+1:std, n:enriched -> n+1:enriched,
    // n:enriched -> n+1: std, n:void ->  n+1:enriched
    else if (((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterstnp != stdnodenp_.end()) or
             (iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end()) or
             (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end()) or
             ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterenp != enrichednodenp_.end()))
    {
      LINALG::Matrix<3,1> bgnodecords(true);
      bgnodecords(0,0) = bgnode->X()[0];
      bgnodecords(1,0) = bgnode->X()[1];
      bgnodecords(2,0) = bgnode->X()[2];

//       std::vector<int> elementsIdsofRelevantBoxes;
//       // loop the patchboxes to find out in which patch element the xfem node is included
//       for ( std::map<int,GEO::CUT::BoundingBox>::const_iterator iter=patchboxes.begin();
//             iter!=patchboxes.end(); ++iter)
//       {
//         const GEO::CUT::BoundingBox & patchbox = iter->second;
//         double norm = 1e-5;
//         bool within = patchbox.Within(norm,bgnode->X());
//         if (within) elementsIdsofRelevantBoxes.push_back(iter->first);
//       }
//       // compute the element coordinates of backgroundflnode due to the patch discretization
//       // loop over all relevant patch boxes
//       bool insideelement;
//       size_t count = 0;
//       for (size_t box=0; box<elementsIdsofRelevantBoxes.size(); ++box)
//       {
//         // get the patch-element for the box
//         DRT::Element* pele = embdis_->gElement(elementsIdsofRelevantBoxes.at(box));
//         cout << " box: " << box << " id: "  << pele->Id() << endl;
//         LINALG::Matrix<4,1> interpolatedvec(true);
//         insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,fluidstate_vector_n);
//         if (insideelement)
//         {
//           // hier set state
//           (*statevnp)[statevnp->Map().LID(bgdis->Dof(bgnode)[0])] = interpolatedvec(0);
//           (*statevnp)[statevnp->Map().LID(bgdis->Dof(bgnode)[1])] = interpolatedvec(1);
//           (*statevnp)[statevnp->Map().LID(bgdis->Dof(bgnode)[2])] = interpolatedvec(2);
//           (*statevnp)[statevnp->Map().LID(bgdis->Dof(bgnode)[3])] = interpolatedvec(3);
//           break;
//         }
//         count ++;

      bool insideelement = false;
      int count = 0;
      // check all embedded elements to find the right one, the patch
      // boxes are not used
      for (int e=0; e<embdis_->NumMyColElements(); e++)
      {
        DRT::Element* pele = embdis_->lColElement(e);
        LINALG::Matrix<4,1> interpolatedvec(true);
        insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,*embstatevn,aledispn,embdis_);
        if (insideelement)
        {
          // hier set state
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[0])] = interpolatedvec(0);
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[1])] = interpolatedvec(1);
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[2])] = interpolatedvec(2);
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[3])] = interpolatedvec(3);
          break;
        }
        count ++;
      }
      if (count == embdis_->NumMyColElements())     // if (count == elementsIdsofRelevantBoxes.size())
      {
        // if there are any enriched values..
        if ((iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end())
            or ((iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end())))
        {
          vector<int> gdofsn = iteren->second;
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[0])] =
            (*bgstatevn)[bgstatevn->Map().LID(gdofsn[0])];
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[1])] =
            (*bgstatevn)[bgstatevn->Map().LID(gdofsn[1])];
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[2])] =
            (*bgstatevn)[bgstatevn->Map().LID(gdofsn[2])];
          (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[3])] =
            (*bgstatevn)[bgstatevn->Map().LID(gdofsn[3])];
        }
        else
        {
          cout << YELLOW_LIGHT << " Warning FFT: No patch element found for the node " << bgnode->Id() ;
          if ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterstnp != stdnodenp_.end())
            cout << YELLOW_LIGHT << " n:void -> n+1:std  " << END_COLOR<< endl;
          else if (iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end())
            cout << YELLOW_LIGHT << " n:enriched -> n+1:enriched " << END_COLOR<< endl;
          else if (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end())
            cout << YELLOW_LIGHT << " n:enriched -> n+1: std " << END_COLOR<< endl;
          else if ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterenp != enrichednodenp_.end())
            cout << YELLOW_LIGHT << " n:void ->  n+1:enriched " << END_COLOR<< endl;
        }

        //dserror("STOP");
      }
    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:std->n+1:enriched, n:enriched->n+1:void
    else if ( (iterstn == stdnoden_.end() and iteren == enrichednoden_.end()
               and iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end() ) or
              (iterstn != stdnoden_.end() and (iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end())) or
              (iterstn != stdnoden_.end() and iterenp != enrichednodenp_.end()) or
              (iteren != enrichednoden_.end() and (iterenp == enrichednodenp_.end() and iterstnp == stdnodenp_.end())) )
    {
      //cout << "do nothing" << bgnode->Id() << endl; ;
    }
    else
      cout << "warum bin ich da?! " <<  bgdis->NumDof(bgnode)   << " " <<    bgnode->Id() <<  endl;
  }
}//SetNewBgStatevectorAndProjectEmbToBg

// -------------------------------------------------------------------
// In this function:
//
// - pele (in): The element where we want to extract our values from
// - x    (in): The coordinates of the node without values (x,y,z)
// - interpolatedvec (out): The vector after interpolation
// - embstate_n (in): Source state vector
// - embeddedddisp (in): The displacement of the embedded fluid at the
//   time that it is considered as source of the interpolation
// - sourcedis (in): the source discretization where we get the values
//   from
//
//   Check whether the unknown node lies in the pele. If yes fill the
//   interpolatedvec and return true.
//
// -------------------------------------------------------------------
bool XFEM::XFluidFluidTimeIntegration::ComputeSpacialToElementCoordAndProject(DRT::Element*                       pele,
                                                                              LINALG::Matrix<3,1>&                x,
                                                                              LINALG::Matrix<4,1>&                interpolatedvec,
                                                                              Epetra_Vector                       embstate_n,
                                                                              Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                                                              Teuchos::RCP<DRT::Discretization>   sourcedis)
{
  const std::size_t numnode = pele->NumNode();
  LINALG::SerialDenseMatrix pxyze(3,numnode);
  DRT::Node** pelenodes = pele->Nodes();

  std::vector<double> myval(4);
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);
  LINALG::Matrix<3,1> xsi(true);

  bool inside = false;

  switch ( pele->Shape() )
  {
    case DRT::Element::hex8:
    {
      const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::Matrix<4,numnodes> veln(true);
      LINALG::Matrix<4,numnodes> disp;

      for (int inode = 0; inode < numnodes; ++inode)
      {
        sourcedis->Dof(pelenodes[inode],0,pgdofs);
        DRT::UTILS::ExtractMyValues(embstate_n,myval,pgdofs);
        veln(0,inode) = myval[0];
        veln(1,inode) = myval[1];
        veln(2,inode) = myval[2];
        veln(3,inode) = myval[3];

        // get the coordinates of patch element
        pxyze(0,inode) = pelenodes[inode]->X()[0];
        pxyze(1,inode) = pelenodes[inode]->X()[1];
        pxyze(2,inode) = pelenodes[inode]->X()[2];

        if (sourcedis->Name() == "xfluid") // embedded fluid
        {
          // we have to take aledispnm_ because the aledispn_ is already updated
          DRT::UTILS::ExtractMyValues(*embeddeddisp,mydisp,pgdofs);
          disp(0,inode) = mydisp[0];
          disp(1,inode) = mydisp[1];
          disp(2,inode) = mydisp[2];

          // add the current displacement to the coordinates
          pxyze(0,inode) += disp(0,inode);
          pxyze(1,inode) += disp(1,inode);
          pxyze(2,inode) += disp(2,inode);
        }
      }


      // check whether the xfemnode (the node with no values) is in the element
      LINALG::Matrix< 3,numnodes > xyzem(pxyze,true);
      GEO::CUT::Position <DRT::Element::hex8> pos(xyzem,x);
      double tol = 1e-10;
      bool insideelement = pos.ComputeTol(tol);
      //bool insideelement = pos.Compute();

      // von ursula
      // bool in  = GEO::currentToVolumeElementCoordinates(DRT::Element::hex8, pxyze, x, xsi);
      // in  = GEO::checkPositionWithinElementParameterSpace(xsi, DRT::Element::hex8);
      // cout << " ursula " << in << endl;

      if (insideelement)
      {
        // get the coordinates of x in element coordinates of patch element pele (xsi)
        xsi = pos.LocalCoordinates();
        // evaluate shape function
        LINALG::SerialDenseVector shp(numnodes);
        DRT::UTILS::shape_function_3D( shp, xsi(0,0), xsi(1,0), xsi(2,0), DRT::Element::hex8 );
        // Interpolate
        for (int inode = 0; inode < numnodes; ++inode){
          for (std::size_t isd = 0; isd < 4; ++isd){
            interpolatedvec(isd) += veln(isd,inode)*shp(inode);
          }
        }
        inside = true;
        break;
      }
      else
      {
        inside = false;
        break;
      }
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      dserror("No support for hex20 and hex27!");
      break;
    }
    default:
    {
      dserror("Element-type not supported here!");
    }
  }
  return inside;
}//ComputeSpacialToElementCoordAndProject

// -------------------------------------------------------------------
// In this function:
//
// - pele (in): The element where we want to extract our values from
// - x    (in): The coordinates of the node without values (x,y,z)
// - interpolatedvec (out): The vector after interpolation
// - embstate_n (in): Source state vector
// - embeddedddisp (in): The displacement of the embedded fluid at the
//   time that it is considered as source of the interpolation
// - sourcedis (in): the source discretization where we get the values
//   from
//
//   Check whether the unknown node lies in the pele. If yes fill the
//   interpolatedvec and return true.
//
// -------------------------------------------------------------------
bool XFEM::XFluidFluidTimeIntegration::ComputeSpacialToElementCoordAndProject2(DRT::Element*                       pele,
                                                                              LINALG::Matrix<3,1>&                x,
                                                                              LINALG::Matrix<4,1>&                interpolatedvec,
                                                                              Epetra_Vector                       embstate_n,
                                                                              Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                                                              Teuchos::RCP<DRT::Discretization>   sourcedis)
{
  const std::size_t numnode = pele->NumNode();
  LINALG::SerialDenseMatrix pxyze(3,numnode);
  DRT::Node** pelenodes = pele->Nodes();

  std::vector<double> myval(4);
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);
  LINALG::Matrix<3,1> xsi(true);

  bool inside = false;

  switch ( pele->Shape() )
  {
    case DRT::Element::hex8:
    {
//      cout <<" bg ele id" << pele->Id() << endl;
      const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::Matrix<4,numnodes> veln(true);
      LINALG::Matrix<4,numnodes> disp;

      for (int inode = 0; inode < numnodes; ++inode)
      {
        // the enriched nodes are also included
        if (sourcedis->NumDof(pelenodes[inode]) > 0)
        {
          sourcedis->Dof(pelenodes[inode],0,pgdofs);
//           cout << sourcedis->NumDof(pelenodes[inode]) << endl;
//           cout << pgdofs.at(0) << " " << pgdofs.at(1) << " " << pgdofs.at(2) << " "<< pgdofs.at(3)<<  endl;
          DRT::UTILS::ExtractMyValues(embstate_n,myval,pgdofs);
          veln(0,inode) = myval[0];
          veln(1,inode) = myval[1];
          veln(2,inode) = myval[2];
          veln(3,inode) = myval[3];

          // get the coordinates of patch element
          pxyze(0,inode) = pelenodes[inode]->X()[0];
          pxyze(1,inode) = pelenodes[inode]->X()[1];
          pxyze(2,inode) = pelenodes[inode]->X()[2];
        }
      }


      // check whether the node with no values is in the element
      LINALG::Matrix< 3,numnodes > xyzem(pxyze,true);
      GEO::CUT::Position <DRT::Element::hex8> pos(xyzem,x);
      double tol = 1e-10;
      bool insideelement = pos.ComputeTol(tol);

      if (insideelement)
      {
        // get the coordinates of x in element coordinates of patch element pele (xsi)
        xsi = pos.LocalCoordinates();
        // evaluate shape function
        LINALG::SerialDenseVector shp(numnodes);
        DRT::UTILS::shape_function_3D( shp, xsi(0,0), xsi(1,0), xsi(2,0), DRT::Element::hex8 );
        // Interpolate
        for (int inode = 0; inode < numnodes; ++inode){
          for (std::size_t isd = 0; isd < 4; ++isd){
            interpolatedvec(isd) += veln(isd,inode)*shp(inode);
          }
        }
        inside = true;
        break;
      }
      else
      {
        inside = false;
        break;
      }
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      dserror("No support for hex20 and hex27!");
      break;
    }
    default:
    {
      dserror("Element-type not supported here!");
    }
  }
  return inside;
}//ComputeSpacialToElementCoordAndProject2

// -------------------------------------------------------------------
// Needed for fluid-fluid-fsi
//
// In this function:
// loop over all embedded fluid nodes and find the embedded element
// where this node was located before we solved the ale again.(the
// displacement is saved at aledispnpoldstate_). If found interpolate
// the values on the new embedded fluid vector. If not found we have
// to look for a matching element in the background fluid. Right now
// we just take the old values of the embedded mesh: TODO
//
// statevbg_n (in): source state vector from gbackground fluid
// statevemb_n (in): source state vector from embedded fluid
// statevembnew_n (out)
// aledispnp (in): displacement of ale at the current time (time of
//               interpolation)
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewEmbStatevector(const RCP<DRT::Discretization>        bgdis,
                                                            Teuchos::RCP<Epetra_Vector>    statevbg_n,
                                                            Epetra_Vector                  statevemb_n,
                                                            Teuchos::RCP<Epetra_Vector>    statevembnew_n,
                                                            Teuchos::RCP<Epetra_Vector>    aledispnp,
                                                            Teuchos::RCP<Epetra_Vector>    aledispnpoldstate)
{
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);

  // Gmsh----------------------------------------------------------
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("emb_element_node_id", 0, 0, 0, bgdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    {
      // draw embedded elements with associated gid at the old  position
      gmshfilecontent << "View \" " << "emb Element(old)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        vector<double> myolddisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnpoldstate, myolddisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + myolddisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + myolddisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + myolddisp[2+(inode*4)];

          mapofnodepos[pelenodes[inode]->Id()] = inodepos;
        }
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, mapofnodepos, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }

    {
      // draw embedded elements with associated gid at the current position
      gmshfilecontent << "View \" " << "emb Element(new)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnp, mydisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + mydisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + mydisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + mydisp[2+(inode*4)];

          mapofnodepos[pelenodes[inode]->Id()] = inodepos;
        }
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, mapofnodepos, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }

    {
      // draw embedded nodes with associated gid at the current position
      gmshfilecontent << "View \" " << "emb Node(new)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnp, mydisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + mydisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + mydisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + mydisp[2+(inode*4)];

          IO::GMSH::cellWithScalarToStream(DRT::Element::point1, pelenodes[inode]->Id(), inodepos, gmshfilecontent);
        }
      };
      gmshfilecontent << "};\n";
    }
  }
  gmshfilecontent.close();
  // --------------------------------------------------------------

  for (int lnid=0; lnid<embdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* embnode = embdis_->lRowNode(lnid);

    DRT::UTILS::ConditionSelector conds(*embdis_, "FSICoupling");
    if (conds.ContainsNode(embnode->Id())==false)//if no fsi node
    {
      //cout << "unknown node id " << embnode->Id() << endl;
      embdis_->Dof(embnode,0,pgdofs);

      DRT::UTILS::ExtractMyValues(*aledispnp,mydisp,pgdofs);

      // the coordinates of the unknown node
      LINALG::Matrix<3,1> embnodecords(true);
      embnodecords(0,0) = embnode->X()[0] + mydisp[0];
      embnodecords(1,0) = embnode->X()[1] + mydisp[1];
      embnodecords(2,0) = embnode->X()[2] + mydisp[2];

      bool insideelement = false;
      int count = 0;
      // check all embedded elements to find the right one

      for (int e=0; e<embdis_->NumMyColElements(); e++)
      {
        DRT::Element* pele = embdis_->lColElement(e);
        LINALG::Matrix<4,1> interpolatedvec(true);
        //cout << "pele id " << pele->Id() << endl;

        insideelement = ComputeSpacialToElementCoordAndProject(pele,embnodecords,interpolatedvec,statevemb_n,aledispnpoldstate,embdis_);
        if (insideelement)
        {
          // here set state
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[0])] = interpolatedvec(0);
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[1])] = interpolatedvec(1);
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[2])] = interpolatedvec(2);
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[3])] = interpolatedvec(3);
          break;
        }
        count ++;
      }
      if (count == embdis_->NumMyColElements()) //if no embedded elements found
      {
        cout << "no embedded element found for: " <<  embnode->Id() << endl;
         count = 0;
//         // check all background elements to find the right one
         for (int e=0; e<bgdis->NumMyColElements(); e++)
         {
            DRT::Element* bgele = bgdis->lColElement(e);
            LINALG::Matrix<4,1> interpolatedvec(true);

            insideelement = ComputeSpacialToElementCoordAndProject2(bgele,embnodecords,interpolatedvec,*statevbg_n,aledispnpoldstate,bgdis);
//            if (insideelement)
//            {
//              // here set state
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[0])] = interpolatedvec(0);
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[1])] = interpolatedvec(1);
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[2])] = interpolatedvec(2);
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[3])] = interpolatedvec(3);
//              break;
//            }
//            count ++;
        }
//          if (count == bgdis->NumMyColElements())
//            cout << "no bg elements found for: " << embnode->Id() << endl;
        }
    }
  }

}//SetNewEmbStatevector

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::GmshOutput(const RCP<DRT::Discretization>        bgdis)
{
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("std_enriched_map", step_, 30, 0, bgdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << "std/enriched/void n\" {\n";
    for (int i=0; i<bgdis->NumMyColNodes(); ++i)
    {
      int kind = 0;
      const DRT::Node* actnode = bgdis->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      map<int, vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
      map<int, vector<int> >::const_iterator iteren = enrichednoden_.find(actnode->Id());
      if (iter != stdnoden_.end()) kind = 1;//std
      if (iteren != enrichednoden_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "std/enriched/void n+1\" {\n";
    for (int i=0; i<bgdis->NumMyColNodes(); ++i)
    {
      int kind = 0;
      const DRT::Node* actnode = bgdis->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      map<int, vector<int> >::const_iterator iter = stdnodenp_.find(actnode->Id());
      map<int, vector<int> >::const_iterator iteren = enrichednodenp_.find(actnode->Id());
      if (iter != stdnodenp_.end()) kind = 1;//std
      if (iteren != enrichednodenp_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
}

// // -------------------------------------------------------------------
// //
// // -------------------------------------------------------------------
// void FLD::XFluidFluid::CreatePatchBoxes(std::map<int, GEO::CUT::BoundingBox> & patchboxes)
// {
//   // get column version of the displacement vector
//   Teuchos::RCP<const Epetra_Vector> col_embfluiddisp =
//     DRT::UTILS::GetColVersionOfRowVector(embdis_, aledispnm_);

//   // Map of all boxes of embedded fluid discretization
//   for (int pele=0; pele<embdis_->NumMyColElements(); ++pele)
//   {
//     const DRT::Element* actpele = embdis_->lColElement(pele);
//     const DRT::Node*const* pelenodes = actpele->Nodes();

//     vector<int> lm;
//     vector<int> lmowner;
//     vector<int> lmstride;
//     actpele->LocationVector(*embdis_, lm, lmowner, lmstride);

//     vector<double> mydisp(lm.size());
//     DRT::UTILS::ExtractMyValues(*col_embfluiddisp, mydisp, lm);

//     GEO::CUT::BoundingBox patchbox;
//     //patchboxes
//     for (int pnode = 0; pnode < actpele->NumNode(); ++pnode)
//     {
//       // the coordinates of the actuall node
//       LINALG::Matrix<3,1> pnodepos(true);
//       pnodepos(0,0) = pelenodes[pnode]->X()[0] + mydisp[0+(pnode*4)];
//       pnodepos(1,0) = pelenodes[pnode]->X()[1] + mydisp[1+(pnode*4)];
//       pnodepos(2,0) = pelenodes[pnode]->X()[2] + mydisp[2+(pnode*4)];

//       // fill the patchbox
//       patchbox.AddPoint(pnodepos);
//     }
// //     cout << actpele->Id() << endl;
// //     patchbox.Print();
//     patchboxes[actpele->Id()] = patchbox;
//   }
// }

#endif // CCADISCRET
