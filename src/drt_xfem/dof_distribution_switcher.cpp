/*!
\file dof_distribution_switcher.cpp

\brief provides the dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include "dof_management.H"
#include "dof_distribution_switcher.H"
#include "dofkey.H"
#include "interfacexfsi.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_element.H"
#include <iostream>


void XFEM::DofDistributionSwitcher::mapVectorToNewDofDistribution(
    RCP<Epetra_Vector>&    vector
) const
{
  // create new vector with new number of dofs
  const RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);

  if (vector == null)
  {
#ifdef DEBUG
    std::cout << "  created new vector with all zeros" << endl;
#endif
  }
  else
  {
#ifdef DEBUG
    bool completely_unchanged = true;
#endif
    const RCP<Epetra_Vector> oldVector = vector;
    const Epetra_BlockMap& oldmap = oldVector->Map();
//    std::cout << "olddofrowmap_" << endl;
//    std::cout << (olddofrowmap_) << endl;
//    std::cout << "newdofrowmap_" << endl;
//    std::cout << (newdofrowmap_) << endl;

    if (not oldmap.SameAs(olddofrowmap_)) dserror("bug!");

    // step 1: find predecessor of new nodal dofkey
    for (map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.begin();
                                                     newdof != newNodalDofDistrib_.end();
                                                     ++newdof)
    {
      const DofKey<onNode> newdofkey = newdof->first;
      const int newdofpos = newdof->second;

      map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
      if (olddof != oldNodalDofDistrib_.end()) // if dofkey has existed before, use old value
      {
        const DofKey<onNode> olddofkey = olddof->first;
        const int olddofpos = olddof->second;
        //cout << newdofkey.toString() << " -> init to old value" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
      }
      else // if dofkey has not been existed before, check for other dofs on the dofkeys node
      {
        //const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();

        // initialize to zero
        (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
      }
    }

    // step 2: find successor of old nodal dofkey to sum up values
    for (map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.begin();
                                                     olddof != oldNodalDofDistrib_.end();
                                                     ++olddof)
    {
      const DofKey<XFEM::onNode> olddofkey = olddof->first;
      const int olddofpos = olddof->second;
      const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();

      // try to find successor
      map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
      if (newdof == newNodalDofDistrib_.end())  // if no successor found (was handled already in step 1)
      {
        // try to find another usefull value
        // current assumption: there is only one type of enrichment per node
        // no overlapping enrichments allowed for now
        const int nodegid = olddofkey.getGid();
        const XFEM::Enrichment oldenr(olddofkey.getFieldEnr().getEnrichment());

        // create alternative dofkey
        XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, oldphysvar, dofman_));

        if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
        {
          // find dof position of alternative key
          const XFEM::FieldEnr altfieldenr(olddofkey.getFieldEnr().getField(), altenr);
          const DofKey<onNode> altdofkey(nodegid, altfieldenr);
          const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;

          //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
          if (newdofpos < 0)
          {
            std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
            std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
            dserror("bug!");
          }

          // add old value to already existing values
          //(*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofrowmap_.LID(olddofpos)];
          (*newVector)[newdofrowmap_.LID(newdofpos)] += (*oldVector)[olddofrowmap_.LID(olddofpos)];
#ifdef DEBUG
          completely_unchanged = false;
#endif
        }
        else // if not alternative is found
        {
          // this can only happen in the void enrichment case and in that case,
          // the dof value is zero anyway, which coincides with the fact that we have no place,
          // where we could store it ;-)
        }
      }
      else
      {
        // do nothing, this case was handled in step 1
      }
    }

#if 0
    // remark: Only if DLM condensation is turned off, there will be element dofs visible in the global
    //         DofManager. Of course, irrespective of the DLM condensation setting, there will be
    //         element dofs on the element level (ElementDofManager), if a stress-based approach is
    //         chosen to handle the boundary terms.
    //         So, if there are element dofs, they are not known here, since 'elementalDofs_' in
    //         fillDofRowDistributionMaps() only knows the different field enrichments for a specific
    //         element, that is
    //
    //         e.g.
    //         element 103 has 28 element dofs ((6 x stess + 1 x disc pressure) x element ansatz (e.g. quad4))
    //
    //         but 'elementalDofs_' contains only 7 "dof-type entries"
    //         <ele 103, Tauxx + quad4>
    //         <ele 103, Tauxy + quad4>
    //         ...
    //         <ele 103, DiscPres + quad4>
    //
    //         That is, element dofs cannot be switched correctly here, since not all element dofs
    //         are known. They only exist on the element level. They are generated in
    //         ElementDofManager::ComputeDependentInfo()
    //                                                                                   henke 03/10
    dserror("This will not work properly! Read comment above!");

    // step 3: find predecessor of new elemental dofkey
    for (map<DofKey<onElem>, DofGID>::const_iterator newdof = newElementalDofDistrib_.begin();
                                                     newdof != newElementalDofDistrib_.end();
                                                     ++newdof)
    {
      const DofKey<XFEM::onElem> newdofkey = newdof->first;
      const int newdofpos = newdof->second;

      map<DofKey<onElem>, DofGID>::const_iterator olddof = oldElementalDofDistrib_.find(newdofkey);
      if (olddof != oldElementalDofDistrib_.end())  // if dofkey has existed before, use old value
      {
        const DofKey<XFEM::onElem> olddofkey = olddof->first;
        const int olddofpos = olddof->second;
        //cout << "init to old value" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
      }
      else // if dofkey has not been existed before, initialize to zero
      {
        //cout << "init to zero" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
        completely_unchanged = false;
      }
    }
#endif

#if 0
    if (completely_unchanged)
      cout << "completely unchanged vector" << endl;
    else
      cout << "modified vector" << endl;
#endif
  }


  // set vector to zero or initialized vector
  vector = newVector;
}


void XFEM::DofDistributionSwitcher::extractDofKeysForInitialization(
    std::map<int, set<XFEM::FieldEnr> >& unknownFieldEnr
) const
{
  unknownFieldEnr.clear();

  // step 1: find predecessor of new nodal dofkey
  for (map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.begin();
                                                   newdof != newNodalDofDistrib_.end();
                                                   ++newdof)
  {
    const DofKey<onNode> newdofkey = newdof->first;
    map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
    if (olddof != oldNodalDofDistrib_.end()) // if dofkey has existed before, use old value
    {

    }
    else // if dofkey has not been existed before
    {
      unknownFieldEnr[newdofkey.getGid()].insert(newdofkey.getFieldEnr());
    }
  }

  // step 2: find successor of old nodal dofkey to sum up values
  for (map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.begin();
                                                   olddof != oldNodalDofDistrib_.end();
                                                   ++olddof)
  {
    const DofKey<XFEM::onNode> olddofkey = olddof->first;

    // try to find successor
    map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
    if (newdof == newNodalDofDistrib_.end())  // if no successor found (was handled already in step 1)
    {
      // try to find another usefull value
      // current assumption: there is only one type of enrichment per node
      // no overlapping enrichments allowed for now
      const int nodegid = olddofkey.getGid();

      // create alternative dofkey
      XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, olddofkey.getFieldEnr().getField(), dofman_));

      if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
      {
        // find dof position of alternative key
        const DofKey<onNode> altdofkey(nodegid, XFEM::FieldEnr(olddofkey.getFieldEnr().getField(), altenr));
        const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;

        //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
        if (newdofpos < 0)
        {
          std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
          std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
          dserror("bug!");
        }
//        unknownDofKeys.erase(altdofkey);
        unknownFieldEnr.erase(altdofkey.getGid());
      }
      else // if not alternative is found
      {
        // this can only happen in the void enrichment case and in that case,
        // the dof value is zero anyway, which coincides with the fact that we have no place,
        // where we could store it ;-)
      }
    }
    else
    {
      // do nothing, this case was handled in step 1
    }
  }

}
//-------------------------------------------------------------------
//-----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::extractDofKeysForInitializationFluidXFluid(
  std::map<int, set<XFEM::FieldEnr> >& unknownFieldEnr,
  std::map<int,int>                    nodelabeln
) const
{
  unknownFieldEnr.clear();

  // step 1: find predecessor of new nodal dofkey
  for (map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.begin();
                                                   newdof != newNodalDofDistrib_.end();
                                                   ++newdof)
  {
    const DofKey<onNode> newdofkey = newdof->first;
    map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
    if (olddof != oldNodalDofDistrib_.end()) // if dofkey has existed before, use old value
    {

    }
    else // if dofkey has not been existed before
    {
      unknownFieldEnr[newdofkey.getGid()].insert(newdofkey.getFieldEnr());
    }
  }

  // step 2: find successor of old nodal dofkey to sum up values
  for (map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.begin();
                                                   olddof != oldNodalDofDistrib_.end();
                                                   ++olddof)
  {
    const DofKey<XFEM::onNode> olddofkey = olddof->first;

    // was in xfluid? yes=0
    map<int,int>::const_iterator iter = nodelabeln.find(olddofkey.getGid());

    // try to find successor. If the node was in xfluid in the last time step skip.
    map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);

    if (newdof == newNodalDofDistrib_.end() and iter->second != 1)  // if no successor found (was handled already in step 1)
    {
      // try to find another usefull value
      // current assumption: there is only one type of enrichment per node
      // no overlapping enrichments allowed for now
      const int nodegid = olddofkey.getGid();

      // create alternative dofkey
      XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, olddofkey.getFieldEnr().getField(), dofman_));

      if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
      {
        // find dof position of alternative key
        const DofKey<onNode> altdofkey(nodegid, XFEM::FieldEnr(olddofkey.getFieldEnr().getField(), altenr));
        const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;

        //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
        if (newdofpos < 0)
        {
          std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
          std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
          dserror("bug!");
        }
//        unknownDofKeys.erase(altdofkey);
        unknownFieldEnr.erase(altdofkey.getGid());
      }
      else // if not alternative is found
      {
        // this can only happen in the void enrichment case and in that case,
        // the dof value is zero anyway, which coincides with the fact that we have no place,
        // where we could store it ;-)
      }
    }
    else
    {
      // do nothing, this case was handled in step 1
    }
  }

}
//-------------------------------------------------------------------
//-----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::extrapolateOldTimeStepValues(
    const RCP<DRT::Discretization>    bdis,
    const std::map<int,LINALG::Matrix<3,1> >&  cutterposn,
    const RCP<const Epetra_Vector>    ivector,
    const RCP<Epetra_Vector>          state_vector
) const
{

  // Achtung: extrapolation should be from fluid material velocities along the interface,
  // not from the interface velocity.
  // this works here, because for moving impermeable walls both velocities are equal.
  // for fluid-fluid coupling, this will fail!!!

  // for now, only convex structures are tested

  for (std::map<int, set<XFEM::FieldEnr> >::const_iterator nodaldofs = unknownFieldEnr_.begin();
      nodaldofs != unknownFieldEnr_.end();
      nodaldofs++)
  {
    const int nodeGid = nodaldofs->first;
    const set<XFEM::FieldEnr> fieldenrset = nodaldofs->second;
    const DRT::Node* node = ih_->xfemdis()->gNode(nodeGid);

    // generalize, such that  nearest object is returned with respect to specific interface
    GEO::NearestObject nearestobject;
    ih_->PositionWithinConditionN(LINALG::Matrix<3,1>(node->X()), nearestobject);

    if (nearestobject.getObjectType() == GEO::SURFACE_OBJECT)
    {
      const DRT::Element* bele = bdis->gElement(nearestobject.getSurfaceId());

      const std::size_t nsd = 3;
      const std::size_t numnode_boundary = bele->NumNode();

      LINALG::Matrix<nsd,9> veln;
      const DRT::Node*const* nodes = bele->Nodes();
      {
        std::vector<double> myval(nsd);
        std::vector<int> gdofs(nsd);
        for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
        {
          bdis->Dof(nodes[inode],0,gdofs);
          DRT::UTILS::ExtractMyValues(*ivector,myval,gdofs);
          veln(0,inode) = myval[0];
          veln(1,inode) = myval[1];
          veln(2,inode) = myval[2];
        }
      }

      LINALG::Matrix<2,1> xsiB(true);
      GEO::CurrentToSurfaceElementCoordinates(
          bele->Shape(),
          GEO::getCurrentNodalPositions(bele, cutterposn),
          nearestobject.getPhysCoord(),
          xsiB);
      //cout << "xsiB: " << xsiB << endl;

      LINALG::SerialDenseVector shp(DRT::UTILS::getNumberOfElementNodes(bele->Shape()));
      //static LINALG::Matrix<numnodefix_boundary,1> funct_boundary;
      DRT::UTILS::shape_function_2D(shp, xsiB(0),xsiB(1),bele->Shape());

      LINALG::Matrix<nsd,1> v(true);
      for (std::size_t iparam = 0; iparam < numnode_boundary; ++iparam)
        for (std::size_t isd = 0; isd < nsd; ++isd)
          v(isd) += veln(isd,iparam)*shp(iparam);

//      cout << "v = " << v << endl;


      for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
          fieldenr != fieldenrset.end();
          ++fieldenr)
      {
        const DofKey<onNode> newdofkey(nodeGid, *fieldenr);
        map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(newdofkey);
        const bool found = (newdof != newNodalDofDistrib_.end());
        if (not found) dserror("bug");
        const int newdofpos = newdof->second;
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(0);
        if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(1);
        if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(2);
      }


    }
    else if (nearestobject.getObjectType() == GEO::LINE_OBJECT)
    {
      // TODO: needs closer attention
      cout << "Line: " << endl;
      cout << nearestobject.getLineId() << endl;

      const DRT::Element* bele = bdis->gElement(nearestobject.getLineId());

      bele->Print(cout);
      cout << endl;

      const std::size_t nsd = 3;
      const std::size_t numnode_boundary = bele->NumNode();

      LINALG::Matrix<nsd,9> veln;
      const DRT::Node*const* nodes = bele->Nodes();
      {
        std::vector<double> myval(nsd);
        std::vector<int> gdofs(nsd);
        for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
        {
          bdis->Dof(nodes[inode],0,gdofs);
          DRT::UTILS::ExtractMyValues(*ivector,myval,gdofs);
          veln(0,inode) = myval[0];
          veln(1,inode) = myval[1];
          veln(2,inode) = myval[2];
        }
      }

      LINALG::Matrix<2,1> xsiB(true);
      GEO::CurrentToSurfaceElementCoordinates(
          bele->Shape(),
          GEO::getCurrentNodalPositions(bele, cutterposn),
          nearestobject.getPhysCoord(),
          xsiB);
      cout << "xsiB: " << xsiB << endl;

      LINALG::SerialDenseVector shp(DRT::UTILS::getNumberOfElementNodes(bele->Shape()));
      //static LINALG::Matrix<numnodefix_boundary,1> funct_boundary;
      DRT::UTILS::shape_function_2D(shp, xsiB(0),xsiB(1),bele->Shape());

      LINALG::Matrix<nsd,1> v(true);
      for (std::size_t iparam = 0; iparam < numnode_boundary; ++iparam)
        for (std::size_t isd = 0; isd < nsd; ++isd)
          v(isd) += veln(isd,iparam)*shp(iparam);

//      cout << "v = " << v << endl;


      for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
          fieldenr != fieldenrset.end();
          ++fieldenr)
      {
        const DofKey<onNode> newdofkey(nodeGid, *fieldenr);
        map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(newdofkey);
        const bool found = (newdof != newNodalDofDistrib_.end());
        if (not found) dserror("bug");
        const int newdofpos = newdof->second;
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(0);
        if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(1);
        if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(2);
      }


    }
    else if (nearestobject.getObjectType() == GEO::NODE_OBJECT)
    {
      cout << "Node: " << endl;
      cout << nearestobject.getNodeId() << endl;

      const DRT::Element* bele = bdis->gElement(nearestobject.getNodeId());

      bele->Print(cout);

      const std::size_t nsd = 3;
      const std::size_t numnode_boundary = bele->NumNode();

      LINALG::Matrix<nsd,9> veln;
      const DRT::Node*const* nodes = bele->Nodes();
      {
        std::vector<double> myval(nsd);
        std::vector<int> gdofs(nsd);
        for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
        {
          bdis->Dof(nodes[inode],0,gdofs);
          DRT::UTILS::ExtractMyValues(*ivector,myval,gdofs);
          veln(0,inode) = myval[0];
          veln(1,inode) = myval[1];
          veln(2,inode) = myval[2];
        }
      }

      LINALG::Matrix<2,1> xsiB(true);
      GEO::CurrentToSurfaceElementCoordinates(
          bele->Shape(),
          GEO::getCurrentNodalPositions(bele, cutterposn),
          nearestobject.getPhysCoord(),
          xsiB);
      //cout << "xsiB: " << xsiB << endl;

      LINALG::SerialDenseVector shp(DRT::UTILS::getNumberOfElementNodes(bele->Shape()));
      //static LINALG::Matrix<numnodefix_boundary,1> funct_boundary;
      DRT::UTILS::shape_function_2D(shp, xsiB(0),xsiB(1),bele->Shape());

      LINALG::Matrix<nsd,1> v(true);
      for (std::size_t iparam = 0; iparam < numnode_boundary; ++iparam)
        for (std::size_t isd = 0; isd < nsd; ++isd)
          v(isd) += veln(isd,iparam)*shp(iparam);

      //      cout << "v = " << v << endl;


      for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
          fieldenr != fieldenrset.end();
          ++fieldenr)
      {
        const DofKey<onNode> newdofkey(nodeGid, *fieldenr);
        map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(newdofkey);
        const bool found = (newdof != newNodalDofDistrib_.end());
        if (not found) dserror("bug");
        const int newdofpos = newdof->second;
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(0);
        if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(1);
        if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          (*state_vector)[newdofrowmap_.LID(newdofpos)] = v(2);
      }


    }
    else
    {
      // create file
      const std::string filename = "extrapolation_error.pos";
      std::ofstream gmshfilecontent(filename.c_str());

      // write fluid and boundary mesh
      IO::GMSH::disToStream("Fluid", 0.0, ih_->xfemdis(), gmshfilecontent);
      IO::GMSH::disToStream("Cutter", 1.0, bdis, cutterposn, gmshfilecontent);

      // translate test position to parameter for gmsh
      LINALG::Matrix<3,1> gmsh_pos(true);
      gmsh_pos(0,0) = node->X()[0];
      gmsh_pos(1,0) = node->X()[1];
      gmsh_pos(2,0) = node->X()[2];

      // write test position
      gmshfilecontent << "View \" " << "Problemkandidat" << " \" {\n";
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, 1.0, gmsh_pos, gmshfilecontent);
      gmshfilecontent << "};\n";
      gmshfilecontent.close();

      dserror("wieso bin ich hier?");
    }

  }

}

//----------------------------------------------------------------------
//-------------------------------------------------------------------
void XFEM::DofDistributionSwitcher::projectEmbeddedtoBackgroundfluid(
    const RCP<DRT::Discretization>                patchdis,
    const std::map<int,GEO::CUT::BoundingBox>&    patchboxes,
    const RCP<const Epetra_Vector>                fluidstate_vector_n,
    const RCP<Epetra_Vector>                      xfluidstate_vector_n,
    const RCP<const Epetra_Vector>                fluiddispnm,
    int                                           step
) const
{
  cout << "number of unknown nodes: " << unknownFieldEnr_.size() << endl;

  //-------------------------------------------
  // round robin

  // get number of processors and the current processors id
  int numproc=ih_->xfemdis()->Comm().NumProc();

#ifdef PARALLEL
  // create an exporter for point to point comunication
  DRT::Exporter exporter(ih_->xfemdis()->Comm());

  // necessary variables
  MPI_Request request;
#endif

  // define send and receive blocks
  vector<char> sblock;
  vector<char> rblock;

  for (int np=0;np<numproc+1;++np)
  {
    // in the first step, we cannot receive anything
    if (np > 0)
    {
#ifdef PARALLEL
      ReceiveBlock(rblock,exporter,request);
#else
      rblock = sblock;
#endif

      // Unpack info from the receive block from the last proc
      //UnpackValues(state_vector,unknownFieldEnr_,rblock);
    }

    if (np < numproc)
    {
      std::map<int,int> mapforoutput;

      // -----------------------
      // do what we wanted to do
      // loop over the xfluid unknown dofs
      for (std::map<int, set<XFEM::FieldEnr> >::const_iterator nodaldofs = unknownFieldEnr_.begin();
           nodaldofs != unknownFieldEnr_.end();
           nodaldofs++)
      {
        std::vector<int> elementsIdsofRelevantBoxes;

        const int nodeGid = nodaldofs->first;
        const DRT::Node* backgroundflnode = ih_->xfemdis()->gNode(nodeGid);
        // coordinates of backgroundfluid node
        LINALG::Matrix<3,1>   bgnodecords(true);
        bgnodecords(0,0) = backgroundflnode->X()[0];
        bgnodecords(1,0) = backgroundflnode->X()[1];
        bgnodecords(2,0) = backgroundflnode->X()[2];

        // loop the patchboxes to find out in which patch element the xfem node is included
        for ( std::map<int,GEO::CUT::BoundingBox>::const_iterator iter=patchboxes.begin();
              iter!=patchboxes.end(); ++iter)
        {
          const GEO::CUT::BoundingBox & patchbox = iter->second;
          // tolerance = boxtolerance(1e-7)*norm
          double norm = 1e-9;
          bool within = patchbox.Within(norm,backgroundflnode->X());

//        const double * x = backgroundflnode->X();
//        cout <<  " within " <<  within << " Box " << " ele " << iter->first << " xnode id" << backgroundflnode->Id()
//             << " (" << x[0] << " " << x[1] << " " << x[2] << ") ";
//         patchbox.Print();

          mapforoutput[iter->first] = within;
          if (within) elementsIdsofRelevantBoxes.push_back(iter->first);
        }

        // compute the element coordinates of backgroundflnode due to the patch discretization
        // loop over all patch boxes
        bool insideelement;
        for (size_t box=0; box<elementsIdsofRelevantBoxes.size(); ++box)
        {
          // get the patch-element for the box
          const DRT::Element* pele = patchdis->gElement(elementsIdsofRelevantBoxes.at(box));

           insideelement = ComputeSpacialToElementCoordAndProject(pele,patchdis,bgnodecords,
                                                                  fluidstate_vector_n,xfluidstate_vector_n,
                                                                  nodaldofs,fluiddispnm);
          if (insideelement)
          {
            break;
          }
        }
      }

      GmshOutput(step,patchdis,mapforoutput);

      // Pack UnknownFieldEnr and state_vector into block to send
      DRT::PackBuffer data;
      PackValuestoSend(xfluidstate_vector_n,unknownFieldEnr_,data);
      data.StartPacking();
      PackValuestoSend(xfluidstate_vector_n,unknownFieldEnr_,data);
      swap( sblock, data() );

#ifdef PARALLEL
      SendBlock(sblock,exporter,request);
#endif
      }
  }
}

//----------------------------------------------------------------
//-----------------------------------------------------------------
bool XFEM::DofDistributionSwitcher::ComputeSpacialToElementCoordAndProject(
  const DRT::Element*                                  pele,  // patch element
  const RCP<DRT::Discretization>                       patchdis, // patch dis
  const LINALG::Matrix<3,1>&                           x,     // background node's coordinates (x,y,z)
  const RCP<const Epetra_Vector>                       fluidstate_vector_n,
  const RCP<Epetra_Vector>                             xfluidstate_vector_n,
  std::map<int, set<XFEM::FieldEnr> >::const_iterator  unknownnodaldof,
  const RCP<const Epetra_Vector>                       fluiddispnm
  ) const
{
  // if (pele->ElementType() != DRT::Element::hex8) dserror("fluidxfluidInterpolation just inplemented for Hex8!");
  const std::size_t nsd = 4;
  LINALG::Matrix<3,1> xsi(true);
  std::vector<int> gdofs(nsd);

  const std::size_t numnode = pele->NumNode();
  LINALG::SerialDenseMatrix xyze(3,numnode);
  const DRT::Node*const* pelenodes = pele->Nodes();

  std::vector<double> myval(nsd);
  std::vector<double> mydispnm(nsd);
  LINALG::Matrix<nsd,8> veln;
  LINALG::Matrix<nsd,8> dispnm;

  for (std::size_t inode = 0; inode < numnode; ++inode)
  {
    patchdis->Dof(pelenodes[inode],0,gdofs);
    DRT::UTILS::ExtractMyValues(*fluidstate_vector_n,myval,gdofs);
    veln(0,inode) = myval[0];
    veln(1,inode) = myval[1];
    veln(2,inode) = myval[2];

    DRT::UTILS::ExtractMyValues(*fluiddispnm,mydispnm,gdofs);
    dispnm(0,inode) = mydispnm[0];
    dispnm(1,inode) = mydispnm[1];
    dispnm(2,inode) = mydispnm[2];

    // get the coordinates of patch element and add the current displacement to it
    xyze(0,inode) = pelenodes[inode]->X()[0] + dispnm(0,inode);
    xyze(1,inode) = pelenodes[inode]->X()[1] + dispnm(1,inode);
    xyze(2,inode) = pelenodes[inode]->X()[2] + dispnm(2,inode);
  }

  // check whether the xfemnode is in the element
  LINALG::Matrix<3,DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement > xyzem(xyze,true);
  GEO::CUT::Position <DRT::Element::hex8> pos(xyzem,x);
  bool insideelement = pos.Compute();
  if (insideelement)
  {
    // get the coordinates of x in element coordinates of patch element pele (xsi)
    xsi = pos.LocalCoordinates();

    LINALG::SerialDenseVector shp(numnode);
    // evaluate shape functions
    DRT::UTILS::shape_function_3D(shp,xsi(0,0),xsi(1,0),xsi(2,0),pele->Shape());

    // Interpolate
    LINALG::Matrix<nsd,1> v(true);
    for (std::size_t inode = 0; inode < numnode; ++inode)
      for (std::size_t isd = 0; isd < nsd; ++isd)
        v(isd) += veln(isd,inode)*shp(inode);


    // get the set of Dofsets for the unknown node
    const set<XFEM::FieldEnr> fieldenrset = unknownnodaldof->second;
    const int nodeGid = unknownnodaldof->first;

    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
         fieldenr != fieldenrset.end();
         ++fieldenr)
    {
      const DofKey<onNode> newdofkey(nodeGid, *fieldenr);
      map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(newdofkey);
      const bool found = (newdof != newNodalDofDistrib_.end());
      if (not found) dserror("bug");

      const int newdofpos = newdof->second;
      if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        (*xfluidstate_vector_n)[newdofrowmap_.LID(newdofpos)] = v(0);
      if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        (*xfluidstate_vector_n)[newdofrowmap_.LID(newdofpos)] = v(1);
      if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        (*xfluidstate_vector_n)[newdofrowmap_.LID(newdofpos)] = v(2);
      if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        (*xfluidstate_vector_n)[newdofrowmap_.LID(newdofpos)] = v(3);
    }
    return true;
  }
  else
    return false;
}
//----------------------------------------------------------------
//----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::PackValuestoSend(
  const RCP<Epetra_Vector>                      state_vector,
  std::map<int, set<XFEM::FieldEnr> >           unknownFieldEnr,
  DRT::PackBuffer                            &  sblock
  ) const
{
  // Pack UnknownFieldEnr
  int numentries = (int)unknownFieldEnr.size();

  // add size  to sendblock
  DRT::ParObject::AddtoPack(sblock,numentries);

  int i=0;
  for (std::map<int, set<XFEM::FieldEnr> >::const_iterator unknowndofs = unknownFieldEnr.begin();
       unknowndofs != unknownFieldEnr.end(); unknowndofs++)
  {
    DRT::ParObject::AddtoPack(sblock,unknowndofs->first);
    //DRT::ParObject::AddtoPack(sblock,unknowndofs->second);
    dserror( "not working" );
    ++i;
  }

  if(i!=numentries) dserror("Something wrong with number of elements");

  // Pack state_vector
  double * state_vector_cp  = state_vector->Values();

  // loop all nodes on the processor
  int numstatevectorentries = 0;
  for(int lnodeid=0;lnodeid<ih_->xfemdis()->NumMyRowNodes();lnodeid++)
  {
    // get the processor local node
    DRT::Node*  lnode  = ih_->xfemdis()->lRowNode(lnodeid);
    int numdofs = ih_->xfemdis()->NumDof(lnode);
    numstatevectorentries += numdofs;
  }

  // add size  to sendblock
  DRT::ParObject::AddtoPack(sblock,numstatevectorentries);

  DRT::ParObject::AddtoPack(sblock,state_vector_cp,numstatevectorentries);
}

#ifdef PARALLEL
//----------------------------------------------------------------
//----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::ReceiveBlock(
  vector<char>         & rblock,
  DRT::Exporter        & exporter,
  MPI_Request          & request
) const
{
  // get number of processors and the current processors id
  int numproc = ih_->xfemdis()->Comm().NumProc();
  int myrank  = ih_->xfemdis()->Comm().MyPID();

  // necessary variables
  int         length  = -1;
  int         frompid = (myrank+numproc-1)%numproc;
  int         tag     = frompid;

  // make sure that you do not think you received something if you didn't
  if(rblock.empty() == false)
  {
    dserror("rblock not empty");
  }

  // receive from predecessor
  exporter.ReceiveAny(frompid,tag,rblock,length);

  if(tag != (myrank+numproc-1)%numproc)
  {
    dserror("received wrong message (ReceiveAny)");
  }

  exporter.Wait(request);

  // for safety
  exporter.Comm().Barrier();
}
#endif

#ifdef PARALLEL
//----------------------------------------------------------------
//----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::SendBlock(
  vector<char>         & sblock,
  DRT::Exporter        & exporter,
  MPI_Request          & request
) const
{
  if(sblock.empty() == true)
  {
    dserror("sblock is empty");
  }

  // get number of processors and the current processors id
  int numproc = ih_->xfemdis()->Comm().NumProc();
  int myrank  = ih_->xfemdis()->Comm().MyPID();

  // Send block to next proc.
  int         tag    =myrank;
  int         frompid=myrank;
  int         topid  =(myrank+1)%numproc;

  exporter.ISend(frompid,topid,
                 &(sblock[0]),sblock.size(),
                 tag,request);

  // for safety
  exporter.Comm().Barrier();
}
#endif

//----------------------------------------------------------------
//----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::UnpackValues(
  RCP<Epetra_Vector>                            state_vector,
  std::map<int, set<XFEM::FieldEnr> >           unknownFieldEnr,
  vector<char>                               &  rblock
  ) const
{

  // unpack into state_vector and unknownFieldEnr
  state_vector = null;
  unknownFieldEnr.clear();

}

//---------------------------------------------------------------
//---------------------------------------------------------------
void XFEM::DofDistributionSwitcher::GmshOutput(
  int                               step,
  const RCP<DRT::Discretization>    patchdis,
  std::map<int,int>                 mapforoutput
  ) const
{
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("Interpolation",step, 5, false, patchdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  gmshfilecontent << "View \" " << "xfem Node \" {\n";


  for (int i=0; i<ih_->xfemdis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->xfemdis()->lColNode(i);
    const LINALG::Matrix<3,1> pos(actnode->X());
    map<int,set<XFEM::FieldEnr> >::const_iterator iter = unknownFieldEnr_.find(actnode->Id());
    int unknown = 0;
    if (iter != unknownFieldEnr_.end()) unknown = 1;
    IO::GMSH::cellWithScalarToStream(DRT::Element::point1, unknown, pos, gmshfilecontent);
  }
  gmshfilecontent << "};\n";

  gmshfilecontent << "View \" " << "xfem Node Id\" {\n";
  for (int i=0; i<ih_->xfemdis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->xfemdis()->lColNode(i);
    const LINALG::Matrix<3,1> pos(actnode->X());
    map<int,set<XFEM::FieldEnr> >::const_iterator iter = unknownFieldEnr_.find(actnode->Id());
    int id = -1;
    if (iter != unknownFieldEnr_.end()) id = actnode->Id();
    IO::GMSH::cellWithScalarToStream(DRT::Element::point1, id, pos, gmshfilecontent);
  }
  gmshfilecontent << "};\n";

  gmshfilecontent << "View \" " << "patch element \" {\n";
  // get the patch-element for the box
  {
    for (int i=0; i<patchdis->NumMyColElements(); ++i)
    {
      const DRT::Element* pele = patchdis->lColElement(i);
      std::map<int,int>::const_iterator patchiter =mapforoutput.find(pele->Id());
      int inside;
      if (patchiter != mapforoutput.end()) inside = patchiter->second;
      const std::size_t numnode = pele->NumNode();
      LINALG::SerialDenseMatrix xyze(3,numnode);
      const DRT::Node*const* pelenodes = pele->Nodes();

      for (std::size_t inode = 0; inode < numnode; ++inode)
      {
        // get the coordinates of patch element
        xyze(0,inode) = pelenodes[inode]->X()[0];
        xyze(1,inode) = pelenodes[inode]->X()[1];
        xyze(2,inode) = pelenodes[inode]->X()[2];
      }
      IO::GMSH::cellWithScalarToStream(pele->Shape(),inside,xyze, gmshfilecontent);
    }
  }
  gmshfilecontent << "};\n";

  gmshfilecontent << "View \" " << "patch element id \" {\n";
  // get the patch-element for the box
  {
    for (int i=0; i<patchdis->NumMyColElements(); ++i)
    {
      const DRT::Element* pele = patchdis->lColElement(i);
      const std::size_t numnode = pele->NumNode();
      LINALG::SerialDenseMatrix xyze(3,numnode);
      const DRT::Node*const* pelenodes = pele->Nodes();

      for (std::size_t inode = 0; inode < numnode; ++inode)
      {
        // get the coordinates of patch element
        xyze(0,inode) = pelenodes[inode]->X()[0];
        xyze(1,inode) = pelenodes[inode]->X()[1];
        xyze(2,inode) = pelenodes[inode]->X()[2];
      }
      IO::GMSH::cellWithScalarToStream(pele->Shape(),pele->Id(),xyze, gmshfilecontent);
    }
  }
  gmshfilecontent << "};\n";

  gmshfilecontent.close();
}
//----------------------------------------------------------------
//----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::mapVectorToNewDofDistributionCombust(
    RCP<Epetra_Vector>&    vector,
    const bool             quasi_static_enr
) const
{
  // create new vector with new number of dofs
  const RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);

  if (vector == null)
  {
#ifdef DEBUG
    std::cout << "  created new vector with all zeros" << endl;
#endif
  }
  else
  {
#ifdef DEBUG
//     bool completely_unchanged = true;
#endif
    const RCP<Epetra_Vector> oldVector = vector;
    const Epetra_BlockMap& oldmap = oldVector->Map();
//    std::cout << "olddofrowmap_" << endl;
//    std::cout << (olddofrowmap_) << endl;
//    std::cout << "newdofrowmap_" << endl;
//    std::cout << (newdofrowmap_) << endl;

    if (not oldmap.SameAs(olddofrowmap_)) dserror("bug!");

    // step 1: find predecessor of new nodal dofkey
    for (map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.begin();
                                                     newdof != newNodalDofDistrib_.end();
                                                     ++newdof)
    {
      const DofKey<onNode> newdofkey = newdof->first;
      const int newdofpos = newdof->second;

      map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
      if (olddof != oldNodalDofDistrib_.end()) // if dofkey has existed before, use old value
      {
        const DofKey<onNode> olddofkey = olddof->first;
        const int olddofpos = olddof->second;
        //cout << newdofkey.toString() << " -> init to old value" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];

        if (quasi_static_enr == true)
        {
          //std::cout << "-------------Warning: enriched dofs reset to zero------------" << std::endl;
          if (newdofkey.getFieldEnr().getEnrichment().Type() != XFEM::Enrichment::typeStandard)
            (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
//        if (newdofkey.getFieldEnr().getEnrichment().Type() == XFEM::Enrichment::typeKink)
//            (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
        }
      }
      else // if dofkey has not been existed before, check for other dofs on the dofkeys node
      {
        //const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();
        // initialize to zero
        (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
      }
    }

//    // step 2: find sucessor of old nodal dofkey to sum up values
//    for (map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.begin();
//                                                     olddof != oldNodalDofDistrib_.end();
//                                                     ++olddof)
//    {
//      const DofKey<XFEM::onNode> olddofkey = olddof->first;
//      const int olddofpos = olddof->second;
//      const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();
//
//      // try to find successor
//      map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
//      if (newdof == newNodalDofDistrib_.end())  // if no successor found (was handled already in step 1)
//      {
//        // try to find another usefull value
//        // current assumption: there is only one type of enrichment per node
//        // no overlapping enrichments allowed for now
//        const int nodegid = olddofkey.getGid();
//        const XFEM::Enrichment oldenr(olddofkey.getFieldEnr().getEnrichment());
//
//        // create alternative dofkey
//        XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, oldphysvar, dofman_));
//
//        if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
//        {
//          // find dof position of alternative key
//          const XFEM::FieldEnr altfieldenr(olddofkey.getFieldEnr().getField(), altenr);
//          const DofKey<onNode> altdofkey(nodegid, altfieldenr);
//          const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;
//
//          //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
//          if (newdofpos < 0)
//          {
//            std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
//            std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
//            dserror("bug!");
//          }
//
//          // add old value to already existing values
//          //(*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofrowmap_.LID(olddofpos)];
//          (*newVector)[newdofrowmap_.LID(newdofpos)] += (*oldVector)[olddofrowmap_.LID(olddofpos)];
//          completely_unchanged = false;
//        }
//        else // if not alternative is found
//        {
//          // this can only happen in the void enrichment case and in that case,
//          // the dof value is zero anyway, which coincides with the fact that we have no place,
//          // where we could store it ;-)
//        }
//      }
//      else
//      {
//        // do nothing, this case was handled in step 1
//      }
//    }

#if 0
    // remark: Only if DLM condensation is turned off, there will be element dofs visible in the global
    //         DofManager. Of course, irrespective of the DLM condensation setting, there will be
    //         element dofs on the element level (ElementDofManager), if a stress-based approach is
    //         chosen to handle the boundary terms.
    //         So, if there are element dofs, they are not known here, since 'elementalDofs_' in
    //         fillDofRowDistributionMaps() only knows the different field enrichments for a specific
    //         element, that is
    //
    //         e.g.
    //         element 103 has 28 element dofs ((6 x stess + 1 x disc pressure) x element ansatz (e.g. quad4))
    //
    //         but 'elementalDofs_' contains only 7 "dof-type entries"
    //         <ele 103, Tauxx + quad4>
    //         <ele 103, Tauxy + quad4>
    //         ...
    //         <ele 103, DiscPres + quad4>
    //
    //         That is, element dofs cannot be switched correctly here, since not all element dofs
    //         are known. They only exist on the element level. They are generated in
    //         ElementDofManager::ComputeDependentInfo()
    dserror("This will not work properly! Read comment above!");

    // step 3: find predecessor of new elemental dofkey
    for (map<DofKey<onElem>, DofGID>::const_iterator newdof = newElementalDofDistrib_.begin();
                                                     newdof != newElementalDofDistrib_.end();
                                                     ++newdof)
    {
      const DofKey<XFEM::onElem> newdofkey = newdof->first;
      const int newdofpos = newdof->second;

      map<DofKey<onElem>, DofGID>::const_iterator olddof = oldElementalDofDistrib_.find(newdofkey);
      if (olddof != oldElementalDofDistrib_.end())  // if dofkey has existed before, use old value
      {
        const DofKey<XFEM::onElem> olddofkey = olddof->first;
        const int olddofpos = olddof->second;
        //cout << "init to old value" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
      }
      else // if dofkey has not been existed before, initialize to zero
      {
        //cout << "init to zero" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
        completely_unchanged = false;
      }
    }
#endif

#if 0
    if (completely_unchanged)
      cout << "completely unchanged vector" << endl;
    else
      cout << "modified vector" << endl;
#endif
  }


  // set vector to zero or initialized vector
  vector = newVector;
}

//----------------------------------------------------------------
//----------------------------------------------------------------
void XFEM::DofDistributionSwitcher::mapVectorToNewDofDistributionFluidXFluid(
  RCP<Epetra_Vector>&    vector,
  std::map<int,int>      nodelabeln
) const
{
  // create new vector with new number of dofs
  const RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);

  if (vector == null)
  {
#ifdef DEBUG
    std::cout << "  created new vector with all zeros" << endl;
#endif
  }
  else
  {
#ifdef DEBUG
    bool completely_unchanged = true;
#endif
    const RCP<Epetra_Vector> oldVector = vector;
    const Epetra_BlockMap& oldmap = oldVector->Map();
//    std::cout << "olddofrowmap_" << endl;
//    std::cout << (olddofrowmap_) << endl;
//    std::cout << "newdofrowmap_" << endl;
//    std::cout << (newdofrowmap_) << endl;


    if (not oldmap.SameAs(olddofrowmap_)) dserror("bug!");

    // step 1: find predecessor of new nodal dofkey
    for (map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.begin();
         newdof != newNodalDofDistrib_.end(); ++newdof)
    {
      const DofKey<onNode> newdofkey = newdof->first;
      const int newdofpos = newdof->second;

      map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
      if (olddof != oldNodalDofDistrib_.end()) // if dofkey has existed before, use old value
      {
        const DofKey<onNode> olddofkey = olddof->first;
        const int olddofpos = olddof->second;
        //cout << newdofkey.toString() << " -> init to old value" << endl;
        (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
      }
      else // if dofkey has not been existed before, check for other dofs on the dofkeys node
      {
        //const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();

        // initialize to zero
        (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
      }
    }

   // step 2: find successor of old nodal dofkey to sum up values
   for (map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.begin();
        olddof != oldNodalDofDistrib_.end(); ++olddof)
   {
      const DofKey<XFEM::onNode> olddofkey = olddof->first;
      const int olddofpos = olddof->second;
      const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();

      // was in xfluid? yes=0
      map<int,int>::const_iterator iter = nodelabeln.find(olddofkey.getGid());

      // try to find successor
      map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);

      // if no successor found (was handled already in step 1) and if node was
      // in xfluid in the last time step
      if (newdof == newNodalDofDistrib_.end() and iter->second != 1)
      {
        // try to find another usefull value
        // current assumption: there is only one type of enrichment per node
        // no overlapping enrichments allowed for now
        const int nodegid = olddofkey.getGid();
        const XFEM::Enrichment oldenr(olddofkey.getFieldEnr().getEnrichment());

        // create alternative dofkey
        XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, oldphysvar, dofman_));

        if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
        {
          // find dof position of alternative key
          const XFEM::FieldEnr altfieldenr(olddofkey.getFieldEnr().getField(), altenr);
          const DofKey<onNode> altdofkey(nodegid, altfieldenr);
          const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;

          //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
          if (newdofpos < 0)
          {
            std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
            std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
            dserror("bug!");
          }

          // add old value to already existing values
          //(*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofrowmap_.LID(olddofpos)];
          (*newVector)[newdofrowmap_.LID(newdofpos)] += (*oldVector)[olddofrowmap_.LID(olddofpos)];
#ifdef DEBUG
          completely_unchanged = false;
#endif
        }
        else // if not alternative is found
        {
          // this can only happen in the void enrichment case and in that case,
          // the dof value is zero anyway, which coincides with the fact that we have no place,
          // where we could store it ;-)
        }
      }
      else
      {
        // do nothing, this case was handled in step 1
      }
    }
  }
  // set vector to zero or initialized vector
  vector = newVector;
}
////////////////////////////////////////////
void XFEM::DofDistributionSwitcher::generateTransferInformation(
    const RCP<Epetra_Vector>&             vector
) const
{

  dserror("not ready, yet!");
  // create new vector with new number of dofs
  RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);

  if (vector == null)
  {
#ifdef DEBUG
    std::cout << "  created new vector with all zeros" << endl;
#endif
  }
  else
  {
    bool completely_unchanged = true;
    const RCP<Epetra_Vector> oldVector = vector;
    const Epetra_BlockMap& oldmap = oldVector->Map();

    if (not oldmap.SameAs(olddofrowmap_))
      dserror("bug!");

    map<DofKey<onNode>,set<DofKey<onNode> > > transferop;

    // step 1: find predecessor of new nodal dofkey
    for (map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.begin();
    newdof != newNodalDofDistrib_.end();
    ++newdof)
    {
      const DofKey<onNode> newdofkey = newdof->first;
      //            const int newdofpos = newdof->second;
      if (newdofkey.getFieldEnr().getField() == XFEM::PHYSICS::Velx)
      {
        transferop[newdofkey] = set<DofKey<onNode> >();
      }
      map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
      if (olddof != oldNodalDofDistrib_.end())  // if dofkey has existed before, use old value
      {
        const DofKey<onNode> olddofkey = olddof->first;
        //                const int olddofpos = olddof->second;
        //cout << newdofkey.toString() << " -> init to old value" << endl;
        //                (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];

        if (newdofkey.getFieldEnr().getField() == XFEM::PHYSICS::Velx)
        {
          transferop[newdofkey].insert(olddofkey);
        }


      }
      else // if dofkey has not been existed before, check for other dofs on the dofkeys node
      {
        //              const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();

        // initialize to zero
        //                (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
      }
    }

    // step 2: find sucessor of old nodal dofkey to sum up values
    for (map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofDistrib_.begin();
    olddof != oldNodalDofDistrib_.end();
    ++olddof)
    {
      const DofKey<onNode> olddofkey = olddof->first;
      const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();

      // try to find successor
      map<DofKey<onNode>, DofGID>::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
      if (newdof == newNodalDofDistrib_.end())  // if no successor found
      {
        // try to find another usefull value
        // current assumption: there is only one type of enrichment per node
        // no overlapping enrichments allowed for now
        const int nodegid = olddofkey.getGid();
        const XFEM::Enrichment oldenr = olddofkey.getFieldEnr().getEnrichment();

        // create alternative dofkey
        XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, oldphysvar, dofman_));

        if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
        {
          // find dof position of alternative key
          const XFEM::FieldEnr altfieldenr(olddofkey.getFieldEnr().getField(), altenr);
          const DofKey<onNode> altdofkey(nodegid, altfieldenr);
          const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;

          //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
          if (newdofpos < 0)
          {
            std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
            std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
            dserror("bug!");
          }

          // add old value to already existing values
          //(*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofrowmap_.LID(olddofpos)];
          //                    (*newVector)[newdofrowmap_.LID(newdofpos)] += (*oldVector)[olddofrowmap_.LID(olddofpos)];
          completely_unchanged = false;
          if (altdofkey.getFieldEnr().getField() == XFEM::PHYSICS::Velx)
          {
            transferop[altdofkey].insert(olddofkey);
          }
        }
        else // if not alternative is found
        {
          // this can only happen in the void enrichment case and in that case,
          // the dof value is zero anyway, which coincides with the fact that we have no place,
          // where we could store it ;-)
        }
      }
      else
      {
        // do nothing, this case was handled in step 1
      }
    }

#if 0
    if (completely_unchanged)
      cout << "completely unchanged vector" << endl;
    else
      cout << "modified vector" << endl;
#endif
    map<DofKey<onNode>,set<DofKey<onNode> > >::const_iterator entry;
    for (entry = transferop.begin(); entry!=transferop.end(); ++entry)
    {
      if (entry->second.size() > 1)
        dserror("only one successor allowed for now");

      if (entry->second.empty())
      {
        cout << RED << entry->first << END_COLOR << endl;
      }
      else
      {
        set<DofKey<onNode> >::const_iterator first_set_entry = entry->second.begin();
        DofKey<onNode> old = *entry->second.begin();
        if (entry->first != old)
        {
          cout << GREEN << entry->first << END_COLOR << " <-- " << GREEN << old << END_COLOR << endl;
        }
      }
    }
    //exit(1);
  }
}



//! try to find another enrichment for this physical field
XFEM::Enrichment XFEM::genAlternativeEnrichment(
    const int                    gnodeid,
    const XFEM::PHYSICS::Field   oldphysvar,
    const RCP<XFEM::DofManager>& dofman
)
{
  const std::set<XFEM::FieldEnr>& fieldset(dofman->getNodeDofSet(gnodeid));
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenriter = fieldset.begin(); fieldenriter != fieldset.end(); ++fieldenriter)
  {
    const XFEM::PHYSICS::Field physvar = fieldenriter->getField();
    if (oldphysvar == physvar)
    {
      return fieldenriter->getEnrichment();
      break;
    }
  }
  return XFEM::Enrichment();
}

#endif  // #ifdef CCADISCRET
