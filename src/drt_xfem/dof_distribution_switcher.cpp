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
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_io/io_gmsh.H"
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
    bool completely_unchanged = true;
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
          completely_unchanged = false;
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
#if 1

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

#ifdef DEBUG
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


void XFEM::DofDistributionSwitcher::extrapolateOldTimeStepValues(
    const RCP<DRT::Discretization>    bdis,
    const std::map<int,LINALG::Matrix<3,1> >&  cutterposn,
    const RCP<const Epetra_Vector>    ivector,
    const RCP<Epetra_Vector>          state_vector
) const
{

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









void XFEM::DofDistributionSwitcher::mapVectorToNewDofDistributionCombust(
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
    bool completely_unchanged = true;
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
#if 1

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

#ifdef DEBUG
    if (completely_unchanged)
      cout << "completely unchanged vector" << endl;
    else
      cout << "modified vector" << endl;
#endif
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

#ifdef DEBUG
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
