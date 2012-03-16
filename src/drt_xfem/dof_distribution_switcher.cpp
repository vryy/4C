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
#include "../drt_combust/combust_interface.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_element.H"
#include <iostream>


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


  for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->FluidDis()->lColNode(i);
    const LINALG::Matrix<3,1> pos(actnode->X());
    map<int,set<XFEM::FieldEnr> >::const_iterator iter = unknownFieldEnr_.find(actnode->Id());
    int unknown = 0;
    if (iter != unknownFieldEnr_.end()) unknown = 1;
    IO::GMSH::cellWithScalarToStream(DRT::Element::point1, unknown, pos, gmshfilecontent);
  }
  gmshfilecontent << "};\n";

  gmshfilecontent << "View \" " << "xfem Node Id\" {\n";
  for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->FluidDis()->lColNode(i);
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
      int inside = -1;
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
